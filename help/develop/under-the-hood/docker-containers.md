---
title: "Docker containers (under the hood)"
sidebar_label: "Docker containers"
keywords:
  - grok_spawner
  - grok_registry
  - grok_registry_proxy
  - docker
  - plugin container
  - ECR
---

This page describes how Datagrok manages plugin Docker images and containers end-to-end:
how a package's `Dockerfile` becomes a running container, what services participate, and
what an administrator has to set up on each deployment path.

For the developer-facing workflow (how to add a `Dockerfile` to a package and call it from
plugin code), see [Creating a docker container](../how-to/packages/docker-containers.md).

## Components

| Service                | Role                                                                                                  |
|------------------------|-------------------------------------------------------------------------------------------------------|
| `grok` CLI             | Builds the image on the developer's workstation and pushes it to the registry during `grok publish`.  |
| Datlas                 | Owns the `DockerImage` and `DockerContainer` entities. Drives the validation and deployment loop.     |
| `grok_spawner`         | Talks to the orchestrator (Docker, Docker Swarm, ECS, Kubernetes), validates images, runs containers. |
| `grok_registry`        | Self-managed Docker registry for on-prem stands. Datagrok-JWT-authenticated.                          |
| `grok_registry_proxy`  | Optional JWT-auth proxy in front of a backing registry (ECR, Docker Hub, any V2 registry).            |
| Backing registry       | Where image layers actually live: ECR, Docker Hub, the internal `grok_registry`, etc.                 |

`grok_spawner` is the only service that talks to the orchestrator. Datlas never calls Docker
or Kubernetes directly; it only calls `grok_spawner` over the spawner REST API.

See [Infrastructure](infrastructure.md) and [Deployment](../../deploy/deploy.md) for where
these services sit in the broader platform.

## End-to-end flow

```
Developer workstation                            Datagrok server
─────────────────────                            ───────────────
grok publish
 1. docker build --platform linux/amd64
 2. docker tag {package}-{container}:{version}
 3. docker login registry.{domain}      ──►  grok_registry_proxy
       (JWT via /api/docker/token)               ▼
 4. docker push                         ──►  ECR / grok_registry / Docker Hub
 5. POST package ZIP                    ──►  Datlas
                                                 │
                                       creates DockerImage row (PENDING_VALIDATION)
                                                 │
                                                 ▼
                                       Datlas timer (processImages)
                                       POST /validate/cvm/{name}/{tag}  ──►  grok_spawner
                                                                          │
                                                              checks ECR / Docker Hub
                                                                          │
                                       image status = READY  ◄────────────┘
                                                 │
                       User toggles "Start"  ────┤
                                                 ▼
                                       Datlas timer (processContainers)
                                       POST /deploy/cvm/{service}/1  ──►  grok_spawner
                                                                       │
                                                          orchestrator (Docker/ECS/K8s)
                                                                       │
                                       container running, address known
                                                 │
Plugin JS code                                   │
fetchProxy(id, '/path')                ──►  Datlas /docker/containers/proxy/{id}/{path}
                                                 │
                                                 ▼
                                       Running container on internal network
```

## Image lifecycle

The `DockerImage` entity tracks one image (one `Dockerfile` from one package). Its lifecycle
is driven by a Datlas timer (`processImages` in `docker_service.dart`) calling
`grok_spawner`'s `POST /validate/{host}/{service}/{tag}` endpoint.

```
PENDING_VALIDATION ──► VALIDATING ──► READY
                                  └─► ERROR
```

| Status               | Meaning                                                                       |
|----------------------|-------------------------------------------------------------------------------|
| `PENDING_VALIDATION` | New image (just pushed) or rebuild requested. Waiting for the next timer tick. |
| `VALIDATING`         | Spawner is checking that the image exists and is pullable.                    |
| `READY`              | Image confirmed in the registry. Containers can now start.                    |
| `ERROR`              | Image not found, registry unreachable, or the image manifest is invalid.      |

Validation checks ECR first when the spawner is running in ECS mode (the typical AWS
deployment), then falls back to Docker Hub. The full URI is resolved at this point —
on ECS, a short tag like `datagrok/bio:2.26.6` is rewritten to the cluster's
`<account>.dkr.ecr.<region>.amazonaws.com/datagrok/bio:2.26.6` using
`ecr:GetAuthorizationToken`.

To re-validate an image after a manual registry change (for example, retagging in ECR
without re-publishing the package), use `grok.dapi.docker.dockerImages.revalidate(id)`.

## Container lifecycle

A `DockerContainer` is one running instance of a `DockerImage`. Each container has a
`desired_instances` counter (0 or 1 in practice). Changing it triggers the next iteration
of Datlas's `processContainers` timer, which calls `grok_spawner` to deploy or destroy.

```
                  ┌──► PENDING_START ──► STARTING ──► STARTED
(desired=1)  ─────┤                                       │
                  └──► ERROR                              │
                                                          │
(desired=0)  ─────┬──► PENDING_USER_STOP ──► STOPPING ──► USER_STOPPED
                  └──► PENDING_SYSTEM_STOP ──► STOPPING ─► SYSTEM_STOPPED

(any state stuck > containerStatusTimeoutMinutes)  ──► ERROR
```

Once `STARTED`, Datlas caches the container's address and port (and DNS name on ECS / K8s)
and starts proxying requests from plugin code. Containers configured with
`shutdown_timeout` move to `SYSTEM_STOPPED` automatically after the idle period and are
re-started on the next request.

### Environment variables injected at deploy time

Every container is started with these variables injected by Datlas before the deploy call:

| Variable | Source |
|----------|--------|
| `DATAGROK_AMPQ_HOST` / `_PORT` / `_USER` / `_PASSWORD` / `_TLS` | Queue settings (call queue / AMQP). |
| `DATAGROK_PIPE_HOST` / `_PORT` / `_KEY` | Pipe settings (binary streaming). |
| `DATAGROK_API_URL` | Server API root, so the container can call back into Datagrok. |
| `USER_API_KEY` | Session token of the user whose action started the container. |
| Custom credentials | Encrypted entity credentials and `env` entries from `container.json`. |

Anything in `container.json`'s `env` map is appended last and can reference entities via
`#{x.Package:Entity}` expressions (rendered server-side). Credentials are only
resolvable for entities owned by the same package.

## How images get to the registry

`grok publish` does both the image build and the registry push. The push goes through
the Datagrok JWT auth layer rather than direct registry credentials — developers never
hold AWS keys or registry passwords.

```
docker login registry.{domain} -u {dev-key} -p {dev-key}
  ▼
grok_registry_proxy returns 401 WWW-Authenticate (realm = /token)
  ▼
docker client calls /token with Basic auth (dev-key)
  ▼
proxy forwards to Datagrok /api/docker/token
  ▼
Datagrok validates dev-key, returns RS256 JWT (30 min, claim lim=docker)
  ▼
docker client retries with Bearer JWT
  ▼
proxy validates JWT signature/issuer/expiry/scope
  ▼
proxy attaches ECR or registry creds, forwards to backing registry
```

The JWT signing key is rotated on the Datagrok server. The proxy fetches the current
public key from `GET /api/docker/keys` at startup and re-fetches when validation
returns an unknown `kid`.

Read operations (GET / HEAD on `/v2/...`) pass through without auth. Only writes
(PUT / POST / PATCH / DELETE) require a JWT. ECR repositories are auto-created on first
push when the proxy is in ECR mode.

## Orchestration backends

`grok_spawner` selects a backend from its configuration. One spawner instance can only
target one backend at a time.

| Backend         | Selected by                                                              | Used on                                |
|-----------------|--------------------------------------------------------------------------|----------------------------------------|
| Docker (host)   | `GROK_SPAWNER_FORCE_DOCKER_TYPE_HOST=true` or a non-Swarm Docker socket. | Compose, single-VM, local development. |
| Docker Swarm    | Docker socket reports a Swarm manager and `FORCE_DOCKER_TYPE_HOST=false`. | Swarm clusters (rare).                |
| AWS ECS         | `DATAGROK_ECS__CLUSTER` set.                                              | AWS CloudFormation (ECS / EKS).        |
| Kubernetes      | `kubernetes` block configured.                                            | Helm chart / EKS / GKE / on-prem K8s. |

### Docker (host)

Runs containers directly via the Docker socket. Simplest setup; assumes the spawner
container has `/var/run/docker.sock` mounted. Containers attach to the Docker bridge
network shared with Datlas so the proxy can reach them by name.

### Docker Swarm

Uses `docker.services.create()` for replicated services on overlay networks. Healthchecks
piggy-back on task state. Node placement constraints map to Swarm node role labels.

### AWS ECS (Fargate / EC2)

The primary cloud backend. Key behaviours:

* **Launch type** — `allocate_instance_or_fargate()` checks EC2 capacity first and falls
  back to Fargate if there are no slots. Set `DATAGROK_ECS__LAUNCH_TYPE=FARGATE` to skip
  the EC2 check.
* **Fargate resource buckets** — requested CPU / memory are mapped to the nearest valid
  Fargate combination (0.25–16 vCPU, 0.5–120 GB). Don't expect `cpu: 3.7` to land
  literally.
* **Image resolution** — ECR-first, Docker Hub fallback. The base ECR URI is discovered
  from the IAM task role via `ecr:GetAuthorizationToken`.
* **Service deployment** — ECS services with the deployment circuit breaker enabled.
  The spawner waits for rollout completion before returning success.
* **Logging** — all plugin containers log to the shared CloudWatch group set by
  `DATAGROK_ECS__LOG_GROUP`.

### Kubernetes

Creates a `Deployment` plus `Service` per plugin container in the configured namespace.
Service DNS follows `{service_name}.{namespace}.svc.cluster.local`. Pulls from private
registries use the `imagePullSecrets` referenced in the spawner's `kubernetes.dockerAuth`
config.

## Registry setup

You only need a registry if the deployment hosts custom plugin images. A stand that
only runs Docker-Hub-public images works without `grok_registry` and without
`grok_registry_proxy`.

### grok_registry (on-prem, self-managed)

Bundles Docker Distribution Registry + Nginx + `joxit/docker-registry-ui` in a single
container. JWT authentication is delegated to Datagrok.

Required configuration:

| Variable        | Description                                          |
|-----------------|------------------------------------------------------|
| `DATAGROK_URL`  | Datagrok API URL **with `/api` suffix**.             |

Ports:

| Port | Purpose                            |
|------|------------------------------------|
| 5000 | Docker V2 API (login / push / pull). |
| 8080 | Web UI (browse pushed images).     |

Storage is a `/var/lib/registry` filesystem volume — size it for image growth (each
image version is kept until manually pruned).

### grok_registry_proxy (cloud, backed by ECR / Hub / another registry)

A small OpenResty/Nginx image that adds the same Datagrok JWT auth in front of an
external registry. It does not store image layers — those live in the backing registry.

| Variable             | Required        | Description                                                    |
|----------------------|-----------------|----------------------------------------------------------------|
| `DATAGROK_URL`       | yes             | Datagrok API URL (with `/api`).                                |
| `REGISTRY_URL`       | yes             | Backing registry endpoint (ECR repository host, Docker Hub, …). |
| `ECR_MODE`           | no              | `true` enables ECR IAM auth via `aws ecr get-login-password`.  |
| `AWS_REGION`         | ECR mode        | Region for ECR.                                                |
| `DATAGROK_ISSUER`    | no              | Overrides JWT issuer; auto-detected from `/api/admin/startup_data`. |
| `REGISTRY_USERNAME`  | non-ECR mode    | Username for the backing registry.                             |
| `REGISTRY_PASSWORD`  | non-ECR mode    | Password for the backing registry.                             |

In ECR mode the proxy:

1. Reads ECR creds with the ECS task role (refreshes every 6 hours; ECR tokens expire at 12).
2. Auto-creates the ECR repository on first push (`ecs:CreateRepository`).
3. Caches the bearer token on disk for the worker Lua processes.

## Datagrok configuration

The Datagrok side of the wiring lives in `PluginsDockerSettings` (server config). The
default values are correct for the Docker Compose dev stack; cloud deployments override
a few of them via Helm values or the CloudFormation parameters.

| Setting                          | Default            | Description                                                                                 |
|----------------------------------|--------------------|---------------------------------------------------------------------------------------------|
| `useGrokSpawner`                 | `true`             | Master toggle for plugin containers. Set false to disable everything.                       |
| `grokSpawnerHost`                | `grok_spawner`     | DNS name / hostname for the spawner.                                                        |
| `grokSpawnerPort`                | `8000`             | Spawner HTTP port.                                                                          |
| `grokSpawnerApiKey`              | `test-x-api-key`   | Shared secret sent in `X-Api-Key`. Always override in prod.                                 |
| `containerStatusTimeoutMinutes`  | `5`                | Max time a container can spend in a transitional state before being marked `ERROR`.         |
| `proxyRequestTimeout`            | `300000` ms        | Timeout for `fetchProxy` / `webSocketProxy`. Must exceed orchestrator cold-start latency.   |
| `registryProxyHost`              | (empty)            | Public hostname of `grok_registry_proxy`.                                                   |
| `registryProxyPort`              | `5000`             | Port for `docker push` / `docker pull`.                                                     |

The `proxyRequestTimeout` default of 5 minutes matches the ECS Fargate cold-start budget.
If you also front Datagrok with a load balancer, raise that load balancer's idle timeout
to match (the AWS internal ALB ships with 300 s for this reason).

## Setup per deployment path

### Docker Compose (single-machine / development)

Bring up the spawner alongside Datlas with the `grok_spawner` profile. The compose file
mounts the Docker socket, so containers are sibling to the spawner on the host.

```bash
cd core/server/datlas/resources/compose
docker compose -f deps.docker-compose.yaml -p datagrok-dev \
  --profile grok_spawner up -d
```

This stand does **not** need a registry — every plugin container is built against
Docker Hub. To host custom images locally, also start `grok_registry`:

```bash
docker compose -f deps.docker-compose.yaml -p datagrok-dev \
  --profile grok_registry up -d
```

…and point Datagrok at it with `registryProxyHost=grok_registry`.

See [Docker Compose deployment](../../deploy/docker-compose/docker-compose.mdx) and the
advanced variant for the full compose configuration.

### Kubernetes / Helm

The [Helm chart](../../deploy/k8s/install-helm-chart.md) ships `spawner.*` and
`grokRegistryProxy.*` value blocks. Both are enabled by default and pull image tags
that match the chart version.

Required pieces:

* `spawner.image.tag` — pin to the version compatible with your Datagrok core (see
  [Images and versions](../../deploy/images.md)).
* `grokRegistryProxy.image.tag` — same pinning.
* `grokRegistryProxy.env.REGISTRY_URL` — backing registry. On EKS this is the ECR base
  URI; on other clusters typically a Docker Hub org or a self-hosted V2 registry.
* `serviceAccount.annotations` — if the backing registry is ECR, the Datagrok service
  account needs IRSA with `AmazonEC2ContainerRegistryFullAccess` (or a narrower
  read/write/CreateRepository policy).

The chart auto-grants `grok_spawner` the in-cluster permissions it needs to create
Deployments and Services in the same namespace.

### AWS CloudFormation (ECS / EKS)

Both [CloudFormation templates](../../deploy/aws/deploy-amazon-eks.mdx) declare:

| Resource                          | What it does                                                                                  |
|-----------------------------------|-----------------------------------------------------------------------------------------------|
| `GrokSpawnerService`              | Spawner Fargate service on the internal ALB, port 8000.                                       |
| `GrokRegistryProxyService`        | Registry proxy Fargate service on the internal ALB (5000) and external ALB via host-header.  |
| `DatagrokTCP443Listener` + rule   | Routes `registry.<host>` traffic on the external ALB to the proxy.                            |
| `GrokRegistryProxyECRPolicy`      | IAM policy for `ecr:Get*`, `ecr:Put*`, `ecr:CreateRepository`, `ecr:GetAuthorizationToken`.   |
| `EnableRegistryProxy` parameter   | Set to `false` to skip the proxy entirely (no custom plugin images).                          |

The internal ALB's idle timeout is set to 300 s to absorb Fargate cold starts. Don't
lower it — the spawner waits on `runTask` synchronously and a shorter timeout returns
spurious `ERROR` states.

### Bare metal / on-prem VM

Run `grok_spawner` on the same host as Datlas, mount the Docker socket, and start
`grok_registry` if you need to host custom images. The spawner is happy with a plain
`docker run`:

```bash
docker run -d --name grok_spawner \
  --user root \
  --restart unless-stopped \
  -e X_API_KEY={shared-secret} \
  -e GROK_SPAWNER_ENVIRONMENT=onprem \
  -e GROK_SPAWNER_FORCE_DOCKER_TYPE_HOST=true \
  -p 8000:8000 \
  -v /var/run/docker.sock:/var/run/docker.sock \
  --network datagrok \
  datagrok/grok_spawner:{tag}
```

Match `X_API_KEY` to Datlas's `grokSpawnerApiKey`, and the network to Datlas's network.

## DNS, ports, and connectivity

Plugin containers are only reachable through Datlas's proxy — clients never talk to a
plugin container directly. The connectivity requirements are therefore:

* Datlas → `grok_spawner` (port 8000) — REST API for deploy / destroy / status.
* `grok_spawner` → orchestrator (Docker socket / ECS API / Kubernetes API).
* `grok_spawner` → plugin containers (HTTP, port from `EXPOSE`) for healthchecks and
  proxying.
* Datlas → plugin containers (same path; the spawner returns address + port and Datlas
  caches them).
* Developer workstation → `registry.<host>` (HTTPS) — only at publish time.

A single `EXPOSE PORT` per Dockerfile is supported. Containers exposing multiple ports
will not work — the spawner picks the first declared port.

## Logs and troubleshooting

| Want to see…                       | Where to look                                                                     |
|------------------------------------|-----------------------------------------------------------------------------------|
| Image validation errors            | **Manage → Dockers** → image card → context panel → **Build logs**.               |
| Container runtime logs             | **Manage → Dockers** → container card → context panel → **Logs**.                 |
| Spawner-side decisions             | `grok_spawner` container logs (deploy / destroy / validate request traces).        |
| Datlas state machine               | Datlas logs, search for `DockerService` / `processImages` / `processContainers`.   |
| Registry auth failures             | `grok_registry_proxy` logs (look for the `/token` and `/v2/` access lines).        |
| ECR pull errors on ECS             | The ECS task event stream — `ResourceInitializationError` is almost always IAM.    |

Common failures and their first checks:

* **Image stuck in `PENDING_VALIDATION`** — spawner not reachable, or `grokSpawnerApiKey`
  mismatch (HTTP 401).
* **Image goes `ERROR` immediately after publish** — registry proxy can't pull from the
  backing registry. On ECR, check IAM; on Docker Hub, check rate limits.
* **Container goes `ERROR` after start** — open the container logs first; the application
  inside the image is the usual culprit. If logs are empty, the orchestrator killed the
  task before it bound the port — bump `memory` in `container.json`.
* **`fetchProxy` returns 504** — `proxyRequestTimeout` is too low or the container is
  cold-starting. Raise the timeout, or pre-warm the container by removing
  `shutdown_timeout`.

## See also

* [Creating a docker container](../how-to/packages/docker-containers.md) — developer how-to.
* [Deployment](../../deploy/deploy.md) — deployment paths and components.
* [Images and versions](../../deploy/images.md) — image tag conventions for `grok_spawner`
  and `grok_registry_proxy`.
* [Install Datagrok with Helm](../../deploy/k8s/install-helm-chart.md) — chart values
  for spawner and registry proxy.
* [AWS CloudFormation (EKS)](../../deploy/aws/deploy-amazon-eks.mdx) — full stack
  including spawner, registry proxy, IAM, and ALB rules.
* [Infrastructure](infrastructure.md) — where `grok_spawner` sits in the platform.
