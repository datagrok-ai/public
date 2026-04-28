---
title: "Images and versions"
sidebar_position: 12
---

The Datagrok services that run on every deployment path
([Components](deploy.md#components)) are pulled from
[Docker Hub](https://hub.docker.com/u/datagrok). Defaults below are pinned to the latest
stable Datagrok release and are tested together. Override only when you need to pin a
specific release or test a candidate.

## Latest defaults

| Service | Image | Version |
|---|---|---|
| Datagrok | [`datagrok/datagrok`](https://hub.docker.com/r/datagrok/datagrok) | `1.27.3` |
| Grok Pipe | [`datagrok/grok_pipe`](https://hub.docker.com/r/datagrok/grok_pipe) | `1.19.0` |
| Grok Spawner | [`datagrok/grok_spawner`](https://hub.docker.com/r/datagrok/grok_spawner) | `2.16.0` |
| Grok Connect | [`datagrok/grok_connect`](https://hub.docker.com/r/datagrok/grok_connect) | `2.6.2` |
| Jupyter Kernel Gateway | [`datagrok/jupyter_kernel_gateway`](https://hub.docker.com/r/datagrok/jupyter_kernel_gateway) | `1.31.0` |
| Grok Registry Proxy | [`datagrok/grok_registry_proxy`](https://hub.docker.com/r/datagrok/grok_registry_proxy) | `1.27.1` |
| RabbitMQ | [`rabbitmq`](https://hub.docker.com/_/rabbitmq) | `4.0.5-management` |

The Helm chart (`oci://registry-1.docker.io/datagrok/datagrok`) is published in the same
repo as the `datagrok` image, with chart tags suffixed `-helm` to keep the namespaces
disjoint — so for Datagrok `1.27.3` use chart tag `1.27.3-helm`.

## Tag conventions

| Tag                            | Published from        | Use for             |
|--------------------------------|-----------------------|---------------------|
| `1.27.3`, `1.27.2`, ...        | Release tags          | Production          |
| `1.27.3-rc`, ...               | `release/*` branches  | Release candidates  |
| `bleeding-edge`                | `master`, nightly     | Dev / evaluation    |
| `latest`                       | Most recent stable    | Convenience alias   |

Helm chart tags follow the same scheme with a `-helm` suffix
(`1.27.3-helm`, `1.27.3-rc-helm`, `bleeding-edge-helm`).

## Independent release cadence

* **RabbitMQ** is an upstream image. Bump its tag only when there's a reason to —
  Datagrok versions don't pin to specific RabbitMQ versions.
* **PostgreSQL** is the [`pgvector/pgvector:pg17`](https://hub.docker.com/r/pgvector/pgvector)
  image when the Helm chart provisions Postgres in-cluster. For managed Postgres (RDS,
  Cloud SQL, Azure Database, on-prem), use any PostgreSQL 17 instance.

## Pinning per-deployment

| Deployment            | How to pin                                                                              |
|-----------------------|-----------------------------------------------------------------------------------------|
| AWS CFN (EKS / ECS)   | `DatagrokVersion`, `GrokPipeVersion`, ... CloudFormation parameters                     |
| Helm chart            | `--set datagrok.image.tag=1.27.3 --set grokPipe.image.tag=1.18.0 ...` or values overlay |
| Local Docker Compose  | `DATAGROK_VERSION=1.27.3` etc. before `docker compose up`                               |
| GCP Terraform module  | `datagrok_version = "1.27.3"` input variable                                            |

## Older releases

See the [release history](releases/release-history.md) for archived versions and breaking
changes. Older Datagrok releases pin different versions for sub-services — bump them as a
set when upgrading.
