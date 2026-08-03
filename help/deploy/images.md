---
title: "Images and versions"
sidebar_position: 12
description: Reference for Datagrok's Docker Hub images, version tags, and per-service pinning conventions across deployment paths.
keywords:
  - docker hub images
  - version tags
  - bleeding-edge tag
  - pin image version
  - helm chart tag
---

The Datagrok services that run on every deployment path
([Components](deploy.md#components)) are pulled from
[Docker Hub](https://hub.docker.com/u/datagrok). Defaults below are pinned to the latest
stable Datagrok release and are tested together. Override only when you need to pin a
specific release or test a candidate.

## Latest defaults

| Service | Image | Version |
|---|---|---|
| Datagrok | [`datagrok/datagrok`](https://hub.docker.com/r/datagrok/datagrok) | `1.27.7` |
| Grok Pipe | [`datagrok/grok_pipe`](https://hub.docker.com/r/datagrok/grok_pipe) | `1.22.2` |
| Grok Spawner | [`datagrok/grok_spawner`](https://hub.docker.com/r/datagrok/grok_spawner) | `2.23.0` |
| Grok Connect | [`datagrok/grok_connect`](https://hub.docker.com/r/datagrok/grok_connect) | `2.7.0` |
| Jupyter Kernel Gateway | [`datagrok/jupyter_kernel_gateway`](https://hub.docker.com/r/datagrok/jupyter_kernel_gateway) | `1.34.1` |
| Grok Registry Proxy | [`datagrok/grok_registry_proxy`](https://hub.docker.com/r/datagrok/grok_registry_proxy) | `1.30.2` |
| RabbitMQ | [`rabbitmq`](https://hub.docker.com/_/rabbitmq) | `4.0.5-management` |

These match what `datagrok/<service>:latest` currently resolves to (verified against
Docker Hub digests on 2026-08-03). Sub-service rows are pinned rather than tracking
`:latest` directly so that an unattended deployment using these exact versions stays
reproducible. See [tag conventions](#tag-conventions) below for what `:latest`
points at and the orphan-tag caution before substituting a different numeric version.

The Helm chart (`oci://registry-1.docker.io/datagrok/datagrok`) is published in the same
repo as the `datagrok` image, with chart tags suffixed `-helm` to keep the namespaces
disjoint — so for Datagrok `1.27.7` use chart tag `1.27.7-helm`.

## Configuration is versioned with the image {#release-binding}

Image tags are only half of a version. The configuration that starts a service is the
other half, and the two are released together — a Compose file, a Helm chart, or a
CloudFormation template written for one release can be invalid for another. Between
1.27.x and 1.28.0, for example, `connectorsSettings` in `GROK_PARAMETERS` changed from a
single object to an array of endpoints; a 1.27.x server given the array form exits at
startup with `Invalid argument(s): useGrokConnect`.

So always take the configuration from the same release as the images:

| Deployment            | Configuration comes from                                        |
|-----------------------|-----------------------------------------------------------------|
| Local Docker Compose  | `docker/*.docker-compose.yaml` on branch `release/<version>`    |
| Helm chart            | chart tag `<version>-helm` (the chart carries its own defaults) |
| AWS CFN (EKS / ECS)   | the template published for that release                         |

The Compose files carry the binding in a top-level `x-datagrok-release` field naming the
release they belong to, and every `datagrok/*` image in them defaults to a tag from that
release. On `master` the field reads `bleeding-edge`, so master's Compose file pulls
`bleeding-edge` images and never a stable release.

Every `datagrok/datagrok` image also records the branch it was built from in a `branch`
label, which is what makes a moving tag safe to use:

```shell
docker image inspect datagrok/datagrok:latest --format '{{.Config.Labels.branch}}'
# release/1.27.7
```

The install scripts read that label, fetch the Compose file from the branch it names, and
refuse to start a Compose file from one release against images from another. So
`DATAGROK_VERSION=latest` is bound just as tightly as an exact version — the pairing comes
from the image, not from the tag's name.

## Tag conventions

| Tag                          | Published from                       | Use for                       |
|------------------------------|--------------------------------------|-------------------------------|
| `1.27.7`, `1.27.6`, ...      | Release tags                         | Production (reproducible pin) |
| `1.27.8-rc`, ...             | `release/*` branches                 | Release candidates            |
| `latest`                     | Patch / minor / major bump builds    | Auto-track latest stable      |
| `bleeding-edge`              | `master` nightly cron                | Dev / evaluation              |

A version having a `release/X.Y.Z` branch does not mean it is published — the branch is
cut at release-candidate time, well before the image is promoted. Check
[release history](releases/release-history.md) or Docker Hub before pinning to one.

`:latest` is repointed each time a `Build-<service>` job runs with `BUILD_VERSION` set
to `patch`, `minor`, or `major` (not `bleeding-edge`). Use it when you want unattended
deployments to track the most recent stable build; pin to a specific semver from
[release history](releases/release-history.md) when you need reproducibility.

Helm chart tags follow the same scheme with a `-helm` suffix
(`1.27.3-helm`, `1.27.3-rc-helm`, `bleeding-edge-helm`). There is no `latest-helm`
chart — pin the chart by version.

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
| Helm chart            | `--version 1.27.7-helm` (pins the chart and its image defaults together)                |
| Local Docker Compose  | `DATAGROK_VERSION=1.27.7 bash datagrok-install-local.sh` — the script fetches the Compose file from the branch that image was built from. `latest` works the same way. |
| GCP Terraform module  | `datagrok_version = "1.27.7"` input variable                                            |

For Compose, setting `DATAGROK_VERSION` on a Compose file you already have is not enough:
it changes the image but not the configuration. Let the install script fetch the matching
file, or check out the `release/<version>` branch yourself — see
[configuration is versioned with the image](#release-binding).

## Older releases

See the [release history](releases/release-history.md) for archived versions and breaking
changes. Older Datagrok releases pin different versions for sub-services — bump them as a
set when upgrading.
