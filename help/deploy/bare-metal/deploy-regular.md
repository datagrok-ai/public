---
title: "Regular machine"
sidebar_position: 1
---

Datagrok runs as a set of Docker containers on top of a PostgreSQL metadata database and
persistent file storage. This page covers manual on-host installs — bare-metal servers,
on-prem VMs, single EC2 / GCE instances, or any other host you manage directly with
Docker Compose. The same Datagrok services run on every deployment — see
[Components](../deploy.md#components) for the canonical list.

:::tip Use container orchestration when you can

For new AWS stands prefer the [CloudFormation (EKS)](../aws/deploy-amazon-eks.mdx) or
[CloudFormation (ECS)](../aws/deploy-amazon-ecs.mdx) templates — they automate everything on
this page and are the AWS path under active development. For Kubernetes (on-prem, GKE, AKS,
existing clusters) use the [Helm chart](../k8s/install-helm-chart.md). This page is for
hosts without an orchestrator.

:::

## Prerequisites

* A Linux host (or a Linux VM on Windows / macOS) with at least 4 CPUs, 8 GB RAM, and
  60 GB free disk for the full stack including server-side scripting.
* [Docker Engine](https://docs.docker.com/engine/install/) and
  [Docker Compose v2](https://docs.docker.com/compose/install/) installed; the user that
  will run the stack added to the `docker` group.
* A PostgreSQL 17 database. The bundled compose stack runs Postgres in-cluster; for
  production prefer a managed instance ([AWS RDS](https://aws.amazon.com/rds/),
  [GCP Cloud SQL](https://cloud.google.com/sql), Azure Database for PostgreSQL, or an
  on-prem cluster).
* Object storage. The bundled compose stack uses a local volume; production stands
  typically use [S3](https://aws.amazon.com/s3/), [GCS](https://cloud.google.com/storage),
  Azure Blob, or an S3-compatible service.
* DNS or load-balancer pointing at the host on port `8080` (Datagrok) — direct port
  exposure works for evaluation, but production stands should sit behind a TLS-terminating
  reverse proxy.

## Install

1. Clone the public repository on the host (it ships the canonical compose file):

    ```bash
    git clone https://github.com/datagrok-ai/public.git
    cd public/docker
    ```

2. Open `localhost.docker-compose.yaml` and edit the `GROK_PARAMETERS` JSON on the
   `datagrok` service. Replace the values inline with your database and storage details
   (drop the `amazonStorage*` block if you're using local file storage):

    ```json
    {
      "dbServer": "<DATABASE_HOST>",
      "dbPort": "5432",
      "db": "datagrok",
      "dbLogin": "datagrok",
      "dbPassword": "<DB_PASSWORD>",
      "dbAdminLogin": "<POSTGRES_ADMIN_USER>",
      "dbAdminPassword": "<POSTGRES_ADMIN_PASSWORD>",
      "amazonStorageRegion": "us-east-2",
      "amazonStorageBucket": "<S3_BUCKET>"
    }
    ```

   See [Server configuration](../configuration.md) for every supported key — including
   GCS, Azure, RDS IAM auth, and TLS options.

3. Pull the images and start the full stack:

    ```bash
    docker compose -f localhost.docker-compose.yaml --project-name datagrok \
      --profile all up -d
    ```

   Use the `--profile` flags from
   [Local machine: advanced](../docker-compose/docker-compose-advanced.mdx#compose-profiles)
   to skip optional services (e.g., drop server-side scripting).

4. After about a minute the server is ready at
   `http://<HOST>:8080`. Sign in as `admin` / `admin` and
   [change the admin password](../complete-setup/configure-auth.md) on first login.

## Multi-host topologies

For multi-host installs (Datagrok services on one host, scripting / Jupyter Kernel
Gateway on another, or larger), use the [Helm chart](../k8s/install-helm-chart.md) on
Kubernetes. A single-node K8s distribution like [k3s](https://k3s.io/) or
[kind](https://kind.sigs.k8s.io/) is enough if you don't already run a cluster.

## On AWS EC2

For a single EC2 instance with [RDS](https://aws.amazon.com/rds/) and
[S3](https://aws.amazon.com/s3/) attached, follow this page and supply the RDS endpoint
and S3 bucket details in `GROK_PARAMETERS` — see
[AWS EC2 specifics](deploy-amazon-ec2.md). For multi-AZ, autoscaling, or load-balanced
production stands on AWS, use the [CFN ECS](../aws/deploy-amazon-ecs.mdx) or
[CFN EKS](../aws/deploy-amazon-eks.mdx) template instead — they provision the host
fleet, RDS, S3, and ALBs end-to-end.

## See also

* [Components](../deploy.md#components) — service list and roles
* [Server configuration](../configuration.md) — full `GROK_PARAMETERS` reference
* [Helm chart](../k8s/install-helm-chart.md) — single- or multi-node Kubernetes installs
