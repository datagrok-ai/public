---
title: "Docker Swarm"
sidebar_position: 3
---

Deploy Datagrok across two or more hosts with [Docker Swarm](https://docs.docker.com/engine/swarm/) —
a lighter-weight alternative to Kubernetes for multi-host installs. The same Datagrok
services run on every deployment — see [Components](../deploy.md#components) for the
canonical list.

:::tip Use Kubernetes when you can

For multi-node clusters under active load, prefer the
[Helm chart](../k8s/install-helm-chart.md) on Kubernetes — it ships rolling updates,
horizontal autoscaling, and proper service-mesh integration. This page is for swarms
where Kubernetes is overkill or unavailable.

:::

## Prerequisites

* Two or more Linux hosts with [Docker Engine](https://docs.docker.com/engine/install/) and
  the user that will run the stack added to the `docker` group on each.
* A PostgreSQL 17 instance reachable from the swarm — managed
  ([RDS](https://aws.amazon.com/rds/), [Cloud SQL](https://cloud.google.com/sql), Azure)
  or self-hosted.
* Object storage (S3, GCS, Azure Blob, or a local volume).
* DNS or a load balancer pointing at the manager node on port `8080`.

Per-node sizing:

* Datagrok manager node — 4 CPUs, 8 GB RAM, 30 GB disk.
* Worker node (scripting / Jupyter Kernel Gateway) — 4 CPUs, 8 GB RAM, 60 GB disk.

## Initialize the swarm

1. On the host that will run Datagrok core, initialize swarm mode:

    ```shell
    docker swarm init
    ```

   The output prints a `docker swarm join` command — copy the token.

2. On each worker host, run the join command from step 1 to add it to the swarm:

    ```shell
    docker swarm join --token <TOKEN> <MANAGER_HOST>:2377
    ```

3. Back on the manager, label nodes by role. The compose file places services on nodes
   based on these labels:

    ```shell
    docker node ls

    docker node update --label-add role=datagrok <DATAGROK_NODE_ID>
    docker node update --label-add role=cvm      <WORKER_NODE_ID>
    ```

## Deploy the stack

1. Download `docker/swarm.docker-compose.yaml` from the
   [public repository](https://github.com/datagrok-ai/public/blob/master/docker/swarm.docker-compose.yaml)
   to the manager node.

2. Edit the `GROK_PARAMETERS` JSON on the `datagrok` service. Point it at your database
   and storage; see [Server configuration](../configuration.md) for every supported key:

    ```json
    {
      "dbServer": "<DATABASE_HOST>",
      "dbPort": "5432",
      "db": "datagrok",
      "dbLogin": "datagrok",
      "dbPassword": "<DB_PASSWORD>",
      "dbAdminLogin": "<POSTGRES_ADMIN_USER>",
      "dbAdminPassword": "<POSTGRES_ADMIN_PASSWORD>",
      "amazonStorageRegion": "<AWS_REGION>",
      "amazonStorageBucket": "<S3_BUCKET>"
    }
    ```

3. Deploy the stack:

    ```shell
    docker stack deploy -c swarm.docker-compose.yaml datagrok
    ```

4. After about a minute the server is ready at `http://<MANAGER_HOST>:8080`. Sign in
   as `admin` / `admin` and
   [change the admin password](../complete-setup/configure-auth.md) on first login.

## User access

Datagrok and the optional Jupyter Kernel Gateway are reachable on `8080`. For
production stands put a TLS-terminating reverse proxy
([nginx](https://www.nginx.com/), [Traefik](https://traefik.io/),
[Caddy](https://caddyserver.com/)) in front of the manager node and route
`datagrok.example` to port `8080` on a public IP or load balancer.

## See also

* [Regular machine](deploy-regular.md) — single-host Docker Compose install
* [Helm chart](../k8s/install-helm-chart.md) — multi-node Kubernetes install
* [Components](../deploy.md#components) — service list and roles
* [Server configuration](../configuration.md) — `GROK_PARAMETERS` reference
