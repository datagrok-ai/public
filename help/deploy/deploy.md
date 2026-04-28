---
title: "Deployment"
---

Datagrok runs as a set of Docker containers on top of a PostgreSQL metadata database and
persistent file storage. The same containers ship across every deployment path; what changes is
where they run and what manages their lifecycle.

## Components

Every Datagrok stand runs the following services. The same images are used on every deployment
path; the [Helm chart](k8s/install-helm-chart.md) is the canonical reference for service
configuration, and [Images and versions](images.md) tracks the latest pinned tags.

| Service                     | Image                                | Role |
|-----------------------------|--------------------------------------|------|
| **Datagrok**                | `datagrok/datagrok`                  | Core REST API, web client, authn/authz, metadata persistence, Nginx reverse proxy. |
| **PostgreSQL**              | `pgvector/pgvector:pg17`             | Metadata store (users, projects, packages, queries, scripts, file index). RDS / Cloud SQL / in-cluster. |
| **grok\_pipe**              | `datagrok/grok_pipe`                 | WebSocket multiplexer for streaming DataFrames and script results between clients and Jupyter workers. |
| **grok\_spawner**           | `datagrok/grok_spawner`              | Manages plugin container lifecycle on Docker / ECS / Kubernetes (selected per deployment). |
| **grok\_connect**           | `datagrok/grok_connect`              | JDBC bridge for 30+ external databases. |
| **JupyterKernelGateway**    | `datagrok/jupyter_kernel_gateway`    | Server-side script execution (Python, R, Julia, JavaScript, Octave). |
| **RabbitMQ**                | `rabbitmq`                           | AMQP broker for the call queue (script and function execution). Independent release cadence. |
| **grok\_registry\_proxy**   | `datagrok/grok_registry_proxy`       | Optional. Proxies plugin image pulls from a backing registry (ECR, Docker Hub) using Datagrok JWT auth, so users never see registry credentials. |

For object storage, use AWS S3, Google Cloud Storage, Azure Blob, or a local volume — see the
chart's `storage.type` value or [File storage](../develop/under-the-hood/architecture.md).

## Deployment paths

Five paths are supported. They share images and configuration parameters; pick by where the
Datagrok stand will live.

| Path | Use when |
|------|----------|
| [Local Docker Compose](docker-compose/docker-compose.mdx) | Single machine — laptop or VM — for evaluation, demos, or development. Self-contained PostgreSQL inside Compose. |
| [Advanced Docker Compose](docker-compose/docker-compose-advanced.mdx) | Single-machine deployments that need separate data volumes, the JS-API debug stack, or other custom topology. |
| [Kubernetes Helm chart](k8s/install-helm-chart.md) | Any Kubernetes cluster: on-prem, GKE, AKS, kind, k3s, MicroK8s, or a pre-existing EKS. The EKS CFN template uses the same chart. |
| [AWS CloudFormation (EKS)](aws/deploy-amazon-eks.mdx) | **Recommended for new AWS stands.** Provisions EKS, RDS, S3, IAM with IRSA, and installs the Helm chart automatically. |
| [AWS CloudFormation (ECS)](aws/deploy-amazon-ecs.mdx) | Existing ECS stacks. Same RDS / S3 logical IDs as the EKS template, so an in-place stack-template swap migrates without re-creating data. Targeted for deprecation. |

The [AWS Marketplace](aws/deploy-marketplace.md) listing wraps the EKS template for one-click,
infrastructure-isolated installs.

[Terraform on AWS](aws/deploy-amazon-terraform.md) and [Terraform on GCP](GCP/deploy-gcp-gke-terraform.md)
are available for teams that integrate Datagrok into existing infrastructure-as-code pipelines.

[Bare-metal / VM](bare-metal/deploy-regular.md) is the manual Docker-on-host path for environments
without container orchestration.

## EKS or Helm directly?

The EKS CloudFormation template **calls** the Helm chart — it is not a separate deployment. Pick
based on what infrastructure you already manage:

* Use **CloudFormation (EKS)** if you want one stack to provision the cluster, RDS, S3, IAM, and the
  application together. The template can also target a pre-existing EKS cluster
  (`UseExistingCluster=true`).
* Use the **Helm chart directly** if your cluster, database, and object storage already exist, or
  if the cluster is not on AWS (GKE, AKS, on-prem). The chart ships ready-made overlays for EKS
  (`values-eks.yaml`) and GKE (`values-gke.yaml`).

## Complete the setup

After the platform is reachable, configure cross-cutting concerns:

1. [Authentication](complete-setup/configure-auth.md) (LDAP, OAuth, SAML, IAP, etc.)
2. [SMTP](complete-setup/configure-smtp.md)
3. [Install packages](complete-setup/install-packages.md)
4. [S3 backups](complete-setup/configure-s3-backup.md) (cloud deployments)
