---
title: "Deployment"
---

To deploy Datagrok services, you can use [Docker containers](https://www.docker.com/resources/what-container/#:~:text=A%20Docker%20container%20image%20is,tools%2C%20system%20libraries%20and%20settings.). Datagrok consists of [core](../develop/under-the-hood/infrastructure.md#1-core-components) containers, compute containers, a PostgreSQL database to store metadata, and [persistent file storage](../develop/under-the-hood/infrastructure.md#1-core-components) to store files.

Using Docker containers, you can deploy Datagrok on many environments, such as container services in the cloud providers, for example, [AWS ECS](#aws-deployment), [Kubernetes](#kubernetes-deployment), [bare-metal machines](#regular-machine-deployment), [virtual machines](#regular-machine-deployment), and so on.

To store data for Datagrok, we recommend using scalable and highly reliable solutions such as [AWS S3](https://aws.amazon.com/s3/) for persistent file storage and [AWS RDS](https://aws.amazon.com/rds/) as a PostgreSQL database.

## Local deployment

[Local deployment](docker-compose/docker-compose.md) is a quick way to see Datagrok in action using [Docker Compose](https://docs.docker.com/compose/). You can use it for local evaluation and development.

<!-- ### Deploy script

The interactive way to deploy the platform is to use
our [deployment script](https://github.com/datagrok-ai/public/blob/master/help/develop/admin/deploy/deploy.sh)

1. Download the script from
   repository: [deploy.sh](https://raw.githubusercontent.com/datagrok-ai/public/master/help/develop/admin/deploy/deploy.sh)
2. For AWS deployment, check that you have
   all [required permissions](https://github.com/datagrok-ai/public/blob/master/help/develop/admin/deploy/iam.list)
   on AWS account and installed [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html) with your credentials
3. Run the script. It will ask questions and deploy a Datagrok stand based on your answers. The supported deployment
   platform:
   ECS, Kubernetes, Virtual Machine.

   * EC2 instance should be treated like Virtual Machine. It is required to create EC2 instances before the script run.
     You can check how to create instances
     in [regular machine example preparations steps](bare-metal/deploy-regular.md#preparations)

```bash
sh deploy.sh
```
--->

## AWS deployment

We strongly recommend using [AWS ECS](https://aws.amazon.com/ecs/) for the Datagrok deployment. It provides a highly
scalable, fast container management service that makes it easy to manage application components.

We prepared three options for effortless and secure deployments to AWS:

* [Marketplace](aws/deploy-marketplace.md). The easiest way to start with Datagrok on AWS. [Marketplace](https://aws.amazon.com/marketplace) deployment scripts create a separate infrastructure for Datagrok from scratch.
* [CloudFormation (EKS)](aws/deploy-amazon-eks.md). Using the [CloudFormation template](https://aws.amazon.com/cloudformation/), you can customize the Datagrok infrastructure with an elaborate template that considers all standard security policies. The legacy [ECS variant](aws/deploy-amazon-ecs.md) is deprecated.
* [Terraform](aws/deploy-amazon-terraform.md). It is the most flexible solution. You can integrate Datagrok into your existing infrastructure with consideration of your security policies. However, [Terraform](https://www.terraform.io/) is also an advanced option that requires additional knowledge in infrastructure as a code area.

## Kubernetes deployment

Deploy Datagrok to any [Kubernetes](https://kubernetes.io/) cluster using the unified Helm chart. The
chart supports two deployment profiles:

* **Internal**: self-contained deployment where PostgreSQL and file storage run inside the cluster.
  Suitable for on-premises, MicroK8s, kind, or any Kubernetes cluster without managed cloud services.
* **Cloud**: PostgreSQL and file storage are offloaded to managed services (AWS RDS + S3, or GCP
  Cloud SQL + GCS). The chart includes value overlays for
  [EKS](https://aws.amazon.com/eks/) and [GKE](https://cloud.google.com/kubernetes-engine).

Any combination is possible (for example, internal database with S3 storage).

### Quick start (internal)

```bash
helm install datagrok ./deploy/k8s/datagrok/ \
  --set postgres.adminPassword=<ADMIN_PASS> \
  --set postgres.password=<DB_PASS> \
  --set ingress.host=datagrok.example.com \
  -n datagrok --create-namespace
```

This deploys all Datagrok services, an in-cluster PostgreSQL instance with persistent storage, and
local PVCs for file storage. RabbitMQ, Grok Connect, Grok Pipe, Jupyter Kernel Gateway, and
Grok Spawner (Kubernetes mode) are included.

### Cloud deployment (EKS)

```bash
helm install datagrok ./deploy/k8s/datagrok/ \
  -f ./deploy/k8s/datagrok/values-eks.yaml \
  --set postgres.external.host=<RDS_ENDPOINT> \
  --set postgres.password=<DB_PASS> \
  --set storage.s3.bucket=<BUCKET_NAME> \
  --set ingress.host=datagrok.example.com \
  -n datagrok --create-namespace
```

The EKS overlay (`values-eks.yaml`) configures external RDS, S3 storage, ALB ingress, IRSA service
accounts, and External Secrets Operator integration. A similar overlay exists for GKE
(`values-gke.yaml`) with Cloud SQL and GCS.

### Container registry

The chart supports three registry modes for plugin container images:

* **Internal** (`registry.type: internal`): a shared Docker Distribution registry for all instances
  in the cluster. Needs an ingress for `grok publish` from outside the cluster.
* **Proxy** (`registry.type: proxy`): deploys `grok_registry_proxy` that translates Datagrok JWT
  authentication to ECR or Google Artifact Registry credentials. Used for cloud deployments.
* **External** (`registry.type: external`): references an existing registry by URL without deploying
  anything.

### Grok Spawner

Grok Spawner runs in Kubernetes mode by default. The chart creates a dedicated ServiceAccount with
RBAC permissions to manage Deployments and Services in its namespace. Plugin containers are spawned
as Kubernetes pods on the same cluster.

### Configuration

All services are configured through a single `values.yaml`. The chart builds the `GROK_PARAMETERS`
JSON automatically from structured Helm values. See
[Server configuration](configuration.md) for the full list of configuration options.

For existing deployments, set `postgres.existingClaim`, `storage.local.existingDataClaim`, and
`storage.local.existingCfgClaim` to reuse existing PersistentVolumeClaims during migration.

## Regular machine deployment

You can deploy Datagrok to a [regular machine](bare-metal/deploy-regular.md): bare-metal servers or virtual machines, including [EC2 instances](https://aws.amazon.com/ec2/). However, this method is less reliable, scalable, and maintainable than others. You need to set up hosts manually and manage the data storage. Consider using other options if possible.

## Complete the setup

After the deployment, open the platform to complete the setup:  

1. [Configure authentication](complete-setup/configure-auth.md)
2. [Configure SMTP](complete-setup/configure-smtp.md)
3. [Install packages](complete-setup/install-packages.md)
