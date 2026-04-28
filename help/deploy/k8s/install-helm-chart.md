---
title: "Install Datagrok with Helm"
sidebar_position: 5
---

The Datagrok Helm chart is published as an OCI artifact on Docker Hub. It supports
in-cluster PostgreSQL with local PVCs (laptop / single-node clusters), as well as
managed cloud databases (RDS / Cloud SQL) with object storage (S3 / GCS) on EKS or
GKE.

:::tip On AWS?

For turnkey EKS deployments, use the [AWS CloudFormation (EKS)](../aws/deploy-amazon-eks.mdx)
template — it provisions the cluster, RDS, S3, and IAM, and installs this chart
automatically. This page covers manual Helm installs for any Kubernetes cluster
(on-prem, GKE, AKS, kind, k3s, MicroK8s).

:::

## Prerequisites

* Kubernetes 1.27+ with a default StorageClass (or an existing PVC for stateful data)
* Helm 3.8+ (for OCI registry support)
* For cloud installs: a managed Postgres instance, an object storage bucket, and an
  IAM role / service account that the cluster can use to reach them

## Chart version

The chart is published as an OCI artifact in the shared `datagrok/datagrok` Docker Hub
repo under tags suffixed with `-helm` to keep the chart and image tag namespaces
disjoint. `--version 1.27.0-helm` pulls a chart that deploys Datagrok `1.27.0`, and
every sub-service (grok-pipe, grok-spawner, grok-connect, Jupyter Kernel Gateway)
defaults to the same app tag. Override individual service image tags via
`--set <service>.image.tag=...` — see [Service versions](#service-versions) below.

Available chart versions:

| Version                            | Published from               | Use for           |
|------------------------------------|------------------------------|-------------------|
| `1.27.0-helm`, `1.26.5-helm`, ...  | release tags                 | Production        |
| `1.27.0-rc-helm`, ...              | `release/*` branches         | Release candidates|
| `bleeding-edge-helm`               | `master`, nightly            | Dev / evaluation  |

## Quick install (in-cluster Postgres + local storage)

```bash
helm install datagrok oci://registry-1.docker.io/datagrok/datagrok \
  --version 1.26.5-helm \
  --namespace datagrok --create-namespace \
  --set postgres.password=$(openssl rand -base64 24) \
  --set postgres.adminPassword=$(openssl rand -base64 24) \
  --set ingress.host=datagrok.example.com
```

This installs PostgreSQL, RabbitMQ, the datagrok app, grok-pipe, grok-spawner,
grok-connect, and JupyterKernelGateway with default resources. Suitable for
evaluation, dev, or single-tenant production.

To track the latest unstable build (rebuilt after every merge to `master`), use
`--version bleeding-edge-helm` instead of a release version.

## Service versions

Every service image tag defaults to the chart version. Override individual tags
when you need to run a newer grok-connect against an older datagrok core, or
to pin a specific service during a rollout:

```bash
helm install datagrok oci://registry-1.docker.io/datagrok/datagrok \
  --version 1.27.0-helm \
  --set datagrok.image.tag=1.27.0 \
  --set grokPipe.image.tag=1.27.0 \
  --set grokConnect.image.tag=1.27.1 \
  --set spawner.image.tag=1.27.0 \
  --set jkg.image.tag=1.27.0 \
  --set rabbitmq.image.tag=4.0.5-management \
  -n datagrok
```

RabbitMQ follows its own upstream release cadence and is not tied to the Datagrok
version.

## Production install on AWS EKS

The easiest path is the [AWS CloudFormation (EKS) template](../aws/deploy-amazon-eks.mdx),
which provisions EKS, RDS, S3, IAM with IRSA, and installs this Helm chart for
you. The manual steps below are useful when you want to manage the cluster
yourself or when deploying into a pre-existing EKS cluster.

The `datagrok-eks-cfn.yaml` CloudFormation template (in the chart repo under
`deploy/k8s/`) is the same source-of-truth template used by the shipped CFN
variants; you can deploy it directly if you prefer a single-file template.

1. Deploy the CFN stack and capture its outputs (RDS endpoint, bucket name, IAM
   role ARN, ACM cert ARN).
2. Configure kubectl:
   ```bash
   aws eks update-kubeconfig --name <cluster-name>
   ```
3. Install the AWS Load Balancer Controller (the CFN stack's `PostDeployCommands`
   output prints the exact `helm install` commands).
4. Save the EKS overlay below as `values-prod.yaml`, filling in the SET fields:

   ```yaml
   postgres:
     internal: false
     external:
       host: your-db.xxxxx.us-east-1.rds.amazonaws.com   # from CFN
       port: 5432
       ssl: true

   storage:
     type: s3
     s3:
       bucket: your-datagrok-bucket                       # from CFN
       region: us-east-1

   ingress:
     enabled: true
     className: alb
     host: datagrok.example.com
     annotations:
       alb.ingress.kubernetes.io/scheme: internet-facing
       alb.ingress.kubernetes.io/target-type: ip
       alb.ingress.kubernetes.io/listen-ports: '[{"HTTPS":443}]'
       alb.ingress.kubernetes.io/certificate-arn: arn:aws:acm:...:certificate/...
     tls:
       enabled: false      # ALB terminates TLS

   serviceAccount:
     create: true
     annotations:
       eks.amazonaws.com/role-arn: arn:aws:iam::ACCOUNT:role/datagrok-role  # from CFN

   registry:
     type: proxy
     proxy:
       backendUrl: https://ACCOUNT.dkr.ecr.REGION.amazonaws.com

   credentials:
     source: externalSecrets
     externalSecrets:
       enabled: true
       remoteKey: datagrok/prod                             # AWS Secrets Manager key
   ```

5. Install:
   ```bash
   helm install datagrok oci://registry-1.docker.io/datagrok/datagrok \
     --version 1.26.5-helm \
     -f values-prod.yaml \
     -n datagrok --create-namespace
   ```

The chart repo also ships a ready-made `values-eks.yaml` skeleton you can copy.

### CFN profiles

The `datagrok-eks-cfn.yaml` template supports two deployment profiles:

| Profile | `UseExistingCluster` | Creates | Use when |
|---------|---------------------|---------|----------|
| **Full stack** (default) | `false` | EKS cluster, node group, RDS, S3, IAM, secrets | Starting from scratch |
| **Existing cluster** | `true` | RDS, S3, IAM/IRSA, secrets only | You already have an EKS cluster |

#### Using an existing EKS cluster

If your organization already manages an EKS cluster, you can deploy only the
Datagrok-specific infrastructure (database, storage, IAM roles) into it:

1. Gather your cluster details:
   ```bash
   CLUSTER=my-cluster
   aws eks describe-cluster --name $CLUSTER --query '{
     OIDCIssuerUrl: cluster.identity.oidc.issuer,
     SecurityGroupId: cluster.resourcesVpcConfig.clusterSecurityGroupId
   }'
   ```

2. Deploy the CFN stack with `UseExistingCluster=true`:
   ```bash
   aws cloudformation create-stack --stack-name datagrok \
     --template-body file://datagrok-eks-cfn.yaml \
     --capabilities CAPABILITY_NAMED_IAM \
     --parameters \
       ParameterKey=UseExistingCluster,ParameterValue=true \
       ParameterKey=ExistingClusterName,ParameterValue=$CLUSTER \
       ParameterKey=ExistingClusterOIDCIssuerUrl,ParameterValue=<oidc-url> \
       ParameterKey=ExistingClusterSecurityGroupId,ParameterValue=<sg-id> \
       ParameterKey=VpcId,ParameterValue=<vpc-id> \
       ParameterKey=PrivateSubnets,ParameterValue='<subnet-1>,<subnet-2>' \
       ParameterKey=PublicSubnets,ParameterValue='<subnet-1>,<subnet-2>' \
       ParameterKey=DBSubnets,ParameterValue='<subnet-1>,<subnet-2>'
   ```

3. Install the Helm chart using the stack outputs (same as the full-stack flow,
   step 4 onward).

> The stack does **not** create an EKS cluster, node group, or Fargate profile.
> Ensure your existing cluster has nodes (or Fargate profiles) capable of
> scheduling the Datagrok pods and that the AWS Load Balancer Controller is
> already installed.

## Production install on GCP GKE

Same flow as EKS, but use the `values-gke.yaml` overlay (Cloud SQL + GCS + GKE
Workload Identity).

## Upgrades

```bash
helm upgrade datagrok oci://registry-1.docker.io/datagrok/datagrok \
  --version 1.26.6-helm \
  -f values-prod.yaml \
  -n datagrok
```

The chart's PostgreSQL StatefulSet, datagrok-data, and datagrok-cfg PVCs are
preserved across upgrades. Database schema migrations run automatically on the
first start of a new app version.

## Reusing existing PVCs

If you're migrating an existing Datagrok install to the chart and want to keep
your data, point the chart at the existing claims:

```yaml
postgres:
  internal: true
  existingClaim: my-existing-postgres-pvc

storage:
  type: local
  local:
    existingDataClaim: my-existing-data-pvc
    existingCfgClaim: my-existing-cfg-pvc
```

When `postgres.existingClaim` is set, the chart skips its own
`volumeClaimTemplates` and mounts the data volume at the PV root (no `pgdata`
subPath), matching the on-disk layout used by Datagrok installs prior to the
chart.

## Uninstall

```bash
helm uninstall datagrok -n datagrok
```

PVCs are NOT deleted by `helm uninstall`. To remove all data permanently:

```bash
kubectl delete pvc -n datagrok -l app.kubernetes.io/instance=datagrok
kubectl delete namespace datagrok
```
