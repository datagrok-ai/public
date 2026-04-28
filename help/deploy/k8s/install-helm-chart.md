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
disjoint. `--version 1.27.3-helm` pulls a chart that deploys Datagrok `1.27.3`, and
every sub-service (grok-pipe, grok-spawner, grok-connect, Jupyter Kernel Gateway)
defaults to the same app tag. Override individual service image tags via
`--set <service>.image.tag=...` — see [Service versions](#service-versions) below.

Chart tags follow the same scheme as the Datagrok image tags with a `-helm` suffix; see
[Images and versions](../images.md#tag-conventions) for the full convention.

## Quick install (in-cluster Postgres + local storage)

```bash
helm install datagrok oci://registry-1.docker.io/datagrok/datagrok \
  --version 1.27.3-helm \
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

Every service image tag defaults to the chart version. The default set is documented on
[Images and versions](../images.md) and matches the
[AWS CloudFormation templates](../aws/deploy-amazon-eks.mdx). Override individual tags when
you need to run a newer `grok_connect` against an older `datagrok` core, or to pin a
specific service during a rollout:

```bash
helm install datagrok oci://registry-1.docker.io/datagrok/datagrok \
  --version 1.27.3-helm \
  --set datagrok.image.tag=1.27.3 \
  --set grokPipe.image.tag=1.18.0 \
  --set grokConnect.image.tag=2.6.2 \
  --set spawner.image.tag=2.15.0 \
  --set jkg.image.tag=1.31.0 \
  --set grokRegistryProxy.image.tag=1.18.0 \
  -n datagrok
```

RabbitMQ follows its own upstream release cadence and is not pinned to the Datagrok version.

## Production install on AWS EKS

For new AWS stands use the [AWS CloudFormation (EKS) template](../aws/deploy-amazon-eks.mdx) —
it provisions EKS, RDS, S3, IAM with IRSA, and installs this chart for you. The steps
below are for installing the chart directly into an EKS cluster you already manage.

1. Configure kubectl:
   ```bash
   aws eks update-kubeconfig --name <cluster-name>
   ```
2. Install the [AWS Load Balancer Controller](https://kubernetes-sigs.github.io/aws-load-balancer-controller/)
   in the cluster — the chart's `ingress.className: alb` annotations require it.
3. Provision RDS, S3, and the IRSA role for the Datagrok ServiceAccount yourself
   (`rds-db:connect`, `s3:Get/Put/List` on the bucket, optionally Secrets Manager read).
4. Save the EKS overlay below as `values-prod.yaml`, filling in the SET fields:

   ```yaml
   postgres:
     internal: false
     external:
       host: your-db.xxxxx.us-east-1.rds.amazonaws.com
       port: 5432
       ssl: true

   storage:
     type: s3
     s3:
       bucket: your-datagrok-bucket
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
       eks.amazonaws.com/role-arn: arn:aws:iam::ACCOUNT:role/datagrok-role

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
     --version 1.27.3-helm \
     -f values-prod.yaml \
     -n datagrok --create-namespace
   ```

The chart repo also ships a ready-made `values-eks.yaml` skeleton you can copy.

## Production install on GCP GKE

Same flow as EKS, but use the `values-gke.yaml` overlay (Cloud SQL + GCS + GKE
Workload Identity).

## Upgrades

```bash
helm upgrade datagrok oci://registry-1.docker.io/datagrok/datagrok \
  --version 1.27.4-helm \
  -f values-prod.yaml \
  -n datagrok
```

Always upgrade through consecutive minor versions for production stands; database schema
migrations run automatically on the first start of each new app version.

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
