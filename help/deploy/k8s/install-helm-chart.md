---
title: "Install Datagrok with Helm"
sidebar_position: 5
---

The Datagrok Helm chart is published as an OCI artifact on Docker Hub. It supports
in-cluster PostgreSQL with local PVCs (laptop / single-node clusters), as well as
managed cloud databases (RDS / Cloud SQL) with object storage (S3 / GCS) on EKS or
GKE.

## Prerequisites

* Kubernetes 1.27+ with a default StorageClass (or an existing PVC for stateful data)
* Helm 3.8+ (for OCI registry support)
* For cloud installs: a managed Postgres instance, an object storage bucket, and an
  IAM role / service account that the cluster can use to reach them

## Quick install (in-cluster Postgres + local storage)

```bash
helm install datagrok oci://registry-1.docker.io/datagrok/datagrok \
  --version 1.26.5 \
  --namespace datagrok --create-namespace \
  --set postgres.password=$(openssl rand -base64 24) \
  --set postgres.adminPassword=$(openssl rand -base64 24) \
  --set ingress.host=datagrok.example.com
```

This installs PostgreSQL, RabbitMQ, the datagrok app, grok-pipe, grok-spawner,
grok-connect, and JupyterKernelGateway with default resources. Suitable for
evaluation, dev, or single-tenant production.

To track the latest unstable build (rebuilt nightly from `master`), use
`--version 0.0.0-bleeding-edge` instead of a release version. (Helm requires
SemVer 2.0, so the bleeding-edge chart is published as a prerelease version.)
Stable releases are
published on a versioned cadence; bleeding-edge is rebuilt by Jenkins after every
merge to master.

## Production install on AWS EKS

The `datagrok-eks-cfn.yaml` CloudFormation template (in the chart repo under
`deploy/k8s/`) provisions EKS, RDS, S3, and the necessary IAM with IRSA.

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
     --version 1.26.5 \
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
  --version 1.26.6 \
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
