---
title: "GCP GKE Terraform"
sidebar_position: 1
---

Deploy Datagrok on Google Cloud Platform using a Terraform module that provisions GKE,
Cloud SQL (PostgreSQL), GCS, and a VPC, then installs the [Helm chart](../k8s/install-helm-chart.md)
into the cluster. This is the GCP analogue to the [AWS CloudFormation (EKS)](../aws/deploy-amazon-eks.mdx)
template.

The same Datagrok services run on every deployment — see
[Components](../deploy.md#components) for the canonical list.

:::tip Already have a GKE cluster?

If your team already manages GKE, install the [Helm chart](../k8s/install-helm-chart.md)
directly with the `values-gke.yaml` overlay (Cloud SQL + GCS + Workload Identity). The
Terraform module is for stands that need GKE, Cloud SQL, GCS, IAM, and the chart
provisioned together.

:::

## Prerequisites

* GCP project with billing enabled and the
  [GKE](https://cloud.google.com/kubernetes-engine/docs/quickstarts/create-cluster),
  [Cloud SQL](https://cloud.google.com/sql/docs/postgres/quickstart),
  [Cloud Storage](https://cloud.google.com/storage/docs/quickstart),
  and Compute Engine APIs enabled.
* Terraform 1.5+ and a GCP service account credential with `Project Owner` (or the
  equivalent narrower permissions: GKE admin, Cloud SQL admin, Storage admin, IAM admin,
  Network admin).
* `gcloud` and `kubectl` installed locally.

## Deploy

The [Datagrok GCP Terraform module](https://github.com/datagrok-ai/tf-module-datagrok-core/tree/main/gcp)
creates the GKE cluster, Cloud SQL instance, GCS bucket, and supporting IAM, then
installs the [Helm chart](https://github.com/datagrok-ai/tf-module-datagrok-core/tree/main/helm/datagrok)
into the cluster.

```hcl
module "datagrok_gcp" {
  source  = "github.com/datagrok-ai/tf-module-datagrok-core//gcp"
  project = "<gcp-project-id>"
  region  = "europe-west3"

  cluster_name      = "datagrok"
  datagrok_version  = "1.27.3"
  datagrok_dns      = "datagrok.example.com"
  acm_cert          = ""        # use Google-managed cert via the chart's ingress instead
}
```

Apply:

```bash
terraform init
terraform apply
```

After Terraform finishes, the chart is already installed and Datagrok is reachable
through the ingress. Configure kubectl to inspect the running stand:

```bash
gcloud container clusters get-credentials datagrok --region europe-west3 \
  --project <gcp-project-id>
kubectl -n datagrok get pods
```

## Service versions

The Terraform module's `datagrok_version` input pins all service image tags to the same
release. Defaults track the latest stable Datagrok release; see the
[release history](../releases/release-history.md) and the equivalent
[Helm chart overrides](../k8s/install-helm-chart.md#service-versions) for granular
per-service control.

## Update Datagrok components

Bump `datagrok_version` (and the related variables, if you pinned services
individually) and re-apply:

```bash
terraform apply -var datagrok_version=1.27.4
```

Cloud SQL and the GCS bucket are not replaced on `apply`. Database schema migrations run
automatically on the first start of the new app version. For a manual upgrade with the
chart only (no Terraform changes), see
[Helm chart upgrades](../k8s/install-helm-chart.md#upgrades).

## See also

* [Install Datagrok with Helm](../k8s/install-helm-chart.md) — direct chart install for
  any cluster (including pre-existing GKE)
* [AWS CloudFormation (EKS)](../aws/deploy-amazon-eks.mdx) — AWS analogue
* [GCP Terraform module README](https://github.com/datagrok-ai/tf-module-datagrok-core/blob/main/gcp/README.md)
