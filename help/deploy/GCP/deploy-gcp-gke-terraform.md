---
title: "GCP GKE Terraform"
sidebar_position: 3
---

The deployment consists of a few docker containers, [database](../../develop/under-the-hood/infrastructure.md#1-core-components) for storing metadata,
and [persistent file storage](../../develop/under-the-hood/infrastructure.md#1-core-components) for storing files

This document contains instructions to deploy Datagrok using [Terraform](https://www.terraform.io/)
on [Google Kubernetes Engine](https://cloud.google.com/kubernetes-engine) with [SQL Server](https://cloud.google.com/sql-server)
and [Cloud Storage](https://cloud.google.com/storage).

We considered a lot of typical security nuances during the Terraform code development. As a result, you will
create a Datagrok infrastructure in GCP that applies to all standard security policies.

More information about Datagrok design and components:

* [Architecture](../../develop/under-the-hood/architecture.md)
* [Infrastructure](../../develop/under-the-hood/infrastructure.md)

## Basic usage

Use [Datagrok GCP Terraform module](https://github.com/datagrok-ai/tf-module-datagrok-core/tree/main/gcp) to deploy 

The use [HELM charts](https://github.com/datagrok-ai/tf-module-datagrok-core/tree/main/helm/datagrok) to deploy Kubernetes cluster

After the Datagrok container starts, the Datagrok server will deploy the database. You can check the status by
   checking the running task log in container output

## Advanced usage

The Terraform code is highly configurable. Feel free to adapt the code and variables to meet your needs and
requirements.
Terraform modules documentation:

* [Datagrok Core README](https://github.com/datagrok-ai/tf-module-datagrok-core/blob/main/gcp/README.md)
