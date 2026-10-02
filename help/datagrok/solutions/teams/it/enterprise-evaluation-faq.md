---
title: "Enterprise evaluation FAQ"
description: FAQ covering architecture, security, scalability, and extensibility for enterprise Datagrok evaluations.
keywords:
  - enterprise evaluation
  - security assessment
  - encryption at rest
  - disaster recovery
  - infrastructure as code
  - high availability
  - vendor evaluation checklist
---

* [Architecture](../../../../develop/under-the-hood/architecture.md)
  * [Infrastructure](../../../../develop/under-the-hood/infrastructure.md)
  * [Deployment](../../../../deploy/deploy.md)

* Security
  * [Security, authentication, and authorization](../../../../govern/access-control/access-control.md)
  * [Encryption at rest](#encryption-at-rest)
  and [encryption in transit](#encryption-in-transit)
  * [Vulnerability remediation](security.md#vulnerability-remediation)

* Enterprise readiness
  * [Managing an enterprise instance](../../../../govern/manage-enterprise.md)
  * [Logging and monitoring](#logging-and-monitoring)
  * [Backup and restore](#backup-and-restore)
  * [Disaster recovery (HA/DR)](#disaster-recovery)
  * [Infrastructure as code](#infrastructure-as-a-code) (ability to deploy using standard DevOps tools)

* Interoperability
  * Calling web services from client and server with proper auth: [OpenAPI](../../../../access/open-api.md)
  * Interacting with Datagrok
    * [server API](#server-api)
    * [client API](../../../../develop/packages/js-api.md)
  * Connecting to common data sources
    * [relational databases](https://youtu.be/YJmSvh3_uCM)
    * [local files](https://datagrok.ai/img/slides/access-file-formats.mp4)
    * [file shares and cloud storage](../../../../access/files/files.md)
  <!--Incorrect GIF* [Embedding a Datagrok visualization into a custom web application](https://datagrok.ai/embed_test.html)-->
  * [Embedding a custom visualization into Datagrok](../../../../visualize/viewers/markup.md)

* Developer experience
  * [Debug environment when developing customizations](https://youtu.be/PDcXLMsu6UM)
  * [Devops process including deployment of packages, their dependencies and versioning](../../../../develop/develop.md)
  * [Concurrent work by team of developers](../../../../develop/develop.md#development)

* Scalability and performance
  * [Maximum dataset sizes and in-memory performance](../../../../develop/under-the-hood/performance.md)
  * [Stability under concurrent user load](stress-testing-results.md)
  * [Scaling and stability under load](../../../../develop/under-the-hood/infrastructure.md#scalability)

* Extensibility
  * [Creating custom visualizations](https://github.com/datagrok-ai/public/tree/master/packages/BiostructureViewer)
  * [Creating custom server-side components](https://github.com/datagrok-ai/public/tree/master/packages/Admetica)
  * [Creating custom scripts](../../../../compute/scripting/scripting.mdx) and utilizing them in other components
  * [Ability to reskin Datagrok to appear as a fit-for-purpose web application](https://public.datagrok.ai/apps/HitTriage/HitTriage?browse=apps)
  * [Ability to build custom application including data entry, workflow, data model, state management, persistence, etc](https://github.com/datagrok-ai/public/tree/master/packages)

* Frontend
  * [Holding and exploring large datasets in the browser](../../../../develop/under-the-hood/performance.md#in-memory-database)
  * [Visualizing datasets with interactive viewers](../../../../visualize/viewers/viewers.md)
  * [2D and 3D structure rendering, sketching, and search](../../../../datagrok/solutions/domains/chem/chem.md)
  * [3D biostructures](../../../../visualize/viewers/biostructure.md)
  * Interactivity
    * [Live data masking](https://youtu.be/67LzPsdNrEc)
    * [Filter by selection](https://youtu.be/67LzPsdNrEc)
  * [Developing custom viewers, including non-native components such as React containers](../../../../develop/how-to/viewers/develop-custom-viewer.md)

## Encryption at rest

For AWS deployment, we rely on Amazon's built-in encryption for
[RDS](https://docs.aws.amazon.com/AmazonRDS/latest/UserGuide/Overview.Encryption.html)
and
[S3 buckets](https://docs.aws.amazon.com/AmazonS3/latest/userguide/bucket-encryption.html).
Credentials for data connections are additionally encrypted by Datagrok with
[server keys](../../../../govern/access-control/server-keys.md).

## Encryption in transit

All client-server communications use the [HTTPS](https://en.wikipedia.org/wiki/HTTPS) protocol, which means it is secure
and encrypted.

## Server API

The Datagrok client uses an HTTP REST API to interact with the server. You must pass an authentication token to access all
features.
[Proof of concept video](https://www.youtube.com/watch?v=TjApCwd_3hw)

## Logging and monitoring

Datagrok works with the monitoring tools you already use, in any cloud or on-premises:

* **Logs and alerts.** **Log sync** pushes logs, audit events, alerts, and a heartbeat to
  Amazon CloudWatch, Google Cloud Logging, or any OpenTelemetry (OTLP) collector. See
  [Export logs](../../../../govern/audit/audit.md#export-logs).
* **Health checks.** The `/api/admin/health` endpoint reports the status of every service and needs
  no sign-in, so load balancers and uptime monitors can probe it.
* **Usage and performance.** [Usage Analysis](../../../../govern/audit/usage-analysis.md) shows
  user activity, errors, and server and database metrics.

For details, see [Monitor the platform](../../../../govern/manage-enterprise.md#monitor-the-platform).

## Backup and restore

Back up these together so that you can restore them to the same point in time:

* The Postgres database that holds metadata, audit logs, and encrypted credentials. You can back it up and restore it
  as a standard Postgres database, for example with scheduled RDS snapshots on AWS.
* The file storage (S3, Google Cloud Storage, or a local volume). On AWS, see
  [S3 backup](../../../../deploy/complete-setup/configure-s3-backup.md).
* [Server key](../../../../govern/access-control/server-keys.md) material, if you keep keys in an external backend.

## Disaster recovery

[Datagrok supports Docker installation, Amazon cluster will immediately restart failed instance](https://www.youtube.com/watch?v=oFs9RShkHT8)

## Infrastructure as a Code

Datagrok Docker containers are built using Jenkins. All software is upgraded and patched on every build.

You can deploy Datagrok with standard DevOps tools:

* [Docker Compose](../../../../deploy/docker-compose/docker-compose.mdx)
* [Helm chart for Kubernetes](../../../../deploy/k8s/install-helm-chart.md)
* [CloudFormation template for Amazon EKS](../../../../deploy/aws/deploy-amazon-eks.mdx)
* [Terraform for AWS](../../../../deploy/aws/deploy-amazon-terraform.md)
* [Terraform for Google Cloud (GKE)](../../../../deploy/GCP/deploy-gcp-gke-terraform.md)
