---
title: "AWS EC2"
sidebar_position: 2
---

Use [Docker Compose on a regular host](deploy-regular.md) on a single EC2 instance,
configured with [RDS](https://aws.amazon.com/rds/) for the metadata database and
[S3](https://aws.amazon.com/s3/) for file storage.

:::tip Prefer the CFN templates for production AWS

For new AWS stands the [CloudFormation (EKS)](../aws/deploy-amazon-eks.mdx) or
[CloudFormation (ECS)](../aws/deploy-amazon-ecs.mdx) templates provision the cluster, RDS,
S3, IAM, and load balancers end-to-end — and EKS is the AWS path under active
development. This page is for evaluation stands or one-off deployments where you manage
the EC2 instance yourself.

:::

## Configuration

Follow the [Regular machine](deploy-regular.md) install instructions, but replace the
`GROK_PARAMETERS` block on the `datagrok` service with one that points at your RDS
instance and S3 bucket. If the EC2 instance has an IAM role with S3 access, omit
`amazonStorageId` and `amazonStorageKey` — Datagrok resolves credentials from the
instance metadata service.

```json
{
  "dbServer": "<RDS_ENDPOINT>",
  "dbPort": "5432",
  "db": "datagrok",
  "dbLogin": "datagrok",
  "dbPassword": "<DB_PASSWORD>",
  "dbAdminLogin": "<POSTGRES_ADMIN_USER>",
  "dbAdminPassword": "<POSTGRES_ADMIN_PASSWORD>",
  "amazonStorageRegion": "<AWS_REGION>",
  "amazonStorageBucket": "<S3_BUCKET>",
  "amazonStorageId": "<AWS_ACCESS_KEY_ID>",
  "amazonStorageKey": "<AWS_SECRET_ACCESS_KEY>"
}
```

For RDS IAM authentication (passwordless DB login), set `dbUseIam: true` and omit
`dbPassword`. Datagrok mints SigV4 tokens via the IAM role attached to the EC2 instance.
See [Server configuration](../configuration.md) for the full list of options.

## See also

* [Regular machine](deploy-regular.md) — full Docker Compose install on any host
* [AWS CloudFormation (EKS)](../aws/deploy-amazon-eks.mdx) — recommended AWS path
* [AWS CloudFormation (ECS)](../aws/deploy-amazon-ecs.mdx) — Fargate-based AWS deployment
* [Components](../deploy.md#components) — service list and roles
