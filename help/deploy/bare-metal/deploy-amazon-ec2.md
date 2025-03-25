---
title: "AWS EC2"
sidebar_position: 2
---

The deployment consists of a few docker containers, [database](../../develop/under-the-hood/infrastructure.md#1-core-components) for storing metadata,
and [persistent file storage](../../develop/under-the-hood/infrastructure.md#1-core-components) for storing files

Like a regular machine, any bare-metal server or virtual machine, including virtual machines in cloud providers, for
example, [AWS EC2](https://aws.amazon.com/ec2/), can be used.

As [database](../../develop/under-the-hood/infrastructure.md#1-core-components) Datagrok supports any PostgreSQL database out-of-the-box, including cloud
solutions for PostgreSQL database, for example [AWS RDS](https://aws.amazon.com/rds/).

For [persistent file storage](../../develop/under-the-hood/infrastructure.md#1-core-components), Datagrok supports a lot of options, including cloud solutions,
for example [AWS S3](https://aws.amazon.com/s3/) and Local File System storage.

To deploy Datagrok on an AWS EC2 instance with RDS as a database and S3 as a storage
use [deployment process for regular machine](deploy-regular.md) with
different [datagrok configuration](../configuration.md):

1. Replace `GROK_PARAMETERS` in [deployment process for regular machine](deploy-regular.md). `amazonStorageId`
   and `amazonStorageKey` can be omitted if you attach the IAM role to EC2 instances.

```json
{
  "amazonStorageRegion": "<S3_BUCKET_REGION>",
  "amazonStorageBucket": "<S3_BUCKET_NAME>",
  "amazonStorageId": "<S3_BUCKET_CREDS_ACCOUNT_ID>",
  "amazonStorageKey": "<S3_BUCKET_CREDS_SECRET_KEY>",
  "dbServer": "<RDS_ENDPOINT>",
  "db": "datagrok",
  "dbLogin": "datagrok",
  "dbPassword": "SoMeVeryCoMpLeXpAsSwOrD",
  "dbAdminLogin": "postgres",
  "dbAdminPassword": "postgres"
}
```
