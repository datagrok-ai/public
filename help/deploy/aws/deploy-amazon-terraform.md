---
title: "AWS Terraform"
sidebar_position: 2
---

The deployment consists of several Docker containers, a [database](../../develop/under-the-hood/infrastructure.md#1-core-components) for storing metadata,
and [persistent file storage](../../develop/under-the-hood/infrastructure.md#1-core-components) for storing files.

This document contains instructions to deploy Datagrok using [Terraform](https://www.terraform.io/)
on an [AWS ECS cluster](https://aws.amazon.com/ecs/) with [AWS RDS](https://aws.amazon.com/rds/)
and [AWS S3](https://aws.amazon.com/s3/).

We addressed many typical security considerations during the development of the Terraform code. As a result, the deployed Datagrok infrastructure adheres to standard AWS security best practices.

More information about Datagrok design and components:

* [Architecture](../../develop/under-the-hood/architecture.md)
* [Infrastructure](../../develop/under-the-hood/infrastructure.md)

## Basic usage

### Prerequisites

1. Check that you
   have the [required permissions](https://github.com/datagrok-ai/public/blob/master/help/develop/admin/deploy/iam.list)
   on the AWS account to perform Terraform deployment to ECS.
2. Check that you
   have the required [S3 permissions](https://www.terraform.io/language/settings/backends/s3#s3-bucket-permissions)
   and [DynamoDB permissions](https://www.terraform.io/language/settings/backends/s3#dynamodb-table-permissions).

3. Create an IAM role for CloudFormation stack deployment with the following policies attached:

   * IAMFullAccess
   * PowerUserAccess
  
  You can name it something like `CloudFormationExecutionRole`. Its ARN will be used as the `iam_role_arn` variable in the Terraform module.

## Deploy Datagrok components

1. You can use our example as a template and modify it to your needs.
   Example: [link](https://github.com/datagrok-ai/tf-module-datagrok-core/blob/main/aws/examples/vpc/main.tf).

2. Make sure to change the module source to `source = "git@github.com:datagrok-ai/tf-module-datagrok-core.git//aws?ref=main"` or any other version you want to use.

3. [Set up AWS credentials to run Terraform code](https://developer.hashicorp.com/terraform/language/settings/backends/s3#credentials-and-shared-configuration)

4. Apply the Terraform code
   ```shell
   terraform init --upgrade
   terraform apply
   ```

5. Once the Datagrok container starts, the server will initialize the database. You can monitor the status by checking the running task logs in [CloudWatch](https://aws.amazon.com/cloudwatch/).

## Advanced usage

The Terraform code is highly configurable. Feel free to adapt the code and variables to meet your needs and requirements.
Terraform modules documentation:

* [Datagrok Core README](https://github.com/datagrok-ai/tf-module-datagrok-core/blob/main/aws/README.md)
