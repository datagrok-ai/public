<!-- TITLE: Deployment on AWS ECS Cluster-->
<!-- SUBTITLE: -->

# Deployment on AWS ECS Cluster

This document contains instructions to deploy Datagrok on [AWS ECS cluster](https://aws.amazon.com/ecs/).

## Prerequisites

1. Install [Docker Compose](https://docs.docker.com/compose/). If you do not have it, follow
   these [installation instructions](https://docs.docker.com/compose/install/) for your operating system.
2. To perform AWS CLI commands provided in the document
    1. [Install AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html)
    2. [Configure authorization for AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-quickstart.html)
3. Check that your default VPC has three subnets with internet gateway routing for internet facing Application Load
   Balancer
4. Create S3 bucket. For security reasons, we recommend to:
    * Disable public access to bucket and objects
    * Enable encryption
    ```shell
    #!/bin/sh

    DATAGROK_S3_BUCKET_NAME=example-unique-name-datagrok-s3
    aws s3api create-bucket \
       --bucket "$DATAGROK_S3_BUCKET_NAME" \
       --create-bucket-configuration LocationConstraint="$AWS_REGION" \
       --region "$AWS_REGION"
    aws s3api put-bucket-encryption \
       --bucket "$DATAGROK_S3_BUCKET_NAME" \
       --server-side-encryption-configuration '{"Rules": [{"ApplyServerSideEncryptionByDefault": {"SSEAlgorithm": "AES256"}}]}'
    aws s3api put-public-access-block \
       --bucket "$DATAGROK_S3_BUCKET_NAME" \
       --public-access-block-configuration "BlockPublicAcls=true,IgnorePublicAcls=true,BlockPublicPolicy=true,RestrictPublicBuckets=true"
    ```
5. Create RDS database for Datagrok
    * PostgreSQL 12
    * No publicly accessible
    ```shell
    #!/bin/sh

    aws rds create-db-instance \
      --db-instance-identifier "datagrok-rds" \
      --db-name "datagrok" \
      --engine 'postgres' \
      --engine-version '12.8' \
      --auto-minor-version-upgrade \
      --allocated-storage 50 \
      --max-allocated-storage 100 \
      --db-instance-class 'db.t3.large' \
      --master-username "postgres" \
      --master-user-password "postgres" \
      --port "5432" \
      --multi-az \
      --no-publicly-accessible \
      --storage-encrypted \
      --deletion-protection \
      --backup-retention-period 3 \
      --tags Key=project,Value=datagrok \
      --region "$AWS_REGION"
    ```
6. Create Docker ECS context: `docker context ecs create --from-env ECS`

## Setup Datagrok components

1. Download Docker Compose yaml
   file: [link](https://github.com/datagrok-ai/public/blob/master/docker/ecs.datagrok.docker-compose.yaml).
2. Replace `GROK_PARAMETERS` variables in `ecs.datagrok.docker-compose.yaml` with actual values or use environment
   variables automatic substitution
    ```shell
    #!/bin/sh

    export DATAGROK_S3_BUCKET_REGION=$AWS_REGION
    export DATAGROK_S3_BUCKET_NAME=example-unique-name-datagrok-s3
    export DATAGROK_RDS_ENDPOINT=$(aws rds describe-db-instances --db-instance-identifier "datagrok-rds" --output text --query 'DBInstances[].Endpoint.Address')
    ```
3. Use docker context for ECS: `docker context use ECS`
4. [Configure authorization for docker ECS context as for AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-quickstart.html)
5. Run docker compose deploy to ECS: `docker compose -p datagrok -f ecs.datagrok.docker-compose.yaml up -d`. It will
   create ECS cluster with tasks in [Fargate](https://aws.amazon.com/fargate/).

6. Allow newly created ECS cluster connect to RDS database for Datagrok: add Datagrok ECS Service Security Group to RDS
   database for Datagrok
    ```shell
    #!/bin/sh

    aws rds modify-db-instance \
      --db-instance-identifier "datagrok-rds" \
      --region "$AWS_REGION" \
      --vpc-security-group-ids "$(aws ecs describe-services \
        --cluster datagrok \
        --services "$(aws ecs list-services --cluster datagrok --output text --query 'serviceArns[0]')" \
        --output text --query 'services[].networkConfiguration.awsvpcConfiguration.securityGroups[]')" \
      --apply-immediately
    ```

7. Check deployment result: `docker compose -p datagrok -f ecs.datagrok.docker-compose.yaml ps`
8. Datagrok starts to deploy the database immediately after startup, you can check the status by checking running task
   log in [CloudWatch](https://aws.amazon.com/cloudwatch/)

8. Switch back to default docker context: `docker context use default`

9. Create CNAME DNS record <DATAGROK_DNS> to the newly created Application Load Balancer

10. Go in the web browser to <DATAGROK_DNS>, login to Datagrok using username "admin" and password "admin"

[//]: # (11. Edit settings in the Datagrok &#40;Tools | Settings...&#41;. Do not forget to click Apply to save new settings.)

[//]: # ()

[//]: # (* Connectors)

[//]: # (    * External Host: `grok_connect`)

## Setup Compute Virtual Machine

1. Download Docker Compose yaml
   file: [link](https://github.com/datagrok-ai/public/blob/master/docker/ecs.cvm.docker-compose.yaml)
   .
2. Use docker context for ECS: `docker context use ECS`
3. [Export environment variables for AWS](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-quickstart.html)

4. Run docker compose deploy to ECS: `docker compose -p cvm -f ecs.cvm.docker-compose.yaml up -d`. It will create ECS
   cluster with tasks in [Fargate](https://aws.amazon.com/fargate/).

5. Check deployment results: `docker compose -p cvm -f ecs.cvm.docker-compose.yaml ps`

6. Switch back to default docker context: `docker context use default`

7. Create CNAME DNS record <CVM_DNS> to the newly created Application Load Balancer

8. Go in the web browser to <CVM_DNS>, login to Datagrok using username "admin" and password "admin"

9. Edit settings in the Datagrok (Tools | Settings...). Do not forget to click Apply to save new settings.

* Scripting:
    * CVM Url: `http://<CVM_DNS>`
    * CVM Url Client: `http://<CVM_DNS>`
    * H2o Url: `http://<CVM_DNS>:54321`
    * Api Url: `http://<DATAGROK_DNS>/api`
    * Cvm Split: `true`
* Dev:
    * CVM Url: `http://<CVM_DNS>`
    * Cvm Split: `true`
    * Api Url: `http://<DATAGROK_DNS>/api`

See also:

* [Architecture](architecture.md)
* [Architecture Details](infrastructure.md)
