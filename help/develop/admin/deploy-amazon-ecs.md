<!-- TITLE: Deployment on AWS ECS using Docker Compose -->
<!-- SUBTITLE: -->

# Deployment on AWS ECS using Docker Compose

Datagrok consist of Docker containers, [database](infrastructure.md#database)
and [persistent file storage](infrastructure.md#storage).

Docker containers allow installing Datagrok on any platform, including container services in cloud providers, for
example [AWS ECS](https://aws.amazon.com/ecs/).

As [database](infrastructure.md#database) Datagrok supports any PostgreSQL database out-of-the-box, including cloud
solutions for PostgreSQL database, for example [AWS RDS](https://aws.amazon.com/rds/).

For [persistent file storage](infrastructure.md#storage) Datagrok supports a lot of options, including cloud solutions,
for example [AWS S3](https://aws.amazon.com/s3/).

This document contains instructions to deploy Datagrok using [Docker Compose](https://docs.docker.com/compose/)
on [AWS ECS cluster](https://aws.amazon.com/ecs/) with [AWS RDS](https://aws.amazon.com/rds/)
and [AWS S3](https://aws.amazon.com/s3/).

More information about Datagrok design and components:

* [Architecture](architecture.md)
* [Infrastructure](infrastructure.md)

## Prerequisites

1. We use native Docker compose commands to run applications in Amazon ECS. It simplifies multi-container application
   development on Amazon ECS using familiar Compose files.
    1. Download and install the latest version of [Docker Desktop](https://docs.docker.com/desktop/), with Docker
       compose included, following installation instructions
       for [Windows](https://docs.docker.com/desktop/windows/install/)
       or [Mac](https://docs.docker.com/desktop/mac/install/). Or install
       the [Docker Compose CLI for Linux](https://docs.docker.com/cloud/ecs-integration/#install-the-docker-compose-cli-on-linux)
       which also enables Docker Compose features
    2. Check that you have [required permissions](https://docs.docker.com/cloud/ecs-integration/#requirements) on AWS
       account to perform Docker containers deployment to ECS
2. Additional components, such as database and storage, can be created using AWS CLI. To perform AWS CLI commands
   provided in the document
    1. [Install AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html)
    2. [Configure authorization for AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-quickstart.html)
3. Check that your default [VPC](https://docs.aws.amazon.com/vpc/latest/userguide/VPC_Subnets.html) has three subnets
   with [internet gateway routing](https://docs.aws.amazon.com/vpc/latest/userguide/VPC_Internet_Gateway.html) for
   internet facing Application Load Balancer. By default, default VPC already has all the required settings.
4. [Create S3 bucket](https://docs.aws.amazon.com/AmazonS3/latest/userguide/create-bucket-overview.html). For security
   reasons, we recommend to:
    * Disable public access to the bucket and objects inside
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

5. [Create RDS database](https://docs.aws.amazon.com/AmazonRDS/latest/UserGuide/USER_CreateDBInstance.html) for Datagrok
   with following parameters:
    * Required engine: PostgreSQL 12
    * Recommended DB instance class: db.t3.large
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
   > **_NOTE:_**  Datagrok provides demo databases with demo data for the full experience.
   > If you want to try datagrok with demo data choose
   > [Docker Compose Demo file](https://github.com/datagrok-ai/public/blob/master/docker/ecs.datagrok.demo.docker-compose.yaml)
   > instead. Rename downloaded file to `ecs.datagrok.docker-compose.yaml` to use commands in the instruction below
   > without any further changes

2. Replace `GROK_PARAMETERS` variables in `ecs.datagrok.docker-compose.yaml` with actual values or use environment
   variables automatic substitution

    ```shell
    #!/bin/sh

    export DATAGROK_S3_BUCKET_REGION=$AWS_REGION
    export DATAGROK_S3_BUCKET_NAME=example-unique-name-datagrok-s3
    export DATAGROK_RDS_ENDPOINT=$(aws rds describe-db-instances --db-instance-identifier "datagrok-rds" --output text --query 'DBInstances[].Endpoint.Address')
    ```

3. Use docker context for ECS: `docker context use ECS`
4. [Configure authorization for Docker ECS context as for AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-quickstart.html)
   and export required credentials as environment variables.

    ```shell
    #!/bin/sh
    export AWS_ACCESS_KEY_ID=<AWS_ACCESS_KEY_ID>
    export AWS_SECRET_ACCESS_KEY=<AWS_SECRET_ACCESS_KEY>
    export AWS_REGION=<AWS_REGION>
    ```

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
8. Datagrok starts to deploy the database immediately after startup. You can check the status by checking running task
   log in [CloudWatch](https://aws.amazon.com/cloudwatch/)

9. Switch back to default docker context: `docker context use default`

10. Create CNAME DNS record <DATAGROK_DNS> to the newly created Application Load Balancer

11. Go in the web browser to <DATAGROK_DNS>, login to Datagrok using username "admin" and password "admin"

## Setup Compute Virtual Machine

1. Download Docker Compose YAML
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
    * CVM URL Client: `http://<CVM_DNS>`
    * H2o Url: `http://<CVM_DNS>:54321`
    * API Url: `http://<DATAGROK_DNS>/api`
    * Cvm Split: `true`
* Dev:
    * CVM Url: `http://<CVM_DNS>`
    * Cvm Split: `true`
    * API Url: `http://<DATAGROK_DNS>/api`

## Additional information

This instruction will create all Datagrok components in the ECS cluster, including all related resources. However, the
resulting Application Load Balancers use HTTP instead of HTTPS, which is less secure. To improve security, create HTTPS
[Application Load Balancer listeners](https://docs.aws.amazon.com/elasticloadbalancing/latest/application/create-https-listener.html)
with SSL certificate, which matches newly created DNS names <CVM_DNS> and <DATAGROK_DNS>. And
then [change the HTTP listeners](https://docs.aws.amazon.com/elasticloadbalancing/latest/application/listener-update-rules.html)
to forward requests to the HTTPS ones.
