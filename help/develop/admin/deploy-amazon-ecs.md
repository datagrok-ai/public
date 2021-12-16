<!-- TITLE: Deployment on AWS ECS Cluster-->
<!-- SUBTITLE: -->

# Deployment on AWS ECS Cluster

This document contains instructions to deploy Datagrok on [AWS ECS cluster](https://aws.amazon.com/ecs/).

## Prerequisites

1. Install [Docker Compose](https://docs.docker.com/compose/). If you do not have it, follow
   these [installation instructions](https://docs.docker.com/compose/install/) for your operating system.
2. Default VPC with three subnets with internet gateway routing for internet facing Application Load Balancer
3. Configure S3 bucket and RDS database, which should be available from default VPC
4. Docker ECS context: `docker context ecs create --from-env ECS`

## Setup Datagrok virtual machine

1. Download Docker Compose yaml
   file: [link](https://github.com/datagrok-ai/public/blob/master/docker/ecs.datagrok.docker-compose.yaml).
2. Replace `GROK_PARAMETERS` placeholders in `ecs.datagrok.docker-compose.yaml` with actual credentials
3. Use docker context for ECS:
   `docker context use ECS`
4. Export environment variables for AWS
   ```bash
   export AWS_REGION=<AWS_REGION>
   export AWS_ACCESS_KEY_ID=<AWS_ACCESS_KEY_ID>
   export AWS_SECRET_ACCESS_KEY=<AWS_SECRET_ACCESS_KEY>
   ```
5. Run docker compose deploy to ECS:
   `docker compose -p datagrok -f ecs.datagrok.docker-compose.yaml up -d`. It will create ECS cluster with tasks
   in [Fargate](https://aws.amazon.com/fargate/).

6. Check deployment result. Datagrok starts to deploy the database immediately after startup, you can check the status
   by checking running task log in
   [CloudWatch](https://aws.amazon.com/cloudwatch/)
   `docker compose -p datagrok -f ecs.datagrok.docker-compose.yaml ps`

7. Switch back to default docker context
   `docker context use default`

8. Create CNAME DNS record <DATAGROK_DNS> to the newly created Application Load Balancer

9. Go in the web browser to <DATAGROK_DNS>, login to Datagrok using username "admin" and password "admin"

10. Edit settings in the Datagrok (Tools | Settings...). Do not forget to click Apply to save new settings.
* Connectors
    * External Host: `grok_connect`

## Setup Compute Virtual Machine

1. Download Docker Compose yaml
   file: [link](https://github.com/datagrok-ai/public/blob/master/docker/ecs.cvm.docker-compose.yaml).
2. Use docker context for ECS:
   `docker context use ECS`
3. Export environment variables for AWS
   ```bash
   export AWS_REGION=<AWS_REGION>
   export AWS_ACCESS_KEY_ID=<AWS_ACCESS_KEY_ID>
   export AWS_SECRET_ACCESS_KEY=<AWS_SECRET_ACCESS_KEY>
   ```
4. Run docker compose deploy to ECS:
   `docker compose -p cvm -f ecs.cvm.docker-compose.yaml up -d`. It will create ECS cluster with tasks
   in [Fargate](https://aws.amazon.com/fargate/).

5. Check deployment results.   
   `docker compose -p cvm -f ecs.cvm.docker-compose.yaml ps`

6. Switch back to default docker context
   `docker context use default`

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

* [Architecture](architecture.md#application)
* [Architecture Details](architecture-details.md)
* [Compute VM](compute-vm.md)
