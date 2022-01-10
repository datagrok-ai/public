<!-- TITLE: Deployment -->
<!-- SUBTITLE: -->

# Deployment

Datagrok consist of Docker containers which can be installed on any platform including but not limited to bare-metal
machine, on-premise virtual machine or virtual machine in cloud provider, for
example [AWS EC2](https://aws.amazon.com/ec2/), on-premise Kubernetes cluster or Kubernetes service in cloud provider,
for example [AWS EKS](https://aws.amazon.com/eks/), and container services in cloud provides, for
example [AWS ECS](https://aws.amazon.com/ecs/).

Also, Datagrok server requires [PostgreSQL database](infrastructure.md#database).
As [database](infrastructure.md#database)
Datagrok supports any PostgreSQL cluster out-of-the-box, including cloud solutions for PostgreSQL database, for
example [AWS RDS](https://aws.amazon.com/rds/). We recommend to use scalable and highly reliable solutions for database
and avoid single database instance setup to prevent datagrok internal information loss such as created users, created
connections, etc. User data won't be affected anyhow on Datagrok database crash.

For [persistent file storage](infrastructure.md#storage) Datagrok supports Local File System, Network shares or cloud
solutions, for example [AWS S3](https://aws.amazon.com/s3/) or [Google Cloud Storage](https://cloud.google.com/storage).
We recommend to use scalable and highly reliable solutions for storage and avoid local file system setup to prevent
datagrok internal information loss, such as projects, settings, etc. User data won't be affected anyhow on Datagrok
storage loss.

This document contains instructions to deploy Datagrok.

More information:

* [What is Datagrok?](../../home.md)
* [Architecture](architecture.md)
* [Infrastructure explanation](infrastructure.md)

## Prerequisites

### 1. Deployment Software

[Install Docker](https://docs.docker.com/get-docker/). Install [Docker Compose](https://docs.docker.com/compose/). If
you do not have it, follow these [installation instructions](https://docs.docker.com/compose/install/) for your
operating system.

### 2. Resources

Provide required [resources for every component](infrastructure.md#resources)

### 3. Database

Create [database](infrastructure.md#database) and [storage](infrastructure.md#storage), which should be available
for [Datlas](infrastructure.md#datlas). If you use local setup with local file storage and local database, skip this
step.

### 4. Storage

Create [database](infrastructure.md#database) and [storage](infrastructure.md#storage), which should be available
for [Datlas](infrastructure.md#datlas). If you use local setup with local file storage and local database, skip this
step.

### 5. Endpoints

Both [Compute](infrastructure.md#compute-components) and [Datagrok](infrastructure.md#datagrok-components) engines
should be accessible by users. The easiest way is to create DNS endpoints pointing to load balancers in front of the
services:
`datagrok.example` and `cvm.example`.

## Cloud providers deployment

Datagrok supports deployment to several cloud providers. To use cloud services for database, storage, load balancing,
DNS create all resources before Datagrok deployment.

For example, to deploy to AWS:

1. Create VPC with three subnets with internet gateway routing for internet facing Application Load Balancer
2. Configure S3 bucket and RDS database, which should be available from this VPC
3. Create Application Load Balancer which will point to instances or ECS cluster.
4. Create DNS records in Route53 which leads to Application Load Balancer
5. Start the Datagrok deployment

There are multiple Infrastructure as a Code examples can be found:

* [CloudFormation template](deploy/cloudformation.yaml) to deploy to AWS ECS
* [Terraform scripts](deploy/terraform.tf) to deploy to AWS ECS
* [Docker Compose for ECS deployment](deploy-amazon-ecs.md)

## Deployment

### Deploy script

The easiest way to deploy application is to use
our [deployment script](https://github.com/datagrok-ai/public/blob/master/help/develop/admin/deploy/deploy.sh)

1) Download script from
   repository: [deploy.sh](https://raw.githubusercontent.com/datagrok-ai/public/master/help/develop/admin/deploy/deploy.sh)
2) Run script, it will ask questions and deploy datagrok stand based on your answers. The supported deployment platform:
   Virtual Machine, Kubernetes, ECS.

```bash
sh deploy.sh
```

### Manual installation

1. Create SSH access to the host using SSH keys
2. Check that your user is in `docker` group
3. Create docker contexts:

```
docker context create --docker 'host=ssh://<DATAGROK_HOST_NAME>:22' datagrok`
docker context create --docker 'host=ssh://<CVM_HOST_NAME>:22' cvm`
```

4. Download Docker Compose yaml
   file: [link](https://github.com/datagrok-ai/public/blob/master/docker/localhost.docker-compose.yaml).
5. Prepare JSON string `GROK_PARAMETERS` and replace it in Docker Compose file. Full list of available options can be
   found in [Configuration](configuration.md).

For example, to use S3 as storage and RDS as database with Datagrok deployed on ec2 the result JSON will be

```
{
  \"amazonStorageRegion\": \"us-east-2\",
  \"amazonStorageBucket\": \"datagrok-test\",
  \"amazonStorageId\": \"ACCOUNTID\",
  \"amazonStorageKey\": \"SECRETKEY\",
  \"dbServer\": \"datagrok-db-1.abc.us-east-2.rds.amazonaws.com\",
  \"db\": \"datagrok\",
  \"dbLogin\": \"datagrok\",
  \"dbPassword\": \"SoMeVeryCoMpLeXpAsSwOrD\",
  \"dbAdminLogin\": \"postgres\",
  \"dbAdminPassword\": \"postgres\"
}
```

6. Switch to the datagrok context `docker context use datagrok`
7. Run Datagrok deploy. Wait for the deployment process to complete.
   `docker-compose --project-name datagrok --profile datagrok up -d`
8. Switch to the CVM context `docker context use cvm`
9. Run Datagrok deploy. Wait for the deployment process to complete.
   `docker-compose --project-name cvm --profile cvm up -d`
10. Check if Datagrok started successfully: `http://<DATAGROK_HOST_NAME>:8080`, login to Datagrok using
    username "`admin`" and password "`admin`".
11. Check settings in the Datagrok (Tools | Settings...). Do not forget to click Apply to save new settings if needed.

* Connectors
    * External Host: `grok_connect`
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
