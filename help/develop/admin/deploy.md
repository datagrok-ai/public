<!-- TITLE: Deployment -->
<!-- SUBTITLE: -->

# Overview

Datagrok consist of Docker containers, [database](infrastructure.md#database)
and [persistent file storage](infrastructure.md#storage).

Docker containers allow installing Datagrok on any platform, including but not limited to bare-metal machines,
on-premise virtual machines or virtual machines in cloud providers, for example [AWS EC2](https://aws.amazon.com/ec2/),
on-premise Kubernetes cluster or Kubernetes service in cloud providers, for
instance [AWS EKS](https://aws.amazon.com/eks/), and container services in cloud provides, for
example [AWS ECS](https://aws.amazon.com/ecs/).

Also, Datagrok server requires [PostgreSQL database](infrastructure.md#database).
As [database](infrastructure.md#database)
Datagrok supports any PostgreSQL database out-of-the-box, including cloud solutions for PostgreSQL database, for
example [AWS RDS](https://aws.amazon.com/rds/). We recommend using scalable and highly reliable solutions for databases
and avoiding single database instance setup to prevent datagrok internal information loss such as created users, created
connections, etc. User data won't be affected anyhow on Datagrok database crash.

For [persistent file storage](infrastructure.md#storage) Datagrok supports Local File System, Network shares or cloud
solutions, for example [AWS S3](https://aws.amazon.com/s3/) or [Google Cloud Storage](https://cloud.google.com/storage).
We recommend using scalable and highly reliable solutions for storage and avoiding local file system setup to prevent
datagrok internal information loss, such as projects, settings, etc. User data won't be affected anyhow on Datagrok
storage loss.

This document contains different deployment options for Datagrok.

More information:

* [What is Datagrok?](../../home.md)
* [Architecture](architecture.md)
* [Infrastructure explanation](infrastructure.md)

## Deployment options

Datagrok supports different deployment scenarios. You can choose the one which suits you best.

### Deploy script

The easiest way to deploy the application is to use
our [deployment script](https://github.com/datagrok-ai/public/blob/master/help/develop/admin/deploy/deploy.sh)

1. Download script from
   repository: [deploy.sh](https://raw.githubusercontent.com/datagrok-ai/public/master/help/develop/admin/deploy/deploy.sh)
2. For AWS deployment, check that you have
   all [required permissions](https://github.com/datagrok-ai/public/blob/master/help/develop/admin/deploy/iam.list)
   on AWS account
3. Run the script. It will ask questions and deploy a datagrok stand based on your answers. The supported deployment
   platform:
   ECS, Kubernetes, Virtual Machine.

   * EC2 instance should be treated like Virtual Machine. It is required to create EC2 instances before the script run.
     You can check how to create instances
     in [regular machine example preparations steps](deploy-regular.md#preparations)

```bash
sh deploy.sh
```

### Local deployment

This method does not require anything except one machine on which docker-compose will deploy all the components for
Datagrok. We recommend this method if you want to jump-start with Datagrok on your local device.

This is a reference implementation that is NOT recommended for production usage.

More information:

* [Try Datagrok Locally](docker-compose.md)

### ECS deployment

We strongly recommend using [AWS ECS](https://aws.amazon.com/ecs/) for the Datagrok deployment. It provides a highly
scalable, fast container management service that makes it easy to manage application components. You can go completely
serverless using [AWS Fargate engine](https://aws.amazon.com/fargate/), which will reduce efforts on service support.

We prepared two options for effortless deployments to ECS:

1. [Docker Compose](deploy-amazon-ecs.md). It deploys the default application without enabled security features such as
   SSL and encryption. The deployment process is simple, and no foreknowledge is necessary. However, this option
   requires separate creation for an RDS database and S3 bucket.
2. [CloudFormation](deploy-amazon-cloudformation.md). It is an elaborate setup that considers all common security
   policies. It uses the CloudFormation template, which creates all required resources with enabled security options
   such as encryption and privileges control.
3. [Terraform](deploy-amazon-terraform.md). It is an elaborate setup that considers all common security
   policies. We developed terraform modules, which create all required resources with enabled security options
   such as encryption and privileges control. The modules are checked with checkov and tfsec to provide you the best
   quality code.

### Kubernetes deployment

[Kubernetes](https://kubernetes.io/) is software for scaling and managing applications in Docker containers. It is an
extensive system with complex architecture. Kubernetes provides you with a framework to run distributed systems
resiliently.

Datagrok requires persistent storage, which can be created using Kubernetes persistent volumes. PostgreSQL database can
be both deployed in Kubernetes or on any other server. Set the database address and credentials
in [Datagrok configuration](configuration.md). Datagrok server will deploy all required schemas and users on the
startup.

To deploy Datagrok in Kubernetes, including Cloud solutions for Kubernetes, such
as [AWS EKS](https://aws.amazon.com/eks/), you can use our prepared deployment scripts and ingress configuration.

### Regular machine deployment

Datagrok can be deployed to a regular machine: bare-metal servers or virtual machines. This method is not as reliable,
scalable, and maintainable as others, so we do not recommend it for production usage. This deployment method will
require separate deployment for Datagrok required resources: PostgreSQL database and persistent storage(local filesystem
can be used).

We use native Docker compose commands to run applications on machines. It simplifies multi-container application
development and deployment.

More information:

* [Deployment on a regular machine](deploy-regular.md)
* [Deployment on AWS EC2](deploy-amazon-ec2.md)
