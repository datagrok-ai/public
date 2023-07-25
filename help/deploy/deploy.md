---
title: "Deployment"
---

The deployment consists of a few docker containers, [database](../develop/under-the-hood/infrastructure.md#database) for storing metadata,
and [persistent file storage](../develop/under-the-hood/infrastructure.md#storage) for storing files.

Using Docker containers, you can deploy Datagrok on many environments, such as bare-metal machines,
on-premise virtual machines or virtual machines in cloud providers, for example [AWS EC2](https://aws.amazon.com/ec2/),
on-premise Kubernetes cluster or Kubernetes service in cloud providers, for
instance [AWS EKS](https://aws.amazon.com/eks/), and container services in the cloud providers, for
example [AWS ECS](https://aws.amazon.com/ecs/).

Datagrok requires PostgreSQL [database](../develop/under-the-hood/infrastructure.md#database) to store metadata.
We recommend using scalable and highly reliable solutions, such as [AWS RDS](https://aws.amazon.com/rds/).

For [persistent file storage](../develop/under-the-hood/infrastructure.md#storage) Datagrok supports Local File System, Network shares or cloud
solutions, for example [AWS S3](https://aws.amazon.com/s3/) or [Google Cloud Storage](https://cloud.google.com/storage).
We recommend using scalable and highly reliable solutions for storage, such as [AWS S3](https://aws.amazon.com/s3/).

This document contains different deployment options for Datagrok.

More information:

* [What is Datagrok?](../home.md)
* [Architecture](../develop/under-the-hood/architecture.md)
* [Infrastructure explanation](../develop/under-the-hood/infrastructure.md)

## Deployment options

Datagrok supports different deployment scenarios. You can choose the one which suits you best.

### Local deployment

Local deployment is a quick way to see Datagrok in action, and should not be used as a production environment.
[Try Datagrok Locally](docker-compose.mdx) for this case.

### Deploy script

The interactive way to deploy the platform is to use
our [deployment script](https://github.com/datagrok-ai/public/blob/master/help/develop/admin/deploy/deploy.sh)

1. Download the script from
   repository: [deploy.sh](https://raw.githubusercontent.com/datagrok-ai/public/master/help/develop/admin/deploy/deploy.sh)
2. For AWS deployment, check that you have
   all [required permissions](https://github.com/datagrok-ai/public/blob/master/help/develop/admin/deploy/iam.list)
   on AWS account and installed [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html) with your credentials
3. Run the script. It will ask questions and deploy a Datagrok stand based on your answers. The supported deployment
   platform:
   ECS, Kubernetes, Virtual Machine.

   * EC2 instance should be treated like Virtual Machine. It is required to create EC2 instances before the script run.
     You can check how to create instances
     in [regular machine example preparations steps](deploy-regular.md#preparations)

```bash
sh deploy.sh
```

### ECS deployment

We strongly recommend using [AWS ECS](https://aws.amazon.com/ecs/) for the Datagrok deployment. It provides a highly
scalable, fast container management service that makes it easy to manage application components. You can go completely
serverless using [AWS Fargate engine](https://aws.amazon.com/fargate/), which will reduce efforts on service support.

We prepared three options for effortless deployments to ECS:

1. [Docker Compose](deploy-amazon-ecs.md). It deploys the default application without enabled security features such as
   SSL and encryption. The deployment process is simple, and no foreknowledge is necessary. However, this option
   requires separate creation for an RDS database and S3 bucket.
2. [CloudFormation](deploy-amazon-cloudformation.md). It is an elaborate setup that considers all common security
   policies. It uses the CloudFormation template, which creates all required resources with enabled security options
   such as encryption and privileges control. Deployment can be provided on EC2 or FARGATE instances.
3. [Terraform](deploy-amazon-terraform.md). It is an elaborate setup that considers all common security
   policies. We developed terraform modules, which create all required resources with enabled security options
   such as encryption and privileges control. The modules are checked with checkov and tfsec to provide you the best-quality code. Deployment can be provided on EC2 or FARGATE instances.

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
scalable, and maintainable as others, so **we do not recommend it for production usage**. This deployment method will
require separate deployment for Datagrok required resources: PostgreSQL database and persistent storage(local filesystem
can be used).

More information:

* [Deployment on a regular machine](deploy-regular.md)
* [Deployment on AWS EC2](deploy-amazon-ec2.md)

Next steps:

* [Configure authentification](configure-auth.md)
* [Configure SMTP](configure-smtp.md)
