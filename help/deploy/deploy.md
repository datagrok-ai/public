---
title: "Deployment"
---

To deploy Datagrok services, you can use [Docker containers](https://www.docker.com/resources/what-container/#:~:text=A%20Docker%20container%20image%20is,tools%2C%20system%20libraries%20and%20settings.). Datagrok consists of [core](../develop/under-the-hood/infrastructure.md#datagrok-components) containers, compute containers, a PostgreSQL [database](../develop/under-the-hood/infrastructure.md#database) to store metadata, and [persistent file storage](../develop/under-the-hood/infrastructure.md#storage) to store files.

Using Docker containers, you can deploy Datagrok on many environments, such as container services in the cloud providers, for example, [AWS ECS](#aws-deployment), [Kubernetes](#kubernetes-deployment), [bare-metal machines](#regular-machine-deployment), [virtual machines](#regular-machine-deployment), and so on.

To store data for Datagrok, we recommend using scalable and highly reliable solutions such as [AWS S3](https://aws.amazon.com/s3/) for persistent file storage and [AWS RDS](https://aws.amazon.com/rds/) as a PostgreSQL database.

## Local deployment

[Local deployment](docker-compose/docker-compose.md) is a quick way to see Datagrok in action using [Docker Compose](https://docs.docker.com/compose/). You can use it for local evaluation and development.

<!-- ### Deploy script

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
     in [regular machine example preparations steps](bare-metal/deploy-regular.md#preparations)

```bash
sh deploy.sh
```
--->

## AWS deployment

We strongly recommend using [AWS ECS](https://aws.amazon.com/ecs/) for the Datagrok deployment. It provides a highly
scalable, fast container management service that makes it easy to manage application components.

We prepared three options for effortless and secure deployments to AWS:

* [Marketplace](aws/deploy-marketplace.md). The easiest way to start with Datagrok on AWS. [Marketplace](https://aws.amazon.com/marketplace) deployment scripts create a separate infrastructure for Datagrok from scratch.
* [CloudFormation](aws/deploy-amazon-cloudformation.md). Using the [CloudFormation template](https://aws.amazon.com/cloudformation/), you can customize the Datagrok infrastructure with an elaborate template that considers all standard security policies.
* [Terraform](aws/deploy-amazon-terraform.md). It is the most flexible solution. You can integrate Datagrok into your existing infrastructure with consideration of your security policies. However, [Terraform](https://www.terraform.io/) is also an advanced option that requires additional knowledge in infrastructure as a code area.

## Kubernetes deployment

To deploy Datagrok to [Kubernetes](https://kubernetes.io/), we prepared [deployment scripts and ingress configuration](https://github.com/datagrok-ai/public/tree/master/help/deploy/k8s). It creates namespace and allocate all the necessary resources.

## Regular machine deployment

You can deploy Datagrok to a [regular machine](bare-metal/deploy-regular.md): bare-metal servers or virtual machines, including [EC2 instances](https://aws.amazon.com/ec2/). However, this method is less reliable, scalable, and maintainable than others. You need to set up hosts manually and manage the data storage. Consider using other options if possible.

## Complete the setup

After the deployment, open the platform to complete the setup:  

1. [Configure authentication](complete-setup/configure-auth.md)
2. [Configure SMTP](complete-setup/configure-smtp.md)
3. [Install packages](complete-setup/install-packages.md)
