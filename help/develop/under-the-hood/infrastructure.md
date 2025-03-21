---
title: "Infrastructure"
---

![Infrastructure diagram](architecture-diagram.drawio.png)

## Datagrok Platform Infrastructure Overview

The Datagrok platform consists of multiple interconnected services that support data analysis, scripting, plugin management, 
and external integrations. The infrastructure includes the following key components:

### 1. Core Components
- **Datlas (Main Web Server)**: The central brain of the Datagrok platform, responsible for managing user interactions, processing requests, and orchestrating services.
- [**Nginx (Reverse Proxy)**](https://www.nginx.com/): Routes incoming web traffic to the server, allowing multiple instances of Datlas to run.
- [**PostgreSQL**](https://www.postgresql.org/): Stores platform metadata, including user accounts, projects, and configurations.
- **[S3](https://docs.aws.amazon.com/AmazonS3/latest/userguide/Welcome.html)/[Google Cloud Storage](https://cloud.google.com/storage/docs)/Local Storage**: Stores user files and tabular data used for analysis.

---

### 2. Authentication & Credential Management
- **Credentials Service**: Manages secure storage and retrieval of user credentials or credentials to external data sources.
- **External Credentials Storage**: A separate service that securely stores enterprise authentication credentials.

---

### 3. External Database Connectivity
- [**Grok Connect**](../../access/access.md#data-connection): A dedicated in-house service for integrating with external databases (PostgreSQL, Oracle, MySQL, MSSQL, etc.).

---

### 4. Scripting and Computation
- **AMQP**: Handles the queuing of script execution requests. [RabbitMQ](https://www.rabbitmq.com/docs) is typically used, but can be replaced with any message queue that implements AMQP protocol.
- **Scripting Workers**: A set of workers that process user scripts:
  - [**Jupyter Kernel-Based Worker**](https://docs.jupyter.org/en/stable/projects/kernels.html): Executes Python, R, Julia, Octave, and Node.js scripts.
  - **Custom Workers**: Additional computation engines with specialized capabilities that can be created on client requests.
- **Grok Pipe (Data Transfer Service)**: In-house service that handles binary data exchange between Datlas and Jupyter workers.

---

### 5. Plugin & Docker Container Management
- **Grok Spawner**: Manages the deployment of Docker Containers that are delivered within [Datagrok Plugins](../packages/extensions).
- **Plugin Docker Containers**: There are various specialized containers for AI, Bioinformatics, Machine Learning, Cheminformatics and Jupyter Notebooks that are supplied within Datagrok Plugins.

---

## Resources

There are different requirements for every component. In general:

* Core Components require 2 vCPU and 4 GB RAM.
  * For the active usage of Datagrok we would recommend: 4 vCPU and 8 GB RAM
* Additional components require at least 4 vCPU and 8 GB RAM.
  * For the active usage of Datagrok we would recommend: 8 vCPU and 32 GB RAM

## Deployment

Datagrok is a web application, which means no deployment efforts per user once the server is set up. All following
administration tasks could be performed via the web interface as well.

Enterprises typically prefer on-premise deployment for multiple reasons, such as security, ability to easily access
internal data, and other features such as integration with the enterprise
[authentication](../../govern/access-control/access-control.md#authentication) methods. Regarding Datagrok infrastructure it can be easily done. For
more information check [Enterprise Evaluation FAQ](../../datagrok/solutions/enterprise/enterprise-evaluation-faq.md)
page.

Datagrok consist of [Docker containers](https://hub.docker.com/u/datagrok) which can be installed on any platform
including but not limited to bare-metal
machine, on-premise virtual machine or virtual machine in cloud provider, for
example [AWS EC2](https://aws.amazon.com/ec2/), on-premise Kubernetes cluster or Kubernetes service in cloud provider,
for example [AWS EKS](https://aws.amazon.com/eks/), and container services in cloud provides, for
example [AWS ECS](https://aws.amazon.com/ecs/).

As [database](#database) Datagrok supports any PostgreSQL database out-of-the-box, including cloud solutions for
PostgreSQL database, for example [AWS RDS](https://aws.amazon.com/rds/). We recommend to use scalable and highly
reliable solutions for database and avoid single database instance setup to prevent datagrok internal information loss
such as created users, created connections, etc. User data won't be affected anyhow on Datagrok database crash.

For [persistent file storage](#storage) Datagrok supports Local File System, Network shares or cloud solutions, for
example [AWS S3](https://aws.amazon.com/s3/) or [Google Cloud Storage](https://cloud.google.com/storage). We recommend
to use scalable and highly reliable solutions for storage and avoid local file system setup to prevent datagrok internal
information loss, such as projects, settings, etc. User data won't be affected anyhow on Datagrok storage loss.

Check [How to deploy datagrok?](../../deploy/deploy.md) for details.

## Scalability

All Docker containers are based on Ubuntu 20.04 and use the latest software available. Datagrok can be deployed in any cloud
or on a regular machine. Components are easily scalable using cloud container services, for example [Amazon ECS](https://aws.amazon.com/ecs/),
or [Kubernetes](https://kubernetes.io/) services. We provide full support to make application work stable under concurrent user load.


## Useful links

* [Deployment](../../deploy/deploy.md)
* [Configuration](../../deploy/configuration.md)
* [Continuous integration](continuous-integration.png)
* [Versioning policy](../dev-process/versioning-policy.md)
* [Try Datagrok locally](../../deploy/docker-compose/docker-compose.md)
