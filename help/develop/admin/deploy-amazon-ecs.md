<!-- TITLE: Deploy Datagrok on AWS ECS cluster-->
<!-- SUBTITLE: -->

# Deploy Datagrok on AWS ECS cluster

This document contains instructions to deploy Datagrok on AWS ECS cluster.

## Prerequisites

1. Configure S3 bucket and RDS database

## Setup Datagrok Virtual Machine

1. Download latest Datagrok Virtual Machine docker image from [dev.datagrok.ai/docker_images](https://dev.datagrok.ai/docker_images) and upload it to your ECR
2. Create "Task Definition" with parameters:
    - Task Definition Name: datagrok
    - Task Role: ecsTaskExecutionRole
    - Network Mode: Host
    - Compatibilities: EC2
    - Task execution role: ecsTaskExecutionRole
    - Task memory (MiB): 2048
    - Task CPU (unit): 2048
    - Container Definitions:
        * Container name: datagrok
        * Image: datagrok:latest
        * Memory Limits (MiB): Soft Limit: 2048
        * Port mappings: 80
        * CPU units: 2048
        * Environment variables (remove comments and make single line):
            - GROK_MODE = deploy
            - GROK_PARAMETERS = see below
 ```GROK_PARAMETERS
{
"amazonStorageRegion": "us-east-2",                             # S3 region
"amazonStorageBucket": "datagrok-test",                         # S3 bucket name
"amazonStorageId": "ACCOUNTID",                                 # S3 credential ID, Datagrok will resolve ECS task role if empty
"amazonStorageKey": "SECRETKEY",                                # S3 credential secret key, Datagrok will resolve ECS tas if empty
"dbServer": "datagrok-db-1.abc.us-east-2.rds.amazonaws.com",    # RDS endpoint
"db": "datagrok_docker",                                        # RDS new database name
"dbLogin": "datagrok_docker",                                   # RDS new user name
"dbPassword": "SoMeCoMpLeXpAsSwOrD",                            # RDS new user password
"dbAdminLogin": "postgres",                                     # RDS admin login
"dbAdminPassword": "postgres"# RDS admin password
}
```

3. Create cluster (EC2 Linux + Networking):
    - Name: datagrok
    - EC2 instance type: t2.medium
    - Number of instances: 1
    - Key pair: <your ssh key pair>
    - Container instance IAM role: ecsInstanceRole
4. Create service:
    - Launch type: EC2
    - Task Definition: datagrok
    - Cluster: datagrok
    - Service name: datagrok
    - Service type: REPLICA
    - Number of tasks: 1
    - Deployment type: Rolling update

    Datagrok starts to deploy database immediately, you can check status by checking running task log

5. Create new revision for task definition
    - Container Definitions:
        * Modify environment variable
            - GROK_MODE = start
6. Update service, select new task definition
7. Go to the web browser, login to Datagrok using username "admin" and password "SM9ekKEkZuBDp5eD", open Tools | Settings.


## Setup Compute Virtual Machine

1. Download latest Compute Virtual Machine docker image from [dev.datagrok.ai/docker_images](https://dev.datagrok.ai/docker_images) and upload it to your ECR
2. Create "Task Definition" with parameters:
    - Task Definition Name: cvm
    - Task Role: ecsTaskExecutionRole
    - Network Mode: Host
    - Compatibilities: EC2
    - Task execution role: ecsTaskExecutionRole
    - Task memory (MiB): 4096
    - Task CPU (unit): 4096
    - Container Definitions:
        * Container name: grok_cvm
        * Image: grok_cvm:1.0.X-XXXXXXX
        * Memory Limits (MiB): Soft Limit: 4096
        * Port mappings: 80 and 54321
        * CPU units: 4096
        * Environment variables:
            - GROK_COMPUTE_NUM_CORES = 4
3. Create cluster (EC2 Linux + Networking):
    - Name: cvm
    - EC2 instance type: t2.medium
    - Number of instances: 2
    - Key pair: <your ssh key pair>
    - Container instance IAM role: ecsInstanceRole
4. Create load balancer (from EC2 console):
    - Load Balancers | Create Load Balancer
    - Application Load Balancer | Create
    - Name: cvm
    - Listeners for ports: 5005, 8004, 8888, 8889, 54321 or (80 and 54321)
    - Availability Zones: add all available zones
5. Create service:
    - Launch type: EC2
    - Task Definition: cvm
    - Cluster: cvm
    - Service name: cvm
    - Service type: DAEMON
    - Number of tasks: 2
    - Deployment type: Rolling update
    - Load Balancer: cvm

Edit settings in the Datagrok (Tools | Settings...):
* Scripting:
    * OpenCPU, OpenCPU Client: http://CVM_HOSTNAME/ocpu
    * Jupyter Notebook, Jupyter Notebook: http://CVM_HOSTNAME
    * Jupyter Gateway, Jupyter Gateway Client: http://CVM_HOSTNAME/jupyter
    * Grok Compute, Grok Compute Client: http://CVM_HOSTNAME/grok_compute
* Machine Learning:
    * H2O: http://CVM_HOSTNAME:54321

See also:
* [Compute VM](../../compute/compute-vm.md)
* [Architecture](architecture.md#application)
