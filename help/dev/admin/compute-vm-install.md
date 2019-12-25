<!-- TITLE: GrokCompute Installation -->
<!-- SUBTITLE: -->

# GrokCompute Installation

This document contains instructions for installation of the [GrokCompute](architecture.md#grokcompute),
platform's advanced computational services. See [Datagrok installation](datagrok-install.md) 
for instructions on deploying platform core.  

## Setup Docker image locally

Minimum requirements:
* CPU: 4x
* RAM: 8GB

For AWS usual case is "c5.xlarge".

Docker commands:

```bash
docker pull grok_cvm:1.0.X-XXXXXXX
docker run -it -e GROK_COMPUTE_NUM_CORES=1 \
        -p 5005:5005 -p 8004:8004 -p 8888:8888 -p 8889:8889 -p 54321:54321 \
        grok_cvm:1.0.X-XXXXXXX
```

Edit settings in the platform (Main menu | Tools | Settings...):
* Dev: 
    * OpenCPU: http://localhost:8004/ocpu
    * Jupyter Notebook: http://localhost:8889
    * Jupyter Gateway: http://localhost:8888
    * Grok Compute: http://localhost:5005
* Machine Learning: 
    * H2O: http://localhost:54321
    
GrokCompute also contain Nginx service to reduce number of opened ports. For this case ports 80 and 54321 
only are required. 

```bash
docker pull grok_cvm:1.0.X-XXXXXXX
docker run -it -e GROK_COMPUTE_NUM_CORES=1 -p 80:80 -p 54321:54321 grok_cvm:1.0.X-XXXXXXX
```

Settings:
* Dev: 
    * OpenCPU: http://localhost/ocpu
    * Jupyter Notebook: http://localhost
    * Jupyter Gateway: http://localhost/jupyter
    * Grok Compute: http://localhost/grok_compute
* Machine Learning: 
    * H2O: http://localhost:54321


## AWS bare EC2 installation

1. Get access to CVM docker images
2. [Install Docker](https://phoenixnap.com/kb/how-to-install-docker-on-ubuntu-18-04)
3. [Install AWS console](https://linuxhint.com/install_aws_cli_ubuntu/) 
4. Login to Amazon `${aws ecr get-login --no-include-email --region us-east-2 --registry-ids 766822877060})`
5. Download CVM image `docker pull 766822877060.dkr.ecr.us-east-2.amazonaws.com/grok_cvm:1.0.X-XXXXXXX`    
6. Run CVM image     
`docker run -it -e GROK_COMPUTE_NUM_CORES=4 -p 80:80 -p 54321:54321 766822877060.dkr.ecr.us-east-2.amazonaws.com/grok_cvm:1.0.X-XXXXXXX`


## AWS installation

Setup: 
    - EC2: c5.xlarge x2 (4 vCPU and 8GiB)
    - 2 instances of Compute VM

1. Get access to Grok console.aws.amazon.com/ecr/repositories/grok_cvm
1. Create "Task Definition" with parameters:
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
        * Port mappings: 5005, 8004, 8888, 8889, 54321 or (80 and 54321)
        * CPU units: 4096
        * Environment variables: 
            - GROK_COMPUTE_NUM_CORES = 4
1. Create cluster (EC2 Linux + Networking):
    - Name: cvm
    - EC2 instance type: t2.medium
    - Number of instances: 2
    - Key pair: <your ssh key pair>
    - Container instance IAM role: ecsInstanceRole
1. Create load balancer (from EC2 console):
    - Load Balancers | Create Load Balancer
    - Application Load Balancer | Create
    - Name: cvm
    - Listeners for ports: 5005, 8004, 8888, 8889, 54321 or (80 and 54321)
    - Availability Zones: add all available zones
1. Create service:
    - Launch type: EC2
    - Task Definition: cvm
    - Cluster: cvm
    - Service name: cvm
    - Service type: DAEMON
    - Number of tasks: 2
    - Deployment type: Rolling update
    - Load Balancer: cvm


See also: 
* [Compute VM](../../features/compute-vm.md)    
