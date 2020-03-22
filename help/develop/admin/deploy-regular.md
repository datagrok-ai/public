
<!-- TITLE: Deploy Datagrok on a regular machine -->
<!-- SUBTITLE: -->

# Deploy Datagrok on a regular machine

This document contains instructions to deploy Datagrok on a regular machine without cloud-based hosting.

## Prerequisites

1. [Install Docker](https://phoenixnap.com/kb/how-to-install-docker-on-ubuntu-18-04)
2. Install PostgreSQL, allow port 5432 in the firewall

## Setup Datagrok Virtual Machine

Requirements: 2 vCPU and 4GiB RAM

1. Download latest Datagrok Virtual Machine docker image from [dev.datagrok.ai/docker_images](https://dev.datagrok.ai/docker_images)
2. Import Docker image `docker load -i FILE.tar`
3. Prepare JSON string GROK_START_PARAMETERS:
 ```
{
\"dbServer\": \"host.docker.internal\",     # Postgres database server, Datagrok will use default 5432 port (use host.docker.internal or 172.17.0.1 to connect to localhost)
\"db\": \"datagrok\",                       # New database name
\"dbLogin\": \"datagrok\",                  # New DB user name, Datagrok will use it to connect to Postgres database
\"dbPassword\": \"SoMeCoMpLeXpAsSwOrD\",    # New DB user password, Datagrok will use it to connect to Postgres database
\"dbAdminLogin\": \"postgres\",             # Postgres admin login
\"dbAdminPassword\": \"postgres\"           # Postgres admin password
}
```
4. Prepare local directory to store data: GROK_DATA_PATH, allow access to write to group
5. Prepare local directory to store config files: GROK_CFG_PATH, allow access to write to group
6. Run Datagrok image. Datagrok will create database automatically.
`docker run -it -v <GROK_DATA_PATH>:/home/grok/data/prod -v <GROK_CFG_PATH>:/home/grok/cfg -e GROK_PARAMETERS="<GROK_START_PARAMETERS>" -p 80:80 <IMAGE_NAME>`,
wait for deploy process to complete
7. Check if Datagrok started successfully: http://HOST_NAME, login to Datagrok using username "admin" and password "SM9ekKEkZuBDp5eD".
 Do not use "localhost", use 127.0.0.1 or machine name instead. 

## Setup Compute Virtual Machine

Requirements: 4 vCPU and 8GiB RAM

1. Download latest Compute Virtual Machine docker image from [dev.datagrok.ai/docker_images](https://dev.datagrok.ai/docker_images)
2. Run CVM image `docker run -it -e GROK_COMPUTE_NUM_CORES=4 -p 5005:5005 -p 8004:8004 -p 8888:8888 -p 8889:8889 -p 54321:54321 <IMAGE_NAME>`

Edit settings in the Datagrok (Tools | Settings...):
* Scripting:
    * OpenCPU, OpenCPU client: http://CVM_HOST_NAME:8004/ocpu
    * Jupyter Notebook, Jupyter Notebook Client: http://CVM_HOST_NAME:8889
    * Jupyter Gateway, Jupyter Gateway Client: http://CVM_HOST_NAME:8888
    * Grok Compute, Grok Compute Client: http://CVM_HOST_NAME:5005
* Machine Learning:
    * H2O: http://HOST_NAME:54321

ComputeVM also contain Nginx service to reduce number of opened ports. For this case ports 80 and 54321
only are required: `docker run -it -e GROK_COMPUTE_NUM_CORES=1 -p 80:80 -p 54321:54321 <IMAGE_NAME>`

Settings:
* Scripting:
    * OpenCPU, OpenCPU client: http://CVM_HOST_NAME/ocpu
    * Jupyter Notebook, Jupyter Notebook Client: http://CVM_HOST_NAME
    * Jupyter Gateway, Jupyter Gateway Client: http://CVM_HOST_NAME/jupyter
    * Grok Compute, Grok Compute Client: http://CVM_HOST_NAME/grok_compute
* Machine Learning:
    * H2O: http://CVM_HOST_NAME:54321

See also:
* [Compute VM](../../compute/compute-vm.md)
* [Architecture](architecture.md#application)
