<!-- TITLE: Deploy Datagrok on a regular machine -->
<!-- SUBTITLE: -->

# Deploy Datagrok on a regular machine

This document contains instructions to deploy Datagrok on a regular machine without cloud-based hosting.

## Prerequisites

1. [Install Docker](https://phoenixnap.com/kb/how-to-install-docker-on-ubuntu-18-04)
2. Install PostgreSQL

## Setup Datagrok Virtual Machine

Requirements: 2 vCPU and 4GiB RAM

1. Download latest Datagrok Virtual Machine docker image from [dev.datagrok.ai/docker_images](https://dev.datagrok.ai/docker_images)
2. Prepare JSON string GROK_START_PARAMETERS:
 ```
{
\"dbServer\": \"host.docker.internal\",     # Postgres database server (use host.docker.internal to connect to localhost)
\"db\": \"datagrok\",                       # New database name
\"dbLogin\": \"datagrok\",                  # New user name
\"dbPassword\": \"SoMeCoMpLeXpAsSwOrD\",    # New user password
\"dbAdminLogin\": \"postgres\",             # Poatgres admin login
\"dbAdminPassword\": \"postgres\"           # Postgres admin password
}
```
3. Prepare local directory to store data: GROK_DATA_PATH
4. Prepare local directory to store config files: GROK_CFG_PATH
5. Run Datagrok image in deploy mode. Datagrok will create database automatically.
`docker run -it -v <GROK_DATA_PATH>:/home/grok/data/prod -v <GROK_CFG_PATH>:/home/grok/cfg -e GROK_PARAMETERS="<GROK_START_PARAMETERS>" -e GROK_MODE=deploy -p 80:80 <IMAGE_NAME>`, will
wait for deploy process to complete
6. Run Datagrok image in regular mode. You can remove dbAdminLogin and dbAdminPassword parameters.
`docker run -it -v <GROK_DATA_PATH>:/home/grok/data/prod -v <GROK_CFG_PATH>:/home/grok/cfg -e GROK_PARAMETERS="<GROK_START_PARAMETERS>" -e GROK_MODE=start -p 80:80 <IMAGE_NAME>`
7. Check if Datagrok started successfully: http://localhost, login to Datagrok using username "admin" and password "SM9ekKEkZuBDp5eD"

## Setup Compute Virtual Machine

Requirements: 4 vCPU and 8GiB RAM

1. Download latest Compute Virtual Machine docker image from [dev.datagrok.ai/docker_images](https://dev.datagrok.ai/docker_images)
2. Run CVM image `docker run -it -e GROK_COMPUTE_NUM_CORES=4 -p 5005:5005 -p 8004:8004 -p 8888:8888 -p 8889:8889 -p 54321:54321 <IMAGE_NAME>`

Edit settings in the Datagrok (Tools | Settings...):
* Scripting:
    * OpenCPU, OpenCPU client: http://localhost:8004/ocpu
    * Jupyter Notebook, Jupyter Notebook Client: http://localhost:8889
    * Jupyter Gateway, Jupyter Gateway Client: http://localhost:8888
    * Grok Compute, Grok Compute Client: http://localhost:5005
* Machine Learning:
    * H2O: http://localhost:54321

ComputeVM also contain Nginx service to reduce number of opened ports. For this case ports 80 and 54321
only are required: `docker run -it -e GROK_COMPUTE_NUM_CORES=1 -p 80:80 -p 54321:54321 <IMAGE_NAME>`

Settings:
* Scripting:
    * OpenCPU, OpenCPU client: http://localhost/ocpu
    * Jupyter Notebook, Jupyter Notebook Client: http://localhost
    * Jupyter Gateway, Jupyter Gateway Client: http://localhost/jupyter
    * Grok Compute, Grok Compute Client: http://localhost/grok_compute
* Machine Learning:
    * H2O: http://localhost:54321

See also:
* [Compute VM](../../compute/compute-vm.md)
* [Architecture](architecture.md#application)
