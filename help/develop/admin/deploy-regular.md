
<!-- TITLE: Deployment on a Regular Machine -->
<!-- SUBTITLE: -->

# Deployment on a Regular Machine

This document contains instructions to deploy Datagrok on a regular machine without cloud-based hosting.

In case you want to jump-start using Datagrok with minimum manual effort on a local machine, check [Deployment with Docker Compose](docker-compose.md). The below method gives more control and includes manual installation of a PostgreSQL instance along with locating Datagrok's working files on a host machine's file system.

## Prerequisites

1. [Install Docker](https://phoenixnap.com/kb/how-to-install-docker-on-ubuntu-18-04).
2. Install PostgreSQL, allow port 5432 in the firewall.

## Setup Datagrok Virtual Machine

Requirements: 2 vCPU and 4 GiB RAM.

1. Get the latest Datagrok Virtual Machine docker image from [Docker Hub](https://hub.docker.com/u/datagrok):
   ```
   docker pull datagrok/datagrok:VERSION
   ```
2. Prepare JSON string `GROK_START_PARAMETERS`:
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
4. Prepare local directory to store data: `GROK_DATA_PATH`, allow write access to the group.
5. Prepare local directory to store config files: `GROK_CFG_PATH`, allow write access to the group.
6. Run Datagrok image. Datagrok will create a database automatically.
   ```
   docker run -it -v <GROK_DATA_PATH>:/home/grok/data -v <GROK_CFG_PATH>:/home/grok/cfg -e GROK_PARAMETERS="<GROK_START_PARAMETERS>" -p 80:80 <IMAGE_NAME>
   ```
   Wait for the deployment process to complete
7. Check if Datagrok started successfully: `http://HOST_NAME`, login to Datagrok using username "`admin`" and password "`SM9ekKEkZuBDp5eD`".
   Do not use "`localhost`", use `127.0.0.1` or machine name instead. 

## Setup Compute Virtual Machine

Requirements: 4 vCPU and 8 GiB RAM.

1. Get the latest Compute Virtual Machine docker image from [Docker Hub](https://hub.docker.com/u/datagrok):
   ```
   docker pull datagrok/cvm:VERSION 
   ```
2. Run CVM image: `docker run -it -e GROK_COMPUTE_NUM_CORES=4 -p 8080:8080 -p 54321:54321 <IMAGE_NAME>`.

Edit settings in the Datagrok (_Tools | Settings..._):

* Scripting:
    * CVM: `http://CVM_DOCKER_URL:8080`
    * CVM Client: `http://CVM_URL:8080`
    * API Url: `http://DATAGROK_DOCKER_URL/api`

See also:

  * [Compute VM](../../compute/compute-vm.md)
  * [Architecture](architecture.md#application)
  * [Deployment with Docker Compose](docker-compose.md)
