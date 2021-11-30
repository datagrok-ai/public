<!-- TITLE: Deployment on a regular machine -->
<!-- SUBTITLE: -->

# Deployment on a regular machine

This document contains instructions to deploy Datagrok on a regular machine without cloud-based hosting.

In case you want to jump-start using Datagrok with minimum manual effort on a local machine,
check [Deployment with Docker Compose](docker-compose.md). The below method gives more control and includes manual
installation of a PostgreSQL instance along with locating Datagrok's working files on a host machine's file system.

## Prerequisites

1. [Install Docker](https://docs.docker.com/get-docker/).
2. [Install PostgreSQL](https://www.postgresql.org/download/), allow port 5432 in the firewall.
3. Create grok user with ID 1001 and group with ID 1001

## Setup Datagrok virtual machine

Requirements: 2 vCPU and 4 GiB RAM.

1. Get the latest Datagrok VM docker images from [Docker Hub](https://hub.docker.com/u/datagrok):

```bash
docker pull datagrok/datagrok:latest
docker pull datagrok/grok_connect:latest
```

2. Prepare JSON string `GROK_PARAMETERS`:

```json
{
   "dbServer": "host.docker.internal",      # Postgres database server (use host.docker.internal or 172.17.0.1 to connect to localhost)
   "dbPort": "5432",                        # Postgres database server port
   "db": "datagrok",                        # New database name
   "dbLogin": "datagrok",                   # New DB user name, Datagrok will use it to connect to Postgres database
   "dbPassword": "SoMeVeRyCoMpLeXpAsSwOrD", # New DB user password, Datagrok will use it to connect to Postgres database
   "dbAdminLogin": "postgres",              # Postgres admin login
   "dbAdminPassword": "postgres"            # Postgres admin password
}
```

3. Prepare local directory to store data: `GROK_DATA_PATH`, allow write access to the group grok.
4. Prepare local directory to store config files: `GROK_CFG_PATH`, allow write access to the group grok.
5. Create docker network for Datagrok: `docker network create datagrok`
6. Run Grok Connect image

```bash
docker run -it -d \
  --network datagrok \
  --network-alias grok_connect \
  --restart unless-stopped \
  datagrok/grok_connect:latest
```

7. Run Datagrok image. Wait for the deployment process to complete.

```bash
docker run -it -d \
  -e GROK_PARAMETERS="<GROK_PARAMETERS>" \
  -v <GROK_DATA_PATH>:/home/grok/data \
  -v <GROK_CFG_PATH>:/home/grok/cfg \
  --network datagrok \
  --network-alias datagrok \
  -p 8080:8080 \
  --restart unless-stopped \
  datagrok/datagrok:latest
```

7. Check if Datagrok started successfully: `http://HOST_NAME:8080`, login to Datagrok using username "`admin`" and
   password "`admin`".

8. Edit settings in the Datagrok (Tools | Settings...). Do not forget to click Apply to save new settings.
    * Connectors
        * External Host: `grok_connect`

## Setup Compute Virtual Machine

Requirements: 4 vCPU and 8 GiB RAM.

1. Get the latest Compute Virtual Machine docker images from [Docker Hub](https://hub.docker.com/u/datagrok):

```bash
docker pull datagrok/jupyter_kernel_gateway:latest
docker pull datagrok/jupyter_notebook:latest
docker pull datagrok/grok_compute:latest
docker pull datagrok/h2o:latest
docker pull datagrok/cvm_nginx:latest
```

2. Create docker network for CVM: `docker network create cvm`
3. Run CVM images:

```bash
# Grok Compute
docker run -it -d \
  -e GROK_COMPUTE_NUM_CORES=4 \
  --network cvm \
  --network-alias grok_compute \
  --restart unless-stopped \
  datagrok/grok_compute:latest
# H2O
docker run -it -d \
  --network cvm \
  --network-alias h2o \
  -p 54321:54321 
  -p 5005:5005 
  --restart unless-stopped \
  datagrok/h2o:latest
# JKG
docker run -it -d \
  --network cvm \
  --network-alias jupyter_kernel_gateway \
  --restart unless-stopped \
  datagrok/jupyter_kernel_gateway:latest
# JN
docker run -it -d \
  --network cvm \
  --network-alias jupyter_notebook \
  --restart unless-stopped \
  datagrok/jupyter_notebook:latest
# CVM Nginx
docker run -it -d \
  --network cvm \
  --network-alias cvm \
  -p 8090:8090 
  --restart unless-stopped \
  datagrok/cvm_nginx:latest
```

4. Edit settings in the Datagrok (Tools | Settings...). Do not forget to click Apply to save new settings.
    * Scripting:
        * CVM Url: `http://host.docker.internal:8090`
        * CVM Url Client: `http://localhost:8090`
        * H2o Url: `http://localhost:54321`
        * Api Url: `http://host.docker.internal:8080/api`
        * Cvm Split: `true`
    * Dev:
        * CVM Url: `http://localhost:8090`
        * Cvm Split: `true`
        * Api Url: `http://host.docker.internal:8080/api`

See also:

* [Architecture](architecture.md#application)
* [Architecture Details](architecture-details.md)
* [Compute VM](compute-vm.md)
* [Deployment with Docker Compose](docker-compose.md)
