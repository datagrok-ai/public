<!-- TITLE: Deployment with Docker Compose -->
<!-- SUBTITLE: -->

# Deployment with Docker Compose

This document contains instructions for running Datagrok on a regular machine
via [Docker Compose](https://docs.docker.com/compose/).

This method doesn't require cloud-based hosting. It automatically fetches, configures, and runs the required Docker
images.

If you want to jump-start with Datagrok on your local machine, we recommend this method. If you need to install manually
PostgreSQL and put Datagrok's working data on a host machine's file system,
check [Deployment on a regular machine](deploy-regular.md).

## Prerequisites

1. [Docker Compose](https://docs.docker.com/compose/). If you do not have it, follow
   these [installation instructions](https://docs.docker.com/compose/install/) for your operating system.
2. Ideally, you should have at least 30 GB of free disk space.

## Instructions

1. Download a Docker Compose YAML
   file: [link](https://github.com/datagrok-ai/public/blob/master/docker/localhost.docker-compose.yaml)
   .

2. To start up Datagrok, run this command:
   ```bash
   docker-compose --project-name datagrok --profile all up
   ```  
   Datagrok will deploy a new database automatically.

   In case you get an error on Windows running `docker compose up` related to a `WriteFile`
   function, try running `cmd`
   in Administrator mode (this is a [known issue](https://github.com/docker/compose/issues/4531) of Docker on some
   computers).

3. Once the server is up and running, the Login page should be available
   at [`http://localhost:8080`](http://localhost:8080). For a quick setup, login to Datagrok using a username `admin`
   and a password `admin`. To change your password, pass a key-value pair `"adminPassword": "yourPassword"` to the JSON
   string `GROK_PARAMETERS`.

4. After Datagrok is deployed for the first time, you can shut it down using `Ctrl+C`. Alternatively, run the command:
   ```bash
   docker-compose --project-name datagrok --profile all down
   ```  
   All the data will be saved in the persistent storage ([Docker volumes](https://docs.docker.com/storage/volumes/)). If
   you want to reset Datagrok to factory settings, run the following command instead:
   ```bash
   docker-compose --project-name datagrok --profile all down --volumes
   ```  
5. You may use the following commands to continue working with the existing containers comfortably:
   ```bash
   docker-compose --project-name datagrok --profile all up -d
   docker-compose --project-name datagrok --profile all stop
   ```
   Starting in a detached mode allows running the containers in the background, leaving out the logs, while the `stop`
   command, as opposed to `docker-compose down`, does not remove the network and the stopped containers after use.

6. Check the settings in the Datagrok (Tools | Settings...).

* Connectors
  * External Host: grok_connect
* Scripting:
  * CVM Url: `http://cvm:8090`
  * CVM Url Client: `http://localhost:8090`
  * H2o Url: `http://h2o:54321`
  * Api Url: `http://datagrok:8080/api`
  * Cvm Split: `true`
* Dev:
  * CVM Url: `http://localhost:8090`
  * Cvm Split: `true`
  * Api Url: `http://datagrok:8080/api`

See also:

* [Docker Compose](https://docs.docker.com/compose/)
* [Deployment on a regular machine](deploy-regular.md)
