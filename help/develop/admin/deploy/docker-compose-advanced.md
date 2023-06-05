---
title: "Local machine: advanced"
slug: /develop/admin/docker-compose-advanced
sidebar_position: 1
---

This document contains instructions for running Datagrok on a local machine
via [Docker Compose](https://docs.docker.com/compose/).
This method doesn't require cloud-based hosting.
It automatically fetches, configures, and runs the required Docker images.
We recommend this method if you want to jump-start with Datagrok on your
local machine.

## Prerequisites

1. [Docker Compose](https://docs.docker.com/compose/).
   If you do not have it, follow
   these [installation instructions](https://docs.docker.com/compose/install/)
   for your operating system.
2. Minimal hardware requirements: 30 GB of free disk space, 2 CPUs, 4 GB RAM.
   Recommended hardware requirements: >30 GB of free disk space,
   4 CPUs, 8 GB RAM, or higher.

## Installing Datagrok

To install Datagrok manually, download a
[Docker Compose YAML file](https://raw.githubusercontent.com/datagrok-ai/public/master/docker/localhost.docker-compose.yaml)
and then run the following commands to spin up Datagrok or get the latest updates:

```shell
docker compose -f localhost.docker-compose.yaml --profile all pull
docker compose -f localhost.docker-compose.yaml --project-name datagrok --profile all up -d
```

Datagrok deploys a new database automatically.

:::note

   If you encounter an error related to a `WriteFile` function when running `docker-compose up` on Windows,
   try running the command prompt (`cmd`) in Administrator mode.
   This is a [known issue](https://github.com/docker/compose/issues/4531)
   with Docker on certain computers.

:::
  
After the docker compose process is completed, wait for approximately 1 minute for the Datagrok server to spin up.
Once the server is up and running, open `http://localhost:8080` in your browser.

### Login to Datagrok

On the login page enter the following credentials:

* username `admin`
* password `admin`

:::note

   If you see the message `Datagrok server is unavaliable` on the login page,
   wait for approximately 1 minute and reload the page.

:::

### Shutting down Datagrok

To shut down Datagrok use this command:

```shell
docker compose -f localhost.docker-compose.yaml --project-name datagrok --profile all stop
```

All the data is saved in the [Docker volumes](https://docs.docker.com/storage/volumes/).

To reset Datagrok to factory settings, including all created users, projects, connections, etc.,
run the following command:

```shell
docker compose -f localhost.docker-compose.yaml --project-name datagrok --profile all down --volumes
```

### CVM features

If you do not need CVM features, you can run only Datagrok application containers:

```shell
docker compose -f localhost.docker compose.yaml --project-name datagrok --profile datagrok --profile db up -d
```

If you need CVM features only, you can run only CVM application containers:

```shell
docker compose -f localhost.docker-compose.yaml --project-name datagrok --profile cvm up -d
```

To run Datagrok with exact CVM features, specify them in the command line using the `--profile` flag

* Cheminformatics

   ```shell
   docker compose -f localhost.docker-compose.yaml --project-name datagrok --profile datagrok --profile db --profile chem up -d
   ```

* Jupyter notebook

   ```shell
   docker compose -f localhost.docker-compose.yaml --project-name datagrok --profile datagrok --profile db --profile jupyter_notebook up -d
   ```

* Scripting

   ```shell
   docker compose -f localhost.docker-compose.yaml --project-name datagrok --profile datagrok --profile db --profile scripting up -d
   ```

* Modeling

   ```shell
   docker compose -f localhost.docker-compose.yaml --project-name datagrok --profile datagrok --profile db --profile modeling up -d
   ```

* Features can be enabled in any combination

   ```shell
   docker compose -f localhost.docker-compose.yaml --project-name datagrok \
     --profile datagrok \
     --profile db \
     --profile chem \
     --profile scripting \
     --profile jupyter_notebook \
     --profile modeling \
     up -d
   ```

* Datagrok container is not required to be started for any feature, so you can do not include in run parameters

   ```shell
   docker compose -f localhost.docker-compose.yaml --project-name datagrok \
     --profile chem \
     --profile scripting \
     --profile jupyter_notebook \
     --profile modeling \
     up -d
   ```

### Multiple stands

It is possible to run multiple stands of Datagrok on one host machine. To do so:

Run the first stand as described in [instruction](#installing-datagrok).

Set the Datagrok image version with the `DATAGROK_VERSION` environment variable. It can be any tag from
   [Docker Hub](https://hub.docker.com/r/datagrok/datagrok/tags). The default value is `latest`.

For Windows users:

 ```shell
 set DATAGROK_VERSION=latest
 ```

For Unix/Linux users:

 ```shell
 export DATAGROK_VERSION='latest'
 ```

Set environment variables for mapped ports:

| Environment variable       | Default value |
|----------------------------|---------------|
| DATAGROK_PORT              | 8080          |
| DATAGROK_DB_PORT           | 5432          |
| DATAGROK_CVM_PORT          | 8090          |
| DATAGROK_H2O_PORT          | 54321         |  
| DATAGROK_H2O_HELPER_PORT   | 5005          |
| DATAGROK_GROK_SPAWNER_PORT | 8000          |

To start the second stand properly the values should
differ from the ports of the existing stands. For example, you can increment every port value by 1.

For Windows users:

 ```shell
 set DATAGROK_PORT=8081
 set DATAGROK_DB_PORT=5433
 set DATAGROK_CVM_PORT=8091
 set DATAGROK_H2O_PORT=54322
 set DATAGROK_H2O_HELPER_PORT=5006
 set DATAGROK_GROK_SPAWNER_PORT=8001
 ```

For Linux/MacOS users:

 ```shell
 export DATAGROK_PORT='8081'
 export DATAGROK_DB_PORT='5433'
 export DATAGROK_CVM_PORT='8091'
 export DATAGROK_H2O_PORT='54322'
 export DATAGROK_H2O_HELPER_PORT='5006'
 export DATAGROK_GROK_SPAWNER_PORT='8001'
 ```

The last step is to run the second stand. It is important to change the project name to start the second Datagrok
stand. The project name in [the standard instructions](#installing-datagrok), which were used for the first stand,
is `datagrok`. For example, you can add an increment to the project name: `datagrok_2`

   ```shell
   docker compose -f localhost.docker-compose.yaml --project-name datagrok_2 --profile all up -d
   ```

### Demo Databases

Datagrok offers demo databases that include sample data.
To install them locally, run the following command:

```shell
docker compose -f localhost.docker-compose.yaml --project-name datagrok --profile all --profile demo up -d
```

## Troubleshooting

1. In case of any issues, check the settings in the Datagrok (Tools -> Settings...).

    * Connectors
      * External Host: `grok_connect`
    * Scripting:
      * CVM Url: `http://cvm:8090`
      * CVM URL Client: `http://localhost:8090`
      * H2o Url: `http://h2o:54321`
      * API Url: `http://datagrok:8080/api`
      * Cvm Split: `true`
    * Dev:
      * CVM Url: `http://localhost:8090`
      * Cvm Split: `true`
      * API Url: `http://datagrok:8080/api`

2. Check containers logs for any possible errors and report the problem if there are any

   ```shell
   docker compose -f localhost.docker-compose.yaml --project-name datagrok --profile all logs
   ```

   You can also watch the logs of the desired service in real-time

   ```shell
   # Get a list of running containers and choose a necessary one
   docker ps
   # Replace datagrok_datagrok_1 with the necessary container name
   # Replace 50 with desired log lines to watch or remove --tail 50 at all, if you need to
   # watch the full log
   docker logs -f --tail 50  datagrok_datagrok_1
   ```

3. Restart Docker compose stand

    ```shell
    docker compose -f localhost.docker-compose.yaml --project-name datagrok --profile all stop
    docker compose -f localhost.docker-compose.yaml --project-name datagrok --profile all up -d
    ```

4. For advanced service troubleshooting, you can access the containers shell

    ```shell
    # Replace <service> with one of the services: db, datagrok, grok_connect, grok_compute, jupyter_notebook, jupyter_kernel_gateway, h2o
    docker compose -f localhost.docker-compose.yaml --project-name datagrok --profile all exec <service> /bin/sh
    ```

5. Docker logs might take up all your free disk space. If such a situation has already taken place,
   add to the Docker daemon
   configuration log properties. For more about configuring Docker log options, you can refer
   to [the official documentation](https://docs.docker.com/config/containers/logging/local/#usage).

## Useful links

* [Docker Compose](https://docs.docker.com/compose/)
* [Infrastructure](../infrastructure.md)
* [Deployment](deploy.md)
* [Configuration](../configuration.md)
