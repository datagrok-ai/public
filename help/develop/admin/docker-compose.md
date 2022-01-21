<!-- TITLE: Try Datagrok Locally -->
<!-- SUBTITLE: -->

# Try Datagrok Locally

This document contains instructions for running Datagrok on a machine
via [Docker Compose](https://docs.docker.com/compose/).

This method doesn't require cloud-based hosting. It automatically fetches, configures, and runs the required Docker
images.

If you want to jump-start with Datagrok on your local machine, we recommend this method.

## Prerequisites

1. [Docker Compose](https://docs.docker.com/compose/). If you do not have it, follow
   these [installation instructions](https://docs.docker.com/compose/install/) for your operating system.
2. Ideally, you should have at least 30 GB of free disk space, 2 CPU, 4 GiB RAM

## Instructions

1. Download a Docker Compose YAML
   file: [link](https://github.com/datagrok-ai/public/blob/master/docker/localhost.docker-compose.yaml).

2. To start up Datagrok, run this command:
   ```shell
   docker-compose -f docker/localhost.docker-compose.yaml --project-name datagrok --profile all up -d
   ```
   Datagrok will deploy a new database automatically.

   In case you get an error on Windows running `docker compose up` related to a `WriteFile`
   function, try running `cmd`
   in Administrator mode (this is a [known issue](https://github.com/docker/compose/issues/4531) of Docker on some
   computers).

3. Once the server is up and running, the Login page should be available
   at [http://localhost:8080](http://localhost:8080). For a quick setup, login to Datagrok using a username `admin`
   and a password `admin`.

4. After Datagrok is deployed for the first time, you can shut it down using command:
   ```shell
   docker-compose -f docker/localhost.docker-compose.yaml --project-name datagrok --profile all stop
   ```
   All the data will be saved in the [Docker volumes](https://docs.docker.com/storage/volumes/). If you want to reset
   Datagrok to factory settings, including all created users, projects, connections, etc., run the following command
   instead:
   ```shell
   docker-compose -f docker/localhost.docker-compose.yaml --project-name datagrok --profile all down --volumes
   ```

## Troubleshooting

1. In case of any issues, check the settings in the Datagrok (Tools | Settings...).
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

2. Check containers logs for any possible errors and report the problem if there is any
   ```shell
   docker-compose -f docker/localhost.docker-compose.yaml --project-name datagrok --profile all logs
   ```
3. Restart docker compose stand
    ```shell
    docker-compose -f docker/localhost.docker-compose.yaml --project-name datagrok --profile all stop
    docker-compose -f docker/localhost.docker-compose.yaml --project-name datagrok --profile all up -d
    ```
4. For advanced service troubleshooting you can access the containers shell
    ```shell
    # Replace <service> with one of the services: db, datagrok, grok_connect, grok_compute, jupyter_notebook, jupyter_kernel_gateway, h2o
    docker-compose -f docker/localhost.docker-compose.yaml --project-name datagrok --profile all exec <service> /bin/sh
    ```

## Advanced usage

If you do not need CVM features you can run only Datagrok application containers:
```shell
docker-compose -f docker/localhost.docker-compose.yaml --project-name datagrok --profile datagrok up -d
```

To run Datagrok with exact CVM features, specify them in the command line using `--profile` flag
* Cheminformatics
   ```shell
   docker-compose -f docker/localhost.docker-compose.yaml --project-name datagrok --profile datagrok --profile chem up -d
   ```
* Jupyter notebook
   ```shell
   docker-compose -f docker/localhost.docker-compose.yaml --project-name datagrok --profile datagrok --profile jupyter_notebook up -d
   ```
* Scripting
   ```shell
   docker-compose -f docker/localhost.docker-compose.yaml --project-name datagrok --profile datagrok --profile scripting up -d
   ```
* Modeling
   ```shell
   docker-compose -f docker/localhost.docker-compose.yaml --project-name datagrok --profile datagrok --profile modeling up -d
   ```
* Features can be enabled in any combination
   ```shell
   docker-compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
     --profile datagrok \
     --profile chem \
     --profile scripting \
     --profile jupyter_notebook \
     --profile modeling \
     up -d
   ```

## Useful links

* [Docker Compose](https://docs.docker.com/compose/)
* [Infrastructure](infrastructure.md)
* [Deployment](_deploy.md)
* [Configuration](configuration.md)
