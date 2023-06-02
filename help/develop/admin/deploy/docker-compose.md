---
title: "Local machine: basic"
slug: /develop/admin/docker-compose
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

To download the installation script, right-click the link below 
and select "Save link as...".
Don't add any file extension.

   1. MacOS/Linux [Download installation script](https://raw.githubusercontent.com/datagrok-ai/public/master/docker/datagrok-install-local.sh). 
   2. Windows [Download installation script](https://raw.githubusercontent.com/datagrok-ai/public/master/docker/datagrok-install-local.cmd).

Once the script is downloaded and run the script to initiate the installation. For MacOS/Linux usage don't forget
make script executable by running `chmod +x localhost.docker-compose.yaml` before you start it.

:::note

   Docker daemon must be started and running before you run the script

:::

Once the installation process is complete,
proceed to the [login to Datagrok](#login-to-datagrok) section.

### Login to Datagrok

Once the server is up and running, open in your browser `http://localhost:8080` to see login page.
On the login page enter the following credentials:

* username `admin`
* password `admin`

:::note

   If you see the message `Datagrok server is unavaliable` on the login page,
   wait for approximately 1 minute and reload the page.

:::

### Shutting down Datagrok

To shut down Datagrok use this command from the directory in which you downloaded the script:

```shell
docker compose -f localhost.docker-compose.yaml --project-name datagrok --profile all stop
```

All the data is saved in the [Docker volumes](https://docs.docker.com/storage/volumes/).

## Useful links

* [Local machine deploy: advanced](docker-compose-advanced.md)
* [Docker Compose](https://docs.docker.com/compose/)
* [Infrastructure](../infrastructure.md)
* [Deployment](deploy.md)
* [Configuration](../configuration.md)
