---
title: "Local machine: advanced"
sidebar_position: 1
format: 'mdx'
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

This article provides additional options for running Datagrok on your local machine
using [Docker Compose](https://docs.docker.com/compose/). By following these instructions, you can customize your local
Datagrok stand or even start the second one.

:::info Hardware requirements

Minimal hardware requirements: 60 GB of free disk space, 4 CPUs, 8 GB RAM.

:::

## Prerequisites

1. Install and launch the latest Docker Desktop application for your operating system:

    - [Docker Desktop for Windows](https://docs.docker.com/desktop/install/windows-install/)
    - [Docker Desktop for MacOS](https://docs.docker.com/desktop/install/mac-install/)
    - [Docker Desktop for Linux](https://docs.docker.com/desktop/install/linux-install/)
2. [Clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)
   public [repository](https://github.com/datagrok-ai/public) with docker-compose file
3. Open the command-line interface and navigate to the directory where you cloned the repository.

## Installing Datagrok

1. Get latest Datagrok docker images

   <Tabs groupId="os" queryString>
   <TabItem value="win" label="Windows" default>

   ```cmd
   docker compose -f docker\localhost.docker-compose.yaml --profile all pull
   ```

   </TabItem>
   <TabItem value="bash" label="MacOS/Linux">

   ```shell
   docker compose -f docker/localhost.docker-compose.yaml --profile all pull
   ```

   </TabItem>
   </Tabs>

2. Run the basic Datagrok stand

   <Tabs groupId="os" queryString>
   <TabItem value="win" label="Windows" default>

   ```cmd
   docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok ^
     --profile all up -d
   ```

   :::note

   If you encounter an error related to a `WriteFile` function when running `docker-compose up` on Windows,
   try running the command prompt (`cmd`) in Administrator mode.
   This is a [known issue](https://github.com/docker/compose/issues/4531)
   with Docker on certain computers.

   :::

   </TabItem>
   <TabItem value="bash" label="MacOS/Linux">

   ```shell
   docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
     --profile all up -d
   ```

   </TabItem>
   </Tabs>

3. After the docker compose process is completed, wait for approximately 1 minute for the Datagrok server to spin up.
   Once the server is up and running, open [http://localhost:8080](http://localhost:8080) in your browser
   and [log in](#log-in-to-datagrok).

### Log in to Datagrok

1. Once the server is up and running, open your browser and go to [http://localhost:8080](http://localhost:8080).
2. On the login page, use the following credentials to login:
    - Login or Email: `admin`
    - Password `admin`

:::note

If you see the message `Datagrok server is unavaliable`, wait for approximately 1 minute for the server to start, and
then reload the page.

:::

### CVM features

If you do not need CVM features, you can run only Datagrok application containers to save space and resources:

<Tabs groupId="os" queryString>
<TabItem value="win" label="Windows" default>

```cmd
docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok ^
  --profile datagrok --profile db up -d
```

</TabItem>
<TabItem value="bash" label="MacOS/Linux">

```shell
docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
  --profile datagrok --profile db up -d
```

</TabItem>
</Tabs>

If you need CVM features only, you can run only CVM application containers:

<Tabs groupId="os" queryString>
<TabItem value="win" label="Windows" default>

```cmd
docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok ^
  --profile cvm up -d
```

</TabItem>
<TabItem value="bash" label="MacOS/Linux">

```shell
docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
  --profile cvm up -d
```

</TabItem>
</Tabs>

To run Datagrok with exact CVM features, specify them in the command line using the `--profile` flag

- Cheminformatics

   <Tabs groupId="os" queryString>
   <TabItem value="win" label="Windows" default>

   ```cmd
   docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok ^
    --profile datagrok --profile db --profile chem up -d
   ```

   </TabItem>
   <TabItem value="bash" label="MacOS/Linux">

   ```shell
   docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
    --profile datagrok --profile db --profile chem up -d
   ```

   </TabItem>
   </Tabs>

- Jupyter notebook

   <Tabs groupId="os" queryString>
   <TabItem value="win" label="Windows" default>

   ```cmd
   docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok ^
    --profile datagrok --profile db --profile jupyter_notebook up -d
   ```

   </TabItem>
   <TabItem value="bash" label="MacOS/Linux">

   ```shell
   docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
    --profile datagrok --profile db --profile jupyter_notebook up -d
   ```

   </TabItem>
   </Tabs>

- Scripting

   <Tabs groupId="os" queryString>
   <TabItem value="win" label="Windows" default>

   ```cmd
   docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok ^
    --profile datagrok --profile db --profile scripting up -d
   ```

   </TabItem>
   <TabItem value="bash" label="MacOS/Linux">

   ```shell
   docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
    --profile datagrok --profile db --profile scripting up -d
   ```

   </TabItem>
   </Tabs>

- Modeling

   <Tabs groupId="os" queryString>
   <TabItem value="win" label="Windows" default>

   ```cmd
   docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok ^
    --profile datagrok --profile db --profile modeling up -d
   ```

   </TabItem>
   <TabItem value="bash" label="MacOS/Linux">

   ```shell
   docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
    --profile datagrok --profile db --profile modeling up -d
   ```

   </TabItem>
   </Tabs>

- Features can be enabled in any combination

   <Tabs groupId="os" queryString>
   <TabItem value="win" label="Windows" default>

   ```cmd
   docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok ^
     --profile datagrok ^
     --profile db ^
     --profile chem ^
     --profile scripting ^
     --profile jupyter_notebook ^
     --profile modeling ^
     up -d
   ```

   </TabItem>
   <TabItem value="bash" label="MacOS/Linux">

   ```shell
   docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
     --profile datagrok \
     --profile db \
     --profile chem \
     --profile scripting \
     --profile jupyter_notebook \
     --profile modeling \
     up -d
   ```

   </TabItem>
   </Tabs>

- Datagrok container is not required to be started for any feature, so you can omit it in run parameters

   <Tabs groupId="os" queryString>
   <TabItem value="win" label="Windows" default>

   ```cmd
   docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok ^
     --profile chem ^
     --profile scripting ^
     --profile jupyter_notebook ^
     --profile modeling ^
     up -d
   ```

   </TabItem>
   <TabItem value="bash" label="MacOS/Linux">

   ```shell
   docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
     --profile chem \
     --profile scripting \
     --profile jupyter_notebook \
     --profile modeling \
     up -d
   ```

   </TabItem>
   </Tabs>

### Multiple stands

It is possible to run multiple stands of Datagrok on one host machine. To do so:

1. Run the first stand as described in [instruction](#installing-datagrok).

2. Set the Docker images versions with the environment variables. It can be any tag from
   [Docker Hub](https://hub.docker.com/r/datagrok/).

   | Environment variable | Default value |
      |----------------------|---------------|
   | DATAGROK_VERSION     | latest        |
   | GROK_SPAWNER_VERSION | latest        |
   | GROK_CONNECT_VERSION | latest        |
   | GROK_COMPUTE_VERSION | latest        |  

   <Tabs groupId="os" queryString>
   <TabItem value="win" label="Windows" default>

   ```cmd
    set DATAGROK_VERSION=latest
    set GROK_SPAWNER_VERSION=latest
    set GROK_CONNECT_VERSION=latest
    set GROK_COMPUTE_VERSION=latest
   ```

   </TabItem>
   <TabItem value="bash" label="MacOS/Linux">

   ```shell
    export DATAGROK_VERSION='latest'
    export GROK_SPAWNER_VERSION='latest'
    export GROK_CONNECT_VERSION='latest'
    export GROK_COMPUTE_VERSION='latest'
   ```

   </TabItem>
   </Tabs>

3. Set environment variables for mapped ports:

   | Environment variable                  | Default value |
   |---------------------------------------|---------------|
   | DATAGROK_PORT                         | 8080          |
   | DATAGROK_DB_PORT                      | 5432          |
   | DATAGROK_CVM_PORT                     | 8090          |
   | DATAGROK_H2O_PORT                     | 54321         |  
   | DATAGROK_H2O_HELPER_PORT              | 5005          |
   | GROK_SPAWNER_PORT                     | 8000          |
   | GROK_CONNECT_PORT                     | 1234          |
   | DATAGROK_DEMO_POSTGRES_NORTHWIND_PORT | 5433          |
   | DATAGROK_DEMO_POSTGRES_CHEMBL_PORT    | 5434          |
   | DATAGROK_DEMO_POSTGRES_UNICHEM_PORT   | 5435          |
   | DATAGROK_DEMO_POSTGRES_STARBUCKS_PORT | 5436          |
   | DATAGROK_DEMO_POSTGRES_WORLD_PORT     | 5437          |

   To start the second stand properly the values should
   differ from the ports of the existing stands. For example, you can increment every port value by 11.

   <Tabs groupId="os" queryString>
   <TabItem value="win" label="Windows" default>

   ```cmd
   set DATAGROK_PORT=8091
   set DATAGROK_DB_PORT=5444
   set DATAGROK_CVM_PORT=8101
   set DATAGROK_H2O_PORT=54332
   set DATAGROK_H2O_HELPER_PORT=5016
   set GROK_SPAWNER_PORT=8011
   set GROK_CONNECT_PORT=1245
   set DATAGROK_DEMO_POSTGRES_NORTHWIND_PORT=5444
   set DATAGROK_DEMO_POSTGRES_CHEMBL_PORT=5445
   set DATAGROK_DEMO_POSTGRES_UNICHEM_PORT=5446
   set DATAGROK_DEMO_POSTGRES_STARBUCKS_PORT=5447
   set DATAGROK_DEMO_POSTGRES_WORLD_PORT=5448
   ```

   </TabItem>
   <TabItem value="bash" label="MacOS/Linux">

   ```shell
   export DATAGROK_PORT=8091
   export DATAGROK_DB_PORT=5444
   export DATAGROK_CVM_PORT=8101
   export DATAGROK_H2O_PORT=54332
   export DATAGROK_H2O_HELPER_PORT=5016
   export GROK_SPAWNER_PORT=8011
   export GROK_CONNECT_PORT=1245
   export DATAGROK_DEMO_POSTGRES_NORTHWIND_PORT=5444
   export DATAGROK_DEMO_POSTGRES_CHEMBL_PORT=5445
   export DATAGROK_DEMO_POSTGRES_UNICHEM_PORT=5446
   export DATAGROK_DEMO_POSTGRES_STARBUCKS_PORT=5447
   export DATAGROK_DEMO_POSTGRES_WORLD_PORT=5448
   ```

   </TabItem>
   </Tabs>

4. The last step is to run the second stand. It is important to change the project name to start the second Datagrok
   stand. The project name in [the standard instructions](#installing-datagrok), which were used for the first stand,
   is `datagrok`. For example, you can add an increment to the project name: `datagrok_2`

   <Tabs groupId="os" queryString>
   <TabItem value="win" label="Windows" default>

   ```cmd
   docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok_2 ^
     --profile all up -d
   ```

   </TabItem>
   <TabItem value="bash" label="MacOS/Linux">

   ```shell
   docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok_2 \
     --profile all up -d
   ```

   </TabItem>
   </Tabs>

### Demo Databases

Datagrok offers demo databases that include sample data.
To install them locally, run the following command:

<Tabs groupId="os" queryString>
<TabItem value="win" label="Windows" default>

```cmd
docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok_2 ^
  --profile all --profile demo up -d
```

</TabItem>
<TabItem value="bash" label="MacOS/Linux">

```shell
docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok_2 \
  --profile all --profile demo up -d
```

</TabItem>
</Tabs>

## Shutting down Datagrok

To shut down Datagrok use this command:

<Tabs groupId="os" queryString>
<TabItem value="win" label="Windows" default>

```cmd
docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok ^
  --profile all --profile demo stop
```

</TabItem>
<TabItem value="bash" label="MacOS/Linux">

```shell
docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
  --profile all --profile demo stop
```

</TabItem>
</Tabs>

All the data is saved in the [Docker volumes](https://docs.docker.com/storage/volumes/).

To reset Datagrok to factory settings and remove all stored data, including all created users, projects, connections,
etc.,
run the following command:

<Tabs groupId="os" queryString>
<TabItem value="win" label="Windows" default>

```cmd
docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok ^
  --profile all --profile demo down --volumes
```

</TabItem>
<TabItem value="bash" label="MacOS/Linux">

```shell
docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
  --profile all --profile demo down --volumes
```

</TabItem>
</Tabs>

## Troubleshooting

1. In case of any issues, check the settings in the Datagrok (Tools -> Settings...).

    - Connectors
        - External Host: `grok_connect`
    - Scripting:
        - CVM Url: `http://cvm:8090`
        - CVM URL Client: `http://localhost:8090`
        - H2o Url: `http://h2o:54321`
        - API Url: `http://datagrok:8080/api`
        - Cvm Split: `true`
    - Dev:
        - CVM Url: `http://localhost:8090`
        - Cvm Split: `true`
        - API Url: `http://datagrok:8080/api`

2. Check containers logs for any possible errors and report the problem if there are any

   <Tabs groupId="os" queryString>
   <TabItem value="win" label="Windows" default>

   ```cmd
   docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok ^
     --profile all --profile demo logs
   ```

   </TabItem>
   <TabItem value="bash" label="MacOS/Linux">

   ```shell
   docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
     --profile all --profile demo logs
   ```

   </TabItem>
   </Tabs>

   You can also watch the logs of the desired service in real-time.
    - Replace `<service>` with the necessary service name
    - Replace `<number>` with desired log lines to watch or remove `--tail <number>` at all, if you want to see the full log

   <Tabs groupId="os" queryString>
   <TabItem value="win" label="Windows" default>

   ```cmd
   docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok ^
     --profile all --profile demo logs -f --tail <number> <service>
   ```

   </TabItem>
   <TabItem value="bash" label="MacOS/Linux">

   ```shell
   docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
     --profile all --profile demo logs -f --tail <number> <service>
   ```

   </TabItem>
   </Tabs>

3. Restart Docker compose stand

   <Tabs groupId="os" queryString>
   <TabItem value="win" label="Windows" default>

   ```cmd
   docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok ^
     --profile all --profile demo down
   docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok ^
     --profile all --profile demo up -d
   ```

   </TabItem>
   <TabItem value="bash" label="MacOS/Linux">

   ```shell
   docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
     --profile all --profile demo down
   docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
     --profile all --profile demo up -d
   ```

   </TabItem>
   </Tabs>

4. For advanced service troubleshooting, you can access the containers shell.
    - Replace `<service>` with one of the services: db, datagrok, grok_connect, grok_compute, grok_spawner,
      jupyter_notebook, jupyter_kernel_gateway, h2o, etc.

   <Tabs groupId="os" queryString>
   <TabItem value="win" label="Windows" default>

   ```cmd
   docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok ^
     --profile all --profile demo exec <service> /bin/sh
   ```

   </TabItem>
   <TabItem value="bash" label="MacOS/Linux">

   ```shell
   docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
     --profile all --profile demo exec <service> /bin/sh
   ```

   </TabItem>
   </Tabs>

5. Docker logs might take up all your free disk space. If such a situation has already taken place:

    - Limit Disk Usage
      using [Docker Desktop Resources Settings](https://docs.docker.com/desktop/settings/windows/#resources)
    - Or add to the Docker Daemon
      configuration [log properties](https://docs.docker.com/config/containers/logging/local/#usage).

6. As your last resort, recreate stand completely

   <Tabs groupId="os" queryString>
   <TabItem value="win" label="Windows" default>

   ```cmd
   docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok ^
     --profile all --profile demo down --volumes
   docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok ^
     --profile all --profile demo pull
   docker compose -f docker\localhost.docker-compose.yaml --project-name datagrok ^
     --profile all --profile demo up -d
   ```

   </TabItem>
   <TabItem value="bash" label="MacOS/Linux">

   ```shell
   docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
     --profile all --profile demo down --volumes
   docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
     --profile all --profile demo pull
   docker compose -f docker/localhost.docker-compose.yaml --project-name datagrok \
     --profile all --profile demo up -d
   ```

   </TabItem>
   </Tabs>

## Useful links

- [Basic local machine deployment](docker-compose.md)
- [All deployment options](../deploy.md)
- [Infrastructure overview](../../develop/under-the-hood/infrastructure.md)
- [Server configuration properties](../configuration.md)