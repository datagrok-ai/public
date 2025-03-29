---
title: "Local machine"
slug: /deploy/docker-compose
sidebar_position: 0
format: 'mdx'
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

This article provides step-by-step instructions for running Datagrok on your local machine
using [Docker Compose](https://docs.docker.com/compose/). By following these instructions, you can quickly set up and
run the necessary Docker images to get started with Datagrok.

:::info Hardware requirements

Minimal hardware requirements: 40 GB of free disk space, 2 CPUs, 4 GB RAM.

:::

## Start Datagrok

1. Install and launch the latest Docker Desktop application for your operating system:

    - [Docker Desktop for Windows](https://docs.docker.com/desktop/install/windows-install/)
    - [Docker Desktop for MacOS](https://docs.docker.com/desktop/install/mac-install/)
    - [Docker Desktop for Linux](https://docs.docker.com/desktop/install/linux-install/)

2. Download the installation script based on your operating system:

   <Tabs groupId="os" queryString>
   <TabItem value="win" label="Windows" default>

   Right-click
   on [this link](https://raw.githubusercontent.com/datagrok-ai/public/master/docker/datagrok-install-local.cmd) and
   select "Save Link As..."

   </TabItem>
   <TabItem value="silicon" label="MacOS - Apple Silicon">

   Right-click
   on [this link](https://raw.githubusercontent.com/datagrok-ai/public/refs/heads/master/docker/datagrok-install-local-macos-silicon.sh) and
   select "Save Link As..."

   </TabItem>
   <TabItem value="intel" label="MacOS/Linux - Intel">

   Right-click
   on [this link](https://raw.githubusercontent.com/datagrok-ai/public/master/docker/datagrok-install-local.sh) and
   select "Save Link As..."

   </TabItem>
   </Tabs>

3. Run the downloaded script to start the installation:

   <Tabs groupId="os" queryString>
   <TabItem value="win" label="Windows" default>

   Double-click the downloaded script to begin the installation.
   
   </TabItem>
   <TabItem value="silicon" label="MacOS - Apple Silicon">

   Open the command-line interface and navigate to the directory where you saved the script. Then, run
   the following command:

   ```shell
   chmod +x datagrok-install-local.sh
   bash datagrok-install-local-macos-silicon.sh
   ```

   </TabItem>
   <TabItem value="intel" label="MacOS/Linux">

   Open the command-line interface and navigate to the directory where you saved the script. Then, run
   the following command:

   ```shell
   chmod +x datagrok-install-local.sh
   bash datagrok-install-local.sh
   ```

   </TabItem>
   </Tabs>

## Log in to Datagrok

1. Once the server is up and running, open your browser and go to [http://localhost:8080](http://localhost:8080).
2. On the login page, use the following credentials to login:
   - Login or Email: `admin`
   - Password `admin`

:::note

If you see the message `Datagrok server is unavaliable`, wait for approximately 1 minute for the server to start, and then reload the page.

:::

## Shut down Datagrok

Use the following command from the directory where you downloaded the installation script:

<Tabs groupId="os" queryString>
<TabItem value="win" label="Windows" default>

```cmd
datagrok-install-local.cmd stop
```

</TabItem>
<TabItem value="bash" label="MacOS/Linux">

```shell
bash datagrok-install-local.sh stop
```

</TabItem>
</Tabs>

All the data is saved in the [Docker volumes](https://docs.docker.com/storage/volumes/).

## Useful links

- [Advanced local machine deployment](docker-compose-advanced.md)
- [All deployment options](../deploy.md)
- [Infrastructure overview](../../develop/under-the-hood/infrastructure.md)
- [Server configuration properties](../configuration.md)