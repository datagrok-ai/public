---
title: "Deployment on Docker Swarm"
---

Datagrok is based on Docker containers, [database](../infrastructure.md#database)
and [persistent file storage](../infrastructure.md#storage).

As [database](../infrastructure.md#database) Datagrok supports any PostgreSQL database out-of-the-box, including cloud
solutions for PostgreSQL database, for example [AWS RDS](https://aws.amazon.com/rds/).

For [persistent file storage](../infrastructure.md#storage), Datagrok supports a lot of options, including cloud solutions,
for example [AWS S3](https://aws.amazon.com/s3/) and Local File System storage.

This document contains instructions to deploy Datagrok using [Docker Swarm](https://docs.docker.com/engine/swarm/)
on virtual machines or dedicated hosts with Local File System for persistent storage.

More information about Datagrok design and components:

* [Architecture](../architecture.md)
* [Infrastructure](../infrastructure.md)

## Prerequisites

1. We use native [Docker compose](https://docs.docker.com/compose/) commands to run applications on machines. It
   simplifies multi-container application development and deployment.
    1. Download and install the latest version of [Docker Compose](https://docs.docker.com/compose/install/) to your
       environment (virtual machines, dedicated hosts, etc.)

## Preparations

1. Activate docker swarm mode on your manager node (it should be a datagrok host):

   ```shell
   docker swarm init
   ```

   The shell will print a token for adding worker nodes something like that:

   ```shell
   docker swarm join \
    --token SWMTKN-1-49nj1cmql0jkz5s954yi3oex3nedyz0fb0xx14ie39trti4wxv-8vxv8rssmk743ojnwacrr2e7c \
    192.168.99.100:2377
   ```

2. Use the received token on worker nodes to join the swarm as described below. If you lost a token, you can get
   it by typing in the manager node shell for workers:

   ```shell
   docker swarm join-token worker
   ```

   Or, if you want to add another manager node:

   ```shell
   docker swarm join-token manager
   ```

3. From the manager node set labels to all nodes in the swarm cluster.

   1. List all nodes:

      ```shell
      docker node ls
      ```

   2. Add labels ```role==datagrok``` and ```role==cvm``` to nodes:

      ```shell
      docker node update --label-add  role==datagrok <prefer_node>
      ```

4. Download a Docker Compose YAML
   file: [link](https://github.com/datagrok-ai/public/blob/master/docker/swarm.docker-compose.yaml)

5. Adjust in the downloaded file `GROK_PARAMETERS` values,
`X_API_KEY`,`DOCKER_REGISTRY`,`DOCKER_REGISTRY_USER`,
`DOCKER_REGISTRY_PASS`;`POSTGRES_USER` and `POSTGRES_PASSWORD`must be the same as
`dbAdminLogin` and `dbAdminPassword` in `GROK_PARAMETERS`.

6. Run the docker stack

   ```shell
   docker stack deploy -c ./swarm,docker-compose.yaml <STACK_NAME>
   ```

7. Check if Datagrok started successfully: `http://<DATAGROK_VM_IP_ADDRESS>:8080`, login to Datagrok using the
   username "`admin`" and password "`admin`".

8. Edit settings in the running Datagrok platform (Tools -> Settings...). Do not forget to click Apply to save new settings.
    * Scripting:
        * CVM Url: `http://<CVM_VM_IP_ADDRESS>:8090`
        * CVM URL Client: `http://<CVM_VM_IP_ADDRESS>:8090`
        * H2o Url: `http://<CVM_VM_IP_ADDRESS>:54321`
        * API Url: `http://<DATAGROK_VM_IP_ADDRESS>:8080/api`
        * Cvm Split: `true`
    * Dev:
        * CVM Url: `http://<CVM_VM_IP_ADDRESS>:8090`
        * Cvm Split: `true`
        * API Url: `http://<DATAGROK_VM_IP_ADDRESS>:8080/api`

## Users access

Both Compute and Datagrok engines should be accessible by users.
The easiest way is to create DNS endpoints pointing to public IPs or load balancers in front of the
services: `datagrok.example`
and `cvm.example`.
