---
title: "Docker Swarm"
sidebar_position: 3
---

The deployment consists of a few docker containers, [database](../../develop/under-the-hood/infrastructure.md#1-core-components) for storing metadata,
and [persistent file storage](../../develop/under-the-hood/infrastructure.md#1-core-components) for storing files

Datagrok requires PostgreSQL [database](../../develop/under-the-hood/infrastructure.md#1-core-components) to store metadata.
We recommend using scalable and highly reliable solutions, such as [AWS RDS](https://aws.amazon.com/rds/).

For [persistent file storage](../../develop/under-the-hood/infrastructure.md#1-core-components), Datagrok supports a lot of options, including cloud solutions,
for example [AWS S3](https://aws.amazon.com/s3/) and Local File System storage.

This document contains instructions to deploy Datagrok using [Docker Swarm](https://docs.docker.com/engine/swarm/)
on virtual machines or dedicated hosts with Local File System for persistent storage.

More information about Datagrok design and components:

* [Architecture](../../develop/under-the-hood/architecture.md)
* [Infrastructure](../../develop/under-the-hood/infrastructure.md)

## Prerequisites

1. We use native [Docker compose](https://docs.docker.com/compose/) commands to run applications on machines. It
   simplifies multi-container application development and deployment.
    1. Platform requires 2 hosts: 1 for Datagrok, 1 for CVM with hardware parameters:
       * Minimal: 30 GB free disk space, 2 CPUs, 4 GB RAM.
       * Recommended: Over 30 GB free disk space, 4 CPUs, 8 GB RAM, or higher.
    2. Download and install the latest version of [Docker Compose](https://docs.docker.com/compose/install/) to your
       environment (virtual machines, dedicated hosts, etc.)

## Preparations

1. Activate docker swarm mode on your manager node (it should be a datagrok host):

   ```shell
   docker swarm init
   ```

   The shell will print a token for adding worker nodes something like that:

   ```shell
   Swarm initialized: current node (bvz81updecsj6wjz393c09vti) is now a manager.

   To add a worker to this swarm, run the following command:

    docker swarm join \
    --token SWMTKN-1-testetsttestetst-qwewsdwqdwdqwddwqdwq XX.XX.XX.XX:2377

   To add a manager to this swarm, run 'docker swarm join-token manager' and follow the instructions.
   ```

2. Login on worker node (it should be a CVM host) and add it to the docker-swarm cluster:

   ```shell
       docker swarm join --token SWMTKN-1-testetsttestetst-qwewsdwqdwdqwddwqdwq XX.XX.XX.XX:2377
   ```

3. From the manager node (datagrok) set labels to all nodes in the swarm cluster.

   1. List all nodes. We will need ids of the datagrok and CVM nodes for the next step:

      ```shell
      docker node ls
      ```

   2. Add labels ```role==datagrok``` to datagrok node which is the manager and ```role==cvm```
      to CVM node which is the worker in the cluster:

      ```shell
      docker node update --label-add  role==datagrok <datagrok_node_id>
      docker node update --label-add  role==cvm <cvm_node_id>
      ```

4. Download a Docker Compose YAML
   file: [link](https://github.com/datagrok-ai/public/blob/master/docker/swarm.docker-compose.yaml)

5. Run the docker stack

   ```shell
   docker stack deploy -c ./swarm,docker-compose.yaml datagrok
   ```

6. Check if Datagrok started successfully: `http://<DATAGROK_VM_IP_ADDRESS>:8080`, login to Datagrok using the
   username "`admin`" and password "`admin`".

7. Edit settings in the running Datagrok platform (Tools -> Settings...). Do not forget to click Apply to save new settings.
    * Scripting:
        * CVM Url: `http://<CVM_VM_DNS_ADDRESS>:8090`
        * CVM URL Client: `http://<CVM_VM_DNS_ADDRESS>:8090`
        * H2o Url: `http://<CVM_VM_DNS_ADDRESS>:54321`
        * API Url: `http://<DATAGROK_VM_DNS_ADDRESS>:8080/api`
        * Cvm Split: `true`
    * Dev:
        * CVM Url: `http://<CVM_VM_DNS_ADDRESS>:8090`
        * Cvm Split: `true`
        * API Url: `http://<DATAGROK_VM_DNS_ADDRESS>:8080/api`

## Users access

Both Compute and Datagrok engines should be accessible by users.
The easiest way is to create DNS endpoints pointing to public IPs or load balancers in front of the
services: `datagrok.example`
and `cvm.example`.

## Useful links

* [Server configuration properties](../configuration.md)
