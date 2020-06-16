
<!-- TITLE: Deploy Datagrok using Docker Compose -->
<!-- SUBTITLE: -->

# Deploy Datagrok using Docker Compose

This document contains instructions to deploy Datagrok on a regular machine without cloud-based hosting.

1. Download latest [docker-compose.yaml](../docker-compose.yaml)
2. Put to new empty folder
3. Run `docker-compose up`

Datagrok will start and automatically deploy a new database. After it deployed first time you can shut it down using `Ctrl+C`.
All data will be saved in persistent storage.

See also:

* [Docker Compose](https://docs.docker.com/compose/)