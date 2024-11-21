---
title: "Connecting to database inside package Docker container"
---

This document explains how to connect to database inside package Docker container.

## Creating connection

A connection to a database inside a package Docker container can only be created programmatically (see [Access section](access-data.md#creating-a-connection)). 
The JSON file with connection parameters should include the real credentials of the database deployed inside the container. Additionally, 
if you want to use a particular schema (or database), you should add the `db` property to the JSON. The port parameter should not be added, 
as Datagrok will resolve it dynamically based on the container address.

To ensure the Datagrok platform recognizes the connection as one that requires Docker, the `server` property of the JSON 
should have the following format: `${<package name>:<container friendly name><DockerContainer>}`. The `friendly name` of 
the container is the same as the package name when you have only one container declared in the package. If there are 
several folders inside the `dockerfiles` of the package with Dockerfile inside, then friendly name follows this format: `<package name>-<folder name>`.

For example, consider the following Dockerfile created inside a package called `Test`:

```shell
FROM postgres
ENV POSTGRES_PASSWORD datagrok
ENV POSTGRES_USER datagrok
ENV POSTGRES_DB world
EXPOSE 5432
```

This Docker container runs a [Postgres](https://www.postgresql.org/) database, exposing its port. It uses environment variables to specify the user, password and database name.

> Note: EXPOSE instruction is crucial and should always be declared.

The connection JSON file could look like this:

```json
{
  "#type": "DataConnection",
  "name": "PostgresDocker",
  "friendlyName": "PostgresDocker",
  "parameters": {
    "server": "${Test:Test<DockerContainer>}",
    "db": "world"
  },
  "credentials" : {
    "parameters" : {
      "login": "datagrok",
      "password": "datagrok"
    }
  },
  "dataSource": "Postgres"
}
```

> Note: Since there is only one Dockerfile inside this package, we use the package name as the container name.

After publishing package you will see created connection under the **Databases** in the **Browse**. Such connections can be used like any other; you can execute database queries, browse schemas, and use them within projects.

For a complete working example, please refer to the source code of our [DbTests package](https://github.com/datagrok-ai/public/tree/master/packages/DBTests).

See also:

- [Packages Docker containers](docker_containers.md)
- [Packages](../develop.md#packages)
- [Data Access](./access-data.md)
