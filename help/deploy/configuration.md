---
title: "Server configuration"
sidebar_position: 13
---

Datagrok supports several deployment schemas which can be configured using `GROK_MODE` and `GROK_PARAMETERS` environment variables.

## Datlas Configuration

`GROK_PARAMETERS` is a JSON-formatted environment variable with the Datlas configuration options:

| Option                   | Required     | Default                                         | Description                                                                                                                                                 |
|--------------------------|--------------|-------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------|
| dbServer                 | **Required** |                                                 | Postgres database server                                                                                                                                    |
| dbPort                   | Optional     | 5432                                            | Postgres database server port                                                                                                                               |
| db                       | **Required** |                                                 | Datagrok database name                                                                                                                                      |
| dbLogin                  | **Required** |                                                 | Username to connect to database                                                                                                                             |
| dbPassword               | **Required** |                                                 | Password to connect to database                                                                                                                             |
| dbSsl                    | Optional     | false                                           | If set to true, TLS connection will be used to connect to database                                                                                          |
| dbAdminLogin             | Optional     |                                                 | Postgres admin username to create user and database for Datagrok                                                                                            |
| dbAdminPassword          | Optional     |                                                 | Postgres admin password to create user and database for Datagrok                                                                                            |
| googleStorageCert        | Optional     |                                                 | Access certificate to Google Cloud Storage. If set, GCS will be used for persistent data storage                                                            |
| amazonStorageRegion      | Optional     |                                                 | S3 region                                                                                                                                                   |
| amazonStorageBucket      | Optional     |                                                 | S3 bucket name                                                                                                                                              |
| amazonStorageId          | Optional     |                                                 | S3 credential ID, Datagrok will resolve EC2 role if empty                                                                                                   |
| amazonStorageKey         | Optional     |                                                 | S3 credential secret key, Datagrok will resolve EC2 role if empty                                                                                           |
| googleStorageBucket      | Optional     |                                                 | Google Cloud Storage bucket name                                                                                                                            |
| googleStorageCredentials | Optional     |                                                 | Google Cloud Storage credentials                                                                                                                            |
| googleStorageProject     | Optional     |                                                 | Google Cloud Storage project ID                                                                                                                             |
| adminPassword            | Optional     |                                                 | Datagrok admin user password which will be created on first start                                                                                           |
| debug                    | Optional     | false                                           | Extended logging and saving stack traces                                                                                                                    |
| useSSL                   | Optional     | false                                           | If set to true, Datlas serves TLS connections                                                                                                               |
| certPath                 | Optional     |                                                 | Path to the TLS certificate                                                                                                                                 |
| certKeyPath              | Optional     |                                                 | Path to the TLS certificate key                                                                                                                             |
| certKeyPwd               | Optional     |                                                 | Password to the TLS certificate key                                                                                                                         |
| queueSettings            | Optional     | See [Queue Settings](#queue-settings)           | Configuration for [Scripting and Computations](../develop/under-the-hood/infrastructure.md#4-scripting-and-computation)                                     |
| dockerSettings           | Optional     | See [Docker Settings](#docker-settings)         | Configuration for [Plugins Docker Management](../develop/under-the-hood/infrastructure.md#5-plugin--docker-container-management)                            |
| connectorsSettings       | Optional     | See [Connectors Settings](#connectors-settings) | Configuration for [External Database Connectivity](../develop/under-the-hood/infrastructure.md#3-external-database-connectivity) and file storage providers |

### Queue Settings

`queueSettings` object:

| Option                 | Default      | Description                                                         |
|------------------------|--------------|---------------------------------------------------------------------|
| useQueue               | true         | Enables usage of queue for calls execution (required for scripting) |
| amqpHost               | rabbitmq     | Host of the AMQP message queue                                      |
| amqpPort               | 5672         | Port of the AMQP message queue                                      |
| amqpUser               | guest        | AMQP username                                                       |
| amqpPassword           | guest        | AMQP password                                                       |
| tls                    | false        | Enables TLS for AMQP connection                                     |
| taskPickupTimeout      | 60000        | Timeout in ms for task pickup by worker                             |
| queueReconnectWaitTime | 1500         | Time in ms to wait before retrying connection                       |
| pipeHost               | grok_pipe    | Host of the grok_pipe service                                       |
| pipePort               | 3000         | Port of the grok_pipe service                                       |
| pipeKey                | datagrok-key | Key for grok_pipe authentication                                    |

### Docker Settings

`dockerSettings` object:

| Option                   | Default        | Description                                                 |
|--------------------------|----------------|-------------------------------------------------------------|
| useGrokSpawner           | true           | Enables Grok Spawner for managing containers                |
| grokSpawnerApiKey        | test-x-api-key | API key for Grok Spawner                                    |
| grokSpawnerHost          | grok_spawner   | Host for Grok Spawner                                       |
| grokSpawnerPort          | 8000           | Port for Grok Spawner                                       |
| imageBuildTimeoutMinutes | 30             | Max wait time in minutes for Docker image build             |
| proxyRequestTimeout      | 60000          | Max wait time in ms for proxy request to a Docker container |

### Connectors Settings

`connectorsSettings` object:

| Option                    | Default      | Description                                                     |
|---------------------------|--------------|-----------------------------------------------------------------|
| useGrokConnect            | true         | Enables Grok Connect for additional data connectors             |
| grokConnectHost           | grok_connect | Host for Grok Connect                                           |
| grokConnectPort           | 1234         | Port for Grok Connect                                           |
| externalDataFrameCompress | true         | Enables external data frame compression                         |
| sambaVersion              | 3.0          | Samba version                                                   |
| sambaSpaceEscape          | none         | Specifies how spaces are escaped in Samba (none, space, quotes) |
| dataframeParsingMode      | New Process  | DataFrame parsing mode (Main Thread, New Thread, New Process)   |
| localFileSystemAccess     | false        | Enables local file system access                                |
| windowsSharesProxy        |              | Proxy for Windows shares                                        |

## Datlas Startup Mode

`GROK_MODE` possible values:

* `start` - Starts the application without database and storage deployment.
* `deploy` - Datlas will perform the full deployment.
* `auto` - Datlas will check the existing database and storage and perform deployment only if needed.

## Useful links

* [Infrastructure](../develop/under-the-hood/infrastructure.md)
* [Architecture](../develop/under-the-hood/architecture.md)
