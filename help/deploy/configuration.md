---
title: "Server configuration"
sidebar_position: 13
description: Reference for GROK_MODE and GROK_PARAMETERS environment variables that configure the Datlas server.
keywords:
  - grok_parameters
  - grok_mode
  - environment variable overrides
  - queue settings
  - docker settings
  - connectors settings
---

Datagrok supports several deployment schemas which can be configured using `GROK_MODE` and `GROK_PARAMETERS` environment variables.

## Datlas Configuration

`GROK_PARAMETERS` is a JSON-formatted environment variable with the Datlas configuration options:

| Option                   | Required     | Default                                         | Description                                                                                                                                                 |
|--------------------------|--------------|-------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------|
| dbServer                 | **Required** |                                                 | Postgres database server                                                                                                                                    |
| dbPort                   | Optional     | 5432                                            | Postgres database server port                                                                                                                               |
| db                       | **Required** |                                                 | Datagrok database name                                                                                                                                      |
| dbLogin                  | **Required** |                                                 | Database user Datagrok uses for everything: runtime, schema migrations, and bootstrap. The user must own the database and have `CREATEROLE` (plus `CREATEDB` for a fresh install) |
| dbPassword               | **Required** |                                                 | Password of `dbLogin` |
| dbSsl                    | Optional     | false                                           | If set to true, TLS connection will be used to connect to database                                                                                          |
| dbAdminLogin             | Optional     |                                                 | Deprecated. If set to a different user, Datagrok uses this user instead of `dbLogin` as its only database credential |
| dbAdminPassword          | Optional     |                                                 | Deprecated. Password of `dbAdminLogin` |
| recreateDb               | Optional     | false                                           | Drops and recreates the database on a full deployment. Destructive. Without it, Datagrok creates a missing database and never drops an existing one |
| googleStorageCert        | Optional     |                                                 | Access certificate to Google Cloud Storage. If set, GCS will be used for persistent data storage                                                            |
| amazonStorageRegion      | Optional     |                                                 | S3 region                                                                                                                                                   |
| amazonStorageBucket      | Optional     |                                                 | S3 bucket name                                                                                                                                              |
| amazonStorageId          | Optional     |                                                 | S3 credential ID, Datagrok will resolve EC2 role if empty                                                                                                   |
| amazonStorageKey         | Optional     |                                                 | S3 credential secret key, Datagrok will resolve EC2 role if empty                                                                                           |
| amazonStorageEndpoint    | Optional     |                                                 | Custom S3-compatible endpoint host (e.g. eu2.contabostorage.com); uses path-style addressing when set                                                       |
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
| connectorsSettings       | Optional     | See [Connectors Settings](#connectors-settings) | List of Grok Connect endpoints for [External Database Connectivity](../develop/under-the-hood/infrastructure.md#3-external-database-connectivity) |

Datagrok also creates a read-only `<db>_reader` database role and gives it a new random
password on every start. You don't configure it. The former `dbReaderLogin` and
`useAdminForMigrations` options are no longer used.

### Queue Settings

`queueSettings` object:

| Option                 | Default      | Description                                                         |
|------------------------|--------------|---------------------------------------------------------------------|
| useQueue               | true         | Enables usage of queue for calls execution (required for scripting) |
| amqpHost               | rabbitmq     | Host of the AMQP message queue                                      |
| amqpPort               | 5672         | Port of the AMQP message queue                                      |
| amqpAuthMode           | password     | AMQP authentication: `password` (static user) or `jwt` (Datagrok-issued tokens) |
| mqAudience             | datagrok-mq  | Audience of broker tokens. Must match the RabbitMQ `resource_server_id` |
| amqpUser               | guest        | AMQP username in `password` mode                                    |
| amqpPassword           | guest        | AMQP password in `password` mode                                    |
| tls                    | false        | Enables TLS for AMQP connection                                     |
| taskPickupTimeout      | 60000        | Timeout in ms for task pickup by worker                             |
| queueReconnectWaitTime | 1500         | Time in ms to wait before retrying connection                       |
| pipeHost               | grok_pipe    | Host of the grok_pipe service                                       |
| pipePort               | 3000         | Port of the grok_pipe service                                       |
| pipeAuthMode           | password     | grok_pipe authentication: `password` (static `pipeKey`) or `jwt`    |
| pipeAudience           | datagrok-pipe | Audience of grok_pipe tokens. Must match the grok_pipe token audience |
| pipeKey                | datagrok-key | Legacy static key for grok_pipe, used in `password` mode           |

In `jwt` mode, Datagrok signs short-lived tokens with its [server keys](../govern/access-control/server-keys.md),
and RabbitMQ, grok_pipe, and Grok Spawner verify them against the Datagrok JWKS endpoint
(`/api/.well-known/jwks.json`). No shared password or key is stored. Static keys still work but
are legacy. The Helm chart enables `jwt` for grok_pipe and Grok Spawner by default, and for
RabbitMQ when `ingress.host` is set.

### Docker Settings

`dockerSettings` object:

| Option                   | Default        | Description                                                 |
|--------------------------|----------------|-------------------------------------------------------------|
| useGrokSpawner           | true           | Enables Grok Spawner for managing containers                |
| spawnerAuthMode          | jwt            | Grok Spawner authentication: `jwt` (Datagrok-issued tokens, see [Queue Settings](#queue-settings)) or `password` (static `grokSpawnerApiKey`) |
| spawnerAudience          | datagrok-spawner | Audience of Grok Spawner tokens. Must match the spawner token audience |
| grokSpawnerApiKey        | test-x-api-key | Legacy static API key for Grok Spawner, used in `password` mode |
| grokSpawnerHost          | grok_spawner   | Host for Grok Spawner                                       |
| grokSpawnerPort          | 8000           | Port for Grok Spawner                                       |
| imageBuildTimeoutMinutes | 30             | Max wait time in minutes for Docker image build             |
| proxyRequestTimeout      | 60000          | Max wait time in ms for proxy request to a Docker container |

### Connectors Settings

`connectorsSettings` is a list of Grok Connect endpoints. The first entry is the default
endpoint. Add the ADBC connector and the optional extended connectors as extra entries. If
the list is empty or missing, Datagrok runs without Grok Connect. A single object, the form
used before 1.28, is still accepted and becomes a one-element list.

```json
"connectorsSettings": [
  {"grokConnectHost": "grok_connect", "grokConnectPort": 1234},
  {"grokConnectHost": "grok_connect_adbc", "grokConnectPort": 1235, "externalDataFrameCompress": false},
  {"grokConnectHost": "grok_connect_extended", "grokConnectPort": 1234}
]
```

Each entry accepts these options:

| Option                    | Default      | Description                                                     |
|---------------------------|--------------|-----------------------------------------------------------------|
| grokConnectHost           | grok_connect | Host of the endpoint                                            |
| grokConnectPort           | 1234         | Port of the endpoint                                            |
| useGrokConnect            | true         | Polls this endpoint. Set to `false` to skip it                  |
| externalDataFrameCompress | true         | Compresses data frames sent to the client for queries served by this endpoint |
| gzipLevel                 | 1            | Compression level (1–9) for these data frames                   |
| dataFrameBatchSize        | 8000000      | Maximum WebSocket frame size, in bytes, for query results       |

File storage and parsing options apply to the whole server and go to the root of `GROK_PARAMETERS`, not
into a `connectorsSettings` entry:

| Option                    | Default      | Description                                                     |
|---------------------------|--------------|-----------------------------------------------------------------|
| sambaVersion              | 3.0          | Samba version                                                   |
| sambaSpaceEscape          | none         | Specifies how spaces are escaped in Samba (none, space, quotes) |
| dataframeParsingMode      | New Process  | DataFrame parsing mode (Main Thread, New Thread, New Process)   |
| localFileSystemAccess     | false        | Enables local file system access                                |
| windowsSharesProxy        |              | Proxy for Windows shares                                        |


## Overriding Datlas Configuration with Environment Variables

In addition to supplying the full JSON, you can override individual values using environment variables. This is useful when deploying in containerized or cloud environments where injecting single parameters is easier than rebuilding the entire configuration object.

### Naming Convention

The environment variable name is derived from the configuration key by applying the following rules:

1. Prefix every variable with `GROK_PARAMETERS_`
2. Flatten nested objects by joining keys with a double underscore `__`.
   Example:
    - `queueSettings.queueReconnectWaitTime` → `QUEUE_SETTINGS__QUEUE_RECONNECT_WAIT_TIME`
3. Convert key names to uppercase.
- `dbLogin` → `DB_LOGIN`
- `grokConnectHost` → `GROK_CONNECT_HOST`

### Examples

| Config Option                          | Environment Variable Name                                   |
|----------------------------------------|-------------------------------------------------------------|
| `dbLogin`                              | `GROK_PARAMETERS_DB_LOGIN`                                  |
| `dbPassword`                           | `GROK_PARAMETERS_DB_PASSWORD`                               |
| `queueSettings.queueReconnectWaitTime` | `GROK_PARAMETERS_QUEUE_SETTINGS__QUEUE_RECONNECT_WAIT_TIME` |
| `connectorsSettings[0].grokConnectHost` | `GROK_PARAMETERS_CONNECTORS_SETTINGS__GROK_CONNECT_HOST`   |
| `connectorsSettings[1].grokConnectPort` | `GROK_PARAMETERS_CONNECTORS_SETTINGS__1__GROK_CONNECT_PORT` |
| `dockerSettings.proxyRequestTimeout`   | `GROK_PARAMETERS_DOCKER_SETTINGS__PROXY_REQUEST_TIMEOUT`    |
| `sambaSpaceEscape`                     | `GROK_PARAMETERS_SAMBA_SPACE_ESCAPE`                        |

To address a `connectorsSettings` entry other than the first, add its zero-based index after
`CONNECTORS_SETTINGS__`. Without an index, the override applies to the first entry. The legacy
`GROK_PARAMETERS_CONNECTORS_SETTINGS__<option>` form of these server-wide options, such as
`CONNECTORS_SETTINGS__SAMBA_SPACE_ESCAPE`, still works.


## Datlas Startup Mode

`GROK_MODE` possible values:

* `start` - Starts the application without database and storage deployment.
* `deploy` - Datlas will perform the full deployment.
* `auto` - Datlas will check the existing database and storage and perform deployment only if needed.

## Settings 

You can set Datargok server settings as a part of `GROK_PARAMETERS`. Get the template JSON in `/settings` view using `{}` button near the server settings section.
Put any parameter to `settings` map of `GROK_PARAMETERS` respecting the hierarchy. 

## Useful links

* [Infrastructure](../develop/under-the-hood/infrastructure.md)
* [Architecture](../develop/under-the-hood/architecture.md)
