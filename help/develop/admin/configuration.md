<!-- TITLE: Datagrok Server Configuration -->
<!-- SUBTITLE: -->

# Configuration

Datagrok supports several deployment schemas which can be configured using `GROK_MODE` and `GROK_PARAMETERS` environment
variables.

## Datlas configuration

`GROK_PARAMETERS` is a JSON-formatted environment variable with the Datlas configuration options :

| Option              | Required        | Default | Description                                                                                      |
|---------------------|-----------------|---------|--------------------------------------------------------------------------------------------------|
| dbServer            | <b>Required</b> |         | Postgres database server.                                                                        |                                                                                                  |
| dbPort              | Optional        | 5432    | Postgres database server port                                                                    |
| db                  | <b>Required</b> |         | Datagrok database name                                                                           |
| dbLogin             | <b>Required</b> |         | Username to connect to database                                                                  |
| dbPassword          | <b>Required</b> |         | Password to connect to database                                                                  |
| dbSsl               | Optional        | false   | If set to true, TLS connection will be used to connect to database                               |
| dbAdminLogin        | Optional        |         | Postgres admin username to create user and database for Datagrok                                 |
| dbAdminPassword     | Optional        |         | Postgres admin password to create user and database for Datagrok                                 |
| googleStorageCert   | Optional        |         | Access certificate to Google Cloud Storage. If set, GCS will be used for persistent data storage |
| amazonStorageRegion | Optional        |         | S3 region                                                                                        |
| amazonStorageBucket | Optional        |         | S3 bucket name                                                                                   |
| amazonStorageId     | Optional        |         | S3 credential ID, Datagrok will resolve EC2 role if empty                                        |
| amazonStorageKey    | Optional        |         | S3 credential secret key, Datagrok will resolve EC2 role if empty                                |
| adminPassword       | Optional        |         | Datagrok admin user password which will be created on first start                                |
| debug               | Optional        | false   | Extended logging and saving stack traces                                                         |
| useSSL              | Optional        | false   | If set to true, Datlas serves TLS connections                                                    |
| certPath            | Optional        |         | Path to the TLS certificate                                                                      |
| certKeyPath         | Optional        |         | Path to the TLS certificate key                                                                  |
| certKeyPwd          | Optional        |         | Password to the TLS certificate key                                                              |

## Datlas startup mode

`GROK_MODE` possible values:

* `start` - starts the application without database and storage deployment
* `deploy` - Datlas will perform the full deploy
* `auto` - Datlas will check the existing database and storage and perform deployment only if needed
