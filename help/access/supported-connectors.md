---
title: "Supported connectors"
---

A connector can work with a database, an Excel file, a CSV file, a web service,
or basically anything that is capable of providing the data. Datagrok currently
supports over 30 different connectors, and the list is quickly growing. Most of
Datagrok’s data connectors are open-sourced and could be found on
[GitHub](https://github.com/datagrok-ai/public/tree/master/connectors) (MIT
license). The supported connectors with their specific parameters are the
following:

| Data Source                                         | Server | Port | DB  |Browse Schema | Cache Schema | Cache Results | SSL | Connection String | Login | Password | Other Parameters |
|-----------------------------------------------------|--------|------|-----|--------------|---------------|-----|-------------------|-------|----------|---------------------------------------------------------------------------|-----|
| [Access](../access/connectors/access.md)            |        |      | ✓   |     |              |               |     | ✓                 | ✓     | ✓        |                                                                           |
| [Athena](../access/connectors/athena.md)            | ✓      | ✓    | ✓   |     |              |               |     | ✓                 |       |          | [See the list](../access/connectors/athena.md)                            |
| [BigQuery](../access/connectors/bigquery.md)        |        |      |     |      |             |               |     | ✓                 | ✓     | ✓        | [See the list](../access/connectors/bigquery.md#connection-parameters)    |
| [Cassandra](../access/connectors/cassandra.md)      | ✓      | ✓    | ✓   |     | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [DB2](../access/connectors/db2.md)                  | ✓      | ✓    | ✓   |     | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [Denodo](../access/connectors/denodo.md)            | ✓      | ✓    | ✓   |     | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [DropBox](../access/connectors/dropbox.md)          |        |      |     |     |              |               |     |                   |       | ✓        | [See the list](../access/connectors/dropbox.md#connection-parameters)     |
| [Files](../access/connectors/files.md)              |        |      |     |     |              |               |     |                   | ✓     | ✓        | [See the list](../access/connectors/files.md#connection-parameters)       |
| [Firebird](../access/connectors/firebird.md)        | ✓      | ✓    | ✓   |     | ✓            | ✓             |     | ✓                 | ✓     | ✓        |                                                                           |
| [Git](../access/connectors/git.md)                  |        |      |     |     |              |               |     |                   |       |          | [See the list](../access/connectors/git.md#connection-parameters)         |
| [Google Cloud](../access/connectors/googlecloud.md) |        |      |     |     |              |               |     |                   |       |          | [See the list](../access/connectors/googlecloud.md#connection-parameters) |
| [HBase](../access/connectors/hbase.md)              | ✓      | ✓    | ✓   |     | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [Hive](../access/connectors/hive.md)                | ✓      | ✓    | ✓   |     | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [Hive2](../access/connectors/hive2.md)              | ✓      | ✓    | ✓   |     | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [Impala](../access/connectors/impala.md)            | ✓      | ✓    | ✓   |     |              |               |     | ✓                 | ✓     | ✓        | [See the list](../access/connectors/impala.md#connection-parameters)      |
| [MariaDB](../access/connectors/mariadb.md)          | ✓      | ✓    | ✓   |   ✓     |✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [MongoDB](../access/connectors/mongodb.md)          | ✓      | ✓    | ✓   |     | ✓            | ✓             |     | ✓                 | ✓     | ✓        |                                                                           |
| [MS SQL](../access/connectors/mssql.md)             | ✓      | ✓    | ✓   |  ✓    | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [MySql](../access/connectors/mysql.md)              | ✓      | ✓    | ✓   |  ✓    | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [Neo4j](../access/connectors/neo4j.md)              | ✓      | ✓    | ✓   |     | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [OData](../access/connectors/odata.md)              |        |      |     |     |              |               |     |                   |       |          | [See the list](../access/connectors/odata.md#connection-parameters)       |
| [Oracle](../access/connectors/oracle.md)            | ✓      | ✓    | ✓   |   ✓    |✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [Postgres](../access/connectors/postgres.md)        | ✓      | ✓    | ✓   |  ✓     | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [PostgresDart](../access/connectors/postgres.md)      | ✓      |      | ✓   |  ✓    |              | ✓             | ✓   |                   | ✓     | ✓        |                                                                           |
| [Redshift](../access/connectors/redshift.md)        | ✓      | ✓    | ✓   |     | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [S3](../access/connectors/s3.md)                    |        |      |     |     |              |               |     |                   |       |          | [See the list](../access/connectors/s3.md#connection-parameters)          |
| [Snowflake](../access/connectors/snowflake.md)      | ✓      | ✓    | ✓   |     | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [Socrata](../access/connectors/socrata.md)          |        |      |     |     |              |               |     |                   |       |          | [See the list](../access/connectors/socrata.md#connection-parameters)     |
| [Sparql](../access/connectors/sparql.md)            |        |      |     |     |              |               |     |                   |       |          | [See the list](../access/connectors/sparql.md#connection-parameters)      |
| [SQLite](../access/connectors/sqlite.md)            |        |      | ✓   |     |              |               |     | ✓                 | ✓     | ✓        |                                                                           |
| [Teradata](../access/connectors/teradata.md)        | ✓      | ✓    | ✓   |     | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [Twitter](../access/connectors/twitter.md)          |        |      |     |     |              |               |     |                   |       |          | [See the list](../access/connectors/twitter.md#connection-parameters)     |
| [Vertica](../access/connectors/vertica.md)          | ✓      | ✓    | ✓   |      |✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [Virtuoso](../access/connectors/virtuoso.md)        | ✓      | ✓    | ✓   |     | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [Web](../access/connectors/web.md)                  |        |      |     |     |              |               |     |                   |       |          | [See the list](../access/connectors/web.md#connection-parameters)         |

See also:

* [Databases](databases.md)
* [File shares](file-shares.md)
