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

| Data Source                                         | Type | Server | Port | DB  | Browse Schema | Cache Schema | Cache Results | SSL | Connection String | Login | Password | Other Parameters |
|-----------------------------------------------------|--------|--------|------|-----|--------------|---------------|-----|-------------------|-------|----------|---------------------------------------------------------------------------|-----|
| [Access](../access/connectors/access.md) | JDBC | | | ✓ | | | | | ✓  | ✓  | ✓  | |
| [Athena](../access/connectors/athena.md) | JDBC | ✓ | ✓ | ✓ | | | | | ✓ | | | [See the list](../access/connectors/athena.md) |
| [BigQuery](../access/connectors/bigquery.md) | JDBC | | | | | | | | ✓ | ✓ | ✓ | [See the list](../access/connectors/bigquery.md#connection-parameters)    |
| [Cassandra](../access/connectors/cassandra.md) | JDBC | ✓ | ✓ | ✓ | | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | |
| [ClickHouse](../access/connectors/clickhouse.md) | JDBC | ✓ | ✓ | ✓ | | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | |
|[CoreWeave](../access/connectors/coreweave.md) | Files | | | | | | | | | | | [See the list](../access/connectors/coreweave.md#connection-parameters) |
| [DB2](../access/connectors/db2.md) | JDBC | ✓ | ✓ | ✓ | | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | |
| [Denodo](../access/connectors/denodo.md) | JDBC  | ✓ | ✓ | ✓ | | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | |
| [DropBox](../access/connectors/dropbox.md) | Files | | | | | | | | | | ✓ | [See the list](../access/connectors/dropbox.md#connection-parameters) |
| [Files](../access/connectors/files.md) | Files | | | | | | | | | ✓ | ✓ | [See the list](../access/connectors/files.md#connection-parameters) |
| [Firebird](../access/connectors/firebird.md) | JDBC  | ✓ | ✓ | ✓ | | ✓ | ✓ | | ✓ | ✓ | ✓ | |
| [Git](../access/connectors/git.md) | Files | | | | | | | | | | | [See the list](../access/connectors/git.md#connection-parameters) |
| [Google Cloud](../access/connectors/googlecloud.md) | Files | | | | | | | | | | | [See the list](../access/connectors/googlecloud.md#connection-parameters) |
| [HBase](../access/connectors/hbase.md) | JDBC | ✓ | ✓ | ✓ | | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | |
| [Hive](../access/connectors/hive.md) | JDBC | ✓ | ✓ | ✓ | | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | |
| [Hive2](../access/connectors/hive2.md) | JDBC  | ✓ | ✓ | ✓ | | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | |
| [Impala](../access/connectors/impala.md) | JDBC | ✓ | ✓ | ✓ | | | | | ✓ | ✓ | ✓ | [See the list](../access/connectors/impala.md#connection-parameters) |
| [MariaDB](../access/connectors/mariadb.md) | JDBC | ✓ | ✓ | ✓ | ✓ |✓ | ✓ | ✓ | ✓ | ✓ | ✓ | |
| [MongoDB](../access/connectors/mongodb.md) | JDBC  | ✓ | ✓ | ✓ | | ✓ | ✓ | | ✓ | ✓ | ✓ ||
| [MS SQL](../access/connectors/mssql.md)| JDBC  | ✓      | ✓    | ✓   |  ✓    | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        ||
| [MySql](../access/connectors/mysql.md)              | JDBC  | ✓      | ✓    | ✓   |  ✓    | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [Neo4j](../access/connectors/neo4j.md)              | JDBC  | ✓      | ✓    | ✓   |     | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [Neptune](../access/connectors/neptune.md)      | JDBC  | ✓      | ✓    |    |     | ✓            | ✓             |    | ✓                 |      |         |                                       [See the list](../access/connectors/neptune.md#connection-parameters) |
| [OData](../access/connectors/odata.md)              |       |        |      |     |     |              |               |     |                   |       |          | [See the list](../access/connectors/odata.md#connection-parameters)       |
| [Oracle](../access/connectors/oracle.md)            | JDBC  | ✓      | ✓    | ✓   |   ✓    |✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [PI](../access/connectors/pi.md)                   | JDBC  | ✓       |      | ✓    | ✓    | ✓             | ✓               |     |  ✓                 |    ✓   | ✓         |    [See the list](../access/connectors/pi.md#connection-parameters)       |
| [Postgres](../access/connectors/postgres.md)        | JDBC  | ✓      | ✓    | ✓   |  ✓     | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [PostgresDart](../access/connectors/postgres.md)      | JDBC  | ✓      |      | ✓   |  ✓    |              | ✓             | ✓   |                   | ✓     | ✓        |                                                                           |
| [Redshift](../access/connectors/redshift.md)        | JDBC  | ✓      | ✓    | ✓   | ✓    | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [S3](../access/connectors/s3.md)                    | Files  |        |      |     |     |              |               |     |                   |       |          | [See the list](../access/connectors/s3.md#connection-parameters)          |
| [SAP HANA](../access/connectors/sap-hana.md)                    | JDBC  |    ✓    |   ✓   |  ✓   | ✓    |              |               |     |                   |   ✓    |  ✓        | [See the list](../access/connectors/sap-hana.md#connection-parameters)          |
| [Snowflake](../access/connectors/snowflake.md)      | JDBC  | ✓      | ✓    | ✓   | ✓    | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [Socrata](../access/connectors/socrata.md)          |       |        |      |     |     |              |               |     |                   |       |          | [See the list](../access/connectors/socrata.md#connection-parameters)     |
| [Sparql](../access/connectors/sparql.md)            |       |        |      |     |     |              |               |     |                   |       |          | [See the list](../access/connectors/sparql.md#connection-parameters)      |
| [SQLite](../access/connectors/sqlite.md)            | JDBC  |        |      | ✓   |     |              |               |     | ✓                 | ✓     | ✓        |                                                                           |
| [Teradata](../access/connectors/teradata.md)        | JDBC  | ✓      | ✓    | ✓   |  ✓   | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [Twitter](../access/connectors/twitter.md)          |       |        |      |     |     |              |               |     |                   |       |          | [See the list](../access/connectors/twitter.md#connection-parameters)     |
| [Vertica](../access/connectors/vertica.md)          | JDBC  | ✓      | ✓    | ✓   |      |✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [Virtuoso](../access/connectors/virtuoso.md)        | JDBC  | ✓      | ✓    | ✓   |     | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                                           |
| [Web](../access/connectors/web.md)                  |       |        |      |     |     |              |               |     |                   |       |          | [See the list](../access/connectors/web.md#connection-parameters)         |

See also:

* [Databases](databases.md)
* [File shares](file-shares.md)
