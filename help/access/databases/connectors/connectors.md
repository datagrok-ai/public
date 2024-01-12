---
title: "Supported connectors"
---


A connector can work with a database, an Excel file, a CSV file, a web service,
or basically anything that is capable of providing the data. Datagrok currently
supports over 30 different connectors, and the list is quickly growing. Most of
Datagrok's data connectors are open-sourced and could be found on
[GitHub](https://github.com/datagrok-ai/public/tree/master/connectors) (MIT
license). The supported connectors with their specific parameters are the
following:

| Data Source                      | Type  | Server | Port | DB  | Browse Schema | Cache Schema | Cache Results | SSL | Connection String | Login | Password | Other Parameters                                       |
|----------------------------------|-------|--------|------|-----|---------------|--------------|---------------|-----|-------------------|-------|----------|--------------------------------------------------------|
| [Access](./access.md)            | JDBC  |        |      | ✓   | ✓             | ✓            |               |     | ✓                 | ✓     | ✓        |                                                        |
| [Athena](./athena.md)            | JDBC  | ✓      | ✓    | ✓   | ✓             | ✓            |               |     | ✓                 |       |          | [See the list](./athena.md)                            |
| [BigQuery](./bigquery.md)        | JDBC  |        |      |     |               |              |               |     | ✓                 | ✓     | ✓        | [See the list](./bigquery.md#connection-parameters)    |
| [Cassandra](./cassandra.md)      | JDBC  | ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [ClickHouse](./clickhouse.md)    | JDBC  | ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [DB2](./db2.md)                  | JDBC  | ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Denodo](./denodo.md)            | JDBC  | ✓      | ✓    | ✓   |               | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Firebird](./firebird.md)        | JDBC  | ✓      | ✓    | ✓   |               | ✓            | ✓             |     | ✓                 | ✓     | ✓        |                                                        |
| [HBase](./hbase.md)              | JDBC  | ✓      | ✓    | ✓   |               | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Hive](./hive.md)                | JDBC  | ✓      | ✓    | ✓   |               | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Hive2](./hive2.md)              | JDBC  | ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Impala](./impala.md)            | JDBC  | ✓      | ✓    | ✓   | ✓             |              |               |     | ✓                 | ✓     | ✓        | [See the list](./impala.md#connection-parameters)      |
| [MariaDB](./mariadb.md)          | JDBC  | ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [MongoDB](./mongodb.md)          | JDBC  | ✓      | ✓    | ✓   |               | ✓            | ✓             |     | ✓                 | ✓     | ✓        ||
| [MS SQL](./mssql.md)             | JDBC  | ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        ||
| [MySql](./mysql.md)              | JDBC  | ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Neo4j](./neo4j.md)              | JDBC  | ✓      | ✓    | ✓   |               | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Neptune](./neptune.md)          | JDBC  | ✓      | ✓    |     |               | ✓            | ✓             |     | ✓                 |       |          | [See the list](./neptune.md#connection-parameters)     |
| [OData](./odata.md)              |       |        |      |     |               |              |               |     |                   |       |          | [See the list](./odata.md#connection-parameters)       |
| [Oracle](./oracle.md)            | JDBC  | ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [PI](./pi.md)                    | JDBC  | ✓      |      | ✓   |               | ✓            | ✓             |     | ✓                 | ✓     | ✓        | [See the list](./pi.md#connection-parameters)          |
| [Postgres](./postgres.md)        | JDBC  | ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [PostgresDart](./postgres.md)    | JDBC  | ✓      |      | ✓   | ✓             |              | ✓             | ✓   |                   | ✓     | ✓        |                                                        |
| [Redshift](./redshift.md)        | JDBC  | ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [SAP HANA](./sap-hana.md)        | JDBC  | ✓      | ✓    | ✓   | ✓             |              |               |     |                   | ✓     | ✓        | [See the list](./sap-hana.md#connection-parameters)    |
| [Snowflake](./snowflake.md)      | JDBC  | ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Socrata](./socrata.md)          |       |        |      |     |               |              |               |     |                   |       |          | [See the list](./socrata.md#connection-parameters)     |
| [Sparql](./sparql.md)            |       |        |      |     |               |              |               |     |                   |       |          | [See the list](./sparql.md#connection-parameters)      |
| [SQLite](./sqlite.md)            | JDBC  |        |      | ✓   | ✓             | ✓            |               |     | ✓                 | ✓     | ✓        |                                                        |
| [Teradata](./teradata.md)        | JDBC  | ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Twitter](./twitter.md)          |       |        |      |     |               |              |               |     |                   |       |          | [See the list](./twitter.md#connection-parameters)     |
| [Vertica](./vertica.md)          | JDBC  | ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Virtuoso](./virtuoso.md)        | JDBC  | ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Web](./web.md)                  |       |        |      |     |               |              |               |     |                   |       |          | [See the list](./web.md#connection-parameters)         |

See also:

* [Databases](../databases.mdx)
* [File shares](../../files/files.mdx)
