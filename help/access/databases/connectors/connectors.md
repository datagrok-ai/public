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

| Data Source                      |  Server | Port | DB  | Browse Schema | Cache Schema | Cache Results | SSL | Connection String | Login | Password | Other Parameters                                       |
|----------------------------------|-------|------|-----|---------------|--------------|---------------|-----|-------------------|-------|----------|--------------------------------------------------------|
| [Access](./access.md)            |         |      | ✓   | ✓             | ✓            |               |     | ✓                 | ✓     | ✓        |                                                        |
| [Athena](./athena.md)            |  ✓      | ✓    | ✓   | ✓             | ✓            |               |     | ✓                 |       |          | [See the list](./athena.md)                            |
| [BigQuery](./bigquery.md)        |         |      |     |               |              |               |     | ✓                 | ✓     | ✓        | [See the list](./bigquery.md#connection-parameters)    |
| [Cassandra](./cassandra.md)      |  ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [ClickHouse](./clickhouse.md)    |  ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [DB2](./db2.md)                  |  ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Denodo](./denodo.md)            |  ✓      | ✓    | ✓   |               | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Firebird](./firebird.md)        |  ✓      | ✓    | ✓   |               | ✓            | ✓             |     | ✓                 | ✓     | ✓        |                                                        |
| [HBase](./hbase.md)              |  ✓      | ✓    | ✓   |               | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Hive](./hive.md)                |  ✓      | ✓    | ✓   |               | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Hive2](./hive2.md)              |  ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Impala](./impala.md)            |  ✓      | ✓    | ✓   | ✓             |              |               |     | ✓                 | ✓     | ✓        | [See the list](./impala.md#connection-parameters)      |
| [MariaDB](./mariadb.md)          |  ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [MongoDB](./mongodb.md)          |  ✓      | ✓    | ✓   |               | ✓            | ✓             |     | ✓                 | ✓     | ✓        ||
| [MS SQL](./mssql.md)             |  ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        ||
| [MySql](./mysql.md)              |  ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Neo4j](./neo4j.md)              |  ✓      | ✓    | ✓   |               | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Neptune](./neptune.md)          |  ✓      | ✓    |     |               | ✓            | ✓             |     | ✓                 |       |          | [See the list](./neptune.md#connection-parameters)     |
| [OData](./odata.md)                |        |      |     |               |              |               |     |                   |       |          | [See the list](./odata.md#connection-parameters)       |
| [Oracle](./oracle.md)            |  ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [PI](./pi.md)                    |  ✓      |      | ✓   |               | ✓            | ✓             |     | ✓                 | ✓     | ✓        | [See the list](./pi.md#connection-parameters)          |
| [Postgres](./postgres.md)        |  ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [PostgresDart](./postgres.md)    |  ✓      |      | ✓   | ✓             |              | ✓             | ✓   |                   | ✓     | ✓        |                                                        |
| [Redshift](./redshift.md)        |  ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [SAP HANA](./sap-hana.md)        |  ✓      | ✓    | ✓   | ✓             |              |               |     |                   | ✓     | ✓        | [See the list](./sap-hana.md#connection-parameters)    |
| [Snowflake](./snowflake.md)      |  ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Socrata](./socrata.md)             |        |      |     |               |              |               |     |                   |       |          | [See the list](./socrata.md#connection-parameters)     |
| [Sparql](./sparql.md)            |          |      |     |               |              |               |     |                   |       |          | [See the list](./sparql.md#connection-parameters)      |
| [SQLite](./sqlite.md)            |         |      | ✓   | ✓             | ✓            |               |     | ✓                 | ✓     | ✓        |                                                        |
| [Teradata](./teradata.md)        |  ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Twitter](./twitter.md)           |        |      |     |               |              |               |     |                   |       |          | [See the list](./twitter.md#connection-parameters)     |
| [Vertica](./vertica.md)          |  ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Virtuoso](./virtuoso.md)        |  ✓      | ✓    | ✓   | ✓             | ✓            | ✓             | ✓   | ✓                 | ✓     | ✓        |                                                        |
| [Web](./web.md)                    |        |      |     |               |              |               |     |                   |       |          | [See the list](./web.md#connection-parameters)         |

See also:

* [Databases](../databases.md)
* [File shares](../../files/files.md)
