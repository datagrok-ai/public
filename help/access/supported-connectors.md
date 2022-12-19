# Supported connectors

A connector can work with a database, an Excel file, a CSV file, a web service,
or basically anything that is capable of providing the data. Datagrok currently
supports over 30 different connectors, and the list is quickly growing. Most of
Datagrok’s data connectors are open-sourced and could be found on
[GitHub](https://github.com/datagrok-ai/public/tree/master/connectors) (MIT
license). The supported connectors with their specific parameters are the
following:

| Data Source                                         | Server  | Port    | DB      | Cache Schema | Cache Results | SSL     | Connection String | Login   | Password | Other Parameters                                                          |
|-----------------------------------------------------|---------|---------|---------|--------------|---------------|---------|-------------------|---------|----------|---------------------------------------------------------------------------|
| [Access](../access/connectors/access.md)            |         |         | &check; |              |               |         | &check;           | &check; | &check;  |                                                                           |
| [Athena](../access/connectors/athena.md)            | &check; | &check; | &check; |              |               |         | &check;           |         |          | [See the list](../access/connectors/athena.md)                            |
| [BigQuery](../access/connectors/bigquery.md)        |         |         |         |              |               |         | &check;           | &check; | &check;  | [See the list](../access/connectors/bigquery.md#connection-parameters)    |
| [Cassandra](../access/connectors/cassandra.md)      | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                           |
| [DB2](../access/connectors/db2.md)                  | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                           |
| [Denodo](../access/connectors/denodo.md)            | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                           |
| [DropBox](../access/connectors/dropbox.md)          |         |         |         |              |               |         |                   |         | &check;  | [See the list](../access/connectors/dropbox.md#connection-parameters)     |
| [Files](../access/connectors/files.md)              |         |         |         |              |               |         |                   | &check; | &check;  | [See the list](../access/connectors/files.md#connection-parameters)       |
| [Firebird](../access/connectors/firebird.md)        | &check; | &check; | &check; | &check;      | &check;       |         | &check;           | &check; | &check;  |                                                                           |
| [Git](../access/connectors/git.md)                  |         |         |         |              |               |         |                   |         |          | [See the list](../access/connectors/git.md#connection-parameters)         |
| [Google Cloud](../access/connectors/googlecloud.md) |         |         |         |              |               |         |                   |         |          | [See the list](../access/connectors/googlecloud.md#connection-parameters) |
| [HBase](../access/connectors/hbase.md)              | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                           |
| [Hive](../access/connectors/hive.md)                | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                           |
| [Hive2](../access/connectors/hive2.md)              | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                           |
| [Impala](../access/connectors/impala.md)            | &check; | &check; | &check; |              |               |         | &check;           | &check; | &check;  | [See the list](../access/connectors/impala.md#connection-parameters)      |
| [MariaDB](../access/connectors/mariadb.md)          | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                           |
| [MongoDB](../access/connectors/mongodb.md)          | &check; | &check; | &check; | &check;      | &check;       |         | &check;           | &check; | &check;  |                                                                           |
| [MS SQL](../access/connectors/mssql.md)             | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                           |
| [MySql](../access/connectors/mysql.md)              | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                           |
| [Neo4j](../access/connectors/neo4j.md)              | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                           |
| [OData](../access/connectors/odata.md)              |         |         |         |              |               |         |                   |         |          | [See the list](../access/connectors/odata.md#connection-parameters)       |
| [Oracle](../access/connectors/oracle.md)            | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                           |
| [Postgres](../access/connectors/postgres.md)     | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                           |
| [PostgreSQL](../access/connectors/postgres.md)      | &check; |         | &check; |              | &check;       | &check; |                   | &check; | &check;  |                                                                           |
| [Redshift](../access/connectors/redshift.md)        | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                           |
| [S3](../access/connectors/s3.md)                    |         |         |         |              |               |         |                   |         |          | [See the list](../access/connectors/s3.md#connection-parameters)          |
| [Snowflake](../access/connectors/snowflake.md)      | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                           |
| [Socrata](../access/connectors/socrata.md)          |         |         |         |              |               |         |                   |         |          | [See the list](../access/connectors/socrata.md#connection-parameters)     |
| [Sparql](../access/connectors/sparql.md)            |         |         |         |              |               |         |                   |         |          | [See the list](../access/connectors/sparql.md#connection-parameters)      |
| [SQLite](../access/connectors/sqlite.md)            |         |         | &check; |              |               |         | &check;           | &check; | &check;  |                                                                           |
| [Teradata](../access/connectors/teradata.md)        | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                           |
| [Twitter](../access/connectors/twitter.md)          |         |         |         |              |               |         |                   |         |          | [See the list](../access/connectors/twitter.md#connection-parameters)     |
| [Vertica](../access/connectors/vertica.md)          | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                           |
| [Virtuoso](../access/connectors/virtuoso.md)        | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                           |
| [Web](../access/connectors/web.md)                  |         |         |         |              |               |         |                   |         |          | [See the list](../access/connectors/web.md#connection-parameters)         |

See also:

* [Databases](databases.md)
* [File shares](file-shares.md)
