# Supported connectors

A connector can work with a database, an Excel file, a CSV file, a web service,
or basically anything that is capable of providing the data. Datagrok currently
supports over 30 different connectors, and the list is quickly growing. Most of
Datagrok’s data connectors are open-sourced and could be found on
[GitHub](https://github.com/datagrok-ai/public/tree/master/connectors) (MIT
license). The supported connectors with their specific parameters are the
following:

| Data Source                                            | Server  | Port    | DB      | Cache Schema | Cache Results | SSL     | Connection String | Login   | Password | Other Parameters                                                             |
|--------------------------------------------------------|---------|---------|---------|--------------|---------------|---------|-------------------|---------|----------|------------------------------------------------------------------------------|
| [Access]( access.md)            |         |         | &check; |              |               |         | &check;           | &check; | &check;  |                                                                              |
| [Athena]( athena.md)            | &check; | &check; | &check; |              |               |         | &check;           |         |          | [See the list]( athena.md)                            |
| [BigQuery]( bigquery.md)        |         |         |         |              |               |         | &check;           | &check; | &check;  | [See the list]( bigquery.md#connection-parameters)    |
| [Cassandra]( cassandra.md)      | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [DB2]( db2.md)                  | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [Denodo]( denodo.md)            | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [DropBox]( dropbox.md)          |         |         |         |              |               |         |                   |         | &check;  | [See the list]( dropbox.md#connection-parameters)     |
| [Files](files.md)              |         |         |         |              |               |         |                   | &check; | &check;  | [See the list]( files.md#connection-parameters)       |
| [Firebird]( firebird.md)        | &check; | &check; | &check; | &check;      | &check;       |         | &check;           | &check; | &check;  |                                                                              |
| [Git]( git.md)                  |         |         |         |              |               |         |                   |         |          | [See the list]( git.md#connection-parameters)         |
| [Google Cloud]( googlecloud.md) |         |         |         |              |               |         |                   |         |          | [See the list]( googlecloud.md#connection-parameters) |
| [HBase]( hbase.md)              | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [Hive]( hive.md)                | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [Hive2]( hive2.md)              | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [Impala]( impala.md)            | &check; | &check; | &check; |              |               |         | &check;           | &check; | &check;  | [See the list]( impala.md#connection-parameters)      |
| [MariaDB]( mariadb.md)          | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [MongoDB]( mongodb.md)          | &check; | &check; | &check; | &check;      | &check;       |         | &check;           | &check; | &check;  |                                                                              |
| [MS SQL]( mssql.md)             | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [MySql]( mysql.md)              | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [Neo4j]( neo4j.md)              | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [OData]( odata.md)              |         |         |         |              |               |         |                   |         |          | [See the list]( odata.md#connection-parameters)       |
| [Oracle]( oracle.md)            | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [PostgresNet]( postgres.md)     | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [PostgreSQL]( postgres.md)      | &check; |         | &check; |              | &check;       | &check; |                   | &check; | &check;  |                                                                              |
| [Redshift]( redshift.md)        | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [S3]( s3.md)                    |         |         |         |              |               |         |                   |         |          | [See the list]( s3.md#connection-parameters)          |
| [Snowflake]( snowflake.md)      | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [Socrata]( socrata.md)          |         |         |         |              |               |         |                   |         |          | [See the list]( socrata.md#connection-parameters)     |
| [Sparql]( sparql.md)            |         |         |         |              |               |         |                   |         |          | [See the list]( sparql.md#connection-parameters)      |
| [SQLite]( sqlite.md)            |         |         | &check; |              |               |         | &check;           | &check; | &check;  |                                                                              |
| [Teradata]( teradata.md)        | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [Twitter]( twitter.md)          |         |         |         |              |               |         |                   |         |          | [See the list]( twitter.md#connection-parameters)     |
| [Vertica]( vertica.md)          | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [Virtuoso]( virtuoso.md)        | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [Web]( web.md)                  |         |         |         |              |               |         |                   |         |          | [See the list]( web.md#connection-parameters)         |

See also:

* [Databases](../databases.md)
* [File shares](../file-shares.md)
