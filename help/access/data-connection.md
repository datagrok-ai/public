<!-- TITLE: Data Connection -->
<!-- SUBTITLE: -->

# Data Connection

Data connection is used for accessing data in a particular data source. 
Connection parameters depend on the data source. Typically, you would need to provide 
server name and login credentials.

## Connectors

A connector could work with a database, an Excel file, a CSV file, a web service, or basically 
anything that is capable of providing the data. We currently support over 20 different 
connectors, and the list is quickly growing. Most of our data connectors are open-sourced and could be found on 
[GitHub](https://github.com/datagrok-ai/public/tree/master/connectors) (MIT license).

| Name                                     | Type  |
|------------------------------------------|-------|
| [Access](connectors/access.md)           | JDBC  |
| [Athena](connectors/athena.md)           | JDBC  |
| [BigQuery](connectors/bigquery.md)       | JDBC  |
| [Cassandra](connectors/cassandra.md)     | JDBC  |
| [DB2](connectors/db2.md)                 | JDBC  |
| [Denodo](connectors/denodo.md)           | JDBC  |
| [DropBox](connectors/dropbox.md)         | Files |
| [Files](connectors/files.md)             | Files |
| [Firebird](connectors/firebird.md)       | JDBC  |
| [Git](connectors/git.md)                 | Files |
| [GitHub](connectors/github.md)           | Files |
| [GoogleCloud](connectors/googlecloud.md) | Files |
| [HBase](connectors/hbase.md)             | JDBC  |
| [Hive](connectors/hive.md)               | JDBC  |
| [Hive2](connectors/hive2.md)             | JDBC  |
| [Impala](connectors/impala.md)           | JDBC  |
| [MariaDB](connectors/mariadb.md)         | JDBC  |
| [MongoDB](connectors/mssql.md)           | JDBC  |
| [MS SQL](connectors/mongodb.md)          | JDBC  |
| [MySql](connectors/mysql.md)             | JDBC  |
| [Neo4j](connectors/neo4j.md)             | JDBC  |
| [OData](connectors/odata.md)             |       |
| [Oracle](connectors/oracle.md)           | JDBC  |
| [Postgres](connectors/postgres.md)       | JDBC  |
| [Redshift](connectors/postgres.md)       | JDBC  |
| [S3](connectors/s3.md)                   | Files |
| [Socrata](connectors/socrata.md)         |       |
| [Sparql](connectors/sparql.md)           |       |
| [SQLite](connectors/sqlite.md)           | JDBC  |
| [Teradata](connectors/teradata.md)       | JDBC  |
| [Twitter](connectors/twitter.md)         |       |
| [Vertica](connectors/vertica.md)         | JDBC  |
| [Virtuoso](connectors/virtuoso.md)       | JDBC  |
| [Web](connectors/web.md)                 |       |


## Creating a new connection

To create a new data connection, open the "Databases" pane (Open | Databases), right-click on the appropriate 
connector in the tree, and choose "Add connection...". Alternatively, click on "New Connection" under the
"Actions" tab, and select the appropriate connector.

![](data-connection-tree.png)

Then, edit the connection attributes, and click on TEST to confirm that you've entered everything
correctly. The set of attributes you can edit depends on the connector. Typically, for JDBC-based connectors
you can provide a custom connection string (but do not enter login and password there, they will still
be picked up from the corresponding fields).   

![](data-connection-create.png)
        
Once a connection is set up, you are ready to start creating queries. There are multiple ways to 
do so: manually, or         
        
## Access control

Connections are first-class entities in the Datagrok platform, and as such are subjects to the standard 
checks and routines performed against them whenever they are used in the specific context. Some of the 
most popular privileges are: `view`, `edit`, `delete`, and `share`. 
Those privileges can be given to individual users, or to [groups](../govern/group.md).
For more information on the access privilege model, check out [privileges](../govern/security.md#privileges).

Another “out of the box” feature that comes with connections being first-class entity is the audit trail for 
every action performed against the connection. For details on that, check out [Audit](../govern/audit.md) page.

## Filtering

You can use these fields to filter connections with [smart search](../overview/smart-search.md):

| Field       | Description                                 |
|-------------|---------------------------------------------|
| id          |                                             |
| name        |                                             |
| server      |                                             |
| port        |                                             |
| db          |                                             |
| login       |                                             |
| dataSource  |                                             |
| description |                                             |
| createdOn   |                                             |
| updatedOn   |                                             |
| author      | [User](../govern/user.md) object            |
| starredBy   | [User](../govern/user.md) object            |
| commentedBy | [User](../govern/user.md) object            |
| usedBy      | [User](../govern/user.md) object            |

## JDBC connection

For some cases connection may require custom JDBC connection string. For this case JDBC-based data connection has 
parameter "Conn. string". If it filled, it will be used for connection, all other parameters will be ignored except 
'login' and 'password'.  

See also:

* [Data Pipeline](data-pipeline.md)
* [Data Query](data-query.md)
* [Data Job](data-job.md)
* [Function Call](../overview/functions/function-call.md)
