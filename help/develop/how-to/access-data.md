<!-- TITLE: Access Data -->
<!-- SUBTITLE: -->

# Data Access

This article will walk you through different ways of data extraction, from managing connections to data sources and running queries to reading files into dataframe objects in your application.

## Connections

In the Datagrok platform, you can retrieve data from a variety of sources, be it a file, database, cloud, or web service. And what's remarkable, this can be accomplished equally well from the user interface and from your program. You can learn more about connections and look at how to create them from the UI in the [dedicated article](../../access/data-connection), while here we go over the technical details for developers, namely:

  * how to [add a data connection](#adding-connections)
  * which [parameters](#parameters) to specify
  * secure ways to [transfer credentials](#managing-credentials)
  * creating and executing [parameterized queries](#queries)
  * [sharing connections](#sharing-connections)

### Adding Connections

As with everything else, development starts with a [package](../develop.md#packages). Packages may contain one or more data connections under the `connections` folder. For each connection, you need to create a separate `json` file with the required parameters. Here is an example:

```json
{
  "name": "ChEMBL",
  "parameters": {
    "server": "$GROK_DB_SERVER",
    "db": "chembl_24"
  },
  "dataSource": "PostgreSQL",
  "description": "CHEMBL db",
  "tags": ["demo", "chem"]
}
```

The field `name` is optional, if omitted, the file name (without the extension) will be used as the connection name. In any case, remember that you should not rely on letter case to distinguish between connections, since their names are not case-sensitive. Giving parameters for the connection in `json` is completely equivalent to what you can do from the platform's interface: you would go to `Data | Databases` and right-click on the data source `PostgreSQL` to add such a connection (or, more generally, perform it from `Actions | Add New Connection`).

Our package utilities provide a similar template on running the `grok add connection <name>`command. To see other examples, open [Chembl](https://github.com/datagrok-ai/public/tree/master/packages/Chembl) or [UsageAnalysis](https://github.com/datagrok-ai/public/tree/master/packages/UsageAnalysis) packages in our public repository.

### Parameters

Connection parameters are specific to a data source. However, most of the data providers require some of these common parameters to be specified:

| Data Source                                            | Server  | Port    | DB      | Cache Schema | Cache Results | SSL     | Connection String | Login | Password |  Other Parameters                                                            |
|--------------------------------------------------------|---------|---------|---------|--------------|---------------|---------|-------------------|-------|----------|------------------------------------------------------------------------------|
| [Access](../../access/connectors/access.md)            |         |         | &check; |              |               |         | &check;           |&check;| &check;  |                                                                              |
| [Athena](../../access/connectors/athena.md)            | &check; | &check; | &check; |              |               |         | &check;           |       |          | [See the list](../../access/connectors/athena.md#parameters)                 |
| [BigQuery](../../access/connectors/bigquery.md)        |         |         |         |              |               |         | &check;           |&check;| &check;  | [See the list](../../access/connectors/bigquery.md#connection-parameters)    |
| [Cassandra](../../access/connectors/cassandra.md)      | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           |&check;| &check;  |                                                                              |
| [DB2](../../access/connectors/db2.md)                  | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           |&check;| &check;  |                                                                              |
| [Denodo](../../access/connectors/denodo.md)            | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           |&check;| &check;  |                                                                              |
| [DropBox](../../access/connectors/dropbox.md)          |         |         |         |              |               |         |                   |       | &check;  | [See the list](../../access/connectors/dropbox.md#connection-parameters)     | 
| [Files](../../access/connectors/files.md)              |         |         |         |              |               |         |                   |&check;| &check;  | [See the list](../../access/connectors/files.md#connection-parameters)       |
| [Firebird](../../access/connectors/firebird.md)        | &check; | &check; | &check; | &check;      | &check;       |         | &check;           |&check;| &check;  |                                                                              |
| [Git](../../access/connectors/git.md)                  |         |         |         |              |               |         |                   |       |          | [See the list](../../access/connectors/git.md#connection-parameters)         | 
| [Google Cloud](../../access/connectors/googlecloud.md) |         |         |         |              |               |         |                   |       |          | [See the list](../../access/connectors/googlecloud.md#connection-parameters) | 
| [HBase](../../access/connectors/hbase.md)              | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           |&check;| &check;  |                                                                              |
| [Hive](../../access/connectors/hive.md)                | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           |&check;| &check;  |                                                                              |
| [Hive2](../../access/connectors/hive2.md)              | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           |&check;| &check;  |                                                                              |
| [Impala](../../access/connectors/impala.md)            | &check; | &check; | &check; |              |               |         | &check;           |&check;| &check;  | [See the list](../../access/connectors/impala.md#connection-parameters)      |
| [MariaDB](../../access/connectors/mariadb.md)          | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           |&check;| &check;  |                                                                              |
| [MongoDB](../../access/connectors/mongodb.md)          | &check; | &check; | &check; | &check;      | &check;       |         | &check;           |&check;| &check;  |                                                                              |
| [MS SQL](../../access/connectors/mssql.md)             | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           |&check;| &check;  |                                                                              |
| [MySql](../../access/connectors/mysql.md)              | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           |&check;| &check;  |                                                                              |
| [Neo4j](../../access/connectors/neo4j.md)              | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           |&check;| &check;  |                                                                              |
| [OData](../../access/connectors/odata.md)              |         |         |         |              |               |         |                   |       |          | [See the list](../../access/connectors/odata.md#connection-parameters)       |
| [Oracle](../../access/connectors/oracle.md)            | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           |&check;| &check;  |                                                                              |
| [PostgresNet](../../access/connectors/postgres.md)     | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           |&check;| &check;  |                                                                              |
| [PostgreSQL](../../access/connectors/postgres.md)      | &check; |         | &check; |              | &check;       | &check; |                   |&check;| &check;  |                                                                              |
| [Redshift](../../access/connectors/redshift.md)        | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           |&check;| &check;  |                                                                              |
| [S3](../../access/connectors/s3.md)                    |         |         |         |              |               |         |                   |       |          | [See the list](../../access/connectors/s3.md#connection-parameters)          |
| [Snowflake](../../access/connectors/snowflake.md)      | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           |&check;| &check;  |                                                                              |
| [Socrata](../../access/connectors/socrata.md)          |         |         |         |              |               |         |                   |       |          | [See the list](../../access/connectors/socrata.md#connection-parameters)     |
| [Sparql](../../access/connectors/sparql.md)            |         |         |         |              |               |         |                   |       |          | [See the list](../../access/connectors/sparql.md#connection-parameters)      |
| [SQLite](../../access/connectors/sqlite.md)            |         |         | &check; |              |               |         | &check;           |&check;| &check;  |                                                                              |
| [Teradata](../../access/connectors/teradata.md)        | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           |&check;| &check;  |                                                                              |
| [Twitter](../../access/connectors/twitter.md)          |         |         |         |              |               |         |                   |       |          | [See the list](../../access/connectors/twitter.md#connection-parameters)     |
| [Vertica](../../access/connectors/vertica.md)          | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           |&check;| &check;  |                                                                              |
| [Virtuoso](../../access/connectors/virtuoso.md)        | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           |&check;| &check;  |                                                                              |
| [Web](../../access/connectors/web.md)                  |         |         |         |              |               |         |                   |       |          | [See the list](../../access/connectors/web.md#connection-parameters)         |

When providing a connection string, you do not have to pass other parameters, as they will not be taken into account. The only exception is credentials: for such parameters as `Login` and `Password`, or `Access Key` and `Secret Key`, the case is different. We will cover that in the [next section](#managing-credentials). For now, let's assume that parameters carrying sensitive data form a distinct category. Given the parameters listed above, in the most general case, you would add a connection with the following contents to pull data from a provider:

```json
{
  "name": "Northwind",
  "parameters": {
    "server": "dev.datagrok.ai",
    "port": 23306,
    "db": "Northwind",
    "cache schema": false,
    "cache results": false,
    "SSL": false,
    "conn. string": ""
  },
  "credentials" : {
    "parameters": {
      "login": "",
      "password": ""
    }
  },
  "dataSource": "Maria DB",
  "description": "Northwind Connection",
  "tags": ["demo"]
}
```

### Managing Credentials

Private information is always a special case. Datagrok has a built-in credentials management system that protects such data (for more details, see the [Security](../../govern/security.md#credentials-storage) article). For this reason, parameters of a connection that regulate access to the resource are processed independently. Therefore, you should not include these parameters to a custom connection string. However, pushing them in a `json` file to your project's repository is not a good idea either:

```json
{
  "credentials" : {
    "parameters": {
      "login": "<login>",
      "password": "<password>"
    }
  }
}
```

What you can do instead is to deploy a connection and send a POST request to `$(GROK_HOST)/api/credentials/for/$(PACKAGE_NAME).$(CONNECTION_NAME)` with raw body containing JSON, such as `{"login": "abc", "password": "123"}`, and headers `{"Authorization": $(API_KEY), "Content-Type": "application/json"}` (the API key should be taken from your profile page in Datagrok, e.g., [https://public.datagrok.ai/u](https://public.datagrok.ai/u)).

### Queries

#### Creating Queries

Once the connection is established, the next step is to extract data. This can be done by sending a [query](../../access/data-query) to the data source. In a package, queries are typically placed in the `queries` folder. Let's start with a simple example for your `queries.sql` file:

```sql
--name: protein classification
--connection: chembl
select * from protein_classification
--end
```

SQL statements are annotated with comments, just like [scripts](../scripting.md), since the underlying mechanism is essentially the same (read more on the concept of [functions](../../overview/functions/function.md)). Here we have two header parameters: the query `name` and the `connection` to use. In fact, this particular query could have been even simpler: there is no need to specify `connection` if the package only has one. Similarly, the tag `end` is not required if there is only one query per file: the parser needs it to understand where the current query ends and the next one begins. So safely omit the name of `connection` and/or the `end` tag if these conditions are met.

#### Executing Queries

There are several ways in which queries can be run in Datagrok. The first and most natural way is to launch a query from the interface, which will be equivalent to the line `$(PACKAGE_NAME):$(QUERY_NAME)()` in the console, for example, `Chembl:ProteinClassification()`. As you might know, it is possible to call any function that can be run in the console through Datagrok's JS API. Thus, the fact that a query behaves like a regular function allows us to use the corresponding methods in JavaScript:

```javascript
grok.functions.call('Chembl:ProteinClassification')
  .then(t => grok.shell.addTableView(t));
```

### Sharing Connections

Data connections can be shared as part of a [project](../../overview/project.md), [package](../develop.md#packages) (and [repository](../../access/connectors/git.md) containing this package), or as an independent [entity](../../overview/objects.md). Access rights of a database connection inherit access rights of a query. However, access rights of the query don't inherit access rights of the database connection. Thus, if one shares a query, the associated database connection shall automatically be shared. At the same time, when you are sharing a connection, your queries aren't going to be shared automatically. As for web queries, they are automatically shared along with sharing the corresponding connection.

See also:
  * [Data Connection](../../access/data-connection)
  * [Data Query](../../access/data-query)
  * [Functions](../../overview/functions/function.md)
  * [Parameterized Queries](../../access/parameterized-queries.md)
