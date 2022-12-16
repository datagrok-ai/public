<!-- TITLE: Access data -->
<!-- SUBTITLE: -->

# Data access

This article will walk you through different ways of data extraction, from [managing connections](#connections) to data
sources and running queries, to [using REST APIs](#rest-endpoints), and to [reading files](#reading-files) into
dataframe objects in your application.

## Connections

In the Datagrok platform, you can retrieve data from a variety of sources, be it a file, database, cloud, or web
service. And what's remarkable, this can be accomplished equally well from the user interface and from your program. You
can learn more about connections and look at how to create them from the UI in
the [dedicated article](../../access/data-connection.md), while here we go over the technical details for developers,
namely:

* how to [add a data connection](#adding-connections)
* which [parameters](#parameters) to specify
* secure ways to [transfer credentials](#managing-credentials)
* [creating](#creating-queries) and [executing parameterized queries](#executing-queries)
* [sharing connections](#sharing-connections)

Let's outline a general workflow for accessing data using connections. You will add a connection to your package, use it
to query the database, and get a dataframe after running the query in JavaScript code. This is how you obtain the data
for further work.

### Adding connections

As with everything else, development starts with a [package](../develop.md#packages). Packages may contain one or more
data connections under the `connections` folder. For each connection, you need to create a separate `json` file with the
required parameters. Here is an example:

```json
{
  "name": "ChEMBL",
  "parameters": {
    "server": "$GROK_DB_SERVER",
    "db": "chembl_24"
  },
  "dataSource": "PostgresDart",
  "description": "CHEMBL db",
  "tags": [
    "demo",
    "chem"
  ]
}
```

The field `name` is optional, if omitted, the filename (without the extension) will be used as the connection name. In
any case, remember that you should not rely on letter case to distinguish between connections, since their names are not
case-sensitive. Giving parameters for the connection in `json` is completely equivalent to what you can do from the
platform's interface: you would go to `Data | Databases` and right-click on the data source `PostgreSQL` to add such a
connection (or, more generally, perform it from `Actions | Add New Connection`).

Our package utilities provide a similar template on running the `grok add connection <name>`
command. To see other examples, open [Chembl](https://github.com/datagrok-ai/public/tree/master/packages/Chembl)
or [UsageAnalysis](https://github.com/datagrok-ai/public/tree/master/packages/UsageAnalysis)
packages in our public repository.

### Parameters

Connection parameters are specific to a data source. However, most of the data providers require some of these common
parameters to be specified:

| Data Source                                            | Server  | Port    | DB      | Cache Schema | Cache Results | SSL     | Connection String | Login   | Password | Other Parameters                                                             |
|--------------------------------------------------------|---------|---------|---------|--------------|---------------|---------|-------------------|---------|----------|------------------------------------------------------------------------------|
| [Access](../../access/connectors/access.md)            |         |         | &check; |              |               |         | &check;           | &check; | &check;  |                                                                              |
| [Athena](../../access/connectors/athena.md)            | &check; | &check; | &check; |              |               |         | &check;           |         |          | [See the list](../../access/connectors/athena.md)                            |
| [BigQuery](../../access/connectors/bigquery.md)        |         |         |         |              |               |         | &check;           | &check; | &check;  | [See the list](../../access/connectors/bigquery.md#connection-parameters)    |
| [Cassandra](../../access/connectors/cassandra.md)      | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [DB2](../../access/connectors/db2.md)                  | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [Denodo](../../access/connectors/denodo.md)            | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [DropBox](../../access/connectors/dropbox.md)          |         |         |         |              |               |         |                   |         | &check;  | [See the list](../../access/connectors/dropbox.md#connection-parameters)     |
| [Files](../../access/connectors/files.md)              |         |         |         |              |               |         |                   | &check; | &check;  | [See the list](../../access/connectors/files.md#connection-parameters)       |
| [Firebird](../../access/connectors/firebird.md)        | &check; | &check; | &check; | &check;      | &check;       |         | &check;           | &check; | &check;  |                                                                              |
| [Git](../../access/connectors/git.md)                  |         |         |         |              |               |         |                   |         |          | [See the list](../../access/connectors/git.md#connection-parameters)         |
| [Google Cloud](../../access/connectors/googlecloud.md) |         |         |         |              |               |         |                   |         |          | [See the list](../../access/connectors/googlecloud.md#connection-parameters) |
| [HBase](../../access/connectors/hbase.md)              | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [Hive](../../access/connectors/hive.md)                | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [Hive2](../../access/connectors/hive2.md)              | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [Impala](../../access/connectors/impala.md)            | &check; | &check; | &check; |              |               |         | &check;           | &check; | &check;  | [See the list](../../access/connectors/impala.md#connection-parameters)      |
| [MariaDB](../../access/connectors/mariadb.md)          | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [MongoDB](../../access/connectors/mongodb.md)          | &check; | &check; | &check; | &check;      | &check;       |         | &check;           | &check; | &check;  |                                                                              |
| [MS SQL](../../access/connectors/mssql.md)             | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [MySql](../../access/connectors/mysql.md)              | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [Neo4j](../../access/connectors/neo4j.md)              | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [OData](../../access/connectors/odata.md)              |         |         |         |              |               |         |                   |         |          | [See the list](../../access/connectors/odata.md#connection-parameters)       |
| [Oracle](../../access/connectors/oracle.md)            | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [Postgres](../../access/connectors/postgres.md)        | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [PostgresDart](../../access/connectors/postgres.md)    | &check; |         | &check; |              | &check;       | &check; |                   | &check; | &check;  |                                                                              |
| [Redshift](../../access/connectors/redshift.md)        | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [S3](../../access/connectors/s3.md)                    |         |         |         |              |               |         |                   |         |          | [See the list](../../access/connectors/s3.md#connection-parameters)          |
| [Snowflake](../../access/connectors/snowflake.md)      | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [Socrata](../../access/connectors/socrata.md)          |         |         |         |              |               |         |                   |         |          | [See the list](../../access/connectors/socrata.md#connection-parameters)     |
| [Sparql](../../access/connectors/sparql.md)            |         |         |         |              |               |         |                   |         |          | [See the list](../../access/connectors/sparql.md#connection-parameters)      |
| [SQLite](../../access/connectors/sqlite.md)            |         |         | &check; |              |               |         | &check;           | &check; | &check;  |                                                                              |
| [Teradata](../../access/connectors/teradata.md)        | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [Twitter](../../access/connectors/twitter.md)          |         |         |         |              |               |         |                   |         |          | [See the list](../../access/connectors/twitter.md#connection-parameters)     |
| [Vertica](../../access/connectors/vertica.md)          | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [Virtuoso](../../access/connectors/virtuoso.md)        | &check; | &check; | &check; | &check;      | &check;       | &check; | &check;           | &check; | &check;  |                                                                              |
| [Web](../../access/connectors/web.md)                  |         |         |         |              |               |         |                   |         |          | [See the list](../../access/connectors/web.md#connection-parameters)         |

When providing a connection string, you do not have to pass other parameters, as they will not be taken into account.
The only exception is credentials: for such parameters as `Login` and `Password`
, or `Access Key` and `Secret Key`, the case is different. We will cover that in
the [next section](#managing-credentials). For now, let's assume that parameters carrying sensitive data form a distinct
category. Given the parameters listed above, in the most general case, you would add a connection with similar contents
to pull data from a provider:

```json
{
  "name": "Northwind",
  "parameters": {
    "server": "dev.datagrok.ai",
    "port": 23306,
    "db": "Northwind",
    "cacheSchema": false,
    "cacheResults": false,
    "ssl": false,
    "connString": ""
  },
  "credentials": {
    "parameters": {
      "login": "",
      "password": ""
    }
  },
  "dataSource": "MariaDB",
  "description": "Northwind Connection",
  "tags": [
    "demo"
  ]
}
```

### Managing credentials

Private information is always a special case. Datagrok has a built-in credentials management system that protects such
data (for more details, see the [Security](../../govern/security.md#credentials-storage) article). For this reason,
parameters of a connection that regulate access to the resource are processed independently. Therefore, you should not
include these parameters to a custom connection string. However, pushing them in a `json`
file to your project's repository is not a good idea either:

```json
{
  "credentials": {
    "parameters": {
      "login": "<login>",
      "password": "<password>"
    }
  }
}
```

What you can do instead is to deploy a connection and send a POST request
to `$(GROK_HOST)/api/credentials/for/$(PACKAGE_NAME).$(CONNECTION_NAME)` with raw body containing JSON, such
as `{"login": "abc", "password": "123"}`, and
headers `{"Authorization": $(API_KEY), "Content-Type": "application/json"}` (the API key should be taken from your
profile page in Datagrok, e.g., [https://public.datagrok.ai/u](https://public.datagrok.ai/u)).

### Queries

#### Creating queries

Once the connection is established, the next step is to extract data. This can be done by sending
a [query](../../access/data-query.md) to the data source. In a package, queries are typically placed in the `queries`
folder. Let's start with a simple example for your `queries.sql` file:

```sql
--name: protein classification
--connection: chembl
select * from protein_classification
--end
```

SQL statements are annotated with comments, just like [scripts](../../compute/scripting.md), since the underlying
mechanism is essentially the same (read more on the concept of [functions](../../datagrok/functions/function.md)). Here
we have two header parameters: the query `name` and the `connection` to use. In fact, this particular query could have
been even simpler: there is no need to specify `connection` if the package only has one. Similarly, the tag `end` is not
required if there is only one query per file: the parser needs it to understand where the current query ends and the
next one begins. So safely omit the name of `connection` and/or the `end` tag if these conditions are met.

To use an existing connection in a query, specify its name along with the namespace in the `connection` parameter. For
example, the above `chembl` connection that lives in the `Chembl`
package has the following path: `chembl:chembl`. When browsing connections on the platform, you can identify such path
by opening `Links` in the tab `Details` of the property panel.

You can find a list of header parameters and other details related to the query annotation
in [this article](../../access/parameterized-queries.md). In addition to this, examples of data queries are available in
the [Chembl](https://github.com/datagrok-ai/public/tree/master/packages/Chembl/queries) package. To quickly insert a
query template into your package, type `grok add query <name>` in the terminal.

#### Executing queries

There are several ways in which queries can be run in Datagrok. The first and most natural way is to launch a query from
the interface, which will be equivalent to the line `$(PACKAGE_NAME):$(QUERY_NAME)()` in the console, for
example, `Chembl:ProteinClassification()`
. As you might know, it is possible to call any function that can be run in the console through Datagrok's JS API. Thus,
the fact that a query behaves like a regular function allows us to use the corresponding method in JavaScript:

```javascript
grok.functions.call('Chembl:ProteinClassification')
  .then(t => grok.shell.addTableView(t));
```

There is also a special method for queries that takes the query name as a required parameter and a few additional ones (
query parameters, whether it is ad hoc or not, and the polling interval):

```javascript
grok.data.query(`${PACKAGE_NAME}:${QUERY_NAME}`, {'parameter': 'value'}, true, 100);
```

See how this method works in the [example](https://public.datagrok.ai/js/samples/data-access/parameterized-query) on
Datagrok.

### Sharing connections

Data connections can be shared as part of a [project](../../datagrok/project.md)
, [package](../develop.md#packages) (
and [repository](../../access/connectors/git.md) containing this package), or as an
independent [entity](../../datagrok/objects.md). Access rights of a database connection inherit access rights of a
query. However, access rights of the query don't inherit access rights of the database connection. Thus, if one shares a
query, the associated database connection shall automatically be shared. At the same time, when you are sharing a
connection, your queries aren't going to be shared automatically. As for web queries, they are automatically shared
along with sharing the corresponding connection.

### Caching results

See [db query caching](../../access/data-connection.md#caching) for details

## Rest endpoints

Web services provide endpoints that you can programmatically connect to. There are two main options for this: the first
is to use [OpenAPI/Swagger](https://swagger.io/docs/specification/about/) format supported by Datagrok, the second one
involves the
[use of the platform's server](https://dev.datagrok.ai/js/samples/dapi/fetch) to send a network request. The details of
Swagger-based connections are further explained in the [dedicated article](../../access/open-api.md). The method used to
proxy requests has an interface similar to the standard `fetch` API and can be applied as follows:

```javascript
const url = 'https://jsonplaceholder.typicode.com/posts';
const data = {name: 'username', password: 'password'};

grok.dapi.fetchProxy(url, {
  method: 'POST',
  headers: {'Content-Type': 'application/json'},
  body: JSON.stringify(data)
}).then(response => grok.shell.info(response.ok));
```

## Reading files

In addition to connections and queries, it is useful to look at methods for reading files. In
our [JavaScript API examples](https://public.datagrok.ai/js), you can find methods that provide data for
demonstration (`grok.data.demo` or `grok.data.getDemoTable`) and testing purposes (`grok.data.testData`). However, here
we will cover the part that you will often use to deliver data to your applications.

### Package files

If you want to access data that resides within your package (e.g. open a `csv` table), load it from URL, as you would do
with other [external files](https://public.datagrok.ai/js/samples/data-access/external/stock-prices). The package root
for client-side can be found with `webRoot` property. The example shown below gets the `test.csv` table from
the `data-samples` subdirectory and opens a table view for it:

```javascript
export let _package = new DG.Package();

grok.data.loadTable(`${_package.webRoot}data-samples/test.csv`)
  .then(t => grok.shell.addTableView(t));
```

### File shares

Connecting to [file shares](../../access/connect-a-file-share.md) offers you more opportunities: files form a hierarchy,
which you can browse naturally from the interface. Let's start with an existing file share — the user's home directory.
In Datagrok's [file browser](https://public.datagrok.ai/files), each user has a special `HOME` folder to store their
files, which makes it a perfect example. Here is how you can work with your files located there:

```javascript
grok.functions.eval(`OpenServerFile("${USER}:Home/data.csv")`)
  .then(t => grok.shell.addTableView(t[0]));
```

But more importantly, file shares let you gain access to data from various locations. So let's find out how to create
such a connection in your package. Again, it's just a `json` file with the corresponding parameters:

```json
{
  "name": "New File Share",
  "parameters": {
    "dir": "/home/www/master/servergrok/data/demo",
    "indexFiles": true
  },
  "credentials": {
    "parameters": {
      "login": "",
      "password": ""
    }
  },
  "dataSource": "Files",
  "tags": [
    "demo"
  ]
}
```

You should specify two parameters for connection: the directory you are going to work with and whether you want to index
its files (if you do, there will be an indexing data job). When referring to a file from your code, put the names of
package and connection before the file. The path to it should be relative to what you previously specified in the `dir`
parameter:

```javascript
grok.functions.eval(`OpenServerFile("${PACKAGE_NAME}:${CONNECTION_NAME}/data.csv")`)
  .then(t => grok.shell.addTableView(t[0]));
```

### Reading and writing files in file shares

The `grok.dapi` provides for fine-grained operations on files from file shares. It's possible to write files, read from
them, check existence, search for presence by a pattern, rename, move and delete. For example, here is how simple it is
to create a text file in your local file share:

`grok.dapi.files.writeAsText('<YOUR_NAME>:Home/testFile.txt', 'Hello, world!');`

All the `dapi.files` methods accept three types of inputs:

* a fully specified file path, as
  in [the examples](https://github.com/datagrok-ai/public/blob/master/packages/ApiSamples/scripts/dapi/files.js)
* a variable of type `file` (
  see [FileInfo](https://github.com/datagrok-ai/public/blob/14eb2acd6e36b33f64c4a0d108e940f7624af479/js-api/src/entities.js#L317))
  which may come, for example, from info panels working on files
* a string with a file share connection GUID

All the `dapi.files` methods are asynchronous.

This [example](https://github.com/datagrok-ai/public/blob/master/packages/ApiSamples/scripts/dapi/files.js)
provides a full understanding of the `files` API.

### Other ways for reading files

Finally, let's walk through other methods that can be used to open files from JavaScript:

* If you define a function that takes an input of `file` type (also
  see [FileInfo](https://github.com/datagrok-ai/public/blob/14eb2acd6e36b33f64c4a0d108e940f7624af479/js-api/src/entities.js#L317))
  , this creates a number of opportunities. First, you can call `file.readAsBytes()`
  or `file.readAsString()` methods on it. For instance, if you pass an obtained string
  to `grok.data.parseCsv(csv, options)`, it is possible to fine-tune the construction of dataframe from comma-separated
  values. Alternatively, you can pass a file to a [script](../../compute/scripting.md) to compute something and get the
  results back to your application's code.
* The method `grok.data.openTable(id)` might come in handy when you are reproducing a process and need to open a
  specific table by its ID. See an [example](https://public.datagrok.ai/js/samples/data-access/open-table-by-id).

See also:

* [JavaScript Development](../develop.md)
* [Data connection](../../access/data-connection.md)
* [Data query](../../access/data-query.md)
* [Data job](../../access/data-job.md)
* [File shares](../../access/connect-a-file-share.md)
* [Functions](../../datagrok/functions/function.md)
* [Parameterized queries](../../access/parameterized-queries.md)
