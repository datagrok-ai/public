<!-- TITLE: Importing data -->
<!-- SUBTITLE: -->

# Importing data

Datagrok lets you easily ingest virtually any kind of data from almost anywhere. It supports dozens of file types and
can connect to [30+ popular databases](data-connection.md#connectors) and other services, including cloud data
providers.

To import a local file, drag-and-drop it into the browser. The important thing to understand is that the data will stay
in the browser. It won't get sent to the server, unless you explicitly do that
by [uploading a project](../overview/project.md#uploading-a-project).

## Supported file types

| Extension     | Description              |
|---------------|--------------------------|
| csv, tsv, txt | Comma-separated file     |
| xml           | XML                      |
| json          | JSON                     |
| HTML          | HTML                     |
| xlsx          | Excel file               |
| edf           | European Data Format     |
| sas7bdat      | SAS                      |
| kml, kmz      | Geographic annotations   |
| rds, rda      | R Data Format            |
| h5            | Hierarchical Data Format |
| nc            | NetCDF                   |
| mat           | MATLAB MAT               |
| d42           | Datagrok project         |
| zip           | ZIP                      |
| gz, gzip      | gzip                     |
| tar           | tar                      |
| ipynb         | Jupyter Notebook         |

## Supported data sources

| Data Source                      |
|----------------------------------|
| Access                           |
| Amazon S3                        |
| Amazon Athena                    |
| Amazon Redshift                  |
| Google BigQuery                  |
| Google Cloud                     |
| Cassandra                        |
| DB2                              |
| Dropbox                          |
| File Network Shares              |
| Firebird                         |
| HBase                            |
| Hive                             |
| MS SQL                           |
| MariaDB                          |
| MongoDB                          |
| MySQL                            |
| ODATA                            |
| Oracle                           |
| PostgreSQL                       |
| SQLite                           |
| [Socrata](edit-socrata-query.md) |
| Sparql                           |
| Teradata                         |
| Twitter                          |
| Vertica                          |

## Videos

[![Notebooks](../uploads/youtube/data_access.png "Open on Youtube")](https://www.youtube.com/watch?v=dKrCk38A1m8&t=336s)

See also:

* [Import text](import-text.md)
* [Data query](data-query.md)
* [JS API: Load CSV](https://public.datagrok.ai/js/samples/data-access/load-csv)
