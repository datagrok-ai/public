---
title: "Data preparation pipeline"
---

DPP is a core component of the Datagrok platform designed to let end users define jobs that would get data from
disparate data sources, clean or merge the data if needed, run transformations, build interactive dashboards based on
the retrieved data, and publish these dashboards. On top of that, it provides users with a web interface for monitoring
and managing of the jobs. That includes creating, editing, deleting, scheduling, browsing log files, as well as the
ability to cancel running jobs.

## Concepts

Datagrok offers [30+ connectors](../databases/connectors/connectors.md) to different data sources. They include all popular
databases, as well as other queryable data sources, such as Twitter or Facebook. To view a list of data sources
currently available to you, run **File | Connect to Data...**.

Before querying data, a [connection](../access.md#data-connection) has to be created. To do that, right-click on a connection and
choose **Add connection**. At this step, you should provide connection parameters. Fields are specific to the data
source, typically you would need to specify server's address and login credentials. Keep in mind that the credentials
are stored on the server in a secure way, and can be shared with others using our privilege management system. To see
all connection at once, open **Admin | Data Connections**

Now we are ready to create a data query. To do that, choose **Add Query** from the connection's context menu. Query
editors are specific to data sources. For instance, for relational databases you would need to enter a SQL query, for
linked data that would be a SPARQL query, and for Twitter it would be something completely different. Try running the
query before saving it to make sure it works as intended. It is also possible to create
[queries that accept parameters](../access/databases/databases.mdx#parameterized-queries), which will be either explicitly
provided by user, or provided by the parent data job. To browse all queries, open `Admin | Data Queries`.

To run a query, double-click on it, or select **Run** from the context menu. To see information associated with the
query, click on it to make it a current object, and then expand panels in the context panel on the right. The same
technique applies to other objects in the platform.

Running a query yields a table. To incorporate custom processing into the query, click **Edit**
under the 'Transformations' panel in the context panel. In addition to the basic data transformation routines, advanced
functions are available as well. Examples of such transformations are R or Python scripts, applying predictive models,
etc.

Sometimes, you want more than a table to be returned. Data jobs let you get results from multiple queries
at once, massage the data using transformations, and apply any visualizations on top of it. Just as data queries, data
jobs can be parameterized as well. The output of the data job is a [project](../datagrok/concepts/project/project.md), which is
essentially a dashboard.

## Security

We take the issue of security and access privileges very seriously, and the concepts of access control and users' roles
and privileges are at the heart at the system. Each system entity (such as `Data Source`, `Data Connection`
, `Data Query`, or `Data Job`)
has a list of possible actions associated with it, such as `view`, `edit`, or `publish`. The platform has a flexible
access control mechanism that lets us define people (or groups of people)
that are allowed to execute actions against different entities, based on the entity attributes. For instance, it is
possible to define a group of people who would be able to open dashboards, but would not have access to the underlying
connection. See
[Access Control](../access.md#sharing-and-access-control) for details.

See also:

* [Data connection](../access.md#data-connection)
* [Data query](access.md#data-query)
* [Function call](../datagrok/functions/function-call.md)
