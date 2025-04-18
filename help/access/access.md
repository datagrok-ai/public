---
title: "Access"
sidebar_position: 0
---

Datagrok unifies data access across your organization, simplifying retrieving data from multiple sources. You can easily visualize, explore, and learn from your data, and use the insights gained to take action. Additionally, the platform provides access controls, security features, caching, and automatic monitoring of connection health.

## Data sources

In addition to local files that you can drag and drop from your computer, Datagrok integrates with various data providers. You can connect to any machine readable source: [a file storage](files/files.md) (like third-party cloud services or an organization-hosted Datagrok server), [databases](databases/databases.md), or [webservices](open-api.md).

Datagrok also supports scripting in various languages, such as R, Julia, and Python, which means you can create custom data sources. For example, you can [load a dataframe](https://public.datagrok.ai/js/samples/data-access/load-csv) from an external website or package, [open a specific table using its ID](https://public.datagrok.ai/js/samples/data-access/open-table-by-id), or write a package to extract data from multiple sources and combine them into one. For more information on getting data using functions and scripts, see [Access data](../develop/how-to/access-data.md#reading-files) section in the developers' documentation.

Datagrok also hosts [public datasets](public-datasets.md) that can be used for analysis, testing, and prototyping. These datasets cover various domains, including cheminformatics, clinical trials, and more.

### Data connection

A data connection is an [entity](../datagrok/concepts/objects.md) that contains the necessary information to connect to a specific data
source, such as its address and credentials. A data connection allows you to work with files and database tables
directly in Datagrok. You can access data either manually through the UI or programmatically
[through an application](../develop/how-to/access-data.md). For manual data access, Datagrok provides a convenient UI
that lets you connect directly to any of the [30+ supported connectors](databases/connectors/connectors.md), retrieve data using
queries, and securely share data with collaborators.

:::note

A connector is a plugin that enables the integration of external data providers into the platform. It can work with a
database, an Excel file, a CSV file, a web service, or any other source capable of providing the data. Most of our data
connectors are [open-sourced and extendable](https://github.com/datagrok-ai/public/tree/master/connectors) (under MIT
license).

:::

:::tip

To see all available data source connections, on the **Sidebar**, select **Manage** > **Connections**. From there, you can search connection by name or by tag.

<details>
<summary> You can also search and filter data connections using the following fields </summary>

| Field       | Description                      |
|-------------|----------------------------------|
| ID          |                                  |
| name        |                                  |
| server      |                                  |
| port        |                                  |
| db          |                                  |
| login       |                                  |
| dataSource  |                                  |
| description |                                  |
| createdOn   |                                  |
| updatedOn   |                                  |
| author      | [User](../govern/access-control/users-and-groups#users) object |
| starredBy   | [User](../govern/access-control/users-and-groups#users) object |
| commentedBy | [User](../govern/access-control/users-and-groups#users) object |
| usedBy      | [User](../govern/access-control/users-and-groups#users) object |

</details>

:::

A data connection is an [entity](../datagrok/concepts/objects.md), which means it can be shared, assigned permissions, annotated, and more.

For instructions on how to add a supported data source, set credentials, share, and manage it from the UI, see documentation for each data source type.

For instructions on how to add a supported data source, set credentials, share, and manage it programmatically, see [developer's documentation](../develop/how-to/access-data.md#connections).

For specific details on the configuration required, see each individual connector's documentation page in the Connectors directory.

### Data query

A data query is a [function](../datagrok/concepts/functions/functions.md) associated with a [data connection](#data-connection)
that typically returns a dataframe. Queries can be executed either manually, or as part of data jobs.
Datagrok has a convenient interface for creating, running, and sharing query results, including aggregation editor,
auto-generated parameter dialogs, and an ability to create dynamic dashboards to visualize query results. All data
governance features — such as data lineage, history, and security — are applicable to queries. For more information about
queries, see documentation for the respective data source type.

In addition to querying databases with SQL, you can query other data storages:

| Data source          | Query      |
|----------------------|------------|
| Relational databases | SQL        |
| File share           | Filename   |
| Excel file           | Sheet name |
| Linked data          | SPARQL     |
| Box                  | Filename   |

The result of a query execution is represented by the [function call](../datagrok/concepts/functions/function-call.md). For instructions on how to run, share queries, and manage queries programmatically, see [developer's documentation](../develop/how-to/access-data.md).

A data query is an [entity](../datagrok/concepts/objects.md), which means it can be shared, assigned permissions, annotated, and more.

## Data formats

With Datagrok, you can retrieve both structured and unstructured data. Datagrok supports [multiple data formats](files/supported-formats.md), including popular formats like CSV, TXT, JSON, and scientific formats like MAT, molecular structure formats (like PDB, MOL, or SDF), geographic annotation, and others. Datagrok also offers a flexible system for extending the platform with organization-specific data formats (see [Extensible framework](#extensible-framework)).

## Browsing and preview

Datagrok offers a rich set of features to help users efficiently browse, manage, and preview data. For more information, see:

* [Database Manager](databases/databases.md#database-manager)
* [File Manager](files/files.md#file-manager)
* [Webservices Manager](open-api.md#webservices-manager).

## Sharing and access control

Datagrok treats data connections, file shares, database tables and columns, and queries as [entities](../datagrok/concepts/objects.md), which means there is a common set of operations that can be applied to them. These entities can be shared with others, assigned access privileges, commented on, versioned, [audited](../govern/audit/audit.md), and so on. Some of the most popular privileges are: `view`, `edit`, `delete`, and `share`. These privileges can be given to individual users, or
to [groups](../govern/access-control/users-and-groups.md#groups). For more information on the access privilege model, see [Access Control](../govern/access-control/access-control.md).

Data connections can be shared as part of a [project](../datagrok/concepts/project/project.md), [package](../develop/develop.md#packages), or as a standalone [entity](../datagrok/concepts/objects.md). Sharing a query automatically shares the associated database connection, since query access depends on connection permissions. However, if you share a database connection with someone, your queries won't be shared automatically. You need to share them separately. 

To learn how to control access for each data source, see the documentation for the corresponding data source.

## Extensible framework

Datagrok is designed as an extensible environment, where extensions can customize or enhance virtually any part of the platform. For example, you can [create custom connectors](databases/create-custom-connectors.md), customize menus, add context actions, customize data preview, and more.

To learn more about extending and customizing Datagrok, see the [Develop](../develop/develop.md) section of our documentation. For examples related to data access, see [File Manager](files/files.md#file-manager).

## Resources

[![Data connection](../uploads/youtube/data_access.png "Open on Youtube")](https://www.youtube.com/watch?v=dKrCk38A1m8\&t=1048s)
