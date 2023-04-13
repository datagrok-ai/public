<!-- TITLE: Access -->
<!-- SUBTITLE: -->

Datagrok provides a single, unified access point for data accross an organization, enforcing data restrictions throughout the platform. With Datagrok, you can centralize data you collect from disparate sources to visualize, explore, learn from it, and use the insights gained to take action.

## Data sources

Besides local files that you can drag and drop from your computer, Datagrok integrates with various data providers. You can connect to any machine readable source: [a remote file storage](file-shares.md) (like third-party cloud services or an organization-hosted Datagrok server), relational and non-relational [databases](databases.md), or [webservices](open-api.md).

Datagrok also supports scripting in various languages, such as R, Julia, and Python, which means you can create custom data sources. For example, you can [load a dataframe](https://public.datagrok.ai/js/samples/data-access/load-csv) from an external website or package, [open a specific table using its ID](https://public.datagrok.ai/js/samples/data-access/open-table-by-id), or write a package to extract data from multiple sources and combine them into one. For more information on getting data using functions and scripts, see [Access data](../develop/how-to/access-data.md/#reading-files) section in the developers' documentation.

Datagrok also hosts [public datasets](public-datasets.md) that can be used for analysis, testing, and prototyping. These datasets cover various domains, including cheminformatics, clinical trials, and more.

### Data connection

A data connection allows you to work with files and database tables directly in Datagrok. When connecting to a data source, you can access data manually from the UI, or programmatically, [through an application](../develop/how-to/access-data.md). For manual data access, Datagrok provides a convenient UI that lets you connect directly to any of the [30+ supported connectors](supported-connectors.md), retrieve data using queries, and securely share data with others.

:::note

A connector can work with a database, an Excel file, a CSV file, a web service, or any other source capable
of providing the data. Most of our data connectors are [open-sourced and extendable](https://github.com/datagrok-ai/public/tree/master/connectors) (under MIT license).

:::

:::tip

To see all available data source connections at once, on the **Sidebar**, select **Manage** > **Connections**. From there, you can search connection by name or by tag.

<details>
<summary> You can also search or filter data connections using these fields </summary>

| Field       | Description                                 |
|-------------|---------------------------------------------|
| ID          |                                             |
| name        |                                             |
| server      |                                             |
| port        |                                             |
| db          |                                             |
| login       |                                             |
| dataSource  |                                             |
| description |                                             |
| createdOn   |                                             |
| updatedOn   |                                             |
| author      | [User](../govern/user.md) object         |
| starredBy   | [User](../govern/user.md) object         |
| commentedBy | [User](../govern/user.md) object         |
| usedBy      | [User](../govern/user.md) object         |

</details>

:::

A data connection is an [entity](../datagrok/objects.md), which means it can be shared, assigned permissions, annotated, and more.

For instructions on how to add a supported data source, set credentials, share, and manage it from the UI, see documentation for each data source type.

For instructions on how to add a supported data source, set credentials, share, and manage it prorammatically, see [developer's documentation](../develop/how-to/access-data.md/#connections).

For specific details on the configuration required, see each individual connector's documentation page in the Connectors directory.

### Data query

Datagrok has a convenient interface for creating, running, and sharing query results, including visual query editors, auto-generated parameter dialogs, and an ability to create dynamic dashboards to visualize query results. All data governance features, such as data lineage, history, and security, are applicable to queries. For more information about queries, see documentation for the respective data source type.

Typically, a query is used against a database, however the same concepts apply for other data sources that are listed
below:

| Data source          | Query      |
|----------------------|------------|
| Relational databases | SQL        |
| File share           | Filename   |
| Excel file           | Sheet name |
| Linked data          | SPARQL     |
| Box                  | Filename   |

The result of a query execution is represented by the [function call](../datagrok/functions/function-call.md). For instructions on how to run, share queries, and manage queries prorammatically, see [developer's documentation](../develop/how-to/access-data.md).

A data query is an [entity](../datagrok/objects.md), which means it can be shared, assigned permissions, annotated, and more.

## Data formats

With Datagrok, you can retrieve both structured and unstructured data. Datagrok supports [multiple data formats](supported-formats.md), including popular formats like CSV, TXT, JSON, and scientific formats like MAT, molecular structure formats (like PDB, MOL, or SDF), geographic annotation, and others. Datagrok also offers a flexible system for extending the platform with organization-specific data formats (see [Extensible framework](#extensible-framework)).

## Browsing and preview

Datagrok offers an array of capabilities and features designed to help users efficiently browse, manage, and preview the content of their data. For more information, see:

* [Database Manager](databases.md/#database-manager)
* [File Manager](file-shares.md/#file-manager)
* [Webservices Manager](open-api.md/#webservices-manager).

## Sharing and access control

Datagrok treats data connections, file shares, database tables and columns, and queries as [entities](../datagrok/objects.md), which means there is a comon set of operations that can be applied to them. These entities can be shared with others, assigned access privileges, commented on, versioned, [audited](../govern/audit.md), and so on. Some of the most popular privileges are: `view`, `edit`, `delete`, and `share`. These privileges can be given to individual users, or
to [groups](../govern/group.md). For more information on the access privilege model, see [Privileges](../govern/security.md#privileges).

Data connections can be shared as part of a [project](../datagrok/project.md), [package](../develop/develop.md#packages) (and [repository](connectors/git.md) containing this package), or as a standalone [entity](../datagrok/objects.md). The access rights of a database connection are inherited from the access rights of a query. However, the reverse is not true: access rights of a query don't inherit the access rights of the database connection. Consequently, when sharing a query, the associated database connection is shared automatically. On the other hand, sharing a database connection does not automatically share your queries. For web queries, they are shared automatically when the corresponding connection is shared.

To learn how to control access for each data source, see the documentation for the corresponding data source.

## Extensible framework

We designed Datagrok as an extensible environment, where extensions can customize or enhance any part
of the platform. For example, you can [create custom connectors](create-custom-connectors.md), add organization-specific data formats, customize menus, add context actions, customize data preview, and more.

To learn more about extending and customizing Datagrok, see the [Develop](../develop/) section of our documentation. For examples related to data access, see [File Manager](file-shares.md/#file-manager).

## Resources

[![Data connection](../uploads/youtube/data_access.png "Open on Youtube")](https://www.youtube.com/watch?v=dKrCk38A1m8\&t=1048s)
