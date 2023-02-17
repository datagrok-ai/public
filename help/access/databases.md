# Databases

```mdx-code-block
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';
```

Datagrok connects to relational databases, allowing each database to function as a standalone data connection. Some of the features include:

* Inspecting schema, tables, and columns.
* Querying a database with the option to use parameterized queries for compatible connectors.
* Previewing and saving query results as a dataframe or as a dynamic dashboard.
* Sharing connections and queries with others.
* Configuring user privileges and managing data access.
* Creating data pipelines.

## Connecting to database

You can connect to a database via a direct connection. Out-of-the-box, Datagrok provides connectors to 30+ databases. To view information about a specific data source, see [the list of supported connectors](supported-connectors.md).

<details>

<summary> Developers </summary>

You can create [custom connectors](https://github.com/datagrok-ai/public/tree/master/connectors). <!--need to create a How-to under Develop-->

</details>

### Adding connection

![Add new database connection](database-connection.gif)

To add a new connection, follow these steps:

1. Go to **Data > Databases**.
2. In the **Database Manager**, right-click the desired connector and select **Add connection…** to open the **Add new connection** dialog.
3. Fill in the connection parameters.
4. Click **TEST** to test the connection, then click **OK** to save it. If the connection fails, verify your connection details and that you added Datagrok's IP addresses to your allowlist.

 > Notes:
 >
 >While many connection parameters are straightforward, some have unique characteristics.
 >
 >_JDBC connection_. For connections that support JDBC, you can use the _Conn. string_ parameter to enter a custom JDBC connection string. Leave all other parameter fields empty. You still need to enter credentials.
 >
 >_Credentials_. You have two options to enter credentials:
 >
 > * Manually (typically, login/password). When entered manually, Datagrok stores secrets in
   a [secure privilege management system](../govern/security.md/#credentials).
 > * Use the [Secrets Manager](data-connection-credentials.md/#secrets-managers).
 >
 >To specify who can change the connection credentials, click the **Gear** icon and then select from the **Credential owner** dropdown.

Upon successful connection, the database appears in the **Database Manager** under the respective data source. By expanding the database, you can view its saved queries. [If connectors support it](supported-connectors.md), you can also inspect the schemas, tables, and columns of relational databases.

>Note: Like other objects in Datagrok, newly created connections are only visible to the user who created them. To let others access the connection, you must [share it](databases.md/#sharing-and-managing-connections).

### Caching data

When working with database connections and queries, you can cache data to improve performance, especially with slow servers or connections. One common use case for caching is storing values used to build the user interface automatically, typically in the form of a `select distinct <name> from <table>`.

Datagrok provides three caching options:

```mdx-code-block
<Tabs>
<TabItem value="schema" label="Cache schema" default>
```

Choosing this option will cache all tables and table columns in the database connection. To cache schema:

1. Right-click a desired connection and select **Edit** from the list.
1. In the dialog that opens, select the checkbox **Cache Schema**.
1. Click **OK** to save the changes.

```mdx-code-block
</TabItem>

<TabItem value="connection" label="Cache connection">
```

Choosing this option will cache the results for all queries executed on this connection. To cache all queries:

1. Right-click a desired connection and select **Edit** from the list.
1. In the dialog that opens, select the checkbox **Cache Results**.
1. In the **Invalidate On** field that appears, type in a [cron expression](https://www.freeformatter.com/cron-expression-generator-quartz.html) that defines when the cache is invalidated.
    <details>

    <summary> Example </summary>

    Here is an example of a cron expression that invalidates the cache at 1 AM every night for a database that refreshes overnight: `0 0 1 * * ?`

    </details>

    >**IMPORTANT**: If you leave the **Invalidate On** field empty, the cache won't be reloaded.
1. Click **OK** to save the changes.

```mdx-code-block
</TabItem>

<TabItem value="query" label="Cache query">
```

Choosing this option will cache the results of an individual query. This is the most specific way to enable caching in Datagrok. To cache a specific query:

1. Right-click a query and select **Edit** from the list.
1. In the dialog that opens, select the checkbox **Cache Results**.
1. In the **Invalidate On** field that appears, type in a [cron expression](https://www.freeformatter.com/cron-expression-generator-quartz.html) that defines when the cache is invalidated.
    <details>

    <summary> Example </summary>

    Here is an example of a cron expression that invalidates the cache at 1 AM every night for a database that refreshes overnight: `0 0 1 * * ?`

    </details>

    >**IMPORTANT**: If you leave the **Invalidate On** field empty, the cache won't be reloaded.
1. Click **OK** to save the changes.

```mdx-code-block
</TabItem>

</Tabs>
```

To set caching for connections that are created automatically as part of a package, you can specify the `cacheResults` parameter in the connection json definition:

```json
{
  "name": "Northwind",
  "parameters": {
    "server": "dev.datagrok.ai",
    "port": 23306,
    "db": "Northwind",
    "cacheSchema": false,
    "cacheResults": true,
    "ssl": false,
    "connString": ""
  }
}
```

To set caching for query results programmatically, you need to edit the corresponding .sql file in the package queries and specify the `meta.cache` and `meta.invalidate` fields. These fields determine whether and how to cache the query results:

* `meta.cache` specifies whether the query is cached or not.
* `meta.invalidate` specifies when the cached result should be invalidated.

By setting these fields, you can control the caching behavior of individual queries and improve overall query performance.

<details>

<summary> Example </summary>

```Grok Script
--name: getProductNames
--input: string department
--meta.cache: true
--meta.invalidate: 0 0 1 * * ?
select distinct name from products p
where p.department = @department
```

</details>

### Modifying connection

1. Right-click a connection and select **Edit...** A dialog opens.
2. In the dialog, change the connection name, parameters, or credentials as needed.
3. Click **TEST** to test the connection, then click **OK** to save the changes.

> Tip: To quickly create a connection similar to an existing one, clone a connection and make changes:
>
> 1. Right-click a connection and select **Clone...** A dialog opens.
> 2. In the dialog, type in a name for the new connection and make other changes as needed.
> 3. Re-enter password or access keys.
> 4. Click **OK** to save the new connection.

## Browsing data

### Database Manager

**Database Manager** is a [hierarchical viewer](../datagrok/workspace.md) that lets you browse and manage database objects. Right-click any object in the **Database Manager** to access a list of available actions. Depending on your permissions, you can edit, delete, share, and perform other actions.

Each object in Datagrok has multiple contexts, or _data properties_. Users get access to these contexts through [info panels](../discover/info-panels.md), which are displayed in the **Context Panel**. The [**Context Panel**](../datagrok/navigation.md#properties) is dynamic and context-sensitive, meaning it adapts to the selected object, presenting only the relevant information and options.

<details>

<summary> Developers </summary>

**Context Panel** can be extended. You can add custom [info panels](../develop/how-to/add-info-panel.md) and [context actions](../develop/how-to/context-actions.md).

</details>

**Database Manager** and **Context Panel** work in tandem. For example, when you click on a table in the **Database Manager**, the **Context Panel** updates, allowing you to view the table's metadata, dynamically preview the table's contents, run queries, and access other relevant information and options.

![DB Hierarchy Browser](../uploads/features/db-hierarchy-browser.gif "DB Hierarchy Browser")

> Tip: You can reposition, resize, detach, or hide any panel. To learn more about customizing your workspace, see [Navigation](../datagrok/navigation.md).

### Schema Browser

**Schema Browser** is a graphical representation of a database schema. It provides a quick and easy way to view and access all tables and their corresponding columns. To open a **Schema Browser**, right-click on the desired connection or table in the **Database Manager** and select **Browse schema**.

Like the **Database Manager**, the **Schema Browser** is designed to work with the **Context Panel**. Use the **Context Panel** to access information and options for the selected object, such as preview table content or run a query.

![Schema browser](browse-schema.gif)

>Tip: The context menu of an object provides quick access to the list of available actions. This functionality is available throughout the platform. For example, to run a query, right-click on the desired table either in the **Schema Browser** or the **Database Manager** and select from the context menu. You can also locate the object name at the top of the **Context Pane** and click the **Down Arrow** control next to it.
>
>If you don't see a certain action, it may be due to insufficient permissions. Contact your Datagrok administrator for assistance.

### Viewing schema as dataframe

Another way to explore a database schema is through an interactive spreadsheet. To access this feature, in the **Database Manager**, right-click on the database connection and select **Open schema as table**. This will open the dataframe in Datagrok, allowing you to interact with and explore the database schema in a tabular format.

## Querying and filtering data

Datagrok provides several tools for querying and filtering data. Choose from the options below depending on the type of your query and the results you want:

* [Query editors](#query-editors): These tools let you manually create or edit a query.
* [Built-in table queries](#built-in-table-queries): These queries return either the first 100 rows, or all rows in a table.
* [Table joiner](#joining-tables): This tool lets you merge multiple tables within a database connection.
* [Filters](#filtering)

> Note: The ability to query tables is connetor-specific and may not be available for all.

When you successfully create a query, it appears in the **Database Manager** under the corresponding connection or table. From there, you can run it, share it with others, manage access, and perform other tasks. To edit an existing query or create a copy of it, locate the query in the **Database Manager** and right-click to choose either **Edit...** or **Clone...**

> Tip: To quickly find the exact query you're looking for, go to **Data** > **Queries** and type in the query's name or tag. To see all querities belonging to a connection, right-click it and select **Browse queries** from the list of options.

### Query editors

Datagrok offers two editors for querying databases:

* [**SQL Query**](#sql-query-editor): This editor is the main interface for executing database queries. Using this tool, you can write and edit SQL statements, add post-processing steps, and preview results in real-time.

   To open, right-click a database connection or a table and select **Add query** (for connections) or **New SQL Query...** (for tables) from the list of options.
  
* [**Aggregation Query**](#aggregation-query-editor): This editor works with tables. Use it to manipulate, summarize, filter, and pivot table data.

   To open, in the **Database Manager** right-click a table and select **Aggregate table data**.

#### SQL Query editor

The **SQL Query** editor has two tabs:

```mdx-code-block
<Tabs>
<TabItem value="query" label="**Query Tab**" default>
```

This tab is where you write and edit SQL statements.

As you work on your query, you can preview the query output any time by clicking the **Run/Refresh** button on the menu ribbon or by pressing the F5 key on the keyboard. The the query result is displayed as an interactive dataframe. You can scroll to view object details, perform actions on columns, and more.

When your query is complete, give it a name and click the **Save** button. If you don't want to save the query, close the editor without saving.

![Create a database query](query-add.gif)

>Note: You can also add the query results to the workspace for further analysis. To do so, click the **Dropdown Arrow** control in the bottom left corner of the **Query** tab and selecting **Add results to workspace**.

```mdx-code-block
</TabItem>
<TabItem value="transformations" label="**Transformations Tab**">
```

Use this tab to apply post-processing operations to your query:

![Query transformations](transformations.gif)

1. Select a function to apply to a query. Use checkboxes next to the operation categories to filter.
1. Set the function parameters in the dialog that opens.
1. Click **OK**.
1. Repeat the steps as needed.

Each transformation is recorded. To see the results of each record, locate the record on the left and click it. The data output for this step is displayed in the **Preview**. You can also edit or remove transformation records as needed. To edit a transformation step, click the **Dropdown Arrow** control next to the transformation record and select from the list of options. To remove the transformation record entirely, click the **Delete** icon.

> Tip: Use the menu ribbon on top to add or remove columns quickly.

<details>

<summary> Developers </summary>

You can create custom transformation functions in R, Python, or any other language. See [Scripting](../compute/scripting.md).

</details>

```mdx-code-block
</TabItem>
</Tabs>
```

#### Aggregation Query editor

Similar to the **SQL Query** editor, the **Aggregation Query** editor has two tabs (**Queries** and **Transformations**), and the **Query Tab** in both editors looks and works the same.

To aggregate data from a table, use the **Transformation** tab. In this tab, you can choose which columns to include in your report and decide how to pivot and group them:

![Visual query](databases-visual-query.gif)

1. First, define the columns that stay the same and add them to the **Rows** section. These columns serve as keys, and their unique values become row identifiers.

1. To group and handle non-unique data, set the aggregation parameters:
   * Select the aggregation column. To do so, click the **Add...** icon next to **Measures** and select from the list.
   * Choose the aggregating function. To do so, clicking the **Add...** icon again, then click **Aggregation...** to select from the list.

      >Note: The list of functions may vary based on the connector and any packages that have been installed. If connectors expose custom functions, such as GIS functions in PostreSQL, Datagrok makes them readily available in the interface.

      <details>

      <summary> Developers </summary>

      You can [create custom aggregation functions](../compute/scripting.md).
      </details>
   * (Optional) Repeat these steps to add additional parameters.

     Each parameter takes its own aggregating function. To modify or remove a parameter, right-click and select the action from the list.

1. To pivot unique values from one column into multiple columns in the query output, add the corresponding columns to **Columns**.

1. (Optional) Add **Filters** to specify which items are returned when you run a query.

### Built-in table queries

Datagrok offers two pre-built table queries available from the context menu:

* **Get All**: This query retrieves all data. Use it with caution.
* **Get TOP 100**: This query retrieves the first 100 rows.

  >Tip: To retrieve specific columns, in the **Database Manager** hold down the Shift key on your keyboard while clicking on the desired columns in the schema. Once selected, use the context actions to run the query just for these columns. ![Get Columns](db-exploration-get-columns.png)

### Joining tables

You can merge data from multiple tables in a schema:

1. Right-click a table and select **Join tables** to open a dialog with a list of columns connected by keys.
1. To choose columns, check the corresponding checkboxes. Datagrok automatically generates the SQL statement based on your choices. You can view the result in the **Preview** section of the dialog.
1. (Optional) Manually edit the SQL statement in the dialog. As you make changes, **Preview** updates to reflect your edits.
1. When you're done, use the **Dropdown Arrow** control in the bottom left corner of the dialog to choose between the two options:

   * **Save as query**: This opens the **SQL Query** editor where you can further edit the query or add additional processing steps before saving.
   * **Add result to workspace**: This opens the query output as a dataframe, which you can inspect and edit further.

![Create a join query](query-builder.gif)

### Filtering

You can use these fields to filter queries with [smart search](../datagrok/smart-search.md):

| Field       | Description                                                       |
|-------------|-------------------------------------------------------------------|
| ID          |                                                                   |
| name        |                                                                   |
| query       |                                                                   |
| runs        | list of [FuncCall](../datagrok/functions/function-call.md) object |
| connection  |                                                                   |
| jobs        | [DataJob](data-job.md) object                                     |
| createdOn   |                                                                   |
| updatedOn   |                                                                   |
| author      | [User](../govern/user.md) object                                  |
| starredBy   | [User](../govern/user.md) object                                  |
| commentedBy | [User](../govern/user.md) object                                  |
| usedBy      | [User](../govern/user.md) object                                  |

## Parameterized queries

A parameterized query uses _parameters_ as placeholders for constant values. The parameters are assigned values only during query execution, so you can run the same query with different values instead of creating distinct queries for each scenario. Parameterized queries are also useful for dynamic data where the values are not determined until the statement is executed. All parameters are optional and are defined at the beginning of a query.

### Adding parameter

The syntax for defining query parameters is based on [scripting](../compute/scripting.md) with additions specific to queries.

To add a parameter, follow these steps:

1. Right-click the query to which you wish to add a parameter. Select **Edit** to open the **SQL Query** editor.
1. In the **Query** tab, annotate parameters. Use `--` for SQL and `#` for Sparql followed by the parameter type.
    >Note: If you declare a parameter without a type, the type is assumed to be `string`.
1. (Optional) Add a default value. The table below lists supported types.
   >Note: If you declare a parameter without a default value, a value must be entered at execution.

   | Parameter  | Description, supported type | Input template |
   |----------------|----------|---------|
   | `name`         |       |      |
   | `friendlyName` | Friendly name  |       |
   | `description`  | Description           |
   | `help`         | Help URL    |          |
   | `tags`         | Tags           |       |
   | `input`        | `int` - integer number<br />`double`  - float number <br />`bool` - boolean<br />`string` - string<br />`datetime`- DateTime<br />`list<T>` - a list of type T | `--input: <type> <name> = <value> {<option>: <value>; ...} [<description>]` |
   | `output`       | Output parameter        |
1. When you're done, click **Save**.

<details>

<summary> Example </summary>

>```sql
>--input: string productName                                               
>select * from products where productname = @productName
>```

</details>

To define input parameters, you have these options:

1. Enter a single value.
1. Use parameter _functions_:
   * [Choices](../compute/scripting.md/#parameter-choices): These functions return a list of strings with no parameters.
      <details>

      <summary> Example </summary>

      ```sql
      --input: string shipCountry = "France" {choices: ['France', 'Italy', 'Germany']}
      ```

      </details>
   * [Editors](../compute/scripting.md/#parameter-editors): These functions generate parameters from the output of another function.
      <details>

      <summary> Example </summary>

      ```sql
      --input: string shipCountry = "France" {choices: Query("SELECT DISTINCT shipCountry FROM Orders")}
      ```

      </details>
   * [Suggestions](../compute/scripting.md/#parameter-suggestions): These functions take a string argument and return a list of matching strings.
      <details>
      <summary> Example </summary>

      ```sql
      --input: string shipCountry = "France" {suggestions: Demo:northwind:countries}
      ```

      </details>
   * [Validators](../compute/scripting.md/#parameter-validators): These functions take a value and validate it against specified criteria. If validation fails, validators return an object and an error message (optional).

1. Reuse other parameters as parameter values.
    <details>

    <summary> Example </summary>

     ```sql
     --input: string firstLetter = "F"
     --input: string shipCountry = "France" {choices: Query("SELECT DISTINCT shipCountry FROM Orders WHERE shipCountry LIKE @firstLetter || '%')}
     SELECT * FROM Orders WHERE (shipCountry = @shipCountry)
     ```

     </details>
1. Use parameter _patterns_ as filters to allow users to enter free text (supported for the `string` and `datetime` only).
    >Note: A _parameter pattern_ is an expression to create a filter condition within the query. Much like a formula, it consists of field references, operators, and constants. To see a list of available patterns, see [Search patterns](../explore/data-search-patterns.md).

    <details>

    <summary> Example </summary>

    Add filtering criteria for the column _freightValue_:

    ```sql
    --input: string freightValue = >= 10.0 {pattern: double}
    select * from Orders where @freightValue(freight)
    ```

    Add filtering criteria for the column _orderDate_:

    ```sql
    --input: string orderDate = "after 1/1/1995" {pattern: datetime}
    select * from orders where @orderDate(orderDate)
    ```

    </details>

> Note: You can persist the collected parameters to use them with more than one query.

## Running queries

To run a query, you have these options:

* Run a query from within a [query editor](#query-editors). Use this option to preview the query results when writing or modifying a query.
* Use the **Run** _context action_ in all other cases. This opens the query output as a dataframe, [giving you several options to handle it](#query-results).

When you run a parameterized query, the system automatically generates a dialog for entering parameters. For each input parameter, there is an input control that prompts you to enter values (for example, a _calendar_ for dates, _checkbox_ for boolean, or a _button_ that launches a sketcher application for filtering substructures.) Datagrok automatically parses the entry
and executes a parameterized, secure, connector-specific SQL query on the backend.

<details>

<summary> Example </summary>

```sql
--input: int employeeId = 5
--input: string shipVia = = 3 {pattern: int}
--input: double freight = 10.0
--input: string shipCountry = "France" {choices: Query("SELECT DISTINCT shipCountry FROM Orders")}
--input: string shipCity = "starts with r" {pattern: string}
--input: bool freightLess1000 = true
--input: datetime requiredDate = "1/1/1995"
--input: string orderDate = "after 1/1/1995" {pattern: datetime}
SELECT * FROM Orders WHERE (employeeId = @employeeId)
    AND (freight >= @freight)
    AND @shipVia(shipVia)
    AND ((freight < 1000) OR NOT @freightLess1000)
    AND (shipCountry = @shipCountry)
    AND @shipCity(shipCity)
    AND @orderDate(orderDate)
    AND (requiredDate >= @requiredDate)
```

![Run a parameterized query](parameterized-query.gif)

</details>

<details>

<summary> Developers </summary>

You can expose the parameter dialog to end-users as an [application](../develop/how-to/create-package.md).<!--Mention?: when the cartdridge is not deployed on that particular database, the query returns an error-->

Similar to functions in JavaScript, queries are [functions](../datagrok/functions/functions.md) in Datagrok. This means you can “call” a query across multiple Datagrok plugins instead of copy-pasting the same query every time.

To run a query programmatically, see this [code snippet](https://public.datagrok.ai/js/samples/data-access/parameterized-query).

</details>

## Query results

When you run a query, you receive the results in a dataframe.

>Tip: To return a value of a different data type, specify it in the output parameter.  
>
><details>
>
><summary> Example </summary>
>
>Return a string with the semantic type `Molecule`:
>
> ```sql
> --output: string smiles {semType: Molecule}
> ```
>
></details>

In the dataframe, you can perform various operations such as data cleansing, transformation, [and more](/visualize/viewers/grid.md). Once you have the data set up to your satisfaction, you can add viewers and create a visualization. You have the option to [persist the layout as a template](/visualize/view-layout.md) for future use.

![Open schema as dataframe](open-schema-as-table.gif)

### Sharing query results

You have two options to share query results:

* [Share the dataframe's URL](#sharing-query-results-as-url)
* [Upload the dataframe to the server](/datagrok/project.md#uploading-a-project) to share it as a project. For parameterized queries, you can create [dynamic dashboards](#using-parameterized-queries-to-create-dynamic-dashboards) that allow users to view different datasets and interact with data in real-time by changing the parameters.

#### Sharing query results as URL

Each query output has a unique URL that encodes the parameters entered. This URL can be used to share the query results. When you change query parameters, the corresponding URL changes.

After you have executed a query, copy the URL from the address bar and share it with others. To access the query results from the link provided, users must have the necessary permissions to execute this query.

#### Using parameterized queries to create dynamic dashboards

To convert a parameterized query output into a dynamic dashboard, do the following:

1. Go to **Data** > **Databases**.
2. In the **Database Manager**, right-click the desired query and select **Run**. A parameter dialog opens.
3. Fill in the query parameters and click **OK** to open a dataframe.
4. (Optional) Customize the dashboard by adding viewers, filters, or applying [layouts](../visualize/view-layout.md).
5. To save the dashboard, you need to upload the project to the server:
   1. On the **Sidebar**, click **Projects** > **Upload**. This opens a dialog.
   2. In the dialog, enter a name and description for the project (optional) in the fields provided.
   3. Select how to store data: (1) save the data as a static snapshot or (2) store the data as a
      generation script by toggling the **Data sync** control. The second option reruns the query each time the project is opened.
      > Note: For more information on dynamic data updates in projects, see [Dynamic data](../datagrok/project.md/#dynamic-data).
6. Click **OK** to upload the project.

After uploading the project, it can be opened or shared with others. To open a project, go to **Data** > **Projects**. Right-click the project to access its context menu or double-click to open it.

![Dynamic dashboards](dynamic-dashboards.gif)

## Access control

In Datagrok, certain classes of objects we call _entities_ have a [standard set of operations](../datagrok/objects.md) that can be applied to them. This includes connections, queries, tables, and table columns, all of which can be shared, assigned permissions, annotated, and more.

When you create a connection or a query, it's private and visible to you only. You can choose to share it, making it accessible to others.

>Note: If you can't share a connection or a query, you may have insufficient permissions. Contact your Datagrok administrator for assistance.

To share:

1. Right-click the connection or a query and select **Share...**. This opens a dialog.
2. Enter a user or a group that needs access and set corresponding permissions:
   * can_create
   * can_edit
   * can_delete
   * can_query

   These privileges can be given to individuals or to [groups](../govern/group.md) (which can be defined via dynamic filters). For more information on the access privilege model, see [Privileges](../govern/security.md#privileges).
3. (Optional) Add a description in the provided text field. If you don't want to notify the recipients, clear the **Send notification** checkbox.
   > Note: To notify via email, enter the user's email in the identity/email field. The email contains a link to the shared item and entered description. When you enter a user or a group name, they are  notified via the Datagrok interface.
4. Click **OK** to share.

![Share a database connection](sharing-database-connections.gif)

> Note: Datagrok query belongs to the database connection for which it's created. It means you can’t share a query without sharing a connection. Deleting a connection also deletes a query.

Subject to your priveleges, you can use the **Context Pane** on the left to inspect and quickly adjust access permissions to your databases, manage queries, view history and activity details, send comments to those you're sharing with, and more.

> Tip: You can access the same list of actions from the context menu. The context menu displays only actions available to you based on your permissions.

## Data reproducibility

Any action performed on datagrok entities is reproducable and can be used in
automation workflows. For example, you can use data preparation pipeline to define jobs for data ingestion, postprocessing, and transformations.

To learn more about automating workflows usung data preparation pipelines, see [Data preparation pipeline](data-pipeline.md).

## Resources

* Videos
  * [Database exploration](https://www.youtube.com/watch?v=YJmSvh3_uCM)
  * [Parameterized queries - Overview](https://www.youtube.com/watch?v=dKrCk38A1m8&t=1980s)
  * [Parameterized queries - Example](https://www.youtube.com/watch?v=sSJp5CXcYKQ)
  * [Using lists in parameterized queries](https://www.youtube.com/watch?v=meRAEF7ogtw)
* Tutorials
  * [Adding parameters to functions](../datagrok/functions/func-params-enhancement.md)
