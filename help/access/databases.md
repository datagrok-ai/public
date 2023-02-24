# Databases

```mdx-code-block
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';
```

Datagrok has native support for relational databases ([30+ databases supported](supported-connectors.md)), allowing to do the following:

* Browse databases schemas, tables, and columns.
* Query databases using parameterized queries.
* Use query as a dynamically refreshed data source for dashboards.
* Share connections and queries with others.
* Configure user privileges.
* Create data pipelines.

:::note Developers

You can create [custom connectors](https://github.com/datagrok-ai/public/tree/master/connectors).

:::

## Connecting to database

### Adding connection

![Add new database connection](database-connection.gif)

To add a new connection, follow these steps:

1. Go to **Data > Databases**.
2. In the **Database Manager**, right-click the desired connector and select **Add connection…** to open the **Add new connection** dialog.
3. Fill in the connection parameters.
4. Click **TEST** to test the connection, then click **OK** to save it. If the connection fails, verify your connection details and that you added Datagrok's IP addresses to your allowlist.

:::note

The connection originates from the Datagrok server, so make sure your database is accessible from there.

:::

While many connection parameters are straightforward, some have unique characteristics:

* _JDBC connection_. For connections that support JDBC, you can use the _Conn. string_ parameter to enter a custom JDBC connection string. Leave all other parameter fields empty. You still need to enter credentials.
* _Credentials_. You have two ways to specify credentials:
  * Manually. When entered manually, Datagrok stores secrets in
   a [secure privilege management system](../govern/security.md/#credentials). To specify who can change the connection credentials, click the **Gear** icon and select from the **Credential owner** dropdown.

  * Use the [Secrets Manager](data-connection-credentials.md/#secrets-managers), such as the AWS Secrets Manager.

Upon successful connection, the database appears in the **Database Manager** under the respective data source. By expanding the database, you can view its saved queries. [If connectors support it](supported-connectors.md), you can also inspect the schemas, tables, and columns of relational databases.

:::note

Like other objects in Datagrok, newly created connections are only visible to the user who created them. To let others access the connection, you must share it (right-click the connection and select **Share...** from the list of options).

:::

### Caching data

You can cache query results to improve query performance. To change caching options, right-click the connection, select **Edit...**, and the select the appropriate checkboxes in the **Edit** dialog:

1. **Cache Schema**: to speed up building the database schema.
1. **Cache Results**: to reuse results for all queries under this connection. To cache results of individual queries, edit that query's properties.

When you cache the query results, use the **Invalidate On** field to specify a [cron expression](https://www.freeformatter.com/cron-expression-generator-quartz.html) that defines when the cache is invalidated. Here is an example of a cron expression that invalidates the cache at 1 AM every night for a database that refreshes overnight: `0 0 1 * * ?`

:::caution

If you leave the **Invalidate On** field empty, the cache won't ever be updated.
This only makes sense when the data in the database never changes. 

:::

<details>

<summary> Caching package-defined connections </summary>

To cache package-defined connections, you can specify the `cacheResults` parameter in the connection json definition:

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

</details>

<details>

<summary> Setting cache options programmatically </summary>

To cache the query results programmatically, edit the corresponding .sql file in the package queries and specify the `meta.cache` and `meta.invalidate` fields. These fields determine whether and how to cache the query results:

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

</details>

### Modifying connection

To modify a connection, right-click it and select **Edit...** from the list of options. To quickly create a connection similar to an existing one, right-click it and select **Clone...**

## Database Manager

**Database Manager** lets you manage database connections and browse database objects. To access an object's context actions, right-click it. Alternatively, left-click the object and expand the **Actions** pane in the **Context Panel** on the left.

If you don't see a certain action, it may be due to insufficient permissions. Contact your Datagrok administrator for assistance.

:::note

When you click on an object, its properties and context actions are shown in the [**Context Panel**](../datagrok/navigation.md#context-panel) on the right. For example, when you click a table in the **Database Manager**, the **Context Panel** updates, allowing you to view the table's metadata, dynamically preview the table's contents, run queries, and access other relevant information and options.

<details>
<summary> Developers </summary>

**Context Panel** can be extended. You can add custom [info panes](../develop/how-to/add-info-panel.md) and [context actions](../develop/how-to/context-actions.md).

</details>

:::

![DB Hierarchy Browser](../uploads/features/db-hierarchy-browser.gif "DB Hierarchy Browser")

:::tip

You can reposition, resize, detach, or hide any panel. To learn more about customizing your workspace, see [Navigation](../datagrok/navigation.md).

:::

### Schema Browser

A **Schema Browser** provides a convenient way to explore all tables and columns belonging to a connection. Right-clicking a table or a column pulls up the object's context menu, clicking the object updates the **Context Panel**.

To open a **Schema Browser**, in the **Database Manager** right-click on a connection and select **Browse schema**.

![Schema browser](browse-schema.gif)

### Viewing schema as dataframe

Another way to explore a database schema is through an [interactive spreadsheet](../visualize/viewers/grid.md). To access this feature, in the **Database Manager**, right-click a connection and select **Open schema as table**.

## Querying and filtering data

Datagrok provides several tools for querying and filtering data:

<Tabs queryString="query-tools">
<TabItem value="editors" label="Editors" default>

Editors let you manually create or edit a query:

* [**Query Editor**](#query-view-editor) is the main interface for executing database queries. Using this tool, you can write and edit statements in SQL or other languages and use functions to add post-processing steps.
* [**Aggregation Editor**](#aggregation-editor) works with tables. Use it to manipulate, summarize, filter, and pivot table data.

</TabItem>

<TabItem value="built-in queries" label="Built-in queries">

Datagork offers two built-in queries:

* **Get All**: retrieves all table data. Use it with caution.
* **Get TOP 100**: retrieves the first 100 rows.

:::tip
  
To retrieve specific columns, hold down the Shift key on your keyboard while clicking on the desired columns in the schema. Once selected, use the context actions to run the query just for these columns. ![Get Columns](db-exploration-get-columns.png)

:::

</TabItem>

<TabItem value="Join tables" label="Join tables">

You can merge multiple tables within a database connection:

1. Right-click a table and select **Join tables** to open a dialog with a list of columns connected by keys.
1. Choose columns by selecting the corresponding checkboxes. Datagrok automatically generates the SQL statement and output preview based on your choices.
1. (Optional) Manually edit the SQL statement. As you make changes, the preview updates to reflect your edits.
1. When you're done, use the **Dropdown Arrow** control in the bottom left corner of the dialog to choose between the two options:
   * **Save as query**: This opens the **Query Editor** where you can edit the query further.
   * **Add result to workspace**: This opens the query output as a dataframe, which you can inspect and edit further.

![Create a join query](query-builder.gif)

</TabItem>
</Tabs>

:::note

The ability to query tables is connector-specific and may not be available for all.

:::

When you successfully create a query, it appears in the **Database Manager** under the corresponding connection or table. From there, you can run it, share it with others, manage access, and perform other tasks.

:::tip

To quickly find the exact query you're looking for, go to **Data** > **Queries** and type in the query's name or tag. To see all queries belonging to a connection, right-click it and select **Browse queries** from the list of options.

:::

### Query Editor

To open, right-click a database connection or a table and select **Add query** (for connections) or **New SQL Query...** (for tables) from the list of options.

The **Query Editor** editor has two tabs:

<Tabs queryString="query-editor">
<TabItem value="query" label="Query Tab" default>

This tab is where you write and edit queries. As you work on your query, you can preview the query output at any time by clicking the **Run/Refresh** button on the menu ribbon or by pressing the F5 key on the keyboard. The the query result is displayed as an interactive dataframe. In the dataframe, you can view object details, perform actions on columns, and more.

When your query is complete, give it a name and click the **Save** button. If you don't want to save the query, close the editor without saving.

![Create a database query](query-add.gif)

:::note

You can also add the query results to the workspace for further analysis. To do so, click the **Dropdown Arrow** control in the bottom left corner of the **Query Tab** and select **Add results to workspace**.

:::

</TabItem>

<TabItem value="transformations" label="Transformations Tab">

Use this tab to apply post-processing operations to your query:

![Query transformations](transformations.gif)

1. Select a function to apply to a query. Use checkboxes next to the operation categories to filter.
1. Set the function parameters in the dialog that opens.
1. Click **OK**.
1. Repeat the steps as needed.

Each transformation is recorded. To see the results of each record, locate the record on the left and click it. The data output for this step is displayed in the **Preview** section. You can also edit or remove transformation records as needed. To edit a transformation step, click the **Dropdown Arrow** control next to the transformation record and select from the list of options. To remove the transformation record entirely, click the **Delete** icon.

:::tip

Use the menu ribbon on top to add or remove columns quickly.

:::

:::note Developers

You can create custom transformation functions in R, Python, or any other language. See [Scripting](../compute/scripting.md).

:::

</TabItem>
</Tabs>

### Aggregation Editor

To open, in the **Database Manager** right-click a table and select **Aggregate data**.

Similar to the **Query Editor**, the **Aggregation Editor** has both the **Queries** and **Transformations** tabs. The **Query Tab** in both editors [looks and works the same](#query-editor).

To aggregate data, use the **Transformation Tab**. In this tab, you can choose which columns to include in your report and decide how to pivot and group them:

![Visual query](databases-visual-query.gif)

* To pivot unique values from one column into multiple columns in the query output, add the corresponding columns to the **Columns** field. To add the columns that stay the same, add them to the **Rows** field. These columns serve as keys, and their unique values become row identifiers.
* To group data:
   1. First, choose the aggregating function: next to **Measures**, click the **Add...** icon, then click **Aggregation...** and select from the list.
   1. Then, select an aggregation column: next to **Measures**, click the **Add...** icon again, and select a column from the list.

     :::note

     The list of functions may vary based on the connector and any packages that have been installed. If connectors expose custom functions, such as GIS functions in PostreSQL, Datagrok makes them readily available in the interface.

     :::

     :::note developers

     You can [create custom aggregation functions](../compute/scripting.md).

     :::

   1. (Optional) Repeat these steps to add additional aggregation parameters. To modify or remove a parameter, right-click and select the action from the list.
* Use **Filters** to specify which items are returned when you run a query.

## Parameterized queries

A query can be parameterized, allowing you to run the same query with different values instead of creating distinct queries for each scenario. 

The syntax for defining query parameters is based on [scripting](../compute/scripting.md) with additions specific to queries. Here's an example of a simple query that accepts one string parameter "productName":


>```sql
>--input: string productName                                               
>select * from products where productname = @productName
>```


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

    :::note
    
    A _parameter pattern_ is an expression to create a filter condition within the query. Much like a formula, it consists of field references, operators, and constants. To see a list of available patterns, see [Search patterns](../explore/data-search-patterns.md).

    :::

In most cases, Datagrok returns the query results as a dataframe. To return a value of a different data type, specify the type of the output parameter.

<details>

<summary> Example </summary>

To return a string with the semantic type `Molecule`:

 ```sql
 --output: string smiles {semType: Molecule}
 ```

</details>

## Running queries

To run a query, you have two options:

* Run a query from within a [query editor](#query-editors). Use this option to preview the query results when writing or modifying a query.
* Use the **Run** context action in all other cases. This opens the query output in a dataframe that you can [modify, save, and share with others](#query-results).

When you run a parameterized query, the system automatically generates a dialog for entering parameters. For each input parameter, there is an input control that prompts you to enter values. Datagrok automatically parses the entry and executes a parameterized, secure, connector-specific SQL query on the backend.

:::note developers

You can also expose the parameter dialog to end-users as an [application](../develop/how-to/create-package.md).<!--Mention?: when the cartdridge is not deployed on that particular database, the query returns an error-->

Similar to functions in JavaScript, queries are [functions](../datagrok/functions/functions.md) in Datagrok. This means you can [_call_ a query](../datagrok/functions/function-call.md) across multiple Datagrok plugins instead of copy-pasting the same query every time.

To run a query programmatically, see this [code snippet](https://public.datagrok.ai/js/samples/data-access/parameterized-query).

:::

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

## Sharing query results

You have two options to share query results:

* _Share the URL_: Each query output has a unique URL that includes query parameters. After you have executed a query, copy the URL from the address bar and share it with others. When this URL is accessed, a query gets re-executed, so you will always be looking at the latest data. To refine the query, change its parameters in the toolbox on the left. To access the query results from the link provided, users must have the necessary permissions to execute this query.
* _Share the project_: To share query results as a project, you need to [upload it](../datagrok/project.md#uploading-a-project). Use this option when you want to persist a specific layout or create a [dynamic dashboard](#sharing-query-output-as-dynamic-dashboard).

### Sharing query output as dynamic dashboard

When you run a query, you typically get back the dataframe. After that, you can perform various operations on that dataframe, such as data cleansing, transformation, [and more](../visualize/viewers/grid.md). All of these actions get recorded as a macros, and will be replayed when the query is re-executed. Once you have the data set up to your satisfaction, you can add viewers and create a visualization. You can save this visualization as a dashboard and share it with others. For parameterized queries, these dashboards can be dynamic. By changing query parameters right inside the dashboard, users can view different datasets and interact with data in real-time.

To save the query output as a dynamic dashboard, do the following:

1. On the **Sidebar**, click **Projects** > **Upload**.
1. In the **Upload project** dialog, enter a name and description (optional) in the fields provided.
1. Select how to store data:
    * Save the data as a static snapshot.
    * Store the data as a generation script by toggling the **Data sync** control. The query re-executes each time the project is opened. To learn more about dynamic data updates in projects, see [Dynamic data](../datagrok/project.md/#dynamic-data).
1. Click **OK** to upload the project.

![Dynamic dashboards](dynamic-dashboards.gif)

## Access control

In Datagrok, certain classes of objects we call _entities_ have a [standard set of operations](../datagrok/objects.md) that can be applied to them. These entities are connections, queries, tables, and table columns, all of which can be shared, assigned permissions, annotated, and more.

When you create an entity such as a connection or a query, it's private and visible to you only. You can share it to make it accessible to others. If you can't share a connection or a query, you may have insufficient permissions. Contact your Datagrok administrator for assistance.

To share:

1. Right-click the connection or a query and select **Share...**.
2. In the **Share** dialog, enter a user or a group that needs access and set corresponding permissions. These privileges can be given to individuals or to [groups](../govern/group.md) (which can be defined via dynamic filters). For more information on the access privilege model, see [Privileges](../govern/security.md#privileges).
3. (Optional) Add a description in the provided text field. If you don't want to notify the recipients, clear the **Send notification** checkbox.
   :::note

   To notify via email, enter the user's email in the identity/email field. The email contains a link to the shared item and entered description. When you enter a user or a group name, they are  notified via the Datagrok interface.

   :::

4. Click **OK** to share.

![Share a database connection](sharing-database-connections.gif)

:::note

Datagrok query belongs to the database connection for which it's created. It means you can’t share a query without sharing a connection. Deleting a connection also deletes a query.

:::

:::tip

Use the **Context Pane** on the right to inspect and quickly adjust access permissions to your databases, manage queries, view history and activity details, send comments to those you're sharing with, and more.

:::

## Data reproducibility

Any action performed on datagrok entities is reproducible and can be used in
automation workflows. For example, you can use data preparation pipeline to define jobs for data ingestion, postprocessing, and transformations.

To learn more about automating workflows using data preparation pipelines, see [Data preparation pipeline](data-pipeline.md).

## Resources

* Videos
  * [Database exploration](https://www.youtube.com/watch?v=YJmSvh3_uCM)
  * [Parameterized queries - Overview](https://www.youtube.com/watch?v=dKrCk38A1m8&t=1980s)
  * [Parameterized queries - Example](https://www.youtube.com/watch?v=sSJp5CXcYKQ)
  * [Using lists in parameterized queries](https://www.youtube.com/watch?v=meRAEF7ogtw)
* Tutorials
  * [Adding parameters to functions](../datagrok/functions/func-params-enhancement.md)
