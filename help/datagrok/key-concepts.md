---
title: "Key concepts"
tags:
  - key concepts
  - dataframe
  - table
  - entities
  - functions
keywords:
  - dataframe
  - table
  - function
  - entities
description: Key concepts
---

```mdx-code-block
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';
```

## Dataframe

A _dataframe_ (also a _table_ or a _data table_) is a table-like structure, composed of rows and named columns. These columns support various data types including `string`, `bool`, `int`, `bigint`, `double`, `qnum`, `float`, `datetime`. Unlike Excel, columns in Datagrok are strongly-typed, which ensures that every cell in a column adheres to the specified data type.

Dataframes are a key component of Datagrok's in-memory data engine. Serving as an in-memory database, the dataframe loads data from the associated [data source](../access/access.md#data-sources) as soon as you open or generate it in Datagrok. Subsequent data operations like [data transformation](../transform/data-wrangling.md) and [statistical computations](../explore) are performed directly on the dataframe. [Visualizations](../visualize), built on top of the data engine, use the same data and cashed statistics and calculations. This direct interaction allows [interactive exploration](../explore/exploratory-data-analysis.md) of big datasets - like the entire 2.7 million molecules from the ChEMBL database - right in the browser.

:::note developers

Datagrok's data engine runs entirely in the browser, but the same code can run on Datagrok servers for tasks like serialization. Learn how Datagrok [stores and processes data](../develop/admin/infrastructure.md#in-memory-database) and how to [optimize performance for large datasets](../develop/advanced/performance.md). For additional technical information, see [JS API: Dataframe](https://datagrok.ai/js-api/classes/dg.DataFrame).

:::

When parsing data, Datagrok automatically annotates dataframes and dataframe columns with basic metadata like its source or import time. It also automatically detects [certain semantic (logical) types](../discover/semantic-types.md#automatic-semantic-type-detection) like a zip code or a stock ticker. In addition, you can add [custom metadata and logical types](../discover/metadata.md).

Dataframes are visualized in a [Table view](#table-view) using a [grid](../visualize/viewers/grid.md) viewer.

## Entities

_Entities_ are classes of objects (such as users or queries) that share common features and operations applicable to them. For any entity, you can:

* Retrieve its URLs.
* Reference it in a [chat](../collaborate/chat.md), [markup](markup.md), or a [dashboard](../visualize/dashboard.md).
* Assign [privileges](../govern/authorization.md) for viewing, editing, sharing, or deleting particular instances.
* Use it as a parameter in [audit](../govern/audit.md) records.

Examples of entities include _queries_, _data connections_, _users_, _functions_, and more. Subject to your privileges, you can explore and interact with specific sets of entities through dedicated _views_ called _entity browsers_. Each entity browser is accessible via a unique URL. 

<details>
<summary>Entity list</summary>

* Entities that have designated entity browsers
  * Applications ([wiki page](../develop/develop.md) | [entity browser](https://public.datagrok.ai/apps))
  * Data connections ([wiki page](../access/access.md#data-connection) | [entity browser](https://public.datagrok.ai/connections))
  * Data jobs ([entity browser](https://public.datagrok.ai/jobs))
  * Data queries ([wiki page](../access/access.md#data-query) | [entity browser](https://public.datagrok.ai/queries))
  * Dockers ([entity browser](https://public.datagrok.ai/dockers))
  * File shares ([wiki page](../access/files/files.mdx#connecting-to-file-storage) | [open entity browser](https://public.datagrok.ai/files))
  * Functions ([wiki page](#functions) | [entity browser](https://public.datagrok.ai/functions))
  * Groups ([wiki page](../govern/group.md) | [entity browser](https://public.datagrok.ai/groups))
  * Models ([wiki page](../learn/predictive-modeling.md) | [entity browser](https://public.datagrok.ai/models))
  * Notebooks ([wiki page](../compute/jupyter-notebook.md) | [entity browser](https://public.datagrok.ai/notebooks))
  * Packages ([wiki page](../develop/develop.md) | [entity browser](https://public.datagrok.ai/packages))
  * Projects ([wiki page](project.md) | [entity browser](https://public.datagrok.ai/projects))
  * Repositories ([wiki page](../develop/develop.md) | [entity browser](https://public.datagrok.ai/repositories))
  * Scripts ([wiki page](../compute/scripting.md) | [entity browser](https://public.datagrok.ai/scripts))
  * Users ([wiki page](../govern/user.md) | [entity browser](https://public.datagrok.ai/users))
  * View Layouts ([wiki page](../visualize/view-layout.md) | [entity browser](https://public.datagrok.ai/view_layouts))
* Other entities
  * [Dataframes](#dataframe)
  * Database column info
  * Database table info
  * Data pipelines
  * [Function calls](#function-call).

</details>

When an entity is created, used, or starred by a user, Datagrok automatically annotates it with the following metadata:

* ID
* name
* author
* createdOn
* updatedOn
* commentedBy
* usedBy
* starredBy.

Datagrok automatically records and stores user interactions with entities, like executed queries or packages used. In addition, you can provide custom [metadata](../discover/metadata.md) with [packages][LINK] using [function annotations][LINK PLACEHOLDER].
 
<!-- not working: indexing. New feature: add metadata via UI-->

Metadata associated with entities is stored in a [Postgres database](../develop/admin/infrastructure.md#database) and can be used for [searching][LINK] within the platform. Metadata also enables auto-recommendations for layouts, models, scripts, and other data-specific functionalities. This functionality is essential for [data augmentation](../discover/data-augmentation.md) and extends the platform's capabilities.

<details>
<summary>Example: Data augmentation using metadata</summary>

In this example a Python-based script from the [Chem package](https://github.com/datagrok-ai/public/blob/master/packages/Chem/README.md) is annotated with [tags](../discover/tags.md) such as `demo`, `chem`, `rdkit`, `panel`. These annotations enable the script to be automatically applied to chemical datasets, diplay the results in the [info pane](../discover/info-panels.md), and show up in searches for entities applicable to cheminformatics (`#chem`).

```mdx-code-block
<Tabs>
<TabItem value="script" label="Script" default>
```

![Gasteiger partial charges script](img/concepts-entities-metadata-0.png)

```mdx-code-block
</TabItem>
<TabItem value="script-output" label="Script output">
```

![Gasteiger partial charges script output](../domains/chem/img/script-output-gasteiger-part-charges-0.png)

```mdx-code-block
</TabItem>
<TabItem value="script-output-in-info-pane" label="Script output in info pane">
```

![Script-based info pane](../domains/chem/img/script-output-info-pane-0.png)

```mdx-code-block
</TabItem>
</Tabs>
```

</details>

Example 2: Rendering

In this example, ... column cell renderer


<!--Becauses most data operations happen within the platform, this allows Datagrok to learn from user behavior. When activated, Datagrok's AI-driven self-learning component identifies usage patterns and provides users with recommendations. These recommendations could be data visualization ideas, predictive model suggestions, and more. Importantly, these recommendations are informed not just by your activity, but also by the actions of other Datagrok users, promoting knowledge sharing across diverse teams and time zones-->

:::note developers

You can define entities manually or programmatically using `grok.data.getEntities`, which fetches entities linked with a specific dataframe.

:::

## Functions

<!--Functions are a key feature of the Datagrok platform. They are used to perform a wide range of operations, such as querying databases, performing scientific computations, displaying user interface elements, or sending emails. Almost any action that can be executed in Datagrok is a function. 

Similar to programming languages, Datagrok functions accept inputs and produce outputs, and can be annotated. They can be written in any language, such as R, Python, Matlab, and others. 

<details>
<summary>Example: Calculating molecule properties</summary>

Consider the following scenario. You have a dataset with molecules and want to calculate the molecular weight for each compound. With just a few lines of code, you can create a function that uses the built-in cheminformatics library for the calculation:

```csharp

function calculateMolecularWeight(compoundTable) {
  var molWeightColumn = compoundTable.columns.addNew('Molecular Weight', DG.TYPE.FLOAT);
  for (var i = 0; i < compoundTable.rowCount; i++) {
    var compound = compoundTable.rows.get(i).compound;
    var molWeight = chem.getDescriptor(compound, 'Molecular Weight');
    molWeightColumn.set(i, molWeight);
  }
}

```

In this example, the script is a function written in C# that adds a new column to the table, iterates through the rows, calculates the molecular weight for each compound, and populates the corresponding 'Molecular Weight' column with the calculated values.

</details>

Each function type is associated with a specific tag that indicates its purpose. 

<details>
<summary>Examples</summary>

* `#app` for applications
* `#dashboard` for dashboards
* `#panel` for info panels
* `#init` for package initialization
* `#autostart` for automatic execution at platform startup
* `#semTypeDetector` for semantic types detectors
* `#cellRenderer` for custom cell renderers
* `#fileViewer` and `#fileExporter` for file viewers and exporters
* `#packageSettingsEditor` for package settings editors.

</details>

You can use these tags to search for particular functions either from the [platform's user interface](https://public.datagrok.ai/functions?q) or programmatically.

:::note developers

<details>
<summary>To search in code, use this code snippet</summary>

```js
const applications = DG.Func.find({tags: [DG.FUNC_TYPES.APP]});
```

</details>

<details>
<summary>To disable all package functions</summary>

To disable all package functions (for debug purposes), use the
`initPackageFunctions=false` flag in the start URL, such as
`https://public.datagrok.ai?initPackageFunctions=false`.

</details>

:::

Functions can be executed either on the client-side (in the browser), server-side (using a compute engine), or both.-->

Datagrok supports functions, which is an incredibly powerful concept. Virtually any action that can be executed within the platform is a function. Functions vary in terms of their specific functionality and execution methods. Some functions execute on the server, while others run in the browser. Some functions perform scientific computations, while others display UI elements.

<details>
<summary>Example: Different function types</summary>

* [Querying](../access/access.md#data-query) external Postgres database.
* Executing a [JavaScript function](../develop/develop.md) in the browser that uses Grok API for integration purposes.
* Calculating Sin(PI).
* Deleting a column from the table.
* Sending an email.
* Applying a [predictive model](../learn/predictive-modeling.md) to a dataset.
* Calculating molecule properties using a Python [script](../compute/scripting.md).
* Displaying a dialog.

</details>

### Function properties

Despite these differences, all functions
share the same mechanism, which means they have the same set of features:

* _Scriptable_: Each function call is represented as a string, allowing execution from the [Console](#console) and integration into larger scripts.

  <details>
  <summary>Example: Using a function's output in UI elements</summary>
  
  In the example below, a [Python script based on RDKit](https://public.datagrok.ai/script/276a5929-6f21-5105-8eec-576845aabae0) is used to calculate Gasteiger partial charges and generate a visual representation of the results. The script is implemented within a UI component known as an [info pane](../discover/info-panels.md), which dynamically updates as you browse the dataset.

  <Tabs>
  <TabItem value="script" label="Script" default>

  ![Gasteiger partial charges script](../domains/chem/img/script-gasteiger-part-charges-0.png)

  </TabItem>
  <TabItem value="script-output" label="Script output">

  ![Gasteiger partial charges script output](../domains/chem/img/script-output-gasteiger-part-charges-0.png)

  </TabItem>
  <TabItem value="script-output-in-info-pane" label="Script output in info pane">

  ![Script-based info pane](../domains/chem/img/script-output-info-pane-0.png)
  </TabItem>
  </Tabs>
  </details>

* _Reproducible_: Every function call initiated from the user interface can be automatically recorded, enabling effortless replication and reuse of individual actions or sequences of actions.

  <details>
  <summary>Example: Data transformation with macros</summary>
  
  In Datagrok, you can manipulate data using visual tools. Each transformation you perform through the user interface, such as clicking the **Delete selected rows** icon, corresponds to a function. When activated, Datagrok automatically records a history of all changes made to a table in the [Console](#console). This feature is especially useful for data wrangling, as it allows you to understand how a specific table was created or replicate the process on a new table with a similar structure. 
  
  To use:
  
  1. Activate auto-recording by doing [**THIS**].
  1. Perform necessary operations on the table from the user interface.
  1. Open the **Console** by pressing the **Tilde** (`~`) key. The Console displays a log of the actions performed. 
  1. You can save the recorded script and reuse it later to clean similar datasets .[NEED MORE DETAILS?]
  
  :::note

  To access the history of changes at a later date, you need to upload the table to the server. To view the history, open the table and expand the **History** pane in the [Context Panel](#context-panel).

  :::

  </details>

* _Findable_: Datagrok provides a dedicated function browser for easy searching based on function names or tags. To access it, click the **Functions** icon on the **Sidebar**, then select **Functions**.
* _Introspectable_: You can access metadata for function parameters programmatically, as well as produce a table with all parameters used for a particular function. 
* _Secure_: Functions are [entities](#entities) and can be assigned [privileges](../govern/authorization.md) for viewing, editing, sharing, and deleting.
* _Auditable_: Information about function execution, including the user, timestamp, and used parameters, is recorded in the audit log.
* _Runnable_: Datagrok dynamically generates UI for editing parameter values.
* _Linkable_: You can create links to functions by dragging and dropping them into [conversations](../collaborate/chat.md) or dashboards. You can then run a function by right-clicking the link and choosing **Run...**
* _Applicable_: Functions can be used in workflow designers and [query result transformations](../transform/recipe-editor.md).

### Data flow

When dealing with sensitive data, it is important to understand its flow. It is usually not a concern if you run an
enterprise version of the Datagrok platform in a virtual private cloud, but it might be a bigger issue when using the
public environment.

While we try to do as many computations on the client-side (in the browser) as possible, certain operations, such as
[training a predictive model](../learn/predictive-modeling.md),
[running an R script](../compute/scripting.md), or
[computing chemical descriptors](../domains/chem/chem.md#descriptors-and-fingerprints)
run on a [server](../develop/admin/architecture.md#compute-engine). When such operation executes, the
relevant part of the input gets sent to the server, where it gets processed, and the result is sent back to the client.
Results of the computations reside on the server until either the client retrieves it, or an automatic cleanup happens.
Neither inputs nor outputs can be accessed by functions executed by other users. Only administrators can access compute
server, where the data is stored or processed. All traffic between client and server
is [TLS-encrypted](https://en.wikipedia.org/wiki/Transport_Layer_Security).

Most of the actions performed by users are logged for the usage analysis and [audit](../govern/audit.md) purposes.
For logging, only the metadata associated with the parameters (such as table and column names) will be additionally
stored in the database. Additionally, if a table passed as a parameter is already residing on a server, the
corresponding audit record will contain a reference to it.

### Filtering

You can use these fields to filter functions with [smart search]:

| Field       | Description                                 |
|-------------|---------------------------------------------|
| ID          |                                             |
| name        |                                             |
| runs        | list of [FuncCall](functions/function-call.md) object |
| createdOn   |                                             |
| updatedOn   |                                             |
| author      | [User](../govern/user.md) object         |
| starredBy   | [User](../govern/user.md) object         |
| commentedBy | [User](../govern/user.md) object         |
| usedBy      | [User](../govern/user.md) object         |

### Videos

[![Functions](../uploads/youtube/functions.png "Open on Youtube")](https://www.youtube.com/watch?v=p7_qOU_IzLM&t=724s)

See also:

* [Console](../datagrok/navigation.md#console)
* [Scripting](../compute/scripting.md)



[SCRIPTING FOR NON-DEVELOPERS]

Datagrok makes it easy for non-developers to script functions using a simple syntax that anyone can learn. You can write and execute functions directly in the Datagrok Console, or create more complex scripts that can be saved and shared with others.

[FUNCTION FEATURES]

Funactions are entities, which means...

You can determine who can execute them and under what conditions. For example, you might want to restrict access to a function that deletes data or performs other sensitive operations. Function roles can be managed through the Datagrok Security Manager. See ... for details

In addition, ....

Despite the differences in functionality and execution, functions share the same features and mechanism. They are:

<details>
<summary>Findable</summary>

You can search for functions (on the **Sidebar**, click **Functions** and search by name or tag ).

</details>

<details>
<summary>Auditable</summary>

You can see who executed a function, when it was executed, and what parameters were used. In addition, you can get a table with all parameters used for a specific function.

</details>

<details>
<summary>Introspectable</summary>

You can programmatically find out the function parameters' metadata.

</details>
<details>
<summary>Runnable</summary>

You can run a function from the UI or programmatically. For parameterized functions, Datagrok automatically generates a dialog for entering parameters, which can be exposed to users as a [standalone application](../develop/how-to/create-package.md).

![Gasteiger partial charges script](../domains/chem/img/script-gasteiger-part-charges-0.png)

</details>

<details>
<summary>Scriptable</summary>

Each function call can be represented as a string and executed from the [Console](#grok-script) or as part of a bigger script.

In the example below, a [Python script based on RDKit](https://public.datagrok.ai/script/276a5929-6f21-5105-8eec-576845aabae0) is used to calculate Gasteiger partial charges and generate a visual representation of the results. The script is implemented within a UI component known as an [info pane](../discover/info-panels.md), which dynamically updates as you browse the dataset.

<Tabs>
<TabItem value="script" label="Script" default>

![Gasteiger partial charges script](../domains/chem/img/script-gasteiger-part-charges-0.png)

</TabItem>
<TabItem value="script-output" label="Script output">

![Gasteiger partial charges script output](../domains/chem/img/script-output-gasteiger-part-charges-0.png)

</TabItem>
<TabItem value="script-output-in-info-pane" label="Script output in info pane">

![Script-based info pane](../domains/chem/img/script-output-info-pane-0.png)
</TabItem>
</Tabs>
</details> 

<details>
<summary>Secure</summary> 

Functions are _entities_ and as such, have access privileges associated with them. For example, a function can be made available only to users in the "Managers" group.

Some of the most popular privileges
are: `view`, `edit`, `delete`, and `share`. Those privileges can be given to individual users, or
to [groups](../govern/group.md). For more information on the access privilege model, check
out [privileges](../govern/security.md#privileges).

</details>

<details>
<summary>Linkable</summary>

You can drag-and-drop a function to a conversation or dashboard and run by right-clicking on the link and choosing "Run..."

</details>

<details>
<summary>Usable</summary>

Functions can be recorded as a macros and used in workflow designers and query result transformations.

</details>


### Function call

Function Call is a result of executing a data job, [data query](../access/access.md#data-query),
[script](../compute/scripting.md), or any other _function_.

Each function call contains the following data:

* Function
* User that triggered job execution
* Started on
* Completed on
* Status
* Runs produced as a result of executing child actions.

:::tip

<details>
<summary>You can use these fields to filter action runs with [smart search](smart-search.md)</summary>

* ID
* name
* action (a `Func` object)
* childRuns (a list of `FuncCall` object)
* parentRun (a `FuncCall` object)
* status
* started
* finished
* createdOn
* updatedOn

</details>

:::

[FUNCTION PARAMETERS]

Functions in Datagrok can accept parameters, which are values that are passed to the function when it is called. Parameter values can be provided in a variety of formats, including numbers, strings, and tables. When creating a function, you can annotate the parameters to provide additional information about their expected data type and default values.


When you call a function in Datagrok, you pass in the required parameters as arguments. The function then executes and returns a result, which can be stored in a variable or displayed to the user.


[PARAMETER ENHANCEMENTS]

Datagrok provides a number of enhancements to function parameters, including the ability to specify default values and metadata about the parameter. This information can be used to provide additional context for users and improve the readability of your scripts.






