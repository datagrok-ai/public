---
title: "Datagrok: Swiss Army Knife for Data"

---

## Why Datagrok?

Datagrok helps you understand data and take action.

It's fast and powerful: you can load the entire ChEMBL database (2.7 million
molecules) in your browser, run substructure searches, apply filters, visualize,
and interactively explore the chemical space.

Datagrok goes beyond standard data analytics. You can access data from any
source, catalog it, analyze and visualize it, run scientific computations, train
and apply models, and do more. Need a specific tool or functionality? Easily
integrate or add your own code. Datagrok's plugin architecture makes it easy to
deliver cohesive, fit-for-purpose solutions. 

#### Access

Get your data from anywhere - databases, web services, file shares, pipelines.
If it's machine-readable, we can work with it!

* [40+ connectors](../access/databases/connectors/connectors.md) to all major
  databases and file shares (or [create your own](../access/databases/create-custom-connectors.md))
* Support for [OpenAPI](access/open-api.md) and access to [public datasets](../access/public-datasets.md)
* [20+ file formats, 30+ molecule structure formats](../access/files/supported-formats.md). Drag-and-drop files to open
* [Browse relational database schemas](../access/databases/databases.md#schema-browser)
* Create, edit, and debug queries with [visual tools](../access/databases/databases.md#working-with-queries)
* [Annotate queries](../access/databases/databases.md#parameterized-queries)
  and save query results as [dynamic dashboards](../access/databases/databases.md#creating-dynamic-dashboards-for-query-results).

[Learn more about data access](../access/access.md).

#### Govern

Use catalogs, data lineage tools, audit, and usage analysis to take control.
Your data is [FAIR](../govern/catalog/fair.md) and secure.
* Control who, what, where, and how: [roles, groups, and privileges](../govern/access-control/access-control.md#authorization), 
flexible [authentication](../govern/access-control/access-control.md#authentication), 
[secrets managers](../govern/access-control/data-connection-credentials.md)
* Centralized [metadata](concepts/objects.md#metadata)-annotated catalog of
[entities](concepts/objects.md). Powerful "everything" browser for
managing data, connections, users, and more
* Built-in [data provenance](../govern/audit/data-provenance.md), data lineage, impact analysis, [usage analysis](../govern/audit/usage-analysis.md), and [audit](../govern/audit/audit.md) tools
* Global search. 

#### Transform

Automatically generate macros from data transformations and use them on new datasets. 

* [Aggregate, join, filter, and edit data](../transform/transform.md), from the UI or programmatically
* Use [500+ available functions](https://public.datagrok.ai/functions?q), or
  write your own in JavaScript, Python, R (or any other language that compiles
  to WASM)
* Record and apply [macros](navigation/panels/panels.md#recording-macros), use in pipelines
* Visually edit [query transformations](../transform/query-transformations.md).

[Learn more about functions](concepts/functions/functions.md).

#### Explore

Slice, dice, and visualize your data. Render millions of data points
interactively and find patterns. Build dynamic dashboards in seconds. Leverage
[metadata](concepts/objects.md#metadata) for automated data enrichment and contextual
suggestions.

* [50+ interactive viewers](../visualize/viewers/viewers.md) for synchronized, dynamic dashboards
* [Integration with visualizations in R, Python, or Julia](../visualize/viewers/scripting-viewer.md)
* Built-in [regression and formula lines](../visualize/viewers/scatter-plot.md#calculations-and-trends),
  confidence intervals, correlations, and statistics
* Automatic detection of outliers, [missing values](../explore/missing-values-imputation.md),
 or incorrect data types
* Adaptive UI and data-specific suggestions.

#### Compute

Write in any language, annotate, publish, and apply scientific models, methods,
and apps. Solve differential equations and run simulations for complex processes<!--(e.g., for [bioreactors](like) and [PKPD])-->.

* [500+ available functions](https://public.datagrok.ai/functions), or write your own in R, Python, or JavaScript
* [Metadata-annotated](../compute/compute.md#metadata) scripts with [cross-language support](../compute/compute.md#functions-and-cross-language-support)
* [Scalable](../compute/compute.md#scalable-computations) and 
[reproducible](../compute/compute.md#reproducible-computations) computations, [model lifecycle management](../compute/compute.md#models)
* [Auto-generated UI](../compute/compute.md#user-interface).

Learn more about [Compute](../compute/compute.md).

#### Learn

No-code modeling. State-of-the-art cheminformatics engines and ML toolkit included.

* [Train, assess, apply, and share models](../learn/learn.md) (or integrate your own)
* Native support for R, Python, Julia, Matlab, and Octave
* Open any dataset with a [Jupyter notebook](../compute/jupyter-notebook.md)
* ML toolkit: [statistical hypothesis testing](solutions/domains/data-science.md#statistical-hypothesis-testing), [multivariate analysis](../explore/multivariate-analysis/pls.md), [dimensionality reduction](../explore/dim-reduction.md), [data clustering](../explore/cluster-data.md), [variance analysis](../explore/anova.md).

#### Collaborate

Share anything with anyone. Collaborate on decision-making. Use an open source
ecosystem to save costs and innovate.

* Share within Datagrok, as a URL link, or integrate: REST API, JS API, 
or embed as an iframe
* [50+ open source plugins](https://github.com/datagrok-ai/public/tree/master/packages),
  including specialized ones for
  [cheminformatics](solutions/domains/chem/chem.md),
  [bioinformatics](solutions/domains/bio/bio.md),
  [NLP](solutions/domains/nlp/nlp.md), and others
* Data annotations, [team discussions](../collaborate/chat.md)
* [Community forum](../collaborate/forum.md) for ideas, support, and feedback.

#### Extend

Customize anything, from context actions to UI elements. Fast development and
deployment time with seamless integration.

* [JavaScript API](../develop/packages/js-api.md) for [extending Datagrok](../develop/packages/extensions.md#what-can-be-extended)
* App marketplace: use or customize [ours](https://public.datagrok.ai/packages), [build your own](../develop/how-to/build-an-app.md), or integrate with third party apps
* [Developer tools](../develop/dev-process/tools/inspector.md), [UI toolkit](../develop/advanced/ui.md)
* Comprehensive help: [wiki](../develop/develop.md), [exercises](../develop/onboarding/exercises.md), [community forum](https://community.datagrok.ai/).

## Who is it for?

**Data**: Datagrok is optimized for structured, tabular data. It automatically
detects the semantics, like zip codes or molecules, and has built-in support
for areas like [cheminformatics](solutions/domains/chem/chem.md),
[bioinformatics](solutions/domains/bio/bio.md), [data science](solutions/domains/data-science.md),
 and others. Need more? Create [your own plugin](../develop/how-to/create-package.md).

**Skillset**: Datagrok is for anyone who works with data: 

* Chemists analyzing SAR tables? [Perfect fit](solutions/domains/chem/chem.md#chemically-aware-spreadsheet).
* Data analysts? Drag and drop your local files to start analyzing. 
* Data scientists mapping new store locations? [Excellent for strategic planning](https://www.youtube.com/watch?v=tVwpRB8fikQ).
* Research scientists running complex simulations? [Absolutely](../compute/compute.md).
* Data engineers? Automatically convert queries to [dynamic dashboards](../access/databases/databases.md#creating-dynamic-dashboards-for-query-results), no coding needed. 
* Developers? Quickly develop and test data-driven applications.   

**Team size**: Datagrok is for individuals and teams of all sizes - from
startups<!--insert link to customer stories--> to large enterprises<!--insert link to customer stories-->. The
platform is [enterprise-ready](solutions/enterprise/enterprise.md),
[scalable](../develop/under-the-hood/scaling.md), and ideal for sharing and collaboration.

## What makes it so flexible?

Our mission is to help anyone understand their data, even in complex scenarios:

* Data that's scattered across various data sources
* Data that needs specialized, domain-specific tools
* Teams that have different data needs and expertise. 

Here's how we do it. 

**JS API**: With JS API, you aren't confined to pre-built features or interfaces. Add new
data formats, connectors, transformations, augmentations, dynamic calculations,
UI elements, full-scale applications, workflows, and
[more](../develop/function-roles.md). The API also provides seamless
integration with data sources and other tools, crucial for large enterprises
combatting data silos and complex data ecosystems.

**Functions**: In Datagrok, every task is a
[function](concepts/functions/functions.md) that can be
[annotated](concepts/functions/func-params-annotation.md). Annotations make
functions versatile, allowing them to work on their own or within larger
scripts, no matter the function's language or role. This means you can use
functions as blocks to build on your team's collective expertise while fully
leveraging Datagrok's capabilities. (See the cheminformatics example below).

<!---

<details>
<summary>Example</summary>

In this example, a [Python script based on RDKit](https://public.datagrok.ai/script/276a5929-6f21-5105-8eec-576845aabae0)
calculates and visualizes Gasteiger partial charges. When you run the script
explicitly, Datagrok shows a dialog for sketching a query molecule and
visualizes the results. In this case, however, the script is also tagged as a
`panel`. This instructs Datagrok to show the results as an interactive UI
element that updates dynamically for the current molecule.

|  | |
|--|--|
|Script with auto-generated UI based on annotation|![Gasteiger partial charges script](solutions/domains/chem/img/script-gasteiger-part-charges-0.png)|
|Script output|![Gasteiger partial charges script output](solutions/domains/chem/img/script-output-gasteiger-part-charges-0.png)|
|Script output in info pane|![Script-based info pane](solutions/domains/chem/img/script-output-info-pane-0.png)|

</details>

-->

**Semantic types**: [Semantic data types](../govern/catalog/semantic-types.md) provide domain-specific customization:

* Automatic detection of domain-specific data types
* Domain-specific menus and context actions
* Custom data rendering, including spreadsheets and visualizations
* Specialized data editing and filtering interfaces
* Domain-specific calculation and data processing functions
* Fit-for-purpose apps built on top of Datagrok.

[See this example for cheminformatics](solutions/domains/chem/chem.md#chemically-aware-spreadsheet).

<!--
<details>
<summary>Example</summary>

For example, when you open a CSV file containing molecules in the SMILES format, the following happens:

* Data is parsed, and the semantic type `molecule` is assigned to the corresponding column.
* Molecules are automatically rendered in 2D in the spreadsheet.
* Column tooltip now shows the most diverse molecules in your dataset.
* Default column filter is now a sketcher-driven substructure search.
* A top menu item labeled **Chem** appears.
* Molecule-specific info panes, such as **Toxicity** or **Drug Likeness**, appear on the right.

![img](solutions/domains/chem/img/chem-exploration.png)

</details>

-->

## What makes it so fast?

Our goal is to let you explore at the speed of thought. To achieve this, we
designed Datagrok from scratch:

* **Data engine**: In-memory columnar database that runs on both server and web
  browser. Fast random access, efficient data storage, aggregation, compression,
  filtering, transformation, and caching.

* **Native viewers**: Access the data engine directly for maximum performance.
  They share statistics, cached calculations, and cooperate on tasks like
  filtering or selection.

* **App server**: Uses the data engine to exchange binary-optimized datasets
  with the client. Custom ORM to efficiently work with metadata in Postgres.

* **Compute engine**: Supports multiple languages working with
  binary-optimized datasets. Scales well. GPU acceleration of ML routines.
  Supports custom Docker containers.

Learn more about [Datagrok's architecture](../develop/under-the-hood/architecture.md) 
and [performance optimization](../develop/under-the-hood/performance.md).

## Solutions

* [Self-service analytics](solutions/domains/use-cases/eda.md)
* [Data science](solutions/domains/data-science.md)
* Life sciences
  * [Chem](solutions/domains/chem/chem.md)
  * [Bio](solutions/domains/bio/bio.md)
* [NLP](solutions/domains/nlp/nlp.md)
* [Enterprise IT](solutions/enterprise/enterprise.md)
* [Plugins](plugins.md)

<!--

### Customer stories

* A big pharma [problem/solution] company by building a custom application for small molecules SAR
* A biotech startup ...

-->