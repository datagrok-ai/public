<!-- TITLE: Datagrok -->
<!-- SUBTITLE: -->

# Datagrok: Swiss Army Knife for Data

Data is one of the most critical elements in today's world. How fast you can convert data to insights, and insights to actions, is paramount. We created Datagrok to help people understand and act on their data.

## Why Datagrok?

An end-to-end solution. An intuitive UI makes it simple for non-technical users to perform all data-related tasks from the main screen: [access](#access), [transform](#transform), [visualize and explore](#visualize-and-explore), [model and test hypotheses](#machine-learning), [compute](#compute), [discover](#discover), and [share](#share) anything with anyone.

Superior web-based performance and scalability. Thanks to our [proprietary in-memory database](develop/advanced/performance.md##in-memory-database), you can interactively explore data sets with _tens of millions of rows_ in the browser. The platform can scale up a billion rows and a million columns and is limited only by the amount of physical memory available on the client computer.

Data-agnostic and domain-intelligent. You can work with any data, including [complex multidimensional and heterogeneous data](#apply), such as life sciences data.

Flexible. The platform has an open architecture and component-based design. You can rapidly develop and deploy new functionality as plugins or cherry-pick and combine data transformations, statistical computations, predictive models, [packages](https://github.com/datagrok-ai/public/tree/master/packages), and custom applications into your fit-for-purpose solution. <!--todo create a wiki page describing fit-for-purpose solution capability-->Of course, you can also use it out-of-the-box as a self-service analytical tool or a BI application.

Secure and meets the needs of today's digital enterprise. Datagrok is easy to deploy and integrate. Seamlessly bring in data from multiple silos, manage workflows, users, data, and more.

Learns from observed user actions. Datagrok uses numerous AI techniques to identify usage patterns and pushes this information as insights via built-in [data augmentation, collaboration, and knowledge dissemination](#discover) capabilities.

Rich and fast growing ecosystem. We partner with industry and academia<!--link to the partnerships section on our website, when ready--> to bring you a [rich ecosystem of open source plugins](https://github.com/datagrok-ai/public/tree/master/packages). By allowing anyone customize cutting-edge solutions or access the latest features, we lower individual ownership costs, help organizations adapt to the ever-changing environment, drive progress through precompetitive collaboration, and democratize data science.

[Launch Datagrok right now](https://public.datagrok.ai/) and see for yourself!

## Access

* [30+ connectors](access/data-connection.md) to all major databases
* 1,000+ services exposed via [OpenAPI](access/open-api.md)
* Drag-and-drop files to open ([10+ formats](access/importing-data.md)), or
  browse [file shares](https://public.datagrok.ai/files)
* [Visually explore](access/db-exploration.md) and manage relational databases
  using [schema browser](access/db-exploration.md#schema-browser)
  and [visual query](access/db-visual-query.md)
* Connect to [thousands of public datasets](access/public-datasets.md)
* Automate via [data preparation pipelines](access/data-pipeline.md)

## Transform

* [Aggregate, join, filter and edit data](transform/data-wrangling.md) right in the browser
* Record and apply [macros](overview/navigation.md#recording-macros)
* Use 500+ available [functions](overview/functions/function.md), or write your own in R, Python, or JavaScript
* Visually edit [pipelines](transform/job-editor.md)
  and [query transformations](transform/recipe-editor.md)

## Visualize and explore

* [Exploratory data analysis](explore/exploratory-data-analysis.md) in the browser
* [30+ high-performance interactive viewers](visualize/viewers.md)
* [Powerful integration with any visualizations available in R, Python, or Julia languages](visualize/viewers/scripting-viewer.md)
* Built-in viewer features: [regression lines](visualize/viewers/scatter-plot.md), confidence intervals,
  correlations, [statistical tests](learn/data-science.md)
* Automatic detection of outliers, [missing values](transform/missing-values-imputation.md), wrong data types
* [Dashboards](/visualize/dashboard.md)

## Machine learning

* [Train, assess, apply, share models](learn/predictive-modeling.md)
* Use in [pipelines](transform/job-editor.md)
* [Seamless integration with Python, R, or any other language](compute/scripting.md)
* [Statistical hypothesis testing](learn/data-science.md)

## Compute

* [A next-generation environment for scientific computing](compute/compute.md)
* [Scripting](compute/scripting.md)
* [Jupyter notebooks](compute/jupyter-notebook.md)

## Share

* Share anything with anyone
* [Rich open source ecosystem](https://github.com/datagrok-ai/public/tree/master/packages)
* Cross-pollinate knowledge via the knowledge base, [discussions](collaborate/chat.md)
  and [forums](collaborate/forum.md)

## Discover

* Push ideas to users via [data augmentation](discover/data-augmentation.md)
* [Self-learning platform](learn/self-learning-platform.md): the more you use it, the better it gets

## Develop

* [Custom applications](develop/how-to/build-an-app.md)
* [User interface](develop/ui.md)

## Deploy and integrate

* Multiple [hosting options](develop/admin/hosting-options.md)

## Manage

* [Roles, groups and privileges](govern/security.md)
* Flexible [authentication](govern/authentication.md)
* Create [pipelines](transform/job-editor.md), schedule [jobs](access/data-job.md), and set up alerts
* [Customizable by IT](develop/admin/it-customizations.md)
* Learning resources: [interactive help](overview/navigation.md#help), documentation<!--todo: add a link-->, [community forum](https://community.datagrok.ai/), and more

## Govern

* Central [metadata](discover/metadata.md)-annotated [catalogue](https://public.datagrok.ai/)
  of [projects](https://public.datagrok.ai/projects), [queries](https://public.datagrok.ai/queries),
  and [connections](https://public.datagrok.ai/connect)
* [FAIR](discover/fair.md): findable, accessible, interoperable, reusable
* [Secure by design](govern/security.md)
* Built-in [data provenance](govern/data-provenance.md), data lineage, impact analysis
  , [usage analysis](govern/usage-analysis.md), and [audit](govern/audit.md) tools

## Apply

* [Cheminformatics](domains/chem/cheminformatics.md)
* [Bioinformatics](domains/bio/peptides.md)
* Text analytics, [natural language processing](https://github.com/datagrok-ai/public/tree/master/packages/NLP)
  with cloud-based machine translation
* [Location Analytics](https://github.com/datagrok-ai/public/tree/master/packages/Leaflet)
* [Digital signal processing](https://github.com/datagrok-ai/public/tree/master/packages/DSP)
* [Biosignal Processing](https://github.com/datagrok-ai/public/tree/master/packages/BioSignals)
