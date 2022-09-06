<!-- TITLE: Exploratory data analysis -->
<!-- SUBTITLE: -->

# Exploratory data analysis

Before we can learn from data, we need to understand
it. [Exploratory data analysis](https://en.wikipedia.org/wiki/Exploratory_data_analysis) (EDA) is a process of
performing initial investigation on data to discover patterns, spot anomalies, test hypothesis, and check assumptions.

By its nature, EDA is visually-driven. Most of today's datasets are too big, too complex, and diverse to be explored in
a tabular format or by statistical means alone. On the other hand, humans evolved to understand complex information
visually and are better than computers at detecting patterns and anomalies.

Interactivity is key. We may not know what we are looking for until we extract knowledge from data and update our
understanding as we go. To uncover insights that otherwise may go unnoticed, we need to be able to quickly change both
what we are viewing and how we are viewing it:

* Look at data from multiple perspectives at once
* Zoom in and filter
* Manipulate, edit, and add data
* Get details on demand
* Select rows of interest, and see how they compare to other row sets.

From the ground up, we designed Datagrok for visually-driven EDA of big, complex datasets. Unlike other tools that use
conventional client-server architecture,
Datagrok's [proprietary in-memory database](../develop/advanced/performance.md#in-memory-database) makes it possible to
analyze _millions of columns_ and
_billions of rows_ at the speed of thought right in your browser.

With Datagrok, you can:

* [Seamlessly load data from any data source](../access/file-browser-and-file-shares.md). Datagrok supports all popular databases,
  multiple [file formats](../access/supported-data-sources.md#supported-file-types) and is both data-agnostic and
  domain-intelligent. <!--TODO link to a section on domains once ready-->

* Visualize the data using domain-specific value renderers (such as molecules on scatter plot axes).

* Analyze big datasets that other tools struggle with (billions of rows, or millions of columns).

* Use multiple interactive tools to [wrangle data](../transform/data-wrangling.md) right from your visualization
  workspace. [Cluster data](cluster-data.md), [impute missing values](../transform/missing-values-imputation.md), find
  and treat duplicates and outliers.

* Visualize data at the click of a button using [30+ native viewers](../visualize/viewers.md). We support all popular
  visualizations (
  like [scatterplots with built-in regression lines](../visualize/viewers/scatter-plot.md#regression-line)
  or [box-plots with built-in statistical tests](../visualize/viewers/box-plot.md#t-test)) and certain domain-specific
  viewers, such as [chemically-aware   viewers](../domains/chem/chemically-aware-viewers.md). The tabular viewer, [_
  grid_](../visualize/viewers/grid.md), is extremely powerful. Some of its features include:

  * Viewing datasets with _millions of columns_ and _billions of rows_
  * Dataset overview, including summary statistics for numerical data columns and distribution for categorical data
    columns
  * Custom cell renderers for molecules, sequences, dose-response curves, and sparklines
  * Editing datasets (for example, adding new molecules using [sketchers](../domains/chem/sketcher.md))

* Filter, zoom in and out, aggregate, pivot, and cross-link data. All our viewers work in tandem and are customizable,
  [high-performant, and interactive](../develop/advanced/performance.md#viewers).

* Perform calculations on data using predefined statistical [functions](../datagrok/functions/function.md).
* Get details on demand using a variety of [widgets](../visualize/widgets.md), including customizable
  [tooltips](../explore/select-tooltip-columns.md#viewer-tooltips) and
  context-driven [info panels](../discover/info-panels.md).

* Build on collective knowledge of Datagrok users. Using
  built-in [data augmentation capabilities](../discover/data-augmentation.md), Datagrok understands the nature of your
  data, and offers actionable insights based on it. For example, the platform
  [automatically suggests visualizations](../visualize/view-layout.md#layout-suggestions) for datasets or predicts
  properties for chemical structures.

You can also leverage Datagrok's component-based architecture to extend or create any element you like. For example, you
can [add custom viewers](../develop/how-to/develop-custom-viewer.md) or develop new functions
in [R, Python, or Julia](../compute/scripting.md).

Each of these actions can be [automated](../datagrok/functions/function.md#macros) and used in
[pipelines](../access/data-pipeline.md). Sharing the results of your analysis is easy and
[secure](../govern/security.md).<!--TODO rewrite for clarity-->

With Datagrok, anyone can use their domain knowledge and perceptive abilities to explore data and uncover its meaning.

## Examples

<!-- markdownlint-capture -->
<!-- markdownlint-disable -->

<div class="card" style={{width:"512px",}}>
<iframe src="https://www.youtube.com/embed/67LzPsdNrEc?vq=hd1080&rel=0&color=white&autohide=0" width="512" height="288" frameborder="0"></iframe>
  <div class="card-body">
    <h3 class="card-title">Interactive Data Visualization</h3>
    <p class="card-text">An overview of some of the visualization capabilities of the Datagrok platform, including the concepts of views, viewers, selection, filter, and layouts.</p>
    </div>
</div>

<div class="card" style={{width:"512px",}}>
<iframe src="https://www.youtube.com/embed/tVwpRB8fikQ?vq=hd1080&rel=0&color=white&autohide=0" width="512" height="288" frameborder="0"></iframe>
  <div class="card-body">
    <h3 class="card-title">Coffee Company</h3>
    <p class="card-text">How do we choose the best location for a new coffee place, given the historical sales data? Datagrok to the rescue! In less than 20 minutes, we achieve the following:<br />
                         • Retrieve historical data from the Postgres database<br />
                         • Explore, visualize, and clean the dataset<br />
                         • Impute missing values<br />
                         • Extract census data from the long/lat coordinates<br />
                         • Perform multivariate analysis<br />
                         • Build multiple predictive models, and assess their performance<br />
                         • Build an interactive map for predicting sales<br />
                         • Deploy the results as an app to all users in our company<br />
    </p>
  </div>
</div>

<!-- markdownlint-restore -->

## See also

* [Viewers](../visualize/viewers.md)
* [Data science](../learn/data-science.md)
* [Predictive modeling](../learn/predictive-modeling.md)
