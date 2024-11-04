---
title: "Exploratory data analysis"
sidebar position: 4
keywords:
 - EDA
 - exploratory data analysis
 - interactive exploration
 - interactive dashboards

---

>We help you explore at the speed of thought. [Learn what makes it possible](../../../datagrok.md#what-makes-it-so-fast).

The goal of data analysis is to derive knowledge. To learn from data, we need to understand it. For large, complex datasets, conventional tools like tables and stats aren't enough. To recognize patterns and anomalies, we must leverage our ability to process information visually.

For this kind of analysis, we need tools that let us quickly test hypotheses and validate assumptions. We need the ability to slice and dice data, switch contexts, zoom, aggregate, focus, and access information as needed. In other words, we need interactivity, flexibility, and speed.

Datagrok delivers exactly that. It works with _millions of columns_ and _billions of rows_ and lets you explore at the speed of thought. Imagine loading the
entire ChEMBL database (2.7 million molecules) in your browser, searching substructures, sketching, filtering, visualizing, and interactively exploring the chemical
space. Datagrok makes it possible.

What's more, Datagrok understands the nature of your
  data, offering actionable insights. It can
  [suggest suitable visualizations](../../../../visualize/view-layout.md#layout-suggestions) for your datasets, automatically render  chemical structures, calculate descriptors, or predict properties. It gives you every tool you need to explore data and uncover its meaning.

![img](../../../../visualize/viewers/img/viewers-interaction-main.gif)

* [Bring data from anywhere](../../../../access/access.md#data-sources). Your data is automatically parsed and rendered in a spreadsheet. The spreadsheet works with _millions of columns_ and _billions of rows_ and has powerful features:
  * Built-in statistics and dataset overview
  * Custom cell renderers, including for domain-specific data (like molecules, sequences, or dose-response curves)
  * [Summary columns](../../../../deploy/releases/platform/1-17.md#summary-columns) and sparklines
  * Editable rows, and [_more_](../../../../visualize/viewers/grid.md).
* [Wrangle data](../../../../transform/transform.md) right from your visualization
  workspace. [Cluster data](../../../../explore/cluster-data.md), [impute missing values](../../../../explore/missing-values-imputation.md), find
  and treat duplicates and outliers.
* Use statistical [functions](../../../concepts/functions/functions.md) to perform calculations.
* Slice and dice data with [50+ interactive viewers](../../../../visualize/viewers/viewers.md). We support all popular
  visualizations (like [scatterplots with built-in regression lines](../../../../visualize/viewers/scatter-plot.md#formula-lines)
  or [box-plots with built-in statistical tests](../../../../visualize/viewers/box-plot.md#t-test)) and certain domain-specific
  viewers. The viewers also support domain-specific value renderers like molecules on scatterplot axes and points.
* Filter, zoom, aggregate, pivot, and cross-link data on the fly. All our viewers are synchronized, [high-performant, and interactive](../../../../develop/under-the-hood/performance.md#viewers).
* Seamlessly access information with [widgets](../../../../visualize/widgets.md)
 and context-driven [info panes](../../../navigation/panels/info-panels.md). 
* Create dashboards in seconds. Share your analysis in easy and
[secure](../../../../govern/access-control/access-control.md) way: send a URL link or integrate: REST API, JS API, or embed as an iframe. 
* Use data annotations and team discussions to collaborate on decision-making. 

Need a specific tool or functionality? Easily [add custom viewers](../../../../develop/how-to/develop-custom-viewer.md) or develop new functions in [R, Python, or Julia](../../../../compute/scripting/scripting.mdx).

[Learn more about capabilities here](../../../datagrok.md).

## Resources

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
    <p class="card-text">
    How do we choose the best location for a new coffee place, given the historical sales data? Datagrok to the rescue! In less than 20 minutes, we achieve the following:<br />
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

* [Viewers](../../../../visualize/viewers/viewers.md)
* [Data science](../data-science.md)
* [Predictive modeling](../../../../learn/learn.md)
