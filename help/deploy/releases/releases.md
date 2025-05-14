---
title: "Release highlights"
format: mdx
keywords:
  - release history
  - release highlights
  - release summary
---

See details: [issues](release-history.md), [plugin changelogs](plugins/plugins.mdx), [JS API compatibility](compatibility/compatibility.mdx)

## 1.25 || 2025-Mar

### Navigation and usability 

|<div style={{ width:220 }}></div>| <div style={{ width:500 }}></div> |
|----------------- | -----------------------------------|
| **Smarter file viewers**: Previously, viewer selection was based solely on file format (e.g., XLSX). Now, Datagrok also analyzes file content to choose the appropriate viewer.<br/><br/>For example, XLSX files with regular data open in the [table view](../../datagrok/navigation/views/table-view.md), while XLSX files with plate-like content open using the [plate reader](https://github.com/datagrok-ai/public/tree/master/packages/Curves#plate-readers) | ![Smarter file viewer](img/1.25-file-viewer-xlsx-format.gif) |

### Developer updates

* **Breaking changes**: Site markup, UI API, and CSS styles have been updated. If your packages or integrations rely on specific UI elements or styling, it is recommended to check for compatibility and make the necessary adjustments.

## 1.24 || 2025-Feb

### Visualization and analysis

|<div style={{ width:220 }}></div>| <div style={{ width:500 }}></div> |
|----------------- | -----------------------------------|
|**Grid**: Rows and columns now resize together when you change row heights. This creates a zoom effect, where cell content adapts for optimal viewing. Zooming out shows the big picture. Zooming in shows details. This works especially well with dynamic content like proteins, user profiles, [forms](../../visualize/viewers/grid.md#summary-columns), and JIRA tickets. |![grid dynamic column width](img/grid-dynamic-column-widths.gif)|
|**Grid**: We now support right-click panning. Hold the right mouse button to drag and move around your data. This feature is particularly useful when exploring extensive datasets in the zoomed-out mode.|![grid mouse pan](img/grid-mouse-pan-content.gif)|
|**Grid**: You can now create [column groups](../../visualize/viewers/grid.md#column-groups)|![grid column groups](img/grid-column-groups.gif)|
|**Grid**: You can now easily [(un)hide columns](../../visualize/viewers/grid.md#hiding-and-unhiding-columns)|![grid unhide columns](img/grid-unhide-columns.gif)|
|**Cell renderers**: [Tags](../../visualize/viewers/grid.md#cell-renderers) cell renderer now works with boolean columns (previously, it required comma-separated values)|![grid tags renderer](img/grid-tags-renderer.gif)|
|**Cell renderers**: [Sparklines](../../visualize/viewers/grid.md#summary-columns) now support multiple normalization options - by row, column, or global values|![sparklines normalization](img/sparklines-normalization-demo.gif)|
|**Cell renderers**: We now support adaptive rendering of files in cells|![grid-files in cells](img/grid-files-in-cells.gif)|

### Developer updates

**JS API changes**:

* New features:
  * Classes: 0
  * Functions: 10
* [Breaking changes](http://localhost:3000/help/deploy/releases/compatibility/?from=1.23.1&to=1.24.0&breakingchanges=true)


## 1.23 || 2024-Dec

### New tools and apps

|<div style={{ width:220 }}></div>| <div style={{ width:500 }}></div> |
|----------------- | -----------------------------------|
| **Pyodide**:<br/><br/> Run Python in the browser | <iframe src="https://www.youtube.com/embed/mphWv6yfEkc?start=1470" width="512" height="288" frameborder="0" allowfullscreen></iframe>  |
| **Admetica**:<br/><br/> Predict, visualize, and analyze 20+ ADMET properties using Datagrok's open-source ADMET predictor<br/><br/>[Learn more](https://github.com/datagrok-ai/admetica)| <iframe src="https://www.youtube.com/embed/mphWv6yfEkc?start=2100" width="512" height="288" frameborder="0" allowfullscreen></iframe>|


### Data access and management

|<div style={{ width:220 }}></div>| <div style={{ width:500 }}></div> |
|----------------- | -----------------------------------|
| New **Visual Query Builder**:<br/><br/>Seamlessly join tables, aggregate, and pivot from one user interface| <iframe src="https://www.youtube.com/embed/mphWv6yfEkc?start=87" width="512" height="288" frameborder="0" allowfullscreen></iframe> |

### Visualization and analysis

|<div style={{ width:220 }}></div>| <div style={{ width:500 }}></div> |
|----------------- | -----------------------------------|
| **Grid**:<br/><br/>Group multiple columns together and apply style for multiple columns at once| ![Grid grouped columns](../../visualize/viewers/img/grid-grouped-columns.gif)   |
| **Scatterplot**:<br/><br/><li>Interactively explore tens of millions of data points</li><li>Customize how labels appear and behave</li><li>Define columns for scatterplot whiskers to highlight data distribution</li> | <iframe src="https://www.youtube.com/embed/mphWv6yfEkc?start=3370" width="512" height="288" frameborder="0" allowfullscreen></iframe>|


### Developer updates

**JS API changes**:

* New features:
  * Classes: 1
  * Functions: 0
* [Breaking changes](https://datagrok.ai/help/deploy/releases/compatibility/?from=1.22.0&to=1.23.0&breakingchanges=true)

## 1.22 || 2024-Nov

### Data access and management

|<div style={{ width:220 }}></div>| <div style={{ width:500 }}></div> |
|----------------- | -----------------------------------|
|**Custom identifier registration**:<br/>Register custom identifiers (e.g., compound IDs) to search, link, and analyze entity data across the platform<br/><br/>[Video demo](https://www.youtube.com/watch?v=4_NS3q7uvjs&t=2932s)<br/>[Developers: Learn more](https://datagrok.ai/help/develop/how-to/grid/register-identifiers)|![custom identifiers](img/register-idenitifiers.gif)|
|**Search integrated functions**:<br/>Annotate queries and other functions with search patterns to display results when users search with matching terms<br/><br/>[Learn more](../../datagrok/concepts/functions/func-params-annotation.md#search-integrated-functions) |![Search-integrated queries](../../datagrok/concepts/functions/search-integrated-queries.gif)|
|**Plugin databases**|You can now ship a Postgres database (such as chemical registration system) with your plugin|

### Visualization and analysis

#### Scatterplot

|<div style={{ width:220 }}></div>| <div style={{ width:500 }}></div> |
|----------------- | -----------------------------------|
|New data series connection lines let you track datapoint relationship over time<br/><br/>[Learn more](../../visualize/viewers/scatter-plot.md#connecting-lines) |![](../../visualize/viewers/img/scatter-plot-lines.png) |
<!---|Set color coding using expressions| [img] |--->


### Platform performance

* **[Pyodide](https://pyodide.org/en/stable/) integration**: you can now run Python functions, including data transformation steps, directly in the browser
* **Improved automated data cleanup**: keeps the server lean by deleting old data and logs (could be configured)  

<!----

### Developer updates

#### JS API changes

-->

## 1.21 || 2024-Sep

:::warning

This version requires a database migration that cannot be rolled back. Once you upgrade to version 1.21.0, you cannot downgrade. If you need assistance, contact Datagrok Support.

:::

### Data transformation

We improved the calculated column feature (requires a [PowerPack](https://github.com/datagrok-ai/public/tree/master/packages/PowerPack) package):
* **Usability**: Drag-and-drop, search, and function suggestions in **Add New Column** dialog
* **Formula assistance**: Real-time validation, syntax highlighting, autocompletion hints
* **Extended functions**: Added package-specific function support

[Video walkthrough](https://www.youtube.com/watch?v=4_NS3q7uvjs&t=1655s).

### ML & modeling

We continue to improve [ML capabilities](../../learn/learn.md) within the platform. We added:

* New interactive model training dashboard
  * Binary classification support with prediction threshold control
  * AUC-ROC visualization
  * Model comparison tools
  * Automatic model selection
* New algorithms ([EDA package](https://github.com/datagrok-ai/public/tree/master/packages/EDA)):
  * Softmax classifier 
  * XGBoost tools
  * Partial least squares regression

### Other

* File caching is now enabled by default. It works automatically without the need for configuration.

<!--

### Developer updates

#### JS API changes

---->

## 1.20 || 2024-Jul 

### New tools and apps

|Tool / App<div style={{ width:100 }}></div> | Purpose| Requirements <div style={{ width:150 }}></div>|Learn more<div style={{ width:150 }}></div>|
|-----------|-------------|---------------|--------|
| Dimensionality reduction using WebGPU | Run dimensionality reduction on 100K molecules in less than 10 seconds | [EDA package](https://github.com/datagrok-ai/public/tree/master/packages/EDA) | [Video walkthrough](https://www.youtube.com/watch?v=RS163zKe7s8&t=2648s) <br/>[Wiki page](../../explore/dim-reduction.md)|
| Diff Studio | App that turns differential equations into interactive visual models |[Diff Studio package](https://github.com/datagrok-ai/public/tree/master/packages/DiffStudio) |[Launch app](https://public.datagrok.ai/browse/apps/DiffStudio/Templates/basic?params:Initial=0.01&Final=15&Step=0.01&Y=0)<br/>[Wiki page](https://github.com/datagrok-ai/public/blob/master/packages/DiffStudio/README.md)|

###  ML & modeling

We continue to improve [ML capabilities](../../learn/learn.md) within the platform:

Improved in-platform model training:
* Automatic input substitution for models
* Interactive hyperparameter tuning
* New interactive visualization widgets for model performance evaluation

[MLflow](https://mlflow.org/) integration ([learn more](../../learn/mlflow.md)):
* Import and apply MLflow predictive models directly to Datagrok datasets
* Streamlined model management:
   * Automatic MLflow model fetching
   * Input annotation via MLflow tags
* MLflow models inference support

### Platform performance

We've introduced caching for [file shares](../../access/files/files.md) (including S3, Azure, and local storage):
* Faster file retrieval and display in [Browse](../../datagrok/navigation/views/browse.md)
* Entirely browser-based cache, eliminating server requests
* Configurable caching for individual files and file shares
* Automatic cache invalidation when files are modified within the platform

[Learn more](../../develop/how-to/functions/cache-function-results.md).


## 1.19 || 2024-Jun

### Visualization and analysis

|<div style={{ width:220 }}></div>|  |
|----------------- | -----------------------------------|
|[Organize and filter hierarchical data](../../visualize/viewers/filters.md#expression-filter):<br/><li>Group and navigate categories using a tree</li><li>Toggle selections with checkboxes or shortcuts</li><li>Customize hierarchy by rearranging, adding, or removing columns</li> |![Hierarchical Filter](../../uploads/gifs/hierarchical-filter.gif) |
|[Create custom filters using expressions](../../visualize/viewers/filters.md#expression-filter):<br/><li>Use column-specific operations like `>`, `contains`, or [regex](https://en.wikipedia.org/wiki/Regular_expression)</li><li>Combine multiple conditions with AND/OR logic</li>|![](../../uploads/gifs/expression-filter.gif) |
|Highlight text matches in cells with [free-text mode](../../visualize/viewers/filters.md#free-text-filter-mode) |![](../../uploads/gifs/free-text-filter.gif) |

<!-----

## Grid improvements

|<div style={{ width:220 }}></div>|  |
|----------------- | -----------------------------------|
|Customize cell styles:<li>content and header</li><li>color and font</li><li>alignment</li><br/>...and more! | |

TODO: I think the realization is buggy

----->

### Other

* Added new connectors: [Azure Blob](../../access/files/shares/azure.md), [SharePoint](../../access/files/shares/sharepoint.md)
* Enhanced CSV export options for dataframes and views:
  * Preserve row and column order when exporting
  * Export molecules as SMILES strings
  * Selectively export visible columns only
  * Configure column qualifiers for improved data formatting

### Developer updates

You can now integrate with Datagrok using Datagrok's [REST API](../../develop/packages/rest-api.md):

* Supports operations on files, tables, dashboards, and functions
* Allows for programmatic data management, dashboard creation, and function calls
* [Python client library](https://github.com/datagrok-ai/public/tree/master/python-api) available for simplified API interaction


## 1.18 || 2024-Mar

### New tools and apps

|Tool / App<div style={{ width:100 }}></div> | Purpose| Requirements <div style={{ width:150 }}></div>|Learn more<div style={{ width:150 }}></div>|
|-----------|-------------|---------------|--------|
| Matched Molecular Pairs | Tool for analyzing how structural changes affect molecular properties and activity | [Chem package](https://github.com/datagrok-ai/public/tree/master/packages/Chem) | [Interactive demo](https://public.datagrok.ai/browse/apps/Tutorials/Demo/Cheminformatics/Matched-Molecular-Pairs) <br/>[Wiki page](../../datagrok/solutions/domains/chem/chem.md#matched-molecular-pairs)|
| Docking| Tool for screening ligand libraries against AutoDock-prepared targets | [Docking package](https://github.com/datagrok-ai/public/tree/master/packages/Docking)|[Interactive demo](https://public.datagrok.ai/browse/apps/Tutorials/Demo/Bioinformatics/Docking)<br/>[GitHub](https://github.com/datagrok-ai/public/tree/master/packages/Docking)<br/><!--[Wiki page]-->|
| Hit Design|App for collaborative hit design | [HitTriage package](https://github.com/datagrok-ai/public/tree/master/packages/HitTriage) |[Launch app](https://public.datagrok.ai/apps/HitTriage/HitDesign)<br/>[GitHub](https://github.com/datagrok-ai/public/blob/master/packages/HitTriage/README_HD.md)<br/><!--[Wiki page]--->|
| Hit Triage |App for collaborative hit assessment and prioritization |[HitTriage package](https://github.com/datagrok-ai/public/tree/master/packages/HitTriage) |[Launch app](https://public.datagrok.ai/apps/HitTriage/HitTriage)<br/>[GitHub](https://github.com/datagrok-ai/public/blob/master/packages/HitTriage/README_HT.md)<br/><!--TODO[Wiki page]--->|

### Navigation and usability

|<div style={{ width:200 }}></div>| <div style={{ width:400 }}></div> |
|----------------- | -----------------------------------|
|New **Browse** view:<br/><br/>Navigate the tree to access, preview, and manage anything in Datagrok, all from one convenient location. <br/><br/>[Learn more](../../datagrok/navigation/views/browse.md) | ![](img/release1.18-browse.gif) |

### Data access and management

|<div style={{ width:100 }}></div>|  |
|----------------- | -----------------------------------|
|Add custom metadata to anything - molecules, experiments, or users <br/><br/>[Learn more](../../govern/catalog/sticky-meta.md) | ![](img/1.18-sticky-meta.png)  |

### Visualization and analysis

#### Cell renderers

|<div style={{ width:200 }}></div>|  |
|----------------- | -----------------------------------|
|New cell renderers:<br/><br/><li>Dropdown</li><li>MultipleChoice</li><li>Tags</li><br/><br/><!--[Learn more](../../visualize/viewers/grid.md#cell-renderers)-->|![image](img/release1.18-cellrend-tags-multichoice-dropdown.gif)  |

#### Trellis plot

|<div style={{ width:200 }}></div>|  |
|----------------- | -----------------------------------|
|[Trellis plot](../../visualize/viewers/trellis-plot.md) now supports sparklines |![](img/release1.18-trellis-sparklines.gif) |

<!-- //TODO after patch

### New text filter

### Ability to clone projects

|<div style={{ width:200 }}></div>|  |
|----------------- | -----------------------------------|
|You can now clone projects:<br/><br/><li>With data sync</li><li>...</li> |[IMG] |

-->

### Developer updates

* Computation queue
* Harmonized inputs and properties
* Major scripting Platform performance improvements
* CVM: CUDA support
* Push notifications
* Refactor package update mechanism
* JS API for viewer events

## 1.17 || 2023-Oct

### Visualization and analysis

#### Grid: Pinned rows

|<div style={{ width:100 }}></div>|  |
|----------------- | -----------------------------------|
|Pin rows and save pinned view to the layout or project<br/><br/>[Learn more](../../visualize/viewers/grid.md#pin-rows-and-columns) |![image](img/pinned-rows.gif)  |

<!--
<table>
  <thead>
    <tr>
      <th width="30%"></th>
      <th width="70%"></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Pin rows and save pinned view to the layout or project.<br/><a href="https://datagrok.ai/help/visualize/viewers/grid#pinned-rows">Learn more</a></td>
      <td><image src="img/pinned-rows.gif" alt="image"></image></td>
    </tr>
  </tbody>
</table>
--->

#### Grid: Column selection

|  |        <div style={{ width:100 }}></div>     |
|-------------------------|-------------------|
|Hold <kbd>Shift</kbd> and drag the mouse over the column headers|![image](img/column-selection.gif)|

#### Grid: Column navigation

|                         | <div style={{ width:100 }}></div> |
|-------------------------------------|-------------------|
|In the **Column Manager**:<br/><li>Hover over a column to automatically scroll to it and see its details</li><li>Filter columns based on data or semantic type</li>|![Column preview and filter by type](img/column-preview-filter-by-type.gif)|

#### Grid: Summary columns

|                         |  <div style={{ width:400 }}></div> |
|-------------------------|-------------------|
|<li>[Smart forms](https://community.datagrok.ai/t/powergrid-smartform/774/1) show values from multiple columns in one cell </li><li>Adaptive rendering to fit the cell size</li><li>Values inherit color-coding from source columns</li><li>Cell renderers for different data and semantic types</li><br/>(From [PowerGrid package](https://github.com/datagrok-ai/public/blob/master/packages/PowerGrid/README.md))|![](img/smart-forms.gif)|

#### Grid: Linked tables in cells

|                         |<div style={{ width:400 }}></div>|
|-------------------------|-------------------|
|Link tables and show data from one linked table in another|![Linked table data in cell](img/grid-nested-linked-tables.gif) |

<!--
### Visualize data behind grouped rows

|                         |                   |
|-------------------------|-------------------|
|**Viewers**:<br/>[Pivot] now supports the visualization of raw data behind grouped rows|![Viewers within pivot](pivot-viewers.gif) |

-->

#### Line chart: Split by category

| <div style={{ width:125 }}></div> |             |
|-------------------------|-------------------|
|Split line charts by multiple categorical columns| ![Splitting line chart by category](img/line-chart-split-by-columns.gif)|

#### New info plane: Plots

|                         |<div style={{ width:150 }}></div>|
|-------------------------|-------------------|
|Automatically visualizes selected columns for quick profiling<br/>(**Context Panel** > **Plots**)<br/><br/>[Learn more](https://community.datagrok.ai/t/ux-updates/544/5?u=oahadzhanian.datagrok.ai)|![Plots info pane](img/plots-info-pane.gif)|
 
#### New info pane: Content

|                         | <div style={{ width:250 }}></div>  |
|-------------------------|-------------------|
|Displays details for the selected rows in a tabular format (**Context Panel** > **Content**)|![Content info pane](img/content-info-pane.gif)|

### Platform performance

We now support client and server-side [caching for function results](../../develop/how-to/functions/cache-function-results.md). This feature is
particularly useful for functions that produce consistent outputs, like queries
and scripts. The client-side cache, limited to 100 MB or 100,000 records, speeds
up data access and improves network efficiency. The unrestricted server-side
cache improves response times and overall server performance. Together, these
caching mechanisms provide a smoother and more responsive platform experience.

### Platform administration

* New installation wizard for platform configuration during deployment
* Authorization with [Okta](https://www.okta.com/)
* Amazon CloudWatch log export improvements

### Developer updates

#### Auto-fill function parameters

|                         | <div style={{ width:550 }}></div>  |
|-------------------------|-------------------|
|Auto populate function input parameters based on a chosen key, such as selecting a car model to instantly fill in details like its mileage and engine size<br/><br/>[Learn more](https://datagrok.ai/help/datagrok/concepts/functions/func-params-annotation#lookup-tables)| ![image](../../datagrok/concepts/functions/default-values-from-file.gif)|