---
title: "Release history"
slug: develop/admin/releases/release-history
position: 6 # float position is supported
---

## Latest version

| Service                                                   | Docker Image                                                                                      |
|-----------------------------------------------------------|---------------------------------------------------------------------------------------------------|
| [Datagrok](../../develop/under-the-hood/infrastructure.md#datagrok-components)      | [datagrok/datagrok:1.18.6](https://hub.docker.com/r/datagrok/datagrok)                            |
| [Grok Connect](../../access/access.md#data-connection) | [datagrok/grok_connect:2.1.10](https://hub.docker.com/r/datagrok/grok_connect)                    |
| Grok Spawner                                              | [datagrok/grok_spawner:1.4.8](https://hub.docker.com/r/datagrok/grok_spawner)                     |
| [Grok Compute](../../develop/under-the-hood/infrastructure.md#grok-compute)         | [datagrok/grok_compute:1.5.5](https://hub.docker.com/r/datagrok/grok_compute)                     |
| [Jupyter Kernel Gateway](../../compute/scripting/scripting.mdx)   | [datagrok/jupyter_kernel_gateway:1.6.2](https://hub.docker.com/r/datagrok/jupyter_kernel_gateway) |
| [Jupyter Notebook](../../compute/jupyter-notebook.md)  | [datagrok/jupyter_notebook:1.1.1](https://hub.docker.com/r/datagrok/jupyter_notebook)             |
| [H2O](../../develop/under-the-hood/infrastructure.md#h2o)                           | [datagrok/h2o:1.1.1](https://hub.docker.com/r/datagrok/h2o)                                       |
| [CVM Nginx](../../develop/under-the-hood/infrastructure.md#load-balancer)           | [datagrok/cvm_nginx:1.10.0](https://hub.docker.com/r/datagrok/cvm_nginx)                          |

See also:
- [Versioning policy](../../develop/dev-process/versioning-policy.md)
- [Docker-Compose](../docker-compose/docker-compose.md)


## 2024-07-23 Datagrok 1.20.0 release 

### Breaking changes:
 
* File caching updates.
* Inputs harmonization.
 
### Visualization and usability improvements:

* GROK-15886: MlFlow registry integration.
* GROK-15881: Lightweight Predictive modeling improvements.
* GROK-15868: Browse hierarchy updates.
* GROK-15906: Harmonized layout editing for Data Query.
* GROK-15908: Forms autosize.
* GROK-15322: Macromolecules updates.
* GROK-15756: Diff Studio: Debugging mode.
* GROK-15885: Predictive models: pretrain on small batch of data. S
* GROK-14605: Context Pane: add visibility of the column names. 
* GROK-16206: Column and table metadata: register commonly used tags and expose them to UI.
* GROK-15866: Projects: restyling dragging an entity into a project.
* GROK-15009: Context Panel: make a tooltip explaining why the property is inactive.
* GROK-15359: Changing logic for docking arrows. 
* GROK-15427: Add the 'Add to favorites' command to objects.  
* GROK-16274: Prevent auto running of functions in browse preview.
* GROK-16278: Avoid sending extra data in files serialization. 
* [#2847](https://github.com/datagrok-ai/public/issues/2847): Export CSV without rounding.
 
### Viewers
#### Improvements:
* GROK-14382: Add full property panel menu to the context menu on right click.
* GROK-15008: Density plot: add functionality.
* GROK-14607: Network diagram: auto layout.
* [#2658](https://github.com/datagrok-ai/public/issues/2658): Scatter Plot: ColorMin, ColorMax, SizeMin, SizeMax properties.
* [#2926](https://github.com/datagrok-ai/public/issues/2926): Heatmap: save scrollbar position for layout.
* [#2932](https://github.com/datagrok-ai/public/issues/2932): Hide the option “Is Grid” from the property panel. 
  
#### Fixes:
* GROK-11678: Tile viewer: switching between tables doesn't work.
* GROK-15518: Pivot table, Pie chart: filtration issues.
* GROK-15811: Pie Chart legend: Fix click on the category black cross. 
* GROK-16226: DG.Filters sets columnNames incorrectly. 
* GROK-15626: Viewer settings do not reset after switching between different options. 
* GROK-16330: Table view: Form style with refresh button is broken. 
* [#2925]( https://github.com/datagrok-ai/public/issues/2925): Viewer title is not saved in layout if it was edited in viewer's header. 
* [#2771]( https://github.com/datagrok-ai/public/issues/2771): Filter out legend with empty categories in the legend. 


### Grid
#### Improvements:
* GROK-16205: Support markup links. 
* GROK-15984: Formula-based rendering.
* [#2789](https://github.com/datagrok-ai/public/issues/2789): Ability to display the table name/description on the Grid.
#### Fixes:
* GROK-15188: Dataframe/grid: renaming column from '~name' -> 'name' doesn’t change visibility of column. 
* [#2934](https://github.com/datagrok-ai/public/issues/2934): Some hotkeys for selection/deselection are not working for columns (MacOS).
 
### [Scatterplot](../../visualize/viewers/scatter-plot.md)
#### Improvements and fixes:
* GROK-14632: Marker legend isn't rendered in the trellis plot. 
* GROK-14376: Scatterplot: unzooming with mouse scroll can go to infinity.
* [#2913](https://github.com/datagrok-ai/public/issues/2913): Improve different context menu style / functionality between Scatterplot in Trellis vs native Scatterplot.
 
 
 
### [Line Chart](../../visualize/viewers/line-chart.md)
#### Fixes:
* GROK-15466: Table > Add View: Formula Lines for Line Chart should be present in added view.
* GROK-15464: Receiving errors when setting the second split in some cases. 
* GROK-15485: Bugs in RowSource = FilteredSelected, SelectedOrCurrent.
* [#2875](https://github.com/datagrok-ai/public/issues/2875): Line chart, trellis plot: number of selected columns is not updated in the properties panel immediately. 
* [#2904](https://github.com/datagrok-ai/public/issues/2904): When multi axis option is used. together with split, line chart is empty. 
 
 
### Bar chart
#### Improvement:
* [#2881](https://github.com/datagrok-ai/public/issues/2881): Bar chart: add an option to rotate / invert axes. 
#### Fixes:
* GROK-15887: Viewer window can't be resized in height when it's too high on the screen. 
* GROK-16234: Wrong label shortening. 
* GROK-16331: Molecules are rendered incorrectly on resizing bar chart. 
 
 
### Box plot
#### Improvements and a fix:
* [#2905](https://github.com/datagrok-ai/public/issues/2905): Add ''Size by" or "Shape by" or "color by'' options. 
* GROK-9748: Violin-state option. 
* GROK-16244: Violin style: category bug. 
 
 ### JS API
#### Improvements:
* [Compatibility tool](https://datagrok.ai/help/deploy/releases/compatibility/) for JS API changes.
* GROK-15771: RFC: Public API. 
* GROK-15827: Expose core viewer events to JS API. 
* GROK-16215: Force DG.FileInfo to have `name`. 
* GROK-16240: Replace metadata set by tags with column.meta.xxx. 
 #### Fixes:
* GROK-15671: `format` option does not go to `DG.Property.format`. 
* GROK-16228: `grok.functions.calls.list` returns incorrect Func's nqName functions. 

### Data Access
#### Improvements:
* GROK-15413: Sharepoint data provider. 
#### Fixes:
* GROK-15979: Swagger: blank screen after loading swagger file. 
* GROK-16003: Queries: New Aggregation Query: PivotGridDbTable. 
* GROK-16248: Proxy storage writeFileStream method not implemented. 
* GROK-15967: Context panel: When you share the connection \- the information doesn`t update.
 
 
### Enhancements in packages
#### Improvements:
* GROK-15922: EDA: Linear Regression. 
* GROK-15759: PowerPack: 'Add new column': Formula edit improvements.
* GROK-11494: Widgets: events.
* [#2902](https://github.com/datagrok-ai/public/issues/2902): Create "Pivot table" tutorial part.
#### Fixes:
* GROK-15324: Widgets: Recent projects: context menu fixing. 
* [#2924](https://github.com/datagrok-ai/public/issues/2924): Curves: connectDots doesn't work if fit is turned off. 

 
### [Bio](https://github.com/datagrok-ai/public/blob/master/packages/Bio/CHANGELOG.md)
#### Fixes:
* GROK-15525: MSA: add checking unsuitable data to avoid running MSA with them. 
* GROK-15996: Fix macromolecule cell renderer fails in long mode. 
* GROK-15798: Fix to Atomic Level for units FASTA, UN alphabet.
* GROK-15793: Calculate Identity, Similarity throw Index out of bounds.
* GROK-15994: Helm: Color missing monomers. 
* GROK-15995: Helm: Colors for libraries monomers. 
* GROK-15796: Fix to Helm cell renderer for conversion to Helm. 

 ### [Chem](https://github.com/datagrok-ai/public/blob/master/packages/Chem/CHANGELOG.md)
#### Fixes:
* GROK-15299: Functions can't be saved to a project with data sync. 
* GROK-15790: Transform to smiles doesn't change units. 
 
### Scripting
#### Improvements and fixes:
* GROK-15424: Cascade parametrized inputs in scripts. 
* GROK-15876: Split logic for blobs and file parameters. 
* GROK-15973: Behavior changes for clicking the 'Save' button. 
* GROK-15968: Usage: rename 'Run count' into 'Successful Runs'. 
* GROK-15829: Unable to add parameter to script: "+" button is missing. 
 
### Fixes:
* GROK-15834: Table Linking: wrong results when multiple duplicate keys are present in both master and details table. 
* GROK-15871: Projects: Project creation date is shown instead of updated one in the dashboards. 
* GROK-15570: Fixes for slider for filtered numerical column where all values are nulls. 
* GROK-12320: "Add new row" button is not rendered properly when no row headers are used.  
* GROK-15672: ML menu: old menu way 'ML > Tools > Missing' is available. 
* GROK-16208: Tooltips: broken layout for a float column with only nulls. 
* GROK-16238: Docker: Failed package deploy with docker config. 
* GROK-16269: History prints wrong timestamps. 
* GROK-15490: Fix UI issues with scaling on different browser sizes.
* GROK-16320: Predictive modeling: interactive model saves without a name.
* GROK-15985: F2 / Column properties: Tags: value is not saved unless you move the cursor out of the editor. 



## 2024-06-21 1.19.1

### Fixes:

* GROK-16000: HashId for package functions doesn't work.


## 2024-06-07 Datagrok 1.19.0 release
 
### Breaking changes:
 
* Updated project`s management.
* Power Search: updated support for entities.
* Usage analysis: user report system improvements.
* Speedup for platform’s startup.

### Visualization and usability improvements

* New hierarchical and expression filters.
* New recent activity view.
* Categorize function improvements.
* Improved experience with sticky meta.
* Browse updates: search, preview, tooltips. 
* Viewers: hooks for custom JS initialization.
* OptiTool: compute inputs providing the specified output of RichFunctionView functions.
* GROK-15361: Updates for Custom formatting.
* GROK-14203: Calculated columns: Update values if the base column values changed.
* [#2673](https://github.com/datagrok-ai/public/issues/2673): PowerGrid: MultiChoice cell renderer.
* [#2681](https://github.com/datagrok-ai/public/issues/2681): PowerGrid: Tags cell renderer.
* [#2765](https://github.com/datagrok-ai/public/issues/2765): PowerPack: integrate search with the semantic value extractors.
* [#2668](https://github.com/datagrok-ai/public/issues/2668): Implement export selected/filtered rows as a CSV.
* [#2710](https://github.com/datagrok-ai/public/issues/2710): Implement export selected dataframe/view as a CSV.
* [#2772](https://github.com/datagrok-ai/public/issues/2772): Pivot table: add a property to hide sparklines.
* [#2847](https://github.com/datagrok-ai/public/issues/2847): Export CSV without rounding.
 
### Viewers
#### Improvements:
* [#2733](https://github.com/datagrok-ai/public/issues/2733): Legend: provide diverse palette and markers' shape.
* [#2755](https://github.com/datagrok-ai/public/issues/2755): Hide the zoom slider for making screenshots.
* [#2742](https://github.com/datagrok-ai/public/issues/2742): Indicate the selected column in the column selector.
* [#2770](https://github.com/datagrok-ai/public/issues/2770): Ability to toggle on / off multiple vs single regression line per category.
* [#2664](https://github.com/datagrok-ai/public/issues/2664): Implement a cascade column renaming used by Formula lines.
* [#2685](https://github.com/datagrok-ai/public/issues/2685): Histogram viewer: add stacked view option when split is set.
#### Fixes:
* [#2642](https://github.com/datagrok-ai/public/issues/2642): Filtering done using viewers is unexpectedly reset on applying other filters.
* [#2771](https://github.com/datagrok-ai/public/issues/2771): Filter out legend with empty categories in the legend.
 
### Grid
#### Improvements:
* GROK-15475: Ability to customize cell styles (content/header, color, font, alignment etc.).
* GROK-15469: Text cells: highlight the search terms.
#### Fixes:
* GROK-15249: Error when pressing Column Properties on a programmatically created Grid.
* GROK-15457: Fix categorical color-coding styles in the property panel.
* [#2682](https://github.com/datagrok-ai/public/issues/2682): Fix the tooltip for Allow Row Reordering.
* [#2663](https://github.com/datagrok-ai/public/issues/2663): Enable easier composition of rows/columns range selections.
 
### Forms
#### Improvements:
* GROK-14951: Clicking on a field should drive to the appropriate column in the table (like 3DX).
* GROK-14985: Ability to edit a form in the form designer.
* GROK-15416: Integration with sparklines.
* [#2687](https://github.com/datagrok-ai/public/issues/2687): Add an option to adjust the size of the compound structure image.
#### Fixes:
* GROK-15209: Error after renaming a column name.
* GROK-15325: In case form viewer is added through layout, grid doesn't jump to corresponding cell when clicking on form viewer field.
* GROK-15677: Fixing "Build a form for selected columns" option.
* [#2644](https://github.com/datagrok-ai/public/issues/2644): Next row button: incorrect behavior after sorting the grid.
 
### [Scatterplot](../../visualize/viewers/scatter-plot.md)
#### Improvements:
* GROK-15459: Ability to show horz axis labels at 45 degrees.
* [#2704](https://github.com/datagrok-ai/public/issues/2704): Add min/max fields for X/Y axes on plot to improve viewer's usability.
#### Fixes:
* GROK-15069: Incorrect hadling of Color and Markers legend in some cases.
* GROK-15190: Scatterplot, line chart: updating formula lines does not work for those specified via a dataframe.
* GROK-13884: The legend is displayed differently depending on the order in which the options are assigned.
* [#2671](https://github.com/datagrok-ai/public/issues/2671): Filtered out data is unexpectedly shown in some cases on filtering the plot further using markers legend.
 
### [Trellis plot](../../visualize/viewers/trellis-plot.md)
#### Improvements:
* GROK-15393: Ability to scroll via mouse wheel.
* [#2636](https://github.com/datagrok-ai/public/issues/2636): Bar chart: improvements.
#### Fixes:
* GROK-15035: Sparklines: incorrect data display.
* GROK-15036: Histogram: error when Row source = MouseOverGroup.                                      	
* GROK-15038: Box plot: error when selecting rows in grid.
* GROK-15013: Bar chart: numbers (zeroes) on chart have different sizes.
* GROK-15033: Line chart: if you set the Split option, the legend remains visible on other inner viewers.
* GROK-15025: Some axes labels overlap the plot.
* GROK-15137: Sparklines: error when setting selected columns.
* GROK-15304: Selecting categories in the legend works incorrectly.
* [#2749](https://github.com/datagrok-ai/public/issues/2749): Line chart and Trellis plot: '+' icon is missing.
 
### [Line Chart](../../visualize/viewers/line-chart.md)
#### Improvements:
* [#2846](https://github.com/datagrok-ai/public/issues/2846): Add distinguishable color legend for line chart when there are multiple Y axes and split key.
* [#2726](https://github.com/datagrok-ai/public/issues/2726): Implement consistent behavior for aggregated and non-aggregated axes.
#### Fixes:
* GROK-15078: Error on changing column name when column used in x\- and y columns.
* GROK-14501: Aggregation Type and Show Split Selector do not work.
* GROK-15048: Incorrect legend display when Row Source = MouseOverGroup and Filter is set.
* [#2782](https://github.com/datagrok-ai/public/issues/2782): Incorrect color-coding while setting filters.
* [#2776](https://github.com/datagrok-ai/public/issues/2776): Page is frozen if too many categories are used in split.
* [#2747](https://github.com/datagrok-ai/public/issues/2747): Formula lines: inconsistent behavior.
* [#2677](https://github.com/datagrok-ai/public/issues/2677): Formula line position is incorrect after switching axis to log scale.
 
### Data Access:
* GROK-14953: Azure blob: directory shown as a file.
* GROK-15397: Storage: add SharePoint provider.
* GROK-15434: Run polytool on molfile v3k libs.
 
### JS API:
* GROK-14512: JS API for entity properties/types/schemas.
* GROK-15671: `format` option does not go to `DG.Property.format`.
* GROK-14674: DG.ViewBase does not trigger UI update on name change.
* GROK-14707: DG.Column.init contructs invalid DataFrame column.
 
### Enhancements in packages
 
### [Bio](https://github.com/datagrok-ai/public/blob/master/packages/Bio/CHANGELOG.md)
#### Improvements:
* GROK-15039: Expose function for PolyTool lib creation.
* GROK-15044: Add basic documentation for monomer lib file manager UI.
* GROK-15045: Monomer lib file manager UI improvements.
* [#2707](https://github.com/datagrok-ai/public/issues/2707): Monomer type as a class.
#### Fixes:
* GROK-14913: Atomic Level for sequences of fasta MSA with gaps fixing.
* GROK-14916: Biosubstructure filter for sequences of Helm fixing.
* GROK-11982: Composition Analysis: duplicates a WebLogo viewer when applying a layout.
* GROK-15086: GetRegion, Notation: doesn't work.
* GROK-15150: hidden/showed inputs incorrect display in forms.
* GROK-15176: An Error after pressing the Molecule on molfile(msa(Sequence)) column.
* GROK-14914: Helm: Fix Properties panel of Context panel for sequences of MSA.
 
### [Chem](https://github.com/datagrok-ai/public/blob/master/packages/Chem/CHANGELOG.md)
#### Improvements:
* GROK-15018: RGroup analysis improvements.
* GROK-15508: MMPA | Generated \- show what cases are actually generated.
* GROK-15501: MMPA | UI/UX and other improvements.
* GROK-15502: MMPA | Ability to switch activity on/off.
* GROK-15454: Similarity, Diversity search viewers | Implement row source.
* GROK-14773: RdKit service as molecule cache.
* [#2729](https://github.com/datagrok-ai/public/issues/2729): Similarity, Diversity searches: add the ability to perform the search for the filtered data.
* [#2525](https://github.com/datagrok-ai/public/issues/2525): Scaffold  tree highlighting improvements.
#### Fixes:
* GROK-14029: Substructure highlight: settings are not saved to the Layout.
* GROK-15063: 2D widget is not resizing.
* GROK-15144: Layout applies incorrectly when having several clones.
* GROK-15127: Python scripts Mutate, Curate, R groups sporadically fail.
* GROK-15171: Filter doesn't apply correctly after removing rows and disabling.
* GROK-15208: MMPA | 'Refreshing pairs...' progress bar is not closing. 
* GROK-15240: Error in balloon after R-group analysis' start (malformed data).
* GROK-15507: MMPA | Control of fragments size (to options).
* GROK-12906: Context Pane: Fix alignment for the *more actions\* icon.
 
### [Curves](https://github.com/datagrok-ai/public/tree/master/packages/Curves/CHANGELOG.md)
#### Improvements:
* [#2754](https://github.com/datagrok-ai/public/issues/2754): Implement capability just to connect the points (without fitting).
* [#2103](https://github.com/datagrok-ai/public/issues/2103): Property panel changes.
* [#2105](https://github.com/datagrok-ai/public/issues/2105): Outliers and curves colors changes.
* [#2101](https://github.com/datagrok-ai/public/issues/2101): Improve curves properties and rendering.
* [#1645](https://github.com/datagrok-ai/public/issues/1645): MultiCurveViewer.
#### Fixes:
* [#2797](https://github.com/datagrok-ai/public/issues/2797): Custom fit function doesn't fit correctly.
* [#2748](https://github.com/datagrok-ai/public/issues/2748): Whole table is broken if the cell with curves contains specific data.
 
### [Charts](https://github.com/datagrok-ai/public/tree/master/packages/Charts/CHANGELOG.md)
#### Improvements:
* GROK-15447: Radar viewer: add a part of Viewers tutorial.
* GROK-15409: Radar: improvements (interactivity, color-code, etc.).
* GROK-15340: Radar: Spider Plotting \- Multiple series.
* GROK-14610: Sunburst: animation and selection improvements.
* [#2500](https://github.com/datagrok-ai/public/issues/2500): Sunburst improvements.
#### Fixes:
* GROK-13831: Tree viewer: error when trying to set Size or Color.
* GROK-13832: Tree viewer: doesn't respond on table switching.
* GROK-13390: Radar viewer: harmonization.
* GROK-15245: Radar: error when saving to project.
* GROK-15166: Sunburst: error when trying to switch the Table property.
* GROK-14899: Tree Visualization: Color Agg fails on null values.
 
### DiffStudio
#### Improvements:
* GROK-14487: UI improvement.
* GROK-15713: Fitting.
* GROK-15112: Added Bioreactor to examples.
* GROK-15116: Added E-notation parsing.
* GROK-15294: New sidebar feature.
* GROK-15364: New app view feature.
 
### Fixes:
* GROK-14994: Conditions with JS functions don’t work.
* GROK-14673: UiBuilder: validators strange behavior.
* GROK-15341: Semantic value extractors.
* GROK-13568: Help | Compute: Sensitivity Analysis.
* GROK-15192: Incorrect behavior of inputs in scripts and queries.
* GROK-15376: Column data disappearing while changing time format.
* GROK-15613: 'Ctrl + Wheel' behavior as 'Zoom In/Zoom Out'.
* GROK-14828: Downloaded files doubled after Extract Selected Rows.
* GROK-15834: Table Linking: wrong results when multiple duplicate keys are present in both master and details table.
* GROK-15051: Pie chart: error when selecting Segment angle and hovering the empty value.
* [#2824](https://github.com/datagrok-ai/public/issues/2824): Creating Dataframe from csv loses qnum precision.
 
#### [Box Plot](../../visualize/viewers/box-plot.md):
* GROK-15752: Box Plot: can't select points outside of whiskers.
* GROK-15365: Box plot: error when switching Value axis to the log scale in some cases.
* [#2764](https://github.com/datagrok-ai/public/issues/2764): Wrong area is selected with shift-drag when categories at the left are filtered out.
 
#### [Bar chart](../../visualize/viewers/bar-chart.md):
* GROK-15156: Values are not shown for just negative numbers.
* [#2696](https://github.com/datagrok-ai/public/issues/2696): Add highlighting for bar chart subsets.
 
#### [Filter Panel](../../visualize/viewers/filters.md):
* GROK-15650: Combined boolean filter.
* GROK-15122: Categorical filters: Tooltip and mouse over group is automatically cleared after 200ms on category hover.
* GROK-15206: Filters: categorical filter is broken for hidden ~name columns.
* [#2647](https://github.com/datagrok-ai/public/issues/2647): File refresh breaks structure filter.
* [#2684](https://github.com/datagrok-ai/public/issues/2684): Wrong filters are added if there are two views in some cases.
* [#2185](https://github.com/datagrok-ai/public/issues/2185): Fuzzy filter enhancements.
 
#### [Tutorials](https://github.com/datagrok-ai/public/tree/master/packages/Tutorials/CHANGELOG.md):
* GROK-15143: Multivariate Analysis and R-Groups Analysis tutorial: Top menu items don't open properly.
* GROK-15132: Timelines: error when using the legend.




## 2024-06-05 1.18.9

### Addressed Issues:

* GROK-15834: Table Linking: wrong results when multiple duplicate keys are present in both master and details table 
* [#2839](https://github.com/datagrok-ai/public/issues/2839): Show Help in "Settings \- Windows" remains checked always 


## 2024-05-19 1.18.8

### Improvements:
* GROK-15681: Projects Management 
* [#2831](https://github.com/datagrok-ai/public/issues/2831): Default open folders in the Browse Menu 

### Fixes:

* [#2832](https://github.com/datagrok-ai/public/issues/2832): Tab Name changes when reopening Project 
* [#2840](https://github.com/datagrok-ai/public/issues/2840): "Global Search" button on bottom-right bar menu improvements
* [#2834](https://github.com/datagrok-ai/public/issues/2834): Loading of DB tables appears to continue forever and never complete 


## 2024-05-15 1.18.7

### Improvements:

* [#2845](https://github.com/datagrok-ai/public/issues/2845): Changing the colors of labels in the properties

### Fixes:

* [#2824](https://github.com/datagrok-ai/public/issues/2824): Creating Dataframe from csv loses qnum precision


## 2024-05-08 1.18.6

### Fixes:

* FormViewer: fixed an issue with dragging column headers
* Scripting: Clear temporary token for input files


## 2024-05-01 1.18.5

### Improvements:

* Error handling: "Report an error" dialog when you click on exclamation point.
* Sticky meta: improve visibility for users:
   * Show preview for entity types
   * Support multi-values edit 
   * Entity type settings with TabControl 
   * Better type selector for schema editor.
* JS API for entity properties/types/schemas.
* [#2733](https://github.com/datagrok-ai/public/issues/2733): Viewers | Legend: provide diverse palette and markers' shape.
* [#2742](https://github.com/datagrok-ai/public/issues/2742): Viewers: Indicate the selected column in the column selector.
* [#2772](https://github.com/datagrok-ai/public/issues/2772): Pivot table: add a property to hide sparklines.
* GROK-14412: Trellis plot | Column selectors: ESC should close the popup.

### Fixes:

* Azure blob: directory shown as a file.
* GROK-15190: Scatterplot, line chart: updating formula lines does not work for those specified via a dataframe.
* GROK-15249: Error when pressing Column Properties on a programmatically created Grid.
* [#2642](https://github.com/datagrok-ai/public/issues/2642): Filtering done using viewers is unexpectedly reset on applying other filters in some cases.
* [#2644](https://github.com/datagrok-ai/public/issues/2644): Forms | Next row button: incorrect behavior after sorting the grid.
* [#2749](https://github.com/datagrok-ai/public/issues/2749): Line chart and Trellis plot: '+' icon is missing.
* [#1092](https://github.com/datagrok-ai/public/issues/1092): Multiple issues with Conditional Color Dialog.
* [#2764](https://github.com/datagrok-ai/public/issues/2764): Box plot: wrong area is selected with shift-drag when categories at the left are filtered out.
* [#2509](https://github.com/datagrok-ai/public/issues/2509): Pinned columns in multiple views for the same table cause performance issues in some cases.
* [#2782](https://github.com/datagrok-ai/public/issues/2782): Line chart: incorrect color-coding while setting filters.
* [#2636](https://github.com/datagrok-ai/public/issues/2636): Trellis plot | Bar chart: contrast color for blue instead of black works when there is no Stack set to the bar chart inner viewer.
* [#2793](https://github.com/datagrok-ai/public/issues/2793): Scatter plot with 'filter by zoom': filtering by scatter plot is shown in filter panel's '?' after resetting in some cases.


## 2024-04-05 1.18.4

### Addressed Issues

* (Bug) GROK-15334: Large file open fails (WIP)


## 2024-03-28 1.18.3

### Addressed Issues

* (Bug) GROK-15300: Table duplicates after changing ID 


## 2024-03-19 1.18.2

### Addressed Issues

* (Bug) GROK-15215: Async package functions don't work in Add new column 
* (Bug) GROK-15216: Views migration didn't run 
* (Bug) GROK-15166: Sunburst: error when trying to switch the Table property 


## 2024-03-14 1.18.1

### Addressed Issues

* (Improvement) GROK-14698: Color coding: improve traffic lights schemas 


## 2024-03-08 1.18.0

Datagrok 1.18 release focuses on usability improvements and stability:

* Improved **Browser** for navigation, preview, and convenient interaction with all available functionality on the Datagrok platform.
* Ability to clone projects.
* New CVM queue employs multiple JKG nodes to ensure scalable and streamlined script execution.
* New MultiChoice cell renderer to show a number of checkboxes within a cell. For details, see [Add multiple check box and dropdown list to cell of table](https://community.datagrok.ai/t/add-multiple-check-box-and-dropdown-list-to-cell-of-table/835/9?u=oahadzhanian.datagrok.ai).
* New Tags cell renderer for comma-separated rows.
* Automatic error handling report: A "Report an Error" dialog box appears when you click the exclamation point icon (Beta version).
* Sticky meta (Beta version).

### Visualization and usability improvements

* Calculated columns:
  * Save calculated columns to the projects with data synchronization.
  * Update values if the base column values changed.
*  [#2634](https://github.com/datagrok-ai/public/issues/2634): Synchronization between Bar chart and Line chart.
* [#2455](https://github.com/datagrok-ai/public/issues/2455): Viewers: Add the ability to choose several categories with ctrl click.
* [#2548](https://github.com/datagrok-ai/public/issues/2548): Viewers | Legends: implement switching between categories using up and down arrows when one category is selected exclusively.
* [#2635](https://github.com/datagrok-ai/public/issues/2635): Customization of actions on the Context Panel.
* [#2664](https://github.com/datagrok-ai/public/issues/2664): Implement a cascade column renaming used by Formula lines.
* [#2467](https://github.com/datagrok-ai/public/issues/2467): Allow to set tab label for stacked viewers programmatically.
* [#2598](https://github.com/datagrok-ai/public/issues/2598): Context menus to hide selectors.
* [#2564](https://github.com/datagrok-ai/public/issues/2564): Column selectors should not overlap the viewer.
* [#2299](https://github.com/datagrok-ai/public/issues/2299): Make viewers aware of surroundings.
* Rename viewer by hotkey **F2**.
* Column selectors: ESC should close the popup.
* Fixed:
   * [#2570](https://github.com/datagrok-ai/public/issues/2570): Viewers: legend does not show all categories after filtering in specific cases.
   * [#2590](https://github.com/datagrok-ai/public/issues/2590): Calculated columns dialog: not all functions are available in list, search is not working as expected.
* [Grid](../../visualize/viewers/grid.md):
   * [#2645](https://github.com/datagrok-ai/public/issues/2645): Support single decimal place for rounding in calculated columns
   * Rendering table cells with any viewer.
   * Custom numeric format.
   * Fixed:
     * [#2468](https://github.com/datagrok-ai/public/issues/2468): Cell values starting with the quote character are parsed incorrectly.
     * [#2483](https://github.com/datagrok-ai/public/issues/2483): Order and hide columns and Tooltip dialogs: cosmetic issues.
     * [#2503](https://github.com/datagrok-ai/public/issues/2503): Order and Hide: Tooltip with null category shouldn’t be that big.
     * [#2523](https://github.com/datagrok-ai/public/issues/2523): Columns added after performing structure search cannot be shown again if hidden.
* [Filter Panel](../../visualize/viewers/filters.md):
   * [#2185](https://github.com/datagrok-ai/public/issues/2185): Fuzzy filter enhancements:
   * The filter header corresponds to the column name.
   * Ability to combine with free-text filter
   * Save to layout typical searches / search patterns
   * [#2502](https://github.com/datagrok-ai/public/issues/2502): buttons on the filter panel appear on hovering only.
   * Fixed:
     * [#2504](https://github.com/datagrok-ai/public/issues/2504): The question mark hovering doesn’t show the filtering of the scaffold tree viewer.
     * [#2647](https://github.com/datagrok-ai/public/issues/2647): #2647: File refresh breaks structure filter
     * [#2684](https://github.com/datagrok-ai/public/issues/2684): #2684: Wrong filters are added if there are two views in some cases.
* [Scatterplot](../../visualize/viewers/scatter-plot.md):
   * [#2456](https://github.com/datagrok-ai/public/issues/2456): ignore negatives and zero values when switching to log scale.
   * [#2480](https://github.com/datagrok-ai/public/issues/2480):default range should take into account only the rows where both axes values are present.
   * Fixed:
     * [#2330](https://github.com/datagrok-ai/public/issues/2330): Scatterplot has to show empty categories like Box plot and Bar chart do.
     * [#2487](https://github.com/datagrok-ai/public/issues/2487): datetime columns are not supported when adding formula lines.
     * [#2488](https://github.com/datagrok-ai/public/issues/2488): Color selector is missing for formula lines dialog.
     * [#2489](https://github.com/datagrok-ai/public/issues/2489): Color selector opened from the 'style' section in properties panel cannot be closed.
     * [#2505](https://github.com/datagrok-ai/public/issues/2505): If categories are filtered out they shouldn't be shown in the marker legend.
     * [#2527](https://github.com/datagrok-ai/public/issues/2527): error occurs when axes are set to log scale in some cases.
     * [#2530](https://github.com/datagrok-ai/public/issues/2530): Errors on hovering scatterplot with formula lines.
     * [#2665](https://github.com/datagrok-ai/public/issues/2665): Ctrl+Shift+Drag doesn't work for deselection.
     * [#2671](https://github.com/datagrok-ai/public/issues/2671): filtered out data is unexpectedly shown in some cases on filtering the plot further using markers legend.
     * [#671](https://github.com/datagrok-ai/public/issues/671): 'Formula lines' dialog preview is not updated if scatterplot configuration is changed
* [Trellis plot](../../visualize/viewers/trellis-plot.md):
   * [#2457](https://github.com/datagrok-ai/public/issues/2457): redesigning:
     * Move X and Y selectors to the corresponding axes.
     * Most relevant options of the inner viewer should be displayed at the top of the viewer
     * Selectors instead of combo-boxes
   * [#2336](https://github.com/datagrok-ai/public/issues/2336): improve the display of outer and inner viewers properties.
   * [#2636](https://github.com/datagrok-ai/public/issues/2636): Trellis plot | Bar chart: improvements.Now they work with global axes, support multiple splits, and display their axis with category names and colors on the left.
   * Better auto-selection of columns for pie chart and bar chart.
   * Integrate sparklines (Beta version). For details, see [Trellis: sparklines](https://community.datagrok.ai/t/trellis-sparklines/818/1)
   
* [Histogram](../../visualize/viewers/histogram.md):
   * [#2549](https://github.com/datagrok-ai/public/issues/2549): add highlighting for split lines as for line chart.
   * [#2685](https://github.com/datagrok-ai/public/issues/2685): add stacked view option when split is set.
   * Fixed:
     * [#2472](https://github.com/datagrok-ai/public/issues/2472): selected column is not synchronised between in-plot controls and properties panel.
     * [#2524](https://github.com/datagrok-ai/public/issues/2524): Split: no tooltip and Y axis when the values are not normalized.
* [Line chart](../../visualize/viewers/line-chart.md):
   * [#2357](https://github.com/datagrok-ai/public/issues/2357): custom tooltip.
   * [#2623](https://github.com/datagrok-ai/public/issues/2623): implement the one-click way to set the Split by columns.
   *  Interactivity on a stacked bar chart.
   * Fixed:
     * [#2485](https://github.com/datagrok-ai/public/issues/2485): too much empty space if Y axis is logarithmic for specific data.
     * [#2577](https://github.com/datagrok-ai/public/issues/2577): colors in the legend do not match color in the plot in some cases.
     * [#2608](https://github.com/datagrok-ai/public/issues/2608): Line chart with row source = selected: all rows are deselected on clicking line in some cases.
* [Box plot](../../visualize/viewers/box-plot.md): [#2484](https://github.com/datagrok-ai/public/issues/2484): no indication what dot is hovered, wrong dot is selected in some cases.
* [Bar chart](../../visualize/viewers/bar-chart.md):
   * [#2584](https://github.com/datagrok-ai/public/issues/2584): add new aggregations for time ranges such as weeks and days.
   * [#2298](https://github.com/datagrok-ai/public/issues/2298): squeeze the white space around the title.
  * Fixed: [#2562](https://github.com/datagrok-ai/public/issues/2562): error when date column is used as category and split function is non-default.
* For pivot table fixed:
    * [#2535](https://github.com/datagrok-ai/public/issues/2535): non-default aggregations are not saved in layout.
    * [#2606](https://github.com/datagrok-ai/public/issues/2606): colours, formatting, shown/hidden columns are reset on filtering.
*  [Tile viewer](../../visualize/viewers/tile-viewer.md): ability to specify custom lane categories.

*  [Correlation plot](../../visualize/viewers/tile-viewer.md): [#2653](https://github.com/datagrok-ai/public/issues/2653): tooltip issue.

*  [Forms viewer](../../visualize/viewers/forms.md): [#2652](https://github.com/datagrok-ai/public/issues/2652): implemented an option to download a viewer as PNG.

### Data Access

*  Azure Blob integration. For details, see [Azure Blob](../../access/files/shares/azure.md)

### Improvements for developers

* Core: Inputs: validation and harmonization 
   * New `UserInput`.
   * Two-way binding for `DG.InputForm`.
   * `bindFuncCall` method for `InputForm`.
   * Computable "visible" and "enabled" properties
   * `TableInput`: Implement capability to get files from Shares.
   * `NumberInput`: Add capability to set custom step.
   * Removed default rounding for float inputs.
*  Table, RowSource, and Filter properties in js viewers.
*  Scripts logging system.
* Script saving via Ctrl+S.
* Scripting | Variables: added the scrollbar to the variables window.
*  API to send an email from DG.

#### JS API

  * ` Column.getNumber()` method.
  * Dynamical dependencies load for `Utils.loadJsCss(files)`.
   * Auto-generation of Dart interop interface.
  * Exposed `FuncCall` fields to JS.
  * Setter for `InputBase` caption.
  * Made `EventData` for `DataFrame` typed.

## 2024-02-08 1.17.13

### Addressed Issues

* (Bug) [#2647](https://github.com/datagrok-ai/public/issues/2647): #2647: File refresh breaks structure filter
* (Bug) [#2639](https://github.com/datagrok-ai/public/issues/2639): Structure in filter is unexpectedly changed when there are two views in some cases
* (Bug) [#2642](https://github.com/datagrok-ai/public/issues/2642): Filtering done using viewers is unexpectedly reset on applying other filters in some cases (WIP)


## 2024-01-26 1.17.12

### Addressed Issues

* (Bug) Windows share mounts with wrong permissions


## 2024-01-22 1.17.11

### Addressed Issues

* (Bug) [#2628](https://github.com/datagrok-ai/public/issues/2628): Structure filter is not applied in some cases when there are two views opened 


## 2023-12-21 1.17.10

### Addressed Issues

* (Bug) [#2535](https://github.com/datagrok-ai/public/issues/2535): Pivot table: non-default aggregations are not saved in layout 
* (Bug) GROK-14381: Grid | Add Summary column: impossible to hide or add a column through the Order or Hide Columns dialog 
* (Bug) [#2577](https://github.com/datagrok-ai/public/issues/2577): Line chart: colours in the legend do not match colour in the plot in some cases 
* (Bug) [#2570](https://github.com/datagrok-ai/public/issues/2570): Viewers: legend does not show all categories after filtering in specific cases (WIP)
* (Bug) [#2574](https://github.com/datagrok-ai/public/issues/2574): Multiple errors in console on hovering line chart with logarithmic axis and no data 


## 2023-12-11 1.17.9

### Addressed Issues

* GROK-14356: Projects: data-sync option is not shown for .sdf files 


## 2023-12-07 1.17.8

### Addressed Issues

* (Improvement) [#2455](https://github.com/datagrok-ai/public/issues/2455): Viewers: Add the ability to choose several categories with ctrl click (WIP)


## 2023-12-06 1.17.7

### Addressed Issues

* (Improvement) [#2455](https://github.com/datagrok-ai/public/issues/2455): Viewers: Add the ability to choose several categories with ctrl click (WIP)
* (Bug) [#2562](https://github.com/datagrok-ai/public/issues/2562): Bar chart: error when date column is used as category and split function is non-default 


## 2023-11-29 1.17.6

### Addressed Issues

* (Improvement) [#2455](https://github.com/datagrok-ai/public/issues/2455): Viewers: Add the ability to choose several categories with ctrl click 
* (Improvement) [#2527](https://github.com/datagrok-ai/public/issues/2527): Scatterplot: error occurs when axes are set to log scale in some cases
* (Improvement) [#2330](https://github.com/datagrok-ai/public/issues/2330): Scatter plot: has to show empty categories like Box plot and Bar chart do


## 2023-11-28 1.17.5

### Addressed Issues



## 2023-11-23 1.17.4

### Addressed Issues

* (Improvement) GROK-14235: JS API: Utils.loadJsCss(files): dynamically load dependencies (WIP)
* (Bug) GROK-14261: HomeView from PowerPack doesn't work as expected (WIP)
* (Bug) GROK-14229: A user cannot add anything to Favorites 
* (Bug) GROK-14246: JS API: `items` setter on ChoiceInput creates empty option 
* (Bug) GROK-14243: Viewers: ScatterPlot ignores initial filtration (WIP)
* (Bug) GROK-14254: JS API: ui.input.forProperty fails on DF property 
* (Bug) GROK-14272: JS API: ui.input.forProperty has no option to set float format 
* (Bug) GROK-14277: Filter Panel:  Unsupported operation: NaN.ceil() error in some cases 
* (Bug) [#2523](https://github.com/datagrok-ai/public/issues/2523): Columns added after performing structure search cannot be shown again if hidden 
* (Improvement) [#2357](https://github.com/datagrok-ai/public/issues/2357): Line chart: Custom tooltip 
* (Improvement) [#1063](https://github.com/datagrok-ai/public/issues/1063): A way to indicate current viewer 
* (Bug) [#2530](https://github.com/datagrok-ai/public/issues/2530): Errors on hovering scatter plot with formula lines 
* (Improvement) [#2455](https://github.com/datagrok-ai/public/issues/2455): Viewers: Add the ability to choose several categories with ctrl click 
* (Improvement) [#2480](https://github.com/datagrok-ai/public/issues/2480): Scatterplot: default range should take into account only the rows where both axes values are present 
* (Bug) GROK-14271: JS API: defaultValue is ignored by `ui.input.forProperty` 
* GROK-7105: Forms viewer (WIP)
* (Bug) [#2527](https://github.com/datagrok-ai/public/issues/2527): Scatterplot: error occurs when axes are set to log scale in some cases

## 2023-11-13 1.17.3

* (Bug) GROK-14182: ui.choiceInput breaks the old code 
* (Bug) GROK-14210: UI: Format is ignored on floatInput 
* (Improvement) [#2298](https://github.com/datagrok-ai/public/issues/2298): Bar Chart: squeeze the white space around the title (WIP)
* (Improvement) [#2456](https://github.com/datagrok-ai/public/issues/2456): Scatterplot: ignore negatives and zero values when switching to log scale 
* (Improvement) [#2480](https://github.com/datagrok-ai/public/issues/2480): Scatterplot: default range should take into account only the rows where both axes values are present 
* (Bug) [#2485](https://github.com/datagrok-ai/public/issues/2485): Line chart: too much empty space if Y axis is logarithmic for specific data 
* (Improvement) [#2484](https://github.com/datagrok-ai/public/issues/2484): Box plot: no indication what dot is hovered, wrong dot is selected in some cases 
* (Bug) [#2483](https://github.com/datagrok-ai/public/issues/2483): Order and hide columns and Tooltip dialogs: cosmetic issues 
* (Improvement) [#2467](https://github.com/datagrok-ai/public/issues/2467): Allow to set tab label for stacked viewers programmatically 
* GROK-14215: Tiles Viewer: Ability to specify custom lane categories 
* (Bug) [#2489](https://github.com/datagrok-ai/public/issues/2489): Colour selector opened from the 'style' section in properties panel cannot be closed 
* (Bug) [#2503](https://github.com/datagrok-ai/public/issues/2503): Order and Hide: Tooltip with null category shouldn't be that big 
* (Bug) GROK-14221: ChoiceInputs are not refreshing in FRAC classification query 
* (Improvement) GROK-14235: JS API: Utils.loadJsCss(files): dynamically load dependencies (WIP)


## 2023-11-03 1.17.2



## 2023-11-02 1.17.1

* (Improvement) GROK-14180: Implement waiting mode in db.runLocked 
* (Bug) [#2330](https://github.com/datagrok-ai/public/issues/2330): Scatter plot: has to show empty categories like Box plot and Bar chart do 
* (Bug) [#2468](https://github.com/datagrok-ai/public/issues/2468): Cell values starting with the quote character are parsed incorrectly 


## 2023-10-31 1.17.0

Datagrok 1.17 release focuses on stability, performance, and usability improvements:

* The ability to configure the platform through the **Settings** wizard.
* Browser designed for navigation, preview, and convenient access to everything available on the platform: features, applications, plugins, models, shared dashboards, and more.
* Function view now shows function signature if parameters are not user-editable.
* Improved client-side caching of function and query results. To learn more, see [Caching function results](https://datagrok.ai/help/develop/function_results_cache#client-side-cache)
* Summary viewer that aggregates the numeric attributes of features.
* Capability to render table cells with any viewer, along with support for linked tables in in-grid dataframes.
* **Content** tab on the **Context Panel**, making it easy to compare selected rows, filters, highlights, and more.
* **Plots** context pane automatically visualizing selected columns, see [UX updates](https://community.datagrok.ai/t/ux-updates/544/5?u=oahadzhanian.datagrok.ai) for details.
* Lookup tables that let you initialize inputs with a set of pre-defined values. For details, see [Lookup tables](https://datagrok.ai/help/datagrok/concepts/functions/func-params-annotation#lookup-tables).
* [Okta](https://www.okta.com/) support for authorization.

### Visualization and usability improvements

* [#2335](https://github.com/datagrok-ai/public/issues/2335): Tooltip harmonization.
* [#1063](https://github.com/datagrok-ai/public/issues/1063): A way to indicate current viewer.
* Auto sizing for viewers' legend.
* New color input for color coding.
* Columns pane: the ability to filter columns by semantic type and preview their content in the grid on mouse-over.
* Viewers: handling of doubles, `infinity`, `-infinity`, `NaN`, and large-scale floating point numbers.
* Slider in number inputs.
* Fixed:
   * [#2368](https://github.com/datagrok-ai/public/issues/2368): Sketcher | structure lookup service: fix the behavior.
   * [#1960](https://github.com/datagrok-ai/public/issues/1960): "Link Tables" does not persist when table data refreshed/reloaded.
   * [#2011](https://github.com/datagrok-ai/public/issues/2011): Viewers: Inconsistencies in saved layouts won't allow creating a viewer when such layout is applied.
   * [#2060](https://github.com/datagrok-ai/public/issues/2060): Projects: the tabs order is changed from the original uploaded project state.
   * [#2272](https://github.com/datagrok-ai/public/issues/2272): Loading page performance issues.
   * [#2284](https://github.com/datagrok-ai/public/issues/2284): Certain scenarios make it impossible to move or stack a viewer.
   * [#2297](https://github.com/datagrok-ai/public/issues/2297): Viewers | Tooltip: a custom tooltip should not be limited in the number of values.
   * [#2315](https://github.com/datagrok-ai/public/issues/2315): Unexpected columns appear in 'Order and hide' columns dialog if the scaffold tree viewer or structure filter is added.
   * [#2344](https://github.com/datagrok-ai/public/issues/2344): Categorical filter hovering does not highlight data in plots.
   * [#2328](https://github.com/datagrok-ai/public/issues/2328): Viewer docked at the left or at the top is sometimes unexpectedly hidden.
* [Grid](../../visualize/viewers/grid.md):
  * [#2271](https://github.com/datagrok-ai/public/issues/2271): Introduced a property to disable "read-only" warning.
  * [#2367](https://github.com/datagrok-ai/public/issues/2367): Pinned rows improvements.
  *  [#2348](https://github.com/datagrok-ai/public/issues/2348): Add the Sync New Columns setting to the grid.
  * Custom column order is now saved after refresh.
  * The ability to clone grid columns not affecting the dataframe.
  * Column selection by holding Shift and mouse-dragging the column header, see [Grid updates](https://community.datagrok.ai/t/grid-updates/616/11).
  * Fixed:
    * [#2293](https://github.com/datagrok-ai/public/issues/2293): Cannot search in 'order and hide columns' if table contains a summary column
    * Sorting isn't saved after the project reopens.
    * Custom cell sizes get reset when you open a project.
    * Order and Hide Columns: the Column List's hamburger menu appears in the Search field.
* [Filter Panel](../../visualize/viewers/filters.md):
  * [#2358](https://github.com/datagrok-ai/public/issues/2358): Layouts: save the persistence of filter item search visibility.
  * [#2165](https://github.com/datagrok-ai/public/issues/2165): Always add new filters to the top.
  * Fixed:
    * [#2316](https://github.com/datagrok-ai/public/issues/2316): Errors on combining scaffold tree with some other filters.
    * [#2359](https://github.com/datagrok-ai/public/issues/2359): For some filters the item search is not working when another filter search is open
    * The min and max values for dates do not fit the entire value.
* For [Scatter plot](../../visualize/viewers/scatter-plot.md) fixed:
  * [#2285](https://github.com/datagrok-ai/public/issues/2285): Default marker size is not changed on moving slider.
  * [#2408](https://github.com/datagrok-ai/public/issues/2408): Incorrect tooltip for 'zoom and filter'.
  * [#2402](https://github.com/datagrok-ai/public/issues/2402): Min/max values are not saved/restored in the layout.
  * [#2330](https://github.com/datagrok-ai/public/issues/2330): Has to show empty categories like Box plot and Bar chart do.
* [Trellis plot](../../visualize/viewers/trellis-plot.md):

  * [#2336](https://github.com/datagrok-ai/public/issues/2336): Improved the display of outer and inner viewers properties.
  * [#2300](https://github.com/datagrok-ai/public/issues/2300): Added the ability to visualize just one category.
       
* [Line Chart](../../visualize/viewers/line-chart.md):
  * [#2302](https://github.com/datagrok-ai/public/issues/2302): Allows to split by multiple columns.
  * [#2357](https://github.com/datagrok-ai/public/issues/2357): Custom tooltip.
  * [#2360](https://github.com/datagrok-ai/public/issues/2360): Hovering and selecting improvements:
    * Hovering over a line highlights data on other viewers.
    * Ctrl+click, Shift +click, and Ctrl+Shift +click selection toggling actions for lines
    * Clicking the line selects the values on other viewers.
  * [#2049](https://github.com/datagrok-ai/public/issues/2049): Added min/max properties.
  * [#2396](https://github.com/datagrok-ai/public/issues/2396): Line chart with split performance issue on large amounts of data.
  * Fixed:
    * [#2163](https://github.com/datagrok-ai/public/issues/2163): Line chart lines are displayed incorrectly in some cases (drop below the X axis).
    * [#2395](https://github.com/datagrok-ai/public/issues/2395): Unaggregated tooltip and hovered row do not match in some cases.
    * [#2388](https://github.com/datagrok-ai/public/issues/2388): Optimize the display of labels on the x-axis.
    * [#2350](https://github.com/datagrok-ai/public/issues/2350): 'Error loading line chart' on applying layout if line chart has a categorical X axis and **Row Source** is set to `Selected`.

* [Histogram](../../visualize/viewers/histogram.md):

  * [#2197](https://github.com/datagrok-ai/public/issues/2197): Add zoom slider to X axis of Histogram.
  * Fixed:
    * [#2301](https://github.com/datagrok-ai/public/issues/2301): Histogram slider's end can be dragged out of valid range.
    * [#2364](https://github.com/datagrok-ai/public/issues/2364): Some rows are still shown if filter's range is set out of data's boundaries.
    * [#2329](https://github.com/datagrok-ai/public/issues/2329): Histogram 'normalise to filter' option: Y axis values are not updated on zooming.
* [Box Plot](../../visualize/viewers/box-plot.md):
  * [#2203](https://github.com/datagrok-ai/public/issues/2203): Improvements for structures rendering.
  * [#2057](https://github.com/datagrok-ai/public/issues/2057): Added structures rendering on the X axis.
  * FIxed:
    * [#2270](https://github.com/datagrok-ai/public/issues/2270): Selected row dots are not brought to the foreground.
    * [#2286](https://github.com/datagrok-ai/public/issues/2286): Markers overlap with X axis and categories in some cases.
* For [bar chart](../../visualize/viewers/bar-chart.md):  [#2298](https://github.com/datagrok-ai/public/issues/2298): Squeezed the white space around the title.
* Pivot table:
   * [#2166](https://github.com/datagrok-ai/public/issues/2166): Improvements: 
     * An option for filtering or reflecting a filter.
     * The **Row Source** setting.
     * Provided a way to collapse the configuration via the **Show Header** option.
  * Fixed [#2198](https://github.com/datagrok-ai/public/issues/2198): Editing and saving issues.

### Data Access

* Improvements:
  * Schema browsing support.
  * Warning balloon if there is no connection with GrokConnect.
* Fixed [#2454](https://github.com/datagrok-ai/public/issues/2454): Categories contain duplicate values using datadriven queries.

### Improvements for developers

* [#2273](https://github.com/datagrok-ai/public/issues/2273): Implemented global permissions for some entities.
* Improved exporting of logs to CloudWatch: send logs to Amazon CloudWatch according to the Log Export Blocks you can create.
* Support of other than CSV formats to TableInput.
* The scrollbar for the variables window.
* Python scripts harmonization:
  * Fixed DateTime parameter.
  * Fixed column_list parameter.
* Julia scripts harmonization: fixed boolean parameter support.
* Ability to set choice items and values separately.
* Option to allow `DG.Accordion` panes dragging.

#### JS API

* [#2118](https://github.com/datagrok-ai/public/issues/2118): API configuration to define layout import settings 
 * Capability to set grid cell value with `onCellValueEdited` event.
* Capability to get stats from values array.

### Enhancements in packages

* [Bio](https://github.com/datagrok-ai/public/blob/master/packages/Bio/CHANGELOG.md): ToAtomicLevel is available for non-linear HELM structures, enhancements for VdRegionsViewer and WebLogo, other improvements and bug fixes.
* [Chem](https://github.com/datagrok-ai/public/blob/master/packages/Chem/CHANGELOG.md): Highlighting multiple substructures with different colors inside one molecule structure; filtering by superstructure, exact structure, and similarity score; scaffold tree integration with color-coded fragments and other enhancements.
* [Dendrogram](https://github.com/datagrok-ai/public/tree/master/packages/Dendrogram/CHANGELOG.md): hierarchical clustering for molecules, the ability to select all leaves from a certain node, a separate loader view to the dendrogram, the ability to switch distance calculation method for macromolecules (Hamming or Levenstein) in case of MSA, and other enhancements.
* [PowerGrid](https://github.com/datagrok-ai/public/tree/master/packages/PowerGrid/CHANGELOG.md): new [smart forms](https://community.datagrok.ai/t/powergrid-smartform/774/1) - a special cell type that renders values from multiple columns in one cell; image rendering on double click on the image URL.
* [Curves](https://github.com/datagrok-ai/public/tree/master/packages/Curves/CHANGELOG.md): aggregations for series statistics, **Context Pane** changes, improved curves demo application, and bug fixes.
* [Peptides](https://github.com/datagrok-ai/public/blob/master/packages/Peptides/CHANGELOG.md): Select All and Deselect All functionality to all viewers, mean activity column for the Most Potent Residues viewer,  WebLogo to Selection table.
* [Charts](https://github.com/datagrok-ai/public/tree/master/packages/Charts/CHANGELOG.md): improvements for sunburst and surface plot.
* [Tutorials](https://github.com/datagrok-ai/public/tree/master/packages/Tutorials/CHANGELOG.md): new chem tutorials on activity cliffs, R-group analysis, similarity and diversity search, substructure search and filtering.

## 2023-10-04 1.16.7

* (Improvement): Connectors: Added gzip compression of socket message
* (Improvement): Connectors: Added offset, fixed caching 
* (Improvement): DataQuery: Added setting of gzip and batchSize through SQL script meta

## 2023-09-14 1.16.6

* (Bug) GROK-13904: Redocking viewer leads to multiple entries in .viewers list 
* (Bug) [#2350](https://github.com/datagrok-ai/public/issues/2350): 'Error loading line chart' on applying layout, if line chart has categorical X axis and row source = Selected 


## 2023-09-06 1.16.5

* (Bug) GROK-13874: Box plot: incorrectly displays values on the X axis 
* (Bug) [#2329](https://github.com/datagrok-ai/public/issues/2329): Histogram 'normalize to filter' option: Y axis values are not updated on zooming 
* (Bug) [#2328](https://github.com/datagrok-ai/public/issues/2328): Viewer docked at the left or at the top is sometimes unexpectedly hidden 
* (Bug) [#2315](https://github.com/datagrok-ai/public/issues/2315): Unexpected columns appear in 'Order and hide' columns dialog if scaffold tree viewer or structure filter is added 

## 2023-08-30 1.16.4

* (Bug) [#2315](https://github.com/datagrok-ai/public/issues/2315): Unexpected columns appear in 'Order and hide' columns dialog if scaffold tree viewer or structure filter is added 
* (Improvement) [#2057](https://github.com/datagrok-ai/public/issues/2057): Box plot has to render structures on the X axis 
* (Bug) [#2293](https://github.com/datagrok-ai/public/issues/2293): Cannot search in 'order and hide columns' if table contains a summary column 
* (Bug) [#2163](https://github.com/datagrok-ai/public/issues/2163): Line chart lines are displayed incorrectly in some cases (drop below the X axis) 
* (Bug) [#2285](https://github.com/datagrok-ai/public/issues/2285): Scatter plot: default marker size is not changed on moving slider 
* (Bug) [#2316](https://github.com/datagrok-ai/public/issues/2316):  Errors on combining scaffold tree with some other filters 
* (Bug) [#2272](https://github.com/datagrok-ai/public/issues/2272): Loading page performance issues (WIP)
* (Bug) [#2301](https://github.com/datagrok-ai/public/issues/2301): Histogram slider's end can be dragged out of valid range 
* (Bug) GROK-13843: ScatterPlot: displaying the legend issues (WIP)
* (Improvement) [#1063](https://github.com/datagrok-ai/public/issues/1063): A way to indicate current viewer 


## 2023-08-17 1.16.3

* (Improvement) GROK-13737: Cleanup package if it's not in NPM and is not installed 

## 2023-08-15 1.16.2

* Improvement:
* (Improvement) [#2166](https://github.com/datagrok-ai/public/issues/2166): Pivot table: improvements 
* (Improvement) [#2165](https://github.com/datagrok-ai/public/issues/2165): Filter Panel: Always add new filters to the top  (WIP)
* (Bug) GROK-13526: Interpret URLs as hyperlinks: incorrect behavior with line wrapping  (WIP)
* (Bug) GROK-13632: when starting an unsaved query with parameters in dataQueryView, the form with editing parameters is not shown on the first run 
* (Bug) [#2060](https://github.com/datagrok-ai/public/issues/2060): Projects:  the tabs order is changed from the original uploaded project state 
* (Bug) GROK-13718: Grid: custom renderers: custom cell sizes get reset when you open a project 
* (Improvement) [#2197](https://github.com/datagrok-ai/public/issues/2197): Add zoom slider to X axis of Histogram 
* (Bug) GROK-13719: After applying Refresh several times, there is an error:  `Resulting message size [x] is too large for configured max of [y]` 
* (Bug) GROK-13448: "Refresh"-ing DB query resets the sizes of Grid cells 
* (Bug) GROK-13717: Unable to use scripts in add new column 
* (Bug) GROK-13750: Can't create custom viewer from JS 
* (Bug) GROK-13651: Databases | Connections: testing connection issue 


## 2023-07-27 1.16.1

* Improvements:
  * Columns pane: column preview in the grid on mouse-over 
  * Columns pane: ability to filter by semantic type 
* Fixes:
  * GROK-13559: Credentials manager for WebServices is not always working 
  * GROK-13609: Incorrect value for scalar parameters in DB 
  * GROK-13617: FuncCall not serialized properly 


## 2023-07-20 1.16.0

Datagrok 1.16 release focuses on performance and usability improvements:

* Improved data streaming, leading to major performance improvements for small queries.
* Client-side caching of function results, allowing to save time and server capacity when executing immutable functions and queries.
* Improved logging system for queries, featuring a debug mode that helps understand time allocation and optimize performance. Additionally, we provide new tools for tuning specific queries.
* Simplified and updated Datagrok local installation, now available with a [one-click script](https://datagrok.ai/help/develop/admin/docker-compose).
* Self-adjustable viewer layouts that keep viewers usable even in a small window.
* [Fuzzy text filter](../../visualize/viewers/filters.md#text-filter), which enables users to create keyword categories and search or filter text using them.
* URLs are now automatically interpreted as hyperlinks in text. For details, see [visualization-related updates](https://community.datagrok.ai/t/visualization-related-updates/521/36).


### Visualization and usability improvements

* [#1972](https://github.com/datagrok-ai/public/issues/1972): Implemented [fuzzy text filter](../../visualize/viewers/filters.md#text-filter).
* [#2056](https://github.com/datagrok-ai/public/issues/2056): Added [URLs rendering](https://community.datagrok.ai/t/visualization-related-updates/521/36) as hyperlinks in text.
* Enhanced [multiline text rendering](https://community.datagrok.ai/t/visualization-related-updates/521/34)
* Enhanced composed cell renderer layout:
   * Subclasses can now dynamically control the rendering of external content. Labels are evenly distributed over the available space when external content is not rendered.
   *  Added file creation and modification timestamps to file cards for extra information.
* Fixed:
  * [#2011](https://github.com/datagrok-ai/public/issues/2011): Viewers: Inconsistencies in saved layouts won't allow creating a viewer when such layout is applied.
  * [#1929](https://github.com/datagrok-ai/public/issues/1929): Viewers do not support hex colors set in categorical coloring.
  * [#1838](https://github.com/datagrok-ai/public/issues/1838): Error on adding a custom viewer.
  * [#1960](https://github.com/datagrok-ai/public/issues/1960): **Link Tables** does not persist when table data is refreshed/reloaded.
  * [#1964](https://github.com/datagrok-ai/public/issues/1964): Adding group membership via user-based interface does not save.
  * [#1983](https://github.com/datagrok-ai/public/issues/1983): Categorical coloring: colors cannot be changed.
  * [#1986](https://github.com/datagrok-ai/public/issues/1986): Dragging a column to filter panel moves it to the beginning of the table.
  * [#1208](https://github.com/datagrok-ai/public/issues/1208): Default tooltip is unexpectedly reset on duplicating a view
* For [Filter Panel](../../visualize/viewers/filters.md) fixed:
  * [#1230](https://github.com/datagrok-ai/public/issues/1230): Filters: custom filters not synced.
  * [#1609](https://github.com/datagrok-ai/public/issues/1609): FilterGroup.add issue, which caused error message when called through `grok.events.onContextMenu` context.
  * [#1984](https://github.com/datagrok-ai/public/issues/1984): Filter's missing values settings are not properly synced between different tabs/views.
  * [#1987](https://github.com/datagrok-ai/public/issues/1987): Tooltip from the wrong column is shown in filters in some cases.
* [Scatter plot](../../visualize/viewers/scatter-plot.md):
  
  * [#2089](https://github.com/datagrok-ai/public/issues/2089): Added an option to change Scatter plot and Trellis plot's mouse drag action from panning to selection. For trellis plot with an inner scatter plot viewer, the default value is `selection`.
  *  [#2108](https://github.com/datagrok-ai/public/issues/2108): Improved behavior of 'Dot' marker type on Scatter plot.
  *  [#2055](https://github.com/datagrok-ai/public/issues/2055): Implemented indicating zoom state for categorical axes.
* [Trellis plot](../../visualize/viewers/trellis-plot.md):
  *  For trellis plots with an inner PC plot, hidden labels and limited the default number of columns to four.
  * Categorical axes now show a tooltip for small structures, providing clear rendering on the axis. 
  * [#1952](https://github.com/datagrok-ai/public/issues/1952): Improved Trellis plot behavior to exclude displaying empty categories.
  *  [#2144](https://github.com/datagrok-ai/public/issues/2144): Added saving of Trellis inner viewers' states when switching between them.
  * [#2054](https://github.com/datagrok-ai/public/issues/2054): Zooming is now disabled for inner scatter plot, bar chart and box plot viewers. For scatter plot, it can be enabled via a new "Allow zoom" setting.
  * Fixed:
    * [#1957](https://github.com/datagrok-ai/public/issues/1957): Trellis plot: cannot be added in specific case, error on setting structure column as X axis.
    * [#2053](https://github.com/datagrok-ai/public/issues/2053): Trellis plot: tooltip in inner viewer cannot be changed.
    * [#2116](https://github.com/datagrok-ai/public/issues/2116): Column settings are connected on two Trellis plots after settings are picked up/applied.
* [Line Chart](../../visualize/viewers/line-chart.md)
  * [#2049](https://github.com/datagrok-ai/public/issues/2049): Added the  min/max properties.
  * Fixed [#2163](https://github.com/datagrok-ai/public/issues/2163): Line chart lines are displayed incorrectly in some cases (drop below the X axis).
* [Box Plot](../../visualize/viewers/box-plot.md)
  * [#2057](https://github.com/datagrok-ai/public/issues/2057): Added structures rendering on the X axis.
  * Fixed [#1985](https://github.com/datagrok-ai/public/issues/1985): Box plot: Y axis zoom slider behavior is inverted.
* For [correlation plot](../../visualize/viewers/correlation-plot.md), fixed [#2009](https://github.com/datagrok-ai/public/issues/2009): Correlation Plot shows wrong axes.

### Data Access

*  Implemented the ability to call query without loading from server, which leads to acceleration of each query by 200-300ms.
* Added the ability to create hierarchies in the query tree using the `friendlyName`.
*  Added the ability to expand complex structures into columns for Neptune and Neo4j databases. With this update, you no longer need to list node properties explicitly, making query writing much simpler and quicker.
* Fixed:
  * [#1963](https://github.com/datagrok-ai/public/issues/1963): Snowflake Schema View random order for columns.
  * [#1750](https://github.com/datagrok-ai/public/issues/1750): File does not delete from Home Directory in File Manager.

### Improvements for developers

*  Added column tags: `.includeInCsvExport` and  `.includeInBinaryExport`.
*  Added `Property.validateValue(x)`.
*  Added Functions: `PadLeft` and `PadRight(s, length, char)`.
*  Implemented the ability to specify input directly via the `forInputType` method.
*  For Inspector, added dynamic package-originated panels with `meta.inspectorPanel: true`.
*  For Inspector, added the  ability to inspect widgets.
*  Moved adHoc property to FuncCall.

#### [JS API](../../develop/packages/js-api.md)

*  [#2118](https://github.com/datagrok-ai/public/issues/2118): Added the API configuration to define layout import settings.
*  [#1592](https://github.com/datagrok-ai/public/issues/1592): For packages, introduced decorators for functions with special roles (viewers, importers, etc.).
*  For JS API, added `ProgressIndicator.canceled`.
* For FileInfo, exposed `isFile` and `isDirectory` methods to JS API.

### Enhancements in packages

* [Bio](https://github.com/datagrok-ai/public/blob/master/packages/Bio/CHANGELOG.md): separator support for **Sequence Space** and **Activity Cliffs**, Helm monomer type is usable for **MSA**, and other enhancements. 
* [Charts](https://github.com/datagrok-ai/public/tree/master/packages/Charts/CHANGELOG.md): improvements for Timelines and Sunburst plot, Radar chart fixes.
* [Chem](https://github.com/datagrok-ai/public/blob/master/packages/Chem/CHANGELOG.md): RDKit rendering for various databases (Chembl, ChemblAPI, PubChem, and DrugBank) when OCL is used, **Scaffold Tree** integration into the **Filters Panel**, and other improvements.
* [Curves](https://github.com/datagrok-ai/public/tree/master/packages/Curves/CHANGELOG.md): rendering for curves and confidence intervals, user-defined Javascript function support for running curves with caching, enhanced rendering in small cells, and more. 
* [Demo](https://github.com/datagrok-ai/public/tree/master/packages/Demo/CHANGELOG.md): new [Bioreactors](https://public.datagrok.ai/apps/Tutorials/Demo/Bioreactors), [heatmap](https://public.datagrok.ai/apps/Tutorials/Demo/Visualization/General/Heatmap) and [Chem](https://public.datagrok.ai/apps/Tutorials/Demo/Cheminformatics) demos.
* [Dendrogram](https://github.com/datagrok-ai/public/tree/master/packages/Dendrogram/CHANGELOG.md): ability to select leaves from specific nodes, separate loader view, and ability to switch distance calculation methods for macromolecules, while also supporting the semType `macromolecules`.
* [Helm](https://github.com/datagrok-ai/public/tree/master/packages/Helm/CHANGELOG.md): addressed issues.
* [Peptides](https://github.com/datagrok-ai/public/blob/master/packages/Peptides/CHANGELOG.md):  **Invariant Map** now selects sequences instead of filtering.
* [PowerGrid](https://github.com/datagrok-ai/public/tree/master/packages/PowerGrid/CHANGELOG.md): structure editing support for pinned columns.
* [Tutorials](https://github.com/datagrok-ai/public/tree/master/packages/Tutorials/CHANGELOG.md): new tutorials.

  
### Enhancements in libraries

* A new [math](https://github.com/datagrok-ai/public/tree/master/libraries/math/CHANGELOG.md) library intended to be used for high efficiency (mostly wasm) calculations.
* See other updates for [utils](https://github.com/datagrok-ai/public/tree/master/libraries/utils/CHANGELOG.md), [ml](https://github.com/datagrok-ai/public/tree/master/libraries/ml/CHANGELOG.md) and [bio](https://github.com/datagrok-ai/public/tree/master/libraries/bio/CHANGELOG.md).

## 2023-06-29 1.15.4

* Improvements: 
  * [#1730](https://github.com/datagrok-ai/public/issues/1730) Scaffold Tree: Filters panel integration.
* Fixed:
  * [#2052](https://github.com/datagrok-ai/public/issues/2052): Coloring is missing on adding new table view.
  * [#2061](https://github.com/datagrok-ai/public/issues/2061): DataFrame in viewer property is null on "viewer added" event.

## 2023-06-19 1.15.3

* Improvements:  
  * [#1988](https://github.com/datagrok-ai/public/issues/1988): Trellis Plot: Enhancements:
     * Added a column filter to the properties for column selection. X and Y axes accept categorical columns only
     * Removed the "Viewer Type" property from the context panel (users can switch it via the dropdown in the viewer itself)
     * Added the ability to resize a legend for Trellis
     * Adjusted labels appearance on the Y axis. Use ellipsis for long names
     * Added Show GridLines property (when enabled, makes the inner viewer borders visible)
  * [#1951](https://github.com/datagrok-ai/public/issues/1951): Add Global scale option to Trellis plot
  * [#2011](https://github.com/datagrok-ai/public/issues/2011): Viewers: Inconsistencies in saved layouts won't allow creating a viewer when such layout is applied
  * [#1378](https://github.com/datagrok-ai/public/issues/1378): PC Plot: context menu harmonization: added "..." for the Columns item 
  * [#1778](https://github.com/datagrok-ai/public/issues/1778): Viewers: Organize properties consistently in the context panel: moved **Col Label Orientation** to the **Styles** category for heatmap and grid 
* Fixed:
  * [#1984](https://github.com/datagrok-ai/public/issues/1984): Filter's missing values settings are not properly synced between different tabs/views 
  * [#1985](https://github.com/datagrok-ai/public/issues/1985): Box plot: Y axis zoom slider behavior is inverted 
  * [#1986](https://github.com/datagrok-ai/public/issues/1986): Dragging a column to filter panel moves it to the beginning of the table 
  * [#1609](https://github.com/datagrok-ai/public/issues/1609): FilterGroup.add issue 
  * [#1983](https://github.com/datagrok-ai/public/issues/1983): Categorical coloring colors cannot be changed 
  * [#1987](https://github.com/datagrok-ai/public/issues/1987): Tooltip from the wrong column is shown in filters in some cases 
  * [#1959](https://github.com/datagrok-ai/public/issues/1959): Cannot adjust legend width for a Trellis plot 
  * [#1230](https://github.com/datagrok-ai/public/issues/1230): Filters: custom filters not synced 
  * [#2009](https://github.com/datagrok-ai/public/issues/2009): Correlation Plot shows wrong axes
  * [#2028](https://github.com/datagrok-ai/public/issues/2028): Handling the layout save with a viewer loading error 
  * GROK-12902: Chem | Context Pane | Actions: Add to favorites doesn't work 
  * GROK-13238: Adding new columns fails in the presence of virtual columns 
  * [#1838](https://github.com/datagrok-ai/public/issues/1838): Error on adding a custom viewer 
  * GROK-13321: Filter Panel: 'Unsupported operation: NaN.round()' error occurs when hovering the Board Id tab
  * GROK-13343: Trellis plot is not loading from the layout in some cases


## 2023-06-02 1.15.2

* (Bug) GROK-13157: Can't load table without initial permissions 
* GROK-13075: Delete old and inactive containers 
* (Improvement) GROK-13188: Datlas: Improve migration process 


## 2023-05-25 1.15.1

* (Bug) GROK-13109: Databases | Mysql, MariaDB: Get top 100 on some tables throws an error.
* (Bug) GROK-13157: Can't load table without initial permissions.
* (Bug) GROK-13132: Fixed scatter plot legend.
* (Enhancement) [#1378](https://github.com/datagrok-ai/public/issues/1378): PC Plot: context menu harmonization.
* (Bug) [#1929](https://github.com/datagrok-ai/public/issues/1929): Viewers do not support hex colors set in categorical coloring.
* (Bug) [#1208](https://github.com/datagrok-ai/public/issues/1208): Default tooltip is unexpectedly reset on duplicating a view.
* (Bug) GROK-13126: Grid | Layouts: color-coding is applied from the previous layout.


## 2023-05-17 1.15.0

We've launched a new version of the Datagrok platform 1.15.0. This update introduces various enhancements in computations, platform governance, stability, and overall usability. Some of the key improvements in this release include:
* Namespaces view for easy access to all data sources and content within the platform. Now you can explore and access all available data in one centralized location.
* [EDA package](release-history.md#eda) using partial least squares regression for the multivariate data analysis.
* [Bioreactors package](#bioreactors) for the simulation of the mechanism of Controlled Fab-Arm Exchange.
* Usage Analysis package for studying usage statistics. It enables you to analyze user activity, package distribution, and function usage to gain valuable insights for statistical analysis. To learn more, see [Usage Analysis](https://datagrok.ai/help/govern/usage-analysis#usage-analysis-application).
* [Demo application](https://public.datagrok.ai/apps/Tutorials/Demo), an interactive educational resource showcasing the diverse capabilities and features of the DataGrok platform. It offers tutorials and demonstrations for hands-on learning of data manipulation, visualization, modeling, and more.
* Multiple improvements in plugins, such as  [Chem](#chem), [Peptides](#peptides), [Dendrogram](#dendrogram).

### Visualization and usability improvements

* Added core viewers support for  +/- Infinity.
* Introduced new [help property](https://community.datagrok.ai/t/visualization-related-updates/521/33?u=oahadzhanian.datagrok.ai) for all core viewers. It could be either a markdown, or a URL.
* Optimized the BigInt parsing.
* [Grid](../../visualize/viewers/grid.md):
  * [#353](https://github.com/datagrok-ai/public/issues/353): Added "Row Source" option to the configuration of the Grid visualization.
  * Added the Ctrl+Shift+UP shortcut to sort the current column.
  * [#1860](https://github.com/datagrok-ai/public/issues/1860): Implemented the support of **Apply Coloring** for selected columns of compatible data type.
  * Fixed:
    * [#1839](https://github.com/datagrok-ai/public/issues/1839): Datagrok is sporadically frozen in drag-and-drop hover state (moving columns around, adding to filters, calculated columns dialog).
    * [#1840](https://github.com/datagrok-ai/public/issues/1840): Calculated columns: errors on missing values.
    * GROK-12918: An exception when dragging columns to the first position.
    * GROK-13054: Error when saving the layout.
* [Filter Panel](../../visualize/viewers/filters.md):
  * Added the ability to sort the default filters by #categories.
  * Fixed [#1837](https://github.com/datagrok-ai/public/issues/1837): Filters cannot be enabled if all filters were disabled in another view
* [Scatter plot](../../visualize/viewers/scatter-plot.md):
   * [#1746](https://github.com/datagrok-ai/public/issues/1746): Added the ability to set date for the Min and Max values on axes.
   * [#1882](https://github.com/datagrok-ai/public/issues/1882): Reset the Min and Max values on column change.
   * Fixed:
     * [#1764](https://github.com/datagrok-ai/public/issues/1764): Error loading scatter plot for specific data and log scale axes on applying layout.
     * [#1761](https://github.com/datagrok-ai/public/issues/1761): Out of memory on adding viewer with specific data.
     * [#1744](https://github.com/datagrok-ai/public/issues/1744): Table is unexpectedly filtered if scatter plot has 'filter by zoom' setting and empty column on log scale axis.
     * [#1855](https://github.com/datagrok-ai/public/issues/1855): Some data is missing after setting **Zoom and Filter** to `pack and zoom by filter`.
     * [#1858](https://github.com/datagrok-ai/public/issues/1858) Scatter plot with specific data and log scale axis failed to load on applying layout.

* [Line Chart](../../visualize/viewers/line-chart.md)
  * Fixed:
    * [#1852](https://github.com/datagrok-ai/public/issues/1852): Line chart with splitting with specific data is making Datagrok slow (row selection, interaction with line chart).
    * [#1671](https://github.com/datagrok-ai/public/issues/1671) Line chart connects the first and last values when resizing the window.
* [PC Plot](../../visualize/viewers/pc-plot.md):
   * [#1378](https://github.com/datagrok-ai/public/issues/1378): Provided the context menu harmonization:
      * Added the **Filter** menu item with the **Show** and **Show Filtered Out Lines** entries.
      * Added the **Columns... ** menu item that opens the **Select column** dialog.
      * Added the **Selection** menu item with the **Show Current Line**, **Show Mouse Over Line**, and **Show Mouse Over Row Group**  entries.
* [Box Plot](../../visualize/viewers/box-plot.md)
  * [#1635](https://github.com/datagrok-ai/public/issues/1635): Added a zoom slider to adjust axes range.
  * [#1377](https://github.com/datagrok-ai/public/issues/1377): Added Show Category Axis to the X-axis menu.

### Data Access

* Improvements:
  * Added the ability to set JS postprocessor for query transformation.
  * Added  `is null` and `is not null` operators for patterns. 
  * Created a new operator for datetime pattern `\- last 30 days`, same for hours, weeks, and months.
  * Improved performance for specific queries.
  * Provided cache and streaming compatibility.
  * Improved context menu options for tables:
    * Renamed **Visual Query** to **Aggregation Query** .
    * Renamed **Build Query** to **Join DB Tables**.
  * Renamed  the fields in **Aggregation Query**: 
    * Columns to Pivot 
    * Rows to Group by 
    * Measures to Aggregate 
    * Filters to Filter
* Fixed:
  * GROK-12778: Snowflake: Browse Schema doesn't work.
  * GROK-12780: Snowflake: error on the Content tab.
  * GROK-12782: Snowflake: error when running the query.
  * GROK-12876: ClickHouse: Browse Schema doesn't work.
  * GROK-12880: NULL in inspect menu.

### Enhancements in packages

#### [Chem](https://github.com/datagrok-ai/public/tree/7c62a0c018ec631d3b23760d538a17aaf4d4ca36/packages/Chem#readme)

* Improvements:
   * Performed **Elemental Analysis** refinement. Now the resulting column name includes the name of the column for which Elemental Analysis is performed.
   * Added the ability to choose whether to search for MCS exact atoms or MCS exact bonds when performing R-Groups Analysis.
   * Added the ability to avoid recalculating the similarity search when the current row is changed. To achieve this, we introduced the **Follow Current Row** checkbox in the settings of the similarity search viewer. 
   * For proper handling of properties and rendering, we now check for smarts and molecular fragments separately.
   * Provided the ability to copy data from the **Descriptors** and **Properties** tabs on the **Context Pane**.
   * Moved **Descriptors** and **Fingerprints** from the  **Context Pane**  to the **Top Menu** ( **Chem** > **Calculate**).
   * Added **Substructure Search** to the **Top Menu** ( **Chem** > **Search**)
   * Modified the tooltip for dialog and drag-n-drop to prevent it from overlapping with the data.
   * Implemented RDKit rendering for Chembl, ChemblAPI, PubChem, and DrugBank databases if OCL is used currently.
   * Performed UI polishing. We've adjusted input field names, repositioned elements in certain dialogs, added tooltips, and organized the top menu items into groups: Calculate, ADME/Tox, Search, Analyze, and Transform.  
* Fixed:
   * [#1492](https://github.com/datagrok-ai/public/issues/1492): Elemental analysis: malformed data handling.
  * GROK-11898: Orientation for smiles in Structural Alerts.
  * GROK-12115: Hamburger menu closing while switching a sketcher.
  * GROK-12905: Sketcher is not opening from the **Filter Panel**.
  * GROK-12929: Scripts don't work if called from the package.
  * GROK-12933: Drug likeness: set score precision.
  * GROK-12961: Elemental Analysis: `Unsupported operation: NaN.round()` error on some data.
  * GROK-12962: Similarity Search: doesn't work on malformed data

#### [Usage Analysis](https://github.com/datagrok-ai/public/tree/master/packages/UsageAnalysis#readme)

We've released Usage Analysis 1.0.0, a tool for comprehensive statistics and insights into usage patterns on the Datagrok platform. Gain a deeper understanding of user interactions, make data-driven decisions, and optimize performance to enhance the user experience. To learn more, see [Usage Analysis](https://datagrok.ai/help/govern/usage-analysis#usage-analysis-application).

#### [EDA](https://github.com/datagrok-ai/public/tree/master/packages/EDA)

We've implemented the multivariate data analysis using partial least squares (PLS) regression in the EDA package. Our solution reduces the predictors to a smaller set of uncorrelated components and performs least squares regression on them. To provide high-performance in-browser computations, we use WebAssembly. For details, see [Multivariate analysis](https://datagrok.ai/help/explore/multivariate-analysis/pls).

#### [Bioreactors](https://github.com/datagrok-ai/public/tree/master/packages/Bioreactors#readme)

We've created the Bioreactors package for the simulation of the mechanism of [Controlled Fab-Arm Exchange](https://www.jbc.org/article/S0021-9258(20)40445-4/fulltext). The unique Datagrok WebAutosolver tool provides an in-browser solution for implementing simulations of the considered phenomena. To ensure high performance, we utilize wasm-computations.

#### [PhyloTreeViewer](https://github.com/datagrok-ai/public/tree/master/packages/PhyloTreeViewer#readme)

* We've retired the PhyloTree viewer and removed it from the PhylotreeViewer package.
* Added the ability to reset zoom for PhylocanvasGL viewer.

#### [BiostructureViewer](https://github.com/datagrok-ai/public/tree/master/packages/BiostructureViewer#readme)
Now, when you click on a cell with a 3D molecule, we utilize the Biostructure viewer to display the structure in a separate viewer.

#### [Bio](https://github.com/datagrok-ai/public/tree/master/packages/Bio#readme)

* We have structured the top menu by organizing the items into groups: SAR, Structure, Atomic level, and Search. 
* Fixed GROK-13048: Activity cliffs identification for macromolecules.

#### [Dendrogram](https://github.com/datagrok-ai/public/tree/master/packages/Dendrogram#readme)

* Improvements:
  * Added a separate loader view to the dendrogram.
  * Added the ability to reset zoom for Dendrogram.
  * Moved hierarchical clustering script to the client side.
  * Allowed to switch distance calculation method for macromolecules (Hamming or Levenstein) in case of MSA.
  * Enabled Dendrogram to work with semType macromolecules.

#### [Peptides](https://github.com/datagrok-ai/public/tree/master/packages/Peptides#readme)

* Fixed: 
    * GROK-12424: Logo Summary Table selection.
    * GROK-12793: Logo Summary Table statistics mismatch.
    * GROK-12794: Wrong statistics on the Distribution Panel. 

#### [Tutorials](https://github.com/datagrok-ai/public/tree/master/packages/Tutorials#readme)

Added new tutorials on [grid customization](https://public.datagrok.ai/apps/tutorials/Tutorials/ExploratoryDataAnalysis/GridCustomization) and how to [add a calculated column](https://public.datagrok.ai/apps/tutorials/Tutorials/Datatransformation/CalculatedColumns).

#### [PowerGrid](https://github.com/datagrok-ai/public/tree/master/packages/PowerGrid#readme)

Added structure editing support for pinned columns.

### Improvements for developers

* Improvements:
  * Added the ability to specify input directly via  `forInputType` method.
  * Added filter property (numerical, categorical) for ColumnInput.
  * Added the ability to create Icon by root.
  * Added functions PadLeft and PadRight (s, length, char).
  * Added Property.validateValue(x).
  * Markup: Added the ability to reference commands by their names.
* Fixed:
  * GROK-12768: DateTimeColumn Timezone inconsistency.
  * GROK-12935: DateTime parsing for "5/21/2022 7:47" pattern.
  * GROK-12971: DateTime parsing for "2-May-2007" pattern.
  * GROK-12936: Float Formatter Error.
  * GROK-12787: UI generation: MultiChoice input binding
#### [JS API](../../develop/packages/js-api.md)

* Improvements:
  * Added ProjectSaving and ProjectClosing events to the existing project events.ts in js-api.
  * Exposed additional menu.item options.
  * Implemented TagsInput control.

#### Enhancements in libraries

[utils](https://github.com/datagrok-ai/public/tree/master/libraries/utils#readme): implemented TypeAhead and DropDown controls.
  
### Other bug fixes

* GROK-12430: Scientific numbers do not convert when converting column type from string to double.
* GROK-12472: The **Link Tables** dialog doesn't work.
* GROK-12805: Viewers: `NullError: method not found: 'F$' on null` error when trying to set Color.
* GROK-13046: Incorrect Min/Max work in AddNewColumnDialog.
* [#1866](https://github.com/datagrok-ai/public/issues/1866) Legends in the stack viewers are not saved properly.
* GROK-12912: Sketcher in core: incorrect dialog display.
* GROK-12954: Viewers: the gear icon does not fit on a 1600x1000 px screen at 2 dpi.
* GROK-13027: Fix Molecule3D grid cell renderer in preview.
* GROK-13068: Projects: `Insufficient privileges to save TableInfo` error when uploading.

## 2023-05-09 1.14.4

* (Bug) GROK-12918: Viewers: Grid: An exception when dragging columns to the first position 
* (Bug) [#1855](https://github.com/datagrok-ai/public/issues/1855): "Pack and zoom by filter": some data is missing 
* (Bug) [#1852](https://github.com/datagrok-ai/public/issues/1852): Line chart with splitting with specific data is making Datagrok slow (row selection, interaction with line chart) 
* (Bug) [#1840](https://github.com/datagrok-ai/public/issues/1840): Calculated columns: errors on missing values 
* null 
* (Bug) [#1837](https://github.com/datagrok-ai/public/issues/1837): Filters cannot be enabled if all filters were disabled in another view 
* (Bug) [#1839](https://github.com/datagrok-ai/public/issues/1839): Datagrok is sporadically frozen in drag-and-drop hover state (moving columns around, adding to filters, calculated columns dialog) 
* (Enhancement) [#1860](https://github.com/datagrok-ai/public/issues/1860): Grid: Support "Apply Coloring" for selected columns of compatible data type 


## 2023-05-03 1.14.3

* (Improvement) GROK-12427: Move Chem descriptors and fingerprints to top-menu -> Calculations (WIP)


## 2023-04-28 1.14.2

* (Improvement) GROK-12922: Cache and streaming compatibility 


## 2023-04-20 1.14.1

* (Improvement) GROK-12721: Tooltip modifications for dialog and drag-n-drop 
* (Improvement) GROK-12788: Cat filter: select one on first select, select all on deselecting last 
* (Enhancement) [#1583](https://github.com/datagrok-ai/public/issues/1583): Scatter Plot: "pack and zoom" option should response to changes when a numerical column is used 
* (Bug) [#1758](https://github.com/datagrok-ai/public/issues/1758): Filtering is slow if there are multiple bar charts on several views 
* (Enhancement) [#1778](https://github.com/datagrok-ai/public/issues/1778): Viewers: Organize properties consistently in the context panel 
* (Enhancement) [#1635](https://github.com/datagrok-ai/public/issues/1635): Box Plot: Add a zoom slider to adjust axes range 
* (Enhancement) [#1746](https://github.com/datagrok-ai/public/issues/1746): Scatter Plot: Date picker for X and Y axes min/max 
* (Bug) [#1671](https://github.com/datagrok-ai/public/issues/1671): Line chart: connects the first and last values when  resizing the window 
* (Bug) GROK-12891: File Manager: scrollbar is missing in the File Tree 


## 2023-04-18 1.14.0

We've launched a new version of the Datagrok platform (1.14.0). This release focuses on improving data access speed and convenience, new visualization and usability features, and ensuring platform stability. Here are some of the biggest improvements made in this release:

* New version of Grok Connect with a completely redesigned data delivery system. Thanks to streaming, you can now work with large data volumes and receive data as the query is being executed.
* [Biostructure Viewer](#biostructureviewer) package for interactive exploration of biological structures. For details, see [Biostructure Viewer](https://community.datagrok.ai/t/macromolecules-updates/661/14).
* New user-designed cell forms feature that allows you to include more details in cells, displaying additional values alongside the cell content. In-grid forms also inherit the color-coding scheme of the grid, enabling you to view even more additional information in the cell.
* [Tile Viewer]( https://public.datagrok.ai/apps/Tutorials/Demo/Viewers/Tile%20Viewer) supporting swimlanes mode, where you can drag and drop cards between lanes.
* Completely redesigned [table link editor](../../transform/link-tables.md).
* Improvement in the function execution efficiency with a significant performance boost of 10x, specifically in formula calculations.
* Multiple improvements in plugins, such as [Chem](#chem), [Bio](#bio), [Helm](#helm), [Peptides](#peptides), [Dendrogram](#dendrogram).

### Visualization and usability improvements

* Core viewers now render molecular structures in the legend marking them with their corresponding colors.
* Placed a hint in the form of three dots (...) in the corner of an element to trigger the context menu when clicked.
* Added a scroll bar in the context menu to handle options that do not fit on the screen.
* For all dialogs that require selecting a molecular column, it is now set as the initial option.
* Ability to select and import multiple files at once using the **Open local file** dialog.
* Consistently reorganized viewers' properties in the **Context pane**:
	 * Moved **Label orientation** to the **Style** tab for PC plot and the other viewers with this property.
   * Moved **Segments Font** to the **Style** tab for line chart.
   * Moved tooltip properties from **Misc** to **Tooltip** and column properties from **Misc** to **Columns** for heatmap and grid.
   * Renamed `Point` to `Line` for PC plot
   * Reordered the properties of the viewers in a way that ensures they use the same order of categories.
* [Grid](../../visualize/viewers/grid.md):
  * Ability to add info panels as columns in the grid.
  * Implemented the column-based edit permissions. Now you can set different editing permissions for each column based on user roles or other criteria.
  * Implemented the ability to view the column's private tags by clicking the **Show private** button (**Context Pane** > **Details**).
  * Added the capability for users to design their own tooltips specific to each column.
  * Implemented a warning that displays the reason why a read-only cell cannot be edited when a user attempts to start editing it.
  *  Implemented the ability to specify columns for in-grid [Default Forms](../../visualize/viewers/grid.md#summary-columns) (**Context Menu** > **Add** > **Forms** > **Default HTML**).
  * Improved default scientific format to show two digits after comma.
  * Fixed:
    * The issue where coloring and format changes were resetting unexpectedly upon interacting with filters, if a group tooltip was set
    * Table is not rendered properly after switching dataframe after sorting data in a column
    * An error when double-clicking on a column header.
* [Filter Panel](../../visualize/viewers/filters.md):
  * Added the ability to drag-and-drop custom filters to the filter panel.
  * Now filters use column format in the tooltip for dates.
  * Consolidated filters application across different tabs so that the behavior is intuitive to users. The **Reset** button resets all filters in the analysis:
      * Turns all filters on.
      * Clears all selected options.
      * Clears visualization filters (e.g. scatterplot color label or `filter by zoom`)
      * Turning on / off a filter makes it turned on /off across all tabs
      * Turning off a filter makes the filter inactive but its selection remains set
  * Fixed:
    * The list of filters becomes empty in **Reorder Filters** dialog
    * Filters synchronization when clicking the **Clear** button
    * Filters re-ordering dialog is not opened if a custom filter was added to filters panel
    * The usability of range filters when min/max values existing in the table were shown instead of actually set range
    * Slow filtering if there are multiple bar charts on several views.
* [Scatter plot](../../visualize/viewers/scatter-plot.md):
  * Now it displays legend labels for conditional coloring.
  * Added date picker for X and Y axes min/max values.
  * Harmonized the scatter plot context menu:
    * Renamed **Drop Lines** to **Show Drop Lines**.
    * Renamed **Regression Line** to **Show Regression Line**.
    * Added a checkbox for toggling **Rectangle/Lasso selection**.
    * Added **Label Color As Marker** to **Context menu** > **Markers**.
    * Now, we've disabled marker type when the marker axis is defined.
  * Made some improvements for **Formula Lines** on the scatter plot:
    * Formula lines labels are now parallel to lines
    * You can hide the regression lines formula by toggling the corresponding checkbox in the context menu (**Tools** > **Show Regression Lines**)
  * Fixed:
    * Table unexpected filtering if scatter plot has `filter by zoom` setting and empty column on log scale axis
    * Error loading scatter plot for specific data and log scale axes on applying layout
    * Incorrect response on enabling `pack and zoom` option when a numerical column is used
    * The visibility of the selected dot on zooming when it is outside the currently shown range
    * Error loading scatter plot for specific data on applying layout
    * Zoom slider behavior for inverted axes
    * Disappearing of the legend when specific filter is applied.
* [Histogram](https://community.datagrok.ai/t/visualization-related-updates/521/32):
  * Added the multi-distribution mode: the ability to show markers and a choice of whether or not distributions should be normalized.
  * Added an option to show markers and straight lines for the split mode.
  * Added an option to normalize distribution in  the split mode.
  * Fixed:
     * incorrect behavior when filtering
     * resetting filter when filtering out missing values.
* [Line Chart](../../visualize/viewers/line-chart.md)
  * Now the line chart renders molecule structures on the X axis.
  * Added the zoom slider on the X axis.
  * Fixed:
    * Wrong data is shown on line chart
    * Unsupported operation: NaN.floor() error
    * Incorrect legend for **Split by** when **Row Source** is set to `Selected`
    * Formula lines error
    * Change columns bug.
* For [PC Plot](../../visualize/viewers/pc-plot.md) we've fixed:
   * Density recalculation on normalization or switching columns to log scale
   * Adding Measures in PC plot transformation.
* [Box Plot](../../visualize/viewers/box-plot.md)
  * Added a zoom slider to adjust axes range
  * Harmonized the box plot context menu:
    * Added new menu entries to **Markers**: **Marker Type**, **Marker Size**, **Marker Opacity**, **Show Inside Values**, and **Show Outside Values**.
    * Added menu **Statistics** on the top level with two entries: **Show Statistics** and **Show p-Value**.
    * Added menu Y axis with **Show Value Axis**, **Axis Type**, **Invert Y axis**, and **Show Value Selector**.
    * Added menu X axis with **Show Category Axis**, **Label Orientation**, and **Show Category Selector**.
* [Bar chart](../../visualize/viewers/bar-chart.md)
  * Harmonized the context menu.
  * Implemented a fix that correctly hides the tooltip when the user selects the **Tooltip** > **Hide** option.
*  [Tile Viewer](../../visualize/viewers/tile-viewer.md) now shows long values that are not fully visible in a tooltip.
* Resolved the issue where the  [Trellis Plot](../../visualize/viewers/trellis-plot.md) displayed not all data until a user clicked a bar chart.

### Data Access

* Improvements:
  * Added the support of [CoreWeave Object Storage](../../access/files/shares/coreweave.md).
  * Conducted Postgres harmonization, which includes complete types support and the addition of regex support in the string pattern for the connector.
  * The MySQL connector now fully supports all types.
  * Added the ability to display SQL views along with tables.
* Fixed:
  * Regex for query annotation.
  * Query annotation: Allow extra spaces.
  * **Browse Schema** action for MySQL now returns the name of the current database.
  * Snowflake: Browse Schema doesn't work.
  * Snowflake: error on the **Content** tab.
  * DataQuery: Query created from connection doesn't return a result.
  * Databases | MS SQL | NorthwindTest | Schemas | dbo: 'com.microsoft.sqlserver.jdbc.SQLServerException:'error in the **Inspect** tab for the column.
  * Error on the **Activity** tab for some connections
  * Snowflake handling of doubles, infinity, -infinity, NaN, and large scale floating point numbers.
  * Exception in parametrized queries with two identical parameter names.
  * Databases | PostgresDart | NorthwindTest: error for the `orders` query.
  * Databases | PostgresDart | ChemblSQL: error on the **Content** tab when clicking the `domains` table.
  * Databases: can't add Visual Query.
  * Databases | Build Query: doesn't show query preview.
  * Databases | Visual Query: error when adding a column.

### Enhancements in packages

<!--#### [ApiSamples](https://github.com/datagrok-ai/public/tree/master/packages/ApiSamples#readme)


ApiSamples now includes new control demos and the "hide-empty-columns" script.-->

#### [BiostructureViewer](https://github.com/datagrok-ai/public/tree/master/packages/BiostructureViewer#readme)

The package provides the functionality of two viewers  that enable the visualization of biological structures on the Datagrok platform.  For details, see [Biostructure Viewer](https://community.datagrok.ai/t/macromolecules-updates/661/14)

#### [Bio](https://github.com/datagrok-ai/public/tree/master/packages/Bio#readme)

* Renamed the dialogs and the items in **Top Menu** > **Bio**.
* Fixed:
  * 'NullError: method not found: 'Q' on null' error when running Sequence Activity Cliffs
  * The aggressiveness of the Macromolecule detector was reduced.

#### [Helm](https://github.com/datagrok-ai/public/tree/master/packages/Helm#readme)

* Fixed:
  * `Neither peptide, nor nucleotide` error when calling **Helm Web Editor**
  * Macromolecule detector for alphabet UN.

#### [Chem](https://github.com/datagrok-ai/public/tree/7c62a0c018ec631d3b23760d538a17aaf4d4ca36/packages/Chem#readme)

* Improvements:
   * Renamed the dialogs and the items in **Top Menu** > **Chem**. Specifically, we've removed the **Sketcher** item and added [**ADME/Tox**](https://community.datagrok.ai/t/cheminformatics-updates/457/23?u=oahadzhanian.datagrok.ai).
   * Updated the column and cell properties and now all chemical information on the **Context Pane** is located under the **Chemistry** info panel. Also, we've removed the diversity Search widget from the **Context Pane** > **Details**.
   * Made tooltip modifications for molecular data. Now, the tooltip in dialogs no longer covers the data and shrinks during drag-n-drop.
   * Improved the handling of malformed and empty data. Now a warning appears with the indexes of the failed molecules in chem analyses.
  * Customization of dimensionality reduction algorithms is now possible, especially in the context of **Chemical Space** and **Activity Cliffs** functions. For details, see [Dimensionality reduction algorithms](https://community.datagrok.ai/t/cheminformatics-updates/457/21).
  * Implemented improvements for Similarity and Diversity Search viewers. Information from any column of the initial dataset can now be added to these viewers, and the color coding applied to the initial dataframe is also displayed there. To learn more, see [New functionality in Similarity/Diversity search](https://community.datagrok.ai/t/cheminformatics-updates/457/19).
  * When a large structure is depicted in the sketch box of the filters panel, it is now rendered in a tooltip.
  * [Scaffold tree](../../datagrok/solutions/domains/chem/chem.md#scaffold-tree-analysis) improvements:
    * Scaffold tree now uses the default sketcher set in the Chem package properties.
    * Implemented logical AND, OR, NOT operations for scaffold tree elements.
* Fixed:
   * Structure display in the structure filter.
   * Scaffold Tree fails on a large dataset.
   * Error occurring on hovering the molecular and macromolecular column's header.
   * Infinite updating when switching to the **Chem Draw** sketcher.
   * Filter synchronization between **Marvin** sketcher and the **Filter Panel** while filtering smarts.
   * Handling empty values in **Fingerprints**.
   * **Group Analysis** viewer: can't open properties by clicking the **Gear** icon.
   * `Result of truncating division is NaN` error when hovering over molecular column's header.
   * **Marvin** filtering via `List/NOT List`.
   * Filter synchronization for different tables if they have the same name of the molecular column.
   * Chem tooltip displaying on the **Context Pane**.
   * `Chem | Possibly a malformed molString at init: 'undefined'` error when hovering over molecular column's header.
   * The issue where the **Cancel** button was not resetting the filter when using **Marvin**.
   * Grid scrolling to the right when hovering over a molecular column.
   * Wrong rendering of malformed data.
   * Molfile copying.
   * The issue with the synchronization of the sketcher. Now, when you switch the sketcher on the **Filter Panel**, it changes the sketcher in the **Hamburger Menu** and on the **Context Pane** accordingly.
   * The  **copy as smiles** action.
   * Color coding for the similarity search viewer.
   * **Elemental Analysis**: false message `Charts package is not installed`.
   * The `Sketcher not found` warning on package init.
   * Sketcher popup clipping in the dialog.
   * Empty values handling in **Descriptors** and **Map Identifiers**.
   * Sketcher duplicating in the **Filter Panel**.
   * Browser tab freezing on hovering a molecule column's header after data was removed.
   * Sketcher tooltip after opening a molecule in Marvin.
   * The molblock native wedging usage. Now it's used where available.

#### [Dendrogram](https://github.com/datagrok-ai/public/tree/master/packages/Dendrogram#readme)

The dendrogram viewer now uses standard default colors.

#### [Peptides](https://github.com/datagrok-ai/public/tree/master/packages/Peptides#readme)

We've significantly improved the Peptides package, including a new feature integrating with dendrograms. For more information, see [Peptides updates](https://community.datagrok.ai/t/macromolecules-updates/661/13).

#### [Charts](https://github.com/datagrok-ai/public/tree/master/packages/Charts#readme)

We've added routing to the viewers' demo.

#### [GIS](https://github.com/datagrok-ai/public/tree/master/packages/GIS#readme)

We've reduced the aggressiveness of detectors and added a column name check.

#### [Tutorials](https://github.com/datagrok-ai/public/tree/master/packages/Tutorials#readme)

 You can now open and share a tutorial via a link.

#### [PowerGrid](https://github.com/datagrok-ai/public/tree/master/packages/PowerGrid#readme)

We've now exposed the pinned columns functionality in the PowerGrid package and fixed the following issues:
  * Unclickable Id link on pinned columns
  * Incorrect displaying after applying the layout with the viewer stacked over the table.

### Improvements for developers

We've separated **Images** and **Docker Containers** in the platform.

#### [PowerPack](https://github.com/datagrok-ai/public/tree/master/packages/PowerPack#readme)

We've added `Qnum` to supported column types and fixed the issue with the execution of some functions in the **Add new column** dialog.

#### [Viewers](../../develop/how-to/develop-custom-viewer.md)

* Added canvas grid cell renderers to the gridext library to render multiple values in a grid cell using a vertical layout.
* Cell renderer now can estimate the desired cell size based on the entire column

#### [JS API](../../develop/packages/js-api.md)

* Improvements:
  * Introduced a new feature for application developers. Now you can [place hints](https://community.datagrok.ai/t/javascript-api-updates/526/21) for users in the form of indicators, and popups, as well as describe a series of visual components in the wizard.  UI methods now come equipped with interactive hints, similar to those found in tutorials, that can be attached to various elements. With the addition of these methods, application authors are able to incorporate interactive hints in their app code, allowing them to introduce new features to their users, among other things.
  * To work with custom filters, we've added `.custom-filter-type` and `.ignore-custom-filter` tags . See [column tags](../../visualize/viewers/filters.md#column-tags) to learn more.
  * Provided `typeAhead`, `dropDown` and `breadcrumbs` controls for the platform.
  * Implemented the capability to export graphics of scatter plot to JS API. You can now render scatter plot viewer to an arbitrary graphics context using JS API. Basically, we've made a viewer that works directly with the cell renderer. This means you can now show cells as viewers and manipulate them in other ways too.
  * Added classes for form viewers in the API: `Form` and `FormViewer`.
  * Added the `onViewRemoving` event.
  * Improved the table join feature by adding the option to merge "all columns" without having to list each individual column.
  * Added `DataQuery.connection`.
  * Added JS renderers: `getDefaultSize(gridColumn)`.
  * Exposed the row tooltip of ScatterPlot viewer to JS API
  * Added the `tags` field to `GridColumn` for storing serializable JS objects.
  * Added the `autoHideTabHeader` property to **Accordion**.
  * Balloon: added the ability to pass elements as content
  * Exposed the `Menu.show` options.
  * Grid: exposed `gridRowToTable` and `tableRowToGrid` methods.
* Fixed:
  * Fixed the TableInput reacting to the programmatic value setting.
  * TableInput does not fire onInput when user uses **Open file** button.
  * The ` whereRowMask()` function to correctly apply the mask to the selection.
  * Signatures of `_render` function for all plots with `CanvasViewerMixin`.
  * Color calculation when the data range is 0.

### Other bug fixes

* Console error on closing viewers.
* Browser tab freezes on hovering a molecule column header after data was removed.
* Floating viewer appears unexpectedly on applying layout in some cases.
* Out of memory on adding viewer with specific data.
* While dragging windows, sometimes text in them might get selected.
* Group Analysis viewer opens twice in the demo.
* Nested context panel's expanded state is not remembered.
* Incorrect bounds of `Datetime` pattern with operator `this week`.
* 3d scatter plot: failed to load if **Show Axes** property is set to `false`.
*  Issue with formula line opacity: opacity set in one formula line object also applies to the others.
* issue with the condensed form, where the icons for options and postfixes overlapped.
* `RegExpExtract` function throwing an exception when no match is found.
* Wrong docker status when making a request.
* Databases: error in the **Sharing** tab of the **Context Pane** after deleting the shared query.

