# 2021-07-29 Build 0.93.0

We've released a [new version](https://public.datagrok.ai/) of the Datagrok platform (0.93.0). It is a large release with many new features and stability improvements.

## Highlights

## Major features and improvements

* Ability to share files 
* Connectors: Redshift: Schema browsing
* Trellis Plot: ability to enlarge individual in-trellis viewers 
* Join tables: Don't include key fields by default 
* Multi-value filters 
* Connectors: add default schema to providers and use as condition for schema browsing 
* Table View: Columns pane: add search 
* JS API: DG.TAGS.FORMAT tag as a preferred way to set the column format 
* JS API: add event onAccordionConstructed and Accordion.context getter  
* Grid: ShowVisibleColumnsInTooltip property 

* Scripting: Ability to set custom name for parameter 
* Scripting: Ability to set postfix for parameter 
* UI: Ability to set dialog background
* JS API: add grok.shell.views and grok.shell.tableViews 
* JS API: ScatterPlot: add viewport 

* Grid: automatically pick up cell type from data frame
* JS API: Stats.histogramsByCategories
* Correlation plot: on cell click, show the corresponding scatter plot in the property panel 

* ScatterPlot custom renderer: disable overlay dot rendering on categories hover  

* JS API: Func.find(package, name, tags, returnType) 
* JS API: Func.appy(params): Promise<TResult> 
* JS API: toDart(): ability for JS classes to define custom conversion 
* JS API: ui.icons for commonly used icons 
* JS API: ui.tools.setHoverVisibility(host, elements) 
* JS API: ui.image 
* Scripting: interactivity: show execution results in the property panel 
* UTC support in the datetime parser
* JS API: dialog.clear() 
* Add string name indexing for columns 

* Conditional color-coding: bring options from the property panel to column context menu 

* Binning Functions: BinBySpecificLimits 
* DateTime extractor functions (year, month, day, etc) 

* Math Functions: Log(arg1, arg2) 
* Math Functions: Log10() 
* Math Functions: Round10() 
* Math Functions: RandBetween() 

* JS API: Added grok.shell.startUri 
* JS API: TabPane: new properties: root, content, parent 
* Math Functions: Median() 
* Text Functions: StrFind() 
* Text Functions: StrLeft() 
* Text Functions: StrRight() 
* Text Functions: StrRepeat() 
* Text Functions: RegExpReplace() 
* Text Functions: RegExpExtract()
* JS API: Viewer.onContextMenu event  

* 'Serialization' info panel for columns 
* Add New Column: minor performance improvements 
* JS API: GridColumn.getVisibleCells() 
* Functions: Parameter categories 
* Viewers: support for styles 
* Viewers: implement 'dashboard' style for all standard viewers 
* Histogram and BarChart: washed-out default colors 
* Line Chart: move hit-testing functionality to series renderers 
* Introduced Color.textColor 
* Bar Chart: adaptive font size 
* Bar Chart: in case of many categories, automatically zoom in to a reasonable number of categories 
* Bar Chart: make scrolling smoother 
* ScatterPlot: changed back marker color 
* ScatterPlot: set default marker size to 10 


## Bugs

* (Bug) Chem: R-Groups Analysis with client RDKit 
* (Bug) Incorrect QNum parsing 
* (Bug) JS API: Columns.byTags does not work 
* (Bug) Queries: "Converting object to an encodable object failed" with dataframe as parameter 
* (Bug) Excel import: an empty column is created as part of the dataframe 
* (Bug) Grid: html cells are not re-rendered after sorting a column
* Prevent socket memory consuming
* JS API: ui.info for a yellow info bar
* (Bug) Core: default tooltip config is no longer saved with layout 
* (Bug) JS Viewer: Incorrect height calculation 
* (Bug) Bar Chart: coloring cannot be disabled 
* (Bug) Filters: adding returns a filter for the column added previously instead of the currently selected one 
* (Bug) Bar Chart: coloring gets only applied when editing colors in a grid column 
* Packages: beta flag 
* (Bug) UI: new view appears below the current 
* Dataframe: Detect column max significant digits in CSV loading
* (Bug) Viewers: the menu item `Viewer` is not visible in uploaded projects 
* (Bug) Grid: drag-and-drop column reordering: provide drop zones 

* fixes regression in 2D layout alignment of unknown origin 
* (Bug) Events: onViewerAdded and onViewerClosed are sent twice 
* (Bug) Filters: if multi-value filters are present in the panel, the reset button doesn't work 
* (Bug) Radiobutton throws exception when created 

* (Bug) Pie Chart: unnecessary datetime aggregation 
* (Bug) JS API: toJs won't work on GrokPackage 
* (Bug) JS API: JsViewerHostCore is returned instead of Viewer instance 

* (Bug) Current user is system while deploying 
* (Bug) Hide function deselects the selected columns from the tooltip  
* (Bug) Bar Chart: Viewer coloring settings should take precedence over the grid coloring settings (WIP)
* (Bug) Data | Unpivot does not work 
* (Bug) Column format changes are not persisted in layouts 

* (Bug) Connectors: Oracle: NullPointerException in DB table content  
* (Bug) IntColumn.fromList(values) does not work with values outside of the int32 range 
* (Bug) Query runs forever 

* (Bug) Grid: empty space on the right 
* (Bug) Filters: filter component state is different across pages 
* (Bug) Connectors: Impala: int32max instead of real values 

* (Bug) Filters: range slider filters out nulls 
* (Bug) OnDialogClosed event fires twice 
* (Bug) File Sharing doesn't work 
* (Bug) JS Editor: exception when ApiSamples package is not there 
* (Bug) Functions: function search won't work on packages 
* (Bug) Share button doesn't work for projects 

* Fixed scrollbars hover 
* (Bug) Grok connect: can't browse scheme of CompoundLookup  
* (Bug) Add new column: Strings in functions cannot be enclosed in single quotes 
* Box Plot: showCategorySelector and showValueSelector properties 
* (Bug) Box Plot: an exception when stdev(value) = 0 
* (Bug) Box Plot: improve initial choice of value column (stdev > 0 if possible) 
* Box Plot: better border color 
* (Bug) Functions View: Click on the selected category does not remove the check mark 
* JS: grok.functions.call should return JS object for multiple output parameters 
* (Bug) Extra space breaks function annotation 
* (Bug) grok.dapi.projects.where('bad filter').first() returns a Project instance with d == null 
* (Bug) Layouts: Grid looses event handlers after layout restore 
* (Bug) Filters: Multi-value filters does not turn off when corresponding checkbox is off 
* (Bug) Filters: multi-value filters have no square indicator on top to toggle category selection 
* (Bug) Filters: There are no icons for sorting and searching on hover for Multi-value filters 
* (Bug) Filters: clicking on "search" should open search field AND focus on it 


# 2021-05-06 Build 0.91.10

We've just released a [new version]([https://public.datagrok.ai](https://public.datagrok.ai/)) of the Datagrok platform (0.91.10). It is a major release with multiple stability and performance improvements, many new features, and new APIs for developers. Here are some of the notable advances:
* Viewers improvements, including [trellis plot](https://datagrok.ai/help/visualize/viewers/trellis-plot) usability and performance, [bar chart](https://datagrok.ai/help/visualize/viewers/bar-chart) with DateTime categorization, [conditional color coding](https://community.datagrok.ai/t/table-visualization-options/333/9)
* [Custom cell renderers](https://community.datagrok.ai/t/cheminformatics-updates/457/8) for molecules and other object types as a first-class citizen in viewers' axis, tooltips, tiles and forms
* New viewers, such as [Chord](https://community.datagrok.ai/t/visualization-related-updates/521/15) and [Timelines](https://community.datagrok.ai/t/visualization-related-updates/521/4)
* Cheminformatics support for [client-side RDKit](http://rdkit.blogspot.com/2019/11/introducing-new-rdkit-javascript.html) searches and molecule rendering
* Usability and stability improvements in [filtering](https://community.datagrok.ai/t/visualization-related-updates/521/12), projects, layouts, files browsing and sharing
* Numerous developers' enhancements, including TypeScript support and new APIs
* [Custom comparer functions](https://community.datagrok.ai/t/javascript-api-updates/526/9)

Check [release notes](https://datagrok.ai/help/develop/release-history) for more details, and give the new version a try at [https://public.datagrok.ai](https://public.datagrok.ai/)!

# 2021-04-14 Build 0.89.36

In this release, we've focused on enriching both the experience of the platform end-users and on-platform developers, as well as on connecting these two groups. For instance, custom-built cell renderers, useful for displaying molecules, nucleotide sequences or experiments results, now expand the platform in many places such as [tooltips and tile viewers](https://community.datagrok.ai/t/cheminformatics-updates/457/8). There are also [new viewers](https://community.datagrok.ai/t/visualization-related-updates/521/4), [several bar chart](https://community.datagrok.ai/t/bar-chart-color-by-category/516) [features](https://community.datagrok.ai/t/visualization-related-updates/521/10), [color coding features](https://dev.datagrok.ai/js/samples/grid/color-coding-conditional), and a few dozen of [JS API]((https://datagrok.ai/help/govern/audit#javascript-api)) [enhancements](https://community.datagrok.ai/t/javascript-api-updates/526/8) and [bug fixes](##bug-fixes).

## Highlights

* Rendering of chemical molecules with RDKit [in tooltips, forms, viewers' axis and tile viewer](https://community.datagrok.ai/t/cheminformatics-updates/457/8)
* Bar Chart: support for DateTime [columns categorization](https://community.datagrok.ai/t/visualization-related-updates/521/10), functional improvements
* [Conditional color coding](https://dev.datagrok.ai/js/samples/grid/color-coding-conditional) for the grid, scatter plot, box plot
* New [Timelines Viewer](https://community.datagrok.ai/t/visualization-related-updates/521/4)
* JS API improvements: Typescript support, new [events](https://community.datagrok.ai/t/javascript-api-updates/526/8) and [methods](https://community.datagrok.ai/t/javascript-api-updates/526/5)
* Better JS editing in the [JS fiddle](https://public.datagrok.ai/js) with IntelliSense and `async/await`

## Major features and improvements

* Grid: Custom cell renderers, including RDKit molecules, [in several contexts](https://community.datagrok.ai/t/cheminformatics-updates/457/8): tooltips, tile viewer, form, other viewers' axes
* Bar Chart: [categorizes DateTime columns using functions of Year, Month, Quarter, Year - Quarter, etc.](https://community.datagrok.ai/t/visualization-related-updates/521/10)
* Bar Chart: [add setting for excluding null values category on bar segments](https://community.datagrok.ai/t/bar-chart-color-by-category/516)
* Heatmap: [Improved adaptive rendering of the column names in the heatmap](https://community.datagrok.ai/t/visualization-related-updates/521/9)
* JS API: Migrate to Typescript
* JS API: [Log API](https://datagrok.ai/help/govern/audit#javascript-api)
* JS API: Grid [onColumnResized / onRowsResized events](https://community.datagrok.ai/t/javascript-api-updates/526/8)
* JS API: [dynamic resolution of object handlers](https://dev.datagrok.ai/js/samples/ui/handlers/dynamic-resolving)
* JS API: Logging: New API for collecting telemetry [logging](https://datagrok.ai/help/govern/audit#javascript-api)
* JS Fiddle: Added IntelliSense [dynamic resolution of object handlers](https://dev.datagrok.ai/js/samples/ui/handlers/dynamic-resolving)
* JS fiddle: Use async/await inside 
* Color coding improvements: categorical color coding with values binning, [conditional color coding](https://dev.datagrok.ai/js/samples/grid/color-coding-conditional) for the grid, scatter plot, box plot
* [`.tags` and `.temp` now support JS-native Map-like iteration and modification](https://community.datagrok.ai/t/javascript-api-updates/526/5)
* [New](https://community.datagrok.ai/t/javascript-api-updates/526/7) [`ValueMatcher`](https://dev.datagrok.ai/js/samples/data-frame/value-matching/value-matcher)
* New [Timelines Viewer](https://community.datagrok.ai/t/visualization-related-updates/521/4) becomes available
* Trellis Plot: show more granular X and Y axis ticks on line charts
* [Selection of data subsets to display on scatter plots](https://community.datagrok.ai/t/visualization-related-updates/521/8)
* Filter state is preserved when closing filter panel
* Chem: 10x faster CoordGen-based alignment, flag to hide isotopic labels in RDKit rendering

## Bug fixes

* Histogram: [duplicated Y axis](https://community.datagrok.ai/t/cannot-easily-change-the-aggregation-of-the-histogram/509/4)
* User can't login with OpenId if someone shared something to him by email
* Detectors load on each table open
* dapi.log.events returns empty array
* Layout apply doesn't work on renamed columns
* Chem: Column width with structures is not correct when using OCL
* `files` method `rename` doesn't work properly
* Saved function annotation doesn't delete old parameters from database
* Settings: apply drops password parameters
* tabControl with ui.wait will concatenate contents of tabs and show loader after second tab click
* Package Properties ignored on repository publish
* Grid: Selection display is drawing incorrect on high res.displays (MacOS)
* Grid: empty cells are colored with conditional color-coding on
* error: NullError: method not found: '_ddt$_name' on null
* Visual Query does not work for all supported providers
* Fixed the bool input debouncing issue
* Heat Map: columns do not resize 
* Exception after adding filters on Table View 
* Column format changes are not persisted in layouts 
* S3 AES: Cant post a file
* Ability to disable routing
* UI boolInput: onChange event executes code twice
* Filters panel opens with an error message

## Miscellaneous

* Newick file viewer
* Improved comments on DataFrameViewer::onFrameAttached(..)
* JS API: pass JsViewer in onContextMenu events
* JS API: Ability to create ViewLayout from ViewState
* JS API: preserve metadata while changing column type (WIP)
* Timelines: Add Reset View item to context menu
* datagrok-tools: add package init function template
* /info/packages route added
* Improved bitset notification logic to respect current level of notifications (aka updateLevel)
* Filters: deletion of rows in grid/df results in incorrect filters
* Grid: Harmonize popup menu (WIP)
* Grid: Custom HTML-based cell renderers: investigate performance issues (WIP)
* "Normalize" action does not work
* Improve error handling mechanism for visualizations
* HeatMap: automatically adjust header font size as column width gets smaller
* Query View: Incorrect display of line numbering
* UX: Flickering of scroll bars fixed
* Added missing react import
* Heat Map: automatically adjust header font size as column width gets smaller 
* JS API: onPackageLoaded event
* Ability to discover and choose JS-based grid cell renderers based on a sem type 
* Optimization: faster BitSet construction based on logical function and two arguments
* Grid: make the color of selected rows lighter
* Bar Chart: ability to select categories by clicking on the category label 
* Bar Chart: tooltips on category labels 
* Bar Chart: disable zoom on the X axis
* [Filters: split categories containing a list of values by pipes](https://community.datagrok.ai/t/visualization-related-updates/521/11)
* Conditional Color Coding: Ranges "20-30" and ">40" are always automatically determined for any numeric columns

*Note:* this summary accumulates updates after 0.89.27.
