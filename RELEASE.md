# 2021-07-29 Build 0.93.0

## Highlights

We've released a [new version](https://public.datagrok.ai/) of the Datagrok platform (0.93.0). It's a large release with many new features and stability improvements, addressing both user experience and flexibility for developers.

Here are some interesting areas in the platform updates:

* Viewers improvements, such as [relative values in Bar Charts] (https://community.datagrok.ai/t/bar-chart-updates/580) or [multi-axis Line Chart](https://community.datagrok.ai/t/visualization-related-updates/521/18)
* Visual improvements like [dragging columns with drop zones](https://community.datagrok.ai/t/visualization-related-updates/521/17) and [cleaner columns summaries in tooltips]( https://community.datagrok.ai/t/visualization-related-updates/521/16)
* Functional features, such as [filtering multivalue cells](https://community.datagrok.ai/t/visualization-related-updates/521/15)
* Improvements in [sharing files and folders](https://www.youtube.com/watch?v=0QxzllnBreI&t=3895s)
* Better programmatic operations, such as [more granular event handling for metadata](https://community.datagrok.ai/t/visualization-related-updates/521/15) and
[new customizations in Scatter Plot](https://community.datagrok.ai/t/javascript-api-updates/526/11)
* Numerous improvements in database connectors stability and performance
* Additional built-in functions and improved [documentation on them](https://datagrok.ai/help/transform/functions/datetime-functions)

We also worked intensively on our [public packages](https://github.com/datagrok-ai/public), presenting a few new of them:

* Biosignals – Physionet annotations file viewer,
dynamic scripting capabilities for constructing signal processing pipelines: [overview](https://github.com/datagrok-ai/public/tree/master/packages/BioSignals#readme), [video](https://www.youtube.com/watch?v=0QxzllnBreI&t=1932s)
* ClinicalCase — working with clinical data in SDTM format: [overview](https://github.com/datagrok-ai/public/tree/master/packages/ClinicalCase#readme), [video](https://www.youtube.com/watch?v=2xuxJjpjXi4&t=95s)
* Sequence Translator — oligonucleotide sequences conversion (BioSpring, Axolabs): [overview](https://github.com/datagrok-ai/public/tree/master/packages/SequenceTranslator#readme), [video](https://www.youtube.com/watch?v=2xuxJjpjXi4&t=3782s)
* Oligo Batch Calculator — extinction coefficient, mass, optical density, Nmole, mol. weight: [overview](https://github.com/datagrok-ai/public/tree/master/packages/OligoBatchCalculator#readme), [video](https://www.youtube.com/watch?v=2xuxJjpjXi4&t=4902s)
* Cheminformatics features: [chemical dataset curation](https://community.datagrok.ai/t/cheminformatics-updates/457/11), [activity cliffs](https://www.youtube.com/watch?v=2xuxJjpjXi4&t=1933s), [save as SDF](https://community.datagrok.ai/t/cheminformatics-updates/457/10)

## Major features and improvements

* Serialization info panel for columns
* Ability to share files
* Ability to skip DF reading on server after Grok Connect request
* Add New Column: ColumnGrid Widget, minor performance improvements
* Add Repository: Add attribute "name" for Source Type selector
* Add string name indexing for columns
* Adde pubspec to speed-up packages resolving
* Bar Chart: adaptive font size; automatically zoom in to a reasonable number of categories in case of many categories; make scrolling smoother
* Binning Functions: BinBySpecificLimits
* Box Plot: better border color; showCategorySelector and showValueSelector properties
* Conditional color-coding: bring options from the property panel to column context menu
* Connectors: Redshift: Schema browsing; default schema to providers and use as condition for schema browsing
* Correlation plot: on cell click, show the corresponding scatter plot in the property panel 
* Grid: ShowVisibleColumnsInTooltip property; automatically pick up cell type from dataframe
* Histogram and BarChart: improved washed-out default colors
* JS API: Added grok.shell.startUri
* JS API: DG.TAGS.FORMAT tag as a preferred way to set the column format
* JS API: Func.appy(params): Promise<TResult>
* JS API: Func.find(package, name, tags, returnType)
* JS API: GridColumn.getVisibleCells()
* JS API: ScatterPlot: add viewport
* JS API: Stats.histogramsByCategories
* JS API: TabPane: new properties: root, content, parent
* JS API: Viewer.onContextMenu event 
* JS API: add event onAccordionConstructed and Accordion.context getter 
* JS API: add grok.shell.views and grok.shell.tableViews
* JS API: dialog.clear()
* JS API: toDart(): ability for JS classes to define custom conversion
* JS API: ui.icons for commonly used icons
* JS API: ui.image
* JS API: ui.tools.setHoverVisibility(host, elements)
* Join tables: Don't include key fields by default
* Line Chart: MultiAxis: double-click (choosing a primary series) should have no effect; move hit-testing functionality to series renderers; support for multiple Y axes
* Math Functions: Log(arg1, arg2), Log10(), Median(), RandBetween(), Round10(), etc.
* Multi-value filters
* Query Builder dialog: Add "name" attributes to checkboxes for tables
* ScatterPlot custom renderer: disable overlay dot rendering on categories hover; changed back marker color; set default marker size to 10
* Scripting: Ability to set custom name for parameter; ability to set postfix for parameter; show execution results in the property panel
* Table View: Columns pane: add search
* Text Functions: RegExpExtract(), RegExpReplace(), StrFind(), etc.
* Trellis Plot: ability to enlarge individual in-trellis viewers
* UI: Ability to set dialog background
* UTC support in the datetime parser
* Viewers: implement 'dashboard' style for all standard viewers; support for styles
* datagrok-tools: allow skipping questions in `grok config`

## Bugs

* Add New Column: Dragging functions opens a drop-area for the table
* Add New Column: History bug \- inputs are not filled with history
* Add new column: Strings in functions cannot be enclosed in single quotes
* Bar Chart: Unexpected bar color change after filtering
* Bar Chart: Viewer coloring settings should take precedence over the grid coloring settings (WIP)
* Bar Chart: coloring cannot be disabled; coloring gets only applied when editing colors in a grid column
* Box Plot: an exception when stdev(value) = 0
* Box Plot: improve initial choice of value column (stdev > 0 if possible)
* Chem: R-Groups Analysis with client RDKit
* Column format changes are not persisted in layouts
* Connectors: Impala: int32max instead of real values
* Connectors: Oracle: NullPointerException in DB table content
* Core: Bitset.falseCount returns the number of set bits
* Core: default tooltip config is no longer saved with layout
* Current user is system while deploying
* Custom ML: apply function from different ML engine
* Custom ML: apply function from different ML engine
* Data | Unpivot does not work
* DataQuery with choices throws an exception
* DataQuery with choices throws an exception
* Dataframe: Detect column max significant digits in CSV loading
* Datagrok to Python skips blank lines
* Events: onViewerAdded and onViewerClosed are sent twice
* Excel import: an empty column is created as part of the dataframe
* Extra space breaks function annotation
* File Sharing doesn't work
* Filters: Multi-value filters does not turn off when corresponding checkbox is off
* Filters: There are no icons for sorting and searching on hover for Multi-value filters
* Filters: adding returns a filter for the column added previously instead of the currently selected one
* Filters: clicking on "search" should open search field AND focus on it
* Filters: clicking on "search" should open search field AND focus on it
* Filters: filter component state is different across pages
* Filters: if multi-value filters are present in the panel, the reset button doesn't work
* Filters: multi-value filters have no square indicator on top to toggle category selection
* Filters: range slider filters out nulls
* Fixed scrollbars hover
* Functions View: Click on the selected category does not remove the check mark
* Functions: function search won't work on packages
* Grid: drag-and-drop column reordering: provide drop zones
* Grid: empty space on the right
* Grid: html cells are not re-rendered after sorting a column
* Grid: switching global coloring on / off removes linear color-coding
* Grok connect: "No suitable driver found for..."
* Grok connect: "The method execute() cannot take arguments" error with query parameters
* Grok connect: DriverManager returns wrong driver
* Grok connect: NullPointerException with meta.cache option
* Grok connect: can't browse scheme of CompoundLookup
* Hide function deselects the selected columns from the tooltip
* Incorrect QNum parsing
* IntColumn.fromList(values) does not work with values outside of the int32 range
* JS API: Columns.byTags does not work
* JS API: JsViewerHostCore is returned instead of Viewer instance
* JS API: Label breaks layout of TextInput with icon
* JS API: Label breaks layout of TextInput with icon
* JS API: OnDialogClosed event fires twice
* JS API: Properties cannot be changed in JsViewer (JsViewer.props and JsViewer.setOptions result in errors)
* JS API: toJs won't work on GrokPackage
* JS API: ui.info for a yellow info bar
* JS Editor: exception when ApiSamples package is not there
* JS Viewer: Incorrect height calculation
* JS Viewers: error: NullError: method not found: 'where$1' on null
* JS: grok.functions.call should return JS object for multiple output parameters
* Layouts: Grid looses event handlers after layout restore
* Line Chart: "axes follow filter" feature does not work
* Line Chart: NullReferenceError when changing X axis column
* Line Chart: point hit-testing doesn't work for points on the right
* Matcher: matching on multiple criteria ignores the "and/or" option
* Packages: beta flag
* Pie Chart: unnecessary datetime aggregation
* Prevent socket memory consuming
* Queries: "Converting object to an encodable object failed" with dataframe as parameter
* Query runs forever
* Query-driven dashboards: query controls do not show up when a project is open
* Radiobutton throws exception when created
* Scatter Plot: "axes follow filter" feature does not work
* Scatter Plot: Regression line appear without activation
* Share button doesn't work for projects
* UI: new view appears below the current
* Updated ui.md (ui.stringInput, ui.searchInput)
* Updated ui.md (ui.stringInput, ui.searchInput)
* View Layouts: Error balloon after deleting saved layout
* Viewers: Inconsistent column selection inside a viewer and its properties panel
* Viewers: textColor property misspelling
* Viewers: the legend colors are not synchronized
* Viewers: the menu item `Viewer` is not visible in uploaded projects
* Viewers: Fixes regression in 2D layout alignment of unknown origin
* grok.dapi.projects.where('bad filter').first() returns a Project instance with d == null

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

# Future release

Multiple marker shapes https://community.datagrok.ai/t/visualization-related-updates/521/20
JS API making drag-and-droppable objects https://community.datagrok.ai/t/visualization-related-updates/521/15?u=dskatov

## Highlights

* A new "Add New Column"
* Custom ML
* Function View and Sensitivity Analysis https://www.youtube.com/watch?v=2xuxJjpjXi4&t=2507s

## Major features and improvements

## Bugs
