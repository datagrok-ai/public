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
