# 2021-10-14 build 0.94.0

We've released a [new version](https://public.datagrok.ai/) of the Datagrok platform (0.94.0). This release focuses on the areas of visualization, usability, rich features for the developers, and traditionally on the platform stability.

## Visualization and usability improvements

* Many updates in Scatter plot: a [lasso tool](https://community.datagrok.ai/t/extensions-to-the-scatter-plot-viewer/481/5), axis sliders, data split by [marker](https://community.datagrok.ai/t/visualization-related-updates/521/20) [shape](https://community.datagrok.ai/t/visualization-related-updates/521/22), custom lines by equations; all summarized in [this video](https://www.youtube.com/watch?v=Q3Dn5NSDSEY&t=3018s)
* New Histogram features: splines and bands: [video](https://www.youtube.com/watch?v=Q3Dn5NSDSEY&t=3490s)
* New Line Chart features: marker shape selection, whiskers: [video](https://www.youtube.com/watch?v=Q3Dn5NSDSEY&t=3562s)
* Using “Relative Values” property in combination with the “Stack” property to analyze the distribution of the stacked values: [overview](https://community.datagrok.ai/t/bar-chart-updates/580)
* Many aspects and issues of the Bar Chart addressed: [overview](https://community.datagrok.ai/t/piechart-issues/569/3)
* Search for category names in filters (useful when there are many categories in the column): more in [this video](https://www.youtube.com/watch?v=Q3Dn5NSDSEY&t=3642s) + an [overview](https://community.datagrok.ai/t/visualization-related-updates/521/21)
* Lists in parameterized database queries: [overview](https://community.datagrok.ai/t/lists-in-parameterized-database-queries/597) + [video](https://www.youtube.com/watch?v=meRAEF7ogtw); also check the video about parameterized queries in general: [link](https://www.youtube.com/watch?v=sSJp5CXcYKQ)
* Database query caching with cron jobs ability: [video](https://www.youtube.com/watch?v=Q3Dn5NSDSEY&t=3984s)

## Improvements for developers

* More granular event handling (an example with intercepting conditional color coding settings change: [link](https://community.datagrok.ai/t/javascript-api-updates/526/12))
* Creating drag-and-drop objects: [overview](https://community.datagrok.ai/t/javascript-api-updates/526/13)
* Getting cell colors: [overview](https://community.datagrok.ai/t/javascript-api-updates/526/14)
* Adding custom machine learning models to Datagrok in R, Python, or from external models deployed in clusters: [video](https://www.youtube.com/watch?v=G66MN30ZPGQ), [package](https://github.com/datagrok-ai/public/tree/master/packages/CustomML), [help](https://datagrok.ai/help/learn/custom-machine-learning-models)

## Enhancements in public packages

* [PowerPack and Universal Search](https://github.com/datagrok-ai/public/tree/master/packages/PowerPack) — commonly used platform enhancements, currently covering Start Page widgets, Power Search (ability to search for anything from the Start Screen), and search templates support: [video](https://www.youtube.com/watch?v=Q3Dn5NSDSEY&t=86s)
* Updates in the [Peptides](https://github.com/datagrok-ai/public/tree/master/packages/Peptides) package, learn more in [the video](https://www.youtube.com/watch?v=Q3Dn5NSDSEY&t=2390s)

We have also redesigned our approach to interactive Tutorials. We have moved them to a standalone [package](https://github.com/datagrok-ai/public/tree/master/packages/Tutorials) and made them available through a dedicated [app](https://public.datagrok.ai/apps). Check [this video](https://www.youtube.com/watch?v=Q3Dn5NSDSEY&t=1920s) to get a grasp. Every tutorial on the platform now is a simple TypeScript code, like [here](https://github.com/datagrok-ai/public/tree/master/packages/Tutorials/src/tracks/chem/tutorials). Now it is possible to equip customers' apps and deployments of Datagrok with tutorials as well!

# 2021-07-29 build 0.93.0

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

# 2021-05-06 build 0.91.10

We've just released a [new version]([https://public.datagrok.ai](https://public.datagrok.ai/)) of the Datagrok platform (0.91.10). It is a major release with multiple stability and performance improvements, many new features, and new APIs for developers. Here are some of the notable advances:
* Viewers improvements, including [trellis plot](https://datagrok.ai/help/visualize/viewers/trellis-plot) usability and performance, [bar chart](https://datagrok.ai/help/visualize/viewers/bar-chart) with DateTime categorization, [conditional color coding](https://community.datagrok.ai/t/table-visualization-options/333/9)
* [Custom cell renderers](https://community.datagrok.ai/t/cheminformatics-updates/457/8) for molecules and other object types as a first-class citizen in viewers' axis, tooltips, tiles and forms
* New viewers, such as [Chord](https://community.datagrok.ai/t/visualization-related-updates/521/15) and [Timelines](https://community.datagrok.ai/t/visualization-related-updates/521/4)
* Cheminformatics support for [client-side RDKit](http://rdkit.blogspot.com/2019/11/introducing-new-rdkit-javascript.html) searches and molecule rendering
* Usability and stability improvements in [filtering](https://community.datagrok.ai/t/visualization-related-updates/521/12), projects, layouts, files browsing and sharing
* Numerous developers' enhancements, including TypeScript support and new APIs
* [Custom comparer functions](https://community.datagrok.ai/t/javascript-api-updates/526/9)

Check [release notes](https://datagrok.ai/help/develop/release-history) for more details, and give the new version a try at [https://public.datagrok.ai](https://public.datagrok.ai/)!

# 2021-04-14 build 0.89.36

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

# 2021 Future Release

