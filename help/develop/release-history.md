<!-- TITLE: Release History -->
<!-- SUBTITLE: -->

# 2021-02-10 Dev build 0.89.17

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.17`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Programming exercises (WIP)
* Updated exercises.md 
* User Data Storage overview 
* Update user-data-storage.md 
* Fixed links in user-data-storage.md 
* Js Packages: Replace GUID hash with build number
Packages: delete debug version when user publishes release 
* Packages: delete debug version when user publishes release 
* Chord Diagram (WIP)
* Wiki: applications - general, building, debugging. wip 
* Wiki: replaced outdated links 
* Wiki: how-to build an app. wip 
* Wiki: debugging packages. 
* Wiki: How to Build Custom Viewers (WIP)
* JS API: added accordion persistence key 
* Simplified phrasing and linked to the main article 
* Better phrasing 
* DataFrame: JSDocs for append and unpivot. 
* JS API: DataFrame.append: ability to specify columns and rows to append (WIP)
* Harmonize ui framework (WIP)
* Wiki: how-to build an app. WIP: dataframe, credentials 
* Wiki: Project types description 
* Domains / NLP (WIP)
* Wiki: Data Access 
* Harmonize ui framework
(Bug) JS API: column name search box disappears from ColumnComboBox (WIP)
* Changed Redo shortcut to Ctrl+Shift+Z 
* Wiki: switch to recommended package naming 
* (Bug) Multiple column selector: mouse-driven selection is broken 
* Property panel title: made the text user-selectable 
* Wiki: URL fix 
* Wiki: how-to build an app. WIP: restructuring 
* Wiki: how-to build an app. WIP: Semantic annotation and metadata 
* Charts package (WIP)
* Wiki: how-to build an app. WIP: Dataframe aggregations 
* Filters: Wrap filter blocks 
* JS API: Column.init(value | indexToValueFunction) 
* (Bug) Add New Column failed to calculate a formula if all output values are null 
* Wiki: how-to build an app. WIP: User data storage 
* datagrok-tools: add `suffix` key for package version hash 
* Fixing grok_compute_test.dart: more descriptors from the latest rdkit. 
* CVM: fixing a missing kNN. 
* Wiki: how-to build an app. WIP: Subscribing to events 
* Wiki: how-to build an app. WIP: Scripting 
* (Bug) Restore CVM to normal (WIP)
* (Bug) @AutoProperties(includeGetters: true) won't work
(Bug) Notebooks: New notebooks are created empty 
* Add Chem package properties 
* Single scaffold alignment mode 'chem-scaffold': simplified 
* (Bug) DockManager: width rounding error 
* [(Bug) Bar Chart: null values break color-coding of stacked bars](https://community.datagrok.ai/t/bar-chart-color-by-category/516|https://community.datagrok.ai/t/bar-chart-color-by-category/516) 
* (Bug) Scaffold highlights have offsets 
* ClinicalCase - a plugin for analyzing SDTM-based clinical studies data (WIP)
* Chem, RDKit-based (WIP)
* [Histogram: add a dropdown for 'Color Aggr Type' property](https://community.datagrok.ai/t/cannot-easily-change-the-aggregation-of-the-histogram/509+|https://community.datagrok.ai/t/cannot-easily-change-the-aggregation-of-the-histogram/509) 
* Filters: improve sizing
Filters: add DatePicker 
* Fixed dev documentation 
* Show package fullName 
* Wiki: how-to build an app. WIP: REST endpoints, groups, sharing 
* Wiki: minor updates 
* added GrokJsObject to Data/TableQuery as Mixin + defined className 
* Wiki: how-to build an app. WIP: Cleaning up 
* TestData.demog: provided a seed value to the random() functions to make results reproducible 
* (Bug) Multiple selection and cloning will corrupt data frame 
* Line chart: hide X and Y axes in a trellis plot when a chart is too small 
* (Bug) Trellis Plot: the Y axis configuration of inner viewer is not shown 
* (Bug) Grid: hiding columns leaves empty space on the right (WIP)
* Biosignals (WIP)
* (Bug) Viewers | Properties: Sliders move to the next row if width of the second properties is not enough
(Bug) Table | Tooltip: After sketching form, fill of the values ​​in tooltip becomes white and their font changes 
* (Bug) Viewers | Properties:  With long table name, the right column of properties is greatly expanded in width 
* Ability to store data in a layout (WIP)
* Clearing grid's background (so that it will stay white even if cell renderer fails) 
* (Bug) Table info panels are duplicated on the PP 
* Add column filter options are partially hidden 
* (Bug) JS API: toJs doesn't seem to work on DataQuery objects 
* (Bug) Opening a dialog resets the property panel 
* Charts: Timelines (WIP)


# 2021-01-27 Dev build 0.89.16

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.16`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Converted to ui.onSizeChanged 
* (Bug) Parser grammar: Power operation goes before multiplication 
* Clean coordinate data by a render-through-smiles tag (WIP)
* Auto release notes 
* Render an envelope instead of an empty cell for failed molecules 
* (Bug) capitalizeWords('showXAxis') should be 'Show X Axis' 
* Viewers: Scatter Plot: axis-specific context menus 
* Programming exercises (WIP)
* Update exercises.md 
* Updated fetchProxy examples 
* Simplified package `test` function 
* Fixed templates in datagrok-tools 
* Chem, RDKit-based (WIP)


# 2021-01-26 Dev build 0.89.15

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.15`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* [Ability to call the table variable from "Add new column"](https://community.datagrok.ai/t/test-topic-for-integration-testing/426) 
* Auto generated release-notes 
* Molecules rendering in TableViews (WIP)


# 2020.12.03 Stable Version

In this release, we have introduced new [filtering functionality](https://www.youtube.com/watch?v=GM3XixUFFUs&t=2688s) 
and exposed plenty of methods to our JavaScript API. In addition to general improvements of the platform's core, 
Datagrok has been extended with packages in various domains, 
from [cheminformatics](https://community.datagrok.ai/t/cheminformatics/457) 
to [biosignal processing](https://community.datagrok.ai/t/biosignals-package/443) and
 [natural language processing](https://www.youtube.com/watch?v=GM3XixUFFUs&t=94s)

## Latest Docker Images

* Datagrok (new): 
  *  `docker pull datagrok/datagrok:1.0.88-451665c`
  *  `docker pull datagrok/datagrok:stable`
  
* CVM (new): 
  *  `docker pull datagrok/cvm:1.0.88-451665c`
  *  `docker pull datagrok/cvm:stable`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* [Form viewer is now available in tooltips](https://community.datagrok.ai/t/using-tile-viewer-for-tooltip/424)
* [Qualified numbers support for forms](https://community.datagrok.ai/t/some-viewers-cant-be-added/470)
* JS API: DataFrame in cells
* Save column format as part of the layout
* ScatterPlot: do not apply "filter by zoom" _while_ zooming
* Adjust columnsInput to return column name
* [JS API: Viewer.type throws an error](https://community.datagrok.ai/t/viewers-type-error/444)
* [Viewer properties: Collapse non-important sections by default](https://community.datagrok.ai/t/organization-of-viewer-property-panels/446)
* Stata file handler
* Chem: use LRU molecule cache to speed up rendering
* JS API: an efficient LRU cache
* [Grid: Not able to render custom html in cell if format tag is set](https://community.datagrok.ai/t/not-able-to-render-custom-html-in-cell-if-format-tag-is-set/451)
* JS API: add showPopup() function
* JS API: add event on viewer axis click
* JS API: CsvImportOptions
* Package deletion issue
* Projects: after opening a project with a table without view, the automatically added view is not displayed in Scratchpad
* Projects: if you open a project that has a table view layout, another empty view will be added to Scratchpad
* Grid: content won't be copied between browser tabs
* Grid: support for pasting to/from Excel
* Saving settings breaks password fields
* URI should correspond to an app, not a package
* Scatter Plot: sometimes the plot gets zoomed out while panning
* [Scatter Plot: log scale tickmarks improvements](https://community.datagrok.ai/t/scatterplot-issues-requests/434)
* JS API: Ability to dock viewers
* Mime-type sporadically comes as plain text
* OpenId OAUTH Iframe support
* Trellis Plot: Scatter plot is always displayed for entire table
* [JS API: add getters and events for GridColumnList, GridColumn, and Events](https://community.datagrok.ai/t/requests-for-js-properties/441)
* Bar Chart: an error when a bool column is selected for split
* JS API: DockManager improvements
* Cell renderers only work in related viewers the second time
* JS API: add dirty flag cleaning method
* Error after opening a sample table from a script
* [Virtual columns are broken](https://community.datagrok.ai/t/virtual-columns/382)
* JS API: add onEvent method for viewers
* JS API: ui library: ability to easily pass css parameters
* JS API: Automatically get an iterator
* JS-based cell renderers: ability to specify default width and height
* JS-based cell renderers only work after the table is opened for the second time
* JS-based semantic type detection: short notation
* JS API: add grok.dapi.logTypes
* JS API: Add ColumnList.numerical(), ColumnList.byTags()
* JS API: Ability to add Viewer as instance
* Package files API endpoint must return correct HTTP codes
* SAS importing: use R's haven instead of sasb7dat
* JS API: ability to save / apply filter state
* JS API: add ColumnList.replace()
* Grid: mouse wheel event is overridden even when you can't scroll up/down anymore
* JS API: add HtmlTable with remove method
* Correlation Plot: Column names are displayed incorrectly
* DbContext leaks connections
* Usage Analysis: add a dashboard with the most important metrics 
* [Grid: setting gridColumn.visible = true does not unhide it](https://community.datagrok.ai/t/visible-property-of-a-gridcolumn/429)
* Queries: Parameterized queries do not return a dataframe
* Queries: After execution, the "Refresh" button does not work in the table view
* DB Tables: Inspecting the content of tables and columns does not work
* Package file caching
* Fire d4-viewer-added event on adding default view
* Incorrect layout serialization
* Python date format
* JS API: Ability to get package version
* Call functions from JS packages in JS clients
* DataQueries with non-table results don't work
* JS API: Ability to order categories in a string column
* Columns | Sort: Exception when calling dialog for columns with int, double, datetime or bool types
* JS API: Ability to retrieve oauth token
* OpenID Auto-login
* Filters: ability to filter missing values on/off
* Filters: ability to select categories by Ctrl, Shift, or Ctrl+Shift-clicking
* Viewer won't update on underlying dataframe changed
* Stochastic Proximity Embedding improvements
* Pie Chart: clicking on "no value" category returns 0 rows even if there are empty values
* Pie Chart: "no value" color in the legend differs from the one on the pie chart
* Grid: column auto-sizing does not work well with multi-line values
* Packages initialization hangs on error
* [Grid: scrollToPixels does not properly work](https://community.datagrok.ai/t/scrolling-of-grid-content/342/8)
* Grid: resizing columns when the grid is horizontally scrolled does not work
* OpenID: Check token audience
* IP inside docker should be 127.0.0.1
* JS API: grok.shell.addTableView is asynchronous inside
* JS API: Ability to create a TableView without adding to shell
* Filters: tooltip harmonization (available actions, counts)
* Load indicator stuck
* Filters: selection indication
* JS API: Add dialog.close()
* Filters: Categorical: use column header for sorting
* Filters: Categorical: a column with the filtered count
* Filters: Histograms: show filtered-out values in gray
* Filters: Categorical: ability to sort by name or filtered count
* Filters: categories indicator: [Select all / Deselect all / Invert all] context menu
* Filters: "active filters" indicator: "Close others" action
* Filters: an icon to expand / collapse all
* Filters: global indicator: tooltip
* Filters: categorical: in filtering mode, show blue checkboxes
* Filters: categories: color-fill based on the filtered ratio
* Filters: check boxes for category selection
* Oracle NUMBER type handling
* Normalization before clustering works well on dev/public, but throws exception locally on some columns
* Generic Stochastic Proximity Embedding: all data types support
* JS API: Add Column.fromBitSet() constructor
* Grid: onCellValueEdited event
* Filters: ability to turn the whole filter group on/off
* Filters: ability to turn individual filters on/off
* JS API: Chem: Sketcher: ability to get molfile representation
* JS API: RangeSlider: add onValuesChanged event
* JS API: map Dart iterable to JS iterables
* Add normalization option in clustering 
* Layouts: match columns by semantic type
* Filters: top panel: ability to search filters
* Filters: top panel: filter count indicators
* Search boxes: ESC should clear input
* Histograms: ignore margins between bins for hit testing
* JS API: Ability to control xp.tables
* Filters: Numerical: "filtering" indicator
* Filters: Categorical: 'selected categories' indicator
* Filters: Categorical Filters: search box
* Filters: headers
* Histograms: changing any visual setting via property panel triggers filter application
* [JS API: add getters and setters for Range Slider](https://community.datagrok.ai/t/scrolling-of-grid-content/342/4)
* JS API: Add first() aggregation
* Compatibility between npm and Datagrok naming conventions
* Scripting: Created scripts do not open
* Add selection of distance metrics for clustering
* Aggregation: Do not add a space if second part of a column name is an empty string
* OpenId: pick key by kid instead of first
* datagrok-tools: do not require dist for old packages (without webpack.config.js)
* Packages: When function condition is set to another function call -- name is wrong
* JS API: unpivoting
* CSV import: column names are not trimmed
* JS API: metadata handling for data frames and columns
* JS API: BitSet.fromBytes(buffer, bitLength)
* JS API: BitSet: similarity functions
* JS API: creating a column from objects
* JS API: BitSet.fromString('100101')
* DataFrame: virtual columns
* JS API: Support FileInfo & file type


# 2020.08.11 Stable version

## Latest docker images

* Datagrok (new): 
  * `766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok:1.0.82-75b821b` 
  *  `docker pull datagrok/datagrok:1.0.82-75b821b`
  *  [download](https://dev.datagrok.ai/docker_images/datagrok-1.0.82-75b821b.tar)
* CVM (new): 
  * `766822877060.dkr.ecr.us-east-2.amazonaws.com/cvm:1.0.82-75b821b` 
  *  `docker pull datagrok/cvm:1.0.82-75b821b`
  *  [download](https://dev.datagrok.ai/docker_images/cvm-1.0.82-75b821b.tar)
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed issues

* Keep scroll positions on view change
* Box Plot: an exception when category column is not specified
* Viewers: ability to specify the table to bind to
* Shape Map: allowPanZoom property not having any effect
* Restart isolate after error
* If all DB connections are closed with socket error they will not be reconnected
* Sharing a package fails second time
* Shape Map: not rendered when initialized from the saved project 
* entity.save adds _1 to script name if there is another version of a package available
* Disable login form after first click
* Empty values in a float column become zeros when converted to int
* JS API: functions/info-panels/info-panels.js: Exception after opening info panel from the example
* Notebooks don't open
* Use npmignore and gitignore to skip files during upload
* JS API: expose methods for changing column type
* View.close() doesn't fire onViewRemoved event
* Unable to run R script from JS if there is DataFrame input parameter
* API to push credentials to DataConnection or Package
* Scripts to return list of Py installed packages
* Box Plot: an exception when selecting the <empty> x category
* Excel import: treat "N/A" as missing values
* Missing "..." after "Edit" in the context menu for packages and repositories
* Harmonize JS-API ViewLayout class
* JS async functions handling
* Support for qualified numbers
* Support package.ts file
* Print additional information on deploy
* Ability to keep package function log
* datagrok-tools: exclude .~ files from uploading
* S3: .~lock.tree.a# filename breaks S3 storage
* Js Function log
* 404 at credentials-storage help
* Scripts do not run from the console (Namespace error)
* Add cash framework
* Notebooks: Notebook to Script 
* Chem | Similarity Search: Similar structures are not displayed after the progress bar disappears. If you click on the smiles, an exception will be thrown
* JS API Examples: "menu.js" does not work
* If you add Tile Viewer and open the "Scatter Plot" panel for "demog" table, then extra balloons will be displayed
* Tile Viewer: Viewer is displayed empty after editing its form
* Rewrite datagrok-upload to nodeJs
* Support typical shortcuts in Jupyter Notebook
* Viewer embedding doesn't work
* Ambiguous column name exception during package creation
* VirtualItemView: horizontal scrolling mode
* RDKitDemo App: Exception after "Close All" when the application is open
* SDTM App does not work
* Discovery App: Exception after opening
* JS API Examples: "virtual-view.js" does not work
* Chemprop: Does not train on "demo/smiles.csv"
* Files: "Image Classification" panel does not work
* SPARQL connections are deployed without sharing
* ODATA provider visible only with admin mode enabled
* Scripting: DateTime issue
* Apps | RDKitDemo: Application does not open by URL https://dev.datagrok.ai/apps/RDKitDemo
* Projects with Notebooks do not open
* Notebooks: Community demo
* Notebooks: Default name as table name
* Selenium: Create project with tests for viewers serialization
* Predictive Modelling: Improve engine naming
* Tile Viewer: In form edit mode, the "Close and Apply" does not close view
* OpenAPI: Open connection in PP after import
* Upgrade Maps viewer
* Ability to execute "models" on server
* H2O: support list parameters
* Predictive Modelling: create diagrams of train and apply process
* Add AM/PM with space case to date time parser
* Investigate possibility of WebAssembly  usage
* Improve table reading in Jupyter Notebooks


# 2020.07.15 Stable version

## Latest docker images
* Datagrok (new): 
  * `766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok:1.0.80-71eb62e` 
  *  `docker pull datagrok/datagrok:1.0.80-71eb62e`
  *  [download](https://dev.datagrok.ai/docker_images/datagrok-1.0.80-71eb62e.tar)
* CVM (new): 
  * `766822877060.dkr.ecr.us-east-2.amazonaws.com/cvm:1.0.80-71eb62e` 
  *  `docker pull datagrok/cvm:1.0.80-71eb62e`
  *  [download](https://dev.datagrok.ai/docker_images/cvm-1.0.80-71eb62e.tar)

## Addressed issues

* RDKitDemo App: Exception after "Close All" when the application is open
* SDTM App does not work
* Create new test with new structure for Visual Query for all supported providers
* Error when adding LeafletViewer: Failed assertion: boolean expression must not be null
* JS API Examples: "virtual-view.js" does not work
* Predictive Modeling: Update Chemprop version
* Ignore .git on package upload
* Predictive Modelling: Chemprop fails
* _1 suffix on Smiles project
* SPARQL connections are deployed without sharing
* Notebooks: Missing "Sharing" tab on PP
* ODATA provider visible only with admin mode enabled
* View switcher: show all views in it, but highlight the current one
* Oracle: "northwind" connection does not work on dev setup
* JS API: Stack overflow after setting ViewBase name
* Selenium: Change tables names after Visual Query, Build Query, Get All (100)
* ChemPredict Demo Application
* Notebooks: Creating a new notebook doesn't work
* Notebooks: Missed toolbar buttons after open from URL
* Selenium: Bind the name of closure icon to the parent element of the ribbon
* Data Queries: After running queries with same names but with different namespaces, "_#" is added to name in URL
* New data connections are not created
* Package is visible only for deployer
* In open file try to find a connection entity, but not all
* Query View: Exception after opening, missing view name, buttons and menu
* Ability to change project.is_root from UI
* Add a number if there is another entity with the same name in the project
* Unable to call JS function from top menu
* Selenium: vertica-query-test.side does not work
* Selenium: notebook-test.side does not work
* JS API Samples: scripting.js example does not work
* Visual Query: There are no sections "Actions" and "Columns" on the toolbox
* JS API Examples do not open by URL
* Async package functions support
* Visual query: After deleting the used column, the preview is not updated
* Query View not displayed under relevant sidebar section when open
* Selenium: data-query-test.side: After saving query, an exception is thrown
* Neo4J queries do not work
* Selenium: Tables export test: Selenium clicks on items with names and nothing happens
* Oracle provider does not work
* Repositories are not published
* AddressToCoordinates does not work
* JS API: Add TreeView
* Register JS function in top menu using function header
* File connector: Indexing doesn't work
* File connector: Indexing run is not displayed for just created file connection
* JS Examples: "parameterized query.js" does not work
* ChemPredict demo
* Projects do not open by URL
* Upload: Data Sync does not turn on
* Reset password link includes 8080 port
* JS Credentials storage example
* Scripts do not run from the console
* Apps: Applications do not work
* Custom npm registry and proxy settings
* Ability to skip webpack building on server
* Empty error message when publishing webpack package
* Use session token hash instead of token on the server side 
* Ability to add service user
* Files: REST method to upload user file
* Data sync: Project does not open if the table was created by "Get All"
* Hide service projects from Projects View
* Demo Notebooks: Notebooks that contain tables do not work
* Packages: Change servers to global variables
* Chem Panels: "Solubility Prediction" panel does not work (Unable to get project asset)
* FilesView: Folder first sorting
* Data Query: Missing parameter list for parameterized queries on PP
* Files: After renaming connection, its namespace remains old
* Data Sync: Share connection with project
* Selenium: Create ui-test for upload project with Notebook
* File Shares: If connection was not indexed, then in the "Index files" field show corresponding job for running
* Visual Query: Exception after opening view
* Ability to run docker for local debugging
* Make namespaces permanent
* "Home" folder should be indexed by default
* Query Transformations: Addition of transformations is not saved
* Reopening a project with a file share-based table source does not work
* Query Transformations: If transformations view is open, an exception is thrown after "Close All"
* Routing: When you open Open | Databases, the URL is dev.datagrok.ai/connect, which when you refresh the page is redirected to Datagrok page
* Table View: Table name is displayed with extra spaces
* Files browser: Exception is if you click on an empty space in tree if you have an extended connection here
* Data Queries: If the "Transformations" tab was opened, the query is not saved
* S3: Edit connection: bucket name is ignored
* Share connection along with query sharing
* Notebook: Return table to UI
* Events documentation, events handler code generation, expose global event bus
* 'Refresh' icon for entity-based views
* Files: New File Share - newly created file share is not shown in the tree
* Scripts: Scripts outputs do not appear in the console if it was not open before running
* Table columns: "quality" tag is not saved if you set it custom manually
* DB Table Queries: After executing "Get All", "Get Top 100" and "Build a query" the table has name "result"
* Benchmarks: "Benchmark Encoders" and "Benchmark Encoders + Archive" do not work
* Upload dialog: Large indent between tables
* Files View Url to connection not handled
* Share dialog: User permissions selector is in the wrong place
* Join Tables: Selecting or Filtering not matching rows selects the wrong rows
* SPGi: After opening project, URL contains /SupData instead of /Main
* Jira web API queries does not work
* Help view: Incorrect tree style on toolbox if you go to the help page from search (Alt + Q)
* File shares: integrate with the URL
* ProjectUploadView: Tooltip on sync switch
* ProjectUploadView: Remove activity section
* Ability to use both indexed and online file browsing
* Table refreshing doesn't work
* Chem | Map identifiers: Can't find "chembl" project
* Routing for uploaded applications no longer works
* Environments: Name of env should use name from YAML file if it is exist there
* Project renaming
* GrokConnect: Oracle issue with "clob" for some cases
* InfoPanel "scatter-plot.grok" does not work
* Scripting: GrokScript: After running "Template" script, table in the selector is duplicated if run again
* Rename user project on user renaming
* Audit: add "start" event
* Project renaming doesn't affect Name, only friendlyName
* DbContext leaks connections
* Data | Aggregate Rows: Missing first() function in measures list (in menu under "+" | Aggregation)
* Grid: If you close PP, then the grid properties do not appear after clicking on the gear
* Credentials Management Service
* Landing: Slides do not change, if you donвЂ™t first switch it by clicking on the title (on mobile)
* Load packages files with ?hash
* Users: Open views remain after logout
* Context help: Wrong page is displayed after clicking on the "?" on entity deletion dialog
* Predictive models: "Train" button is available if features are not selected
* Chem | Diversity Search: If opened during processing in Chem | Similarity Search, the structures will not be displayed
* Amazon: Chem | Map Identifiers does not work (not connected unichem db)
* Uploading projects: picture is not generated
* Impact analysis
* Bootstrap-backed package that demonstrates grok features embedded into the custom UI
* Test case: job-editor-test
* JupyterLab integration


# 2020.06.16 Stable version

## Latest docker images
* Datagrok (new): 
  * `766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok:1.0.79-93dc7bd` 
  *  `docker pull datagrok/datagrok:1.0.79-93dc7bd`
  *  [download](https://dev.datagrok.ai/docker_images/datagrok-1.0.79-93dc7bd.tar)
* CVM (new): 
  * `766822877060.dkr.ecr.us-east-2.amazonaws.com/cvm:1.0.77-c1d42b4` 
  *  `docker pull datagrok/cvm:1.0.79-93dc7bd`
  *  [download](https://dev.datagrok.ai/docker_images/cvm-1.0.79-93dc7bd.tar)

# 2020.05.27 Stable version

## What's new

### JavaScript API changed. Use `datagrok-api` NPM package as a code reference.

Now API have 3 entry points: 
* `DG` contains complete API code
* `grok` is the main entry point to start using API for most of the tasks
* `ui` is for constructing UI

Here is how code changes:
```
grok.t => grok.shell.t
grok.v => grok.shell.v
grok.tables => grok.shell.tables
grok.presentationMode => grok.shell.presentationMode
grok.parseCsv => grok.data.parseCsv
GrokPackage => DG.Package
TYPE_FLOAT => DG.TYPE.FLOAT
```

Please, refer to the [JavaScript examples](https://public.datagrok.ai/js) to migrate your code.
See also:
* [Developing packages](https://datagrok.ai/help/develop/develop) 

### New scripting demos

  * [image classification](https://github.com/datagrok-ai/public/blob/master/packages/PythonScripts/scripts/image_classification.py)
  * [cell imaging segmentation](https://github.com/datagrok-ai/public/blob/master/packages/PythonScripts/scripts/cell_imaging_segmentation.py)

###File metadata extractors

  * [Tika](https://github.com/datagrok-ai/public/tree/master/packages/Tika)
  * [EXIF](https://github.com/datagrok-ai/public/blob/master/packages/PythonScripts/scripts/exif.py)
  
### First class support of Command Line Interface tools for Linux in scripting:
  * [CLI](https://github.com/datagrok-ai/public/tree/master/packages/CLI)
  * [Tika](https://github.com/datagrok-ai/public/tree/master/packages/Tika)
  
### [ChemSpace integration](https://github.com/datagrok-ai/public/tree/master/packages/Chemspace)

## Latest docker images
* Datagrok (new): `766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok:1.0.77-c1d42b4` [download](https://dev.datagrok.ai/docker_images/datagrok-1.0.77-c1d42b4.tar)
* CVM (new): `766822877060.dkr.ecr.us-east-2.amazonaws.com/grok_cvm:1.0.77-c1d42b4` [download](https://dev.datagrok.ai/docker_images/grok_cvm-1.0.77-c1d42b4.tar)

## Addressed issues

* Scripting/Notebooks: change "localhost" on docker to right datlas address
* JS API: provide JSDoc-style comments for all properties and methods
* Optimize upload utility
* Unable to load JS example by URL
* Models trained on Caret server do not apply to the table
* Ability to download grok API
* Swagger queries deleted on second package upload
* ChemSpace integration
* JS API: manipulate.js sample does not work
* Provide a registry for canvas-based grid renderers
* JS API: Grid: canvas-based event-driven rendering
* JS API: convert all event management to streams / RxJS
* JS API: expose internal storage buffers
* JS API: support Property.choices
* Ability to get function progress from server
* Convert repository publish command to function
* Check database connection before start
* Amazon S3 Storage, optimize directory listing
* Connection Tree: Exception is thrown if you expand broken DB-connection
* GrokConnect: Add default aggregations for SQL-based providers
* Package upload utility: Upload doesn't affect sources list
* Retired EventSource in favor of EventBusProvider.
* Ability to autorefresh file from share on project open
* Generic view persistence mechanism
* Models on h2o server do not train (on dev setup)
* Modeling: AUC estimation failed on Caret
* Files View: null reference error
* Command-line utilities demo
* File Editors: support "sas" format
* File Browser: preview for empty folders does not work
* Check if DB is populated before deploy
* File Shares: indexing doesn't work
* Tune DataSourceCardView elements limit
* DataSourceCardView: Refresh button
* Wrong UTF8 handling in UUID generation library
* Proper exceptions in minified code
* File Browser: show URL for the currently selected file
* Swagger: Support request content type
* Home folder shows path
* When adding new file connector -- wrong connector opens
* Entity search: "matches" operator
* Entity search: support for units
* Files: Errors in browser console after clicking on the "texts" folder from the "demo" connection
* JS API: Examples from the "misc/audit" directory do not work
* File Browser: double-clicking on a folder does not open it in certain cases
* File Browser: Index search shows deleted files
* File Browser: filter examples
* File Browser: Confused AuthorMixin.createdOn and created
* File Browser: double-clicking on a .doc file throws an error
* Source maps don't work when URL Hash is present in script URLs
* File Browser: error in the "Sharing" panel
* JS Examples: "predictive-model.js" samples does not work
* JS API samples: viewers.js: Exception is thrown after closing the view of the script result
* File: show file name in the tooltip
* Package PP: make URL clickable
* Predictive Models: Some models do not train on Caret server (demog, race(weight, height, age)
* Add the "name" attribute for span with the table name (for the first, if there are several)
* JS: Event "d4-current-view-changed" is fired twice the first time a new view is changed
* JS Examples: "link-tables.js" samples does not work (tables are not linked)
* JS Examples: Exception after running the "events.js" sample
* Markup Viewer: If error was made in the text, then further editing does not apply
* Selenium: Add “name” attributes to file export items
* Ability to override semantic type detection
* TableView: show 'Chem' top menu when table contains molecules
* Allow to modify calculated columns
* TableView: add "Tooltip..." and "Reset filter" menu items to the new UI
* SPGi App: Exception is thrown after opening project
* Scripting: File metadata extractors
* String->int converter doesn't convert negative values
* Project upload doesn't upload table data
* Modeling: H2O does not work on HTTPS
* Data Sync: Ambiguous entity name error
* Plugin development environment
* Tika integration
* Demo: Cell imaging segmentation
* JS API samples: manipulating views
* JS API: docking
* Property panel: ability to inspect JS objects
* Layouts: matching by column id if it is present
* Tables view added from DB Table | Build Query is duplicated in the "Open" section under "Databases"
* Demo: Image classification
* Files: Exception in the panels "EXIF" and "Tika" for some pictures
* MVA: Layout is broken


# 2020.05.06 Stable Version

* HTTP ports changed. Datagrok is now listening on 8080 port and CVM port is 8090

## Latest docker images
* Datagrok (new): `766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok:1.0.75-e7d5cbd` [download](https://dev.datagrok.ai/docker_images/datagrok-1.0.75-e7d5cbd.tar)
* CVM (new): `766822877060.dkr.ecr.us-east-2.amazonaws.com/grok_cvm:1.0.75-e7d5cbd` [download](https://dev.datagrok.ai/docker_images/grok_cvm-1.0.75-e7d5cbd.tar)

## Addressed issues

* Navigation: a second Files View is created when you browse files, open a file, and press Back
* Simplified UX - Keep views under corresponding list views
* Home View 
* Double-click on the "Open" icon on a toolbox to open a local file
* GrokCvm: Reduce image size
* Markup rendering: override links to /help to open help pages within the platform
* TileViewer: Remove does not work on macOS
* Harmonized UI: Table View: "save as" combo popup
* Markup: ability to omit ".md" when referencing help pages
* Viewers forget name after deserialization
* Viewers forget name after deserialization
* GrokApi: View Menu
* OpenAPI: Hardcoded "Content-Type" for queries without body
* Notebooks: Add "Apply to" as function
* Notebook: Renaming does not influence to Notebook name resolving
* GrokApi: DataFrame "currentCol" and "currentCell"
* ScriptingViewer: Viewer continues handling events after dataframe was closed
* Editing SQL queries: tags are not recognized
* ScriptingViewer: Viewer continues handling events after dataframe was closed
* Query-driven info panels do not show up
* Files: API for uploading, renaming, moving, and deleting
* File Shares: abilty to delete, move, rename files, and create new folders - delete, move and rename
* File Browser: ability to sort by name, size, author, created, modified
* Events for client-originated file operations
* File browser: move files by drag-and-drop
* CSS support for drag-and-drop
* Support for drag-and-dropping multiple objects
* GrokApi: Reset layout for TableView
* Notebook: Download as PDF, HTML buttons
* JS API: use Rx.js for eventing
* FilesClient: Improve client-side API for working with file shares
* Docker: Datagrok and CVM without sudo
* UI: Toolbar button misplaced
* UI: Emply views list
* Files View: Folders disappear after dragging more then 1 element
* Check DB connection before creating data folder
* UI: context menus
* Moved ClinBrowser methods to the ClinBrowserPlugin class.
* Side Bar: Scrollbar does not appear in "Windows" section if many views are open
* Deprecate OpenCPU
* Bind query, script and notebook to View.entity in corresponding views
* Files View
* Clicking on profile on sidebar opens two profile views
* UI: two separators on the toolbar
* File Browser: Connections in tree are not updated after deletion or renaming
* File preview displayed incorrectly on PP
* Inspector: integrate event viewer with the property panel
* Inspector: integrate event viewer with the JS Editor (create code for intercepting and inspecting the event)
* DB deployment: Error when deploying swaggers on dev server
* User personal folder
* New UI: "favorites" pane
* Drag-and-drop: ability to custom-style drop zones with CSS
* Favorites pane: context menu for removing individual items and for removing all items
* Selenium: Add attributes "name" to controls in the "Windows" section
* Unable to create new S3 bucket (UI defaults to a regular file share)
* Data | Categorize is not available in new UI
* UI: Do not collapse toggle buttons list if there are less then 4 items or toolbox is able to fit
* UI: Convert menu button to pane switch
* A drop-down for quickly switching between views
* File Browser: Files are not displayed in the browser 
* New UI: renamed Home to Datagrok, and moved it under Help.
* Display "exception" icon in the new UI
* Exception not displayed as element in UI
* Set proper tab active on sidebar on load
* User Activity harmonization
* Horizontal scroll appears on the toolbox
* Layout gallery causes an exception
* Error when retrieving one row from Oracle with empty string value
* Teradata: "Orders" query does not work
* JS Examples: stock-prices.js does not work
* Web connections: Adding a query is available from UI
* JS Examples: stock-prices.js does not work
* Scripting: Hash Blob input and use as cache ID to avoid resending data to CVM
* Source Tree: Empty context menu for "Tables" and DB columns
* Swagger: Add manual filling of "host" and "schemes" sections for v2.0 as were added for 3.0
* UI: Workspace
* GrokCvm: Reduce image size
* Change type string -> double select by default "extract numbers"
* Saved Visual Query: Name of the resulting table should match the query name
* Forbid to change file connection root
* File indexer doesn't work with windows shares
* Files View makes few list requests after receiving an exception
* Notebooks: Optimize PDF conversion
* Optimize Windows-shares connector
* New UI: "ML" top-level menu
* Refactor Accordion and TabControl, introduce common subclasses, and expose both to JS API
* Add 'Favorites' drop down to the Property Panel header
* Multivariate Analysis does not work
* Missing Values Imputation does not work
* OpenAPI: String list parameter escaping issue
* Prevent browser caching for html files
* "Time series" test dataset
* EntityView: view name is empty
* Bar chart: Alt + Clicking on an empty space throws an exception
* Bar Chart: Starting alt-zooming on a bar selects it
* Bar Chart: mouse wheel for vertical scroll
* Modelling: Chemprop does not work
* ChemSpace integration
* Scripts: "Test" script on R: Exception after execution
* Scripts: R samples: "IIR Filter" does not work
* S3 file shares: the first letter of the file name disappears in some cases
* Scripts: R samples: "LDA" does not work
* Scripts: R samples: "Contour Plot" does not work
* Notebooks: New Notebook throws an exception
* Semantic Type Detectors: latitude/longitude
* Optimize upload utility
* Help: Internal links treated as external
* Split grok_api.js
* Fix invited users namespace permission
* Share to unregistered user is broken
* Move CVM to 8090 port.

# 2020.04.07 Stable Version


## Latest docker images
* Datagrok (new): `766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok-1.0.73-9f58f01` [download](https://dev.datagrok.ai/docker_images/datagrok-1.0.73-9f58f01.tar)
* CVM (old): `766822877060.dkr.ecr.us-east-2.amazonaws.com/grok_cvm-1.0.72-23064a3` [download](https://dev.datagrok.ai/docker_images/grok_cvm-1.0.72-23064a3.tar)

## Addressed issues

* Simplified UX
* Entities parameters 
* OpenTextView
* Files View Url to connection not handled
* DataSourcesView
* Ability to turn off new UI for selenium mode
* Added 'JS' view under Help (new location to be determined)
* WebConnector: Query with not required path param fails
* SSO Support for Community
* Web connector: Result table of the query is called "result" instead of the query name
* Handle Javascript output parameters
* Correct dapi.root resolution
* GrokApi: Add column and columns inputs
* Phaedra demo
* Unable to publish package
* Ability to override default drag-and-drop behavior for files
* Deploy fails
* Files: uploading files to shares
* FilesClient: Improve client-side API for working with file shares
* Files View should automatically open last visited folder
* Add 'Web Services' view
* Files View: 'new file share' command
* Upload: If data sync is enabled, then re-saving project with synchronization turned off does not turn it off
* SetValues works incorrectly with optional parameters


# 2020.03.30 Stable Version


## Latest docker images
* Datagrok (new): `766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok-1.0.72-23064a3` [download](https://dev.datagrok.ai/docker_images/datagrok-1.0.72-23064a3.tar)
* CVM (new): `766822877060.dkr.ecr.us-east-2.amazonaws.com/grok_cvm-1.0.72-23064a3` [download](https://dev.datagrok.ai/docker_images/grok_cvm-1.0.72-23064a3.tar)

## Addressed issues

* Auto GrokMode 
* Fixed long web api response issue
* PP: Actions tab: added both ParamCommands and Funcs to suggestions. Harmonizing handling of special function types - WIP
* ObjectView: ability to select all / none / invert
* Parsing ~500MB smiles file failed
* Files View, routing, file editors
* FileInfo property panel: Preview pane
* WebQuery: Avoid opening "variables" dialog
* Web connections: Extra scroll bar in the dialog, if there are no credentials
* GrokConnect: Neo4j support
* ORM: Ability to save objects as array
* JDBC providers: Exception after opening "Tables" tab in the tree
* MS SQL | Visual query does not work
* "Open file" file scripting feature
* Public version: AddressToCoordinates does not work
* Scripting: Custom meta parameters do not work
* Help: "pca.md" not found when open Tools | Data Science | Principal Component Analysis...
* Admin | View Layouts does not open
* Class 'ColumnInfo' has no instance getter 'tagsKeys'. Receiver: Instance of 'ColumnInfo' error in suggestions section
* Help: "query-buider.md" not found when open "Build Query..."
* Context help not displayed after opening "Visual query ..."
* Athena | northwind | Orders deployed to dev-setup with error in query text
* Help: Corresponding page does not appear when you click on the application in Admin | Apps
* Chem: Wrong R groups order and column names
* File editors: support for pictures
* File Editors: preview not working for files with spaces in the name
* File Editors: support built-in file handlers
* Help: Links to help pages on datagrok.ai/cheminformatics do not work
* Upload: Table synchronization does not work with web queries
* Unable to add View to a project
* File preview: zero height
* Improve TagsMixin.tags performance
* TreeMap: add context menu
* TreeMap: size-coding
* TreeMap: improve cell text rendering (multiline text, clipping)
* Toolbox: Viewers pane: '...' icon to show the rest of available viewers
* Info Panels: After updating browser page unblocked panels are still blocked
* Connectors: Teradata: table browsing doesn't work
* Unable to save repository
* Do not show multiple "Client is outdated message"
* Show error on login form when server is unavailable
* Selenium mode: Add name to "Unpivot" dialog selectors
* FilesView: Show directory size from index
* GrokConnect: Neo4j does not work for http and bolt+routing protocols
* Add New Column: replace +-Infinity with None
* Swagger: Import without servers
* Phaedra demo
* WebQuery from URL does not work
* Ability to select port with Docker
* Db Exception on signup: insert or update on table "projects" violates foreign key constraint "projects_author_id_fkey"

# 2020.03.13 Stable version

## Latest docker images
* Datagrok (new): `766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok:1.0.71-d0b9043` [download](https://dev.datagrok.ai/docker_images/datagrok-1.0.71-d0b9043.tar)
* CVM (old): `766822877060.dkr.ecr.us-east-2.amazonaws.com/grok_cvm:1.0.69-83c6f08` [download](https://dev.datagrok.ai/docker_images/grok_cvm-1.0.69-83c6f08.tar)

## Addressed issues

* UI: Balloon.error: defaulted autoHide to false
* DB: Database optimization
* QA: Extended query cache test
* GrokConnect: Fixed windows share mounting
* Data Query: If you edit the query, and run it without saving, the old query is called
* SPGi: Open Main tab by default
* Save as Excel: File does not download if there is an open tab with query view
* Import Excel file larger than 10MB: Balloon with error disappears with time
* Pedometer: Line chart is not added at restart
* WebProvider: Add "credentialsTemplate"
* Sparql: Add "credentialsTemplate"
* Help: Add to function auto-help parameter options information
* Ability to use both indexed and online file browsing
* Modeling: Column errors/warning should be updated after change
* Visual Query: quoted string not properly terminated
* Data Query: If all parameters are patterns, "Run" from PP opens dialog
* Viewers | Embed does not work
* GrokCvm: "cvm2" crashes sometimes by no free ram memory
* Table Algorithms Section Improvement
* Improve pattern validator
* GrokConnect: Support SSL for all popular providers
* Parameters and runtime are not updated in "General" tab for the table from the query after changing parameters and refreshing
* JS Examples view doesn't work
* Unable to set global permissions
* Postgres: List<String> parameters binding, uuid explicit conversion
* Chem | To InChi (To InChi Key): Dialog opens when called for column
* Chem | Info Panels: Panel hiding does not work
* Settings View
* Replace dataFrame should keep xp.tableInfos in consistent state
* GrokCompute: Too short proxy timeout
* Settings View improvements (APPLY button, notifications, etc)
* Execute subqueries in parallel
* Add date and git path to app and package meta
* Add author to package.json
* GrokCompute: Failed environment setup is not handled
* Visual Queries do not run
* Query Table: If you change parameters and then click on table tab, then "General" displays the old parameter values
* Domino integration demo
* Failed to get project by ID url
* Ability to specify samba version
* Algorithms section records duplicated
* Balloon: copy message to clipboard icon
* Data Jobs view: Wrong style for function list
* UI: depending on the current view, hide non-relevant top-level menus
* ORM: Rename dirName to dir
* UI: Files View
* Web Connector: Credentials were not distributed after dev version client and DB deploying
* EnamineStore integration demo
* GrokConnect: Add custom connection string feature
* FileInfo property panel: Preview pane
* Help: Redesign the wiki structure
* Ability to specify domain in SMB shares
* Teradata provider: support for schema
* Molecules are not rendered in case of view table update
* Improved semantic type detection
* Missing Values Imputation: added help reference
* Modified "files" table "size" column to bigint
* Renamed "AmazonS3Cloud" to "S3"


# 2020.02.25 Stable version

## Latest docker images
* Datagrok (new): `766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok:1.0.69-83c6f08` [download](https://dev.datagrok.ai/docker_images/datagrok-1.0.69-83c6f08.tar)
* CVM (new): `766822877060.dkr.ecr.us-east-2.amazonaws.com/grok_cvm:1.0.69-83c6f08` [download](https://dev.datagrok.ai/docker_images/grok_cvm-1.0.69-83c6f08.tar)

## Addressed issues
* /p/project.smiles/smiles causes "Ambiguous entity" error
* Ability to disable welcome email
* Ability to insert rows (Edit | Add Rows)
* Ability to store FuncCall results in local instance storage
* Ability to turn off security for current session
* Ability to use both indexed and online file browsing
* About section
* Add some new default python libraries
* Added Developers edit permissions on Datagrok-public repository
* Added options to project
* After re-publishing package with swaggers, the connections in the web connector are duplicated
* Amazon connector doesn't work
* Amazon S3 file browsing doesn't work
* Annotate db query results with the corresponding db, schema, table, and column names
* Autostart for functions
* BoxPlot: p-value setting is not serialized to layout
* Cache DataQueries results
* Chem | Solubility Prediction panel does not work
* Data | Link Tables
* Data connections: Sparql connections testing does not work correctly
* Data Query Editor: "+" icon to add results to workspace
* Data Query: "Error" in the nested parameter (The method 'contains' was called on null)
* Data Query: Queries created as Visual Query are displayed without an icon and popup with actions (down arrow) on PP
* DataQuery: DateTime pattern add "anytime" case
* DataQuery: If all parameters are patterns, "Run" from PP opens dialog
* DataQuery: Parameter values casting is missing on server and java
* DataQuey: Pattern Matcher "IN"
* Date-related functions: weekday, today, year, month, day, weekday
* DB Driver: Transactions "race condition" error
* DB Schema: Flag to cache schema
* DbTable: "Get top 100" context command
* Demo packages deploy locally is broken
* Demo repository not fully published
* DemoData: Clin data deploy as DB
* DemoData: Improve "clin" DB import
* Disabled welcome email
* Do not create index jobs for repositories
* Entity PP: "Sharing" panel
* Entity sharing doesn't work
* Executing a query from PP just pops up a dialog
* File connector: "Download ZIP" does not work for files and folders
* File connector: Files do not open by double click on Source Tree
* Files connector: Connection creation dialog does not close after clicking OK
* Func.apply, applySync: support void functions as well
* Func.setValues: performance improvements, code cleanup
* FuncCall.toUrl / fromUrl - done for existing functions; need to find a way to support functions that are not yet registered
* Funcs.byName now accepts namespace-qualified names
* Get All context command for multiple db columns
* Grid: ability to define custom cell editors
* Grid: ability to edit chemical structures by double-clicking on it
* Grok Connect: Db Table PP: Content: error
* grok.query(query, params): make params optional
* GrokConnect: Athena parametrized queries does not work
* GrokConnect: Denodo integration
* GrokConnect: MSSQL collision on usage @<> syntax with query parameters
* GrokConnect: MSSQL provider - "the value is not set for the parameter #2"
* GrokConnect: Order schemas, tables and columns in tree
* GrokConnect: Postgres SSL support
* GrokConnect: Schema table content errors
* GrokConnect: Test program
* Group property panel: move "belongs to" to "details", show administrators
* Hide chat groups in credentials owner combo
* Improved "Debugging packages" message
* Improved error message formatting when publishing packages
* Incorrect audit records
* InputBase: Validators called four times for every check
* JS API: Parametrized adHoc query does not work
* JS API: Query add AdHoc and polling interval
* JSON in which there are columns "n" and "N" does not open
* Layouts: after saving layout, immediately add it on top
* Link FuncCall to server instance
* Log functions timestamps
* Models Browser | Card view: If you select several elements, and then one, then the border around the first selected ones does not disappear (also for Notebooks)
* More date extraction functions: "Align to minute / hour / day / week / month / year"
* Move swaggers to a public repository
* Notebooks: Rewrite by second process issue is back
* Notebooks: Warnings in notebook "Chemical Space Using tSNE" after applying to the table
* OpenAPI: Support 3.0 version
* Optimized DB queries
* Packages don't deploy from repository if there is a package deployed manually
* Packages security
* Pattern inputs: invalidate an input if the entered value is not a valid pattern
* Pipe table file directly to response from local storage
* Predictive Modeling: Chemprop train fails
* Predictive Modeling: Train: Table field should be before the column selection fields (H2O)
* Predictive Modelling: Suggest models only that match by column name and type
* Project PP: Add "?" button with tooltip and help
* Project property panel: "You are the owner" string instead of own namespace
* Project property panel: add "Via" to project
* Project property panel: add icon to group
* Project property panel: Remove context actions buttonGROK-6192: Source Tree: Schemas have same "name" attribute as connection (for PostgreSQL)
* Project property panel: Share button
* Projects: After opening project, an additional view is created, even if one was already in the project
* ProjectUploadView: "Sync" switch instead of icon
* Queries get duplicated when saved
* Query / View integration
* Query Parameters: Wrong style of array parameters in web api queries
* Query: Improve parameter editor
* Query-based info panels
* Redesign project upload view
* Rename disableLogging: true to saveLog: false
* Renamed AmazonS3Cloud to S3
* Replace news.html with about.html
* Repository publishing error
* Repository publishing kills other packages
* Routing: Project url like "/p/namespace:name/viewname"
* Saving package drops connections information
* Scripting fails with https request with self-signed certificate
* Select Missing Values: Dialog has fixed size
* Semantic type detector: text - added a very simplistic detector based on a number of words and max word length
* Set SSL attribute in datagrok connection on deploy
* Setting view.path should modify currently shown URL
* Store users_data_storage in DB
* Support schemas for PostgreSql provider
* Swagger for our internal elastic search
* Table refreshing doesn't work
* Table serialization: Stream bytes instead of sending blob from memory
* TableMeta: Add more query parameters
* TableQuery: limit is not always propagated to the server
* TableView: Changing table does not change table in platform: added VIEW_TABLE_REPLACED event, adjusting xp.tables but not xp.tableInfos yet
* Treat Amazon secret key as a password
* TreeView: Item context menu exception
* Unpivot / Merge - dialog
* upload.dart fails on second time on package that contains query linked to connection outside the package
* upload.dart fails with self-signed SSL certificates
* Upload: Project picture: Closed views shows as empty in the list
* URL Sharing from JS API / Packages
* URL-encoded queries: error when parsing datetime
* User can't see package if other user has debug version
* UsersDataStorage can't read value with currentUser: false
* Visual Query: Filters
* Web connector: Balloon-error text always with API key for all types of security
* Web: Add connection causes error
* When uploading package - check if package with same name was uploaded and permissions were not granted


# 2020.01.15 Stable version

## Latest docker images
* Datagrok (new): `7766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok:1.0.65-8dc0463`
* CVM (old): `766822877060.dkr.ecr.us-east-2.amazonaws.com/grok_cvm:1.0.31-6261445`

## Addressed issues
* Scripting/Notebooks: Ability to create separate environment for custom packages
* AmazonCloudProvider: FileExists raises signing error
* Grok API: Ability to set Grok variables
* DbContext: connections leaking fixed
* Credentials Management System
* JS API: Grid: cell-based tooltip
* Pressing Ctrl+Z in a text field should not initiate platform's undo command
* Grok API: Rename gr to grok
* Grid: "Ctlr+A" should only select filtered rows
* Grid: Column filter: Stops working after it's used once
* Scatter plot: If you zoom with the touchbar (two finger pinch), then the browser window also zooms
* Hotkeys: Ctrl+Page Up (Page Down) does not go to the first (last) row in Grid
* Grid: Shift+PgUp/PgDown should select rows
* Grid: Ctrl+Shift + ↑↓ does not deselect rows
* Filter Rows: returns no rows if used on filtered values
* GrokConnect: Virtuoso support
* GrokConnect: ORACLE num with int precision should be handled as integer data type
* GrokConnect: Schema does not work for ORACLE
* Oracle connector: error when returned dataset contains columns of the same name
* GrokConnect: Parametrize Query timeout
* Connections list: Show only used datasources on first page
* UI: Compare client and server versions on startup
* Dialogs: Columns Selectors: Ability to mark all columns found after searching
* Grok API: Moved NGL viewer to common libraries
* Git data provider
* HashId: ability to overwrite ID
* Credentials: ability to save bool and int parameters
* Mask passwords in settings
* LDAP authentication
* GrokConnect: Oracle timestamp and clob column types cause errors
* Incorrect datetime parsing (12/01/19 01:00:00 gets imported as Jan 12th)
* LDAP authentication. Invitations support
* Missing name field in connection editor
* Projects: Set default project picture
* GrokConnect: Add Schemas and Views
* GrokConnect: Connection test doesn't work
* GrokConnect: Oracle error on datetime pattern
* UI: "Sharing" entity panel
* Data connectors: Oracle: an issue with table browsing
* Connection Tree: Tables do not load in some providers (MS SQL, MySQL, MariaDB, PostgreSQL)
* Swagger: Unable to import FastAPI
* Scripting: Replace column should join column, if it does not exists
* Data connections: "Test" does not work when creating new connection
* Query View: Subquery handled as result
* Query parameters don't fit in property panel
* GrokConnect: Add Schemas and Views
* GrokConnect: Datetime pattern "anytime"
* Query: Pattern parameters should not required by default
* Show welcome string on login form
* Show logo on the login form
* Share entity: show previous shares
* CSV import: FormatException: Invalid integer
* Repository publishing error
* GrokConnect: Impala support
* GrokConnect: Empty columns except last row in case of unknown column data type
* GrokConnect: Prevent of data losses in case of unknown column datatype
* GrokConnect: Improvements: Float, Bool, DateTime and Int columns
* Exclude system schemas in Oracle provider


# 2020.01.07: Service Update

## Latest docker images
* Datagrok (new): `766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok:1.0.X-0dc9ba6`
* CVM (old): `766822877060.dkr.ecr.us-east-2.amazonaws.com/grok_cvm:1.0.31-6261445`

## Addressed issues
* Scripting/Notebooks: Ability to create separate environment for custom packages
* AmazonCloudProvider: FileExists raises signing error
* Grok API: Ability to set Grok variables
* DbContext: connections leaking fixed
* Credentials Management System
* JS API: Grid: cell-based tooltip
* Pressing Ctrl+Z in a text field should not initiate platform's undo command
* Grok API: Rename gr to grok
* Grid: "Ctlr+A" should only select filtered rows
* Grid: Column filter: Stops working after it's used once
* Scatter plot: If you zoom with the touchbar (two finger pinch), then the browser window also zooms
* Hotkeys: Ctrl+Page Up (Page Down) does not go to the first (last) row in Grid
* Grid: Shift+PgUp/PgDown should select rows
* Grid: Ctrl+Shift + ↑↓ does not deselect rows
* Filter Rows: returns no rows if used on filtered values
* GrokConnect: Virtuoso support
* GrokConnect: ORACLE num with int precision should be handled as integer data type
* GrokConnect: Schema does not work for ORACLE
* Oracle connector: error when returned dataset contains columns of the same name
* GrokConnect: Parametrize Query timeout
* Connections list: Show only used datasources on first page
* UI: Compare client and server versions on startup
* Dialogs: Columns Selectors: Ability to mark all columns found after searching
* Grok API: Moved NGL viewer to common libraries
* Git data provider
* HashId: ability to overwrite ID
* Credentials: ability to save bool and int parameters
* Mask passwords in settings
* LDAP authentication
* GrokConnect: Oracle timestamp and clob column types cause errors

# 2019.12.23: Stable version
## Latest docker images
* Datagrok (new): `766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok:1.0.48-5791556`
* CVM (new): `766822877060.dkr.ecr.us-east-2.amazonaws.com/grok_cvm:1.0.31-6261445`

## Addressed issues
* Data connections: Error when connecting to Oracle database
* Scripting/Notebooks: Ability to create separate environment for custom packages
* Show placeholder on login form while server is deploying
* Scripts Samples: Few bugs fixed
* Support anonymous SMTP server
* Layout issues fixed
* Grok JS API: Added TableView.dataFrame setter
* Grok JS API: Fixed parametrized queries
* Grok JS API: init() now returns Promise
* Grok Connect: Fix javax.servlet old version issue
* SSL support for Postgres provider
* Grid: problems with touchpad-initiated scrolling
* Grid: custom cell editors
* Amazon S3 Adapter: Fixed list() response, AES256 support
* Packages: Made GrokPackage.init method asynchronous
* Packages: An example of using RDKit (compiled to WebAssembly) on the client side
* Postgres provider: bugs fixed
* Notebooks: bugs fixed
* Package upload bug
* Table Manager: After renaming the table, its name is not updated in the Table Manager (Alt+T)
* Toolbox: After renaming the table, its name is not updated on the Toolbox
* Copying a cell when block is selected leads to opening developer's console in the browser (Ctrl+Shift+C)
* View | Tooltip: 'all on' doesn't set all on after the first click


# 2019.12.18: Service Update

Addresses a number of issues identified during the technical evaluation. 

## Latest docker images
* Datagrok (new): `766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok:1.0.42-7866749`
* CVM (old): `766822877060.dkr.ecr.us-east-2.amazonaws.com/grok_cvm:1.0.25-54e90ea`

## Addressed issues
* Grid: problems with touchpad-initiated scrolling
* Grid: custom cell editors
* Amazon S3 Adapter: Fixed list() response, AES256 support
* Packages: Made GrokPackage.init method asynchronous
* Packages: An example of using RDKit (compiled to WebAssembly) on the client side
* Postgres provider: bugs fixed
* Notebooks: bugs fixed

# 2019.12.05: Service Update

Addresses a number of issues identified during the technical evaluation.

## Latest docker images
* Datagrok (new): `766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok:1.0.33-629a22f`
* CVM  (old): `766822877060.dkr.ecr.us-east-2.amazonaws.com/grok_cvm:1.0.25-54e90ea`

## Addressed issues
* Package upload bug
* Table Manager: After renaming the table, its name is not updated in the Table Manager (Alt+T)
* Toolbox: After renaming the table, its name is not updated on the Toolbox
* Copying a cell when block is selected leads to opening developer's console in the browser (Ctrl+Shift+C)
* View | Tooltip: 'all on' doesn't set all on after the first click