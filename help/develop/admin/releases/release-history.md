---
title: "Release History"
---

## Latest version

| Service                                                  | Docker Image                                                                                      |
|----------------------------------------------------------|---------------------------------------------------------------------------------------------------|
| [Datagrok](../infrastructure.md#datagrok-components)     | [datagrok/datagrok:1.13.8](https://hub.docker.com/r/datagrok/datagrok)                            |
| [Grok Connect](../../../access/data-connection.md)       | [datagrok/grok_connect:1.3.2](https://hub.docker.com/r/datagrok/grok_connect)                     |
| Grok Spawner                                             | [datagrok/grok_spawner:1.3.2](https://hub.docker.com/r/datagrok/grok_spawner)                     |
| [Grok Compute](../infrastructure.md#grok-compute)        | [datagrok/grok_compute:1.4.1](https://hub.docker.com/r/datagrok/grok_compute)                     |
| [Jupyter Kernel Gateway](../../../compute/scripting.md)  | [datagrok/jupyter_kernel_gateway:1.5.2](https://hub.docker.com/r/datagrok/jupyter_kernel_gateway) |
| [Jupyter Notebook](../../../compute/jupyter-notebook.md) | [datagrok/jupyter_notebook:1.1.0](https://hub.docker.com/r/datagrok/jupyter_notebook)             |
| [H2O](../infrastructure.md#h2o)                          | [datagrok/h2o:1.1.0](https://hub.docker.com/r/datagrok/h2o)                                       |
| [CVM Nginx](../infrastructure.md#load-balancer)          | [datagrok/cvm_nginx:1.10.0](https://hub.docker.com/r/datagrok/cvm_nginx)                          |

See also:
- [Versioning policy](versioning-policy.md)
- [Docker-Compose](../docker-compose.md)


## 2023-03-22 1.13.8

### Addressed Issues

* (Enhancement) [#1552](https://github.com/datagrok-ai/public/issues/1564713406): Scatter Plot: line labels enhancements 
* (Bug) [#1613](https://github.com/datagrok-ai/public/issues/1598588170): 'Error loading scatter plot' for apecific data on applying layout 
* (Bug) GROK-12382: Websocket reconnects after it closed (WIP)
* Datagrok Docker: Allow grok user use sudo without password 
* (Enhancement) [#1111](https://github.com/datagrok-ai/public/issues/1406443198): Viewers: PivotGrid: transformation result preview 
* histogram imitates filter look 
* (Bug) [#1628](https://github.com/datagrok-ai/public/issues/1611012446): Wrong data is shown on line chart 
* (Bug) GROK-12479: Project deserialization leads to 50MB response size 
* (Bug) GROK-12408: Chem: 'Chem | Possibly a malformed molString at init: `undefined`' error when hovering over mol column's header 
* Help: New release history 
* fix for column tooltip 
* Git: Fixed conflicts 
* fix scatterplot 
* (Bug) [#1651](https://github.com/datagrok-ai/public/issues/1623650870): Line chart: incorrect legend for Split by when Row Source=selected 
* (Bug) [#1650](https://github.com/datagrok-ai/public/issues/1623611567): Line chart: Unsupported operation: NaN.floor() error 
* (Bug) GROK-12498: Chem: grid scrolls to the right when you hover over a molecular column 
* connected not supported 
* updated version of grok_connect-1.2.1.jar 
* added snowflake dependency 
* updated snowflake provider 
* made fixes 
* updated version 
* made smal fix 
* GROK-12104: Schema browsing support 
* added equals and hashcode methods 
* added tests dependency, updated psql version 
* created sql scripts for testing psql provider 
* created base class for testing 
* created test class for postgres provider 
* created tests infrastructure 
* made field params public 
* added access modifiers 
* made slight changes 
* made small changes 
* removed equals&hashcode; added getter for data 
* updated sql scripts 
* added comparator for data frames 
* (Improvement) GROK-12090: GrokConect: Postgres Harmonization 
* GROK-12188: Add regex support in string pattern for Postgres Provider 
* added method order annotation 
* added clear display of parametrizied test purpose when running tests 
* GROK-12100: GrokConnect: Full types support for Postgres Provider 
* added oracle dependency to pom 
* Grok Connect: Bump version to 1.3.3 
* GROK-12089: Grok Connect Harmonization (WIP)
* Update grok_connect.yaml 
* Grok Connect: Updated build script 
* GiHub Actions: Push Grok Connect from release branches 
* Grok Connect: Remove WIP tests 
* created test scripts for oracle provider 
* GROK12089 
* GROK12210: Edited test queries 
* GROK12211: intermediate commit 
* GROK-12211: Create test class and test cases for Oracle Provider 
* GROK-12212: GrokConnect: Full types support for Oracle Provider 
* GROK-12210: GrokConnect: Test Queries for Oracle provider 
* (Improvement) GROK-12091: GrokConect: Oracle Harmonization 
* (Improvement) GROK-12103: GrokConect: Snowflake Harmonization (WIP)
* GROK-12286: GrokConnect: Test Queries for Snowflake provider 
* (Improvement) GROK-12346: Add the ability to display sql views along with tables 
* GROK-12344: GrokConnect: Test Queries for MSSQL provider 
* GROK-12352: Create test class and test cases for MSSQL Provider 
* (Improvement) GROK-12092: GrokConect: MSSQL Harmonization 
* GROK-12372: GrokConnect: Test Queries for Athena provider 
* (Improvement) GROK-12095: GrokConect: Athena Harmonization (WIP)
* GROK-12373: Create test class and test cases for Athena Provider 
* Grok Connect: Bump version 
* (Bug) [#1637](https://github.com/datagrok-ai/public/issues/1616819912): Filters re-ordering dialog is not opened if a custom filter was added to filters panel 
* (Bug) [#1724](https://github.com/datagrok-ai/public/issues/1632100036): Table is not rendered properly after switching dataframe if it was sorted 


## 2023-03-14 1.13.7

### Addressed Issues

* fixes #1609 FilterGroup.add issue


## 2023-02-22 1.13.3

### Addressed Issues

* (Bug) Connecting to local file storage results in error (WIP)


## 2023-02-21 1.13.2

### Addressed Issues

* 1583 manual sketcher resizing (#1598)
* closes #1599 Colouring and format changes are reset unexpectedly on interacting with filters if group tooltip was set
* line chart broken

## 2023-02-15 1.13.1

### Addressed Issues

* Standard color schemes : update "Traffic Lights" palette
* Grid: Enable getCellEditor for meta.columnTags
* Grid: Fix for getGridCell
* fix grid demo dir in its package.json

## 2023-02-13 1.13.0

### Addressed Issues

* update GrokConnect to version 1.1.0
* Charts: syntax fixes in tree-viewer.ts
* update GrokConnect version to 1.1.0
* Charts: fixed filter problem in sankey and chord
* Add null values support to ui.dateInput
* Core: add radioGroup identification to menu options where relevant
* Test Manager: Sort categories alphabetically
* Charts: fixed tree-utils mapRowsToObjects() method
* Charts: removed data duplication in chords and sankey
* bio lib: PDB interfaces
* BiostructureViewer: NGL typings, NglViewer open file
* Dendrogram: Package icon, minor fixes
* Harmonize Demo package
* NglViewer: Move NGL to BiostructureViewer
* Utils: history panel engine refactored
* Utils: removed Cloning step from ComputationView
* Utils: added scriptsCache for HistoryPanel
* Utils: fixed getPackageUrls fails for old runs
* Utils: minor version bump
* datagrok-tools: fix exit codes in grok test
* Fix empty state icon style
* Viewers: added sankey viewer redrawing on filltering
* ClinicalCase: Fix connection to S3 bucket with files
* Test manager: Mark all tests yellow if Unhandled exception occurred
* #863: moved pepsea to Bio package
* (Bug) Table content queries fire events
* Removed StackedBarChart
* Chem: menu: sketcher types as radio buttons
* Viewers: fixed labels which don't move with selected graph node
* #1480: Demo App: WIP
* PIDB: Connection string generation
* bio lib: Fix TreeCutOptions type for min, max optional
* Dendrogram: Fixes allowing empty trees
* Libraries/ml: activity cliffs refactoring
* Utils: added updated indicators to history panel
* PI DB icon
* Compute: deps update
* Bump grok_connect versions
* Compute: updated maintainer email
* update grok connect version to 1.2.0
* Reverted serialization version
* Charts: removed exportation of sankey and chord
* Chem: fix UI tests
* Tutorials: library refactoring
* ScaffoldTree: removing script needed to install ScaffoldGraphDG module + all its occurrences
* Tutorials: add missing property to the library
* bump datagrok tiils version to 4.7.6
* Chem: fixing top menu script based.mutate
* closes #1467 date values should be shown in tooltips.
* Tutorials: reuse code from the tutorials lib
* #1465 Improve usability of range filters
* Chem: map identifiers (increase wait)
* dpetrov/1364: Support for non-standard fragments in structure search
* dpetrov/1217: Scaffold-Tree: fixed environments, fixed aromatics, various ui-changes
* closes #1465 Improve usability of range filters
* #1472 Reuse 'Pick up style' functionality for JsViewers
* Utils: parentView made optional
* chen-meta: bump version
* dpetrov/1217: Scaffold-Tree: fixed autogenerate & deserialize
* Tree view outline
* ML: fix and bump
* Chem: fix and bump
* Chem: skip map identifiers test
* Viewers: changed JS to TS
* Viewers: added Circos for chord viewer
* Viewers: syntax fixes
* Viewers: rewrote code base from JS to TS
* Viewers: fixed word-cloud onPropertyChanged()
* Libraries/ml: fix and version up
* Chem: updated ml dependency, bump version
* Ketcher: updated utils dependence, bump version
* Charts: restructured the package
* Adjusted the type for Entity.getProperties() to {[index: string]: any}
* Issue #1492: Chem | Elemental analysis: malformed data in dataset handling
* Chem: bump new version
* Tooltip.show now accepts String | HtmlElement
* Minor help updates
* Packages Charts and Viewers merged into one package
* closes #1481 Chem | Cell with molecule: ‘Use as Flter ‘ action causes error
* JKG: Disable Datlas autostart in container
* BiostructureViewer: add 1pdbq.sdf
* #1336: creating new cluster
* #1336: custom clusters serialization
* #1336: enable cluster removing
* Charts: added documentation about new viewers
* #1465 Improve usability of range filters, layout problem
* (Bug) GrokConnect: Neptune list parameters don't work
* #1465 filters improvements default size
* (Bug) Viewers | Map: layer data export doesn't work
* Tutorials: items to fix (WIP)
* #1333: multiple views
* update grok connect version to 1.2.1
* Proper parameter visibility
* Peptides: Disabled Get Peptides Structure
* utils lib: fix errorToConsole for number type
* bio lib: Fix for handling errors in toAtomicLevel, fix for gaps in MSA in toAtomicLevel bump version
* Bio: bump version for dependencies on utils liband bio lib
* Peptides: throw error if barchart rendering fails
* #1077: UI fixes
* Peptides: fixed aggregated columns not showing in LST
* Peptides: fixed errors in LST
* Peptides: lint fixes
* Peptides: 1.6.0 release
* GIS: fixed odd adding of Markers GL Selection to layers list
* (Bug) Viewers | Map: after setting the 'Render Type', the Layers menu update is one step late
* bio lib: Fix parse newick for no name root, reduce ITreeHelper.newickToDf args, bump version
* Dendrogram: Add tests for TreeHelper, reduce TreeHelper.toNewickDf args, fix for no name root
* Dendrogram: bump version
* ScaffoldTree: script minor fixes
* Dendrogram: Fix using nodePrefix in TreeHelper.toNewickDf correction no name root node
* Utils: HistoryPanel redesigned
* Utils: added ribbonMenu clearing
* Utils: patch version bump
* MLB: Add VRs count to antigen list, fix PTM filter, checkVrForTree, VR space, fix tree for grid with cellType='html',
  fix merged clone trees on network diagram, calc distance on JUNCTION region
* Charts: fixed GIFs size in md file to 800x500
* Compute: migrated default scripts to JS
* closes #1483 Event.path removed from chrome
* Libraries/chem-meta: updated rdKit Reaction api, version up
* Chem: added chemical reaction semantic type to js-api
* Chem: added rdkit reaction cell renderer, version up
* Detectors: one test for all (WIP)
* Docs: Create script to migrate help to Docusaurus
* Chem: #1497 aromatic query substructure search
* bio lib: Fix export type IDendrogramService, TreeCutOptions
* Dendrogram: bump version for dependency bio lib new version
* Charts: replaced jQuery with cash-dom
* bump node version in action packages
* bump npm version to 9 in actions package
* Charts: moved tree-map to deprecated
* dpetrov/1497: Package Chem: Added support for non-standard fragments to structure highlighter
* Compute: renamed default scripts
* Compute: patch version bump
* Utils: moved all history actions to HistoryPanel
* (Bug) JS API: DG.toJs won't work for predictive models
* Chem: reactions renderer \- split row into reactants and products in case of long reaction string
* Charts: added filtering to radar
* ApiTests: added tests for returning FuncCall options
* DevTools: version bump
* js-api: Add semType Molecule3D, PDB_ID
* Charts: added selection to surface-plot
* bio lib: No re export type IDendrogramServer and TreeCutOptions, bump version
* Charts: fixed radar tests
* Dendrogram: Fix import IDendrogramService, bump version
* Dendrogram: Move NglViewer README, PDB_ID semType detector, PDB_ID semType widget panel to BsV, fix Molecule3D file (
  pre)viewers, fix PDB file handler
* bio lib: Remove index.ts
* (Bug) Failed to share a connection
* Fixed NPE
* ApiSamples: version bump
* bio lib: Reduce vd-regions modules, bump version
* Chem: mol2 importer, fix bug with lacking bonds
* Charts: updated substituent analysis and renamed it to group analysis
* Bio: Fix imports from bio lib, bump version
* Dendrogram: Fix imports from bio lib
* PTV: Fix imports from bio lib
* GIS: UI fixes (header in property panel)
* MLB: Selector for VR space dimensionality reduction method UMAP or t-SNE, fix imports from bio lib
* Fix splitter resizing
* GIS: added title and html table for coordinates in property panel
* (Bug) API: `ui.tableInput` shows invalid DF names if tables are opened before its' creation
* Help: DocSearch token
* Chem: ScaffoldTree: removing the highlighting when necessary
* Chem: ScaffoldTree: added the lines (should still be improved though)
* TreeViewNode: added a generics type to represent node.value
* Chem: Scaffold Tree: major code cleanup; made most of it strongly-typed as well
* Chem: Scaffold Tree: added 'ringCutoff' and 'dischargeAndDeradicalize' parameters (to both script and viewer)
* Chem: Scaffold Tree: minor code cleanup
* Chem: Scaffold Tree: saving and opening \- WIP
* Chem: Scaffold Tree: moved UI code to CSS
* Chem: Scaffold Tree: beautified the lines
* Chem: Scaffold Tree: minor fix
* Minor visual improvements
* Fixed #1502: fixed missing custom cluster name in distribution widget
* #1077: WebLogo fixes
* Chem: Scaffold Tree: support for loading and saving trees
* Chem: Scaffold Tree: ability to invoke it from the top menu
* Fixed #1503: added custom cluster name property in logo summary table
* Fixed #1504: fixed aggregated column not appearing in LST
* Fixed #1505: added aggregated columns to distribution widget
* Peptides version 1.6.1
* Chem: Scaffold Tree: Fixed the rendering of non-standard nodes fragments
* Chem: Scaffold Tree: cosmetic improvements
* Chem: Scaffold Tree: Returned accidently removed version of the aromatic processor
* Chem: Scaffold Tree: help
* Chem: bumped up the version
* (Bug) Scatterplot: range sliders issue when a legend is present
* Chem: Scaffold Tree: Fixed structure highliter bug on node selection
* Chem: Scaffold Tree: Fixed saving scaffold tree as part of layout
* NglViewer: Remove the package
* Chem: added tests for substructure search with aromatic bond/explicit hydrogen
* Packages: Tutorials: version bump
* Widgets: obsolete smiles widget removal
* Tutorials: clarifications for the scripting tutorial
* Charts: added filtering to group-analysis
* (Bug) FuncCall.options can't be saved to DB
* closes #1499 Legend/viewer inconsistencies when table is filtered
* #1354: fixed RGroup failing with OCL sketcher
* Dendrogram: Fix hierarchicalClustering for datetime type column as float
* JS UI Test: SPE
* closes 1464 Change order of filters
* JS API: add showNextTo parameter to Dialog.show
* Bio: fixed typo
* Charts: refactored sankey
* Issue #1417: Chem | Column Actions | Descriptors: error when empty values are present
* Chem: Scaffold Tree: Added support to cancel tree auto generation
* MLB: Fix VR space scatter plot X, Y data column names on changing data frame
* GIS: removed the showing of custom property panel of GisObject
* (Bug) Add missing DB index
* Help: Fix linter errors
* Charts: refactored chord
* Docker: load image to system
* Build: Universal build script for all OS
* Charts: refactored globe
* Build: Windows compatibility
* Revert "Build: Windows compatibility"
* Widgets: bumped up the version
* #928 recheck properties order
* Build: Docker cache for files in Windows
* Build: Docker files copy
* The ability to construct byte array column using fromList constructor
* (Bug) File Shares: Mount Windows shares
* Charts: added aggregation types to group-analysis
* Cloudfimation:add 80 port to cvm_int_lb SG
* Charts: fixed typos in radar
* udate node to 18.12.0 in puppeteer
* Charts: fixed typos in word-cloud
* Fixed #1517: fixed cluster not removing
* add webpack to DSP-package
* Fixes #1521: don't create empty new views
* #1518: replaced settings button with wrench icon
* #1518: removed view name field
* #1518: moved New view button to property panel on selection
* Added Dmytro to beta users
* Chem: add tests for mol2 importer
* Bio: substructure filters refactoring
* Alation: updated api dependency
* GIS: fixed zipcode detector
* bio lib: Move WebLogo viewer to Bio package
* Charts: version up
* (Bug) Color coding: dialog size issue when coloring a molecular column
* Lstolbov/curves (#1526)
* Charts: skipped timelines properties test, version up
* Bio: Move WebLogo viewer to Bio package
* Peptides: Remove reexporting index.ts from bio lib
* Charts: skipped timelines viwer creation tests, version up
* #1314: Chem \- fixed properties widget
* #1487: fixed Structural Alerts panel mol orientation
* Chem: updated Group analysis fucntion, bump version
* Chem: fixed typo in panels
* #1454: use the same instance of rdkit module in structural alerts panel
* Chem: bump version
* FottingTools: bum version
* Revert "GitHub Actions: Bump Meta package version for"
* Peptides: refactored views
* Build: Remove Windows symbols
* Chembl: New server for databases
* #1518: moved custom cluster actions to property panel
* Charts: changed mode to development
* #1518: Monomer-Position viewer hints in icon
* Chembl: Add Chem as grok dependency
* #1518: nullable clusters column input
* #1524: nullable clusters field
* #1524: autogenerate aligned column name
* #1524: added inputs tooltips
* (Bug) Viewers | Trellis plot: errors occur when a trellis plot is displayed in R-Groups Analysis
* (Bug) Package can't be installed with npm proxy specified
* Helm: fixing datagrok-libraries/bio version
* Close #1466 MLB: Fix WebLogo tooltip invalid count rows/sequences with specific monomer at particular position, add
  test, add WebLogo.filterSource property
* Dendrogram: Dendrogram viewer fix props category Style
* JS API: fix Dialog.showModal
* Dendrogram: documentation
* Dendrogram: documentation fix gifs
* Dendrogram: Scroll grid on wheel on tree
* MLB: Fix binding embed columns with VR space scatter plot
* (Bug) Charts | GroupAnalysisViewer: NullError: method not found: 'root' on null
* (Bug) AddNewColumn without passing column type causes the exeption (WIP)
* Viewers: Histogram: An option to render the distribution as spline when no split is defined
* dpetrov/1532: Chem: Fixed structure highlight with explicit hydrogen
* Dpetrov/1533 chem structurre search makeover (#1535)
* Grok Compute: Default configuration
* Grok Compute: downgrade numpy Fix 'AttributeError: module 'numpy' has no attribute 'bool'' error Numpy removed method
  bool in 1.24.0 More infromation: https://github.com/numpy/numpy/releases/tag/v1.24.0
* #1529: Chem \- added check if package has been initiated
* DSP: package.json indent
* Update Meta package
* 1476 chem similarity search improvements (#1536)
* Peptides: minor refactoring
* Chem: described similarity search in README
* Build: Multiple actions for docker build
* Git: ignore files generated by datlas during start
* Chem: version up
* Docs: Release documentation
* Docs: Docusaurus documentation
* Grok Shared: Default values for CVM using ipv4 localhost address
* Chem: fixed typos in readme
* Chem: updated similarity search in Wiki
* Chem: similarity search, updated Wiki
* Grok Spawner: optional registry credentials
* Build: Docker output on build action
* JKG: Octave libspqr2 library
* (Bug) Auto init method of Chem package is called multiple times
* JKG: Octave libopenblas library
* Demo DB: MSSQL healthcheck
* ApiTests: New Demo DB servers
* Fixed some tests
* Error handling
* Bump package version
* (Bug) Some inputs cannot be nullified (WIP)
* DevTools: TM fixes
* (Bug) Viewers | Scatterplot : 'zoom by filter' doesn't work, when Filter is set first
* closes #1498 Line chart coloring legend does not match the right colors in some cases
* 1465 Improving usability of range filters
* closes 1475 Chem | Filtering: Filter Panel has to be reopened to display the molecule after Current value>Use as
  filter action
* BiostructureViewer: Change ownership to atanas@datagrok.ai
* Chem: fixing linter errors in similarity search wiki
* Chem: similarity serach, linter fixes
* (Bug) Timelines viewer: properties return null dateFormat
* (Bug) Timelines viewer cannot be created (error in onFrameAttached)
* ApiTests: Viewer.props.setDefault/resetDefault
* MLB: Check data consistency VRs of the antigen in tree2 and antbody2antigen
* JS API: fix build error
* Grok Spawner: Debug in log
* Grok Spawner: Fail on docker push errors
* MLB: Check data consistency VRs of the antigen in tree2 and antibody2antigen fix
* Logger: format stacktrace for test engine
* Grok Compute: num cores based on cpu limit
* DevTools: testFunctions fix
* Add team page
* Update css styles
* #1307: fixed PubChem api panels fixed width in window
* Top panel refimnement (WIP)
* #1030: limited drug likeness panel vertically
* #1030: limited structural alerts panel vertically
* closes #1500 Line chart coloring should be able to scale to visible range of data
* dpetrov/1545: Fixed appending orphans folders after scaffold tree is loaded from a local file
* Grok Spawner: Deploy containers to datarok host if host label is missing
* BiostructureViewer: Fix imports from bio lib, enable ts strict mode, molstar viewer initial
* JS API: add windows.showContextPane
* Cheminformatics | Stability and Polishings (WIP)
* Grok Spawner: Fail if hosts with required labels are missing
* Release notes 1.12.0
* #1465 improve usability of range filters
* Fixed #1543: MPR row ordering fix
* (Bug) Predictive Modeling | H2O: Input fields are not shown in training view (on https://public.datagrok.ai/)
* trellis plot issue
* added type-ahead control to utils
* Help: Documentation for docusaurus
* Docker Compose: Additional profiles
* Docusaurus: Team page with custom images and css
* JS API: package-lock.json upgrade
* Docusaurus: unique location for images
* set julia version to 1.6.7LTS
* Docusaurus: fix
* MLB: Fix VR space zooming on current igphyml clone/tree
* Peptides: changed category to Bioinformatics
* Help: Release notes headings fix
* Chem: columns on null (fixed for To Inchi and To Inchi Keys)
* Chem: fix tests due to panels refinement
* MLB: Bump bio/Bio dependency version for updated WebLogoViewer, VdRegionsViewer
* Rename the property panel to context panel (#1554)
* Chem: fix test for similarity search
* Rename the property panel to context panel
* Chem: skipped UI tests, version up
* Chem: fixed tests, version up
* (Bug) Test engine: enable script testing mechanism for test manager and all packages
* Docs: Fix link in help documentation
* DevTools: bump version
* Usage Analysis: changed category to General
* closes 1482 Chem | Cell with molecule: ‘Sort by similarity ‘ action causes an error
* Docs: Add Browse Schema to Supported Connectors documentation
* MLB: Fix for network diagram data columns
* Meta: Remove ngl-viewer package, bump version
* Chem: #1306 chem item in panel
* Packages: adjust package categories
* Bio: fixed substructure filters tests
* Bio: Enable tests for substructure search, fix grokDependencies versions
* Helm: Fix imports from bio lib, bump version
* fix julia scripts working
* #1523: mutation cliffs now respond to filtering
* #1523: fixed distribution widget for custom clusters
* #1523: fixed cluster creation
* #1523: fixed wrong tooltip stats in LST if filter is applied
* #1557: fixed filter not applied on project load
* #1523: made logo summary table react to filtering
* #1523: store cluster names in selection
* #1557: fixed wrong object is shown on project load
* #1556: added aggregated column to cluster tooltip
* Peptides: fixed Monomer-Position hints
* #274: fixed save and load test
* Library GridExt: Adding repaint call on adding pinned columns
* Removed extra menu items
* Package PowerGrid: Bumping Up version
* #1498 reset filter legend refresh
* 1500 readjusting X axis
* grid core height by visible cols
* scatterplot invalidate method
* scatterplot invalidate canvas method
* Bio: Downgrade platform version dependency
* Add .columns property to DG.ViewLayout
* Chem: diversity search \- additional molecule properties
* added comments to TypeAhead
* added DropDown control
* Build: Moved grok connect to github
* Avoiding using dapi on server to upload tables
* Build: Grok Connect in public
* Docs: how to set the default sketcher + link to the Release Notes
* Peptides: version 1.7.0 release
* Docs: fix link in help documentation
* Peptides: autotests fix
* Docs: fix links in help documentation
* typeahead styles
* add ribbonPanel div
* update ribbonPanel
* Chem: detectors small fix
* closes 1553 slider in log mode
* closes 1551 Line chart: X categorical axis is not aligned with value
* Mdolotov/sketchers refactoring (#1349)
* (Bug) Databases: error in Sharing tab after deleting the shared query
* grid should not throw errors

## 2023-01-27 Dev build 1.12.1

### Addressed Issues

* (Bug) File Shares: Mount Windows shares (WIP)
* (Bug) FuncCall.options can't be saved to DB
* (Bug) Package can't be installed with npm proxy specified
* (Bug) AddNewColumn without passing column type causes the exeption

## 2023-01-24 Dev build 1.12.0

We've released a new version of the Datagrok platform (1.12.0). This release focuses on new features for visualization
and usability and traditionally on platform stability.
Here are some of the biggest improvements:

* New [Scaffold Tree](release-history.md/#enhancements-in-packages ) visualization that organizes molecular data sets by
  arranging molecules into a tree hierarchy based on their scaffolds. For details,
  see [Scaffold tree](../../../domains/chem/scaffold-tree.md).
* [Dendrogram](release-history.md/#enhancements-in-packages) for interactive exploration of the hierarchical clustering.
  For details, see [Dendrogram](https://community.datagrok.ai/t/dendrogram/721).
* Brand new [Map](release-history.md/#enhancements-in-packages). To learn more,
  see [Map viewer](https://community.datagrok.ai/t/visualization-related-updates/521/27)
* [Peptides](release-history.md/#enhancements-in-packages) new features and capabilities, such as custom clustering and
  multiple views, as well as a slightly redesigned user interface and improved application stability. To learn more,
  see [Macromolecules updates](https://community.datagrok.ai/t/macromolecules-updates/661/11)

### Visualization and usability improvements

* Color coding
  * Scatter plot now has a legend for continuous color coding.
  * In a grid, you can apply [color coding](../../../visualize/viewers/grid.md/#column-color-coding) to the text or
    background. This option is available for all linear, categorical, and conditional schemas.
  * When inheriting the color coding from the grid, adjustments to the min/max made in the plot are reflected in the
    grid. However, if the original column is not color-coded and selected for the color column in the scatter plot, the
    configuration isn't applied to the column in the grid.
  * Fixed changing the linear color scheme issue. Now the null values don’t get colored.
  * Fixed inconsistent behavior of color coding checkboxes in the columns context menu.
* [Scatter plot](../../../visualize/viewers/scatter-plot.md):
  * Added the exact min/max of the column on the axes' ticks.
  * Context menu: marker section doesn’t close on click.
  * Fixed axis buffer and filters interaction.
* [Bar chart](../../../visualize/viewers/bar-chart.md): reordered properties under the Order and Data submenus.
* Formula lines: regarding the line labels, the line title is on the plot, and both title and description are in the
  tooltip.
* We’ve added radio group identification to menu options where relevant.
* Linked tables:
  * Added the UI to specify multiple columns as key columns.
  * Fixed saving linked tables to the project.

### Enhancements in packages

#### [Chem](https://github.com/datagrok-ai/public/tree/7c62a0c018ec631d3b23760d538a17aaf4d4ca36/packages/Chem#readme)

We've added new **Scaffold Tree** visualization that organizes molecular data sets by arranging molecules into a tree
hierarchy based on their scaffolds. For details, see [Scaffold tree](../../../domains/chem/scaffold-tree.md).

* Improvements:
  * Added the package property to [set the default **Sketcher
    **](https://github.com/datagrok-ai/public/tree/master/packages/Chem#sketcher) so that users won’t have to switch on
    the first use manually.
  * Сhanged the result output for **Chem | Find MCS**. Now it returns a variable instead of a column.
  * Improved the handling of invalid molecules and empty inputs.
  * Added support for aromatic bonds when importing MOL2 files.
* Bug fixes:
  * Excluded blank data from **Diversity Search**.
  * Fixed filter: incorrect behavior after **Reset**.
  * Added progress indicator to **Similarity** and **Diversity Search**.
  * Fixed handling tables with a small number of data points when building **Chemical Space**.
  * Fixed 2d molecule orientation in the widget.
  * Fixed **Structural Alerts** detection memory leak.
  * Fixed Elemental Analysis working on MOLV3000.
  * Provided permanent fix for shared canvas re-allocation in the `rdkit-cell-renderer`.
  * Synchronized the selection and deselection of activity cliffs on scatter plot and grid.
  * Fixed an error after **Use as Filter** action under the cell with a molecule.

#### [Dendrogram](https://github.com/datagrok-ai/public/tree/master/packages/Dendrogram#readme)

We've separated **Dendrogram** from [**PhyloTreeViewer
**](https://github.com/datagrok-ai/public/tree/master/packages/PhyloTreeViewer). And now it’s
a [package](https://github.com/datagrok-ai/public/tree/master/packages/Dendrogram) for the Datagrok platform for
phylogenetic tree visualization.

Use Dendrogram viewer to:

* Display the Newick tree format files (NWK, NEWICK).
* Inject a dendrogram into the grid for interactive exploration of hierarchical clustering.

For details, see [Dendrogram](https://community.datagrok.ai/t/dendrogram/721).

#### [GIS](https://github.com/datagrok-ai/public/tree/master/packages/GIS)

We've retired Google Map viewer and implemented a Map viewer in
the GIS package. It shows geospatial data on a map as either markers or a heatmap. It displays data in geographic
formats, like GEOJSON, TOPOJSON, KML, and KMZ. You can also add a map viewer to your custom table. To learn more,
see [Map viewer](https://community.datagrok.ai/t/visualization-related-updates/521/27).

#### [Peptides](https://github.com/datagrok-ai/public/tree/master/packages/Peptides#readme)

This release introduces new features and capabilities, such as custom clustering and multiple views, as well as a
slightly redesigned user interface and improved application stability. To learn more,
see [Macromolecules updates](https://community.datagrok.ai/t/macromolecules-updates/661/11)

Improvement:

* Added aggregated columns to the Monomer-Positon viewer tooltips. This viewer is used
  in [Mutation Cliffs / Invariant Map and Most Potent Residues](../../../domains/bio/macromolecules#quantitative-sequence-based-activity-relationship-analyses-to-enable-design-and-optimization-of-polymer-modalities)

#### [SequenceTranslator](https://github.com/datagrok-ai/public/tree/master/packages/SequenceTranslator#readme)

* Improvements:
  * Added adaptive input fields for chains. They resize automatically upon adding longer sequences.
  * Relocated _direction input fields_ to the right of chain input fields.
  * Improved handling of the malformed data. The malformed string is now filled with red color.
  * Added tooltips.
  * Added routing to open specific tabs from links.
* Bug fixes:
  * Fixed rendering _resulting molecule_ with empty AS, AS2.
  * Fixed rendering for separate strands.

#### [Charts](https://github.com/datagrok-ai/public/tree/master/packages/Charts#readme)

We've merged package **Viewers** to **Charts**, removed obsolete viewers, and refactored the others:

* [Sankey viewer](https://github.com/datagrok-ai/public/tree/master/packages/Charts#sankey):
  * Added a tooltip.
  * Added the selection functionality.
  * Added viewer redrawing on filtering.
  * Fixed focus to highlight the current.object on a viewer
  * Fixed selection on filtered rows.
* [Chord viewer](https://github.com/datagrok-ai/public/tree/master/packages/Charts#chord):
  * Added viewer redrawing on filtering.
  * Fixed the selection functionality.
* [Radar viewer](../../../visualize/viewers/radar-viewer.md):
  * Fixed changing color issue: synchronized the color of percentiles on the legend and **Context Pane**.

#### [MLB](https://github.com/datagrok-ai/public/tree/master/packages/MolecularLiabilityBrowser#readme)

* VRs tree for grid.
* Revert database connection dataSource `PostgresNet`.
* Fix for routing tree in url.

### Improvements for developers

#### [Viewers](../../../develop/how-to/develop-custom-viewer.md)

* Added the ability to specify default viewer settings for the dataframe,
  see [annotation](https://github.com/datagrok-ai/public/issues/1395#issuecomment-1364325511).
* Added the ability to show a custom viewer in the Viewers section of the toolbox. The viewer should have an icon (set
  via `meta.icon tag`). Toolbox visibility can be specified as `meta.toolbox: true`.

#### [JS API](../../../develop/js-api.md)

* Improvements:
  * Added optional parameters to `Column.meta.colors.setLinear()` to set min and max values for linear color coding
    scale.You can also use the `COLOR_CODING_SCHEME_MIN` and `COLOR_CODING_SCHEME_MAX` tags.
  * Added radioGroup to Menu.item.
  * Added the ability to set `friendlyName` for `DataConnections` created from JS API.
  * Added hints to UI methods. The purpose is to attach interactive hints, close to what happens in tutorials, to
    elements. Application developers can then use these methods in the app code to introduce a new feature to the user,
    etc.
  * JS Viewers: disabled the removal of the default values (this likely fixes the issue with the JS properties not
    serialized).
  * Exposed schema property for `DataConnection`.
* Bug fixes:
  * `ui.choiceInput` does not work with strings that contain `\r` symbol.
  * Failure to get meta data from `SemanticValue` in JS.
  * Fixed date in `getValues()` in stock-broker.js.
  * Fixed `Project.removeChild` method.

#### [Widgets](https://github.com/datagrok-ai/public/tree/7c62a0c018ec631d3b23760d538a17aaf4d4ca36/packages/Widgets#readme)

* Added the ability to create a custom package settings widget.
* `Widget.getType()` method major code cleanup.

#### Enhancements in libraries

* [utils](https://github.com/datagrok-ai/public/tree/master/libraries/utils#readme):
  * Input categories added in the input form.
  * Added basic validation, btn disabling and UI fixes to **RichFunctionView**.
  * **RichFunctionView**: run button disabling added.
  * Move distance methods and constants to **ML lib** (utils), remove `index.ts`.
  * Exclude blank data from diversity search.
  * Utils for color.
* [ml](https://github.com/datagrok-ai/public/tree/master/libraries/ml#readme):
  * Moved distance methods from **Chem**, **Bio**, **MLB** to **ML lib**.
  * Fixed **Activity Cliffs** page freezing bug.
  * Fixed **Stochastic Proximity Embedding** dialog does not work.
  * Handled errors from web workers.
* [bio](https://github.com/datagrok-ai/public/tree/master/libraries/bio#readme):
  * Clear **Bio** lib dependency on **Phylocanvas**.
  * Separate **Dendrogram** package.
  * Fix dendrogram interfaces to allow emp.
  * Added `ITreeHelper.getNodesByLeaves`.

### [Compute](../../../compute/compute.md)

* Improvements:
  * `RichFunctionView` renamed and stabilized.
  * Added script generator app example.
  * Added sample of high-order function.
* Bug fixes:
  * Fixed linking to `RichFunctionView`.

### Other bug fixes

* Actions in packages don't support semantic values.
* Connections: Fixed anonymous mode.
* Databases | Queries: error occurs in `Inchi Key To Chembl`  query when setting an input dataframe.
* DataFrame: CSV parser: when there are no rows, space is auto-selected as a delimiter instead of comma.
* Datagrok page freezes after some time of inactivity.
* Grid does not clear the `mouseOver` status if the cursor is below the last row.
* `grok.shell.project` : `newProject.open()` doesn't change active project.
* Jupyter gateway falls on unrecognized message.
* Predictive Models deploy ID issue.
* Scatterplot: error when setting equal min and max axis ranges.
* Script loses package when user saves it.
* Search pane: textual search matches numerical columns.
* Sending request to docker container results in error.
* Table | Normalize Column does not work.
* Table content queries fire events.
* Tutorials: Scripting: opening a sample table is not counted as a completed step.
* Viewers: Statistics: Histograms: an error when trying to add histogram for the second time.
* Viewers: Statistics: Use in the tooltip does not work.
* Fixed molecules orientation for smiles in activity cliffs tooltip and **Context Pane**.
* Remove docked grid if viewer is closed.
* PowerGrid: Summary columns: normalization of bars.
* PowerGrid: fixed tooltip area.
* Fixed using DataConnection dynamic properties.
* Fixed handling empty values in panels.
* Fixed column selector dialog position in prop panel.
* Handling cases when sketcher friendly name is saved in user storage.
* Fixed **Use as Filter** with multiple columns.
* Fixed trellis plot full screen button.
* Box Plot: **Show Outside Values** does not work.
* Line chart `X` category axis is not aligned with value after filtering.
* `Event.path` removed from chrome.
* Fixed SAR grid didn't react to clicks.
* Fixed `deleteFile` method.
* Fixed possible null pointer exception.
* Fixed `readBinaryDataFrame` method.
* Fixed `DataConnectionParams` type.
* Fixed Filters: custom filters added via API fail during initialisation.
* Wrapped DataConnection parameters with MapProxy.
* Menu: do not close menu when checking/unchecking the item.
* Menu: indicate unchecked checkboxes.
* Minimised scatter plot issues: some rows are unexpectedly filtered out, scatter plot shown incorrectly when brought
  back.

### Misc

* Added an improved version of link-all script (grok link / grok unlink) to build and link dependencies of a package for
  local development.
* Added property description for property-bound menu items.
* Avoid keeping DB connection while sending email.

## 2023-01-06 Dev build 1.11.2

### Addressed Issues

* Save table linking to project
* #1412 column selector dialog position in prop panel
* closes #1344 Legend for continuous colorColumn
* Scatter plot coloured by numerical or date column: some data is missing after re-applying saved layout #1411
  Workaround
* Fixed connections leakage
* Merge
* Fixed shapes test
* Datlas: Test: Fix paths to tests in test.all.dart
* Revert "closes #1344 Legend for continuous colorColumn"

## 2022-12-23 Dev build 1.11.1

### Addressed Issues

## 2022-12-22 Dev build 1.11.0

### Addressed Issues

* #1236 Color coding: ability to invert colors for linear color coding
* Elemental analysis: minor fixes (columns hadn't been added if dataset already had the columns with such names)
* Environments: Remove top-menu from demo scripts
* #1282 Chem: fixed molecule size when drawing on canvas
* SequenceTranslator: Handle errors for rows with error while Save SDF
* Environments: Small demo script fix
* ST: Fix error message, bump version
* Implement Logger
* #1289: Chem \- place recent molecules on top of the list, saving the most recent coordinates
* Utils: fixed HistoryInput CSS
* Environments: Test for inline environments
* closes #1253 Formula lines not parrallel in log axes
* Chem: fixed sketcher tests
* Chem: #1293 short smiles detector fix
* Elemental analysis: provided an opportunity to add column with molecular formula
* Created index on name and namespace
* Fixed #1281: revert gasteiger charges panel to script-based
* #1297: temporary disable Info panel
* #1261: Chem \- clear sketcher button style fix
* Fixed package notification
* Chem: filtering incorrect recent/favorite molecules
* Chem: returned saving recent/favorites as molfiles
* closes #1292 Line chart is showing less data than expected if table is filtered
* dpetrov/GITGUB-1231: Fixed canvas resize handler by taking into account the element's offset
* ApiTests version update
* Adding auto check children for TreeView nodes
* Package JS-API: Added JS API for TreeViewNode
* GIS: gis-area handler improvement
* GIS: geocoding update
* Chem: fixed bug in addFavorite
* Chem: returned popup menu to similarity/diversity, structure 2d
* #1236 Color coding: ability to invert colors for linear color coding (WIP)
* Chem: handle cases when previous local storage included malformed molecules
* #1300: handling malformed target molecules in similarity search
* Remove delays in JS UI Tests  (WIP)
* Color settings, fix the view styles
* Fixed #1304: rendering only 3 symbols of monomer
* bio lib: Move newickToDf to ITreeHelper api and extend args, remove INewickHelper
* #1221 PhyloTreeViewer: Extend newickToDf args, move it to TreeHelper, remove NewickHelper, bump version
* Elemental analysis: removing option that adds a column with molecular formulas
* ScaffoldTree: changing function names
* Chem: activity cliffs \- remove duplicate code
* Librarise/ml: activity cliffs \- code cleaning
* Butina Cluster: make the algorithm work on the molblock input data
* Chem: similarity/diversity search \- remove ts-ignore
* Chem: rGRoupAnalysis: check conditions at function start
* Minor code cleanup
* Chem: minor code cleanup
* Preventing changing of the current object by expanding a property panel ()
* Chem: constants for units
* Chem: bump version
* #1221 MLB: show all trees / clones of current antigen with NetworkDiagram viewer
* bio lib: bump version to publish
* PhyloTreeViewer: Bump version to publish
* Chem: inchi sketcher bug fix
* closes #1279 axes limits, zooming issues
* Chem: fixed substructure filters test (for Ketcher)
* Set new 1.10.2 version
* PhyloTreeViewer: Downgrade datagrok-api dependency to 1.8.2, bump version to publish
* Butina Cluster: Handling NoneTypes that may occur if invalid data is present in dataset
* Chem: fixing substructure filters tests
* Adding a document with the instructions on how to run docker container on the datagrok instance
* Adding gif for docker_container doc
* Adding gif to the docker_containers.md
* Closes #1321: Chem: R-Group Analysis: unable to choose a column when there are two columns with molecules
* Libraries/chem-meta: added Reaction class to rdkit api
* Bio: updated ml dependency
* Chem: updated ml dependency
* ScaffoldTree: removing one parameter from the script
* ScaffoldTree: minor fixes
* Closes #1329: JS API: ability to highlight rows
* #1252 Color coding: min/max properties for continuous color schemes
* Peptides: settings subject type
* Peptides: selection rendering fix
* Peptides: fixed barchart not rendering on project load
* Peptides: fixed grid row height
* #1077: SAR viewer switch input interaction fix
* Peptides: code cleanup
* Peptides: visible columns fix
* Peptides: fixed settings not triggering changes
* Peptides: fixed settings not updating columns
* Peptides: stability enhancement
* #1119: refactoring WIP
* #1119: made viewers more standalone
* #1119: using raw data for clusters df creation
* #1119: fixed subtitution info redundant calculations
* Peptides: fixed monomer-position viewer not rendering
* Peptides: don't add columns pane if no columns present
* #1043: Invariant Map color coding
* #1117: Logo Summary Table enhancements
* Peptides: removed constant clusters column name
* Peptides: cluster tooltip fix
* Peptides: removed activity columns from 'Columns to include'
* Wiki: JS API: Custom viewers: documented filter, selection, and highlighting
* Ddt: unit tests: moved to folders, deleted obsolete
* Ddt: unit tests: func benchmarks
* Ddt: unit tests: minor fix
* Func: lazy initialization of the aux property
* Fixed a bug in the Balloon.error
* JS API: Added Balloon.closeAll()
* Functions: minor optimizations
* Help: Map Viewer \- WIP
* Help: Tree Map
* Help: Word cloud
* Help: Heatmap
* fix UI tests (awaitCheck)
* #1236, #1252
* #1330: temporary disable include columns feature
* Peptides: fixed missing value statistics calculation
* Avoid saving functions when called from socket
* Fixed analyzer warning
* Fixed #1332: devicePixelRation for WebLogo header cell renderer
* #1230 add as value does not work when tags are saved in df
* closes #1237 Scatter plot tooltip
* Peptides: release version 1.5.0
* #1330: enable column choice for Logo Summary table
* Peptides: release version 1.5.1
* #1310: Chem \- fixed similarity/diversity limit
* ST: highlight cells with invalid Type upon registration, initial version
* bio lib: Add generic parameter TNode for getLeafList and getNodeList
* Help | Visualize | Viewers | Heatmap: Edit documentation (WIP)
* PhyloTreeViewer: TreeHelper add generic param TNode for getLeafList and getNodeList according to ITreeHelper
* recloses #1209 Viewer legend is not synchronised with in-viewer filter in a specific case
* Datlas: WIP: Debug make configuration script for startup
* ST: minor fixes
* Datlas: Feat: Print debug messages for make configuration script on startup
* Range sliders hover area incorrect when legend is added
* Chem: Fix test top menu chem space/UMAP logic
* (Bug) Color coding: All does not work for boolean columns
* Jenkins: CI: Set UID and GID for jenkins user in puppeteer container
* #1236 Scatter plot support
* Rename PostgreSQL to dev and hide, rename PostgresNet to Postgres
* #1236 Update the invert icon
* JS API: add the special tags for continuous color schemes
* #1077: added shadow to Mutation Cliff cell renderer value
* Peptides: fixed selection not showing immediately
* #1236 Add the 'Scheme' label
* Update the generated files
* Fixes #1237 Scatter plot tooltip data values
* Chem: #1325 fix queries alignement (#1347)
* Fixed links
* Fixed links and trailing spaces
* ST: autostart fix
* DB: Test: HealthCheck start period
* #1236 Update the editor layout
* (Bug) First package publication on fresh database fails
* SequenceTranslator: Fix package-lock.json for @luma.gl verson 8.5.17
* #1345 Formula lines issues
* Chem: remove fix#1325
* GIS: nominatim docker file added
* Help: minor fix
* Utils: added feature request item to ComputationView
* Utils: feature request method override
* MSSQL: Test: Healthcheck start period
* MySQL: Build: microdnf instead of apt-get
* Oracle: Build: user to install
* (Bug) Manage: Users: current object does not change if you make a few consecutive clicks
* #1338: using float for activity delta
* (Bug) d4: Fix TreeMap test
* DSTK: Test: Skip tests until we figure out the new infrastructure
* Datagrok: Build: New version of Grok Connect
* Docker: Build: Use short commit IDs in image tag
* Datagrok: Build: Nginx proxy to 127.0.0.1
* Xamgle: Docs: Fix Help links
* Xample: Docs: Fix Help links
* Datagrok: CI: Olena Ahadzhanian in beta_users.csv
* Fixed possible null pointer exception
* Datlas: Test: Ignore Setup instructutions for chem tests
* Datlas: Test: Test File Connections on demo data path
* Datlas: Test: Add localFileSystemAccess option to connctors settings
* Datlas: Test: Increase timeout for inline scripts tests
* (Bug) Color Coding: When changing the linear color scheme, the null values get colored
* Grid: ability to apply color-coding to text instead of background
* Revert some changes to the color coding
* closes #1348 Minimised scatter plot issues
* Public submodule update
* (Bug) Scatterplot: error when setting equal min and max axis ranges
* Fixed deleteFile method
* Fixed readBinaryDataFrame method

## 2022-11-22 Dev build 1.9.0

### Addressed Issues

* bio lib: Fix properly naming for packages/GUIDE.MD
* Bio: Add VdRegionsViewer.positionHeight with default 'Entropy', fix args of exported package functions, fix
  currentView handling for WebLogo-positions tests
* JS UI Test: Predictive modeling (Caret) (WIP)
* JS UI Test: Chem | To Inchi
* JS UI Test: Fingerprints
* bio lib: Fixed phylocanvas.gl and deck.gl versions
* Issue #1123: RadarViewer. Add tests for the updated functionality WIP
* Library gridext: Fixed an issue the the duplicated Grid row header
* Packages PowerGrid: Bumped up the version of the dependency library gridext
* Bio: cell.renderer #1120: sequence Add property to highlight common
* JS UI Test: Groups
* Investigate/improve test manager performance
* DevTools: updated utills version, bump version
* Bio: updated utils library version, bump package version
* Chem: updated utils library version, bump package version
* ClinicalCase: updated utils library version, package version up
* Alation: updated utils lib version and package-test.ts
* ApiSamples: updated utils lib version and package-test.ts
* ApiTests, Arrow: updated utils lib version and package-test.ts
* extended tests with example of datetime column type and join (#1106)
* skip package filter test
* #1078: revert using virtual columns
* DevTools: test-manager \- ability to copy error text from property panel, version up
* Peptides: selection state saving fix
* PTV: Fixed phylocanvas.gl and deck.gl versions
* GIS: selection of rows-points
* JS UI Test: Tags
* Simplified icon methods
* #943 3dscatterplot onload error sometimes
* Utils: #466: Observable naming fix & patch version bump
* Utils: #466: added file input validators
* Utils: #466: fixed input DF uploading
* Utils: Validation class added
* Utils: minor version bump
* Utils: fixed validation error text
* Utils: patch version bump
* Chem: version up
* bio lib: Add interface property NodeType.isLeaf for PhyloTreeViewer, bump version
* SequenceTranslator: duplex registration
* (Bug) ColorPicker blinks in Chrome (WIP)
* DevTools: test manager \- fixed progress bar bug, version up
* Chem: remove substructure search from panel
* achopovsky/#933 to atomic level (#1133)
* bio lib: bump version
* closes 1110 Viewers: legend cannot be resized again
* SequenceTranslator: Axo Labs Pattern enhancement
* Bio, Peptides, MolecularLiabilityBrowser: raise version of lib bio in dependencies
* closes #1131 MultiColumn Selector: selects only columns filtered by search
* (Bug) Nginx: Add additional logging to investigate timeouts
* Docs: LE certificates remewal
* #1103 layout not saving
* (Bug) Grok Connect: the connection visibility
* Settings: mark instance as production (WIP)
* Removed obsolete packages
* Alation: enable Open With handling
* #943 3dscatterplot fix try
* Button alignment for form
* Issue #1034: Charts. RadarViewer removing unused options
* Issue #1123: RadarViewer. Adding tests for the case when options are changed
* Charts: Increase the version
* GIS: selection map points/grid rows (update)
* #609 aggregation mode for pc plot
* Bio: #264 monomer colours by properties
* closes #1134 trellis plot viewer icon
* Docs: Scenarios for release manual testing (WIP)
* Peptides build fix
* Introduce Viewer Serialization Context for saving layouts and projects
* (Bug) JS: DateTimeColumn allows to set string value
* (Bug) InfoPanels: Returning map or widget doesn't work
* (Bug) Package Entities Naming issues (WIP)
* Implement Logger (WIP)
* DevTools: test manager \- fixed tests, moved buttons from ribbon to view
* bio lib: Fix getStats for single line single monomer dataset, add test
* Issue #976: Add packages check
* Issue #976: Check whether elementalAnalysis had been done before in order to avoid error messages
* Chem: Increase version
* #1047: Bio \- substructure filters (helm filtering performance, modified linear search)
* #1118: distribution controls fix
* #1118: SAR mode switch fix
* #1118: collaborative filtering fix
* PhyloTreeViewer: PhylocanvasGL wrapper, injectTreeToGrid
* Converted const isPlainObject to function
* (Bug) Packages: Debug version duplicates release
* Careers section updates
* Issue #1052: Molstar: enable mol\* at the platform fixes
* GIS: autofocus on grid-row select
* Issue #1052: Adding loadStructureFromData
* datagrok-tools: version bump
* Bio: Add test bio.getStats for one line single monomer dataset, attach ticket link to test setRendererManually
* Fixed ORM test
* Fixed param name
* Bio: #1139 \- refactor library, get all data on loaded libraries
* column manager checkwrong columns after reordering
* closes #776 Viewer legend is not synchronized with in-viewer filter
* Bio: #1139 \- refactoir library, notify of monomers changed
* Charts: Surface Plot
* update deploy.md
* (Bug) TreeViewer: duplicated property "Animation Duration"
* SequenceTranslator: Axo Labs Pattern enhancement (WIP)
* Bio: substructure filters tests
* JS: Icons support for JS Viewers
* #943 scrollbars should not appear
* closes 932 Vertical scrollbar is not working in Order or Hide Columns menu
* Bio: #1139 refactor bio, export of capped monomers
* Fixed analyzer warnings
* HELM: add library as global
* Utils: history utils are extracted into separate namespace
* bio lib: Add getMonomerLib():IMonomerLib, fixes WebLogo for empty data
* Bio: Bump version for bio lib 5.5.0
* Helm: Add getMonomerLibObj package function for IMonomerLib implementation
* Fixed ambiguous column names
* Package functions: set helpUrl
* Onuf/grok 11298 external tutorials (#1116)
* Tutorials: update docs
* Issue #946: Adding ui for admet
* Tutorials: allow defining external tutorials and integrate them into the app
* Bio: Script admet-run.py to read ADMET models
* Peptides: Use IMonomerLib through bio lib
* bio: small fix
* Utils: added text filtering on historical runs
* Help: Docs: Terraform deploy instructions
* Help: Docs: Fix Terraform Template
* Help: Docs: Fix lonter errors
* GIS: box-selection bug fixed
* GIS: drag-box selecting issues fixing
* bio lib: Export MonomerWorks change getCappedMonomer(type, name), remove async init()
* Peptides: Use bio.MonomerWorks.getCappedMonomer()
* Help: Docs: Fix Terraform instruction
* Removed debug code
* Peptides: lib/bio version fix
* (Bug) Grid rollover disappears after applying filters or sort
* (Bug) Abnormal behavior when navigating grid by means of upwards arrow
* bio lib: Arrange exports, node types
* Bio: Fix import DG
* Chem: speed up chem space, activity cliffs, similarity search tests
* Tutorials: version bump
* 1128 Order and hide columns dialog isues part 1
* Adding pinned rows to Grid viewer
* Bio: modified test datasets for sequence space/activity cliffs
* PhyloTreeViewer: Fixes for PhylocanvasGlViewer, fix tree-cut-tests (as leaf lists)
* Utils: fixed bug with sharing runs' inputs/outputs
* Chem: fix test R-groups
* Bio: modified activity cliffs in ml library to be able to show modified lines grid
* Tutorials: docs update
* Peptides: analysis loading fix
* 709 filters sync in different views case 1, 2
* ORM: FuncCall visibility should depend on Func visibility (WIP)
* Packages tests
* Fixed syntax error
* Fixed analyzer warning
* Bump ApiTests version
* Bio: UI improvements, version up
* add images to provider conneectors
* fix image sizes for added connectoin images
* Utils: Ability to skip tests
* #1077: settings dialog UI
* (Bug) Simplify the chain for ViewerSerializationContext
* 1128 order or hide dialog(move top move bottom)
* closes 709 filterpanel sync (case with closing filter in different layout)
* Bio: ADMET python venv requirements
* Charts: version bump
* Helm: fixed freezing in substructure filter
* Chem: modified Chem space to call dimentionality reduce results from Bio
* Bio: sequence space on fingerprints
* GIS: detectors sampling of cathegories
* GIS: code optimised, some TODO's removed
* #709 dont close filters menu, just reset them
* Context menu commands: update conditions for JsViewers
* Charts: update readme
* added viewers test (layout)
* #1110 legend cats count incorrect
* closes #643 regression lines in log mode, r2 value
* Issue #1052: Adding pviz-bundle.min.js file
* #1077: Peptides settings
* DrugBank: version bump
* JS UI Test: Chem Sketcher
* Alation: tests fix
* #1152, Bio: work wth pdbs
* SQLite: version bump
* PubChemApi: version bump
* #1077: analysis start from top menu
* bio lib: NodeType, NodeCuttedType interfaces, fixes imports, PhylocanvasGL service interfaces and types
* Bio: Fix test substructureFilters/helm
* #1047: sequence activity cliffs using fingerprints
* GIS: color coding (linear/conditional). Filtering attempts
* #1047: Bio \- sequence activity cliffs using fingerprints
* #1047: Chem: added functions used in Bio to calculate activity cliffs
* Utils: #466: Historical runs with package isolation
* Add vertical alignment
* utils lib: Fix expectedObject for float tolerance
* PhyloTreeViewer: PhylocanvasGlService, fix for grid scroll events, skip render tasks by key, PhylocanvasGlViewer open
  test
* #1152
* Merged in onuf/1148-color-coding-extended-ui (pull request #282)
* #1077: fixed analysis dialog UI
* #1077: replaced SAR viewer title with inputs
* Wiki: update the color coding section
* closes #1156 filter bitset is reset after adding new column
* bio lib: Add IPhylocanvasGlViewer .onAfterRender and onHover, add ITreeHelper
* PhyloTreeViewer: Fix treeCutAsTree with tests, add TreeHelper, PhylocanvasGlViewer events onAfterRender and onHover
* MLB: main view fix routing, TreeBrowser render clone/tree with PhylocanvasGlService, TreeBrowser prevent recreating
  WebGL context, TreeBrowser show tooltip with long id, TwinPViewer NoSchemeItem
* added bool properties check to viewers test
* Added the tree viewer documentation
* Helm: Fix ES2021.String
* closes #1157 show visible only mode in Order and hide dialog
* Wiki: correct a misspelling
* Packages: NLP: add test files
* gridext: Fix remove tsc products
* SequenceTranslator: use sequence types map in registration
* Packages: NLP: add test skip reason
* bio lib: typings for @phylocanvas/phylocanvas.gl
* SequenceTranslator: fix of ribbon panel
* SequenceTranslator: helpers
* Chem: descriptors dlg fix
* SequenceTranslator: download function
* SequenceTranslator: registration constants
* GIS: filtering optimization, conditional color coding
* SequenceTranslator: separation of calculations from registration
* Lib/bio & HELM: MonomerWorks fix
* help: Add gif for PhylocanvasGlViewer
* PhyloTreeViewer: newickToDf set semType for columns node and parent to prevent detectors, add props nodeColumnName and
  parentColumnName to PhylocanvasGlViewer, sync selection between PhylocanvasGL and nwkDf in PhylocanvasGlViewer
* Bio: #264 \- refactoring
* JS UI Test: Form Viewer
* SequenceTranslator: code clarifications
* Grok Connect: CI: Version upgrade in scripts
* created missing indexes
* SequenceTranslator: minor change
* gridext lib: Up datagrok-api dependency version to 1.8.2, bump version
* #1155: Chem: substituent analysis in progress
* Chem: removed incorrect import from package.ts
* #709 #643 minor improvments
* PhyloTreeViewer: Fix GridWthTree and TreeHelper.setGridOrder for use leafColName and DataFrame tag '
  .newickLeafColumn', add tree generator, fix TreeHelper.cutTreeGrid to clean previous cluster marks, up used gridext
  dependency version
* Utils: historical runs API enhancement
* Added the SPGI demo dataset
* HELM: updated dependency
* Bio: #264 \- monomers read and calculated palettte creation
* #1125 add methods for getting value color into ColumnColorHelper (#1137)
* #1125 add methods for getting value color into ColumnColorHelper
* Release of Pinned Rows
* JS-API: Exposing pinned rows API
* Bio: fix init
* bio: version up
* Reorder Columns Dialog: improved a tooltip
* HomeView: got rid of unused code
* Got rid of debug printout
* Added missing files
* Resolves #1071 Color coding menu and dialog is not properly initialised after re-applying layout
* closes #1126 scatterplot packed cats with sparse columns
* Issue #1034: RadarViewer fixing problem with indicator
* JS UI Test: Leaflet Viewer
* Bio: Fix tests for similarity/diversity and activityCliffs, up Chem version in grokDependencies
* Library gridext: Added support for pinned rows.
* Added AppEvents.isProjectPublishing property
* TrellisPlot: code cleanup
* Grid: fix #1130: No tooltip on column header if column is too narrow
* DataFrameViewer: add onDataFrameChanged event
* Library gridext: Bumping up version
* #1148 Color coding: Disable the context menu command instead of error message when coloring type is not applicable
* Utils: fixed run double-load in ComputationView
* fix #1083: Grid: Popup Dialog from hamburger menu stays after the column is deleted
* fix #805: Grid: Change column type \- some strange appearance
* bio: add capPeptideMonomer
* Utils: favourites updated when filter is changed
* Utils: fixed order of the runs
* bio: bump version to 5.9.5
* bio: toAtomicLevel: fix
* bio: bump to 5.9.6
* bio: fix and raise version to 5.9.7
* closes #1161 Transformation Editor: use column selector instead of popup menus
* bio: refactor
* bio: bump version
* add property changed stream to interop
* add property changed stream to interop (#1168)
* (Bug) Files: Windows shares don't work
* Bio: refactor
* (Bug) Files: Harmonize credentials editing (WIP)
* Peptides: monomer tooltips
* closes 1167 zoom slider with specific data
* OligoBatchCalculator: separate constants
* OligoBatchCalculator: user group name and other constants
* GIS: synchronysing heat map with filteration
* added part about tests skipping
* GitHub 1159 Enable Grid column context menu in Pinned columns
* Help | Govern | User group: Add a description of the new functionality
* Chem: remove molecule without sketcher opening
* Charts: substituent analysis viewer in progress
* Peptides: import polishing
* 1157 Order & hide dialog checkbox replaced woth select
* 1161 transformation editor improvemrnts
* OligoBatchCalculator: calculations speedup
* OligoBatchCalculator: constants in order
* bio: monomer library enhances
* Bio: using monomer library
* Oauth: Bind custom google application
* gridext lib: add .gitignore .npmignore, bump version
* OligoBatchCalculator: adiitional modification validation start
* #1155: Chem: substituent analysis viewer in progress
* bio: eror fixed
* Fixed null pointer exception
* CSS fixes
* Helm: bio related refactor
* #1155: Chem: substituent analysis \- fixing bugs,
* Fixed attirbute exception
* #1077: replaced barchart with weblogo
* Peptides: settings sliders fix
* #1077: include columnds
* GitHub 1159 Multiple feature requests for Pin Column sort
* JS-API: notify for copyFrom and invert methods
* #1166: collaborative filtering fix
* Library gridext: Adressing sorting issues in pinned columns
* Library gridext: Updating the version of JS-API
* Package PowerGrid: Bumping up package version, JS-API, gridext
* fix #1173: Significant figures should be rounded properly
* Utils: report to audit about failing tests
* GIS: new colorcoding type, code optimization
* OligoBatchCalculator: validation of changes to additional modifications
* Closes #1163: PhyloTreeViewer TreeForGrid renderer (optimized for visible rows)
* OligoBatcchCalculator: show NaN if sequence is not valid
* Chem: tests
* PhyloTreeViewer: Bump version to publish npm
* Chem: bump version
* Bio: Bump version for dependencies and to publish npm
* Oauth setup
* #1161 transformation editor improvemrnts
* #1157 Order & hide dialog checkbox alignment
* Bio: fix manage libraries
* closes #1174 Can not select empty value colored data points in scatterplot
* #1119: moved out mutation cliffs calculation
* add clickhouse connector image
* Peptides: fixed multiple analysis starts
* GIS: colorcoding with data stored in cells
* Help | Visualize | Viewers | TreeViewer: Create documentation
* Help | Visualize | Viewers | Radar Chart: Create documentation
* closes #1186 Calculated columns dialog is too large in some cases (#1187)
* Dkryvchuk/add new col 2 (#1188)
* #1128 Order or hide dialog not working correctly in chrome
* (Bug) DG.Logger.audit fails with entity params
* JS: Introduce ObjectColumn class
* Packages: PowerPack: update viewer categories
* add to trellis look to js api
* Dkryvchuk/look trellisable (#1189)
* Packages: Charts: update help URLs
* Bio: fix if 0 libraries chosen
* Ability to set custom url for app
* Chem: test fix
* #1117: filtering out unimprotant clusters (naive)
* CI: Build dockerfiles for packages in GitHub Actions
* (Bug) Grid's Hamburger menu jumps around cursor when resizing columns
* Library gridext: Enabled pinned column resize functionality
* GIS: detectors optimization, colorcoding issue fixed
* (Bug) Parameterized Queries: Parameters with {choices} not working
* bio lib: Calculate .alphabetSize tag on demand
* Closes #1171 Bio: enhance detectMacromolecule speed
* Library gridext: Bumping Up the library's version
* Select tooltip columns: ability to reorder columns (WIP)
* (Bug) Line chart | Overview: Exception if you select a column different from the one selected for display
* Package PowerGrid: Bumping Up the package's version
* (Bug) Queries: scalar queries do not work
* PowerPack: improved error message
* (Bug) Chem: molecules are not rendered if you open a file immediately after starting the platform
* ML: Missing Values Imputation: only show columns with missing values
* JS API: Added SurfacePlot to the list of known viewers
* ApiSamples: add SurfacePlot sample
* (Bug) Packages: if a package's detector.js is not loaded, the whole semantic type detection system stops working
* Packages: deactivate semantic type detector if it throws an error
* SequenceTranslator: remove errorExist flag
* #1077: Header WebLogo selection rendering
* Surface plot: datetime fix
* #1077: header WebLogo interactions
* #1077: fixed most potent residues viewer title
* #1117: rendering mode switch
* #1077: monomer-position viewer cell value color fix (WIP)
* Chem: refinements
* Peptides: settings dialog fix
* #1117: columns aggregation
* JS UI Test: Project Upload
* Remove delays in JS UI Tests  (WIP)
* Bio: Fix import from lib utils, bump version for dependencies and publish to npm
* (Bug) JS: prepare() takes linear time when column is a parameter
* Grok Spawner: Docker Swarm support
* Various issues with PC-PLOT (WIP)
* (Bug) Swagger: Incorrect default value parsing
* Peptides: stats calculation fix
* Peptides: logo monomer ordering fix
* Peptides: version 1.4.0
* Docs: Create Help pages for every provider
* #1155: substituent analysis viewer \- fixing bugs
* #1155: substituent analysis \- fixing bugs
* (Bug) Charts | Radar viewer: the wrong axis occurs
* (Bug) Charts | Radar viewer: doesn’t respond to the column deletion in the grid
* (Bug) JS: check permissions method always returns true
* JS API Tests: connections: create, edit, delete, share
* Summary columns don't show after refreshing the page #930
* Issue #534: exposed TableQuery and DbTableQueryBuilder
* Issue #534: docs for TableQuery and DbTableQueryBuilder
* #534: TableQuery & DbTableQueryBuilder tests WIP
* #534: TableQuery & DbTableQueryBuilder tests & fixes WIP
* (Bug) Charts | Radar chart: 'Show All Rows' property issue
* (Bug) Charts | Radar viewer: toggling 'Show Min' and 'Show Max' issues
* RadarViewer: Updating init function
* RadarViewer: removing div id
* Datagrok Docker fail on errors during build
* GIS: filtering bug fixed, filtering optimization, few unit tests added
* JS-API: removed duplicate TableQuery class
* Meta: added ChemDraw
* Jenkins: CI: Add Skipped characteristic to junit tests converter
* Jenkins: CI: Upgrade tools in puppeteer image
* ST: triplex
* ST: reg fixes
* (Bug) Charts | Radar viewer: incorrectly displays negative values (WIP)
* Compute: renamed catalog feature
* (Bug) Charts | Radar viewer: the number of axes should be limited
* Charts: Added radar-viewer-test to package-test
* Resolves #1204 JS API: add Column.meta.format helper
* #1204 JS API: add Column.meta.format helper
* (Bug) Multiple Issues with Chem functionality in Grid
* JS API Tests: TableQuery: create, run (on different query types) (WIP)
* SequenceTranslator: parse dimer/triplex
* Semantic type detectors performance test
* SequenceTranslator: registration: icons instead of buttons, set row height
* UI Tests fixes
* GIS: code optimization, unit tests
* SequenceTranslator: refactoring with isOverhang function
* SequenceTranslator: move functions used in draw-svg to helpers
* bio lib: Add getAlphabet by name, fix NotationConverter.convertHelm, fix UnitsHandler separator property allowing
  undefined, bump version
* #1204 JS API: update Column.meta.format
* Lstolbov/biostructure/render (#1205)
* Meta: updated chemDraw package version
* Closes #1192 Bio: detectMacromolecule benchmark tests
* Bio: detectMacromolecule benchmark tests skipReason for FASTA notation tests
* Bio: skipReason for Bio.detectorsBenchmark.separatorDnaLong1e6Few50
* Fixed tests
* (Bug) JS: Functions.eval output skips toJs
* (Bug) JS: Fix undefined and null confusion in toJs
* JS: Dapi count() method
* Skipped tests
* add new form class
* Bio: Fix splitters and renderers tests, remove force Bio package init
* Test: CI: Skip External provider: small stress test
* Tests: Skip failing group tests
* Fixed categories
* (Bug) PowerPack: initTemplates() test does not work
* Utils: #466: added About menu item
* JKG: Feat: Upgrade components to the latest versions
* Heathchecks for basic images
* generate python requirements file
* Grok Compute: uograde python
* Exposed TableQuery and DbTableQueryBuilder methods
* Added missing extension
* SequenceTranslator: draw-svg: object with Y coordinates
* SequenceTranslator: filter instead of counter using for loop
* (Bug) Chem: Filters: when there are multiple molecule columns, all filters default to the first one
* #709 trying to preserve both state and column name
* Utils: #466: getAbout func updated
* Docker: CI: KyotoCabinet version based on python 3.8

## 2022-11-08 Dev build 1.8.3

### Addressed Issues

* com
* Fix default value for select (choice) (remove '' and "" from string when chosen default value)
* Implement setRawData for bool columns
* render mols vertically on scatterplot x axis
* lrellis plot not working in LS
* Fixed the method for default choices
* delete useless prints
* bar chart slider positioning
* Grok Spawner: fix default value for GrokSpawnerAddress
* Make modals resizable
* Legends Resizing
* closes #810 Reordered selected columns
* add started/finished/author to funccall
* Layouts | Save to gallery: newly saved layout should appear in the property panel
* Layouts pane: ability to put the layout on the property panel
* CI: Deploy process documentation
* closes #798 selecting dots in jittered scatterplot
* #828 formula lines modal fixed
* #828 order/hide dialog rewworked
* #828 make modal not resizable by default
* JS API: ui.fileBrowser
* Work in progress.
* #831: HitTriage: work in progress
* Moved InputBase descendants to separate files.
* AppEvents.onInputCreated global event
* Ability to reuse database browser, and specify how results are interpreted
* Table Input: ability to get table by executing a database query
* Docs: Publish version from dev to Docker Hub
* Datlas: queue interaction for Grok Connect (WIP)
* Docker image for packages tests
* add compact view to scatterplot 1.0
* Scatterplot compact view
* JS API: View.statusBarPanels (WIP)
* (Bug) Grok connect: decimal, real, and double tests failed
* Grok connect: refactoring of numeric type tests
* JS: UI Forms harmonization (WIP)
* (Bug) Failed to delete package
* DB Indexes
* Fixed filename
* Fixed CSS
* Added a demo dataset with image urls
* PowerGrid: help on linked images
* JS UI Test: Scripting in Add new column (WIP)
* Functions: show function class instead of null in status bar
* Grok connect: refactoring of numeric type tests
* DataFrame.setTag(...) now returns this.
* Renderers: ability to specify column tags as conditions
* adds menu click method
* modal click simulation method adjustment
* Datlas connections leaking investigation
* Datlas: remove trailing ; in the row checking query
* JS API: Menu.click method
* #831: HitTriage: work in progress \- refactoring
* closes #800 Form does not work after filtering
* Introduced Entity.package Implemented safe delete New permissions check Derived properties support in ORM
* Introduced Entity.package
* GH-686: TreeView events
* modal placement, dont leave page on backspace
* GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Fixed deadlock in exchangeDevKey method
* CSS fixes
* Introduced Entity.package, reviewed Entities list
* Introduced Entity.package, introduced PackageEntityMixin
* Autogenerated code
* Speed-up FilesView start
* Fixed analyzer warnings
* Fixed package.functions section exception
* closes #791 Can not uncheck "Filter out missing values" after saving
* Fixed wrong naming on package deploy
* Removed vis.js from autostart
* Async external connectors fetch
* Speed-up DG start
* FilesView: fix file appearance bug in the card view
* Introduced Entity.package, bugfixes
* Fixed possible NPE
* JS API: BigInt support
* closes #857 Deselection in the scatter plot
* Fixed author parsing
* Package content validation: improve the up-to-date package bundle check
* Fixed tableInput icons
* closes #859 3D scatter plot: implement saving to .png
* Updated beta users
* Tested renderer selection by multiple tags (works fine)
* (Bug) Integration tests: Fix all northwinds for 'External provider: test all Northwinds' test
* refactor: configuration files
* Docs: fix documentation links
* Public submodule
* Add delete FuncCall api
* PackagesView installed filter
* CSS harmonization
* Introduced EntityType.isPackageEntity field
* Fixes #923: Grid: Order or hide columns: provide "move to the top" / "move to the bottom" buttons
* Viewer-specific tooltips (WIP)
* Min() and Max() are not compatible with date and datetime columns #689
* Fixed tests
* removed debug lines
* #closes 858 Cell as an input for widget
* (Bug) FilesView: dock manager exception on opening
* Setup fixes
* Update doc for proxy
* Better transaction management
* Added debug lines
* Some packages not publishing fix
* Removed debug printout
* Add save as png to all viewers
* (Bug) Connection is cloned on sharing
* Reverted GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Dialog: fixed a bug with the resizing
* Grid: Fixed an issue with bool cell renderer
* dapi.files: fixed a bug by reverting to the lazy initialization. Code cleanup.
* fixed: ColumnsInput: an exception when checkedColumns is not specified
* js-api update
* (Bug) Server push is sent twice
* (Bug) Row.setValues doesn't fire onDataChanged event
* Removed debug lines
* Add exceljs to sources
* DateTime scalar functions harmonization
* Fixed shell test
* PowerPack: build fix
* Input fixes, and new disabled styles
* Fixed scalar functions null behavior
* Fixed logical functions in Chrome
* Introduce DateInput with dayjs support
* Fixed deploy exception
* Fixed label css
* (Bug) ScatterPlot: Top and Bottom legend position cuts labels
* (Bug) Selection doesn't work when jitter activated
* ML | Random Data doesn't work ML | Missing Values Imputation doesn't work
* CSS fix
* Fixed dapi and ml tests
* Possible NPE fixed
* Scatterplot 3d detach method
* dapi.root
* View close fixed
* Updated public
* Test: setUpAll and tearDownAll with groups for integration tests
* Test: Skip tests which fail everytime
* Added JS wrapper for EntityRecord
* Packed categories zoom by filter
* package_id migration tuning
* Fixed user login error
* (Bug) FilesUploader: ?layout parameter doesn't work
* CI: Sofiia Podolskaia in beta users
* Test: Use custom session token for integration tests Load session token for integration tests from environment
  variable GROK_SESSION_TOKEN As a fallback value the old session token unit_test_token will be used
* Refactor: Complete refactoring of deploy folder Deprecate old scripts
* Test: File with environment variables for local tests
* GH-954: expanded accessor for TreeViewGroup
* (Bug) Integration tests: External provider (WIP)
* Docs: Instructions to run datlas tests locally on local custom stand
* Test: Add group for all integration tests Fix integration tests which were failing because of session
* closes #936 pie chart improvments
* Test: Add testDevKey config option for docker stands
* Update background for macromolecules
* Test: fix Missing Values Imputation: k-Nearest Neighbour Imputation
* Fixed Datlas tests
* Reverted incorrect fixes
* GH-965: TreeViewGroup.group now also accepts Element
* closes #961 gray out unused categories instead of not showing them
* (Bug) ScatterPlot: Categories rendering issue
* closes #962 Save as PNG: not all elements are added to image
* closes 790 Reset filter doesn't work on "Filter out missing values"
* add button style
* closes #963 Continuous selection on shift+click
* CI: Add short commit to docker image
* Test: Fix Datlas tests
* Load package before loading tests
* Balloons: add the "X" icon
* closes interfaces generation
* (Bug) When user shares an entity from package package is shared too
* Chem harmonization
* closes 935 treemap improvments
* Chem: fix fingerprints
* PC Plot: collaborative filtering
* Test: Increase timeout for datlas tests
* (Bug) Packages permission check doesn't work
* closes 935 with all rowSource modes
* closes 945 Line chart: ability to hide and make semi-transparent lines
* FilesView: add the "Add new share" ribbon icon and link below the tree
* FilesView: 'refresh' icon
* CI: Docker Buildx services auto list
* CI: Docker Demo Chembl DB
* CI: Docker Buildx target
* CI: Create test_numeric_types table for external_provider_test
* Test: Increase timeout for Dataset Client: Projects: CRUD: Remove
* Test: Increase timeout for Datlas tests
* Integration tests: datetime test
* Line chart: add tooltips in multiaxis mode
* New FontAwesome version
* Fixed css
* Fix bakcground styles, for macromolecules page
* Chem: restoring sim search, add indeces
* code fix for #935 treemap improvments
* Pipe file downloads
* Folder sharing
* Ask credentials on file share open
* Minor style changes
* Sketch View: add the reset button
* JS API: Grid.sortByColumns, Grid.sortTypes
* Chem: fix for similarity search
* Grok connect: performance test (select 1)
* Integration tests: refactor external_provider_perf_test
* Fixed AD login
* GH-955: subDir and share methods
* GH-955: js api update
* (Bug) Files: Unable to open context menu for files in explorer view
* Credentials Manager
* (Bug) connection.shares.connection == null condition doesn't work
* BarChart: support additive functions 'unique', 'missing value count' and 'value count' in stacked mode (WIP)
* Docker: Datagrok ca-certficates
* Fixed files API
* CI: Add Alexey Chopovsky to beta users
* Docker: Remove elastic from datagrok
* Docs: Upgrade docker instructions with certificates
* PowerGrid: bumped up version to 1.1.0
* #969 forrula-lines vert lines
* Reversed wrong commit
* Public token
* Chem #995: clean up -Similarity
* Chem #995: open phacts clean up
* Enable saving multiple credentials in credential manager
* Fixed logout
* Docs: Alllow askalkin@datagrok.ai to view encrypted files
* (Bug) JS: grok.functions.call doesn't handle error
* JS API: correct return type for viewer.getOptions
* #953: Filters turned off status is not saved in layout
* #966: Deleted column filter should not be filtered on
* #969 new classes system ,aybe should be reworked
* Change order of interop type checks
* #969 linechart hittest
* closes #996 Tooltip not working for partly truncated values
* Funcs: JS-based context functions for multiple columns
* (Bug) ColumnsInput: number of selected columns is not shown upon construction
* JS API: GridColumn.move method
* closes #1015: PowerGrid: Sparklines: Ability to add sparklines from the "Actions" pane for selected columns
* JS API: GridColumn.scrollIntoView, RangeSlider.scrollBy, RangeSlider.scrollTo
* PowerGrid: Sparklines: accessing settings from the hamburger menu
* #930: PowerGrid: sparklines: ability to rename and remove summary columns
* #930: PowerGrid: sparklines: managing column order
* #930: PowerGrid: sparklines: setting the sparkline column as a current object after it is created
* Possible NPE fix
* Deprecated old tests
* Chem: #454 descriptors
* Chem: small fix
* Fixed package init
* Datlas: fixed project info retrieval
* #1000 viewport exists
* #1000 viewport apeears properly, style changed
* dialog cliptoscreen onresize
* #1000 bands
* Remove labs submodule
* (Bug) Custom context menu items for JS widgets and viewers cannot be registered via API
* (Bug) Dropbox provider doesn't work
* Chem: fix descriptors test
* Chem: #995 removing r-groups
* Chem: #454 compute side polishing
* Update API and public token
* (Bug) Proxy to docker container doesn't work
* Wiki: replace links
* (Bug) Package Credentials Editor doesn't work
* JS API: code cleanup and minor refactoring
* Histogram: reduced the marker size to 4 pixels
* Better wording
* PowerPack: testing \- WIP
* closes #1022 Tooltip: do not show value while using custom renderers
* closes 1021 matrix plot from layout
* JS API: Point.distanceTo, Rect.containsPoint
* dialog resizing before dom created fixwd
* Update beta users list
* (Bug) grok.shell.registerViewer fails
* Minor code cleanup
* (Bug) Ability to clone the connection (WIP)
* Grok Spawner: more information for every request in return
* Docs: devops docs
* closes #1001 matrix plot inner viewer settings
* JS API: update param type in Column.aggregate
* Grok Spawner: docker entrypoint script
* Chem: rgroups exclusion
* #1001 syncing tables in wrong place
* PowerGrid: version bump up
* Sketch View: update UI
* closes #657 line chart selectors with special symbols
* (Bug) Substructure Filter: Disabled if you drag the column to the filter
* closes #856 Chem: add molecule filter to a filter-group by default
* Packages: disable certain sources checks for non-webpack packages
* Chem: descriptors fixes
* (Bug) Method click() throws an exception when called for top menu item
* (Bug) JS: disable Uint8List convertion
* Implement a mechanism to test functions
* Hierarchical pipe-delimited info panels (example: Chembl | Substructure)
* Fixed #1029: File Manager: clicking the file in the Folder Tree doesn't trigger file preview in the main view
* JS: DataFrame.toByteArray method
* Harmonize membership editor
* Package content validation: check package properties' types
* (Bug) Scripts: Sample table does not open
* (Bug) __jsThen error during some packages init
* Docking: highlighting the separator
* Add a script to find internal links in help
* Add the function namespace to tests
* Grid: automatic format for small floats
* Chem: api fix
* (Bug) Grok Spawner: After image creation the container does not start
* CI: Add time spent to junit convertion script
* (Bug) Make loadKeyPair method more robust
* Package permissions fix
* Updated beta users list
* Merge 'master'
* CI: add jq to puppeteer image
* Dockerfiles: fixed container address retrieval
* update img
* Update img
* Docs: onboarding docs refactoring
* Grok Spawner: publish ports on local docker
* Test: mute failing tests
* Refactor Secret Manager connection
* (Bug) Core: MapInput cannot be added to a dialog for function parameters
* Grok Spawner: network properties in responses
* CI: Junit test results converter
* 1002 adds axes to matrix plot
* #1002 inner viewwer setting icon not merged
* Improve start mode selection
* Dockerfiles: get container port
* ##1042: Onboarding documentation: moved things around, minor cleanup
* CI: Add tests manually to junit reports
* Dashboard Projects (WIP)
* #1002 matrix plot axes code refactored
* CI: Convert JSON files to Junit reports
* CI: Add Rsync to pepputeer base image
* PowerGrid: version 1.1.2
* CI: Correct messages to jest tests failures
* Provided a hook to borrow Hamburger menu by external components
* #943 basic axes
* Removed grid column reset for Hamburger menu
* CI: Append timeout and evaluation errors to junit reports
* closes #953 Filters turned off status is not saved in layout
* CI: Include multiline stdout to Junit reports
* CI: Add author to Junit tests output
* JKG-Add-libmamba-to-conda-resolver-mechism
* Sketch View: remove the icon for adding viewers when the view is opened to edit a tooltip
* Docs: Access to servers using SSH
* Added exportable indices bounds for visible rows and columns
* Added xangle registration for export of visible rows and columns
* #943 3D scatter plot enhancements (log coordinates, some ui moments)
* Add new method (RangeSlider_SetShowHandles) \- need for on/off resize in slider
* fixes choice unput value not screened
* #closes 1075
* PropertyGrid editors: moved each editor to a separate file.
* Viewers: ability to pre-aggregate source DataFrame (WIP)
* (Bug) Scripting: JS Script fails, when contains more than 6 params
* (Bug) Can't share a folder from Home connection
* (Bug) Project duplicates Data Query name
* PC Plot: normalization mode
* PC Plot: auto-hide unused filters
* PC Plot: Show value axis when the scale is global
* RangeSlider: an option to not show the bar when the slider is expanded
* closes #1002 axes for trellis plot
* closes #1087 Use custom renderers with scatterplot label
* Core: add html2canvas library
* closes #1072 when context menu and column selector popup menus overlap, the click event triggers on the wrong menu
  item
* dpetrov/GitHub-1083 Made Actions closable in the Hamburger menu
* ColumnComboBox: not opening the drop-down part when right-clicked
* #943 use existing methods for axes tickmarks
* labels refactoring
* Dockerfile methods WIP
* Dockerfiles containerAddress fix
* DockerfilesClient proxyRequest method fix
* GH-863: Dockerfiles JS API update
* #1096 x axis fixed
* #1096 y labels fix
* Proper packages names
* (Bug) Package Entities Naming issues (WIP)
* (Bug) Script icon incorrect
* dpetrov/GitHub-1093 Added virtual columns support for grid's columns manager
* #1096 fixes multiple columns labels
* dpetrov/GitHub-1093 Minor fixes for the grid's columns manager
* Jenkins: Docs: Add Jenkins setup information
* Jenkins: Docs: Add password information for Jenkins
* (Bug) Query Transformation name issue
* #1098 add hittesting and selection to pc plot
* JS API: Custom inputs (WIP)
* dpetrov/GitHub-1104: Fixed returning null view field in ViewLayout object
* Merged in onuf/color-coding-update-bug (pull request #270)
* Group passwords (WIP)
* #601 PC plot improvements(add interactivity)
* #928 context menus
* 601 PC plot improvements(pc plot dies when transformation reset)
* #943 axes improvments
* Merged in onuf/color-coding-update-bug (pull request #273)
* update public token
* #943 issue with axes on otherinstances
* add Dmytro Nahovskyi to list
* (Bug) Right Click on the Grid's Row Header generates exception
* 1103 fixes filters in layout behavior
* change grok conneect version 1.5
* Update beta users table
* Charts: Radar: clarified property description
* Fix button styles
* Ability to load entitty from inactive package version
* MultiColumn Selector \- use subsets of columns when needed  (WIP)
* Error when serealizing df with virtual column
* Viewers: "General | Open in workspace" option
* Databases: groom the providers (WIP)
* User can't log in is login contains hyphens
* (Bug) Pinned Columns: row number column appears twice after layout application
* requirements.R default value for plsdepot has been returned
* buildx_docker_image.sh returned default value for cache_push (true)
* deploy/demo_db_clickhouse/demo_db_clickhouse:
  -f772192786c7f01c4f7afacfe68affa1dd2fe9a9-f772192786c7f01c4f7afacfe68affa1dd2fe9a9.docker.tar edited online with
  Bitbucket
* closes 1126 scatterPlot packed cats with no values
* #609 pc plot filter in aggr mode
* skip package filter test
* Simplified icon methods
* #943 3dscatterplot onload error sometimes
* (Bug) ColorPicker blinks in Chrome (WIP)
* Chem: remove substructure search from panel
* closes 1110 Viewers: legend cannot be resized again
* closes #1131 MultiColumn Selector: selects only columns filtered by search
* (Bug) Nginx: Add additional logging to investigate timeouts
* Docs: LE certificates remewal
* #1103 layout not saving
* (Bug) Grok Connect: the connection visibility
* Settings: mark instance as production (WIP)
* #943 3dscatterplot fix try
* Button alignment for form
* #609 aggregation mode for pc plot
* closes #1134 trellis plot viewer icon
* Introduce Viewer Serialization Context for saving layouts and projects
* (Bug) InfoPanels: Returning map or widget doesn't work
* Implement Logger (WIP)
* (Bug) Packages: Debug version duplicates release
* Careers section updates
* Fixed ORM test
* column manager checkwrong columns after reordering
* closes #776 Viewer legend is not synchronized with in-viewer filter
* update deploy.md
* JS: Icons support for JS Viewers
* #943 scrollbars should not appear
* closes 932 Vertical scrollbar is not working in Order or Hide Columns menu
* Fixed ambiguous column names
* Package functions: set helpUrl
* Removed debug code
* (Bug) Grid rollover disappears after applying filters or sort
* (Bug) Abnormal behavior when navigating grid by means of upwards arrow
* 1128 Order and hide columns dialog isues part 1
* Adding pinned rows to Grid viewer (WIP)
* 709 filters sync in different views case 1, 2
* ORM: FuncCall visibility should depend on Func visibility (WIP)
* Fixed syntax error
* Fixed analyzer warning
* add images to provider conneectors
* fix image sizes for added connectoin images
* (Bug) Simplify the chain for ViewerSerializationContext
* 1128 order or hide dialog(move top move bottom)
* closes 709 filterpanel sync (case with closing filter in different layout)
* #709 dont close filters menu, just reset them
* Context menu commands: update conditions for JsViewers
* #1110 legend cats count incorrect
* closes #643 regression lines in log mode, r2 value
* Merged in onuf/1148-color-coding-extended-ui (pull request #282)
* closes #1156 filter bitset is reset after adding new column
* closes #1157 show visible only mode in Order and hide dialog
* Chem: descriptors dlg fix
* created missing indexes
* #709 #643 minor improvments
* Added the SPGI demo dataset
* #1125 add methods for getting value color into ColumnColorHelper
* Release of Pinned Rows
* Reorder Columns Dialog: improved a tooltip
* HomeView: got rid of unused code
* Got rid of debug printout
* Added missing files
* Resolves #1071 Color coding menu and dialog is not properly initialised after re-applying layout
* closes #1126 scatterplot packed cats with sparse columns
* Added AppEvents.isProjectPublishing property
* TrellisPlot: code cleanup
* Grid: fix #1130: No tooltip on column header if column is too narrow
* DataFrameViewer: add onDataFrameChanged event
* #1148 Color coding: Disable the context menu command instead of error message when coloring type is not applicable
* fix #1083: Grid: Popup Dialog from hamburger menu stays after the column is deleted
* fix #805: Grid: Change column type \- some strange appearance
* closes #1161 Transformation Editor: use column selector instead of popup menus
* add property changed stream to interop
* (Bug) Files: Windows shares don't work
* (Bug) Files: Harmonize credentials editing (WIP)
* closes 1167 zoom slider with specific data
* GitHub 1159 Enable Grid column context menu in Pinned columns
* 1157 Order & hide dialog checkbox replaced woth select
* 1161 transformation editor improvemrnts
* Fixed null pointer exception
* Fixed attirbute exception

## 2022-11-02 Dev build 1.8.2

### Addressed Issues

* com
* Fix default value for select (choice) (remove '' and "" from string when chosen default value)
* Implement setRawData for bool columns
* render mols vertically on scatterplot x axis
* lrellis plot not working in LS
* Fixed the method for default choices
* delete useless prints
* bar chart slider positioning
* Grok Spawner: fix default value for GrokSpawnerAddress
* Make modals resizable
* Legends Resizing
* closes #810 Reordered selected columns
* add started/finished/author to funccall
* Layouts | Save to gallery: newly saved layout should appear in the property panel
* Layouts pane: ability to put the layout on the property panel
* CI: Deploy process documentation
* closes #798 selecting dots in jittered scatterplot
* #828 formula lines modal fixed
* #828 order/hide dialog rewworked
* #828 make modal not resizable by default
* JS API: ui.fileBrowser
* Work in progress.
* #831: HitTriage: work in progress
* Moved InputBase descendants to separate files.
* AppEvents.onInputCreated global event
* Ability to reuse database browser, and specify how results are interpreted
* Table Input: ability to get table by executing a database query
* Docs: Publish version from dev to Docker Hub
* Datlas: queue interaction for Grok Connect (WIP)
* Docker image for packages tests
* add compact view to scatterplot 1.0
* Scatterplot compact view
* JS API: View.statusBarPanels (WIP)
* (Bug) Grok connect: decimal, real, and double tests failed
* Grok connect: refactoring of numeric type tests
* JS: UI Forms harmonization (WIP)
* (Bug) Failed to delete package
* DB Indexes
* Fixed filename
* Fixed CSS
* Added a demo dataset with image urls
* PowerGrid: help on linked images
* JS UI Test: Scripting in Add new column (WIP)
* Functions: show function class instead of null in status bar
* Grok connect: refactoring of numeric type tests
* DataFrame.setTag(...) now returns this.
* Renderers: ability to specify column tags as conditions
* adds menu click method
* modal click simulation method adjustment
* Datlas connections leaking investigation
* Datlas: remove trailing ; in the row checking query
* JS API: Menu.click method
* #831: HitTriage: work in progress \- refactoring
* closes #800 Form does not work after filtering
* Introduced Entity.package Implemented safe delete New permissions check Derived properties support in ORM
* Introduced Entity.package
* GH-686: TreeView events
* modal placement, dont leave page on backspace
* GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Fixed deadlock in exchangeDevKey method
* CSS fixes
* Introduced Entity.package, reviewed Entities list
* Introduced Entity.package, introduced PackageEntityMixin
* Autogenerated code
* Speed-up FilesView start
* Fixed analyzer warnings
* Fixed package.functions section exception
* closes #791 Can not uncheck "Filter out missing values" after saving
* Fixed wrong naming on package deploy
* Removed vis.js from autostart
* Async external connectors fetch
* Speed-up DG start
* FilesView: fix file appearance bug in the card view
* Introduced Entity.package, bugfixes
* Fixed possible NPE
* JS API: BigInt support
* closes #857 Deselection in the scatter plot
* Fixed author parsing
* Package content validation: improve the up-to-date package bundle check
* Fixed tableInput icons
* closes #859 3D scatter plot: implement saving to .png
* Updated beta users
* Tested renderer selection by multiple tags (works fine)
* (Bug) Integration tests: Fix all northwinds for 'External provider: test all Northwinds' test
* refactor: configuration files
* Docs: fix documentation links
* Public submodule
* Add delete FuncCall api
* PackagesView installed filter
* CSS harmonization
* Introduced EntityType.isPackageEntity field
* Fixes #923: Grid: Order or hide columns: provide "move to the top" / "move to the bottom" buttons
* Viewer-specific tooltips (WIP)
* Min() and Max() are not compatible with date and datetime columns #689
* Fixed tests
* removed debug lines
* #closes 858 Cell as an input for widget
* (Bug) FilesView: dock manager exception on opening
* Setup fixes
* Update doc for proxy
* Better transaction management
* Added debug lines
* Some packages not publishing fix
* Removed debug printout
* Add save as png to all viewers
* (Bug) Connection is cloned on sharing
* Reverted GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Dialog: fixed a bug with the resizing
* Grid: Fixed an issue with bool cell renderer
* dapi.files: fixed a bug by reverting to the lazy initialization. Code cleanup.
* fixed: ColumnsInput: an exception when checkedColumns is not specified
* js-api update
* (Bug) Server push is sent twice
* (Bug) Row.setValues doesn't fire onDataChanged event
* Removed debug lines
* Add exceljs to sources
* DateTime scalar functions harmonization
* Fixed shell test
* PowerPack: build fix
* Input fixes, and new disabled styles
* Fixed scalar functions null behavior
* Fixed logical functions in Chrome
* Introduce DateInput with dayjs support
* Fixed deploy exception
* Fixed label css
* (Bug) ScatterPlot: Top and Bottom legend position cuts labels
* (Bug) Selection doesn't work when jitter activated
* ML | Random Data doesn't work ML | Missing Values Imputation doesn't work
* CSS fix
* Fixed dapi and ml tests
* Possible NPE fixed
* Scatterplot 3d detach method
* dapi.root
* View close fixed
* Updated public
* Test: setUpAll and tearDownAll with groups for integration tests
* Test: Skip tests which fail everytime
* Added JS wrapper for EntityRecord
* Packed categories zoom by filter
* package_id migration tuning
* Fixed user login error
* (Bug) FilesUploader: ?layout parameter doesn't work
* CI: Sofiia Podolskaia in beta users
* Test: Use custom session token for integration tests Load session token for integration tests from environment
  variable GROK_SESSION_TOKEN As a fallback value the old session token unit_test_token will be used
* Refactor: Complete refactoring of deploy folder Deprecate old scripts
* Test: File with environment variables for local tests
* GH-954: expanded accessor for TreeViewGroup
* (Bug) Integration tests: External provider (WIP)
* Docs: Instructions to run datlas tests locally on local custom stand
* Test: Add group for all integration tests Fix integration tests which were failing because of session
* closes #936 pie chart improvments
* Test: Add testDevKey config option for docker stands
* Update background for macromolecules
* Test: fix Missing Values Imputation: k-Nearest Neighbour Imputation
* Fixed Datlas tests
* Reverted incorrect fixes
* GH-965: TreeViewGroup.group now also accepts Element
* closes #961 gray out unused categories instead of not showing them
* (Bug) ScatterPlot: Categories rendering issue
* closes #962 Save as PNG: not all elements are added to image
* closes 790 Reset filter doesn't work on "Filter out missing values"
* add button style
* closes #963 Continuous selection on shift+click
* CI: Add short commit to docker image
* Test: Fix Datlas tests
* Load package before loading tests
* Balloons: add the "X" icon
* closes interfaces generation
* (Bug) When user shares an entity from package package is shared too
* Chem harmonization
* closes 935 treemap improvments
* Chem: fix fingerprints
* PC Plot: collaborative filtering
* Test: Increase timeout for datlas tests
* (Bug) Packages permission check doesn't work
* closes 935 with all rowSource modes
* closes 945 Line chart: ability to hide and make semi-transparent lines
* FilesView: add the "Add new share" ribbon icon and link below the tree
* FilesView: 'refresh' icon
* CI: Docker Buildx services auto list
* CI: Docker Demo Chembl DB
* CI: Docker Buildx target
* CI: Create test_numeric_types table for external_provider_test
* Test: Increase timeout for Dataset Client: Projects: CRUD: Remove
* Test: Increase timeout for Datlas tests
* Integration tests: datetime test
* Line chart: add tooltips in multiaxis mode
* New FontAwesome version
* Fixed css
* Fix bakcground styles, for macromolecules page
* Chem: restoring sim search, add indeces
* code fix for #935 treemap improvments
* Pipe file downloads
* Folder sharing
* Ask credentials on file share open
* Minor style changes
* Sketch View: add the reset button
* JS API: Grid.sortByColumns, Grid.sortTypes
* Chem: fix for similarity search
* Grok connect: performance test (select 1)
* Integration tests: refactor external_provider_perf_test
* Fixed AD login
* GH-955: subDir and share methods
* GH-955: js api update
* (Bug) Files: Unable to open context menu for files in explorer view
* Credentials Manager
* (Bug) connection.shares.connection == null condition doesn't work
* BarChart: support additive functions 'unique', 'missing value count' and 'value count' in stacked mode (WIP)
* Docker: Datagrok ca-certficates
* Fixed files API
* CI: Add Alexey Chopovsky to beta users
* Docker: Remove elastic from datagrok
* Docs: Upgrade docker instructions with certificates
* PowerGrid: bumped up version to 1.1.0
* #969 forrula-lines vert lines
* Reversed wrong commit
* Public token
* Chem #995: clean up -Similarity
* Chem #995: open phacts clean up
* Enable saving multiple credentials in credential manager
* Fixed logout
* Docs: Alllow askalkin@datagrok.ai to view encrypted files
* (Bug) JS: grok.functions.call doesn't handle error
* JS API: correct return type for viewer.getOptions
* #953: Filters turned off status is not saved in layout
* #966: Deleted column filter should not be filtered on
* #969 new classes system ,aybe should be reworked
* Change order of interop type checks
* #969 linechart hittest
* closes #996 Tooltip not working for partly truncated values
* Funcs: JS-based context functions for multiple columns
* (Bug) ColumnsInput: number of selected columns is not shown upon construction
* JS API: GridColumn.move method
* closes #1015: PowerGrid: Sparklines: Ability to add sparklines from the "Actions" pane for selected columns
* JS API: GridColumn.scrollIntoView, RangeSlider.scrollBy, RangeSlider.scrollTo
* PowerGrid: Sparklines: accessing settings from the hamburger menu
* #930: PowerGrid: sparklines: ability to rename and remove summary columns
* #930: PowerGrid: sparklines: managing column order
* #930: PowerGrid: sparklines: setting the sparkline column as a current object after it is created
* Possible NPE fix
* Deprecated old tests
* Chem: #454 descriptors
* Chem: small fix
* Fixed package init
* Datlas: fixed project info retrieval
* #1000 viewport exists
* #1000 viewport apeears properly, style changed
* dialog cliptoscreen onresize
* #1000 bands
* Remove labs submodule
* (Bug) Custom context menu items for JS widgets and viewers cannot be registered via API
* (Bug) Dropbox provider doesn't work
* Chem: fix descriptors test
* Chem: #995 removing r-groups
* Chem: #454 compute side polishing
* Update API and public token
* (Bug) Proxy to docker container doesn't work
* Wiki: replace links
* (Bug) Package Credentials Editor doesn't work
* JS API: code cleanup and minor refactoring
* Histogram: reduced the marker size to 4 pixels
* Better wording
* PowerPack: testing \- WIP
* closes #1022 Tooltip: do not show value while using custom renderers
* closes 1021 matrix plot from layout
* JS API: Point.distanceTo, Rect.containsPoint
* dialog resizing before dom created fixwd
* Update beta users list
* (Bug) grok.shell.registerViewer fails
* Minor code cleanup
* (Bug) Ability to clone the connection (WIP)
* Grok Spawner: more information for every request in return
* Docs: devops docs
* closes #1001 matrix plot inner viewer settings
* JS API: update param type in Column.aggregate
* Grok Spawner: docker entrypoint script
* Chem: rgroups exclusion
* #1001 syncing tables in wrong place
* PowerGrid: version bump up
* Sketch View: update UI
* closes #657 line chart selectors with special symbols
* (Bug) Substructure Filter: Disabled if you drag the column to the filter
* closes #856 Chem: add molecule filter to a filter-group by default
* Packages: disable certain sources checks for non-webpack packages
* Chem: descriptors fixes
* (Bug) Method click() throws an exception when called for top menu item
* (Bug) JS: disable Uint8List convertion
* Implement a mechanism to test functions
* Hierarchical pipe-delimited info panels (example: Chembl | Substructure)
* Fixed #1029: File Manager: clicking the file in the Folder Tree doesn't trigger file preview in the main view
* JS: DataFrame.toByteArray method
* Harmonize membership editor
* Package content validation: check package properties' types
* (Bug) Scripts: Sample table does not open
* (Bug) __jsThen error during some packages init
* Docking: highlighting the separator
* Add a script to find internal links in help
* Add the function namespace to tests
* Grid: automatic format for small floats
* Chem: api fix
* (Bug) Grok Spawner: After image creation the container does not start
* CI: Add time spent to junit convertion script
* (Bug) Make loadKeyPair method more robust
* Package permissions fix
* Updated beta users list
* Merge 'master'
* CI: add jq to puppeteer image
* Dockerfiles: fixed container address retrieval
* update img
* Update img
* Docs: onboarding docs refactoring
* Grok Spawner: publish ports on local docker
* Test: mute failing tests
* Refactor Secret Manager connection
* (Bug) Core: MapInput cannot be added to a dialog for function parameters
* Grok Spawner: network properties in responses
* CI: Junit test results converter
* 1002 adds axes to matrix plot
* #1002 inner viewwer setting icon not merged
* Improve start mode selection
* Dockerfiles: get container port
* ##1042: Onboarding documentation: moved things around, minor cleanup
* CI: Add tests manually to junit reports
* Dashboard Projects (WIP)
* #1002 matrix plot axes code refactored
* CI: Convert JSON files to Junit reports
* CI: Add Rsync to pepputeer base image
* PowerGrid: version 1.1.2
* CI: Correct messages to jest tests failures
* Provided a hook to borrow Hamburger menu by external components
* #943 basic axes
* Removed grid column reset for Hamburger menu
* CI: Append timeout and evaluation errors to junit reports
* closes #953 Filters turned off status is not saved in layout
* CI: Include multiline stdout to Junit reports
* CI: Add author to Junit tests output
* JKG-Add-libmamba-to-conda-resolver-mechism
* Sketch View: remove the icon for adding viewers when the view is opened to edit a tooltip
* Docs: Access to servers using SSH
* Added exportable indices bounds for visible rows and columns
* Added xangle registration for export of visible rows and columns
* #943 3D scatter plot enhancements (log coordinates, some ui moments)
* Add new method (RangeSlider_SetShowHandles) \- need for on/off resize in slider
* fixes choice unput value not screened
* #closes 1075
* PropertyGrid editors: moved each editor to a separate file.
* Viewers: ability to pre-aggregate source DataFrame (WIP)
* (Bug) Scripting: JS Script fails, when contains more than 6 params
* (Bug) Can't share a folder from Home connection
* (Bug) Project duplicates Data Query name
* PC Plot: normalization mode
* PC Plot: auto-hide unused filters
* PC Plot: Show value axis when the scale is global
* RangeSlider: an option to not show the bar when the slider is expanded
* closes #1002 axes for trellis plot
* closes #1087 Use custom renderers with scatterplot label
* Core: add html2canvas library
* closes #1072 when context menu and column selector popup menus overlap, the click event triggers on the wrong menu
  item
* dpetrov/GitHub-1083 Made Actions closable in the Hamburger menu
* ColumnComboBox: not opening the drop-down part when right-clicked
* #943 use existing methods for axes tickmarks
* labels refactoring
* Dockerfile methods WIP
* Dockerfiles containerAddress fix
* DockerfilesClient proxyRequest method fix
* GH-863: Dockerfiles JS API update
* #1096 x axis fixed
* #1096 y labels fix
* Proper packages names
* (Bug) Package Entities Naming issues (WIP)
* (Bug) Script icon incorrect
* dpetrov/GitHub-1093 Added virtual columns support for grid's columns manager
* #1096 fixes multiple columns labels
* dpetrov/GitHub-1093 Minor fixes for the grid's columns manager
* Jenkins: Docs: Add Jenkins setup information
* Jenkins: Docs: Add password information for Jenkins
* (Bug) Query Transformation name issue
* #1098 add hittesting and selection to pc plot
* JS API: Custom inputs (WIP)
* dpetrov/GitHub-1104: Fixed returning null view field in ViewLayout object
* Merged in onuf/color-coding-update-bug (pull request #270)
* Group passwords (WIP)
* #601 PC plot improvements(add interactivity)
* #928 context menus
* 601 PC plot improvements(pc plot dies when transformation reset)
* #943 axes improvments
* Merged in onuf/color-coding-update-bug (pull request #273)
* update public token
* #943 issue with axes on otherinstances
* add Dmytro Nahovskyi to list
* (Bug) Right Click on the Grid's Row Header generates exception
* 1103 fixes filters in layout behavior
* change grok conneect version 1.5
* Update beta users table
* Charts: Radar: clarified property description
* Fix button styles
* Ability to load entitty from inactive package version
* MultiColumn Selector \- use subsets of columns when needed  (WIP)
* Error when serealizing df with virtual column
* Viewers: "General | Open in workspace" option
* Databases: groom the providers (WIP)
* User can't log in is login contains hyphens
* (Bug) Pinned Columns: row number column appears twice after layout application
* requirements.R default value for plsdepot has been returned
* buildx_docker_image.sh returned default value for cache_push (true)
* deploy/demo_db_clickhouse/demo_db_clickhouse:
  -f772192786c7f01c4f7afacfe68affa1dd2fe9a9-f772192786c7f01c4f7afacfe68affa1dd2fe9a9.docker.tar edited online with
  Bitbucket
* closes 1126 scatterPlot packed cats with no values
* #609 pc plot filter in aggr mode
* skip package filter test
* Simplified icon methods
* #943 3dscatterplot onload error sometimes
* (Bug) ColorPicker blinks in Chrome (WIP)
* Chem: remove substructure search from panel
* closes 1110 Viewers: legend cannot be resized again
* closes #1131 MultiColumn Selector: selects only columns filtered by search
* (Bug) Nginx: Add additional logging to investigate timeouts
* Docs: LE certificates remewal
* #1103 layout not saving
* (Bug) Grok Connect: the connection visibility
* Settings: mark instance as production
* #943 3dscatterplot fix try
* Button alignment for form
* #609 aggregation mode for pc plot
* closes #1134 trellis plot viewer icon
* Introduce Viewer Serialization Context for saving layouts and projects
* (Bug) InfoPanels: Returning map or widget doesn't work
* Implement Logger (WIP)
* (Bug) Packages: Debug version duplicates release
* Careers section updates
* Fixed ORM test
* column manager checkwrong columns after reordering
* closes #776 Viewer legend is not synchronized with in-viewer filter
* update deploy.md
* JS: Icons support for JS Viewers
* #943 scrollbars should not appear
* closes 932 Vertical scrollbar is not working in Order or Hide Columns menu
* Fixed ambiguous column names
* Package functions: set helpUrl
* Removed debug code
* ORM: FuncCall visibility should depend on Func visibility (WIP)
* #1110 legend cats count incorrect
* add images for connector

## 2022-10-19 Dev build 1.8.0

### Addressed Issues

* com
* Fix default value for select (choice) (remove '' and "" from string when chosen default value)
* Implement setRawData for bool columns
* render mols vertically on scatterplot x axis
* lrellis plot not working in LS
* Fixed the method for default choices
* delete useless prints
* bar chart slider positioning
* Grok Spawner: fix default value for GrokSpawnerAddress
* Make modals resizable
* Legends Resizing
* closes #810 Reordered selected columns
* add started/finished/author to funccall
* Layouts | Save to gallery: newly saved layout should appear in the property panel
* Layouts pane: ability to put the layout on the property panel
* CI: Deploy process documentation
* closes #798 selecting dots in jittered scatterplot
* #828 formula lines modal fixed
* #828 order/hide dialog rewworked
* #828 make modal not resizable by default
* JS API: ui.fileBrowser
* Work in progress.
* #831: HitTriage: work in progress
* Moved InputBase descendants to separate files.
* AppEvents.onInputCreated global event
* Ability to reuse database browser, and specify how results are interpreted
* Table Input: ability to get table by executing a database query
* Docs: Publish version from dev to Docker Hub
* Datlas: queue interaction for Grok Connect (WIP)
* Docker image for packages tests
* add compact view to scatterplot 1.0
* Scatterplot compact view
* JS API: View.statusBarPanels (WIP)
* (Bug) Grok connect: decimal, real, and double tests failed
* Grok connect: refactoring of numeric type tests
* JS: UI Forms harmonization (WIP)
* (Bug) Failed to delete package
* DB Indexes
* Fixed filename
* Fixed CSS
* Added a demo dataset with image urls
* PowerGrid: help on linked images
* JS UI Test: Scripting in Add new column (WIP)
* Functions: show function class instead of null in status bar
* Grok connect: refactoring of numeric type tests
* DataFrame.setTag(...) now returns this.
* Renderers: ability to specify column tags as conditions
* adds menu click method
* modal click simulation method adjustment
* Datlas connections leaking investigation
* Datlas: remove trailing ; in the row checking query
* JS API: Menu.click method
* #831: HitTriage: work in progress \- refactoring
* closes #800 Form does not work after filtering
* Introduced Entity.package Implemented safe delete New permissions check Derived properties support in ORM
* Introduced Entity.package
* GH-686: TreeView events
* modal placement, dont leave page on backspace
* GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Fixed deadlock in exchangeDevKey method
* CSS fixes
* Introduced Entity.package, reviewed Entities list
* Introduced Entity.package, introduced PackageEntityMixin
* Autogenerated code
* Speed-up FilesView start
* Fixed analyzer warnings
* Fixed package.functions section exception
* closes #791 Can not uncheck "Filter out missing values" after saving
* Fixed wrong naming on package deploy
* Removed vis.js from autostart
* Async external connectors fetch
* Speed-up DG start
* FilesView: fix file appearance bug in the card view
* Introduced Entity.package, bugfixes
* Fixed possible NPE
* JS API: BigInt support
* closes #857 Deselection in the scatter plot
* Fixed author parsing
* Package content validation: improve the up-to-date package bundle check
* Fixed tableInput icons
* closes #859 3D scatter plot: implement saving to .png
* Updated beta users
* Tested renderer selection by multiple tags (works fine)
* (Bug) Integration tests: Fix all northwinds for 'External provider: test all Northwinds' test
* refactor: configuration files
* Docs: fix documentation links
* Public submodule
* Add delete FuncCall api
* PackagesView installed filter
* CSS harmonization
* Introduced EntityType.isPackageEntity field
* Fixes #923: Grid: Order or hide columns: provide "move to the top" / "move to the bottom" buttons
* Viewer-specific tooltips (WIP)
* Min() and Max() are not compatible with date and datetime columns #689
* Fixed tests
* removed debug lines
* #closes 858 Cell as an input for widget
* (Bug) FilesView: dock manager exception on opening
* Setup fixes
* Update doc for proxy
* Better transaction management
* Added debug lines
* Some packages not publishing fix
* Removed debug printout
* Add save as png to all viewers
* (Bug) Connection is cloned on sharing
* Reverted GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Dialog: fixed a bug with the resizing
* Grid: Fixed an issue with bool cell renderer
* dapi.files: fixed a bug by reverting to the lazy initialization. Code cleanup.
* fixed: ColumnsInput: an exception when checkedColumns is not specified
* js-api update
* (Bug) Server push is sent twice
* (Bug) Row.setValues doesn't fire onDataChanged event
* Removed debug lines
* Add exceljs to sources
* DateTime scalar functions harmonization
* Fixed shell test
* PowerPack: build fix
* Input fixes, and new disabled styles
* Fixed scalar functions null behavior
* Fixed logical functions in Chrome
* Introduce DateInput with dayjs support
* Fixed deploy exception
* Fixed label css
* (Bug) ScatterPlot: Top and Bottom legend position cuts labels
* (Bug) Selection doesn't work when jitter activated
* ML | Random Data doesn't work ML | Missing Values Imputation doesn't work
* CSS fix
* Fixed dapi and ml tests
* Possible NPE fixed
* Scatterplot 3d detach method
* dapi.root
* View close fixed
* Updated public
* Test: setUpAll and tearDownAll with groups for integration tests
* Test: Skip tests which fail everytime
* Added JS wrapper for EntityRecord
* Packed categories zoom by filter
* package_id migration tuning
* Fixed user login error
* (Bug) FilesUploader: ?layout parameter doesn't work
* CI: Sofiia Podolskaia in beta users
* Test: Use custom session token for integration tests Load session token for integration tests from environment
  variable GROK_SESSION_TOKEN As a fallback value the old session token unit_test_token will be used
* Refactor: Complete refactoring of deploy folder Deprecate old scripts
* Test: File with environment variables for local tests
* GH-954: expanded accessor for TreeViewGroup
* (Bug) Integration tests: External provider (WIP)
* Docs: Instructions to run datlas tests locally on local custom stand
* Test: Add group for all integration tests Fix integration tests which were failing because of session
* closes #936 pie chart improvments
* Test: Add testDevKey config option for docker stands
* Update background for macromolecules
* Test: fix Missing Values Imputation: k-Nearest Neighbour Imputation
* Fixed Datlas tests
* Reverted incorrect fixes
* GH-965: TreeViewGroup.group now also accepts Element
* closes #961 gray out unused categories instead of not showing them
* (Bug) ScatterPlot: Categories rendering issue
* closes #962 Save as PNG: not all elements are added to image
* closes 790 Reset filter doesn't work on "Filter out missing values"
* add button style
* closes #963 Continuous selection on shift+click
* CI: Add short commit to docker image
* Test: Fix Datlas tests
* Load package before loading tests
* Balloons: add the "X" icon
* closes interfaces generation
* (Bug) When user shares an entity from package package is shared too
* Chem harmonization
* closes 935 treemap improvments
* Chem: fix fingerprints
* PC Plot: collaborative filtering
* Test: Increase timeout for datlas tests
* (Bug) Packages permission check doesn't work
* closes 935 with all rowSource modes
* closes 945 Line chart: ability to hide and make semi-transparent lines
* FilesView: add the "Add new share" ribbon icon and link below the tree
* FilesView: 'refresh' icon
* CI: Docker Buildx services auto list
* CI: Docker Demo Chembl DB
* CI: Docker Buildx target
* CI: Create test_numeric_types table for external_provider_test
* Test: Increase timeout for Dataset Client: Projects: CRUD: Remove
* Test: Increase timeout for Datlas tests
* Integration tests: datetime test
* Line chart: add tooltips in multiaxis mode
* New FontAwesome version
* Fixed css
* Fix bakcground styles, for macromolecules page
* Chem: restoring sim search, add indeces
* code fix for #935 treemap improvments
* Pipe file downloads
* Folder sharing
* Ask credentials on file share open
* Minor style changes
* Sketch View: add the reset button
* JS API: Grid.sortByColumns, Grid.sortTypes
* Chem: fix for similarity search
* Grok connect: performance test (select 1)
* Integration tests: refactor external_provider_perf_test
* Fixed AD login
* GH-955: subDir and share methods
* GH-955: js api update
* (Bug) Files: Unable to open context menu for files in explorer view
* Credentials Manager
* (Bug) connection.shares.connection == null condition doesn't work
* BarChart: support additive functions 'unique', 'missing value count' and 'value count' in stacked mode (WIP)
* Docker: Datagrok ca-certficates
* Fixed files API
* CI: Add Alexey Chopovsky to beta users
* Docker: Remove elastic from datagrok
* Docs: Upgrade docker instructions with certificates
* PowerGrid: bumped up version to 1.1.0
* #969 forrula-lines vert lines
* Reversed wrong commit
* Public token
* Chem #995: clean up -Similarity
* Chem #995: open phacts clean up
* Enable saving multiple credentials in credential manager
* Fixed logout
* Docs: Alllow askalkin@datagrok.ai to view encrypted files
* (Bug) JS: grok.functions.call doesn't handle error
* JS API: correct return type for viewer.getOptions
* #953: Filters turned off status is not saved in layout
* #966: Deleted column filter should not be filtered on
* #969 new classes system ,aybe should be reworked
* Change order of interop type checks
* #969 linechart hittest
* closes #996 Tooltip not working for partly truncated values
* Funcs: JS-based context functions for multiple columns
* (Bug) ColumnsInput: number of selected columns is not shown upon construction
* JS API: GridColumn.move method
* closes #1015: PowerGrid: Sparklines: Ability to add sparklines from the "Actions" pane for selected columns
* JS API: GridColumn.scrollIntoView, RangeSlider.scrollBy, RangeSlider.scrollTo
* PowerGrid: Sparklines: accessing settings from the hamburger menu
* #930: PowerGrid: sparklines: ability to rename and remove summary columns
* #930: PowerGrid: sparklines: managing column order
* #930: PowerGrid: sparklines: setting the sparkline column as a current object after it is created
* Possible NPE fix
* Deprecated old tests
* Chem: #454 descriptors
* Chem: small fix
* Fixed package init
* Datlas: fixed project info retrieval
* #1000 viewport exists
* #1000 viewport apeears properly, style changed
* dialog cliptoscreen onresize
* #1000 bands
* Remove labs submodule
* (Bug) Custom context menu items for JS widgets and viewers cannot be registered via API
* (Bug) Dropbox provider doesn't work
* Chem: fix descriptors test
* Chem: #995 removing r-groups
* Chem: #454 compute side polishing
* Update API and public token
* (Bug) Proxy to docker container doesn't work (WIP)
* Wiki: replace links
* (Bug) Package Credentials Editor doesn't work
* JS API: code cleanup and minor refactoring
* Histogram: reduced the marker size to 4 pixels
* Better wording
* PowerPack: testing \- WIP
* closes #1022 Tooltip: do not show value while using custom renderers
* closes 1021 matrix plot from layout
* JS API: Point.distanceTo, Rect.containsPoint
* dialog resizing before dom created fixwd
* Update beta users list
* (Bug) grok.shell.registerViewer fails
* Minor code cleanup
* (Bug) Ability to clone the connection (WIP)
* Grok Spawner: more information for every request in return
* Docs: devops docs
* closes #1001 matrix plot inner viewer settings
* JS API: update param type in Column.aggregate
* Grok Spawner: docker entrypoint script
* Chem: rgroups exclusion
* #1001 syncing tables in wrong place
* PowerGrid: version bump up
* Sketch View: update UI
* closes #657 line chart selectors with special symbols
* (Bug) Substructure Filter: Disabled if you drag the column to the filter
* closes #856 Chem: add molecule filter to a filter-group by default
* Packages: disable certain sources checks for non-webpack packages
* Chem: descriptors fixes
* (Bug) Method click() throws an exception when called for top menu item
* (Bug) JS: disable Uint8List convertion
* Implement a mechanism to test functions
* Hierarchical pipe-delimited info panels (example: Chembl | Substructure)
* Fixed #1029: File Manager: clicking the file in the Folder Tree doesn't trigger file preview in the main view
* JS: DataFrame.toByteArray method
* Harmonize membership editor
* Package content validation: check package properties' types
* (Bug) Scripts: Sample table does not open
* (Bug) __jsThen error during some packages init
* Docking: highlighting the separator
* Add a script to find internal links in help
* Add the function namespace to tests
* Grid: automatic format for small floats
* Chem: api fix
* (Bug) Grok Spawner: After image creation the container does not start
* CI: Add time spent to junit convertion script
* (Bug) Make loadKeyPair method more robust
* Package permissions fix
* Updated beta users list
* Merge 'master'
* CI: add jq to puppeteer image
* Dockerfiles: fixed container address retrieval
* update img
* Update img
* Docs: onboarding docs refactoring
* Grok Spawner: publish ports on local docker
* Test: mute failing tests
* Refactor Secret Manager connection
* (Bug) Core: MapInput cannot be added to a dialog for function parameters
* Grok Spawner: network properties in responses
* CI: Junit test results converter
* 1002 adds axes to matrix plot
* #1002 inner viewwer setting icon not merged
* Improve start mode selection
* Dockerfiles: get container port (WIP)
* ##1042: Onboarding documentation: moved things around, minor cleanup
* CI: Add tests manually to junit reports
* Dashboard Projects (WIP)
* #1002 matrix plot axes code refactored
* CI: Convert JSON files to Junit reports
* CI: Add Rsync to pepputeer base image
* PowerGrid: version 1.1.2
* CI: Correct messages to jest tests failures
* Provided a hook to borrow Hamburger menu by external components
* #943 basic axes
* Removed grid column reset for Hamburger menu
* CI: Append timeout and evaluation errors to junit reports
* closes #953 Filters turned off status is not saved in layout
* CI: Include multiline stdout to Junit reports
* CI: Add author to Junit tests output
* JKG-Add-libmamba-to-conda-resolver-mechism
* Sketch View: remove the icon for adding viewers when the view is opened to edit a tooltip
* Docs: Access to servers using SSH
* Added exportable indices bounds for visible rows and columns
* Added xangle registration for export of visible rows and columns
* #943 3D scatter plot enhancements (log coordinates, some ui moments)
* Add new method (RangeSlider_SetShowHandles) \- need for on/off resize in slider
* fixes choice unput value not screened
* #closes 1075
* PropertyGrid editors: moved each editor to a separate file.
* Viewers: ability to pre-aggregate source DataFrame (WIP)
* (Bug) Scripting: JS Script fails, when contains more than 6 params
* (Bug) Can't share a folder from Home connection
* (Bug) Project duplicates Data Query name
* PC Plot: normalization mode
* PC Plot: auto-hide unused filters
* PC Plot: Show value axis when the scale is global
* RangeSlider: an option to not show the bar when the slider is expanded
* closes #1002 axes for trellis plot
* closes #1087 Use custom renderers with scatterplot label
* Core: add html2canvas library
* closes #1072 when context menu and column selector popup menus overlap, the click event triggers on the wrong menu
  item
* dpetrov/GitHub-1083 Made Actions closable in the Hamburger menu
* ColumnComboBox: not opening the drop-down part when right-clicked
* #943 use existing methods for axes tickmarks
* labels refactoring
* Dockerfile methods WIP
* Dockerfiles containerAddress fix
* DockerfilesClient proxyRequest method fix
* GH-863: Dockerfiles JS API update
* #1096 x axis fixed
* #1096 y labels fix
* Proper packages names
* (Bug) Package Entities Naming issues (WIP)
* (Bug) Script icon incorrect
* dpetrov/GitHub-1093 Added virtual columns support for grid's columns manager
* #1096 fixes multiple columns labels
* dpetrov/GitHub-1093 Minor fixes for the grid's columns manager
* Jenkins: Docs: Add Jenkins setup information
* Jenkins: Docs: Add password information for Jenkins
* (Bug) Query Transformation name issue
* #1098 add hittesting and selection to pc plot
* JS API: Custom inputs (WIP)
* dpetrov/GitHub-1104: Fixed returning null view field in ViewLayout object
* Merged in onuf/color-coding-update-bug (pull request #270)
* Group passwords (WIP)
* #601 PC plot improvements(add interactivity)
* #928 context menus
* 601 PC plot improvements(pc plot dies when transformation reset)
* #943 axes improvments
* Merged in onuf/color-coding-update-bug (pull request #273)
* update public token
* #943 issue with axes on otherinstances
* add Dmytro Nahovskyi to list
* (Bug) Right Click on the Grid's Row Header generates exception (WIP)
* 1103 fixes filters in layout behavior
* change grok conneect version 1.5
* Update beta users table
* Charts: Radar: clarified property description
* Fix button styles
* Ability to load entitty from inactive package version
* MultiColumn Selector \- use subsets of columns when needed  (WIP)
* Error when serealizing df with virtual column
* Viewers: "General | Open in workspace" option
* Databases: groom the providers (WIP)
* User can't log in is login contains hyphens
* (Bug) Pinned Columns: row number column appears twice after layout application (WIP)
* requirements.R default value for plsdepot has been returned
* buildx_docker_image.sh returned default value for cache_push (true)
* deploy/demo_db_clickhouse/demo_db_clickhouse:
  -f772192786c7f01c4f7afacfe68affa1dd2fe9a9-f772192786c7f01c4f7afacfe68affa1dd2fe9a9.docker.tar edited online with
  Bitbucket
* closes 1126 scatterPlot packed cats with no values
* #609 pc plot filter in aggr mode
* skip package filter test
* (Bug) ColorPicker blinks in Chrome (WIP)

## 2022-10-06 Dev build 1.6.12

### Addressed Issues

* closes #953 Filters turned off status is not saved in layout

## 2022-09-15 Dev build 1.6.11

### Addressed Issues

* com
* Fix default value for select (choice) (remove '' and "" from string when chosen default value)
* Implement setRawData for bool columns
* render mols vertically on scatterplot x axis
* lrellis plot not working in LS
* Fixed the method for default choices
* delete useless prints
* bar chart slider positioning
* Grok Spawner: fix default value for GrokSpawnerAddress
* Make modals resizable
* Legends Resizing
* closes #810 Reordered selected columns
* add started/finished/author to funccall
* Layouts | Save to gallery: newly saved layout should appear in the property panel
* Layouts pane: ability to put the layout on the property panel
* CI: Deploy process documentation
* closes #798 selecting dots in jittered scatterplot
* #828 formula lines modal fixed
* #828 order/hide dialog rewworked
* #828 make modal not resizable by default
* JS API: ui.fileBrowser
* Work in progress.
* #831: HitTriage: work in progress
* Moved InputBase descendants to separate files.
* AppEvents.onInputCreated global event
* Ability to reuse database browser, and specify how results are interpreted
* Table Input: ability to get table by executing a database query
* Docs: Publish version from dev to Docker Hub
* Datlas: queue interaction for Grok Connect (WIP)
* Docker image for packages tests
* add compact view to scatterplot 1.0
* Scatterplot compact view
* JS API: View.statusBarPanels (WIP)
* (Bug) Grok connect: decimal, real, and double tests failed
* Grok connect: refactoring of numeric type tests
* JS: UI Forms harmonization (WIP)
* (Bug) Failed to delete package
* DB Indexes
* Fixed filename
* Fixed CSS
* Added a demo dataset with image urls
* PowerGrid: help on linked images
* JS UI Test: Scripting in Add new column (WIP)
* Functions: show function class instead of null in status bar
* Grok connect: refactoring of numeric type tests
* DataFrame.setTag(...) now returns this.
* Renderers: ability to specify column tags as conditions
* adds menu click method
* modal click simulation method adjustment
* Datlas connections leaking investigation (WIP)
* Datlas: remove trailing ; in the row checking query
* JS API: Menu.click method
* #831: HitTriage: work in progress \- refactoring
* closes #800 Form does not work after filtering
* Introduced Entity.package Implemented safe delete New permissions check Derived properties support in ORM
* Introduced Entity.package
* GH-686: TreeView events
* modal placement, dont leave page on backspace
* GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Fixed deadlock in exchangeDevKey method
* CSS fixes
* Introduced Entity.package, reviewed Entities list
* Introduced Entity.package, introduced PackageEntityMixin
* Autogenerated code
* Speed-up FilesView start
* Fixed analyzer warnings
* Fixed package.functions section exception
* closes #791 Can not uncheck "Filter out missing values" after saving
* Fixed wrong naming on package deploy
* Removed vis.js from autostart
* Async external connectors fetch
* Speed-up DG start
* FilesView: fix file appearance bug in the card view
* Introduced Entity.package, bugfixes
* Fixed possible NPE
* JS API: BigInt support
* closes #857 Deselection in the scatter plot
* Fixed author parsing
* Package content validation: improve the up-to-date package bundle check
* Fixed tableInput icons
* closes #859 3D scatter plot: implement saving to .png
* Updated beta users
* Tested renderer selection by multiple tags (works fine)
* (Bug) Integration tests: Fix all northwinds for 'External provider: test all Northwinds' test
* refactor: configuration files
* Docs: fix documentation links
* Public submodule
* Add delete FuncCall api
* PackagesView installed filter
* CSS harmonization
* Introduced EntityType.isPackageEntity field
* Fixes #923: Grid: Order or hide columns: provide "move to the top" / "move to the bottom" buttons
* Viewer-specific tooltips (WIP)
* Min() and Max() are not compatible with date and datetime columns #689
* Fixed tests
* removed debug lines
* #closes 858 Cell as an input for widget
* (Bug) FilesView: dock manager exception on opening
* Setup fixes
* Update doc for proxy
* Better transaction management
* Added debug lines
* Some packages not publishing fix
* Removed debug printout
* Add save as png to all viewers
* (Bug) Connection is cloned on sharing
* Reverted GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Dialog: fixed a bug with the resizing
* Grid: Fixed an issue with bool cell renderer
* dapi.files: fixed a bug by reverting to the lazy initialization. Code cleanup.
* fixed: ColumnsInput: an exception when checkedColumns is not specified
* js-api update
* (Bug) Server push is sent twice
* (Bug) Row.setValues doesn't fire onDataChanged event
* Removed debug lines
* Add exceljs to sources
* DateTime scalar functions harmonization
* Fixed shell test
* PowerPack: build fix
* Input fixes, and new disabled styles
* Fixed scalar functions null behavior
* Fixed logical functions in Chrome
* Introduce DateInput with dayjs support
* Fixed deploy exception
* Fixed label css
* (Bug) ScatterPlot: Top and Bottom legend position cuts labels
* (Bug) Selection doesn't work when jitter activated
* ML | Random Data doesn't work ML | Missing Values Imputation doesn't work
* CSS fix
* Fixed dapi and ml tests
* Possible NPE fixed
* Scatterplot 3d detach method
* dapi.root
* View close fixed
* Updated public
* Test: setUpAll and tearDownAll with groups for integration tests
* Test: Skip tests which fail everytime
* Added JS wrapper for EntityRecord
* Packed categories zoom by filter
* package_id migration tuning
* Fixed user login error
* (Bug) FilesUploader: ?layout parameter doesn't work
* CI: Sofiia Podolskaia in beta users
* Test: Use custom session token for integration tests Load session token for integration tests from environment
  variable GROK_SESSION_TOKEN As a fallback value the old session token unit_test_token will be used
* Refactor: Complete refactoring of deploy folder Deprecate old scripts
* Test: File with environment variables for local tests
* GH-954: expanded accessor for TreeViewGroup
* (Bug) Integration tests: External provider (WIP)
* Docs: Instructions to run datlas tests locally on local custom stand
* Test: Add group for all integration tests Fix integration tests which were failing because of session
* closes #936 pie chart improvments
* Test: Add testDevKey config option for docker stands
* Update background for macromolecules
* Test: fix Missing Values Imputation: k-Nearest Neighbour Imputation
* Fixed Datlas tests
* Reverted incorrect fixes
* GH-965: TreeViewGroup.group now also accepts Element
* closes #961 gray out unused categories instead of not showing them
* (Bug) ScatterPlot: Categories rendering issue
* closes #962 Save as PNG: not all elements are added to image
* closes 790 Reset filter doesn't work on "Filter out missing values"
* add button style
* closes #963 Continuous selection on shift+click
* CI: Add short commit to docker image
* Test: Fix Datlas tests
* Load package before loading tests
* Balloons: add the "X" icon
* closes interfaces generation
* (Bug) When user shares an entity from package package is shared too
* Chem harmonization (WIP)
* closes 935 treemap improvments
* Chem: fix fingerprints
* PC Plot: collaborative filtering
* Test: Increase timeout for datlas tests
* (Bug) Packages permission check doesn't work
* closes 935 with all rowSource modes
* closes 945 Line chart: ability to hide and make semi-transparent lines
* FilesView: add the "Add new share" ribbon icon and link below the tree
* FilesView: 'refresh' icon
* CI: Docker Buildx services auto list
* CI: Docker Demo Chembl DB
* CI: Docker Buildx target
* CI: Create test_numeric_types table for external_provider_test
* Test: Increase timeout for Dataset Client: Projects: CRUD: Remove
* Test: Increase timeout for Datlas tests
* Integration tests: datetime test
* Line chart: add tooltips in multiaxis mode
* New FontAwesome version
* Fixed css
* Fix bakcground styles, for macromolecules page
* Chem: restoring sim search, add indeces
* code fix for #935 treemap improvments
* Pipe file downloads
* Folder sharing
* Ask credentials on file share open
* Minor style changes
* Sketch View: add the reset button
* JS API: Grid.sortByColumns, Grid.sortTypes
* Chem: fix for similarity search
* Grok connect: performance test (select 1)
* Integration tests: refactor external_provider_perf_test
* Fixed AD login
* GH-955: subDir and share methods
* GH-955: js api update
* (Bug) Files: Unable to open context menu for files in explorer view
* Credentials Manager
* (Bug) connection.shares.connection == null condition doesn't work
* BarChart: support additive functions 'unique', 'missing value count' and 'value count' in stacked mode (WIP)
* Docker: Datagrok ca-certficates
* Fixed files API
* CI: Add Alexey Chopovsky to beta users
* Docker: Remove elastic from datagrok
* Docs: Upgrade docker instructions with certificates
* PowerGrid: bumped up version to 1.1.0
* #969 forrula-lines vert lines
* Reversed wrong commit
* Public token
* Chem #995: clean up -Similarity
* Chem #995: open phacts clean up
* Enable saving multiple credentials in credential manager (WIP)
* Fixed logout
* Docs: Alllow askalkin@datagrok.ai to view encrypted files
* (Bug) JS: grok.functions.call doesn't handle error
* JS API: correct return type for viewer.getOptions
* #953: Filters turned off status is not saved in layout
* #966: Deleted column filter should not be filtered on
* #969 new classes system ,aybe should be reworked
* Change order of interop type checks
* #969 linechart hittest
* closes #996 Tooltip not working for partly truncated values
* Funcs: JS-based context functions for multiple columns (WIP)
* (Bug) ColumnsInput: number of selected columns is not shown upon construction
* JS API: GridColumn.move method
* closes #1015: PowerGrid: Sparklines: Ability to add sparklines from the "Actions" pane for selected columns
* JS API: GridColumn.scrollIntoView, RangeSlider.scrollBy, RangeSlider.scrollTo
* PowerGrid: Sparklines: accessing settings from the hamburger menu
* #930: PowerGrid: sparklines: ability to rename and remove summary columns
* #930: PowerGrid: sparklines: managing column order
* #930: PowerGrid: sparklines: setting the sparkline column as a current object after it is created
* Possible NPE fix
* Deprecated old tests
* Chem: #454 descriptors
* Chem: small fix
* Fixed package init
* Datlas: fixed project info retrieval
* #1000 viewport exists
* #1000 viewport apeears properly, style changed
* dialog cliptoscreen onresize
* #1000 bands
* Remove labs submodule
* (Bug) Custom context menu items for JS widgets and viewers cannot be registered via API
* (Bug) Dropbox provider doesn't work
* Chem: fix descriptors test
* Chem: #995 removing r-groups
* Chem: #454 compute side polishing
* Update API and public token
* (Bug) Proxy to docker container doesn't work (WIP)
* Wiki: replace links
* (Bug) Package Credentials Editor doesn't work
* JS API: code cleanup and minor refactoring
* Histogram: reduced the marker size to 4 pixels
* Better wording
* PowerPack: testing \- WIP
* closes #1022 Tooltip: do not show value while using custom renderers
* closes 1021 matrix plot from layout
* JS API: Point.distanceTo, Rect.containsPoint
* dialog resizing before dom created fixwd
* Update beta users list
* (Bug) grok.shell.registerViewer fails
* Minor code cleanup
* (Bug) Ability to clone the connection (WIP)
* Grok Spawner: more information for every request in return
* Docs: devops docs
* closes #1001 matrix plot inner viewer settings
* JS API: update param type in Column.aggregate
* Grok Spawner: docker entrypoint script
* Chem: rgroups exclusion
* #1001 syncing tables in wrong place
* PowerGrid: version bump up
* Sketch View: update UI
* closes #657 line chart selectors with special symbols
* (Bug) Substructure Filter: Disabled if you drag the column to the filter
* closes #856 Chem: add molecule filter to a filter-group by default
* Packages: disable certain sources checks for non-webpack packages
* Chem: descriptors fixes
* (Bug) Method click() throws an exception when called for top menu item
* (Bug) JS: disable Uint8List convertion
* Hierarchical pipe-delimited info panels (example: Chembl | Substructure)
* Fixed #1029: File Manager: clicking the file in the Folder Tree doesn't trigger file preview in the main view
* JS: DataFrame.toByteArray method
* Harmonize membership editor
* Package content validation: check package properties' types
* (Bug) Scripts: Sample table does not open
* Migrate Chem tests from ApiTests to Chem package
* SequenceTranslator: fix validation corner case
* (Bug) __jsThen error during some packages init
* #1014: Monomer cell renderer fix
* Peptides: cleanup
* Peptides: fixed Mutation Cliffs cell colors
* Bio: using grok.shell.tv for splitToMonomers
* Peptides: renamed filterting-statistics to statistics
* Peptides: renamed constant
* Peptides: moved gridCellValidation function
* Peptides: renamed grid cell validation function
* #731: Chem \- sketcher tests
* Docking: highlighting the separator
* Bio: #970 fingerprints
* Add a script to find internal links in help

## 2022-09-08 Dev build 1.6.9

### Addressed Issues

* com
* Fix default value for select (choice) (remove '' and "" from string when chosen default value)
* Implement setRawData for bool columns
* render mols vertically on scatterplot x axis
* lrellis plot not working in LS
* Fixed the method for default choices
* delete useless prints
* bar chart slider positioning
* Grok Spawner: fix default value for GrokSpawnerAddress
* Make modals resizable
* Legends Resizing
* closes #810 Reordered selected columns
* add started/finished/author to funccall
* Layouts | Save to gallery: newly saved layout should appear in the property panel
* Layouts pane: ability to put the layout on the property panel
* CI: Deploy process documentation
* closes #798 selecting dots in jittered scatterplot
* #828 formula lines modal fixed
* #828 order/hide dialog rewworked
* #828 make modal not resizable by default
* JS API: ui.fileBrowser
* Work in progress.
* #831: HitTriage: work in progress
* Moved InputBase descendants to separate files.
* AppEvents.onInputCreated global event
* Ability to reuse database browser, and specify how results are interpreted
* Table Input: ability to get table by executing a database query
* Docs: Publish version from dev to Docker Hub
* Datlas: queue interaction for Grok Connect (WIP)
* Docker image for packages tests
* Scatterplot compact view
* add compact view to scatterplot 1.0
* JS API: View.statusBarPanels (WIP)
* (Bug) Grok connect: decimal, real, and double tests failed
* Grok connect: refactoring of numeric type tests
* JS: UI Forms harmonization (WIP)
* (Bug) Failed to delete package
* DB Indexes
* Fixed filename
* Fixed CSS
* Added a demo dataset with image urls
* PowerGrid: help on linked images
* JS UI Test: Scripting in Add new column (WIP)
* Functions: show function class instead of null in status bar
* Grok connect: refactoring of numeric type tests
* DataFrame.setTag(...) now returns this.
* Renderers: ability to specify column tags as conditions (WIP)
* adds menu click method
* modal click simulation method adjustment
* Datlas connections leaking investigation (WIP)
* Datlas: remove trailing ; in the row checking query
* JS API: Menu.click method
* #831: HitTriage: work in progress \- refactoring
* closes #800 Form does not work after filtering
* Introduced Entity.package Implemented safe delete New permissions check Derived properties support in ORM
* Introduced Entity.package
* GH-686: TreeView events
* modal placement, dont leave page on backspace
* GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Fixed deadlock in exchangeDevKey method
* CSS fixes
* Introduced Entity.package, reviewed Entities list
* Introduced Entity.package, introduced PackageEntityMixin
* Autogenerated code
* Speed-up FilesView start
* Fixed analyzer warnings
* Fixed package.functions section exception
* closes #791 Can not uncheck "Filter out missing values" after saving
* Fixed wrong naming on package deploy
* Removed vis.js from autostart
* Async external connectors fetch
* Speed-up DG start
* FilesView: fix file appearance bug in the card view
* Introduced Entity.package, bugfixes
* Fixed possible NPE
* JS API: BigInt support
* closes #857 Deselection in the scatter plot
* Fixed author parsing
* Package content validation: improve the up-to-date package bundle check
* Fixed tableInput icons
* closes #859 3D scatter plot: implement saving to .png
* Updated beta users
* Tested renderer selection by multiple tags (works fine)
* (Bug) Integration tests: Fix all northwinds for 'External provider: test all Northwinds' test
* refactor: configuration files
* Docs: fix documentation links
* Public submodule
* Add delete FuncCall api
* PackagesView installed filter
* CSS harmonization
* Introduced EntityType.isPackageEntity field
* Fixes #923: Grid: Order or hide columns: provide "move to the top" / "move to the bottom" buttons
* Viewer-specific tooltips (WIP)
* Min() and Max() are not compatible with date and datetime columns #689
* Fixed tests
* removed debug lines
* #closes 858 Cell as an input for widget
* (Bug) FilesView: dock manager exception on opening
* Setup fixes
* Update doc for proxy
* Better transaction management
* Added debug lines
* Some packages not publishing fix
* Removed debug printout
* Add save as png to all viewers
* (Bug) Connection is cloned on sharing
* Reverted GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Dialog: fixed a bug with the resizing
* Grid: Fixed an issue with bool cell renderer
* dapi.files: fixed a bug by reverting to the lazy initialization. Code cleanup.
* fixed: ColumnsInput: an exception when checkedColumns is not specified
* js-api update
* (Bug) Server push is sent twice
* (Bug) Row.setValues doesn't fire onDataChanged event
* Removed debug lines
* Add exceljs to sources
* DateTime scalar functions harmonization
* Fixed shell test
* PowerPack: build fix
* Input fixes, and new disabled styles
* Fixed scalar functions null behavior
* Fixed logical functions in Chrome
* Introduce DateInput with dayjs support
* Fixed deploy exception
* Fixed label css
* (Bug) ScatterPlot: Top and Bottom legend position cuts labels
* (Bug) Selection doesn't work when jitter activated
* ML | Random Data doesn't work ML | Missing Values Imputation doesn't work
* CSS fix
* Fixed dapi and ml tests
* Possible NPE fixed
* Scatterplot 3d detach method
* dapi.root
* View close fixed
* Updated public
* Test: setUpAll and tearDownAll with groups for integration tests
* Test: Skip tests which fail everytime
* Added JS wrapper for EntityRecord
* Packed categories zoom by filter
* package_id migration tuning
* Fixed user login error
* (Bug) FilesUploader: ?layout parameter doesn't work
* CI: Sofiia Podolskaia in beta users
* Test: Use custom session token for integration tests Load session token for integration tests from environment
  variable GROK_SESSION_TOKEN As a fallback value the old session token unit_test_token will be used
* Refactor: Complete refactoring of deploy folder Deprecate old scripts
* Test: File with environment variables for local tests
* GH-954: expanded accessor for TreeViewGroup
* (Bug) Integration tests: External provider (WIP)
* Docs: Instructions to run datlas tests locally on local custom stand
* Test: Add group for all integration tests Fix integration tests which were failing because of session
* closes #936 pie chart improvments
* Test: Add testDevKey config option for docker stands
* Update background for macromolecules
* Test: fix Missing Values Imputation: k-Nearest Neighbour Imputation
* Fixed Datlas tests
* Reverted incorrect fixes
* GH-965: TreeViewGroup.group now also accepts Element
* closes #961 gray out unused categories instead of not showing them
* (Bug) ScatterPlot: Categories rendering issue
* closes #962 Save as PNG: not all elements are added to image
* closes 790 Reset filter doesn't work on "Filter out missing values"
* add button style
* closes #963 Continuous selection on shift+click
* CI: Add short commit to docker image
* Test: Fix Datlas tests
* Load package before loading tests
* Balloons: add the "X" icon
* closes interfaces generation
* (Bug) When user shares an entity from package package is shared too
* Chem harmonization (WIP)
* closes 935 treemap improvments
* Chem: fix fingerprints
* PC Plot: collaborative filtering
* Test: Increase timeout for datlas tests
* (Bug) Packages permission check doesn't work
* closes 935 with all rowSource modes
* closes 945 Line chart: ability to hide and make semi-transparent lines
* FilesView: add the "Add new share" ribbon icon and link below the tree
* FilesView: 'refresh' icon
* CI: Docker Buildx services auto list
* CI: Docker Demo Chembl DB
* CI: Docker Buildx target
* CI: Create test_numeric_types table for external_provider_test
* Test: Increase timeout for Dataset Client: Projects: CRUD: Remove
* Test: Increase timeout for Datlas tests
* Integration tests: datetime test
* Line chart: add tooltips in multiaxis mode
* New FontAwesome version
* Fixed css
* Fix bakcground styles, for macromolecules page
* Chem: restoring sim search, add indeces
* code fix for #935 treemap improvments
* Pipe file downloads
* Folder sharing (WIP)
* Ask credentials on file share open
* Minor style changes
* Sketch View: add the reset button
* JS API: Grid.sortByColumns, Grid.sortTypes
* Chem: fix for similarity search
* Grok connect: performance test (select 1)
* Integration tests: refactor external_provider_perf_test
* Fixed AD login
* GH-955: subDir and share methods
* GH-955: js api update
* (Bug) Files: Unable to open context menu for files in explorer view
* Credentials Manager (WIP)
* (Bug) connection.shares.connection == null condition doesn't work (WIP)
* BarChart: support additive functions 'unique', 'missing value count' and 'value count' in stacked mode (WIP)
* Docker: Datagrok ca-certficates
* Fixed files API
* CI: Add Alexey Chopovsky to beta users
* Docker: Remove elastic from datagrok
* Docs: Upgrade docker instructions with certificates
* PowerGrid: bumped up version to 1.1.0
* #969 forrula-lines vert lines
* Reversed wrong commit
* Public token
* Chem #995: clean up -Similarity
* Chem #995: open phacts clean up
* Enable saving multiple credentials in credential manager (WIP)
* Fixed logout
* Docs: Alllow askalkin@datagrok.ai to view encrypted files
* (Bug) JS: grok.functions.call doesn't handle error
* JS API: correct return type for viewer.getOptions
* #953: Filters turned off status is not saved in layout
* #966: Deleted column filter should not be filtered on
* #969 new classes system ,aybe should be reworked
* Change order of interop type checks
* #969 linechart hittest
* closes #996 Tooltip not working for partly truncated values
* Funcs: JS-based context functions for multiple columns (WIP)
* (Bug) ColumnsInput: number of selected columns is not shown upon construction
* JS API: GridColumn.move method
* closes #1015: PowerGrid: Sparklines: Ability to add sparklines from the "Actions" pane for selected columns
* JS API: GridColumn.scrollIntoView, RangeSlider.scrollBy, RangeSlider.scrollTo
* PowerGrid: Sparklines: accessing settings from the hamburger menu
* #930: PowerGrid: sparklines: ability to rename and remove summary columns
* #930: PowerGrid: sparklines: managing column order
* #930: PowerGrid: sparklines: setting the sparkline column as a current object after it is created
* Possible NPE fix
* Deprecated old tests
* Fixed package init
* (Bug) Package Credentials Editor doesn't work (WIP)
* Better wording
* (Bug) grok.shell.registerViewer fails
* Test: Skip Grok Compute tests because of compatility

## 2022-09-07 Dev build 1.6.8

### Addressed Issues

* com
* Fix default value for select (choice) (remove '' and "" from string when chosen default value)
* Implement setRawData for bool columns
* render mols vertically on scatterplot x axis
* lrellis plot not working in LS
* Fixed the method for default choices
* delete useless prints
* bar chart slider positioning
* Grok Spawner: fix default value for GrokSpawnerAddress
* Make modals resizable
* Legends Resizing
* closes #810 Reordered selected columns
* add started/finished/author to funccall
* Layouts | Save to gallery: newly saved layout should appear in the property panel
* Layouts pane: ability to put the layout on the property panel
* CI: Deploy process documentation
* closes #798 selecting dots in jittered scatterplot
* #828 formula lines modal fixed
* #828 order/hide dialog rewworked
* #828 make modal not resizable by default
* JS API: ui.fileBrowser
* Work in progress.
* #831: HitTriage: work in progress
* Moved InputBase descendants to separate files.
* AppEvents.onInputCreated global event
* Ability to reuse database browser, and specify how results are interpreted
* Table Input: ability to get table by executing a database query
* Docs: Publish version from dev to Docker Hub
* Datlas: queue interaction for Grok Connect (WIP)
* Docker image for packages tests
* add compact view to scatterplot 1.0
* Scatterplot compact view
* JS API: View.statusBarPanels (WIP)
* (Bug) Grok connect: decimal, real, and double tests failed
* Grok connect: refactoring of numeric type tests
* JS: UI Forms harmonization (WIP)
* (Bug) Failed to delete package
* DB Indexes
* Fixed filename
* Fixed CSS
* Added a demo dataset with image urls
* PowerGrid: help on linked images
* JS UI Test: Scripting in Add new column (WIP)
* Functions: show function class instead of null in status bar
* Grok connect: refactoring of numeric type tests
* DataFrame.setTag(...) now returns this.
* Renderers: ability to specify column tags as conditions (WIP)
* adds menu click method
* modal click simulation method adjustment
* Datlas connections leaking investigation (WIP)
* Datlas: remove trailing ; in the row checking query
* JS API: Menu.click method
* #831: HitTriage: work in progress \- refactoring
* closes #800 Form does not work after filtering
* Introduced Entity.package Implemented safe delete New permissions check Derived properties support in ORM
* Introduced Entity.package
* GH-686: TreeView events
* modal placement, dont leave page on backspace
* GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Fixed deadlock in exchangeDevKey method
* CSS fixes
* Introduced Entity.package, reviewed Entities list
* Introduced Entity.package, introduced PackageEntityMixin
* Autogenerated code
* Speed-up FilesView start
* Fixed analyzer warnings
* Fixed package.functions section exception
* closes #791 Can not uncheck "Filter out missing values" after saving
* Fixed wrong naming on package deploy
* Removed vis.js from autostart
* Async external connectors fetch
* Speed-up DG start
* FilesView: fix file appearance bug in the card view
* Introduced Entity.package, bugfixes
* Fixed possible NPE
* JS API: BigInt support
* closes #857 Deselection in the scatter plot
* Fixed author parsing
* Package content validation: improve the up-to-date package bundle check
* Fixed tableInput icons
* closes #859 3D scatter plot: implement saving to .png
* Updated beta users
* Tested renderer selection by multiple tags (works fine)
* (Bug) Integration tests: Fix all northwinds for 'External provider: test all Northwinds' test
* refactor: configuration files
* Docs: fix documentation links
* Public submodule
* Add delete FuncCall api
* PackagesView installed filter
* CSS harmonization
* Introduced EntityType.isPackageEntity field
* Fixes #923: Grid: Order or hide columns: provide "move to the top" / "move to the bottom" buttons
* Viewer-specific tooltips (WIP)
* Min() and Max() are not compatible with date and datetime columns #689
* Fixed tests
* removed debug lines
* #closes 858 Cell as an input for widget
* (Bug) FilesView: dock manager exception on opening
* Setup fixes
* Update doc for proxy
* Better transaction management
* Added debug lines
* Some packages not publishing fix
* Removed debug printout
* Add save as png to all viewers
* (Bug) Connection is cloned on sharing
* Reverted GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Dialog: fixed a bug with the resizing
* Grid: Fixed an issue with bool cell renderer
* dapi.files: fixed a bug by reverting to the lazy initialization. Code cleanup.
* fixed: ColumnsInput: an exception when checkedColumns is not specified
* js-api update
* (Bug) Server push is sent twice
* (Bug) Row.setValues doesn't fire onDataChanged event
* Removed debug lines
* Add exceljs to sources
* DateTime scalar functions harmonization
* Fixed shell test
* PowerPack: build fix
* Input fixes, and new disabled styles
* Fixed scalar functions null behavior
* Fixed logical functions in Chrome
* Introduce DateInput with dayjs support
* Fixed deploy exception
* Fixed label css
* (Bug) ScatterPlot: Top and Bottom legend position cuts labels
* (Bug) Selection doesn't work when jitter activated
* ML | Random Data doesn't work ML | Missing Values Imputation doesn't work
* CSS fix
* Fixed dapi and ml tests
* Possible NPE fixed
* Scatterplot 3d detach method
* dapi.root
* View close fixed
* Updated public
* Test: setUpAll and tearDownAll with groups for integration tests
* Test: Skip tests which fail everytime
* Added JS wrapper for EntityRecord
* Packed categories zoom by filter
* package_id migration tuning
* Fixed user login error
* (Bug) FilesUploader: ?layout parameter doesn't work
* CI: Sofiia Podolskaia in beta users
* Test: Use custom session token for integration tests Load session token for integration tests from environment
  variable GROK_SESSION_TOKEN As a fallback value the old session token unit_test_token will be used
* Refactor: Complete refactoring of deploy folder Deprecate old scripts
* Test: File with environment variables for local tests
* GH-954: expanded accessor for TreeViewGroup
* (Bug) Integration tests: External provider (WIP)
* Docs: Instructions to run datlas tests locally on local custom stand
* Test: Add group for all integration tests Fix integration tests which were failing because of session
* closes #936 pie chart improvments
* Test: Add testDevKey config option for docker stands
* Update background for macromolecules
* Test: fix Missing Values Imputation: k-Nearest Neighbour Imputation
* Fixed Datlas tests
* Reverted incorrect fixes
* GH-965: TreeViewGroup.group now also accepts Element
* closes #961 gray out unused categories instead of not showing them
* (Bug) ScatterPlot: Categories rendering issue
* closes #962 Save as PNG: not all elements are added to image
* closes 790 Reset filter doesn't work on "Filter out missing values"
* add button style
* closes #963 Continuous selection on shift+click
* CI: Add short commit to docker image
* Test: Fix Datlas tests
* Load package before loading tests
* Balloons: add the "X" icon
* closes interfaces generation
* (Bug) When user shares an entity from package package is shared too
* Chem harmonization (WIP)
* closes 935 treemap improvments
* Chem: fix fingerprints
* PC Plot: collaborative filtering
* Test: Increase timeout for datlas tests
* (Bug) Packages permission check doesn't work
* closes 935 with all rowSource modes
* closes 945 Line chart: ability to hide and make semi-transparent lines
* FilesView: add the "Add new share" ribbon icon and link below the tree
* FilesView: 'refresh' icon
* CI: Docker Buildx services auto list
* CI: Docker Demo Chembl DB
* CI: Docker Buildx target
* CI: Create test_numeric_types table for external_provider_test
* Test: Increase timeout for Dataset Client: Projects: CRUD: Remove
* Test: Increase timeout for Datlas tests
* Integration tests: datetime test
* Line chart: add tooltips in multiaxis mode
* New FontAwesome version
* Fixed css
* Fix bakcground styles, for macromolecules page
* Chem: restoring sim search, add indeces
* code fix for #935 treemap improvments
* Pipe file downloads
* Folder sharing (WIP)
* Ask credentials on file share open
* Minor style changes
* Sketch View: add the reset button
* JS API: Grid.sortByColumns, Grid.sortTypes
* Chem: fix for similarity search
* Grok connect: performance test (select 1)
* Integration tests: refactor external_provider_perf_test
* Fixed AD login
* GH-955: subDir and share methods
* GH-955: js api update
* (Bug) Files: Unable to open context menu for files in explorer view
* Credentials Manager (WIP)
* (Bug) connection.shares.connection == null condition doesn't work (WIP)
* BarChart: support additive functions 'unique', 'missing value count' and 'value count' in stacked mode (WIP)
* Docker: Datagrok ca-certficates
* Fixed files API
* CI: Add Alexey Chopovsky to beta users
* Docker: Remove elastic from datagrok
* Docs: Upgrade docker instructions with certificates
* PowerGrid: bumped up version to 1.1.0
* #969 forrula-lines vert lines
* Reversed wrong commit
* Public token
* Chem #995: clean up -Similarity
* Chem #995: open phacts clean up
* Enable saving multiple credentials in credential manager (WIP)
* Fixed logout
* Docs: Alllow askalkin@datagrok.ai to view encrypted files
* (Bug) JS: grok.functions.call doesn't handle error
* JS API: correct return type for viewer.getOptions
* #953: Filters turned off status is not saved in layout
* #966: Deleted column filter should not be filtered on
* #969 new classes system ,aybe should be reworked
* Change order of interop type checks
* #969 linechart hittest
* closes #996 Tooltip not working for partly truncated values
* Funcs: JS-based context functions for multiple columns (WIP)
* (Bug) ColumnsInput: number of selected columns is not shown upon construction
* JS API: GridColumn.move method
* closes #1015: PowerGrid: Sparklines: Ability to add sparklines from the "Actions" pane for selected columns
* JS API: GridColumn.scrollIntoView, RangeSlider.scrollBy, RangeSlider.scrollTo
* PowerGrid: Sparklines: accessing settings from the hamburger menu
* #930: PowerGrid: sparklines: ability to rename and remove summary columns
* #930: PowerGrid: sparklines: managing column order
* #930: PowerGrid: sparklines: setting the sparkline column as a current object after it is created
* Possible NPE fix
* Deprecated old tests
* Fixed package init
* added share-the-folder.gif
* JS API: add stats enums
* JS API: document Column.aggregate
* bio lib: move ntseq from removed Sequence
* Sequence: Remove
* Wiki: replace links
* Moved Markdown.md to Datagrok folder
* Chem: code cleanup and minor refactoring
* JS API: code cleanup and minor refactoring
* PowerPack: testing \- WIP
* PowerGrid: a sample that demonstrates adding sparklines programmatically
* GIS: midterm commit \- detectors, working with census
* PowerGrid: small fix, moved constants to separate variables
* (Bug) Package Credentials Editor doesn't work (WIP)
* Better wording

## 2022-09-06 Dev build 1.6.7

### Addressed Issues

* com
* Fix default value for select (choice) (remove '' and "" from string when chosen default value)
* Implement setRawData for bool columns
* render mols vertically on scatterplot x axis
* lrellis plot not working in LS
* Fixed the method for default choices
* delete useless prints
* bar chart slider positioning
* Grok Spawner: fix default value for GrokSpawnerAddress
* Make modals resizable
* Legends Resizing
* closes #810 Reordered selected columns
* add started/finished/author to funccall
* Layouts | Save to gallery: newly saved layout should appear in the property panel
* Layouts pane: ability to put the layout on the property panel
* CI: Deploy process documentation
* closes #798 selecting dots in jittered scatterplot
* #828 formula lines modal fixed
* #828 order/hide dialog rewworked
* #828 make modal not resizable by default
* JS API: ui.fileBrowser
* Work in progress.
* #831: HitTriage: work in progress
* Moved InputBase descendants to separate files.
* AppEvents.onInputCreated global event
* Ability to reuse database browser, and specify how results are interpreted
* Table Input: ability to get table by executing a database query
* Docs: Publish version from dev to Docker Hub
* Datlas: queue interaction for Grok Connect (WIP)
* Docker image for packages tests
* add compact view to scatterplot 1.0
* Scatterplot compact view
* JS API: View.statusBarPanels (WIP)
* (Bug) Grok connect: decimal, real, and double tests failed
* Grok connect: refactoring of numeric type tests
* JS: UI Forms harmonization (WIP)
* (Bug) Failed to delete package
* DB Indexes
* Fixed filename
* Fixed CSS
* Added a demo dataset with image urls
* PowerGrid: help on linked images
* JS UI Test: Scripting in Add new column (WIP)
* Functions: show function class instead of null in status bar
* Grok connect: refactoring of numeric type tests
* DataFrame.setTag(...) now returns this.
* Renderers: ability to specify column tags as conditions (WIP)
* adds menu click method
* modal click simulation method adjustment
* Datlas connections leaking investigation (WIP)
* Datlas: remove trailing ; in the row checking query
* JS API: Menu.click method
* #831: HitTriage: work in progress \- refactoring
* closes #800 Form does not work after filtering
* Introduced Entity.package Implemented safe delete New permissions check Derived properties support in ORM
* Introduced Entity.package
* GH-686: TreeView events
* modal placement, dont leave page on backspace
* GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Fixed deadlock in exchangeDevKey method
* CSS fixes
* Introduced Entity.package, reviewed Entities list
* Introduced Entity.package, introduced PackageEntityMixin
* Autogenerated code
* Speed-up FilesView start
* Fixed analyzer warnings
* Fixed package.functions section exception
* closes #791 Can not uncheck "Filter out missing values" after saving
* Fixed wrong naming on package deploy
* Removed vis.js from autostart
* Async external connectors fetch
* Speed-up DG start
* FilesView: fix file appearance bug in the card view
* Introduced Entity.package, bugfixes
* Fixed possible NPE
* JS API: BigInt support
* closes #857 Deselection in the scatter plot
* Fixed author parsing
* Package content validation: improve the up-to-date package bundle check
* Fixed tableInput icons
* closes #859 3D scatter plot: implement saving to .png
* Updated beta users
* Tested renderer selection by multiple tags (works fine)
* (Bug) Integration tests: Fix all northwinds for 'External provider: test all Northwinds' test
* refactor: configuration files
* Docs: fix documentation links
* Public submodule
* Add delete FuncCall api
* PackagesView installed filter
* CSS harmonization
* Introduced EntityType.isPackageEntity field
* Fixes #923: Grid: Order or hide columns: provide "move to the top" / "move to the bottom" buttons
* Viewer-specific tooltips (WIP)
* Min() and Max() are not compatible with date and datetime columns #689
* Fixed tests
* removed debug lines
* #closes 858 Cell as an input for widget
* (Bug) FilesView: dock manager exception on opening
* Setup fixes
* Update doc for proxy
* Better transaction management
* Added debug lines
* Some packages not publishing fix
* Removed debug printout
* Add save as png to all viewers
* (Bug) Connection is cloned on sharing
* Reverted GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Dialog: fixed a bug with the resizing
* Grid: Fixed an issue with bool cell renderer
* dapi.files: fixed a bug by reverting to the lazy initialization. Code cleanup.
* fixed: ColumnsInput: an exception when checkedColumns is not specified
* js-api update
* (Bug) Server push is sent twice
* (Bug) Row.setValues doesn't fire onDataChanged event
* Removed debug lines
* Add exceljs to sources
* DateTime scalar functions harmonization
* Fixed shell test
* PowerPack: build fix
* Input fixes, and new disabled styles
* Fixed scalar functions null behavior
* Fixed logical functions in Chrome
* Introduce DateInput with dayjs support
* Fixed deploy exception
* Fixed label css
* (Bug) ScatterPlot: Top and Bottom legend position cuts labels
* (Bug) Selection doesn't work when jitter activated
* ML | Random Data doesn't work ML | Missing Values Imputation doesn't work
* CSS fix
* Fixed dapi and ml tests
* Possible NPE fixed
* Scatterplot 3d detach method
* dapi.root
* View close fixed
* Updated public
* Test: setUpAll and tearDownAll with groups for integration tests
* Test: Skip tests which fail everytime
* Added JS wrapper for EntityRecord
* Packed categories zoom by filter
* package_id migration tuning
* Fixed user login error
* (Bug) FilesUploader: ?layout parameter doesn't work
* CI: Sofiia Podolskaia in beta users
* Test: Use custom session token for integration tests Load session token for integration tests from environment
  variable GROK_SESSION_TOKEN As a fallback value the old session token unit_test_token will be used
* Refactor: Complete refactoring of deploy folder Deprecate old scripts
* Test: File with environment variables for local tests
* GH-954: expanded accessor for TreeViewGroup
* (Bug) Integration tests: External provider (WIP)
* Docs: Instructions to run datlas tests locally on local custom stand
* Test: Add group for all integration tests Fix integration tests which were failing because of session
* closes #936 pie chart improvments
* Test: Add testDevKey config option for docker stands
* Update background for macromolecules
* Test: fix Missing Values Imputation: k-Nearest Neighbour Imputation
* Fixed Datlas tests
* Reverted incorrect fixes
* GH-965: TreeViewGroup.group now also accepts Element
* closes #961 gray out unused categories instead of not showing them
* (Bug) ScatterPlot: Categories rendering issue
* closes #962 Save as PNG: not all elements are added to image
* closes 790 Reset filter doesn't work on "Filter out missing values"
* add button style
* closes #963 Continuous selection on shift+click
* CI: Add short commit to docker image
* Test: Fix Datlas tests
* Load package before loading tests
* Balloons: add the "X" icon
* closes interfaces generation
* (Bug) When user shares an entity from package package is shared too
* Chem harmonization (WIP)
* closes 935 treemap improvments
* Chem: fix fingerprints
* PC Plot: collaborative filtering
* Test: Increase timeout for datlas tests
* (Bug) Packages permission check doesn't work
* closes 935 with all rowSource modes
* closes 945 Line chart: ability to hide and make semi-transparent lines
* FilesView: add the "Add new share" ribbon icon and link below the tree
* FilesView: 'refresh' icon
* CI: Docker Buildx services auto list
* CI: Docker Demo Chembl DB
* CI: Docker Buildx target
* CI: Create test_numeric_types table for external_provider_test
* Test: Increase timeout for Dataset Client: Projects: CRUD: Remove
* Test: Increase timeout for Datlas tests
* Integration tests: datetime test
* Line chart: add tooltips in multiaxis mode
* New FontAwesome version
* Fixed css
* Fix bakcground styles, for macromolecules page
* Chem: restoring sim search, add indeces
* code fix for #935 treemap improvments
* Pipe file downloads
* Folder sharing (WIP)
* Ask credentials on file share open
* Minor style changes
* Sketch View: add the reset button
* JS API: Grid.sortByColumns, Grid.sortTypes
* Chem: fix for similarity search
* Grok connect: performance test (select 1)
* Integration tests: refactor external_provider_perf_test
* Fixed AD login
* GH-955: subDir and share methods
* GH-955: js api update
* (Bug) Files: Unable to open context menu for files in explorer view
* Credentials Manager (WIP)
* (Bug) connection.shares.connection == null condition doesn't work (WIP)
* BarChart: support additive functions 'unique', 'missing value count' and 'value count' in stacked mode (WIP)
* Docker: Datagrok ca-certficates
* Fixed files API
* CI: Add Alexey Chopovsky to beta users
* Docker: Remove elastic from datagrok
* Docs: Upgrade docker instructions with certificates
* PowerGrid: bumped up version to 1.1.0
* #969 forrula-lines vert lines
* Reversed wrong commit
* Public token
* Chem #995: clean up -Similarity
* Chem #995: open phacts clean up
* Enable saving multiple credentials in credential manager (WIP)
* Fixed logout
* Docs: Alllow askalkin@datagrok.ai to view encrypted files
* (Bug) JS: grok.functions.call doesn't handle error
* JS API: correct return type for viewer.getOptions
* #953: Filters turned off status is not saved in layout
* #966: Deleted column filter should not be filtered on
* #969 new classes system ,aybe should be reworked
* Change order of interop type checks
* #969 linechart hittest
* closes #996 Tooltip not working for partly truncated values
* #864: splitToMonomers function
* Issue #976 Chem: elemental analysis code refactoring
* Issue #976: Chem: elemental analysis (add conversion from smiles to molblock)
* Funcs: JS-based context functions for multiple columns (WIP)
* JS API: Point.distanceTo(p)
* JS API: More documentation
* JS API: InputBase: addCaption, addPostfix, addOptions, setTooltip now return InputBase (useful for call chaining)
* (Bug) ColumnsInput: number of selected columns is not shown upon construction
* JS API: GridColumn.move method
* JS API: documentation
* closes #1015: PowerGrid: Sparklines: Ability to add sparklines from the "Actions" pane for selected columns
* JS API: GridColumn.scrollIntoView, RangeSlider.scrollBy, RangeSlider.scrollTo
* PowerGrid: Sparklines: accessing settings from the hamburger menu
* #930: PowerGrid: sparklines: ability to rename and remove summary columns
* #930: PowerGrid: sparklines: not including datatime columns by default
* #930: PowerGrid: sparklines: managing column order
* #930: PowerGrid: sparklines: setting the sparkline column as a current object after it is created
* Test utils
* Possible NPE fix
* Deprecated old tests
* Fixed package init

## 2022-09-01 Dev build 1.6.6

### Addressed Issues

* com
* Fix default value for select (choice) (remove '' and "" from string when chosen default value)
* Implement setRawData for bool columns
* render mols vertically on scatterplot x axis
* lrellis plot not working in LS
* Fixed the method for default choices
* delete useless prints
* bar chart slider positioning
* Grok Spawner: fix default value for GrokSpawnerAddress
* Make modals resizable
* Legends Resizing
* closes #810 Reordered selected columns
* add started/finished/author to funccall
* Layouts | Save to gallery: newly saved layout should appear in the property panel
* Layouts pane: ability to put the layout on the property panel
* CI: Deploy process documentation
* closes #798 selecting dots in jittered scatterplot
* #828 formula lines modal fixed
* #828 order/hide dialog rewworked
* #828 make modal not resizable by default
* JS API: ui.fileBrowser
* Work in progress.
* #831: HitTriage: work in progress
* Moved InputBase descendants to separate files.
* AppEvents.onInputCreated global event
* Ability to reuse database browser, and specify how results are interpreted
* Table Input: ability to get table by executing a database query
* Docs: Publish version from dev to Docker Hub
* Datlas: queue interaction for Grok Connect (WIP)
* Docker image for packages tests
* add compact view to scatterplot 1.0
* Scatterplot compact view
* JS API: View.statusBarPanels (WIP)
* (Bug) Grok connect: decimal, real, and double tests failed
* Grok connect: refactoring of numeric type tests
* JS: UI Forms harmonization (WIP)
* (Bug) Failed to delete package
* DB Indexes
* Fixed filename
* Fixed CSS
* Added a demo dataset with image urls
* PowerGrid: help on linked images
* JS UI Test: Scripting in Add new column (WIP)
* Functions: show function class instead of null in status bar
* Grok connect: refactoring of numeric type tests
* DataFrame.setTag(...) now returns this.
* Renderers: ability to specify column tags as conditions (WIP)
* adds menu click method
* modal click simulation method adjustment
* Datlas connections leaking investigation (WIP)
* Datlas: remove trailing ; in the row checking query
* JS API: Menu.click method
* #831: HitTriage: work in progress \- refactoring
* closes #800 Form does not work after filtering
* Introduced Entity.package Implemented safe delete New permissions check Derived properties support in ORM
* Introduced Entity.package
* GH-686: TreeView events
* modal placement, dont leave page on backspace
* GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Fixed deadlock in exchangeDevKey method
* CSS fixes
* Introduced Entity.package, reviewed Entities list
* Introduced Entity.package, introduced PackageEntityMixin
* Autogenerated code
* Speed-up FilesView start
* Fixed analyzer warnings
* Fixed package.functions section exception
* closes #791 Can not uncheck "Filter out missing values" after saving
* Fixed wrong naming on package deploy
* Removed vis.js from autostart
* Async external connectors fetch
* Speed-up DG start
* FilesView: fix file appearance bug in the card view
* Introduced Entity.package, bugfixes
* Fixed possible NPE
* JS API: BigInt support
* closes #857 Deselection in the scatter plot
* Fixed author parsing
* Package content validation: improve the up-to-date package bundle check
* Fixed tableInput icons
* closes #859 3D scatter plot: implement saving to .png
* Updated beta users
* Tested renderer selection by multiple tags (works fine)
* (Bug) Integration tests: Fix all northwinds for 'External provider: test all Northwinds' test
* refactor: configuration files
* Docs: fix documentation links
* Public submodule
* Add delete FuncCall api
* PackagesView installed filter
* CSS harmonization
* Introduced EntityType.isPackageEntity field
* Fixes #923: Grid: Order or hide columns: provide "move to the top" / "move to the bottom" buttons
* Viewer-specific tooltips (WIP)
* Min() and Max() are not compatible with date and datetime columns #689
* Fixed tests
* removed debug lines
* #closes 858 Cell as an input for widget
* (Bug) FilesView: dock manager exception on opening
* Setup fixes
* Update doc for proxy
* Better transaction management
* Added debug lines
* Some packages not publishing fix
* Removed debug printout
* Add save as png to all viewers
* (Bug) Connection is cloned on sharing
* Reverted GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Dialog: fixed a bug with the resizing
* Grid: Fixed an issue with bool cell renderer
* dapi.files: fixed a bug by reverting to the lazy initialization. Code cleanup.
* fixed: ColumnsInput: an exception when checkedColumns is not specified
* js-api update
* (Bug) Server push is sent twice
* (Bug) Row.setValues doesn't fire onDataChanged event
* Removed debug lines
* Add exceljs to sources
* DateTime scalar functions harmonization
* Fixed shell test
* PowerPack: build fix
* Input fixes, and new disabled styles
* Fixed scalar functions null behavior
* Fixed logical functions in Chrome
* Introduce DateInput with dayjs support
* Fixed deploy exception
* Fixed label css
* (Bug) ScatterPlot: Top and Bottom legend position cuts labels
* (Bug) Selection doesn't work when jitter activated
* ML | Random Data doesn't work ML | Missing Values Imputation doesn't work
* CSS fix
* Fixed dapi and ml tests
* Possible NPE fixed
* Scatterplot 3d detach method
* dapi.root
* View close fixed
* Updated public
* Test: setUpAll and tearDownAll with groups for integration tests
* Test: Skip tests which fail everytime
* Added JS wrapper for EntityRecord
* Packed categories zoom by filter
* package_id migration tuning
* Fixed user login error
* (Bug) FilesUploader: ?layout parameter doesn't work
* CI: Sofiia Podolskaia in beta users
* Test: Use custom session token for integration tests Load session token for integration tests from environment
  variable GROK_SESSION_TOKEN As a fallback value the old session token unit_test_token will be used
* Refactor: Complete refactoring of deploy folder Deprecate old scripts
* Test: File with environment variables for local tests
* GH-954: expanded accessor for TreeViewGroup
* (Bug) Integration tests: External provider (WIP)
* Docs: Instructions to run datlas tests locally on local custom stand
* Test: Add group for all integration tests Fix integration tests which were failing because of session
* closes #936 pie chart improvments
* Test: Add testDevKey config option for docker stands
* Update background for macromolecules
* Test: fix Missing Values Imputation: k-Nearest Neighbour Imputation
* Fixed Datlas tests
* Reverted incorrect fixes
* GH-965: TreeViewGroup.group now also accepts Element
* closes #961 gray out unused categories instead of not showing them
* (Bug) ScatterPlot: Categories rendering issue
* closes #962 Save as PNG: not all elements are added to image
* closes 790 Reset filter doesn't work on "Filter out missing values"
* add button style
* closes #963 Continuous selection on shift+click
* CI: Add short commit to docker image
* Test: Fix Datlas tests
* Load package before loading tests
* Balloons: add the "X" icon
* closes interfaces generation
* (Bug) When user shares an entity from package package is shared too
* Chem harmonization (WIP)
* closes 935 treemap improvments
* Chem: fix fingerprints
* PC Plot: collaborative filtering
* Test: Increase timeout for datlas tests
* (Bug) Packages permission check doesn't work
* closes 935 with all rowSource modes
* closes 945 Line chart: ability to hide and make semi-transparent lines
* FilesView: add the "Add new share" ribbon icon and link below the tree
* FilesView: 'refresh' icon
* CI: Docker Buildx services auto list
* CI: Docker Demo Chembl DB
* CI: Docker Buildx target
* CI: Create test_numeric_types table for external_provider_test
* Test: Increase timeout for Dataset Client: Projects: CRUD: Remove
* Test: Increase timeout for Datlas tests
* Integration tests: datetime test
* Line chart: add tooltips in multiaxis mode
* New FontAwesome version
* Fixed css
* Fix bakcground styles, for macromolecules page
* Chem: restoring sim search, add indeces
* code fix for #935 treemap improvments
* Pipe file downloads
* Folder sharing (WIP)
* Ask credentials on file share open
* Minor style changes
* Sketch View: add the reset button
* JS API: Grid.sortByColumns, Grid.sortTypes
* Chem: fix for similarity search
* Grok connect: performance test (select 1)
* Integration tests: refactor external_provider_perf_test
* Fixed AD login
* GH-955: subDir and share methods
* GH-955: js api update
* (Bug) Files: Unable to open context menu for files in explorer view
* Credentials Manager (WIP)
* (Bug) connection.shares.connection == null condition doesn't work (WIP)
* BarChart: support additive functions 'unique', 'missing value count' and 'value count' in stacked mode (WIP)
* Docker: Datagrok ca-certficates
* Fixed files API
* CI: Add Alexey Chopovsky to beta users
* Docker: Remove elastic from datagrok
* Docs: Upgrade docker instructions with certificates
* PowerGrid: bumped up version to 1.1.0
* Reversed wrong commit
* Public token
* Chem #995: clean up -Similarity
* Chem #995: open phacts clean up
* Enable saving multiple credentials in credential manager (WIP)
* Fixed logout
* Docs: Alllow askalkin@datagrok.ai to view encrypted files
* (Bug) JS: grok.functions.call doesn't handle error
* Bio: Fix typing for grok.functions and RDKit interfaces, fix WebLogo-positions tests

## 2022-08-30 Dev build 1.6.5

### Addressed Issues

* com
* Fix default value for select (choice) (remove '' and "" from string when chosen default value)
* Implement setRawData for bool columns
* render mols vertically on scatterplot x axis
* lrellis plot not working in LS
* Fixed the method for default choices
* delete useless prints
* bar chart slider positioning
* Grok Spawner: fix default value for GrokSpawnerAddress
* Make modals resizable
* Legends Resizing
* closes #810 Reordered selected columns
* add started/finished/author to funccall
* Layouts | Save to gallery: newly saved layout should appear in the property panel
* Layouts pane: ability to put the layout on the property panel
* CI: Deploy process documentation
* closes #798 selecting dots in jittered scatterplot
* #828 formula lines modal fixed
* #828 order/hide dialog rewworked
* #828 make modal not resizable by default
* JS API: ui.fileBrowser
* Work in progress.
* #831: HitTriage: work in progress
* Moved InputBase descendants to separate files.
* AppEvents.onInputCreated global event
* Ability to reuse database browser, and specify how results are interpreted
* Table Input: ability to get table by executing a database query
* Docs: Publish version from dev to Docker Hub
* Datlas: queue interaction for Grok Connect (WIP)
* Docker image for packages tests
* add compact view to scatterplot 1.0
* Scatterplot compact view
* JS API: View.statusBarPanels (WIP)
* (Bug) Grok connect: decimal, real, and double tests failed
* Grok connect: refactoring of numeric type tests
* JS: UI Forms harmonization (WIP)
* (Bug) Failed to delete package
* DB Indexes
* Fixed filename
* Fixed CSS
* Added a demo dataset with image urls
* PowerGrid: help on linked images
* JS UI Test: Scripting in Add new column (WIP)
* Functions: show function class instead of null in status bar
* Grok connect: refactoring of numeric type tests
* DataFrame.setTag(...) now returns this.
* Renderers: ability to specify column tags as conditions (WIP)
* adds menu click method
* modal click simulation method adjustment
* Datlas connections leaking investigation (WIP)
* Datlas: remove trailing ; in the row checking query
* JS API: Menu.click method
* #831: HitTriage: work in progress \- refactoring
* closes #800 Form does not work after filtering
* Introduced Entity.package Implemented safe delete New permissions check Derived properties support in ORM
* Introduced Entity.package
* GH-686: TreeView events
* modal placement, dont leave page on backspace
* GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Fixed deadlock in exchangeDevKey method
* CSS fixes
* Introduced Entity.package, reviewed Entities list
* Introduced Entity.package, introduced PackageEntityMixin
* Autogenerated code
* Speed-up FilesView start
* Fixed analyzer warnings
* Fixed package.functions section exception
* closes #791 Can not uncheck "Filter out missing values" after saving
* Fixed wrong naming on package deploy
* Removed vis.js from autostart
* Async external connectors fetch
* Speed-up DG start
* FilesView: fix file appearance bug in the card view
* Introduced Entity.package, bugfixes
* Fixed possible NPE
* JS API: BigInt support
* closes #857 Deselection in the scatter plot
* Fixed author parsing
* Package content validation: improve the up-to-date package bundle check
* Fixed tableInput icons
* closes #859 3D scatter plot: implement saving to .png
* Updated beta users
* Tested renderer selection by multiple tags (works fine)
* (Bug) Integration tests: Fix all northwinds for 'External provider: test all Northwinds' test
* refactor: configuration files
* Docs: fix documentation links
* Public submodule
* Add delete FuncCall api
* PackagesView installed filter
* CSS harmonization
* Introduced EntityType.isPackageEntity field
* Fixes #923: Grid: Order or hide columns: provide "move to the top" / "move to the bottom" buttons
* Viewer-specific tooltips (WIP)
* Min() and Max() are not compatible with date and datetime columns #689
* Fixed tests
* removed debug lines
* #closes 858 Cell as an input for widget
* (Bug) FilesView: dock manager exception on opening
* Setup fixes
* Update doc for proxy
* Better transaction management
* Added debug lines
* Some packages not publishing fix
* Removed debug printout
* Add save as png to all viewers
* (Bug) Connection is cloned on sharing
* Reverted GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Dialog: fixed a bug with the resizing
* Grid: Fixed an issue with bool cell renderer
* dapi.files: fixed a bug by reverting to the lazy initialization. Code cleanup.
* fixed: ColumnsInput: an exception when checkedColumns is not specified
* js-api update
* (Bug) Server push is sent twice
* (Bug) Row.setValues doesn't fire onDataChanged event
* Removed debug lines
* Add exceljs to sources
* DateTime scalar functions harmonization
* Fixed shell test
* PowerPack: build fix
* Input fixes, and new disabled styles
* Fixed scalar functions null behavior
* Fixed logical functions in Chrome
* Introduce DateInput with dayjs support
* Fixed deploy exception
* Fixed label css
* (Bug) ScatterPlot: Top and Bottom legend position cuts labels
* (Bug) Selection doesn't work when jitter activated
* ML | Random Data doesn't work ML | Missing Values Imputation doesn't work
* CSS fix
* Fixed dapi and ml tests
* Possible NPE fixed
* Scatterplot 3d detach method
* dapi.root
* View close fixed
* Updated public
* Test: setUpAll and tearDownAll with groups for integration tests
* Test: Skip tests which fail everytime
* Added JS wrapper for EntityRecord
* Packed categories zoom by filter
* package_id migration tuning
* Fixed user login error
* (Bug) FilesUploader: ?layout parameter doesn't work
* CI: Sofiia Podolskaia in beta users
* Test: Use custom session token for integration tests Load session token for integration tests from environment
  variable GROK_SESSION_TOKEN As a fallback value the old session token unit_test_token will be used
* Refactor: Complete refactoring of deploy folder Deprecate old scripts
* Test: File with environment variables for local tests
* GH-954: expanded accessor for TreeViewGroup
* (Bug) Integration tests: External provider (WIP)
* Docs: Instructions to run datlas tests locally on local custom stand
* Test: Add group for all integration tests Fix integration tests which were failing because of session
* closes #936 pie chart improvments
* Test: Add testDevKey config option for docker stands
* Update background for macromolecules
* Test: fix Missing Values Imputation: k-Nearest Neighbour Imputation
* Fixed Datlas tests
* Reverted incorrect fixes
* GH-965: TreeViewGroup.group now also accepts Element
* closes #961 gray out unused categories instead of not showing them
* (Bug) ScatterPlot: Categories rendering issue
* closes #962 Save as PNG: not all elements are added to image
* closes 790 Reset filter doesn't work on "Filter out missing values"
* add button style
* closes #963 Continuous selection on shift+click
* CI: Add short commit to docker image
* Test: Fix Datlas tests
* Load package before loading tests
* Balloons: add the "X" icon
* closes interfaces generation
* (Bug) When user shares an entity from package package is shared too
* Chem harmonization (WIP)
* closes 935 treemap improvments
* Chem: fix fingerprints
* PC Plot: collaborative filtering
* Test: Increase timeout for datlas tests
* (Bug) Packages permission check doesn't work
* closes 935 with all rowSource modes
* closes 945 Line chart: ability to hide and make semi-transparent lines
* FilesView: add the "Add new share" ribbon icon and link below the tree
* FilesView: 'refresh' icon
* CI: Docker Buildx services auto list
* CI: Docker Demo Chembl DB
* CI: Docker Buildx target
* CI: Create test_numeric_types table for external_provider_test
* Test: Increase timeout for Dataset Client: Projects: CRUD: Remove
* Test: Increase timeout for Datlas tests
* Integration tests: datetime test
* Line chart: add tooltips in multiaxis mode
* New FontAwesome version
* Fixed css
* Fix bakcground styles, for macromolecules page
* Chem: restoring sim search, add indeces
* code fix for #935 treemap improvments
* Pipe file downloads
* Folder sharing (WIP)
* Ask credentials on file share open
* Minor style changes
* Sketch View: add the reset button
* JS API: Grid.sortByColumns, Grid.sortTypes
* Chem: fix for similarity search
* Grok connect: performance test (select 1)
* Integration tests: refactor external_provider_perf_test
* Fixed AD login
* GH-955: subDir and share methods
* GH-955: js api update
* (Bug) Files: Unable to open context menu for files in explorer view
* Credentials Manager (WIP)
* (Bug) connection.shares.connection == null condition doesn't work (WIP)
* BarChart: support additive functions 'unique', 'missing value count' and 'value count' in stacked mode (WIP)
* Docker: Datagrok ca-certficates
* Fixed files API
* CI: Add Alexey Chopovsky to beta users
* Docker: Remove elastic from datagrok
* Docs: Upgrade docker instructions with certificates
* PowerGrid: bumped up version to 1.1.0
* Wiki: edited test manager documentation
* SequenceTranslator: outdated SMILES tests commented
* Bio: Add performance tests
* OligoBatchCalculator: outdated validation tests commented
* DevTools: update dev panel tests
* ClinicalCase: updated README
* ClinicalCase: corected typo
* ClinicalCase: resized image in README
* Reversed wrong commit
* Public token

## 2022-08-29 Dev build 1.6.4

### Addressed Issues

* com
* Fix default value for select (choice) (remove '' and "" from string when chosen default value)
* Implement setRawData for bool columns
* render mols vertically on scatterplot x axis
* lrellis plot not working in LS
* Fixed the method for default choices
* delete useless prints
* bar chart slider positioning
* Grok Spawner: fix default value for GrokSpawnerAddress
* Make modals resizable
* Legends Resizing
* closes #810 Reordered selected columns
* add started/finished/author to funccall
* Layouts | Save to gallery: newly saved layout should appear in the property panel
* Layouts pane: ability to put the layout on the property panel
* CI: Deploy process documentation
* closes #798 selecting dots in jittered scatterplot
* #828 formula lines modal fixed
* #828 order/hide dialog rewworked
* #828 make modal not resizable by default
* JS API: ui.fileBrowser
* Work in progress.
* #831: HitTriage: work in progress
* Moved InputBase descendants to separate files.
* AppEvents.onInputCreated global event
* Ability to reuse database browser, and specify how results are interpreted
* Table Input: ability to get table by executing a database query
* Docs: Publish version from dev to Docker Hub
* Datlas: queue interaction for Grok Connect (WIP)
* Docker image for packages tests
* add compact view to scatterplot 1.0
* Scatterplot compact view
* JS API: View.statusBarPanels (WIP)
* (Bug) Grok connect: decimal, real, and double tests failed
* Grok connect: refactoring of numeric type tests
* JS: UI Forms harmonization (WIP)
* (Bug) Failed to delete package
* DB Indexes
* Fixed filename
* Fixed CSS
* Added a demo dataset with image urls
* PowerGrid: help on linked images
* JS UI Test: Scripting in Add new column (WIP)
* Functions: show function class instead of null in status bar
* Grok connect: refactoring of numeric type tests
* DataFrame.setTag(...) now returns this.
* Renderers: ability to specify column tags as conditions (WIP)
* adds menu click method
* modal click simulation method adjustment
* Datlas connections leaking investigation (WIP)
* Datlas: remove trailing ; in the row checking query
* JS API: Menu.click method
* #831: HitTriage: work in progress \- refactoring
* closes #800 Form does not work after filtering
* Introduced Entity.package Implemented safe delete New permissions check Derived properties support in ORM
* Introduced Entity.package
* GH-686: TreeView events
* modal placement, dont leave page on backspace
* GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Fixed deadlock in exchangeDevKey method
* CSS fixes
* Introduced Entity.package, reviewed Entities list
* Introduced Entity.package, introduced PackageEntityMixin
* Autogenerated code
* Speed-up FilesView start
* Fixed analyzer warnings
* Fixed package.functions section exception
* closes #791 Can not uncheck "Filter out missing values" after saving
* Fixed wrong naming on package deploy
* Removed vis.js from autostart
* Async external connectors fetch
* Speed-up DG start
* FilesView: fix file appearance bug in the card view
* Introduced Entity.package, bugfixes
* Fixed possible NPE
* JS API: BigInt support
* closes #857 Deselection in the scatter plot
* Fixed author parsing
* Package content validation: improve the up-to-date package bundle check
* Fixed tableInput icons
* closes #859 3D scatter plot: implement saving to .png
* Updated beta users
* Tested renderer selection by multiple tags (works fine)
* (Bug) Integration tests: Fix all northwinds for 'External provider: test all Northwinds' test
* refactor: configuration files
* Docs: fix documentation links
* Public submodule
* Add delete FuncCall api
* PackagesView installed filter
* CSS harmonization
* Introduced EntityType.isPackageEntity field
* Fixes #923: Grid: Order or hide columns: provide "move to the top" / "move to the bottom" buttons
* Viewer-specific tooltips (WIP)
* Min() and Max() are not compatible with date and datetime columns #689
* Fixed tests
* removed debug lines
* #closes 858 Cell as an input for widget
* (Bug) FilesView: dock manager exception on opening
* Setup fixes
* Update doc for proxy
* Better transaction management
* Added debug lines
* Some packages not publishing fix
* Removed debug printout
* Add save as png to all viewers
* (Bug) Connection is cloned on sharing
* Reverted GitHub-854 Provided a fix for Modal dialog to consume mouse events
* Dialog: fixed a bug with the resizing
* Grid: Fixed an issue with bool cell renderer
* dapi.files: fixed a bug by reverting to the lazy initialization. Code cleanup.
* fixed: ColumnsInput: an exception when checkedColumns is not specified
* js-api update
* (Bug) Server push is sent twice
* (Bug) Row.setValues doesn't fire onDataChanged event
* Removed debug lines
* Add exceljs to sources
* DateTime scalar functions harmonization
* Fixed shell test
* PowerPack: build fix
* Input fixes, and new disabled styles
* Fixed scalar functions null behavior
* Fixed logical functions in Chrome
* Introduce DateInput with dayjs support
* Fixed deploy exception
* Fixed label css
* (Bug) ScatterPlot: Top and Bottom legend position cuts labels
* (Bug) Selection doesn't work when jitter activated
* ML | Random Data doesn't work ML | Missing Values Imputation doesn't work
* CSS fix
* Fixed dapi and ml tests
* Possible NPE fixed
* Scatterplot 3d detach method
* dapi.root
* View close fixed
* Updated public
* Test: setUpAll and tearDownAll with groups for integration tests
* Test: Skip tests which fail everytime
* Added JS wrapper for EntityRecord
* Packed categories zoom by filter
* package_id migration tuning
* Fixed user login error
* (Bug) FilesUploader: ?layout parameter doesn't work
* CI: Sofiia Podolskaia in beta users
* Test: Use custom session token for integration tests Load session token for integration tests from environment
  variable GROK_SESSION_TOKEN As a fallback value the old session token unit_test_token will be used
* Refactor: Complete refactoring of deploy folder Deprecate old scripts
* Test: File with environment variables for local tests
* GH-954: expanded accessor for TreeViewGroup
* (Bug) Integration tests: External provider (WIP)
* Docs: Instructions to run datlas tests locally on local custom stand
* Test: Add group for all integration tests Fix integration tests which were failing because of session
* closes #936 pie chart improvments
* Test: Add testDevKey config option for docker stands
* Update background for macromolecules
* Test: fix Missing Values Imputation: k-Nearest Neighbour Imputation
* Fixed Datlas tests
* Reverted incorrect fixes
* GH-965: TreeViewGroup.group now also accepts Element
* closes #961 gray out unused categories instead of not showing them
* (Bug) ScatterPlot: Categories rendering issue
* closes #962 Save as PNG: not all elements are added to image
* closes 790 Reset filter doesn't work on "Filter out missing values"
* add button style
* closes #963 Continuous selection on shift+click
* CI: Add short commit to docker image
* Test: Fix Datlas tests
* Load package before loading tests
* Balloons: add the "X" icon
* closes interfaces generation
* (Bug) When user shares an entity from package package is shared too
* Chem harmonization (WIP)
* closes 935 treemap improvments
* Chem: fix fingerprints
* PC Plot: collaborative filtering
* Test: Increase timeout for datlas tests
* (Bug) Packages permission check doesn't work
* closes 935 with all rowSource modes
* closes 945 Line chart: ability to hide and make semi-transparent lines
* FilesView: add the "Add new share" ribbon icon and link below the tree
* FilesView: 'refresh' icon
* CI: Docker Buildx services auto list
* CI: Docker Demo Chembl DB
* CI: Docker Buildx target
* CI: Create test_numeric_types table for external_provider_test
* Test: Increase timeout for Dataset Client: Projects: CRUD: Remove
* Test: Increase timeout for Datlas tests
* Integration tests: datetime test
* Line chart: add tooltips in multiaxis mode
* New FontAwesome version
* Fixed css
* Fix bakcground styles, for macromolecules page
* Chem: restoring sim search, add indeces
* code fix for #935 treemap improvments
* Pipe file downloads
* Folder sharing (WIP)
* Ask credentials on file share open
* Minor style changes
* Sketch View: add the reset button
* JS API: Grid.sortByColumns, Grid.sortTypes
* Chem: fix for similarity search
* Grok connect: performance test (select 1)
* Integration tests: refactor external_provider_perf_test
* Fixed AD login
* GH-955: subDir and share methods
* GH-955: js api update
* (Bug) Files: Unable to open context menu for files in explorer view
* Credentials Manager (WIP)
* (Bug) connection.shares.connection == null condition doesn't work (WIP)
* BarChart: support additive functions 'unique', 'missing value count' and 'value count' in stacked mode (WIP)
* Bio: Optimized renderer, added sampling
* Docker: Datagrok ca-certficates

## 2022-07-18 Dev build 1.5.1

### Addressed Issues

* resizing on right side
* Fixed admin session
* Ability to save favorites for groups
* Grok connect: stress test
* Wiki: macromolecules: harmonization & cleanup
* #10933legends almost work in any position
* #10933 resdrawing legend while dragging
* closes #10933 delete comments and redraw list on setData
* vertical legend maxheight
* generate comments for ts interfaces
* JS API: add optional tooltip check to InputBase.setTooltip
* (Bug) Packages: settings editor doesn't update package properties
* #718 Custom categorical columns sort order in not reflected on plots after page reload
* GitHub #777: exposed DataConnection query method
* closes #10938 Make modals resizable
* StatusBar improvements
* (Bug) JS: ViewBase.basePath doesn't work
* xamgle.lock
* closes #10935 FlexboxLayout should not add extra divs
* recloses #756 bar chart zoom slider visibility
* Grok Spawner port
* CI: Fix all integration tests for any environment
* #424 fixed x axis log availability
* #756 dont show slider asap viewer created
* Deploy docker images
* Cleanup in beta users
* Documentation
* Public token
* Docker buildx release condition
* Docker Datagrok: Set Grok Connect directory
* Scripting: Support async in JS
* Better styling
* Harmonize functions processing
* Ability to disable NPM
* Fixed analyzer warnings
* Updated public token
* Speedup packages view
* Code format
* Lib bio: notation-converter and web-logo minor corrections

## 2022-07-15 Dev build 1.5.0

### Addressed Issues

* resizing on right side
* Fixed admin session
* Ability to save favorites for groups
* Grok connect: stress test
* Wiki: macromolecules: harmonization & cleanup
* #10933legends almost work in any position
* #10933 resdrawing legend while dragging
* closes #10933 delete comments and redraw list on setData
* vertical legend maxheight
* generate comments for ts interfaces
* JS API: add optional tooltip check to InputBase.setTooltip
* (Bug) Packages: settings editor doesn't update package properties
* #718 Custom categorical columns sort order in not reflected on plots after page reload
* GitHub #777: exposed DataConnection query method
* closes #10938 Make modals resizable
* StatusBar improvements
* (Bug) JS: ViewBase.basePath doesn't work
* xamgle.lock
* closes #10935 FlexboxLayout should not add extra divs
* recloses #756 bar chart zoom slider visibility
* Grok Spawner port
* CI: Fix all integration tests for any environment
* #424 fixed x axis log availability
* #756 dont show slider asap viewer created
* Deploy docker images
* Cleanup in beta users
* Documentation
* Public token
* Docker buildx release condition
* Docker Datagrok: Set Grok Connect directory
* debug

## 2022-07-07 Dev build 1.4.14

### Addressed Issues

* Ability to save favorites for groups
* Packages: improve package sources check
* JS UI Test: Density Plot
* Bio: modified getMolfilesFromSeq function
* CSS adjustments
* #740: barchart tooltips
* Model Property Panel layout
* #740: Peptide Space tooltip

## 2022-07-07 Dev build 1.4.13

### Addressed Issues

* (Bug) Integration tests: ChemMapIdentifiers test fails
* CardView: Ribbon shortcuts (WIP)
* Document viewer properties (WIP)
* Fixed warnings
* Public token
* Macromolecules page
* Ability to save favorites for groups
* Fixed NPE
* OnCurrentObjectChangedEvent
* Updated public token
* CSS improvements
* Wiki: Upload data: improvements
* (Bug) Find and replace: an exception when replacing values
* Fixed NPE on admin login
* Added initial permissions for Admin user
* Fix docker datagrok image
* Added missing migration
* #763: distribution panel label fix
* Bio: getting monomer objects from HELM lib

## 2022-07-05 Dev build 1.4.12

### Addressed Issues

* Document viewer properties (WIP)
* Status Bar
* Popup: ability to show context menu next to the element
* Fix link to dart for Windows
* #735 and #736 property panel search improvments, expand all button
* #758 [Tooltip] \- using the option to sketch form for the tooltip completely breaks its responsiveness for further
  modifications (WIP)
* Ability to set layout on third-party system exported data
* Implement JWT tokens (WIP)
* closes #732 tooltip for large text in small cells
* Properties: "depends on" property
* CSS fix

## 2022-07-04 Dev build 1.4.11

### Addressed Issues

* Document viewer properties (WIP)
* string fields
* JS API: support a selection mode via modifier keys only in Bitset.handleClick
* property grid null dependencies
* Compute: Model Catalog app to appear on the Left Panel #411
* JWT test
* RSA KeyPair: Sign and verify methods
* Implement JWT tokens (WIP)
* Bio: top menu added
* fix semtype
* Fixed #745: Compute: several custom viewers are in conflict
* Fixed analyzer warning

## 2022-07-01 Dev build 1.4.10

### Addressed Issues

* JS UI Test: Pie Chart
* Better units CSS
* fix $ sign

## 2022-07-01 Dev build 1.4.9

### Addressed Issues

* dont explicitly set description with several lines

## 2022-07-01 Dev build 1.4.8

### Addressed Issues

* Properties: an option to reorder rows in column grid for multi-column property editor
* Document viewer properties (WIP)
* Fixed exception on AccordionPane.toJs()
* Compute: Model Catalog app to appear on the Left Panel #411
* property generator reads comments from files
* optimize iteration throulines in same files
* open file from other package
* Fixed parameter type in detectSemanticTypes
* Func subtype fixed
* Wiki: fixed broken links
* Wiki: moved dashboards.md
* Add smiles-to-mol scketcher tests
* Added other conversion tests except molV3
* Added value for converted smarts to ketcher_tests

## 2022-06-29 Dev build 1.4.7

### Addressed Issues

* #741: Peptides Space fix
* Peptide space tests fix and performance benchmark
* Helm: Updating the sizes of the HelmWebEditor
* Chem: activity cliffs tests
* Compute: Model Catalog Groups editor #658
* Exercises: lint fixes

## 2022-06-29 Dev build 1.4.6

### Addressed Issues

* OligoBatchCalculator: fix of using new modification
* Compute: Model Catalog Groups editor #658

## 2022-06-29 Dev build 1.4.5

### Addressed Issues

* Test manager: ability to run tests using url (WIP)
* Exercises: fixes
* (Bug) IgnoreCloseAll doesn't work
* Compute: Script icons support
* JS-API: package-test.js fails to load on autostart
* Compute: Model Catalog app to appear on the Left Panel #411
* Docs: Help documentation
* Docusaurus documentation
* (Bug) Parametrized DataQueries don't work (WIP)
* Fixed typo
* Compute: Model Catalog Groups editor #658
* Fixed routes order

## 2022-06-28 Dev build 1.4.4

### Addressed Issues

* [Calculated columns] \- Min() and Max() are not compatible with date and datetime columns #689
* JS API: add initial sync option to grok.data.linkTables
* JS API: Property: expose category, format, nullable, editable properties
* Viewers: ability to use JS-based viewers in Trellis (WIP)
* Fixed packages view css
* PackagesView: Fixed categories refresh
* closes #639 Scaling Structures on the ScatterPlot
* trellis plot with one cat
* Added more beta users
* Grid: updated some properties
* Work in progress
* Properties: "depends on" property (WIP)
* Retired the obsolete `unitTest` tag
* Chem: optimization of dependencies \- WIP
* Chem: optimization of dependencies \- got rid of the built-in js-api
* Chem: removing obsolete stuff from the core
* D4: rebuilt auto-generated files
* Chem: moved openchemlib-full.js under /common for Chem to load it via "sources"
* Delayed execution for autostart functions by default (+ meta.autostartImmediate)
* Optimized favorites requests
* Package Manager -\- audit fix
* Disable package functions logging
* Optimized DataSources list load
* Package Manager: fixed the styles
* Packages: context menu improvements (WIP)
* Core: Apps won't start from Favourites #632
* Css fixes
* Promote MenuItem to JS
* Compute: fix model URL #574
* Compute: Model Catalog Groups editor #658
* Compute: Model Catalog app to appear on the Left Panel #411
* Closes #733: Grid: allow to resize column header by dragging the border between the header and the first row
* CI: Fix links checker
* HELM: issue #700 Smiles + Molifle representation

## 2022-06-23 Dev build 1.4.3

### Addressed Issues

* Simplified JsConvertible
* Closes #709: Filter panel is not synchronised across different views
* (Bug) Correlation plot: histogram cells are not shown
* (Bug) Correlation plot: in-cell scatter plots are empty
* Correlation Plot: cell format set to 0.00
* Package content validation: up-to-date bundle check
* Packages: notify of successful settings saving
* Dialog: showHeader and showFooter properties
* (Bug) Semantic type detection not performed in the file preview pane
* Bio pckg: Fix VdRegions viewer registering

## 2022-06-21 Dev build 1.4.2

### Addressed Issues

* Grok Spawner: create endpoint for docker build from Dockerfile
* Packages Manager (WIP)
* Packages View: visual improvements
* Grid: html rendering: automatically annotating host div with semantic type (WIP)
* Closes #699: Valid formula line on scatter plot is not shown if one of the scatter plot lines is invalid
* Work in progress.
* Merged in spodolskaya/testtrack (pull request #179)
* Changing the behavior of ESC button in the filter panel #644
* (Bug) Semantic type detection: column with manually assigned sem types are being auto-detected
* Closes #642: Autodetection of Structure columns imported from external files
* Package content validation: up-to-date bundle check
* Closes #719: Viewers: current object won't change after viewer properties are edited
* Chem: v1.0.2
* Chem: #697 Close preserves changes, x\- button recovers initial sketch
* Closes #718: Custom categorical columns sort order in not reflected on plots after page reload
* MLBData: set category to "Bioinformatics"
* Closes #605: [Filter Panel] Search doesn't work in the column list for adding column filter ("+" icon)
* MLB: Buildup regions layout for cdr definition/numbering scheme combinations
* MLB: Buildup regions layout for cdr definition/numbering scheme combinations fix

## 2022-06-17 Dev build 1.4.1

### Addressed Issues

* Packages Manager (WIP)

## 2022-06-17 Dev build 1.4.0

### Addressed Issues

* KeyPair::encryptString/decryptString methods fixed for UTF8 symbols
* Grid: add cell context action
* (Bug) Functions: invalid date format in `Time()`  in some acceptable cases
* Packages Manager (WIP)
* Fixed test
* # 700: Helm fixing trouble with rendering
* F2 no longer edits a cell as it is still used as a column editor hotkey
* (Bug) Data | Text: Semantic type is not detected after table loaded
* JS API: Exposed Column.Aggregate
* Closes #642: Autodetection of Structure columns imported from external files
* Opening editor on double click

## 2022-06-15 Dev build 1.3.5

### Addressed Issues

* MLB: Fix camelCase in DataLoader and setup.cmd

## 2022-06-15 Dev build 1.3.4

### Addressed Issues

* Row number is not clickable after applying saved layout #679
* JS API: convert datetime to dayjs object in DG.toJs

## 2022-06-15 Dev build 1.3.3

### Addressed Issues

* Fixed missing taskbar

## 2022-06-15 Dev build 1.3.2

### Addressed Issues

* Packages: Debug packages permissions fix (WIP)
* Closes #701: ositioning of Hamburger menus with inserted pinned columns
* Packages: update the details section
* fixes batch
* Docker Datagrok new image
* tab fix
* tab
* Packages: add url
* Row number is not clickable after applying saved layout #679
* Added the "showCharts" property that apparently was still used
* Closes #703: Dialog: onCancel does not get invoked if you close it by pressing ESC or by clicking on the 'x' icon
* improve enter behavior
* JS API: convert datetime to dayjs object in DG.toJs
* # 700: Fixing bugs in opening the editor

## 2022-06-13 Dev build 1.3.1

### Addressed Issues

* Formatting: add a missing format alias for floats
* Packages: change package card view for debug versions
* Compute: UX harmonization #180 URL fix
* Docs: Dev deployment
* Run Integration tests locally
* JS API Tests: Functions: tests for all built-in arithmetic/logical/date functions (WIP)
* Packages: Debug packages permissions fix (WIP)
* cheminformatics: refactored up to descriptor-based tools
* Packages: add categories

## 2022-06-12 Dev build 1.3.0

### Addressed Issues

* Packages manager WIP
* Packages manager WIP Some styling
* Closes #655: Grid: No auto new line when editing and '+' available
* Exclude broken integration tests
* Datagrok incapsulated grok connect
* Packages: context menu improvements (WIP)
* # 695: style paramter for range slider
* Packages: Package install fail timeout
* Packages: fixed exception on packages without repository
* PackagesClient.deletePublished
* Scripting: understand extra white spaces in the script header (WIP)
* Packages: Ability to delete package (WIP)
* SelectAll: only selecting rows, not columns
* Filters: fixed the indication of the active/inactive state
* Packages: Debug packages permissions fix (WIP)
* # 698: custom property panel
* Peptides: p-value column format fix
* Viewers: harmonizing properties \- WIP
* Peptides: accordion panels width fix
* Update icons.
* Helm: initial update (test files + detectors)
* Helm: updated README.md
* # 700: Helm: implemented renderer, updated wiki
* Compute: UX harmonization #180
* Peptides: visual fixes
* # 700: Helm: help update
* Peptides: analysis save-load fix
* Peptides: lint fix
* Compute: UX harmonization #180 URL fix
* Updated public index

## 2022-06-08 Dev build 1.2.1

### Addressed Issues

* (Bug) Dev: Beta users and test track (WIP)
* SequenceTranslator: update JS API version
* Fixed Datlas tests
* Changes made to grid row header are not saved with layout #528
* Sidebar icons harmonization
* Error searching non-existent layout #670
* Closes #694: Selected/filtered rows: "Save as column" action
* New viewer icons
* Closes #678: Memory consumption is constantly increasing if a view with filters panel is duplicated
* Closes #87 Charts: migrate to TypeScript

## 2022-06-08 Dev build 1.2.0

### Addressed Issues

* Closes #405: Switching table for a viewer: column selection dropdown/dialog in properties panel is not updated
* Row number is not clickable after applying saved layout #679
* Closes #687: Filters: ESC should toggle all filters on/off
* Closes #608: [Filer panel] Substructure search improvements
* Fix new datagrok image
* GtHub actions: Fix category for Mets
* Chem: fixed empty molecule check
* Make CloseAll aware of ignoreCloseAll
* Wiki: fix documentation using markdownlint
* JS-API: fixed check marks in sketcher menu
* # 691: getColumnSeparator function intial

## 2022-06-07 Dev build 1.1.9

### Addressed Issues

* Chem: changed logic of empty molecule check
* OligoBatchCalculator: fix of not valid data in additional modifications table
* SequenceTranslator: fix of conversion
* SequenceTranslator: fix of sequenceToSmiles tests

## 2022-06-07 Dev build 1.1.8

### Addressed Issues

* Removed debug printout.
* Closes #576: Changing the order of values of categorical columns is not reflected in the filter panel
* Remove COMPOSE_PROFILES

## 2022-06-06 Dev build 1.1.7

### Addressed Issues

* Closes #654: Chem: ability to import mol files V3000
* Closes #640: How to switch Structure Search to categorical filter
* Code cleanup.
* Closes #637: [Filter Panel] \- the history of Structure Search fields
* # 604: Chem: "Use as Substructure Filter" action for molecules \- WIP
* Ability to specify actions for semantic types as JS functions
* Closes #604: Chem: "Use as Substructure Filter" action for molecules
* Dockerfile datagrok multiarch
* The modification of the formula is not propagated to other calculated columns #666
* Buildx Docker images
* Docker buildx multiplatform
* Heatmap: add max columns property instead of hardcoded value
* Fix buildx script
* Resolves #669 JS API: add View.temp object to store auxiliary information
* Docker Buildx disable cache
* Docker Buildx multiplatform: one platform at a time
* Docker datagrok new image
* Docker buildx scan images
* Docker buildx log the tar archive
* Closes #604: Chem: "Use as Substructure Filter" action for molecules \- plenty of fixes and improvements
* Closes #604: Chem: "Use as Substructure Filter" action for molecules \- progress indicator, better interactivity
* Closes #604: Chem: "Use as Substructure Filter" action for molecules \- CSS improvements
* Docker Buildx cache push by flag
* JS API: grok.dapi.groups.getGroupsLookup
* Closes #681: Friendly names for columns
* Added a beta user
* Do not tag latest without push
* Datagrok elasticsearch configuration
* Fixing issues with drugNameMolecule function
* Issue #673: InChi-based sketcher structure query integration
* Chem: filtering using Ketcher sketcher

## 2022-05-27 Dev build 1.1.5

### Addressed Issues

* New Datagrok Image with caching layers
* Create release notes for every release
* Js-Api: fixed bug with function name in chem sketcher
* Js-Api: fixed bug with sketcher function

## 2022-05-27 Dev build 1.1.4

### Addressed Issues

* JS API: add Property options

## 2022-05-26 Dev build 1.1.3

### Addressed Issues

* Closes #648: Viewers: Column selectors: Interactivity: previous column is not set back on mouse leave
* Closes #535: Chem: Structure filter (MarvinJS) cannot be removed from filter panel once added
* Charts: add webpack-cli to dev dependencies
* Simple mode: correct misspellings
* # 618 Charts: Timelines: ability to define multiple "Events" columns (WIP)
* Test for dateTime pattern added
* JS API: JsViewer.columnList property
* Meta fix
* (Bug) Old version layout restoring hides row numbers column
* Extending UNDO mechanism (views closing, axes changes) #600
* Fix Meta package format
* Removed debug lines
* Added missing await
* Meta-package in NPM (WIP)
* Fixed S3 file read
* JS API Tests: Functions: tests for all built-in arithmetic/logical/date functions (WIP)
* Oligo Batch Calculator: editing option on the base modification logic (WIP)
* Add uploading ANARCI (aligned to numbering scheme) antibodies sequences (to filter for an antigen) within position
  names columns to separate tables containing scheme and chain within table name.
* JS-Api: implemented sketcher inplace/external modes
* Usage analysis: use groups instead of users in filters
* # 646: replaced objects with Maps & storing indexes
* JS API Tests: Stats: test each exposed method and property (WIP)

## 2022-05-24 Dev build 1.1.2

### Addressed Issues

* Closes #641: Grid.autoSize(maxWidth, maxHeight)
* Grok Scaler
* Unable to set value to datetime column from JS (WIP)
* Removed unused code
* Moved grok_api to separate library
* Grok Scaler documentation
* ApiTests: fix the build error
* ApiTests: add eslint plugin
* ApiTests: linting

## 2022-05-23 Dev build 1.1.1

### Addressed Issues

* Sidebar integration improvements
* PowerGrid: code cleanup, added a description to the "Global Scale", gave friendly names to renderers.
* Removed sparkline renderers from core
* Property description is rendered as markup. Code cleanup.
* (Bug) File Browser: "Bad state: no element" error when clicking on the empty folder in the tree
* Chem: Code cleanup
* Closes #630: Structure queries
* Closes #523: Sketcher: ability to paste molblock
* Better function signatures
* Complete CI Flow
* Moved CI Flow images to Help
* Ignore drawio backups
* New CI Flow Documentation
* Resolves #617 Charts: Timelines: setting "Color by" to empty results in an exception
* CDM: sample project
* Datlas split to libraries (WIP)

## 2022-05-20 Stable version 1.0.0

### Addressed Issues

* Packages: get compatible npm package versions
* (Bug) connection.close in JdbcDataProvider.execute didn't account for exceptions other than SQLException
* Additional tests for calculated columns
* Push version script: releases branch
* Scatter Plot: Scroll bars visibility
* Docker build script: releases branch
* Lines by equations: Ability to create dashed and dotted lines
* Lines by equations: Ability to remove lines from storage
* Datlas: keep trying to connect to Grok connect
* MultiView css fix
* Viewers: harmonize property names for shapes
* Filters: add wildcard support in search over categories
* Formula lines: Update help
* Work in progress
* Minor code cleanup
* (Bug) In-Vis Filter: Does not affect the regression line
* JS: MultiView: DockView support
* Close button on TabControl
* Convert FunctionView to DockView
* JS: MultiView: DockView support
* Peptides: code cleanup and minor improvements
* Model Catalog improvements
* (Bug) Filters: Uncheck/Check a field in the filter panel that does not update the record properly
* Viewers: In-vis filtering: Implement for Bar Chart
* JKG: lost RxODE
* Fixed analyzer warning
* iframe-embeddable views
* Compute: documentation \- WIP
* Added documentation on caching
* OnViewAdding event
* Compute: documentation
* NPM repositories: new package manager UI
* Scatter Plot: Hide property "formula-lines" from PP
* Remove a utility function
* Fixed css
* (Bug) Multiple FuncCall.onCancel subscriptions
* Additional debug lines
* (Bug) Socket fails if there is ID in DataFrame Param value instead of TableInfo
* Formula Lines: Hide '...' button from Filter field in PP
* Change routing for CVM components in AWS
* (Bug) Packages: dapi.packages.find ignores the include string
* (Bug) Packages: URL field is empty for packages deployed from npm repositories
* README
* Datlas: Add logging to ConnectorsPlugin
* Datals: add tryToReconnectToConnectorServer flag
* In-vis Filter: Implement for Line Chart
* Fix concurrent modification during iteration exception
* CSS Fix
* (Bug) Formula Lines: Lines are clipped on some data
* Packages: Add tests for versions
* JS API: Script.sample, Script.reference, Script.tags, Script.environment
* Chem: similarity analysis and SPE removal
* (Bug) Packages: Infinite detectors loop
* Docker versioning
* Chem: panels removal
* Commit version to the according branch
* Fix condition to create branches in public repository
* Git information inside docker containers
* (Bug) In-vis Filter: Take into account property ShowFilteredOutPoints
* Ability to show input parameters in Function View
* (Bug) MariaDB and MySQL databases throw an exception when trying to open a list of tables
* Separate source maps
* Scatter Plot: Make scrollbars like a histogram and place them on the axes
* More error handling for getNpmPackageVersions
* Check js-api build before docker image build
* Push version: get latest datagrok-api version from npm registry
* Packages: optimised search for compatible versions
* Fixed addNewColumn layout
* Unit tests: exclude and include components from command line
* (Bug) Grok connect: NPE after dev deploy
* Updated dev environment v2 for Linux users
* Packages: allow subfolders in the directories reserved for entities
* (Bug) Grok Connect: Athena: java.lang.NullPointerException at java.util.Hashtable.put(Hashtable.java:460)
* JS API: add Property options
* Datlas: remove timer for reconnecting to Grok connect
* Grok connect: improve SettingsManager for tests and GrokConnectShell
* (Bug) Formula Lines: Sometimes tooltips doesn't appear
* Wiki: Access \- WIP
* Package files caching
* Scatter Plot: Marker coding for Split-By-Marker
* JS API: Improved code style (marker coding, formula lines)
* (Bug) Scatter Plot: Y axis shows weird values
* Small help addition
* Scatter Plot: Get rid of code duplication (axes and grid lines)
* (Bug) Bar chart: When "Axis Type" = "logarithmic", numbers on the axis get merged
* (Bug) Grid: categorical color-coding cannot be switched on/off without opening 'Edit' menu in some cases
* Functions: help-url support for scripts in Function View
* Ability to disable events in InputBase
* (Bug) Calculated column with numeric type and empty formula leaves no formula tag
* Ability to show multiple viewers in Function View
* Function View: Show table inputs if there are viewers set
* small css adjustment
* Function View: Adjust table heights
* Grid: add changePropertyPanel property
* Fixed analyzer warnings
* Update public
* Scatter plot: regression line: improve precision for the regression coefficient
* Scatter plot: improve axis formatting (".00" postfixes are often unnecessary)
* Scatter Plot: Rename and move property for MinMax labels visibility
* (Bug) Trellis Plot: full screen mode for inner viewers doesn't work when trellis is zoomed
* Renamed .d to .dart to fix minified code
* CSS fix
* Grok Connect: maven cache in docker image
* Viewers: Unification of properties (alpha/opacity/opaque)
* (Bug) JS API: null reference in JsViewer.dataframe
* Scatter Plot: an option to follow global filter
* Added an SQLite demo file
* SQLite: improved documentation
* Bar chart: improve and refactor collaborative filtering
* Grok connect: NPE when settings in schemas, schema, and query_table_sql are null
* Formula lines: Fixed name conflict
* NPM repositories: add check for deprecated versions
* NPM repositories: improve package search for given scope
* Terraform configuration for CVM
* Function View: Viewer positioning tags
* Function View: Ability to show grid
* Release notes: ignore GitHub Actions commits
* Added isolatesCount to GROK_PARAMETERS
* Added pool settings to GROK_PARAMETERS
* (Bug) Grid: "Show row header" property does not have any effect
* Skip test: ddt.calculated_column.metadata presence after formula application
* (Bug) Tree Map: After selecting column for "Size" property, column selectors on the viewer are duplicated
* Chem: Support filters from packages
* Push version using the version in package.json file
* (Bug) Formula Lines: "=" in the column names breaks the formula lines
* NPM repositories: add registry property
* Actions panel function names fix
* Added build all batch script
* Fixed name conflict (opacity)
* (Bug) Formula Lines: Browser crashes om some analysis
* ElasticSearch log4j vulnerability CVE-2021-44228
* Push version: generate release notes for minor release
* (Bug) ListInput: caption is not shown
* Filters: ability to add multiple filters (for the specified selected package filter) for selected columns
* Filters: ability to remove all filters (via the popup menu)
* Scripting: in-place definition for environments
* (Bug) Viewer Filter: Scatter Plot doesn't refresh after Viewer filter changed
* H2O: log4j CVE-2021-44228 remediation
* Css fixes
* CSS Harmonization
* Changed default run section to condensed mode
* Formula Lines: Name harmonization \- core side (doesn't affect users)
* Formula Lines: Name harmonization \- js api side (affects users). Change property name: equation -> formula. Both
  options are valid, but equation is obsolete.
* (Bug) Filters: closing the viewer does not properly detach it (event subscriptions are leaking)
* (Bug) ColumnComboBox: allowDrop field is ignored
* (Bug) Widget.die() does not kill all descendants
* Filters: ability to drag-and-drop columns to the filter group
* (Bug) Formula Lines: No lines if plot was created from the script annotation (Pmax)
* Propertly detaching JS-based filters
* Minor fixes and improvements
* (Bug) Scatter Plot: Wrong scrollbar length when canvas resized
* JS API: Filters (WIP)
* Bar Chart: Wrong scrollbar length when canvas resized
* JS API: ability to interop with Dart lists (no deep copy)
* Harmonize interop
* Minor code cleanup.
* JS API: getRawData should return real buffers (currently creates a copy)
* (Bug) Bar Chart: Exception when filtered table rows = 0
* SketcherFilter: provided filterType (will be retired soon, anyway)
* Seamless data loading to Function View
* Function View improvements
* Markers: Fixed bug with last marker in List
* Markers: Add new marker for outliers
* Python base image: remove boost artifacts from image
* Package repositories: unique scope check
* Datagrok: disable update certificates by environment variable
* Function View: markup improvements
* Formula Lines: Property harmonization (opacity)
* JS API: Moved properties (column.colors => column.meta.colors, column.markers => column.meta.markers)
* Scripting: global named conda environments
* JS API: Users API
* Formula Lines: Endless lines
* Result Header CSS Fix
* Formula Lines: Minor changes
* Minor code cleanup, switched to the thin barbell look
* Fixed a silly division by zero error
* Fixed filtering and add minor fixes
* _JsBasedHandlers resulting tables names fix
* Bar Chart: Fixed "compact()" call on nullable "stackColAggr"
* Removed parentheses
* (Bug) Line Chart: Fixed bug with "Y Global Scale"
* (Bug) Core: Sensitivity analysis won't work
* Line Chart: Added global Y-axis
* Line Chart: Added title for global Y-axis
* Octave: vectorization
* (Bug) Line Chart: X-selector overlaps the chart when the x-axis is auto-hidden
* (Bug) Line Chart: No X-axis in Multi Axis mode (if there are many lines)
* Line Chart: Show chevrons only on mouse hover
* Line Chart: Ability to set selectors legend position
* JS API: improve the registration of async functions (WIP)
* (Bug) Grid: column popups: popup container left handing in the DOM tree after a popup is closed
* Filters API: passing the requester filter up the event chain
* (Bug) ui.wait() adds default 400x300 size to the container
* Filters in the column property panel
* Grid: column quick panel: ability to add as filter (WIP)
* Bar chart: checkboxes "select" and "filter" in the context menu
* (Bug) Property panel: columns: clicking on the column does not change current object
* Added a comment on comments
* ClientPackageFunc: ability to run function synchronously once the package is loaded
* (Bug) Scatter Plot: programmatic zoom causes an error
* Untyped events since they clearly were not `Stream<String>`
* Bar chart: add categoryValueWidth property (WIP)
* Chem: excludin duplicated features
* (Bug) Scatter Plot: onZoomed returns objects of different data types when the plot is zoomed manually and
  programmaticaly
* Provide Area \- interface for Rect, Polygon, etc.
* Line Chart: Auto margin for yAxisTitle (get rid of yGlobalMarginAxisLeft)
* (Bug) Line Chart: The top marker of each line is cut off
* Line Chart: Changed close-icon style of lines
* Line Chart: Place close buttons to the right of the y-selectors in MultiAxis mode
* Minor CSS improvements
* Wiki harmonization \- WIP
* Name harmonization: API Samples, API Tests
* Set some packages as beta.
* Packages: remove unneeded package detectors
* Dialog.create: made options really optionable :)
* Axes: Ability to display vertical axis on the right side of the chart
* Line Chart: 2 different y-axes when lines count = 2
* Documentation improvements
* Chem: move semantic type detectors to the Chem package
* (Bug) Chem: some panels (drug likeness, toxicity, etc) do not work with molblocks-encoded molecules
* Unit conversion \- WIP
* Chem: get rid of the Molecule-specific renderer tricks in the core
* Line Chart: Special y-selector visual mode for 1 and 2 lines
* Axes: Add more y-tickmarks (need for Line Chart)
* Reusing DruglikenessPredictor
* (Bug) Opening an SDF file results in writing a big binary array to Datagrok's console
* Made the accordion pane's left margin smaller
* (Bug) Info panel invocations end up cluttering Datagrok console
* Build fix
* JS API: Cell.value setter
* (Bug) Chem: double-clicking on a structure: move to the Chem package
* (Bug) Chem: double-clicking on a structure opens an empty sketcher
* (Bug) Chem: Molfile widget does not work with molblocks
* Implemented #243: Chem: sketcher: "Copy as SMILES" and "Copy as MOLBLOCK" context commands
* CSS fixes
* Code cleanup
* Chem: fixed two memory leaks (molecules not being disposed)
* Chem: SDF importer: fixed a bug with 100 first rows being ignored
* Chem: removed unnecessary imports
* JS API: Func.find: ability to search by metadata
* Dart Sketcher: fixed an issue with smiles = null
* Chem: refactoring cell rendering \- WIP
* Formula Lines: Editor
* grok.shell tests fix
* (Bug) Octave: graphics won't work
* DevOps documentation
* Modified field length
* Update the dialog design for developer keys
* Help failed tests troubleshooting
* (Bug) Bar chart: scroll the page when BC scrolls are invisible
* H2O: Grok Helper URL prefix
* (Bug) DataFrame.onCellValue event not fired when a value is deleted in a spreadsheet by pressing Del or Backspace
* Viewer in accordion has zero height
* Bar chart: pass information about is click was on the header
* Bar Chart: maxBarHeight property
* Bar Chart: verticalAlign property
* (Bug) Bar Chart: initial render cuts category names
* (Bug) JS API: Legend.create doesn't return an instance of DG.Legend class
* Add CodeMirror to sources
* Change release-history.md during release
* Updated beta_users.csv
* Bump elasticsearch from 6.4.2 to 7.16.2 Elasticsearch release to upgrade Apache Log4j2
* JS-API: Packages test framework
* Bar chart: add filtering by clicking a category label in the stacked BC
* Bar chart: fixed bar border line
* (Bug) JS API: "Legend.column =" causes an exception
* JKG: Allow gnuplot run on different platforms
* Fix docker build for datagrok
* (Bug) Grid: virtual column of type bool renders as text
* 9913 displaying scalar query outputparams
* 9913 space by arrow
* DataFrame: editable virtual columns
* JS-API: Exposed DataFrame.fromByteArray and FileSource.readAsDataFrame methods
* JS API: updated js-api-common files
* (Bug) Pivot table icon is missing
* JS API: readAsDataFrames renamed to readBinaryDataFrames
* Bar chart: click not on the bar to clear the filter
* PowerPack: added tests for widgets
* ApiTests: added a test for package files
* (Bug) Trellis Plot | Inner Viewer Settings throws an exception
* ParamCommand: default item comands ("open" for project, "details" for users, etc). Support for double-click
* Property panel: decreased the margin for nested accordions
* (Bug) Grid: renderer type is not automatically changed when it is set via column.tags['cell-renderer']
* Resolves #266: Chem: molecule column property panel
* Dialogs: add missing wiki links
* JS API: add Viewer.helpUrl parameter
* Introduce Func.source constants
* (Bug) grok.dapi.files.delete(): NoSuchMethodError
* 9913 dataFame
* API: added getting and saving event types
* Updated demo/bio peptides data
* Added ability to set code in ScriptView constructor
* ViewFactories ScriptView construction enhancement
* Move all demo assets to packages
* (Bug)  JS-based file handlers are not invoked when you double-click on a file in file browser
* # 277 Chem: Substructure filter: fixed collaborative filtering
* # 277 Chem: Substructure filter: picking up current filter
* (Bug) DAPI: Filter returns incorrect query
* Fixed test
* (Bug) Grid: column tools: filter: add: there should be no other filters in the filter group
* Simplified the creation of the popup filter
* # 277 Chem: Substructure filter: better integration for the filter group / popup filter / property panel filter
* Formula Lines: Fixed typo
* Update the beta users list
* Fix a css rule
* NPM repositories: package data caching
* Amazon S3 Adapter: Anonymous access
* Packages: Ability to store content in S3
* Code cleanup and harmonization
* Moved files to folders
* Removed default handlers for Function and Script
* JS API: GridCellRenderer: ability to render values (without GridCell / GridCellStyle)
* Chem: Similarity & Diversity: switched to the renderMolecule() method
* Add a grayed-out text color constant
* DataQueries: Support scalar results in setResult method (WIP)
* UI: Condensed form, units after labels
* added ScriptView.fromParams constructor
* Always create NPM repository
* JS API: Projects API
* Updated public token
* Remove redundant deployDataPath
* Enabled Chem | Descriptors
* Embedding viewer should embed Guest User session
* CI: Vulnerability scan
* (Bug) Filters: the panel is reset when it is closed and re-opened
* JS API: Get last error state (WIP)
* (Bug) Histogram: Opacity for bands has no effect
* (Bug) Line Chart: Exception while activating multiAxis & overviewType at the same time
* (Bug) Line Chart: Previous renders of the chart are visible on the canvas
* (Bug) Line Chart: X axis overlays on lines when { overviewType: "Line Chart" }
* AppEvents.lastError: simplified the code and fixed a bug
* Updated package-locks.
* (Bug) Filters: show header property does not work
* JS API: Filters: fine-grained API (WIP)
* Filters: ability to specify filtering criteria during construction
* JS API: JS Filters: ability to provide custom caption
* JS API: Help & Example for ui.colorInput()
* (Bug) Line Chart: Selectors indents (when count = 2)
* (Bug) Scatter plot: Erroneous ticks when zooming out
* Fixed analyzer error
* JS API: Added additional useful FLines methods
* Usage analysis: add Function errors tab
* Chem: minor fixes
* Func.find: fixed a bug with metadata search
* Func.find: ability to search for the result type/semantic type
* Implemented #138: Chem: Text search for molecules in the sketcher
* Chem: refactoring
* JS API: Added ColumnList.getUnusedName()
* (Bug) FLines: Constant lines do not take min-max into account
* credential manager prototype
* 10269 secret manager prototype
* JnJ: monitoring resources: RAM, CPU, HDD
* Implemented #304: Chem: Sketcher: Add to favorites
* (Bug) Filters: free-text input: emptying the input box throws an error
* Tooltip: ability to show multiple tooltips at once
* Filters: Histograms: Tooltips while dragging
* Css style fix for the checkbox in filters
* JS API: add vertical RangeSlider
* (Bug) Accessing .type in onViewerAdded event breaks filter
* Canvas viewers: fixed fps calculation
* Scatter Plot: selected rows should be rendered on top
* Line Chart: logarithmic axes (WIP)
* (Bug) Scatter Plot: area selector: null values are getting selected
* A quick way to move columns to the beginning/end of the spreadsheet
* Box Plot: support for log scale
* Chem: Filter won't appear
* Fixed analyzer warings
* Fixed wrong method
* (Bug) Demo files connection is missing
* Minor CSS fixes
* Fixed typo
* Upgrade H2O
* (Bug) Filters are not applied in some cases
* (Bug) Function View: Add to Workspace button won't work
* (Bug) Function View: Dataframes selection from a file is buggy
* (Bug) Function View: Dataframes selection from a file is buggy
* (Bug) Scatter Plot: column selectors should be enabled in the invalid state
* (Bug) Form: Structure in form viewer is not rendered after page reload
* (Bug) MoleculeInput: caption is not shown
* JS API: Inputs: validation (WIP)
* Public token
* Fixed grok_shared build
* (Bug) Filters: After reopening filters, all columns appear deactivated in the "Select columns" dialog
* (Bug) Filters: Filters for multi-value columns are not saved after reopening the Filters
* (Bug) Filters: ​double-click on checkbox: weird characters are visible
* (Bug) Scatter Plot: Custom linear coloring is not propagated
* NPM repositories: add an option to skip installation of deployed packages during repository publication
* Push version from package.json for master branch
* syft local image
* DockManager_Dock now attaches viewer to view
* Fixed #320: Filtering: 'Select columns' dialog for filter panel is showing currently selected columns incorrectly in
  some cases
* Fixed #319: Disabling/not disabling filtering behaviour on closing filter panel is inconsistent in some cases
* (Bug) JS API: incorrect constant in DG.DOCK_TYPE.TOP
* Fix datagro build
* Added Viewer_Remove_From_View function
* js-api-common update
* Scripting: Better handle missing output variables (R, Python)
* Package description: add markdown rendering
* JS API: Dialog.getOpenDialogs()
* JS API: Dialog.input(caption): InputBase
* JS API: Grid: ability to get cell back color
* (Bug) Color coding: the editing dialog shows the black color for categories which color is set via hex codes
* (Bug) Line Chart: rightmost line segment is not always drawn
* New demo dataset: dose-response
* New sample script: charts-in-cells
* GridColumnMeta
* (Bug) Chem: double-clicking on an empty structure does not open a sketcher
* (Bug) Grid: an exception when resizing the window
* Wiki: Doc on running Dart tests locally
* Test file fixed precision
* (Bug) Grid scrolls to first row after adding/deleting row
* Grid: ability to easily add rows to the end
* # 323: MultiForm: WIP
* Fixed null exception
* (Bug) Parameter subquery is empty
* Harmonize viewer title
* Property panel: "Distributions" pane (WIP)
* Core: delete the service jobs before the repository deletion
* bug fixed fro smiles string in example
* Formula lines improvements
* Update repository jobs: change user permissions (WIP)
* Error status style update (WIP)
* (Bug) BoxPlot Y axis labels overlapping with axis column name
* Harmonize PackagesView (WIP)
* Filters: header should stay on top when you scroll down
* (Bug) Input too large for RSA cipher exception during package deploy
* Avoid extra data copy
* update ref to head on public
* Porting Ketcher reading molfile in React to package (WIP)
* DB: migrate test and demo DB to docker
* VirtualItemView.refreshItem
* Docking Manager: improve auto-docking algorithm when the initial split is vertical
* Improved the description of the connection string in the "add connection" dialog.
* Security: WIP
* Docs: add local installation issue.
* (Bug) Color coding: missing values lose the reserved color when the rules for categorical color-coding are set
* Add details about this can be applied only for Debian-based systems.
* Added Neptune logo
* Docs: create documentation for CVM connection from Local Dart stand
* Strong mode static checking (WIP)
* (Bug) Scripting: Reset variables before script run (Python, R)
* Fix analyser warnings
* (Bug) Renv creates environments directories with nested UUIDs
* (Bug) Several functions broken in calculated columns
* Checkbox styling improvement
* (Bug) Chem: moleculeInput won't set smiles
* Expose additional map data structure to save meta information along with project entity #331
* Add additional options to configure the view loading behaviour to project.open() method #332
* Viewers: ability to specify default axis type (linear / logarithmic)
* PC Plot: log scaling for Y axis
* Grid: spacebar to toggle row selection
* Made analysis options less strict.
* Scatter Plot: ability to invert axes
* Bar Chart: fixed an issue with the initial viewport position
* Expose id property for tableView entitity #328
* Core: update OpenChemLib
* (Bug) Exception on reading credentials for package
* Webstorm configs
* Ability to set Function as a default value
* WebStrom default configurations
* Viewers: legend: ability to reset"filter by category"
* (Bug) Newly created user failed to login
* Fixed #344: Miss configured in-visualisation filter kills the browser session
* Fixed the initial position of the filter
* Grid: Excel-style column resizing
* Grid: column rearrangement: fixed the "'offsetParent' was called on null" bug
* Grid: ability to reorder multiple columns at once
* Renamed images
* Implement '!=' and 'NOT IN' operators for matching rows
* Fixed #348: HTML cells \- selection & scrolling issues
* EacapeSequences: R
* Updated help
* Line Chart: "Invert X Axis" property
* Exposed grok.events.onResetFilterRequest
* Scatter Plot: highlight the category on axis on mouse-over
* Sharing file connection must share credentials (WIP)
* Made debug configuration the default one
* Grid tooltip should be more discreet
* Add IN operator to AddNewColumn
* Box Plot: "inverse Y axis" property
* Ability to rename tab of docked viewer
* (Bug) applyLayout drops columns if table property set
* (Bug) Unable to switch scatterplot axis to logarithmic mode
* Chem: working with mol
* (Bug) AddNewColumn: Vectorization doesn't work
* Chem: hiding old realizations of similarity
* (Bug) Unable to delete project
* Grok Compute gunicorn workers
* PowerGrid: initial update
* Chem: hide fasr descriptors calculation
* Closes #372: JS API: expose CsvExportOptions
* Closes #350: Add option on toCsv() function to split QNum into two columns
* (Bug) Unable to call a function from AddNewColumn formula
* Missing table attribute in layout for "active table" #345
* DG.Utils.download(filename, content, contentType)
* Closes #369: Persisting tooltip form state
* Fixed #357: Free-text filter doesn't work for columns which names contain more than one word
* (Bug) Viewers: 'Edit Viewer Tooltip..." disregards the table a viewer is bound to
* Closes #374: Grid: rendering dataframe values as HTML #
* (Bug) toCsv() throws an exception
* User status interop
* (Bug) Filter indicator is hidden after restoring the layout
* Closes #340: Expose Grid renderers for DG native data types
* Closes #383: Items are unexpectedly filtered out after new rows were added
* (Bug) Query View: "Add results to workspace" button does not work
* Closes #385: Filters saved on older DG versions are not restored
* Fixed the type annotation
* Fixed the friendly name
* Closes #393: GridCellRenderer: ability to handle mouse input
* Closes #393: GridCellRenderer: ability to handle mouse input \- added onMouseLeave
* VPN access to dev resources
* Function View convert to TS WIP
* Closes #393: GridCellRenderer: ability to handle mouse input \- fixed onMouseLeave
* Closes #395: GridCellRenderer: ability to render HTML elements
* VPN with credentials
* Merged in kdoncov/grok_connect_vulnerability_fixes (pull request #163)
* Remove deploy folder
* Unable to set value to datetime column from JS
* # 396: Chem: support SMARTS for db querying
* Elasticsearch config
* Chem: polishing
* Closes #401: Vertical Axis: incorrect label density for custom renderers (i.e., Molecules)
* Chem: bumped version up to 0.51
* `LruCache<Key, Value>`: added generic arguments and strongly-typed it
* Simplified code to avoid warnings
* Chem rendering: strongly-typed, simplified
* Chem: work in progress
* Fixed positions descriptions
* Elasticsearch commit
* Closes #423: Grid: JS-based sparklines
* Closes #422: PowerGrid: Sparklines: Sparklines
* Closes #435: ui.columnsInput: a way to specify a subset of columns to choose from
* Closes #436: Grid JS API: ability to add virtual column
* Added new functions: Utils.openFile and Utils.openFileBytes
* BinaryImageCellRenderer: work in progress
* Fixed #403: Scatter plot: wrong category name is shown on highlighting
* Fixed #404: Box plot: axis labels are not updated when axis is inverted or changed to log scale
* public
* Ability to get all packages tests
* Grid: better renderer isolation: an exception in the JS renderer should not affect the whole grid
* Logging: ability to report a given exception only once per session
* GridColumn: fixed serialization of JS-based settings
* Closes #442: PowerGrid: Sparklines: an option to normalize values
* Closes $443: PowerGrid: PercentCompleted cell renderer
* Allow minor builds in push_version script
* PGAdmin connect through SSH tunnel
* # 424: Blocking option (graying it out) to switch to log scale for columns with negative values \- WIP
* Allow user to rename the viewers from the arranged tab #341
* ignore IDE files
* removed duplicate
* # 424: Blocking option (graying it out) to switch to log scale for columns with negative values
* New method: ListUtils.move(list, from, to)
* Grid: manual row reordering
* Closes #455: Ability to define order of categories
* Grid: fixed an issue with the exception when row header is not visible
* Default value computing error handling
* Closes #464: ability to rearrange columns
* Closes #465: Column selectors: search field is invisible (above the screen) when the number of columns is big
* Grid: Order columns: ability to show/hide all selected columns at once
* Grid: Order or Hide: ability to reorder all selected columns at once
* # 466: ComputationView \- work in progress
* Added banner
* (Bug) DataFrame: Column.init(...) does not increment dataframe version
* (Bug) Grid: html cells: rendered cells should be invalidated when table content changes
* # 470: JS API: JS-based custom inputs \- work in progress
* Miss configured in-visualisation filter kills the browser session #344
* New VPN server
* Temporarily rolled back JsView.name-related changes
* Dataframe formula lines are not shown on scatter plot, if scatterplot's table was changed in properties panel #375
* (Bug) Saving Connection doesn't save credentials
* APP_DATA variable in package deploy
* (Bug) NPM repositories: package manager opens with an exception
* Add missing import statement
* Method to check current user permissions
* Fixed typing
* (Bug) Registry is undefined when a package repository update job is launched
* Added Alexander to the list of users.
* # 480: Ordering categories: SMILES should be rendered as molecules
* JS-API: added notify parameter to BitSet.init
* Revert grok_server and connectors changes
* HttpDataSource -\- dapi.functions support
* Fixed wrong parameter name
* # 466 Computation View \- Integration with FunctionView
* JS-API: exposed setBufferSilent method
* JS-API: setBuffer method
* Closes #488: Positioning of Grid cells in HTML rendering mode
* Menu capitalization fixes.
* Closes #489: Grid: "Current Column | Sort | Custom..." should open the "Order or Hide Columns" dialog
* bug: duplicated event firing
* Closes #490: Custom order of categories should be saved in the view layout
* Wrong data is shown in viewers using calculated column after calculated column's formula was updated #434
* bug: changing cell-double-click features
* (Bug) Wrong output of getCategoryIndex with caseSensitive = false
* Datlas: use reader/writer login in Datagrok system connection, if reader login isn't available
* (Bug) Previously deployed packages sometimes get deleted during repository publication
* Ability to use Func options for grouping
* JS: Ability to disable events in InputBase
* Exposed DG.Color.fromHtml method in JsApi
* JS-API: add a default parameter value for FileSource.list()
* (Bug) NPM repositories: missing proxy settings
* Closes #433: Data exported to csv is incorrect in 'sign' columns if there is more than one such column and there are
  some nulls
* Closes #526: [LineChart] \- axes do not scale to plotted points but to max values of plotted columns
* Closes #367: Persisting visualizations zoom when saving analysis \- proper handling of the "no action" zoom type
* Updated help links
* Closes #367: Persisting visualizations zoom when saving analysis
* Packages: skip the build step for packages that were published along with build outputs
* Fixed an issue with the barchart.
* Closes #530: Grid: do not include units in the column header
* Scatter Plot: Animated zoom to filter
* Ability to save new user without permissions
* Top Menu: change "mousedown" event to "click" for clicks on menu items
* (Bug) Package content validation: one-line functions get annotation warnings
* Closes #532: Show only filtered category items on legends
* # 342: (Only) show filtered category items on scatterplot axes
* Caching Column.getUsedCategories for performance reasons
* JS: bySemTypeAll method
* DataSourceCardView categories (WIP)
* setup for connect refined
* Sketching: size up
* (Bug) Can't set and get ViewBase URL
* Chem: sketching issues with substructure search
* (Bug) Color coding: Selected palette in "Colors" panel differs from the actual colors
* (Bug) Modal window positioning bugs
* (Bug) Color Coding: Color does not change for boolean columns
* Closes #547: Chem: SubstructureFilter: substructure is not saved properly
* Update the linear color-coding section
* Add .item to grid context menu args
* Added Widget.isDetached property
* Closes #524: Filters: highlighting stays after you close the filter group
* Chem: substrucure search queries fixed
* Molecule input: code cleanup
* JS API: additional options for Menu.items
* Documentation \- WIP
* Zeno integration (WIP)
* (Bug) Grid: ctrl+click no longer inverts selection
* (Bug) Error when closing a view
* ScatterPlot: made the "viewport" property hidden
* Functions: dateTimeFromUnixTimestamp (WIP)
* Closes #367: Persisting visualizations zoom when saving analysis Closes #552: Zoom Sliders get reset on the Scatter
  Plot (Zoom And Filter: no action, Axes Follow Filter: false) on opening the filter panel with some filtering applied
* Color-coding: add a check for the currently used color-coding type to the column context menu
* git-secret for credentials
* git-secret documentation
* git-secret
* (Bug) Fix Demo deploy
* GitHub #566: CsvImportOptions
* (Bug) Packages: missing version check during package upload
* (Bug) JS API: JS-based viewers that are instantiated directly are not linked to the Dart viewer
* Closes #545: Including columns headers when copying data from the table (ctrl + c / cmd +c )
* ddt folder
* restore lines
* Closes #546: [ScatterPlot] \- the config panel is being closed when trying to change some values in dropdowns
* Closes #507: [Filter panel \- Structures] Turning on / off doesn't work for structure columns in the filter panel
* [Formula Lines] Columns with parentheses in its names used in Bands' formulas breaks Formula Lines rendering #542
* HierarchyView: minor UI improvements
* PowerGrid help \- WIP
* Datagrok build context
* (Bug) Color Coding: After enabling Conditional for numeric column, the colors are displayed as Linear
* Closes #588: DataFrame: strongly-typed columns #588
* Closes #589: Strongly-typed InputBase.value
* Share git secret to aparamonov
* #594: the 'Features' dialog in the predictive modeling view shows incorrect count of checked columns when opened for the second time
* TreeView harmonization
* (Bug) Connections deploy doesn't deploy connection credentials
* Folder content preview wiki: fixed the title
* FilterGroup.add: fixed the signature.
* Added shell.tv: TableView
* Grok connect: test numeric types
* Datlas split to libraries (WIP)
* Closes #527: Grid options don't return default values
* Viewers: Transformations (WIP)
* JS API: fix type annotation for DataFrame.getSortedOrder
* Meta-package in NPM (WIP)
* JS API: better typing \- WIP
* Package guidelines (WIP)
* JS API: GridColumn.renderer
* Sparklines: renamed "Normalize" to "Global Scale"
* #612: Compute: Add "open" option for DataFrameInput
