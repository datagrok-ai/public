<!-- TITLE: Release History -->
<!-- SUBTITLE: -->

# 2021-12-22 Dev build 0.107.0

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.107.0`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) NPM repositories: add registry property 
* Formula Lines: Property harmonization (opacity) 
* JS API: Moved properties (column.colors => column.meta.colors, column.markers => column.meta.markers) 
* Scripting: global named conda environments (WIP)
* JS API: Users API (WIP)
* CSS Fix 


# 2021-12-21 Dev build 0.106.0

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.106.0`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Function View: markup improvements 
* Datagrok: disable update certificates by environment variable 


# 2021-12-21 Dev build 0.105.0

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.105.0`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) ListInput: caption is not shown 
* Filters: ability to add multiple filters (for the specified selected package filter) for selected columns 
* Filters: ability to remove all filters (via the popup menu) 
* (Bug) Formula Lines: Browser crashes om some analysis 
* Scripting: in-place definition for environments 
* (Bug) Viewer Filter: Scatter Plot doesn't refresh after Viewer filter changed 
* H2O: log4j CVE-2021-44228 remediation 
* Css fixes 
* CSS Harmonization 
* Changed default run section to condensed mode 
* Fixed analyzer warning 
* Formula Lines: Name harmonization \- core side (doesn't affect users) 
* Formula Lines: Name harmonization \- js api side (affects users). Change property name: equation -> formula. Both options are valid, but equation is obsolete. 
* (Bug) ColumnComboBox: allowDrop field is ignored 
* (Bug) Filters: closing the viewer does not properly detach it (event subscriptions are leaking) 
* (Bug) Widget.die() does not kill all descendants 
* Filters: ability to drag-and-drop columns to the filter group 
* (Bug) Formula Lines: No lines if plot was created from the script annotation (Pmax) 
* Minor fixes and improvements 
* Propertly detaching JS-based filters 
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


# 2021-12-16 Dev build 0.104.0

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.104.0`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Push version: generate release notes for minor release 
* ElasticSearch log4j vulnerability CVE-2021-44228 


# 2021-12-16 Dev build 0.103.0

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.103.0`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Skip test: ddt.calculated_column.metadata presence after formula application 
* (Bug) Tree Map: After selecting column for "Size" property, column selectors on the viewer are duplicated 
* Chem: Support filters from packages (WIP)
* (Bug) Grid: categorical color-coding cannot be switched on/off without opening 'Edit' menu in some cases 
* Push version using the version in package.json file 
* (Bug) Formula Lines: "=" in the column names breaks the formula lines 
* (Bug) NPM repositories: add registry property (WIP)
* Actions panel function names fix 
* Added build all batch script 
* Fixed name conflict (opacity) 
* SQLite closes #47: fixed output error 
* (Bug) Formula Lines: Browser crashes om some analysis 


# 2021-12-14 Dev build 0.102.0

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.102.0`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Function View: Ability to show grid 
* NPM repositories: new package manager UI (WIP)
* Compute: Vmax: Pmax: rounding values, convenient naming objects, transposing tables 
* Compute: export improved 
* Compute: made OutSel buttons flexible 
* (Bug) Grid: "Show row header" property does not have any effect 
* Compute: removed buttons from OutSel 
* Compute: Vmax: remove titles 


# 2021-12-14 Dev build 0.101.0

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.101.0`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Function View: Ability to show grid 
* Compute: Vmax: transposing tables(WIP) 
* Release notes: ignore GitHub Actions commits 
* Added isolatesCount to GROK_PARAMETERS 
* Added pool settings to GROK_PARAMETERS 
* Peptides #93: Add spiral projection viewer. 


# 2021-12-14 Dev build 0.100.0

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.100.0`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* GitHub Actions: Update js-api/package-lock.json 
* Compute: Vmax: consistent capitalization, Update regression line in Input tab after selecting new outliers 
* Timelines: improve column selection heuristics 
* SQLite #47: added missing file 
* Compute: Vmax: flux decay multiaxis linechart 
* Timelines: selection color is missing 
* Scatter Plot: an option to follow global filter 
* Viewers: Unification of properties (alpha/opacity/opaque) 
* Added an SQLite demo file 
* SQLite: improved documentation 
* (Bug) datagrok-tools: `grok api` creates invalid function declarations 
* GitHub Actions: Update tools/package-lock.json 
* Move scripts into DemoPackages 
* Compute: Vmax and Pmax enhancements 
* GitHub Actions: Update packages/Compute/package-lock.json 
* Grok connect: NPE when settings in schemas, schema, and query_table_sql are null  
* Formula lines: Fixed name conflict 
* NPM repositories: add check for deprecated versions 
* NPM repositories: improve package search for given scope (WIP)
* Compute: Vmax: flux decay plot outliers bug fix 
* Terraform configuration for CVM 
* Function View: Viewer positioning tags 
* Function View: Ability to show grid 


# 2021-12-09 Dev build 0.99.0

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.99.0`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* GitHub Actions: Update js-api/package-lock.json 
* Dev-tools: version bump 
* GitHub Actions: Update packages/DevTools/package-lock.json 
* 228 Chem: descriptors, incho, inchi keys, mcs refactoring and polishing, bug fixes 
* Compute: Vmax: editable input tab and transposed output tables 
* Compute: removed margin after info block in SelOut 
* Compute: centered the buttons in OutSel 
* Compute: removed padding on OutSel viewer 
* Compute: added padding to FE preview 
* Compute: added styling to FE preview 
* Compute: fixed export after update 
* Compute: Vmax update 
* Compute: patch version bump 
* GitHub Actions: Update packages/Compute/package-lock.json 
* Grok Connect: rename ClickHouse provider 
* Update public 
* Scatter plot: regression line: improve precision for the regression coefficient  
* Scatter plot: improve axis formatting (".00" postfixes are often unnecessary) 
* Scatter Plot: Rename and move property for MinMax labels visibility 
* Minor code cleanup 
* (Bug) Trellis Plot: full screen mode for inner viewers doesn't work when trellis is zoomed 
* Peptides #194: Add correlating positions highlighting. 
* Peptides #194: Release correlation analysis. 
* Peptides #93: Small adjustments. 
* Peptides #93: Release. 
* Peptides #194: Remove debug logging. 
* Peptides #194: Version up. 
* RepertoireBrowser: minor fixes 
* Compute: improved typing on function params 
* Compute: fixed function propss labels in FE 
* Compute: code preview generation based on prop type 
* Renamed .d to .dart to fix minified code 
* CSS fix 
* Grok Connect: maven cache in docker image 
* Compute: bumped minor version and JS API version 
* Viewers: Unification of properties (alpha/opacity/opaque) 
* Compute: Vmax and Pmax tables rearrangements 
* SQLite #47: initial commit (WIP) 
* (Bug) JS API: null reference in JsViewer.dataframe 


# 2021-12-07 Dev build 0.98.0

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.98.0`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* GitHub Actions: Update js-api/package-lock.json 
* Function View: Adjust table heights 
* Grid: add changePropertyPanel property   
* Fixed analyzer warnings 


# 2021-12-07 Dev build 0.97.0

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.97.0`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* GitHub Actions: Update js-api/package-lock.json 
* Chem: Set defaults in chem_common.ts 
* Scatter Plot: Marker coding for Split-By-Marker 
* Infrastructure Documentation (WIP)
* Peptides: removed unused code 
* GitHub Actions: Update packages/Peptides/package-lock.json 
* Charts: Package documentation (WIP)
* JS API: Improved code style (marker coding, formula lines) 
* (Bug) Dev panel doesn't work for any entities 
* Compute: export timeout increased 
* (Bug) Scatter Plot: Y axis shows weird values 
* Compute: fixed bug on script call from Model Hub 
* Chem: added Identifiers panel 
* Small help addition 
* Compute: export function works with single tab 
* Scatter Plot: Get rid of code duplication (axes and grid lines) 
* (Bug) Bar chart: When "Axis Type" = "logarithmic", numbers on the axis get merged 
* Chem: R Group Analysis fix 
* NPM repositories: new package manager UI (WIP)
* hem: MCS panel 
* Compute: Pmax: add polynomial regression lines (WIP)
* JS API: add Property options 
* Chem: panels Inchi, Inchi Keys 
* Compute: added outliers selection viewer 
* Compute: marked the outliers selection dialog as deprecated 
* Compute: minor version bump 
* GitHub Actions: Update packages/Compute/package-lock.json 
* Functions: help-url support for scripts in Function View (WIP)
* Compute: outliers selection additional column fixed 
* Ability to disable events in InputBase 
* Compute: improved the code generation 
* (Bug) Calculated column with numeric type and empty formula leaves no formula tag 
* Ability to show multiple viewers in Function View 
* Function View: Show table inputs if there are viewers set 
* small css adjustment 


# 2021-12-02 Dev build 0.96.0

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.96.0`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* GitHub Actions: Update js-api/package-lock.json 
* CI: publish libraries packages to npm 
* (Bug) Sequence Translator: translation from GCRS to MM12 is not displayed 
* Chem: a curious case where tests passes, but the code shouldn't work 
* SequenceTranslator: fix of SMILES concatenation of GCRS sequence 
* Chem: Setup routine (npm unlink/link) 
* Chem: Version bump 0.36.0 (no locks) 
* SequenceTranslator: fix of ps linkage concatenation 
* Fix concurrent modification during iteration exception 
* #181 Chem: Workers. wip: Passing ArrayBuffers to the main thread 
* Chem: Fixing external package-lock.json by updates 
* Tutorials: package.json update 
* Script editor: update 
* GitHub Actions: Update packages/Tutorials/package-lock.json 
* CSS Fix 
* Chem: Refactoring, fingerprint-centric methods 
* Chem: Refactoring, small renames 
* Peptides: some fixes and improvemets 
* (Bug) Formula Lines: Lines are clipped on some data 
* NPM repositories: new package manager UI (WIP)
* Packages: Add tests for versions 
* Chem: RdKit module external availability 
* JS API: Script.sample, Script.reference, Script.tags, Script.environment 
* Peptides: viewers are docked down 
* Peptides: Add simple correlation analysis. 
* add tutorial previews 
* Chem: Proper return of fingerprints with small datasets 
* Chem: Incorrect cache invalidation on re-use; common fixed setting for Morgan Fingerprints 
* Chem: Compute Morgan fingerprints on small datasets on a main thread 
* Chem: Remove workers' getSimilarities in favor of the main thread 
* Chem: Handling a case in caching when an input column is modified (on PR #39) 
* Chem: Redundant async calls in rdkit_service.ts 
* Chem: Remove workers' structuralAlerts in favor of the main thread 
* Documentation: Correlation plot can also show Spearman's correlation. 
* corrPlot: Add options description. 
* Chem: similarity analysis and SPE removal 
* Chem: Version bump 0.38.0 (no locks) 
* (Bug) Packages: Infinite detectors loop 
* Docker versioning 
* Chem: Moved service_worker part with dictionary to TypeScript 
* JS Api: added options parameter to chem.svgMol 
* Oligo Batch Calculator: rearrange object with codes 
* Fixed class names 
* PubChem: added PubChem panel 
* GitHub Actions: Update packages/PubChemApi/package-lock.json 
* Peptides: updated build script 
* Sequence Translator: rearrange object with codes 
* Chem: panels removal 
* Commit version to the according branch 
* GitHub Actions: Publish js-api from any branch 
* Fix condition to create branches in public repository 
* Git information inside docker containers 
* (Bug) In-vis Filter: Take into account property ShowFilteredOutPoints 
* Ability to show input parameters in Function View 
* Lib statistics: updated package.json 
* GitHub Actions: Update libraries/statistics/package-lock.json 
* Peptides: updated package.json 
* (Bug) MariaDB and MySQL databases throw an exception when trying to open a list of tables 
* Peptides: chem palette improvement 
* Fixed build errors 
* GitHub Actions: Update packages/Peptides/package-lock.json 
* New build script 
* Production js build setting 
* Separate source maps 
* Scatter Plot: Make scrollbars like a histogram and place them on the axes 
* Peptides #194: Add a box plot. 
* Peptides #194: Add missing changes. 
* More error handling for getNpmPackageVersions 
* Chem #197: descriptors 
* Peptides #194: Switch to Kendall's correlation. 
* Lib statistics: added increment after changes 
* GitHub Actions: Update libraries/utils/package-lock.json 
* Check js-api build before docker image build 
* Grok Connect: Upgrade shade plugin with relocated http 
* JS API: add Property options (WIP)
* Chem #197: descriptors polishang 
* Push version: get latest datagrok-api version from npm registry 
* Peptides: AAR grouping 
* Chem: polishing 
* Compute: fixed export on multiple plots 
* Packages: optimised search for compatible versions 
* Compute: path version bump 
* Compute: Vmax: Pmax: WIP 
* Compute: js-api version bump and patch version bump 
* Compute: fixed buttons height 
* Compute: patch version bump 
* GitHub Actions: Update packages/Compute/package-lock.json 
* Fixed addNewColumn layout 
* Chem: returning locks to the package functions 
* Unit tests: exclude and include components from command line 
* Compute: Vmax: Pmax: output parameters joined in tables 
* (Bug) Grok connect: NPE after dev deploy  
* Updated dev environment v2 for Linux users 
* Clinical case: added warning in case demo fils or SDTM data not loaded 
* Compute: fixed grids export 
* Peptides #194, #93: Fix Peptide Space. 
* SequenceTranslator: fix of SMILES concatenation 
* GitHub Actions: Update packages/SequenceTranslator/package-lock.json 
* (Bug) Filters: Uncheck/Check a field in the filter panel that does not update the record properly 
* Libs and Peptides versions bump 
* Compute: added graphics export feature 
* Packages: allow subfolders in the directories reserved for entities 
* (Bug) Grok Connect: Athena: java.lang.NullPointerException at java.util.Hashtable.put(Hashtable.java:460) 
* Datlas: remove timer for reconnecting to Grok connect 
* Peptides: dataset restoration 
* Grok connect: improve SettingsManager for tests and GrokConnectShell 
* (Bug) Formula Lines: Sometimes tooltips doesn't appear 
* Wiki: Access \- WIP 
* Charts: Package documentation (WIP)
* Peptides #161, #193: improvements 
* Peptides: version bump 
* Compute: added Function Editor 
* SequenceTranslator: version bump 
* Package files caching 
* RepertoireBrowser:  updates, new PTM obs tracks 
* GitHub Actions: Update packages/RepertoireBrowser/package-lock.json 


# 2021-11-25 Dev build 0.95.9

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.95.9`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* GitHub Actions: Update js-api/package-lock.json 
* Jest: Moving things around for Chem and others 
* NPM repositories: new package manager UI (WIP)
* #181 Chem: Workers. wip: Access locks on grok API chem functions 
* Peptides: Add dimensionality reducer worker. 
* Peptides: Remove console.log. 
* Formula Lines: Hide '...' button from Filter field in PP 
* Change routing for CVM components in AWS 
* Peptides: Add peptide space as a widget. 
* Chem: Descriptors app 
* (Bug) Packages: dapi.packages.find ignores the include string 
* (Bug) Packages: URL field is empty for packages deployed from npm repositories 
* (Bug) Sequence Translator: molecule SVG and MOL file don't match (WIP)
* README 
* Chem: RdKit module external availability, refactoring 
* Chem: R-Group Analysis fix 
* add ClickHouse connector 
* Datlas: Add logging to ConnectorsPlugin 
* Datals: add tryToReconnectToConnectorServer flag 
* Chem: version bump, for further evaluation 
* In-vis Filter: Implement for Line Chart 
* #181 Chem: Workers. wip: Simplify Similarity to an array of BitArray-s 
* #181 Chem: Workers. wip: Prepared Substructure Search for pattern fingerprints 
* Peptides: Add peptide space in property panel. 
* Packages: version bump 
* GitHub Actions: Update packages/ChemblBrowser/package-lock.json 
* GitHub Actions: Update packages/BioSignals/package-lock.json 
* GitHub Actions: Update packages/Bio/package-lock.json 
* GitHub Actions: Update packages/CustomML/package-lock.json 
* GitHub Actions: Update packages/DSP/package-lock.json 
* GitHub Actions: Update packages/ChaRPy/package-lock.json 
* GitHub Actions: Update packages/ClinicalCase/package-lock.json 
* GitHub Actions: Update packages/Discovery/package-lock.json 
* GitHub Actions: Update packages/Charts/package-lock.json 
* GitHub Actions: Update packages/DevTools/package-lock.json 
* GitHub Actions: Update packages/DrugBank/package-lock.json 
* GitHub Actions: Update packages/PhyloTreeViewer/package-lock.json 
* GitHub Actions: Update packages/NLP/package-lock.json 
* GitHub Actions: Update packages/MultiPlot/package-lock.json 
* GitHub Actions: Update packages/ScatterPlot3D/package-lock.json 
* GitHub Actions: Update packages/TensorFlow.js/package-lock.json 
* GitHub Actions: Update packages/Tutorials/package-lock.json 
* GitHub Actions: Update packages/Ketcher/package-lock.json 
* GitHub Actions: Update packages/Impute/package-lock.json 
* GitHub Actions: Update packages/Notebooks/package-lock.json 
* GitHub Actions: Update packages/VDJtools/package-lock.json 
* GitHub Actions: Update packages/Viewers/package-lock.json 
* GitHub Actions: Update packages/UsageAnalysis/package-lock.json 
* Chem: Preparing a transport for fingerprints 
* Tutorials: fixed ribbon hint indicator position 
* DevTools, Script-editor: Fix ribbon panels 
* Chem: Version bump 0.31.0 (no locks) 
* Chem: Version bump 0.32.0 (no locks) 
* null 
* Fixed build errors 
* Chem: Version bump 0.33.0 (no locks) 


# 2021-11-23 Dev build 0.95.8

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.95.8`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* GitHub Actions: Update js-api/package-lock.json 
* (Bug) Socket fails if there is ID in DataFrame Param value instead of TableInfo (WIP)


# 2021-11-23 Dev build 0.95.7

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.95.7`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* GitHub Actions: Update js-api/package-lock.json 
* Packages: get compatible npm package versions 
* NPM repositories: new package manager UI (WIP)
* Moved types to deps 
* Fixed a bug with p-values sort. 
* Update compute.md 
* Model Catalog improvements 
* Compute: scripts cleanup 
* Clinical case: fixed width of choice inputs in Time profile, added ability to select values for ALT, AST, BLN when creating Hy's law plot 
* Cimpute: fix popPK 
* Compute: IntervalsFromECG script 
* fixed library 
* Scatter Plot: Hide property "formula-lines" from PP 
* Chem: similarity search port 
* Chem: fixes 
* Peptides: package restructuring 
* Chem: 3d coordinates script optimization 
* Compute: added additional column of flux 
* Remove a utility function 
* (Bug) datagrok-tools: `grok add` inserts invalid import statements 
* Utils: \* Update operations to in-house implementation. \* Add JSDoc. 
* Jest: Moving things around for Chem and others 
* #53 Tutorials package: markup adjustments to the tutorial page 
* Added trick with .bind for callback distance. 
* Fixed build issues 
* Peptides: SAR viewers are split apart 
* Added missing dependency 
* Clinical case: fixed bugs with missing domains validation for ae browser view 
* Compute: fixed column name 
* Compute: moved methods under eslint-ignore line 
* Compute: fixed the column computation 
* Compute: patch version bump 
* Fix in assert import. 
* Added docstrings missing. Implements #93 
* Clinical case: updated Readme 
* Clinical case: moved Laboratory view parameters to property panel 
* Compute: changed reason to rationale term 
* Tutorials: Add sticky header 
* Tutorials: add chevron icon for step 
* Tutorials: update css 
* Fixed css 
* (Bug) Multiple FuncCall.onCancel subscriptions (WIP)
* Additional debug lines 


# 2021-11-22 Dev build 0.95.6

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.95.6`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* GitHub Actions: Update js-api/package-lock.json 
* Model Catalog improvements 
* new version 
* Added documentation on caching 
* Added metadata to demo models 
* Compute: minor UI improvements 
* Compute: help: embedding as iframe 
* OnViewAdding event 
* Added a gif on how to share a connection 
* Compute: documentation 


# 2021-11-21 Dev build 0.95.5

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.95.5`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* GitHub Actions: Update js-api/package-lock.json 
* Compute: documentation \- WIP 
* Model Catalog improvements 


# 2021-11-21 Dev build 0.95.4

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.95.4`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* GitHub Actions: Update js-api/package-lock.json 
* JKG: lost RxODE 
* Fixed analyzer warning 
* Compute: documentation \- moved to the core help 
* Moved help on scripting help to /compute 
* Moved help on jupyter help to /compute 
* iframe-embeddable views 


# 2021-11-21 Dev build 0.95.3

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.95.3`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* GitHub Actions: Update js-api/package-lock.json 
* Compute: documentation \- WIP 
* (Bug) Filters: Uncheck/Check a field in the filter panel that does not update the record properly 
* Viewers: In-vis filtering: Implement for Bar Chart 


# 2021-11-21 Dev build 0.95.2

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.95.2`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* GitHub Actions: Update js-api/package-lock.json 
* Scatter Plot: Scroll bars visibility 
* RepertoireBrowsser: dynamical viewer logo 
* Filters: add wildcard support in search over categories 
* Formula lines: Update help 
* Clinical case: updated missing domains and columns validation 
* Clinical case: added vital signs domain to distributions, correlation matrix and time profile views 
* Peptides: ngl in property panel 
* Peptides: code cleanup and minor improvements 
* Work in progress 
* Minor code cleanup 
* Statistics lib: removed padjust 
* (Bug) In-Vis Filter: Does not affect the regression line 
* JS: MultiView: DockView support 
* Close button on TabControl 
* Convert FunctionView to DockView 
* JS: MultiView: DockView support 
* Disabled MultiView 
* Samples scripts: naming harmonization 
* Model Catalog improvements 


# 2021-11-19 Dev build 0.95.1

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.95.1`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* REST API for functions 
* Core: Progress indicator for server functions 
* Connections: add connection tests (WIP)
* JS API: ui.input.forProperty 
* JS API: 'slider' editor type (WIP)
* JS API: ui.tableFromProperties 
* #180: Compute: UX harmonization \- WIP 
* (Bug) Trellis Plot: full screen mode for inner viewers doesn't work when trellis is zoomed (WIP)
* Viewers: Box Plot: In-vis filtering 
* Wiki improvements 
* Lines by equations: Extending lines so they can be drawn on different visualization independently 
* Update help & js examples 
* Datlas: Add logging to ConnectorsPlugin 
* Settings documentation improvement 
* Viewers: harmonize property names for shapes 
* Fixed misspelling 
* Packages: get compatible npm package versions (WIP)
* (Bug) connection.close in JdbcDataProvider.execute didn't account for exceptions other than SQLException 
* Additional tests for calculated columns 
* Push version script: releases branch 
* GitHub Actions: Update js-api/package-lock.json 
* Scatter Plot: Scroll bars visibility 
* Docker build script: releases branch 
* \* Add molecular weight calculation for a peptide sequence. 
* Lines by equations: Ability to create dashed and dotted lines 
* Peptides: various fixes and improvements 
* Lines by equations: Ability to remove lines from storage 
* Datlas: keep trying to connect to Grok connect 
* Removed relative paths to datagrok-api 
* MultiView css fix 
* GitHub Actions: Update packages/Ketcher/package-lock.json 
* GitHub Actions: Update packages/ScatterPlot3D/package-lock.json 
* Compute: documentation \- WIP 
* Fixed datagrok-api dependency following the concept of use of 'npm link'. 
* Fixes following linter hints. 
* Add correction for multiple hypotheses testing. 
* \* Fix imports. \* Add required packages. 
* Fixed import and warning 
* #53 Tutorials package: markup adjustments to the tutorial page 


# 2021-11-18 Stable version 0.95.0

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.95.0`
  *  `docker pull datagrok/datagrok:stable`
* CVM: 
  *  `docker pull datagrok/cvm:0.95.0`
  *  `docker pull datagrok/cvm:stable`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* updated clinical demo files project 
* Script Editor, snippets menu 
* #53 Tutorials package: EDA track (WIP) 
* Selenium: Add attributes "name" to elements of the new AddNewColumn dialog 
* JS API: onDialogShown event 
* Rename d4-input-editor class 
* Add superclass to the Dialog widget 
* (Bug) Grok connect: non-parameterized query is not logged 
* Grok connect: log params in parameterized query 
* Update README.md 
* Wiki: Peptides 
* Update script-editor.ts 
* Fiexed misinformation in data access tutorial 
* JS API: add dialog's title 
* Connectors: Support Array input parameters #90 
* JS API: Dialog.inputs 
* Can't execute .applyFormula on column without formula tag #95 
* chem:svgMol function -> molfile format support added 
* DrugBank package -> search functions and widgets added 
* Data on demand 
* #53 Tutorials package: Data connectors tutorial 
* DrugBank package -> excessive header removed 
* Added cache checkbox 
* resolves #97 datagrok-tools: package content validation: check scripts location 
* Fix UI tests (WIP)
* Fixed js-api links 
* sdf file-reader 
* Fix a typo 
* #53 Tutorials package: lock `rxjs` version 
* sdf reader fix 
* Update package.js 
* #53 Tutorials package: Multivariate analysis tutorial 
* added survival analisys view 
* #53 Tutorials package: temporary ts-ignore for new functionality 
* Wiki: misc small fixes. 
* ApiSamples: a case typo in a sketcher sample. 
* Wiki: initial Style Guide for Help contributors. 
* Release: adding a versioned hand-crafted Release Notes. 
* updated survival analysis view and r scripts 
* moved reading of validation rules table to clinicalCaseApp() function 
* Release Notes: 2021-05-10, accumulating. wip 
* Release Notes: 2021-05-21, accumulating. wip 
* Release Notes: 2021-06-10, accumulating. wip 
* Update render-items.js 
* Update custom-cell-rendering-indexes.js 
* Release Notes: 2021-07-14, accumulating. wip 
* Release Notes: 2021-07-29, accumulating. wip 
* Release Notes: 2021-07-29, short annotation. wip 
* Release Notes: 2021-07-29 Build 0.93.0 
* #53 Tutorials package: adjustments to new UI, code clean-up 
* Release Notes: cleanup. 
* Wiki: Update filter help 
* corrections of ui for survival view 
* Peptides: generate a dataset for further analysis 
* (Bug) Grok connect: query isn't logged on exception 
* Filters wiki page update 
* Wiki: Dataframe (WIP)
* updated survival analysis view 
* Initial commit 
* Fixed training 
* Fixed typos 
* #104 Wiki: Library tour 
* Add ui.switch to examples and docs 
* added creation of survival dataset (for first SAE occurance), updated survival analisys view 
* UI improvements in survival analysis view 
* ApiSamples: fixed a typo. 
* Docker Compose: a new CVM yaml file. 
* Docker Compose: adding an old yaml for history. 
* Docker Compose: removing the old yaml. 
* Docker Compose: linkified the Docker yaml. 
* Wiki: Custom machine learning models. Cleanup 
* Wiki: Custom machine learning models. Links 
* Charts: fix a package publication error 
* JS API: update MapChangeArgs interface 
* #53 Tutorials package: change node placement 
* #53 Tutorials package: EDA: viewers (WIP) 
* JS API: add event args types to onViewerAdded/onViewerClosed 
* Oligo Batch Calculator: detect errors when code table contains one-char and two-chars codes 
* Sequence Translator: Creating a mol file from the landing page 
* Pairwise sequences alignment 
* JS API: Create percentile calculation function 
* Clarified what the percentile function is. 
* added hospitalization and drug related ae to Kaplan-Meier endpoints, added ability to choose lab values on Patient profile 
* fized tooltips, yAxesLabels and linechart title in multiplot 
* Update package.ts 
* Stacked barchart added 
* Wiki: Capitalization 
* Wiki: Capitalization. Adjustments for source codes with # 
* Decompose the calculated columns example 
* API Samples: a few additional links 
* #110: Chem: dynamic sketchers 
* GROK-GitHub-109 JS API: DataFrame.insert wrong order of optional argument 
* GROK-GitHub-109 JS API: DataFrame.insert wrong order of optional argument. Type fix 
* API Samples: simplify an info panel example 
* Wiki: Small rendering adjustments 
* JS API: grok.shell.windows.showRibbon flag 
* JS API: dataframe.ts cleanup. 
* null 
* aligned sequence split function 
* fixed splitAlignedPeptides function 
* added splitAlignedSequence panel function 
* closes #113 datagrok-tools: package template: add eslint config file 
* Update CONTRIB.md 
* Peptides: compile the initial dataset, detect sequence, split   
* rewrite Peptides package in TypeScript 
* Add SAR viewer 
* JS API: fix type warnings 
* Peptides: refactor 
* Peptides: simplified the aligned sequence detector 
* Stacked barchart: Selection style changed, datagrok palette added, axis display simplified, unreadable bar labeles removed, trying to fix escaping bars bug. 
* #114: Bio: grid cell renderer for sequences 
* Sequence: made it work again; code cleanup. 
* Wiki: Restructure. Access (WIP)
* escaping bar bug is fixed 
* JS API: BitSet.and, or, xor, andNot 
* Peptides: drawing circles in a grid instead of numbers 
* Peptides: viewer rework 
* (Bug) JS API: Table.rows.filter() doesn't work if called after opening the table 
* Minor improvements 
* JS API: a type for property options 
* Peptides: added activity scaling options 
* Peptides: source table selection based on viewer select 
* Peptides: refactor according to coding guidlines 
* JS API: Minor harmonization 
* Peptides: dataset fix 
* added ability to select values for y axe in multiplot, updated patient profile and timelines views 
* Peptide Cell renderer 
* #53 Tutorials package: ML (predictive modeling tutorial): missing values imputation 
* Clinical Cases UI update 
* (Bug) Apps: App button goes to common heap instead of separate tab 
* Rollback 
* Peptides: added help to property panel, minor improvements 
* Fixed views creation 
* Fixed typing 
* Peptide cell renderer 
* Peptides: color scheme for sequences 
* update clinical cases UI 
* Packages: change the API source 
* DevTools: a polyfill for replaceAll 
* Peptides: scaling method replace ln -> lg 
* Fixed Analyzer warning 
* Removed old dependency 
* lock files 
* Barchart: Font changed to serif-sans, column numeration changed, current row highlighting added 
* Peptides: renderer fix 
* Peptides: added histogram and updated stats, scaling fix 
* (Bug) Chem: External sketcher, filtering and SMARTS 
* JS API: add missing toJs in package.getProperties() 
* Wiki: package settings 
* exposed JSMol::get_smarts 
* Peptides: AA sem type detector 
* Peptides: property panel WIP 
* Add "iconTool" to devTools 
* Peptides: property panel fix 
* Peptides: finished converting to TS, code cleanup 
* Linter rules update 
* Removed an unused dependency 
* UsageWidget: code cleanup 
* #117: PowerPack: PowerSearch: search templates \- added server-based collections of templates 
* implemented ability to set min and max values in multiplot 
* #117: PowerPack: PowerSearch: build fix 
* #117: PowerPack: PowerSearch: search templates \- added package settings for templates paths 
* implmented multiple lines on linechart using eChart graphs 
* fixed multiplot tooltips bug 
* implemented changing height of scatter plot with categories multiselect in multiplot 
* #53 Tutorials package: ML (predictive modeling tutorial): actions for the model training view 
* updating survival plots on tab click in case filters changed 
* Peptides: Aligned sequence cell renderer text properly aligned, readability increased 
* Peptides: analyzePeptides fix 
* Peptides: color coding circles based on MAD 
* Peptides: help fix 
* DevTools: add data connection examples 
* #53 Tutorials package: code clean-up 
* #53 Tutorials package: helpers for working with view inputs 
* Update icons for packages 
* Update package.png 
* Update Laboratory view 
* Peptides: color coding optimization 
* #53 Tutorials package: openViewByType method 
* JS API: add examples with the linear color-coding 
* Peptides: color-coding fix, split sequence join 
* Closes #118 JS API: color-coding methods (Column.colors) 
* Peptides: Stacked barchart transfered to ts and added to package 
* logo Viewer for Peptides 
* Peptides: Fixed cell renderer column width 
* Update about-widget.ts 
* Fix property initialization warnings 
* DG.Package: methods for working with files 
* SMARTS: test possible solutions for aromatization/kekulization 
* Tutorials package: give rights for editing 
* JS API: DG.Color.linear() 
* Peptides: ui fixes 
* JS API: menu order parameter 
* JBIO: prepare color schemes for NGL viewer 
* Peptides: various fixes and improvements 
* Peptides: refactoring 
* datagrok-tools: adjust lint scripts for recursive directory check 
* changed enrollent linechart to cumulative sum, implemented request to clinicaltrials.gov for study info, added split of AEs barcharts by treatment arm 
* JS API: Add Legend widget 
* tutorials update 
* Peptides: typo fixed 
* Peptides: Added AA cell-renderer, changed sequence cell renderer 
* update icon 
* #53 Tutorials package: change the storage structure (refactoring, WIP) 
* #53 Tutorials package: a proper way to access a tutorial's track (refactoring, WIP) 
* Closes #119 Chem: Alignment differs for same MolBlocks with different offsets 
* #53 Tutorials package: css files (refactoring, WIP) 
* Closes #120 Chem: Rotation flicker when using alignment 
* Peptides: AA render in sar viewer 
* #117: PowerPack: PowerSearch: search templates 
* update tutorials. Add css style 
* tutorials images 
* update tutorials 
* Closes #121 (+ #77). Chem: Support SMARTS 
* Peptides: histogram rework WIP 
* #53 Tutorials package: code clean-up, hint styling 
* Peptides: histogram fix WIP 
* update imgages for tutorials 
* #122 JS API: RowGroup class 
* #53 Tutorials package: openDialog method 
* Peptides: histogram rework and fix 
* Closes #123 Chem: Substructure search to intercept SMARTS if MolBlock fails 
* tutorials updates 
* #53 Tutorials package: cheminformatics track: descriptors 
* JS API: onContextMenuItemClick 
* fixes #119 
* Chem: property-panel WIP 
* #53 Tutorials package: ability to pass a cusom promise for an action, helper for selecting an item from the context menu, code clean-up 
* Chem: Substructure filter 
* Docker Compose: Remove profiles 
* (Bug) Filter resets when user recalculates column formula 
* Peptides: t-test, u-test and histogram changes WIP 
* Compute package, initial commit. 
* Compute: Readme.md, continued. 
* JS API: add more specific types 
* update tu 
* Revert "Merge branch 'master' of https://github.com/datagrok-ai/public" 
* Pubchem: Search panels 
* implemented legend and switching type of linechart in multiplot. Implemented study summary property panel 
* HitTriage: initial commit 
* HitTriage: documentation: work in progress 
* Grok Connect: Update script 
* Peptides: center aar 
* Peptides: various improvements WIP 
* Peptides: polishing 
* Peptides: add a real-like dataset for analysis 
* updated tooltips in timelines view, added dynamic change of lines width depending on zoom to timelines chart 
* Peptides: Cell renderer split sequence fix & Header barchart 
* #53 Tutorials package: tutorial cancelation (WIP) 
* Peptides: verious fixes and improvements 
* #53 Tutorials package: set initial value in `clearRoot` 
* #53 Tutorials package: remove unnecessary status update 
* Substructure filter:undefined column & searchSubstructure changes fix 
* JBIO: prepare tooltips and update ux 
* Peptides: fix header renderer 
* Peptides: clipped all cells 
* Peptides: minor fixes and improvements 
* Peptides: BC viewer is now rendered in onCellRender 
* Peptides: Snake case function name fixed 
* Peptides: Header BC selection and current column highlighting & labels added 
* Peptides: more minor fixes and improvements 
* #53 Tutorials package: make step description optional 
* #53 Tutorials package: show completed step descriptions on click 
* Peptides:Cell renderers for aa and aligned sequences for advanced dataset & some colors for barchart 
* Peptides:Fixed 
* Peptides:Fixed Errors 
* Peptides: moved logo viewer 
* Peptides: u-test rework 
* Peptides: small fixes and improvements 
* Form: onFormCreating event for adjusting the default column selection 
* Commented out BitSet.setFast() 
* #117: PowerPack: PowerSearch: code cleanup and refactoring 
* Added a picture for dynamic molecular sketchers 
* JS API: "file edited" event (WIP)
* Resolves #134 
* DataFrame.onSorted event 
* DataFrame.onRowsFiltered event 
* #53 Tutorials package: Predicitve Modeling (WIP) 
* #53 Tutorials package: Multivariate Analysis 
* Peptides: Added indent & small fixes 
* JS API: Functions: Added help for toString() function 
* Peptides: Long modifications hidden & tooltip added 
* #53 Tutorials package: allow text selection in steps 
* #53 Tutorials package: track restructuring 
* #53 Tutorials package: track progress formatting 
* property panel in Adverse events view, updated box plot view(added default box plots based on p-value, added strata choice) 
* clinical case fix box-plots view 
* JS API: TabPane.header 
* added check for number of stratas in box plots view 
* JS API: Ability to access FuncCall options and aux 
* Strict always on 
* Better events in AddNewColumn 
* #53 Tutorials package: add missing help URL 
* update tutorial widget 
* #53 Tutorials package: data connectors tutorial 
* Fixed links 
* #53 Tutorials package: handle edge cases in openViewByType 
* Peptides:Basic responsive header barchart 
* Peptides:Fixed color coding rules for AAR && AA are now rendered with short names 
* update tutorial images 
* #53 Tutorials package: Predicitve Modeling 
* Removed debug printout. 
* update widgets 
* #53 Tutorials package: viewer basics (WIP) 
* Peptides:Fixed cell renderer 
* Peptides: analyse and select the most appropriate subset in the "real@ dataset for tool performance 
* Peptides: fixed L color and added full AA names 
* JBIO: removing test files 
* #53 Tutorials package: viewer basics 
* Filters: proper detaching from events 
* (Bug) Chem Filter: After reopening the filter, browser tab with Datagrok freezes (WIP)
* #53 Tutorials package: Filters (EDA track) 
* Peptides: sketcher-tooltip for aar & refactoring chem palette 
* Peptides: Barchart fixes & tooltip added 
* Peptides: Sar aar tooltip added 
* Peptides: Aligned sequens are now aligned 
* HitTriage WIP 
* Peptides: Full sequence mol-graph widget 
* Peptides: polishing 1 for user meeting 
* Peptides: Aligned sequence are aligned for sure 
* Better typing 
* Fixed CSS issue 
* (Bug) JS: Property setter doesn't work on FuncCall 
* Package files: Remove path prefix 
* update learning widget 
* Peptides: Initial table columns are hidden& Full AAname added to tooltip 
* Revert "Better typing" 
* Peptides: Activity column unhidden 
* Peptides: Barchart fix 
* Peptides: help 
* Peptides: fix the number of significant digits in statistics 
* Peptides: minor polishing 
* (Bug) Add New Column: Exception after opening dialog for editing an already added column 
* #135: DevTools: a lightweight package unit-testing framework 
* PowerPack: Ability to hide welcomeView 
* Added "molecules" demo dataset that could be used synchronously (grok.data.demo.molecules) 
* (Bug) Chem: Pass SMARTS properly to searchSubstructure (WIP)
* Widget.temp property (WIP)
* added time profile and adverse event browser views 
* Programming exercises (WIP)
* #53 Tutorials package: Viewers (EDA track) 
* Sequence Translator: enhancements to AxoLabs display  
* Update tutorials widget 
* Peptides: various fixes and imporvements 
* JS API: Add method to easily customize lines & bands in a Scatter Plot 
* SequenceTranslator: OligoBatchCalculator: added image and links to YouTube presentations 
* #53 Tutorials package: chem: Activity prediction (WIP) 
* Help & JS Examples: Add info about lines & band in a Scatter Plot 
* Wiki: Simplify Develop and Getting Started (WIP)
* updated AE browser property panel, fixed bug with property getter in Timelines viewer, fixed bug in Survival analysis view 
* removed dm from additional domains in AE browser property panel 
* Fixes #136: Client-side SDF file reader 
* Removed the sdf file reader from Chem (since it is now part of the OpenChemLib) 
* Peptides: minor improvements 
* Help: Change youtube video to screenshots (WIP)
* #53 Tutorials package: chem: update datasets 
* #132 Chem: Move Substructure Filter panel to the package 
* Removed debug printout 
* #132 Chem: Move Substructure Filter panel to the package. Grooming 
* #138 Chem: Text search for molecules in the sketcher. wip: SMILES only 
* Update CONTRIB.md, grooming 
* Updated UI in Time profile, Box plots, Matrix, Summary views. Fixed bugs in Survival analysis view 
* #132 Chem: Move Substructure Filter panel to the package. Collaborative filtering 
* #132 Chem: Move Substructure Filter panel to the package. Substruct library 
* #139 Chem: Move from .tags to .temp 
* Expose DataSourceCardView to TS-API (WIP)
* ModelHub initial commit 
* Peptides: various fixes and improvements WIP 
* added summary view description to readme 
* summary view screenshots 
* updated images links 
* added gif of summary view 
* added summary gif 
* set gif width and height 
* timelines view description 
* corrected typo 
* patient profile view description 
* correced markdown 
* All dependent columns should be recalculated after new formula was applied (WIP)
* Molecular Liability Browser: refactor to handle common actions in NGL and pViz 
* adverse events view description 
* Set version to 0.94.0 
* Grok connect: log query execution time by steps  
* Seldon: initial package commit 
* added laboratory view description 
* #93 Peptides: fixed 2d view 
* Peptides: Barchart fix again && sar-list-viewer && little fixes 
* Peptides: filter fix and other improvements 
* Peptides: Sequence width 
* biomarkers distribution description, updated biomarkers distribution view(created ribbon panel) 
* #53 Tutorials package: base class improvements 
* added correlation matix view description 
* added time profile description 
* Peptides: viewer row sorting and colors fix 
* added ae browser description 
* Peptides: Sar-list update 
* Revert "Peptides: Sar-list update" 
* Revert "Merge remote-tracking branch 'origin/master'" 
* Update Usage analysis package 
* Peptides: deleted commented code 
* updated boxplots view 
* Update icon 
* JS-API: ContextActions method 
* JS-API: star() method 
* JS-API: Accordion pane counts 
* JS-API: author getter 
* JS-API: Accordion.addTitle method 
* Renamed package 
* add chemSimilarityAnalysis 
* added validation view description, fixed bug with time profile view 
* #53 Tutorials package: add missing titles 
* updated ae browser gif 
* Connection pool test 
* Clinical Case: update property panel 
* View: helpUrl setter 
* #53 Tutorials package: chem: renaming 
* added p-values to box plots, updated survival analysis view 
* Added example of setUpdateIndicator. 
* Peptides: minor refactoring 
* Model Catalog stubs 
* Support Ad-Hoc queries in JS 
* Update localhost.docker-compose.yaml 
* Fixed typos and replaced oudated resources in Help/develop and Help/exercises sections. Added illustrations. 
* Update package icons 
* JBIO: all PTMs tracks syncing 
* Wiki: add missing heading for chem articles 
* Sequence Translator: Enhancements to SVG 
* (Bug) Usage analysis: UA widget doesn't work for JnJ 
* JS API: Menu.root 
* #53 Tutorials package: getMenuItem helper for visual hints 
* #53 Tutorials package: entry function fix 
* #53 Tutorials package: add visual hints for the top menu 
* Locked postgres image version 
* JBIO: PTMs syncing with paratopes and cdr 
* JBIO: motif PTMs in own track and group 
* Peptides: major refactoring, viewers merging 
* Peptides: Refactor & nasty bugs fixed & selection mode for SARV 
* ChemblApi 
* #150 Chem: Move to webpack and TS. wip 
* added laboratory linechart with changes within normalized normal ranges, changed AE browser view to TableView, updated survival analysis view UI 
* added survival analysis description 
* Code cleanjup 
* JS-API: Packages test framework (WIP)
* Fixed naming 
* JS-API: AccordionPane.root property 
* JS API: Filters (WIP)
* JS API: Accordion.header 
* Pubchem API 
* Wiki: Using lists in parameterized queries 
* #150 Chem: Move to webpack and TS. Moved to webpack 
* add web workers 
* #53 Tutorials package: add visual hints for the sidebar panes 
* Accordion: more documentation 
* Minor code cleanup. 
* Removed accidentaly added functions 
* Resolves #165: Folder Browser Preview Functions 
* Formatting in Compute fixed 
* Added manual outliers detection 
* Peptides: removing lock 
* Sequence Translator: option to be able to  view just users’ pattern  
* Added eslint config and formatted all existing files 
* #53 Tutorials package: chem: virtual screening 
* SequenceTranslator: minor code simplification 
* Oligo Batch Calculator: assume sequence representation by it's first valid code 
* datagrok-tools: `grok api` command 
* corrected typos in readme, updated ae browser property panel, adde labels to reference lines in laboratory linechart 
* Added generated dataframes return 
* Removed puppeteer dep 
* implemented lazy view loading 
* PubChemApi: export PubChem class 
* ChemblApi: Iupac name function 
* Added promise resolve on dialog closing 
* Added protection from multiple dialog opening 
* Added column existence check before the column insertion 
* Added function abort on cancel event 
* Added RemoveOutlier button 
* Fixed backgroundColor option to be optional 
* Changed "detection" tem to "selection" term 
* Lasso tool is enabled by default 
* Fixed typo in ui.list example 
* added help for summary view 
* added help urls for views 
* moved seting of help urls to constructors 
* Lesser fixes for tutorials 
* use byteArray type as fingerprint in similarityAnalysis function 
* Fixed scatterPlot width 
* Changed non-mutating outlier selection on mutating 
* updated laboratory help 
* Timelines: remove a hard-coded column name 
* (Bug) Package Properties: default values are ignored 
* Replaced ScatterPlot legend by Grid 
* (Bug) Oligo Batch Calculator: grid doesn't appear if sequence was recognized as valid 
* Added group delete button to the group row 
* #150 Chem: Move to webpack and TS. Fix wrong refactoring 
* new release history 
* Molecular Liability Browser: intremidiate commit 
* Chembl: added few more queries 
* Added remove group button to group grid 
* add minor improvements 
* Grok connect: connection pool 
* #169 Chem: Molecule rendering to a canvas with an optional SMARTS scaffold 
* #150: Chem: Move to webpack and TS. Moving to TS, wip 
* #150: Chem: Move to webpack and TS. Moving to TS, wip. Compiles 
* #150: Chem: Move to webpack and TS. Moving to TS, wip. Workers work 
* #150: Chem: Move to webpack and TS. Moving to TS, wip. Caches work 
* add class BitArray 
* #53 Tutorials package: delete data available via the demo files connection 
* PubChemApi: make IUPAC name function async 
* Prototype of exportCallFunc added 
* Moved exportFuncCall to the distinct file 
* Refactored exportFuncCall 
* JS-API: property.category 
* Refactored to keep interface clean 
* Changed mark/unmark buttons behaviour 
* Closes #166: Chem: moved panels from core 
* Update parameterized-queries.md 
* Linked to the YouTube video 
* Fixed return type on propertyType getter/setter 
* Removed type casts on propertyType calls 
* Added automatic detection dialog 
* Export function refactoring 
* Updated Datagrok/Grok wording, fixed missing links 
* CI: publish packages to npm 
* improved performance of views loading, added check for missing tables 
* Added dynamic autodetection dialog render 
* Chanched functions search principle 
* fix some bug 
* Chem: chembl & pubchem id/link to report card added 
* Grok connect: added Hikari library 
* Chembl: fix 
* ChemblBrowser: fix 
* Added MVP of automatic outliers detection 
* Oligo Batch Calculator: get codes from input object 
* Closes #171 Chem: implemented OCL renderer 
* Fixed the case of multiple autodetections 
* Fixed imports in outliers selection 
* Change return type 
* datagrok-tools: code restructuring 
* Fix package names 
* GitHub Actions: Update packages/PhyloTreeViewer/package-lock.json 
* js-api right version [skip actions] 
* Grok connect: No suitable driver found for Hive query #173 
* Grok connect: connection pool safety test  
* (Bug) Grok connect: connection pool exceptions aren't shown 
* Release Notes: 2021-10-14 Build 0.94.0 
* Column.init: added boolean value initializer 
* #106 Chem: Add R-group analysis to the package. wip 
* GitHub actions: git push 
* GitHub actions: workflow_dispatch 
* tools package-lock.json 
* GitHub Actions: pull 
* GitHub Actions: Update tools/package-lock.json 
* GitHub Actions: upgrade npm 
* add unsigned bit movement 
* #53 Tutorials package: replace the scroll method 
* add default values 
* Added hints for selection usage 
* Changed the dialog layout to be more descriptive 
* GitHub Actions: Update js-api/package-lock.json 
* (Bug) Grok connect: correctly return connections to Hikari  
* datagrok-tools: add ignore files 
* updated validation for missing domains, updated timeines property panel 
* Added Grok connect log string 
* Non-interactive mode for outliers detector 
* Grok Connect: Fail in case of any error during grok_connect.sh 
* Grok connect: Added HikariCP to Maven dependency 
* Grok connect: removed test 
* SimPKPD: RxODE example 
* datagrok-tools: refactoring 
* datagrok-tools: run config test 
* JS API: expose GridCellDataArgs and corresponding events 
* JS API: Func parameters 
* DataFrame.fromProperties(properties, rowCount) constructor 
* View.fromRoot constructor 
* Grok Connect: Fix "Invalid signature file digest for Manifest" error Excluded manifest signature files in maven-shade-plugin 
* Docker Compose: valid YAML file 
* create a new file for SPE 
* Added tooltips to buttons 
* Added help section to the top of the dialog 
* Added context help docs 
* add code style fix 
* #53 Tutorials package: dashboards, scripting (WIP) 
* #53 Tutorials package: delete Molecular Descriptors (covered in Virtual Screening tutorial) 
* Change package naming 
* GitHub Actions: Update packages/BioSignals/package-lock.json 
* Added augmentation for the case with no dialog callled 
* (Bug) “Outliers…” button cannot be pressed after a single run 
* (Bug) Functions: After opening Function View, "Scripts" look different 
* Running calculators from “ModelCatalog” with double-click 
* (Bug) Functions: Function view doesn't reflect the changed script 
* add BitArray tests 
* Added inputs export to file 
* Fixed lists' names 
* Missing package-lock files. 
* Added a dummy simulation example 
* add a bug fix 
* add another code style fix 
* Wiki: minor corrections 
* Release Notes: 2021-10-14 Build 0.94.0 (cleanup) 
* change the function name and add style fix 
* Minor fixes 
* #178: Compute: Function analyzers: Function parameters grid \- WIP 
* Added missing lock files 
* Grok connect: improve Hikari idle connections timeout 
* Chem: SA optimization and other minor improvements 
* Removed obsolete packages. 
* An examples of the JS function that returns multiple values (cool visuals for sensitivity analysis) 
* Compute: Lotka-Volterra model script 
* Seldon: Fetch metadata from the deployment's model (v >= 1.3) 
* Seldon: Fetch metadata from the deployment's model. Stringify outputs 
* Work in progress 
* Compute: Lorenz attractor script 
* #179 Chem: Optimize package.js, externalize OCL 
* #179 Chem: Optimize package.js, externalize OCL. Cleanup 
* Compute: Logistic map script 
* Useful additions 
* Marked demo scripts as models. 
* #180: Compute: UX harmonization \- created ModelCatalogView class, added a menu and a toolbox 
* Adjusted parameters 
* Reflect TypeScript errors caught by the build server 
* #181 Workers: Helpers and Wiki 
* Added 'onClick' and 'handleNode' functions to the ElementOptions type 
* #180: Compute: UX harmonization \- WIP 
* #53 Tutorials package: Virtual Screening 
* (Bug) Grok connect: "Detected both log4j-over-slf4j.jar AND slf4j-log4j12.jar on the class path" 
* Chem: OCL cell renderer speedup 
* Chem: structural alerts fix 
* Update linter config in package template 
* JS API: added DG.InputBase.forColumn(c) 
* (Bug) Grok connect: Connection doesn't close on exception 
* (Bug) Grok connect: Connection doesn't close on test 
* (Bug) JS API: DataFrame.onColumnNameChanged won't work 
* JS API: View.forObject(x) 
* First commit 
* Renamed to follow our file naming convention 
* Code cleanup and documentation. 
* Add #183: Common: Wizard 
* (Bug) Compute: Outliers. Delete with a "bin" does something else 
* Implemented custom blocks and styles 
* Added the dialog resolve on dialog close 
* Seldon: support CATEGORICAL, move to metadata 
* Compute: changed the appearance of buttons on outliers selection 
* Clinical Case: added check for missing sdtm columns for each view, removed hardcoded sdtm columns names \- put them into constants 
* Compute: renamed the files 
* Seldon: support CATEGORICAL, move to metadata. wip 
* A snippet that demonstrates different ways to create a viewer 
* Added #189: JS API: TreeView.fromItemCategories(items, tags[]) 
* Minor code cleanup 
* Clinical Case: updated Timelines property panel 
* Compute: changed the buttons appearance 
* (Bug) Compute: Outliers. Hit "Esc" corrupts the state, dataframe disappears 
* Preparation for model upload 
* Clinical Case: foxed legend in adverse events view 
* Removing redundant code 
* Compute: changed behavior on outliers selection dialog closing 
* Clinical Case: corrected typos in Readme and views help 
* (Bug) Core: Sensitivity analysis won't work (WIP)
* Compute: outliers selection fix on zoom 
* Compute: added changes reverting on non-OK exit 
* Compute: fixed texts 
* Compute: moved buttons to top 
* Compute: removed function search on autodetection 
* Compute: changed extractRows function on stddev rule 
* Compute: Lotka-Volterra Model fix 
* JS API: TabControl.onTabChanged 
* Peptides: minor fixes 
* JS API: View.getRibbonPanels 
* JS API: MultiView 
* JS API: PropertyGrid 
* Clinical Case: renamed biomarkers distribution and correlation matrix views 
* ModelCatalog: pkpd, popPK models added 
* Compute: relevant signatures for Maximal Volume and Fab Arm Exchange scripts 
* Compute: removed mathjs package 
* Clinical Case: updated summary property panel for case when trial is not registeed on clinicaltrials.gov 
* Compute: disabled ESLint for several files 
* Compute: added the editing feature to outliers groups 
* move chem space to analysis folder 
* Grok connect: add logging of open and idle connections 
* update default icon 
* Compute: removed inner dialog 
* Compute: added cell selection to outliers 
* Grok connect: global settings  
* Grok connect: connection closure test 
* Chem: commit to check porting 
* ModelCatalog: PKPD model implementation 
* (Bug) Semantic type detection takes too long sometimes 
* Ability to not load package functions (for debug purposes) 
* JS API: ui.input.form 
* remove null entry exception 
* datagrok-tools: version bump 
* #53 Tutorials package: show completed steps description in a tooltip 
* Chem: replaced usage of PubChemApi function 
* ChemScripts: fixes 
* JS API: Add Qnum.sign() \- returns the string representation of the qualifier 
* JS API: Rename Qnum.sign() -> Qnum.qualifier() 
* Histogram: Add ability to set bands transparency 
* #53 Tutorials package: show the app function as a toolbox command 
* Chem: added margins to structural alerts blocks 
* Minor wiki updates 
* Compute: ModelCatalog separate package is redundant 
* GitHub Actions: package-lock.json and package.json 
* GitHub Actions: Update packages/ChemScripts/package-lock.json 
* GitHub Actions: Update packages/ChemblAPI/package-lock.json 
* GitHub Actions: Update packages/ChemblBrowser/package-lock.json 
* GitHub Actions: Update packages/ChaRPy/package-lock.json 
* GitHub Actions: Update packages/Charts/package-lock.json 
* GitHub Actions: Update packages/CustomML/package-lock.json 
* GitHub Actions: Update packages/DSP/package-lock.json 
* GitHub Actions: Update packages/DrugBank/package-lock.json 
* GitHub Actions: Update packages/Impute/package-lock.json 
* GitHub Actions: Update packages/NLP/package-lock.json 
* GitHub Actions: Update packages/PubChemApi/package-lock.json 
* GitHub Actions: Update packages/Seldon/package-lock.json 
* GitHub Actions: Update packages/Notebooks/package-lock.json 
* GitHub Actions: Update packages/SimPKPD/package-lock.json 
* GitHub Actions: Update packages/TensorFlow.js/package-lock.json 
* GitHub Actions: Update packages/VDJtools/package-lock.json 
* GitHub Actions: install js-api dependencies 
* GitHub Actions: Update packages/Bio/package-lock.json 
* GitHub Actions: Update packages/Compute/package-lock.json 
* GitHub Actions: Update packages/MultiPlot/package-lock.json 
* GitHub Actions: Update packages/HitTriage/package-lock.json 
* GitHub Actions: Update packages/DevTools/package-lock.json 
* GitHub Actions: Update packages/ClinicalCase/package-lock.json 
* GitHub Actions: Update packages/Chem/package-lock.json 
* GitHub Actions: Update packages/OpenChemLib/package-lock.json 
* GitHub Actions: Update packages/Peptides/package-lock.json 
* GitHub Actions: Update packages/SequenceTranslator/package-lock.json 
* GitHub Actions: Update packages/Tutorials/package-lock.json 
* Chem: 3d structure panel fix 
* GitHub Actions: Discovery and Forms packages 
* GitHub Actions: Update packages/Discovery/package-lock.json 
* GitHub Actions: Update packages/Forms/package-lock.json 
* Repertoire Broeser: commit for build 
* RepertoireBrowser package release 
* GitHub Actions: Update packages/RepertoireBrowser/package-lock.json 
* Fixed analyzer warnings 
* Clinical case: improved missing columns validation (in case column is optional view will be created) 
* Remove unnecessary ignore statements 
* Fix package build errors 
* GitHub Actions: Update packages/Viewers/package-lock.json 
* Add wu library #191 
* Chem: R-Group Analysis fix 
* Chem: R-Group columns categories fix 
* GitHub Actions: fix packages build 
* GitHub Actions: Update packages/UsageAnalysis/package-lock.json 
* js-api package-lock.json new version 
* Removing an extra line 
* Chem: To #PR 39. Invalidate categories 
* Peptides: RegEx catastrophic backtracking fix 
* Peptides: string trimming befor regexp testing 
* Chem: version bump 
* Compute: changed 'reason' term to 'rationale' 
* Update .gitignore. 
* Fix to pass linter rules of string length. 
* Add prototype of peptide space functionality. 
* Fix linter suggestions. 
* Chem version bump 
* Compute: npm update after wu package 
* Compute: removed dialog centering 
* RepertoireBrowser: initial refactor 
* RepertoireBrowser: fix 
* Enhance SPE class to be more general. Add UMAP method. 
* Peptides: detector fix 
* Peptides: followup cleanup 
* Docker, Useful additions 
* Compute: outliers detection fixes 
* Compute: version bump \- stable 
* JS API: DataConnection.test() 
* #53 Tutorials package: inspecting script results, script editing, minor adjustments 
* JS API: ability to "take a picture" of a viewer 
* RepertoireBrowser: port to ts initiated 
* Compute: documentation \- WIP 
* Chem: Switching to own fingerprints dictionary 
* Chem: Switching to own fingerprints dictionary. Switch implementations by imports 
* Oligo Batch Calculator: join representations 
* Add t-SNE embedding method. Add Jaro-Winkler string metric. 
* OligoBatchCalculator: spellcheck=false for text area 
* \* Add universal string measures. \* Add universal embedding methods. \* Add auxiliary types & mathematical  operations. 
* Added missing convertion 
* Compute: outlier selection UX improvement 
* Compute: patch version bump 
* RepertoireBrowser: porting works 
* Function View Harmonization 
* OligoBatchCalculator: fix of sequence validation bug 
* GitHub Actions: Update packages/OligoBatchCalculator/package-lock.json 
* #155 Chem: Package autotests. wip: npm run test 
* RepertoireBrowser: porting main to ts, excliding hartdcode 
* JS API: ui.input.forProperty 
* JS API: DG.Grid.fromProperties(items, properties) 
* JS API: ui.tableFromProperties 
* RepertoireBrowser: properties added 
* Undo changes in utils. 
* Add utils used in Peptides package. 
* Lines by equations: Extending lines so they can be drawn on different visualization independently 
* Clinical case: fixed sorting p-values array 
* RepertoireBrowser: ts porting 
* #155 Chem: Package autotests: wip. Local dataset 
* Compute: images export MVP 
* Compute: fixed imports 
* Compute: minor version bump 
* Compute: changed the plot's sheet 
* Modified import wu 
* REST API for functions 
* RepertoireBrowser: class MolecularLiabilityBrowser refactored and ported ts 
* Resolves #202 datagrok-tools: update package template 
* RepertoireBrowser: porting to ts 
* #155 Chem: Package autotests. wip: jest-worker-loader 
* Chem Benchmark. Add more test molecules 
* #181 Chem: Workers. wip: Restructure service and service workers 
* Viewers: harmonize property names for shapes 
* #155 Chem: Package autotests: wip. Substructure search 
* Peptides: split function rework 
* #155 Chem: Package autotests: wip. findSimilar 
* (Bug) connection.close in JdbcDataProvider.execute didn't account for exceptions other than SQLException 
* Jbio: add composition analysis view (WIP)
* #181 Chem: Workers. wip: Structural alerts, architecture 
* #181 Chem: Workers. wip: Structural alerts on an RdKit Service Worker 


# 2021-10-26 Dev build 0.94.1

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.94.1`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* JsEditor: show the result of the script execution, if possible 
* Code cleanup 
* (Possibly) fixed a bug when no network interfaces are present 
* Widget.temp property (WIP)
* \#135: DevTools: a lightweight package unit-testing framework 
* (Bug) Chem Filter: After reopening the filter, browser tab with Datagrok freezes (WIP)
* added time profile and adverse event browser views 
* Dart Utils: Minor code cleanup 
* Programming exercises (WIP)
* (Bug) Trellis Plot: Expanded viewer is only using a third of the horizontal space 
* \#53 Tutorials package: Viewers (EDA track) 
* Wiki: Restructure. Access (WIP)
* Sequence Translator: enhancements to AxoLabs display  
* Update tutorials widget 
* Peptides: various fixes and imporvements 
* Scatter Plot: Add tooltips for lines and bands 
* \#53 Tutorials package: code clean-up 
* JS API: Add method to easily customize lines & bands in a Scatter Plot 
* SequenceTranslator: OligoBatchCalculator: added image and links to YouTube presentations 
* \#53 Tutorials package: chem: Activity prediction (WIP) 
* Help & JS Examples: Add info about lines & band in a Scatter Plot (WIP)
* Wiki: Simplify Develop and Getting Started (WIP)
* updated AE browser property panel, fixed bug with property getter in Timelines viewer, fixed bug in Survival analysis view 
* removed dm from additional domains in AE browser property panel 
* (Bug) Chem: sdf file handler isn't invoked (WIP)
* Fixes #136: Client-side SDF file reader 
* Removed the sdf file reader from Chem (since it is now part of the OpenChemLib) 
* Peptides: minor improvements 
* Help: Change youtube video to screenshots (WIP)
* \#53 Tutorials package: chem: update datasets 
* \#132 Chem: Move Substructure Filter panel to the package 
* Removed debug printout 
* \#132 Chem: Move Substructure Filter panel to the package. Grooming 
* (Bug) Add New Column: Exception after opening dialog for editing an already added column 
* Caret: predict as factor 
* \#138 Chem: Text search for molecules in the sketcher. wip: SMILES only 
* Update CONTRIB.md, grooming 
* Updated UI in Time profile, Box plots, Matrix, Summary views. Fixed bugs in Survival analysis view 
* \#132 Chem: Move Substructure Filter panel to the package. Collaborative filtering 
* \#132 Chem: Move Substructure Filter panel to the package. Substruct library 
* \#139 Chem: Move from .tags to .temp 
* (Bug) Chem | Similarity Analysis: Execution process does not end. Errors in the console 
* Expose DataSourceCardView to TS-API (WIP)
* ModelHub initial commit 
* Peptides: various fixes and improvements WIP 
* added summary view description to readme 
* summary view screenshots 
* updated images links 
* added gif of summary view 
* added summary gif 
* set gif width and height 
* timelines view description 
* corrected typo 
* Fixed JS-API deploy 
* patient profile view description 
* correced markdown 
* All dependent columns should be recalculated after new formula was applied (WIP)
* Molecular Liability Browser: refactor to handle common actions in NGL and pViz 
* adverse events view description 
* Set version to 0.94.0 
* Grok connect: log query execution time by steps  
* Seldon: initial package commit 
* added laboratory view description 
* \#93 Peptides: fixed 2d view 
* Peptides: Barchart fix again && sar-list-viewer && little fixes 
* Peptides: filter fix and other improvements 
* Peptides: Sequence width 
* biomarkers distribution description, updated biomarkers distribution view(created ribbon panel) 
* \#53 Tutorials package: base class improvements 
* added correlation matix view description 
* added time profile description 
* Peptides: viewer row sorting and colors fix 
* added ae browser description 
* Peptides: Sar-list update 
* Revert "Peptides: Sar-list update" 
* Revert "Merge remote-tracking branch 'origin/master'" 
* Update Usage analysis package 
* Peptides: deleted commented code 
* updated boxplots view 
* Update icon 
* JS-API: ContextActions method 
* JS-API: star() method 
* JS-API: Accordion pane counts 
* JS-API: author getter 
* JS-API: Accordion.addTitle method 
* Renamed package 
* JS-API: ContextActions method (WIP)
* JS-API: star() method (WIP)
* JS-API: Accordion pane counts (WIP)
* JS-API: author getter (WIP)
* JS-API: Accordion.addTitle method (WIP)
* Expose DataSourceCardView to TS-API (WIP)
* (Bug) Grid: column reordering throws an error 
* added validation view description, fixed bug with time profile view 
* \#53 Tutorials package: add missing titles 
* updated ae browser gif 
* Clinical Case: update property panel 
* View: helpUrl setter 
* \#53 Tutorials package: chem: renaming 
* added p-values to box plots, updated survival analysis view 
* Added example of setUpdateIndicator. 
* Peptides: minor refactoring 
* Model Catalog stubs 
* Support Ad-Hoc queries in JS 
* (Bug) JS-API deploy ignores deleted files 
* Update localhost.docker-compose.yaml 
* Fixed typos and replaced oudated resources in Help/develop and Help/exercises sections. Added illustrations. 
* Update package icons 
* JBIO: all PTMs tracks syncing (WIP)
* Wiki: add missing heading for chem articles 
* Sequence Translator: Enhancements to SVG 
* JS API: Menu.root 
* \#53 Tutorials package: getMenuItem helper for visual hints 
* \#53 Tutorials package: entry function fix 
* \#53 Tutorials package: add visual hints for the top menu 
* Locked postgres image version 
* JBIO: PTMs syncing with paratopes and cdr 
* Remove form editor toolbar 
* JBIO: motif PTMs in own track and group 
* Peptides: major refactoring, viewers merging 
* Peptides: Refactor & nasty bugs fixed & selection mode for SARV 
* ChemblApi 
* \#150 Chem: Move to webpack and TS. wip 
* Peptides: minor fixes and improvements 
* added laboratory linechart with changes within normalized normal ranges, changed AE browser view to TableView, updated survival analysis view UI 
* added survival analysis description 
* Code cleanjup 
* JS-API: Packages test framework (WIP)
* Fixed naming 
* JS-API: AccordionPane.root property 
* JS API: Filters (WIP)
* JS API: Accordion.header 
* Pubchem API 
* Wiki: Using lists in parameterized queries 
* \#150 Chem: Move to webpack and TS. Moved to webpack 
* \#53 Tutorials package: add visual hints for the sidebar panes 
* Accordion: more documentation 
* Minor code cleanup. 
* Removed accidentaly added functions 
* Resolves #165: Folder Browser Preview Functions 
* Sidebar menu button (WIP)
* Delete the tutorials plugin 
* Formatting in Compute fixed 
* Added manual outliers detection 
* Peptides: removing lock 
* Fixed analyzer warnings 
* Sequence Translator: option to be able to  view just users’ pattern  
* Added eslint config and formatted all existing files 
* \#53 Tutorials package: chem: virtual screening 
* (Bug) Filters: date format is not correctly displayed in a tooltip 
* SequenceTranslator: minor code simplification 
* Oligo Batch Calculator: assume sequence representation by it's first valid code 
* datagrok-tools: `grok api` command 
* corrected typos in readme, updated ae browser property panel, adde labels to reference lines in laboratory linechart 
* Added generated dataframes return 
* Deleted tutorials plugin 
* Removed puppeteer dep 
* implemented lazy view loading 
* Context Help: allow loading external resources 
* PubChemApi: export PubChem class 
* Scripting: Support "id" annotation 
* (Bug) Chem: ChEMBL search throws an error for unknown molecules 
* ChemblApi: Iupac name function 
* Added function abort on cancel event 
* Added protection from multiple dialog opening 
* Added promise resolve on dialog closing 
* Added RemoveOutlier button 
* Added column existence check before the column insertion 
* Removed debug line, Fixed css 
* Harmonize scalar outputs in function view (WIP)
* Removed debug 
* Fixed backgroundColor option to be optional 
* Changed "detection" tem to "selection" term 
* Lasso tool is enabled by default 
* Fixed typo in ui.list example 
* Add demo files 
* added help for summary view 
* added help urls for views 
* moved seting of help urls to constructors 
* Lesser fixes for tutorials 
* Fixed scatterPlot width 
* Changed non-mutating outlier selection on mutating 
* updated laboratory help 
* Timelines: remove a hard-coded column name 
* (Bug) Package Properties: default values are ignored 

# 2021-10-14 build 0.94.0

## Major features and improvements

* 150 Chem: Move to webpack and TS. Fix wrong refactoring
* Ability to access FuncCall from function method
* Accordion: more documentation
* Add "iconTool" to devTools
* Add "Switch" control to JS API
* Add file handlers to function types
* Add NC: Open column selector when type $ in expression
* Add SAR viewer
* Add superclass to the Dialog widget
* Add ui.switch to examples and docs
* AppsViews and PackagesView: make cards smaller
* Bar chart: add collaborative filtering #85
* Bar chart: add int column as split and stack column
* Bar chart: don't handle a click on _catNamesBox in _inlineCategories mode, handle a click to the left of the bar as full bar click
* Barchart: Font changed to serif-sans, column numeration changed, current row highlighting added
* Bar chart: merge two adjacent filter rectangles
* Bar chart: remove barChartFilter bitset
* biomarkers distribution description, updated biomarkers distribution view(created ribbon panel)
* Can't execute .applyFormula on column without formula tag #95
* Caret: predict as factor
* Changed "detection" tem to "selection" term
* changed enrollent linechart to cumulative sum, implemented request to clinicaltrials.gov for study info, added split of AEs barcharts by treatment arm
* Changed non-mutating outlier selection on mutating
* Change endpoints to use new divided cvm images
* Charts: fix a package publication error
* Charts on viewers: Universal ability to show them with interactive legend
* Color-coding: support multiple color representations (hex, rgb) in tags for categorical and linear coloring (WIP)
* Column.ReplaceData method
* Column context actions: take semantic types into account
* Compute: Readme.md, continued.
* Compute package, initial commit.
* Connectors: Support Array input parameters #90
* Context Help: allow loading external resources
* Core: property for ordering menu items
* corrected typos in readme, updated ae browser property panel, adde labels to reference lines in laboratory linechart
* corrections of ui for survival view
* Cron: support @reboot schedule
* Dart Utils: Minor code cleanup
* DataFrame.onRowsFiltered event
* DataFrame.onSorted event
* datagrok-tools: adjust lint scripts for recursive directory check
* datagrok-tools: update package template
* datagrok-tools: `grok api` command
* Datagrok API and system samples, add cross linking to improve connectivity (WIP)
* DevTools: add data connection examples
* DevTools: annotate API examples
* DevTools: a polyfill for replaceAll
* DevTools: integration with Inspector
* exposed JSMol::get_smarts
* Fiexed misinformation in data access tutorial
* Files: package-specific AppData
* Filters: graceful handling of column removal
* Filters: proper detaching from events
* Filters wiki page update
* Fix Dart Analysis warnings
* Fixes #136: Client-side SDF file reader
* Fix grok_connect deploy after upgrade
* Fix misspelled class names
* Fix property initialization warnings
* Fix typos in function registration
* Fix UI tests (WIP)
* fized tooltips, yAxesLabels and linechart title in multiplot
* Form: move setReadOnly implementation into SketchHtmlElementHandler
* Form: onFormCreating event for adjusting the default column selection
* Formatting in Compute fixed
* Frontend track: WIP
* Funcs: Function cache for scripts and data sync
* Function parameter editor (WIP)
* Functions: Tabs in functions markup
* JS API: DataFrame.insert wrong order of optional argument
* Grok connect: add logging of memory and column type
* Grok connect: allow calls without mainCallId
* Grok Connect: cancelable queries
* Grok connect: debug queries using inspector and property panel
* Grok connect: don't show canceled query as an error
* Grok connect: log params in parameterized query
* Grok connect: log query execution time by steps
* Grok Connect: Update script
* Grok script: add timing to log
* Harmonize scalar outputs in function view (WIP)
* Help & JS Examples: Add info about lines & band in a Scatter Plot (WIP)
* Help: Change youtube video to screenshots (WIP)
* JS-API: Accordion.addTitle method
* JS-API: Accordion.addTitle method (WIP)
* JS-API: AccordionPane.root property
* JS-API: Accordion pane counts
* JS-API: Accordion pane counts (WIP)
* JS-API: author getter
* JS-API: author getter (WIP)
* JS-API: ContextActions method
* JS-API: ContextActions method (WIP)
* JS-API: Packages test framework (WIP)
* JS-API: star() method
* JS-API: star() method (WIP)
* JS API: "file edited" event (WIP)
* JS API: Ability to access FuncCall options and aux
* JS API: ability to set current object to SemanticValue
* JS API: Accordion.header
* JS API: add dialog's title
* JS API: add event args types to onViewerAdded/onViewerClosed
* JS API: add examples with the linear color-coding
* JS API: Add Legend widget
* JS API: Add method to easily customize lines & bands in a Scatter Plot
* JS API: add missing toJs in package.getProperties()
* JS API: add more specific types
* JS API: a type for property options
* JS API: BitSet.and, or, xor, andNot
* JS API: correct event args types
* JS API: Create percentile calculation function
* JS API: Create ui.makeDraggable()
* JS API: dataframe.ts cleanup.
* JS API: DG.Color.linear()
* JS API: Dialog.inputs
* JS API: file importers
* JS API: Filters (WIP)
* JS API: fix type warnings
* JS API: Functions: Added help for toString() function
* JS API: Functions: Added toString() function (https://community.datagrok.ai/t/type-detection-in-if-statement/589)
* JS API: grok.shell.view should look for a view instead of a table view
* JS API: grok.shell.windows.showRibbon flag
* JS API: Menu.root
* JS API: menu order parameter
* JS API: Minor harmonization
* JS API: Moved toJs from grok_api to ui.ts
* JS API: Move function initFormulaAccelerators from top level ui to tools or misc
* JS API: onContextMenuItemClick
* JS API: onDialogShown event
* JS API: optional parameter in ui.image()
* JS API: Pass InputBase to Dialog() insted of InpuBase.root
* JS API: ScatterPlot.onResetView
* JS API: Some improvements of info-bar
* JS API: Special logic of passing input to the ui.dialog()
* JS API: TabPane.header
* JS API: textual descriptions of currently applied filters (RowList.filters)
* JS API: ui.image, redirects on the new page when the link doesn't provided
* JS API: update MapChangeArgs interface
* JS API: update ScatterPlot event args types
* JS API Examples: Drag and Drop
* JsEditor: show the result of the script execution, if possible
* Lasso tool is enabled by default
* Legend: Ability to add elements to the legend associated with additional charts (splines)
* Line Chart: Ability to show extreme values (min/max, q1/q3. etc.) of each marker (WIP)
* Line Chart: Marker choice
* Modified Patient Profele View, added AE Risk Assessment View
* Molecular Liability Browser: refactor to handle common actions in NGL and pViz
* moved reading of validation rules table to clinicalCaseApp() function
* moved seting of help urls to constructors
* Multi linechart with variable number of charts
* MultiPlot: add category selection to the mulit linechart
* MultiPlot: added file extension
* MultiPlot: change color of marker of the status chart according to value and min and max columns
* MultiPlot: make status chart to display same data as muli linechart
* MultiPlot: manual categories setting in multi linechart
* MultiPlot: minor bugs fix and add advanced example to Readme.md
* MultiPlot: refactor and documenting (WIP)
* MultiPlot: switch to several tables
* New PowerPack grid styling
* PowerPack: Ability to hide welcomeView
* PowerPack: Described Power Widgets and Power Search
* Replaced ScatterPlot legend by Grid
* Scatter Plot: ability to show lines by equations
* Scatter Plot: Add tooltips for lines and bands
* Scatter Plot: Legend for markers
* Scatter Plot: Remove "Asterisk" marker (https://community.datagrok.ai/t/new-marker-by-scatterplot-feature-fails-rendering-with-many-data-points/586)
* Selenium: Add attributes "name" to elements of the new AddNewColumn dialog
* Set new AddNC Dialog as standart add/edit dialog
* SMARTS: test possible solutions for aromatization/kekulization
* Stacked barchart: Selection style changed, datagrok palette added, axis display simplified, unreadable bar labeles removed, trying to fix escaping bars bug.
* Substructure filter:undefined column & searchSubstructure changes fix
* Support Ad-Hoc queries in JS
* Support npm repositories
* Text-based file viewers: ability to edit and save files
* Timelines: remove a hard-coded column name
* updated boxplots view
* updated clinical demo files project
* updated images links
* updated laboratory help
* Updated public token
* updated survival analysis view
* updated survival analysis view and r scripts
* updated tooltips in timelines view, added dynamic change of lines width depending on zoom to timelines chart
* Updated UI in Time profile, Box plots, Matrix, Summary views. Fixed bugs in Survival analysis view
* updating survival plots on tab click in case filters changed
* UsageWidget: code cleanup
* Usage widget and new for UsageAnalysis queries
* View: helpUrl setter
* Viewers: Legend usage refactoring
* Viewers: Organize context menu entries consistently
* Viewers: Scatter Plot: Bands

## Bugs

* Add New Column: Column names in formula are case sensitive
* Add New Column: Exception after opening dialog for editing an already added column
* Add new package: impossible to add a package without description
* Apps: App button goes to common heap instead of separate tab
* Bar Chart: "by category" sorting throws an exception
* Bar chart: After switching on the "Relative Values" property, bars are displayed completely zoomed
* Bar chart: After switching on the "Relative Values" property, bars are incorrectly colored
* Bar chart: click with ctrl doesn't disable a selected bar
* Bar chart: Viewer view is not refreshed after changing aggregation
* Chem: ChEMBL search throws an error for unknown molecules
* Chem: External sketcher, filtering and SMARTS
* Chem: Pass SMARTS properly to searchSubstructure (WIP)
* Chem: sdf file handler isn't invoked (WIP)
* Chem Filter: After reopening the filter, browser tab with Datagrok freezes (WIP)
* Chem | Similarity Analysis: Execution process does not end. Errors in the console
* Core: FileInfo.readAsString returns undefined
* DataSync: Do not set .query tag if no query was loaded
* Exception on publishing NPM repo
* Filter resets when user recalculates column formula
* Filters: date format is not correctly displayed in a tooltip
* Filters: Histogram: synchronization issue
* Filters: Range sliders don't always appear on hover
* Fix CVM python environments through conda
* Form, Tile: rendering error
* Form: Molecules overlap each other after switching to another row
* Form: wrong representation of missing values for qualified numbers
* Grid: column reordering throws an error
* Grok connect: batch mode sometimes returns the first data frame instead of the last
* Grok connect: log min, max, mean values for column
* Grok connect: non-parameterized query is not logged
* Grok connect: query isn't logged on exception
* Histogram: Change position of the checkbox "Filter out missing values"
* Histogram: Exception if remove column when filter pane (with histogram) open
* JS-API deploy ignores deleted files
* JS: Property setter doesn't work on FuncCall
* JS API: "combo-popup.js" example doesn't work
* JS API: GroupByBuilder.where won't work with an object pattern
* JS API: Table.rows.filter() doesn't work if called after opening the table
* Legend: Does not appear when selecting a boolean column
* Line Chart: There is no Median in the list of aggregations on the chart
* Oligo Batch Calculator: grid doesn't appear if sequence was recognized as valid
* Package environment file deploy does not work
* Package Properties: default values are ignored
* Package publication deletes properties on second run
* Property Panel history is broken
* QNum Column: isNone check is not correct
* Scatter Plot: Wrong color definition for boolean column (blue and red).
* Sequence Translator: additional modifications are on the wrong side of antisense strand
* Tools: False breakpoints in VS Code debugging (WIP)
* Trellis Plot: Expanded viewer is only using a third of the horizontal space
* View.grid throws an exception if there is no grid in the layout
* Viewers: error when adding to a view with the restricted width

# 2021-10-27 Stable version 0.94.0 (autogenerated)

## Latest Docker Images

* Datagrok:
  *  `docker pull datagrok/datagrok:0.94.0`
  *  `docker pull datagrok/datagrok:stable`
* CVM:
  *  `docker pull datagrok/cvm:0.94.0`
  *  `docker pull datagrok/cvm:stable`

* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* JS API: Some improvements of info-bar
* Wiki: File Browsing and Sharing (WIP)
* Wiki: Simplify Develop and Getting Started (WIP)
* Wiki: Fix headings, links, and types; remove duplicates
* Wiki: Add missing links
* DevTools: integration with Inspector
* MultiPlot: switch to several tables
* Tooling: Actualize datagrok-tools for the onboarding
* JS API: Moved toJs from grok_api to ui.ts
* JS API: Create ui.makeDraggable()
* Multi linechart with variable number of charts
* Wiki: a docker compose possible issue.
* Wiki: "Defining semantic type detectors" reformatting and links
* Wiki: video-contents.md adjustments (cases, links, Dev Session 10, 11)
* Wiki: video-contents.md adjustments
* Wiki: links adjustments
* JS API Examples: Drag and Drop
* MultiPlot: manual categories setting in multi linechart
* JS API Examples update
* Update ui.md
* NLP readme.md update
* Set new AddNC Dialog as standart add/edit dialog
* MultiPlot: refactor and documenting (WIP)
* MultiPlot: make status chart to display same data as muli linechart
* Add new icons
* Create property-panel.js
* MultiPlot: change color of marker of the status chart according to value and min and max columns
* JS API: ui.image, redirects on the new page when the link doesn't provided
* (Bug) AddNC: Column names in formula are case sensitive
* datagrok-tools: update package template
* (Bug) View.grid throws an exception if there is no grid in the layout
* Datagrok API and system samples, add cross linking to improve connectivity (WIP)
* Documented the Color class.
* MultiPlot: add category selection to the mulit linechart
* Simplify info panels example
* Wiki: Dataframe (WIP)
* MultiPlot: minor bugs fix and add advanced example to Readme.md
* MultiPlot: added file extension
* Exercises: Simplify Exercises (WIP)
* Oligo Batch Calculator: add compatibility with Axolabs sequences
* Usage widget and new for UsageAnalysis queries
* \#49: Ketcher \- initial update
* Bar Chart \- minor code cleanup.
* Fixed type annotations
* JS API: ability to set current object to SemanticValue
* \#49: Ketcher \- work in progress
* Chem: ability to switch sketchers (OpenChemLib, Marvin, Ketcher, etc)
* Update the docs
* Update the links
* \#51: OpenChemLib package: initial update
* \#49, #51: Ketcher and OpenChemLib sketchers: API harmonization
* \#49, #51: Ketcher and OpenChemLib sketchers: code cleanup and harmonization
* \#49 Ketcher: import fix, add emitted js files to .gitignore
* \#51 OpenChemLib sketcher: fix external modules imports
* API: Column's .toString should return a boolean value
* Grok Connect: cancelable queries
* Update learning-widget.ts
* New PowerPack grid styling
* Function parameter editor (WIP)
* DevTools package
* Update scripting.md
* JS API: Pass InputBase to Dialog() insted of InpuBase.root
* Add "Switch" control to JS API
* Add NC: Open column selector when type $ in expression
* update ui.icons
* Create toolbar.js
* Update icons.js
* New widgets styling
* Update power-pack.css
* toolbox example
* Update toolbox.js
* \- if RdKitSubstructLibrary.search() is passed an invalid query, return no matches rather than crashing \- the renderer should not attempt to call is_valid() on a null mol \- the renderer should check if rdkitScaffoldMol is non-null and valid before passing it to generate_aligned_coords()   to avoid throwing unnecessary exceptions \- mols which are invalid should still be deleted to avoid leaking memory
* JS API: Move function initFormulaAccelerators from top level ui to tools or misc
* Removed the "shift-drag on axes" action description
* \#49, #51: Ketcher and OpenChemLib sketchers: first working version
* JS API: grok.shell.view should look for a view instead of a table view
* JS API: Special logic of passing input to the ui.dialog()
* Packages: JavaScript formatting
* Grok connect: allow calls without mainCallId
* CLI: documentation improvements
* Fixed the typings.
* (Bug) Sequence Translator: additional modifications are on the wrong side of antisense strand
* Introduction to cheminformatics
* DevTools: annotate API examples
* Wiki: Parameterized queries
* Widgets update
* Cheminformatics: prepare an intro for developers
* Minor cleanup
* Update README.md
* Frontend track: WIP
* added widget ordering
* image thumbnails for learn widget
* df.dialogs.addNewColumn fix
* update learn widget, add new css
* \#53 Tutorials package
* (Bug) Tools: False breakpoints in VS Code debugging (WIP)
* (Bug) Grok connect: batch mode sometimes returns the first data frame instead of the last
* Viewers: Legend usage refactoring
* Wiki: Query View: Updated documentation
* Wiki: Switch to a relevant CVM before the new version
* Tutorials package (WIP)
* Grok connect: don't show canceled query as an error
* Grok connect: debug queries using inspector and property panel
* (Bug) Package Properties: default values are ignored
* Update tooltip-events.js
* Fixed markup for some links to videos
* Fixed a bug that created two views of the same object
* PowerPack: Described Power Widgets and Power Search
* Top menu and Popup menu. Add ".endgroup()" method for closing the group of items
* (Bug) JS API: "combo-popup.js" example doesn't work
* Update Tabs 2 column.js
* JS API: ScatterPlot.onResetView
* \#53 Tutorials package: Visualization
* \#53 Tutorials package: code clean-up
* JS API: optional parameter in ui.image()
* \#53 Tutorials package: tutorial runner
* JS API: update ScatterPlot event args types
* (Bug) Grok connect: log min, max, mean values for column
* Grok connect: add logging of memory and column type
* Grok connect: log query execution time by steps
* Wiki: Extending and customizing Datagrok, simplifications
* \#53 Tutorials package: scripting tutorial (WIP)
* JS API: correct event args types
* \#53 Tutorials package: EDA track (WIP)
* [Scatter Plot: Add lasso info into help and to community](https://community.datagrok.ai/t/extensions-to-the-scatter-plot-viewer/481)
* JS API: textual descriptions of currently applied filters (RowList.filters)
* Create Custom-Machine-Learning.md
* Fix UI tests (WIP)
* added function for creation of AE risk assessment dataframe(relative risk, risk difference, odds ratio)
* Update and rename Custom-Machine-Learning.md to custom-ml-models.md
* JS API: file importers
* Add file handlers to function types
* resolves #80 Charts: WordCloud
* OligoBatchCalculator: README.md
* \#53 Tutorials package: Predictive Modeling tutorial (WIP)
* Add files via upload
* Wiki: JS API table of contents update
* Typedoc installed
* Merge commits
* Added customization options to wordcloud
* Enabled the word cloud viewer in Charts
* Fixed output variable names
* Modified Patient Profele View, added AE Risk Assessment View
* modified multiplot to use custom tables instead of tables frm grok.shell, added ability to use multiple timelins plots
* Connectors: Support Array input parameters #90
* updated clinical demo files project
* Script Editor, snippets menu
* Selenium: Add attributes "name" to elements of the new AddNewColumn dialog
* JS API: onDialogShown event
* Rename d4-input-editor class
* Add superclass to the Dialog widget
* (Bug) Grok connect: non-parameterized query is not logged
* Grok connect: log params in parameterized query
* Wiki: Peptides
* Update script-editor.ts
* Fiexed misinformation in data access tutorial
* JS API: add dialog's title
* JS API: Dialog.inputs
* Can't execute .applyFormula on column without formula tag #95
* chem:svgMol function -> molfile format support added
* DrugBank package -> search functions and widgets added
* Data on demand
* \#53 Tutorials package: Data connectors tutorial
* DrugBank package -> excessive header removed
* Added cache checkbox
* resolves #97 datagrok-tools: package content validation: check scripts location
* Fixed js-api links
* sdf file-reader
* Fix a typo
* \#53 Tutorials package: lock `rxjs` version
* sdf reader fix
* Update package.js
* \#53 Tutorials package: Multivariate analysis tutorial
* added survival analisys view
* \#53 Tutorials package: temporary ts-ignore for new functionality
* Wiki: misc small fixes.
* ApiSamples: a case typo in a sketcher sample.
* Wiki: initial Style Guide for Help contributors.
* Release: adding a versioned hand-crafted Release Notes.
* updated survival analysis view and r scripts
* moved reading of validation rules table to clinicalCaseApp() function
* Release Notes: 2021-05-10, accumulating. wip
* Release Notes: 2021-05-21, accumulating. wip
* Release Notes: 2021-06-10, accumulating. wip
* Update render-items.js
* Update custom-cell-rendering-indexes.js
* Release Notes: 2021-07-14, accumulating. wip
* Release Notes: 2021-07-29, accumulating. wip
* Release Notes: 2021-07-29, short annotation. wip
* Release Notes: 2021-07-29 Build 0.93.0
* \#53 Tutorials package: adjustments to new UI, code clean-up
* Release Notes: cleanup.
* Wiki: Update filter help
* corrections of ui for survival view
* Peptides: generate a dataset for further analysis
* (Bug) Grok connect: query isn't logged on exception
* Filters wiki page update
* updated survival analysis view
* Initial commit
* Fixed training
* Fixed typos
* \#104 Wiki: Library tour
* Add ui.switch to examples and docs
* added creation of survival dataset (for first SAE occurance), updated survival analisys view
* UI improvements in survival analysis view
* ApiSamples: fixed a typo.
* Docker Compose: a new CVM yaml file.
* Docker Compose: adding an old yaml for history.
* Docker Compose: removing the old yaml.
* Docker Compose: linkified the Docker yaml.
* Wiki: Custom machine learning models. Cleanup
* Wiki: Custom machine learning models. Links
* Charts: fix a package publication error
* JS API: update MapChangeArgs interface
* \#53 Tutorials package: change node placement
* \#53 Tutorials package: EDA: viewers (WIP)
* JS API: add event args types to onViewerAdded/onViewerClosed
* Oligo Batch Calculator: detect errors when code table contains one-char and two-chars codes
* Sequence Translator: Creating a mol file from the landing page
* Pairwise sequences alignment
* JS API: Create percentile calculation function
* Clarified what the percentile function is.
* added hospitalization and drug related ae to Kaplan-Meier endpoints, added ability to choose lab values on Patient profile
* fized tooltips, yAxesLabels and linechart title in multiplot
* Update package.ts
* Stacked barchart added
* Wiki: Capitalization
* Wiki: Capitalization. Adjustments for source codes with #
* Functions: Tabs in functions markup
* Cron: support @reboot schedule
* Decompose the calculated columns example
* API Samples: a few additional links
* \#110: Chem: dynamic sketchers
* Legend: Ability to add elements to the legend associated with additional charts (splines)
* GROK-GitHub-109 JS API: DataFrame.insert wrong order of optional argument
* Fix typos in function registration
* Fix misspelled class names
* Check package content at the publishing step (WIP)
* GROK-GitHub-109 JS API: DataFrame.insert wrong order of optional argument. Type fix
* API Samples: simplify an info panel example
* Wiki: Small rendering adjustments
* JS API: grok.shell.windows.showRibbon flag
* JS API: dataframe.ts cleanup.
* null
* Bar chart: add collaborative filtering #85
* (Bug) Filters: Histogram: synchronization issue
* Scatter Plot: Legend for markers
* Updated public token
* aligned sequence split function
* fixed splitAlignedPeptides function
* added splitAlignedSequence panel function
* closes #113 datagrok-tools: package template: add eslint config file
* Update CONTRIB.md
* Peptides: compile the initial dataset, detect sequence, split
* Added a demo file with peptide sequences.
* (Bug) Form, Tile: rendering error
* rewrite Peptides package in TypeScript
* Add SAR viewer
* JS API: fix type warnings
* Peptides: refactor
* Peptides: simplified the aligned sequence detector
* Stacked barchart: Selection style changed, datagrok palette added, axis display simplified, unreadable bar labeles removed, trying to fix escaping bars bug.
* \#114: Bio: grid cell renderer for sequences
* Sequence: made it work again; code cleanup.
* Wiki: Restructure. Access (WIP)
* escaping bar bug is fixed
* JS API: BitSet.and, or, xor, andNot
* Peptides: drawing circles in a grid instead of numbers
* Peptides: viewer rework
* (Bug) JS API: Table.rows.filter() doesn't work if called after opening the table
* Form: move setReadOnly implementation into SketchHtmlElementHandler
* Minor improvements
* JS API: a type for property options
* Peptides: added activity scaling options
* Peptides: source table selection based on viewer select
* Peptides: refactor according to coding guidlines
* (Bug) Filters: Range sliders don't always appear on hover
* JS API: Minor harmonization
* Peptides: dataset fix
* AppsViews and PackagesView: make cards smaller
* (Bug) Histogram: Change position of the checkbox "Filter out missing values"
* added ability to select values for y axe in multiplot, updated patient profile and timelines views
* Peptide Cell renderer
* \#53 Tutorials package: ML (predictive modeling tutorial): missing values imputation
* PackageView and Package PP: show github URL
* Clinical Cases UI update
* (Bug) Add new package: impossible to add a package without description
* (Bug) Apps: App button goes to common heap instead of separate tab
* Rollback
* Peptides: added help to property panel, minor improvements
* Fixed views creation
* Fixed typing
* Peptide cell renderer
* Peptides: color scheme for sequences
* (Bug) Bar Chart: "by category" sorting throws an exception
* update clinical cases UI
* Column context actions: take semantic types into account
* Packages: change the API source
* DevTools: a polyfill for replaceAll
* Peptides: scaling method replace ln -> lg
* Fixed Analyzer warning
* Fixed typo
* Removed old dependency
* lock files
* Barchart: Font changed to serif-sans, column numeration changed, current row highlighting added
* Peptides: renderer fix
* Package Settings: don't show the settings editor when there are no properties to edit
* Peptides: added histogram and updated stats, scaling fix
* (Bug) Chem: External sketcher, filtering and SMARTS
* JS API: add missing toJs in package.getProperties()
* Wiki: package settings
* exposed JSMol::get_smarts
* Peptides: AA sem type detector
* (Bug) DataSync: Do not set .query tag if no query was loaded
* Peptides: property panel WIP
* Add "iconTool" to devTools
* Peptides: property panel fix
* Peptides: finished converting to TS, code cleanup
* Linter rules update
* Removed an unused dependency
* UsageWidget: code cleanup
* (Bug) Core: FileInfo.readAsString returns undefined
* \#117: PowerPack: PowerSearch: search templates \- added server-based collections of templates
* implemented ability to set min and max values in multiplot
* \#117: PowerPack: PowerSearch: build fix
* \#117: PowerPack: PowerSearch: search templates \- added package settings for templates paths
* implmented multiple lines on linechart using eChart graphs
* fixed multiplot tooltips bug
* Charts on viewers: Universal ability to show them with interactive legend
* implemented changing height of scatter plot with categories multiselect in multiplot
* \#53 Tutorials package: ML (predictive modeling tutorial): actions for the model training view
* updating survival plots on tab click in case filters changed
* Peptides: Aligned sequence cell renderer text properly aligned, readability increased
* Peptides: analyzePeptides fix
* Funcs: Function cache for scripts and data sync
* Peptides: color coding circles based on MAD
* Peptides: help fix
* DevTools: add data connection examples
* Bar chart: add int column as split and stack column
* \#53 Tutorials package: helpers for working with view inputs
* Update icons for packages
* Update package.png
* Update Laboratory view
* Peptides: color coding optimization
* \#53 Tutorials package: openViewByType method
* JS API: add examples with the linear color-coding
* Peptides: color-coding fix, split sequence join
* Closes #118 JS API: color-coding methods (Column.colors)
* \#118 JS API: color-coding methods (Column.colors)
* Peptides: Stacked barchart transfered to ts and added to package
* logo Viewer for Peptides
* Peptides: Fixed cell renderer column width
* Update about-widget.ts
* Fix property initialization warnings
* (Bug) Fix CVM python environments through conda
* DG.Package: methods for working with files
* Files: package-specific AppData
* (Bug) Bar chart: After switching on the "Relative Values" property, bars are displayed completely zoomed
* SMARTS: test possible solutions for aromatization/kekulization
* Tutorials package: give rights for editing
* Datlas: Build js-api on startup
* \#84: Scatter Plot: Ability to show min/max on axes where appropriate
* Minor code cleanup
* Core: property for ordering menu items
* JS API: DG.Color.linear()
* Bar chart: merge two adjacent filter rectangles
* Peptides: ui fixes
* Renamed the file to reflect the feature name
* Text-based file viewers: ability to edit and save files
* Viewers: Organize context menu entries consistently
* JS API: menu order parameter
* Bar chart: remove barChartFilter bitset
* JBIO: prepare color schemes for NGL viewer
* Peptides: various fixes and improvements
* (Bug) Bar chart: After switching on the "Relative Values" property, bars are incorrectly colored
* Peptides: refactoring
* datagrok-tools: adjust lint scripts for recursive directory check
* changed enrollent linechart to cumulative sum, implemented request to clinicaltrials.gov for study info, added split of AEs barcharts by treatment arm
* JS API: Add Legend widget
* tutorials update
* Peptides: typo fixed
* Peptides: Added AA cell-renderer, changed sequence cell renderer
* update icon
* \#53 Tutorials package: change the storage structure (refactoring, WIP)
* \#53 Tutorials package: a proper way to access a tutorial's track (refactoring, WIP)
* Closes #119 Chem: Alignment differs for same MolBlocks with different offsets
* \#53 Tutorials package: css files (refactoring, WIP)
* Closes #120 Chem: Rotation flicker when using alignment
* Peptides: AA render in sar viewer
* (Bug) Package environment file deploy does not work
* \#117: PowerPack: PowerSearch: search templates
* Support npm repositories
* update tutorials. Add css style
* Packages: Ignore JS files if TS present
* tutorials images
* update tutorials
* Closes #121 (+ #77). Chem: Support SMARTS
* \#121 Chem: Support SMARTS
* \#84: Scatter Plot: Property to show/hide extreme labels
* Fix grok_connect deploy after upgrade
* Peptides: histogram rework WIP
* \#53 Tutorials package: code clean-up, hint styling
* Peptides: histogram fix WIP
* Added user
* update imgages for tutorials
* \#122 JS API: RowGroup class
* Line Chart: Marker choice
* \#53 Tutorials package: openDialog method
* (Bug) Line Chart: There is no Median in the list of aggregations on the chart
* (Bug) Bar chart: Viewer view is not refreshed after changing aggregation
* (Bug) Legend: Does not appear when selecting a boolean column
* Peptides: histogram rework and fix
* (Bug) Scatter Plot: Wrong color definition for boolean column (blue and red).
* Bar chart: don't handle a click on _catNamesBox in _inlineCategories mode, handle a click to the left of the bar as full bar click
* Closes #123 Chem: Substructure search to intercept SMARTS if MolBlock fails
* \#123 Chem: Substructure search to intercept SMARTS if MolBlock fails
* (Bug) Histogram: Exception if remove column when filter pane (with histogram) open
* tutorials updates
* \#53 Tutorials package: cheminformatics track: descriptors
* Line Chart: Ability to show extreme values (min/max, q1/q3. etc.) of each marker (WIP)
* (Bug) Bar chart: click with ctrl doesn't disable a selected bar
* JS API: onContextMenuItemClick
* fixes #119
* Chem: property-panel WIP
* \#53 Tutorials package: ability to pass a cusom promise for an action, helper for selecting an item from the context menu, code clean-up
* Chem: Substructure filter
* (Bug) JS API: GroupByBuilder.where won't work with an object pattern
* Color-coding: support multiple color representations (hex, rgb) in tags for categorical and linear coloring (WIP)
* Docker Compose: Remove profiles
* (Bug) Filter resets when user recalculates column formula
* Change endpoints to use new divided cvm images
* Peptides: t-test, u-test and histogram changes WIP
* Compute package, initial commit.
* Compute: Readme.md, continued.
* JS API: add more specific types
* update tu
* Filters: graceful handling of column removal
* Revert "Merge branch 'master' of https://github.com/datagrok-ai/public"
* Pubchem: Search panels
* implemented legend and switching type of linechart in multiplot. Implemented study summary property panel
* HitTriage: initial commit
* HitTriage: documentation: work in progress
* Grok Connect: Update script
* Peptides: center aar
* Fixed analyzer warnings
* Peptides: various improvements WIP
* Peptides: polishing
* Peptides: add a real-like dataset for analysis
* updated tooltips in timelines view, added dynamic change of lines width depending on zoom to timelines chart
* Peptides: Cell renderer split sequence fix & Header barchart
* \#53 Tutorials package: tutorial cancelation (WIP)
* Peptides: verious fixes and improvements
* \#53 Tutorials package: set initial value in `clearRoot`
* \#53 Tutorials package: remove unnecessary status update
* Substructure filter:undefined column & searchSubstructure changes fix
* JBIO: prepare tooltips and update ux
* Peptides: fix header renderer
* Peptides: clipped all cells
* Peptides: minor fixes and improvements
* Grok script: add timing to log
* Peptides: BC viewer is now rendered in onCellRender
* Peptides: Snake case function name fixed
* Peptides: Header BC selection and current column highlighting & labels added
* Peptides: more minor fixes and improvements
* \#53 Tutorials package: make step description optional
* \#53 Tutorials package: show completed step descriptions on click
* Peptides:Cell renderers for aa and aligned sequences for advanced dataset & some colors for barchart
* Peptides:Fixed
* Peptides:Fixed Errors
* Peptides: moved logo viewer
* Peptides: u-test rework
* Peptides: small fixes and improvements
* Form: onFormCreating event for adjusting the default column selection
* Commented out BitSet.setFast()
* \#117: PowerPack: PowerSearch: code cleanup and refactoring
* Added a picture for dynamic molecular sketchers
* JS API: "file edited" event (WIP)
* Resolves #134
* DataFrame.onSorted event
* DataFrame.onRowsFiltered event
* \#53 Tutorials package: Predicitve Modeling (WIP)
* \#53 Tutorials package: Multivariate Analysis
* Peptides: Added indent & small fixes
* (Bug) Viewers: error when adding to a view with the restricted width
* Scatter Plot: Remove "Asterisk" marker (https://community.datagrok.ai/t/new-marker-by-scatterplot-feature-fails-rendering-with-many-data-points/586)
* JS API: Functions: Added help for toString() function
* JS API: Functions: Added toString() function (https://community.datagrok.ai/t/type-detection-in-if-statement/589)
* Peptides: Long modifications hidden & tooltip added
* \#53 Tutorials package: allow text selection in steps
* \#53 Tutorials package: track restructuring
* \#53 Tutorials package: track progress formatting
* property panel in Adverse events view, updated box plot view(added default box plots based on p-value, added strata choice)
* clinical case fix box-plots view
* JS API: TabPane.header
* added check for number of stratas in box plots view
* Ability to access FuncCall from function method
* JS API: Ability to access FuncCall options and aux
* Strict always on
* Column.ReplaceData method
* Better events in AddNewColumn
* \#53 Tutorials package: add missing help URL
* update tutorial widget
* \#53 Tutorials package: data connectors tutorial
* Fixed links
* \#53 Tutorials package: handle edge cases in openViewByType
* Peptides:Basic responsive header barchart
* Peptides:Fixed color coding rules for AAR && AA are now rendered with short names
* update tutorial images
* Fixed $130: Core: histogram lines are not colored according to category colors
* \#53 Tutorials package: Predicitve Modeling
* Removed debug printout.
* update widgets
* \#53 Tutorials package: viewer basics (WIP)
* Peptides:Fixed cell renderer
* Peptides: analyse and select the most appropriate subset in the "real@ dataset for tool performance
* Peptides: fixed L color and added full AA names
* JBIO: removing test files
* (Bug) Property Panel history is broken
* (Bug) Form: wrong representation of missing values for qualified numbers
* (Bug) QNum Column: isNone check is not correct
* Fix Dart Analysis warnings
* \#53 Tutorials package: viewer basics
* Filters: proper detaching from events
* (Bug) Chem Filter: After reopening the filter, browser tab with Datagrok freezes (WIP)
* \#53 Tutorials package: Filters (EDA track)
* Scatter Plot: ability to show lines by equations
* Peptides: sketcher-tooltip for aar & refactoring chem palette
* Peptides: Barchart fixes & tooltip added
* Peptides: Sar aar tooltip added
* Peptides: Aligned sequens are now aligned
* HitTriage WIP
* Viewers: Scatter Plot: Bands
* Peptides: Full sequence mol-graph widget
* Peptides: polishing 1 for user meeting
* Peptides: Aligned sequence are aligned for sure
* Better typing
* Fixed CSS issue
* (Bug) JS: Property setter doesn't work on FuncCall
* Package files: Remove path prefix
* CSS fix
* (Bug) Exception on publishing NPM repo
* update learning widget
* Peptides: Initial table columns are hidden& Full AAname added to tooltip
* Revert "Better typing"
* Peptides: Activity column unhidden
* Peptides: Barchart fix
* Peptides: help
* Peptides: fix the number of significant digits in statistics
* Peptides: minor polishing
* (Bug) Form: Molecules overlap each other after switching to another row
* (Bug) Add New Column: Exception after opening dialog for editing an already added column
* \#135: DevTools: a lightweight package unit-testing framework
* (Bug) Package publication deletes properties on second run
* PowerPack: Ability to hide welcomeView
* Added "molecules" demo dataset that could be used synchronously (grok.data.demo.molecules)
* (Bug) Chem: Pass SMARTS properly to searchSubstructure (WIP)
* JsEditor: show the result of the script execution, if possible
* Code cleanup
* (Possibly) fixed a bug when no network interfaces are present
* Widget.temp property (WIP)
* added time profile and adverse event browser views
* Dart Utils: Minor code cleanup
* Programming exercises (WIP)
* (Bug) Trellis Plot: Expanded viewer is only using a third of the horizontal space
* \#53 Tutorials package: Viewers (EDA track)
* Sequence Translator: enhancements to AxoLabs display
* Update tutorials widget
* Peptides: various fixes and imporvements
* Scatter Plot: Add tooltips for lines and bands
* JS API: Add method to easily customize lines & bands in a Scatter Plot
* SequenceTranslator: OligoBatchCalculator: added image and links to YouTube presentations
* \#53 Tutorials package: chem: Activity prediction (WIP)
* Help & JS Examples: Add info about lines & band in a Scatter Plot (WIP)
* updated AE browser property panel, fixed bug with property getter in Timelines viewer, fixed bug in Survival analysis view
* removed dm from additional domains in AE browser property panel
* (Bug) Chem: sdf file handler isn't invoked (WIP)
* Fixes #136: Client-side SDF file reader
* Removed the sdf file reader from Chem (since it is now part of the OpenChemLib)
* Peptides: minor improvements
* Help: Change youtube video to screenshots (WIP)
* \#53 Tutorials package: chem: update datasets
* \#132 Chem: Move Substructure Filter panel to the package
* Removed debug printout
* \#132 Chem: Move Substructure Filter panel to the package. Grooming
* Caret: predict as factor
* \#138 Chem: Text search for molecules in the sketcher. wip: SMILES only
* Update CONTRIB.md, grooming
* Updated UI in Time profile, Box plots, Matrix, Summary views. Fixed bugs in Survival analysis view
* \#132 Chem: Move Substructure Filter panel to the package. Collaborative filtering
* \#132 Chem: Move Substructure Filter panel to the package. Substruct library
* \#139 Chem: Move from .tags to .temp
* (Bug) Chem | Similarity Analysis: Execution process does not end. Errors in the console
* Expose DataSourceCardView to TS-API (WIP)
* ModelHub initial commit
* Peptides: various fixes and improvements WIP
* added summary view description to readme
* summary view screenshots
* updated images links
* added gif of summary view
* added summary gif
* set gif width and height
* timelines view description
* corrected typo
* Fixed JS-API deploy
* patient profile view description
* correced markdown
* All dependent columns should be recalculated after new formula was applied (WIP)
* Molecular Liability Browser: refactor to handle common actions in NGL and pViz
* adverse events view description
* Set version to 0.94.0
* Seldon: initial package commit
* added laboratory view description
* \#93 Peptides: fixed 2d view
* Peptides: Barchart fix again && sar-list-viewer && little fixes
* Peptides: filter fix and other improvements
* Peptides: Sequence width
* biomarkers distribution description, updated biomarkers distribution view(created ribbon panel)
* \#53 Tutorials package: base class improvements
* added correlation matix view description
* added time profile description
* Peptides: viewer row sorting and colors fix
* added ae browser description
* Peptides: Sar-list update
* Revert "Peptides: Sar-list update"
* Revert "Merge remote-tracking branch 'origin/master'"
* Update Usage analysis package
* Peptides: deleted commented code
* updated boxplots view
* Update icon
* JS-API: ContextActions method
* JS-API: star() method
* JS-API: Accordion pane counts
* JS-API: author getter
* JS-API: Accordion.addTitle method
* Renamed package
* JS-API: ContextActions method (WIP)
* JS-API: star() method (WIP)
* JS-API: Accordion pane counts (WIP)
* JS-API: author getter (WIP)
* JS-API: Accordion.addTitle method (WIP)
* Expose DataSourceCardView to TS-API (WIP)
* (Bug) Grid: column reordering throws an error
* added validation view description, fixed bug with time profile view
* \#53 Tutorials package: add missing titles
* updated ae browser gif
* Clinical Case: update property panel
* View: helpUrl setter
* \#53 Tutorials package: chem: renaming
* added p-values to box plots, updated survival analysis view
* Added example of setUpdateIndicator.
* Peptides: minor refactoring
* Model Catalog stubs
* Support Ad-Hoc queries in JS
* (Bug) JS-API deploy ignores deleted files
* Update localhost.docker-compose.yaml
* Fixed typos and replaced oudated resources in Help/develop and Help/exercises sections. Added illustrations.
* Update package icons
* JBIO: all PTMs tracks syncing (WIP)
* Wiki: add missing heading for chem articles
* Sequence Translator: Enhancements to SVG
* JS API: Menu.root
* \#53 Tutorials package: getMenuItem helper for visual hints
* \#53 Tutorials package: entry function fix
* \#53 Tutorials package: add visual hints for the top menu
* Locked postgres image version
* JBIO: PTMs syncing with paratopes and cdr
* Remove form editor toolbar
* JBIO: motif PTMs in own track and group
* Peptides: major refactoring, viewers merging
* Peptides: Refactor & nasty bugs fixed & selection mode for SARV
* ChemblApi
* \#150 Chem: Move to webpack and TS. wip
* added laboratory linechart with changes within normalized normal ranges, changed AE browser view to TableView, updated survival analysis view UI
* added survival analysis description
* Code cleanjup
* JS-API: Packages test framework (WIP)
* Fixed naming
* JS-API: AccordionPane.root property
* JS API: Filters (WIP)
* JS API: Accordion.header
* Pubchem API
* Wiki: Using lists in parameterized queries
* \#150 Chem: Move to webpack and TS. Moved to webpack
* \#53 Tutorials package: add visual hints for the sidebar panes
* Accordion: more documentation
* Minor code cleanup.
* Removed accidentaly added functions
* Resolves #165: Folder Browser Preview Functions
* Sidebar menu button (WIP)
* Delete the tutorials plugin
* Formatting in Compute fixed
* Added manual outliers detection
* Peptides: removing lock
* Sequence Translator: option to be able to  view just users’ pattern
* Added eslint config and formatted all existing files
* \#53 Tutorials package: chem: virtual screening
* (Bug) Filters: date format is not correctly displayed in a tooltip
* SequenceTranslator: minor code simplification
* Oligo Batch Calculator: assume sequence representation by it's first valid code
* datagrok-tools: `grok api` command
* corrected typos in readme, updated ae browser property panel, adde labels to reference lines in laboratory linechart
* Added generated dataframes return
* Deleted tutorials plugin
* Removed puppeteer dep
* implemented lazy view loading
* Context Help: allow loading external resources
* PubChemApi: export PubChem class
* Scripting: Support "id" annotation
* (Bug) Chem: ChEMBL search throws an error for unknown molecules
* ChemblApi: Iupac name function
* Added RemoveOutlier button
* Added function abort on cancel event
* Added column existence check before the column insertion
* Added protection from multiple dialog opening
* Added promise resolve on dialog closing
* Removed debug line, Fixed css
* Harmonize scalar outputs in function view (WIP)
* Removed debug
* Fixed backgroundColor option to be optional
* Changed "detection" tem to "selection" term
* Lasso tool is enabled by default
* Fixed typo in ui.list example
* added help for summary view
* added help urls for views
* moved seting of help urls to constructors
* Lesser fixes for tutorials
* Fixed scatterPlot width
* Changed non-mutating outlier selection on mutating
* updated laboratory help
* Timelines: remove a hard-coded column name
* Replaced ScatterPlot legend by Grid
* (Bug) Oligo Batch Calculator: grid doesn't appear if sequence was recognized as valid
* Added group delete button to the group row
* 150 Chem: Move to webpack and TS. Fix wrong refactoring


# 2021-10-10 Dev build 0.93.35

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.35`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* \#53 Tutorials package: make step description optional 
* \#53 Tutorials package: show completed step descriptions on click 
* Peptides:Cell renderers for aa and aligned sequences for advanced dataset & some colors for barchart 
* Peptides:Fixed 
* Peptides:Fixed Errors 
* Peptides: moved logo viewer 
* Peptides: polishing 
* Peptides: u-test rework 
* Peptides: small fixes and improvements 
* Form: onFormCreating event for adjusting the default column selection 
* Commented out BitSet.setFast() 
* \#117: PowerPack: PowerSearch: code cleanup and refactoring 
* Added a picture for dynamic molecular sketchers 
* JS API: "file edited" event (WIP)
* Resolves #134 
* DataFrame.onSorted event 
* DataFrame.onRowsFiltered event 
* \#53 Tutorials package: Predicitve Modeling (WIP) 
* \#53 Tutorials package: Multivariate Analysis 
* Line Chart: Ability to show extreme values (min/max, q1/q3. etc.) of each marker (WIP)
* Peptides: Added indent & small fixes 
* (Bug) Viewers: error when adding to a view with the restricted width 
* Scatter Plot: Remove "Asterisk" marker (https://community.datagrok.ai/t/new-marker-by-scatterplot-feature-fails-rendering-with-many-data-points/586) 
* JS API: Functions: Added help for toString() function 
* JS API: Functions: Added toString() function (https://community.datagrok.ai/t/type-detection-in-if-statement/589) 
* update tutorials 
* Peptides: Long modifications hidden & tooltip added 
* \#53 Tutorials package: allow text selection in steps 
* \#53 Tutorials package: track restructuring 
* \#53 Tutorials package: track progress formatting 
* property panel in Adverse events view, updated box plot view(added default box plots based on p-value, added strata choice) 
* clinical case fix box-plots view 
* JS API: TabPane.header 
* added check for number of stratas in box plots view 
* Ability to access FuncCall from function method 
* JS API: Ability to access FuncCall options and aux 
* Strict always on 
* Column.ReplaceData method 
* Better events in AddNewColumn 
* \#53 Tutorials package: code clean-up 
* \#53 Tutorials package: add missing help URL 
* update tutorial widget 
* \#53 Tutorials package: data connectors tutorial 
* Fixed links 
* \#53 Tutorials package: handle edge cases in openViewByType 
* Peptides:Basic responsive header barchart 
* Peptides:Fixed color coding rules for AAR && AA are now rendered with short names 
* update tutorial images 
* Fixed $130: Core: histogram lines are not colored according to category colors 
* \#53 Tutorials package: Predicitve Modeling 
* Removed debug printout. 
* JBIO: prepare tooltips and update ux 
* update widgets 
* \#53 Tutorials package: viewer basics (WIP) 
* Peptides:Fixed cell renderer 
* Peptides: analyse and select the most appropriate subset in the "real@ dataset for tool performance (WIP)
* Peptides: fixed L color and added full AA names 
* JBIO: removing test files 
* (Bug) Property Panel history is broken 
* (Bug) Form: wrong representation of missing values for qualified numbers 
* (Bug) QNum Column: isNone check is not correct 
* Fix Dart Analysis warnings 
* \#53 Tutorials package: viewer basics 
* Filters: proper detaching from events 
* (Bug) Chem Filter: After reopening the filter, browser tab with Datagrok freezes 
* \#53 Tutorials package: Filters (EDA track) 
* Scatter Plot: ability to show lines by equations 
* Peptides: sketcher-tooltip for aar & refactoring chem palette 
* Peptides: Barchart fixes & tooltip added 
* Peptides: Sar aar tooltip added 
* Peptides: Aligned sequens are now aligned 
* HitTriage WIP 
* Viewers: Scatter Plot: Bands 
* Peptides: Full sequence mol-graph widget 
* Peptides: polishing 1 for user meeting 
* Peptides: Aligned sequence are aligned for sure 
* Better typing 
* Fixed CSS issue 
* (Bug) JS: Property setter doesn't work on FuncCall 
* Package files: Remove path prefix 
* CSS fix 
* (Bug) Exception on publishing NPM repo 
* update learning widget 
* Peptides: Initial table columns are hidden& Full AAname added to tooltip 
* Revert "Better typing" 
* Peptides: Activity column unhidden 
* Peptides: Barchart fix 
* Peptides: help 
* Peptides: fix the number of significant digits in statistics 
* Peptides: minor polishing 
* (Bug) Form: Molecules overlap each other after switching to another row 
* (Bug) Add New Column: Exception after opening dialog for editing an already added column 
* \#135: DevTools: a lightweight package unit-testing framework 
* (Bug) Package publication deletes properties on second run 
* PowerPack: Ability to hide welcomeView 
* Added "molecules" demo dataset that could be used synchronously (grok.data.demo.molecules) 
* (Bug) Chem: Pass SMARTS properly to searchSubstructure (WIP)


# 2021-10-01 Dev build 0.93.34

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.34`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) Filter resets when user recalculates column formula 
* Compute package, initial commit. 
* Compute: Readme.md, continued. 
* JS API: add more specific types 
* update tutorials 
* update tu 
* Updated public token 
* Filters: graceful handling of column removal 
* Revert "Merge branch 'master' of https://github.com/datagrok-ai/public" 
* Pubchem: Search panels 
* implemented legend and switching type of linechart in multiplot. Implemented study summary property panel 
* Color-coding: support multiple color representations (hex, rgb) in tags for categorical and linear coloring (WIP)
* HitTriage: initial commit 
* HitTriage: documentation: work in progress 
* (Bug) Chem: External sketcher, filtering and SMARTS (WIP)
* \#53 Tutorials package: cheminformatics track: descriptors 
* Grok Connect: Update script 
* Peptides: center aar 
* Fixed analyzer warnings 
* Peptides: various improvements WIP 
* Peptides: polishing 
* Peptides: add a real-like dataset for analysis (WIP)
* updated tooltips in timelines view, added dynamic change of lines width depending on zoom to timelines chart 
* Peptides: Cell renderer split sequence fix & Header barchart 
* \#53 Tutorials package: tutorial cancelation (WIP) 
* Peptides: verious fixes and improvements 
* \#53 Tutorials package: set initial value in `clearRoot` 
* \#53 Tutorials package: remove unnecessary status update 
* Substructure filter:undefined column & searchSubstructure changes fix 
* JBIO: prepare tooltips and update ux (WIP)
* Peptides: fix header renderer 
* Peptides: clipped all cells 
* Peptides: minor fixes and improvements 
* Grok script: add timing to log 
* Peptides: BC viewer is now rendered in onCellRender 
* Peptides: Snake case function name fixed 
* Peptides: Header BC selection and current column highlighting & labels added 
* Peptides: more minor fixes and improvements 


# 2021-09-27 Dev build 0.93.33

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.33`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) Package environment file deploy does not work 
* \#117: PowerPack: PowerSearch: search templates 
* Support npm repositories 
* Packages: Ignore JS files if TS present 
* \#121 Chem: Support SMARTS 
* \#84: Scatter Plot: Property to show/hide extreme labels 
* Fix grok_connect deploy after upgrade 
* Added user 
* \#122 JS API: RowGroup class 
* Line Chart: Marker choice 
* (Bug) Line Chart: There is no Median in the list of aggregations on the chart 
* (Bug) Bar chart: Viewer view is not refreshed after changing aggregation 
* (Bug) Legend: Does not appear when selecting a boolean column 
* (Bug) Scatter Plot: Wrong color definition for boolean column (blue and red). 
* Bar chart: don't handle a click on _catNamesBox in _inlineCategories mode, handle a click to the left of the bar as full bar click 
* \#123 Chem: Substructure search to intercept SMARTS if MolBlock fails 
* (Bug) Histogram: Exception if remove column when filter pane (with histogram) open 
* Bar chart: remove barChartFilter bitset 
* Line Chart: Ability to show extreme values (min/max, q1/q3. etc.) of each marker (WIP)
* (Bug) Bar chart: click with ctrl doesn't disable a selected bar 
* JS API: onContextMenuItemClick 
* (Bug) JS API: GroupByBuilder.where won't work with an object pattern 
* Color-coding: support multiple color representations (hex, rgb) in tags for categorical and linear coloring (WIP)
* (Bug) Filter resets when user recalculates column formula (WIP)
* DG.Package: methods for working with files 
* Change endpoints to use new divided cvm images 


# 2021-09-27 Dev build 0.93.32

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.32`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Peptides: Added AA cell-renderer, changed sequence cell renderer 
* update icon 
* \#53 Tutorials package: change the storage structure (refactoring, WIP) 
* \#53 Tutorials package: a proper way to access a tutorial's track (refactoring, WIP) 
* Closes #119 Chem: Alignment differs for same MolBlocks with different offsets 
* \#53 Tutorials package: css files (refactoring, WIP) 
* Closes #120 Chem: Rotation flicker when using alignment 
* Peptides: AA render in sar viewer 
* \#117: PowerPack: PowerSearch: search templates 
* update tutorials. Add css style 
* tutorials images 
* update tutorials 
* Closes #121 (+ #77). Chem: Support SMARTS 
* (Bug) Chem: External sketcher, filtering and SMARTS (WIP)
* Peptides: histogram rework WIP 
* \#53 Tutorials package: code clean-up, hint styling 
* Peptides: histogram fix WIP 
* update imgages for tutorials 
* \#122 JS API: RowGroup class 
* \#53 Tutorials package: openDialog method 
* Peptides: histogram rework and fix 
* Closes #123 Chem: Substructure search to intercept SMARTS if MolBlock fails 
* tutorials updates 
* \#53 Tutorials package: cheminformatics track: descriptors 
* JS API: onContextMenuItemClick 
* datagrok-tools: adjust lint scripts for recursive directory check 
* \#53 Tutorials package: ability to pass a cusom promise for an action, helper for selecting an item from the context menu, code clean-up 
* \#53 Tutorials package: code clean-up 
* Chem: Substructure filter 
* Docker Compose: Remove profiles 
* (Bug) Filter resets when user recalculates column formula (WIP)
* Updated public token 


# 2021-09-23 Dev build 0.93.31

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.31`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues



# 2021-09-23 Dev build 0.93.30

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.30`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Peptides: refactoring 
* datagrok-tools: adjust lint scripts for recursive directory check 
* changed enrollent linechart to cumulative sum, implemented request to clinicaltrials.gov for study info, added split of AEs barcharts by treatment arm 
* JS API: menu order parameter 
* Viewers: Organize context menu entries consistently 
* JS API: Add Legend widget (WIP)
* tutorials update 
* Peptides: typo fixed 


# 2021-09-22 Dev build 0.93.29

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.29`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Renamed the file to reflect the feature name 
* JS API: DG.Color.linear() 
* Text-based file viewers: ability to edit and save files 
* Viewers: Organize context menu entries consistently 
* JS API: menu order parameter 
* Bar chart: remove barChartFilter bitset 
* JBIO: prepare color schemes for NGL viewer (WIP)
* Peptides: various fixes and improvements 
* (Bug) Bar chart: After switching on the "Relative Values" property, bars are incorrectly colored (WIP)


# 2021-09-22 Dev build 0.93.28

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.28`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Peptides: color-coding fix, split sequence join 
* Closes #118 JS API: color-coding methods (Column.colors) 
* \#118 JS API: color-coding methods (Column.colors) 
* Peptides: Stacked barchart transfered to ts and added to package 
* logo Viewer for Peptides 
* Peptides: Fixed cell renderer column width 
* Update about-widget.ts 
* Fix property initialization warnings 
* (Bug) Fix CVM python environments through conda  
* DG.Package: methods for working with files 
* Files: package-specific AppData 
* (Bug) Bar chart: After switching on the "Relative Values" property, bars are displayed completely zoomed 
* Fix UI tests (WIP)
* SMARTS: test possible solutions for aromatization/kekulization 
* Tutorials package: give rights for editing 
* Datlas: Build js-api on startup 
* \#84: Scatter Plot: Ability to show min/max on axes where appropriate 
* Minor code cleanup 
* Charts on viewers: Universal ability to show them with interactive legend 
* Core: property for ordering menu items (WIP)
* JS API: DG.Color.linear() 
* Bar chart: merge two adjacent filter rectangles 
* Peptides: ui fixes 


# 2021-09-21 Dev build 0.93.27

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.27`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Funcs: Function cache for scripts and data sync 
* \#53 Tutorials package: code clean-up 
* \#53 Tutorials package: helpers for working with view inputs 
* Update icons for packages 
* Update package.png 
* Update Laboratory view 
* Peptides: color coding optimization 
* \#53 Tutorials package: openViewByType method 
* JS API: add examples with the linear color-coding 
* Bar chart: add int column as split and stack column 


# 2021-09-20 Dev build 0.93.26

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.26`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Peptides: property panel WIP 
* Add "iconTool" to devTools 
* Peptides: property panel fix 
* Peptides: finished converting to TS, code cleanup 
* Linter rules update 
* Removed an unused dependency 
* UsageWidget: code cleanup 
* (Bug) Core: FileInfo.readAsString returns undefined 
* \#117: PowerPack: PowerSearch: search templates \- added server-based collections of templates 
* implemented ability to set min and max values in multiplot 
* \#117: PowerPack: PowerSearch: build fix 
* \#117: PowerPack: PowerSearch: search templates \- added package settings for templates paths 
* implmented multiple lines on linechart using eChart graphs 
* fixed multiplot tooltips bug 
* Charts on viewers: Universal ability to show them with interactive legend 
* implemented changing height of scatter plot with categories multiselect in multiplot 
* \#53 Tutorials package: ML (predictive modeling tutorial): actions for the model training view 
* updating survival plots on tab click in case filters changed 
* Peptides: Aligned sequence cell renderer text properly aligned, readability increased 
* Peptides: analyzePeptides fix 
* Funcs: Function cache for scripts and data sync 
* Peptides: color coding circles based on MAD 
* Peptides: help fix 
* DevTools: add data connection examples 
* Bar chart: add int column as split column 


# 2021-09-17 Dev build 0.93.25

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.25`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Peptide cell renderer 
* Peptides: color scheme for sequences (WIP)
* Cron: support @reboot schedule 
* (Bug) Bar Chart: "by category" sorting throws an exception 
* update clinical cases UI 
* Column context actions: take semantic types into account 
* Packages: change the API source 
* DevTools: a polyfill for replaceAll 
* Peptides: scaling method replace ln -> lg 
* Fixed Analyzer warning 
* Fixed typo 
* Removed old dependency 
* lock files 
* Barchart: Font changed to serif-sans, column numeration changed, current row highlighting added 
* Peptides: renderer fix 
* Updated public token 
* Package Settings: don't show the settings editor when there are no properties to edit 
* Peptides: added histogram and updated stats, scaling fix 
* (Bug) Chem: External sketcher, filtering and SMARTS (WIP)
* JS API: add missing toJs in package.getProperties() 
* Wiki: package settings 
* Peptides: AA sem type detector 
* (Bug) DataSync: Do not set .query tag if no query was loaded 


# 2021-09-16 Dev build 0.93.24

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.24`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) Form, Tile: rendering error 
* rewrite Peptides package in TypeScript 
* Add SAR viewer 
* JS API: fix type warnings 
* Peptides: refactor 
* Peptides: simplified the aligned sequence detector 
* Stacked barchart: Selection style changed, datagrok palette added, axis display simplified, unreadable bar labeles removed, trying to fix escaping bars bug. 
* \#114: Bio: grid cell renderer for sequences 
* Sequence: made it work again; code cleanup. 
* Wiki: Restructure. Access (WIP)
* escaping bar bug is fixed 
* JS API: BitSet.and, or, xor, andNot 
* Peptides: drawing circles in a grid instead of numbers 
* Peptides: viewer rework 
* (Bug) JS API: Table.rows.filter() doesn't work if called after opening the table 
* Form: move setReadOnly implementation into SketchHtmlElementHandler 
* Minor improvements 
* JS API: a type for property options 
* Peptides: added activity scaling options 
* Peptides: source table selection based on viewer select 
* Peptides: refactor according to coding guidlines 
* (Bug) Filters: Range sliders don't always appear on hover 
* JS API: Minor harmonization 
* Peptides: dataset fix 
* AppsViews and PackagesView: make cards smaller 
* (Bug) Histogram: Change position of the checkbox "Filter out missing values" 
* added ability to select values for y axe in multiplot, updated patient profile and timelines views 
* Peptide Cell renderer 
* \#53 Tutorials package: ML (predictive modeling tutorial): missing values imputation 
* PackageView and Package PP: show github URL (WIP)
* Clinical Cases UI update 
* Update package.ts 
* (Bug) Add new package: impossible to add a package without description 
* (Bug) Apps: App button goes to common heap instead of separate tab (WIP)
* Rollback 
* Peptides: added help to property panel, minor improvements 
* Fixed views creation 
* Fixed typing 


# 2021-09-14 Dev build 0.93.23

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.23`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Stacked barchart added 
* Histogram: Add legend feature 
* Wiki: Capitalization 
* Wiki: Capitalization. Adjustments for source codes with # 
* (Bug) Color coding: custom colors are not persisted (WIP)
* Functions: Tabs in functions markup 
* Cron: support @reboot schedule (WIP)
* Decompose the calculated columns example 
* API Samples: a few additional links 
* \#110: Chem: dynamic sketchers 
* Legend: Ability to add elements to the legend associated with additional charts (splines) 
* GROK-GitHub-109 JS API: DataFrame.insert wrong order of optional argument 
* Fix typos in function registration 
* Fix misspelled class names 
* Check package content at the publishing step (WIP)
* GROK-GitHub-109 JS API: DataFrame.insert wrong order of optional argument. Type fix 
* API Samples: simplify an info panel example 
* Wiki: Dataframe (WIP)
* Wiki: Small rendering adjustments 
* JS API: grok.shell.windows.showRibbon flag 
* JS API: dataframe.ts cleanup. 
* Update package.ts 
* null 
* Bar chart: add collaborative filtering #85 
* (Bug) Filters: Histogram: synchronization issue 
* Scatter Plot: Legend for markers 
* Updated public token 
* aligned sequence split function 
* fixed splitAlignedPeptides function 
* added splitAlignedSequence panel function 
* closes #113 datagrok-tools: package template: add eslint config file 
* Update CONTRIB.md 
* Peptides: compile the initial dataset, detect sequence, split   
* Added a demo file with peptide sequences. 
* (Bug) Form, Tile: rendering error 


# 2021-09-09 Dev build 0.93.22

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.22`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* added hospitalization and drug related ae to Kaplan-Meier endpoints, added ability to choose lab values on Patient profile 
* fized tooltips, yAxesLabels and linechart title in multiplot 
* Update package.ts 
* (Bug) Functions cannot be cached if called from server-side 
* Fixed tests 


# 2021-09-09 Dev build 0.93.21

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.21`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* JS API: Create percentile calculation function 
* Clarified what the percentile function is. 
* Histogram: Ability to add x-values of splines 


# 2021-09-08 Dev build 0.93.20

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.20`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Charts: fix a package publication error 
* JS API: update MapChangeArgs interface 
* \#53 Tutorials package: change node placement 
* \#53 Tutorials package: EDA: viewers (WIP) 
* JS API: add event args types to onViewerAdded/onViewerClosed 
* Selenium: Add attributes "name" to elements of the new AddNewColumn dialog 
* Oligo Batch Calculator: detect errors when code table contains one-char and two-chars codes 
* Sequence Translator: Creating a mol file from the landing page 
* Pairwise sequences alignment 
* Change endpoints to use new divided cvm images (WIP)
* Fixed warning 
* (Bug) Functions cannot be cached if called from server-side 
* Histogram: Ability to add x-values of splines 


# 2021-09-07 Dev build 0.93.19

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.19`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* corrections of ui for survival view 
* Form: Allow editing labels 
* Peptides: generate a dataset for further analysis 
* (Bug) Grok connect: query isn't logged on exception 
* \#101: Charts: Trellis plot: category slider width is suboptimal sometimes when rendering molecules 
* \#102: Remember window settings visibility (toolbox, property panel, etc) between sessions 
* Providing function name when throwing a function call validation error 
* Filters wiki page update 
* Wiki: Dataframe (WIP)
* updated survival analysis view 
* Defaulted showPropeties to true. 
* Css fix 
* (Bug) New “marker-by” scatterplot feature fails rendering with many data points 
* Fixed typos 
* \#104 Wiki: Library tour 
* (Bug) Job: CallFunc processing error 
* Wiki: Update filter help 
* Add ui.switch to examples and docs 
* added creation of survival dataset (for first SAE occurance), updated survival analisys view 
* UI improvements in survival analysis view 
* ApiSamples: fixed a typo. 
* Popup Menu: Ability to add icon groups 
* Scatter Plot: Icons for Marker Shapes in Popup menu 
* Docker Compose: a new CVM yaml file. 
* Docker Compose: adding an old yaml for history. 
* Docker Compose: removing the old yaml. 
* Docker Compose: linkified the Docker yaml. 
* Wiki: Custom machine learning models. Cleanup 
* Wiki: Custom machine learning models. Links 


# 2021-09-03 Dev build 0.93.18

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.18`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Datlas: add deploy logging 
* sdf reader fix 
* Update package.js 
* \#53 Tutorials package: Multivariate analysis tutorial 
* (Bug) Scatter Plot | Markers: "Gradient" type does not work 
* (Bug) Scatter Plot | Markers: When "Pie Bar Chart" type selected, an exception is thrown and the Scatter Plot starts to slow down 
* (Bug) Form: Ellipsis truncation is missing after applying a layout 
* (Bug) App's git URL is not well-formed 
* added survival analisys view 
* Form: Allow editing labels (WIP)
* Minor code cleanup 
* Filters: Allow to paste a list of values to a filter 
* \#53 Tutorials package: temporary ts-ignore for new functionality 
* Histogram: Add ability to set bands transparency 
* Wiki: misc small fixes. 
* ApiSamples: a case typo in a sketcher sample. 
* (Bug) Jobs are not saved correctly 
* Histogram: ability to using column stats in band rules 
* (Bug) Query Builder: Exception when deselecting table in the dialog (WIP)
* Wiki: initial Style Guide for Help contributors. 
* Release: adding a versioned hand-crafted Release Notes. 
* updated survival analysis view and r scripts 
* moved reading of validation rules table to clinicalCaseApp() function 
* Ability to define job output parameters 
* Release Notes: 2021-05-10, accumulating. wip 
* Release Notes: 2021-05-21, accumulating. wip 
* Release Notes: 2021-06-10, accumulating. wip 
* Update render-items.js 
* Update custom-cell-rendering-indexes.js 
* Release Notes: 2021-07-14, accumulating. wip 
* Release Notes: 2021-07-29, accumulating. wip 
* Release Notes: 2021-07-29, short annotation. wip 
* Bar chart: add collaborative filtering #85 
* Release Notes: 2021-07-29 Build 0.93.0 
* Multivariate Analysis -\- always show plots 
* \#53 Tutorials package: adjustments to new UI, code clean-up 
* Release Notes: cleanup. 
* Wiki: Update filter help 
* Revert "Bar chart: add collaborative filtering #85" 


# 2021-08-31 Dev build 0.93.17

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.17`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Package content validation: refine missing export statement check 
* Custom ML: add blob support for JS engines 
* Histogram: Default values for filter slider 
* Grid: horizontal scroll bar jumps to the left as the data frame is filtered 
* resolves #97 datagrok-tools: package content validation: check scripts location 
* DrugBank package -> excessive header removed 
* Histogram: reference distributions (GitHub #91) (WIP)
* Fix UI tests (WIP)
* Fixed js-api links 
* sdf file-reader 
* Fix a typo 
* \#53 Tutorials package: lock `rxjs` version 
* Histogram: Display bands (WIP)


# 2021-08-30 Dev build 0.93.16

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.16`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Scatter Plot: Change caption for Markers in popup menu (Type -> Shape) 
* \#53 Tutorials package: Data connectors tutorial 
* DrugBank package -> excessive header removed 
* Added cache checkbox 
* Grid: horizontal scroll bar jumps to the left as the data frame is filtered 
* Updatd public token 


# 2021-08-30 Dev build 0.93.15

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.15`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* JS API: Dialog.inputs 
* Can't execute .applyFormula on column without formula tag #95 
* chem:svgMol function -> molfile format support added 
* DrugBank package -> search functions and widgets added 
* Data on demand (WIP)


# 2021-08-30 Dev build 0.93.14

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.14`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Script Editor, snippets menu 
* Form: Column Selector 
* \#53 Tutorials package: EDA track (WIP) 
* (Bug) Grid: Sometimes the horizontal scroll bar does not show 
* Selenium: Add attributes "name" to elements of the new AddNewColumn dialog 
* JS API: corrected layout of ui.switchInput 
* (Bug) Core: JavaScript wrapper for Dialog cannot be constructed 
* JS API: onDialogShown event 
* Rename d4-input-editor class 
* Add superclass to the Dialog widget 
* (Bug) Grok connect: non-parameterized query is not logged 
* Grok connect: log params in parameterized query 
* Grok connect: tests for logging 
* Update README.md 
* Wiki: Peptides 
* Update script-editor.ts 
* Fiexed misinformation in data access tutorial 
* Fixe unit test URL 
* Modal: update named parameters 
* (Bug) Connectors: 'Add new connection' returns a dialog for editing 
* JS API: add dialog's title 
* (Bug) Scatter Plot: marker bugs \- dont work caches, borders width, borders display, borders on other shapes (not circles), choices in popup menu, shapes are cut off, wrong hittests, wrong mouseover shapes, wrong sizes, bugs in Markercol functionality, 
* Connectors: Support Array input parameters #90 
* Data on demand (WIP)


# 2021-08-26 Dev build 0.93.13

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.13`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) Scatter Plot: Dot style chart stay on the top left corner of the canvas 
* AddNewColumn editor: allow changing the column name during editing 
* (Bug) Custom ML: task bar disappears when model applied 
* \#53 Tutorials package: Predictive Modeling tutorial (WIP) 
* Custom ML: apply the newly trained model 
* Wiki: JS API table of contents update 
* Typedoc installed 
* Merge commits 
* Typedoc instead of JsDoc 
* Custom ML: add custom engines from JS functions 
* Change endpoints to use new divided cvm images (WIP)
* Typedoc styling 
* (Bug) Custom ML: duplicate accordion parameters 
* Custom ML: use column name provided by apply function 
* Added customization options to wordcloud 
* Enabled the word cloud viewer in Charts 
* Fixed output variable names 
* Modified Patient Profele View, added AE Risk Assessment View 
* modified multiplot to use custom tables instead of tables frm grok.shell, added ability to use multiple timelins plots 
* Connectors: Support Array input parameters #90 
* updated clinical demo files project 


# 2021-08-25 Dev build 0.93.12

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.12`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) Filters: Empty value filter option is missing 
* Scatter Plot: allow using multiple marker shapes 
* (Bug) Heat Map: Zooming out changes custom colors 
* (Bug) Grok connect: log min, max, mean values for column 
* (Bug) Filters: Clicking the reset button does not remove filtering in the table 
* Grok connect: add logging of memory and column type 
* Grok connect: log query execution time by steps  
* Wiki: Extending and customizing Datagrok, simplifications 
* Scatter Plot: lasso selector 
* \#53 Tutorials package: scripting tutorial (WIP) 
* JS API: correct event args types 
* \#53 Tutorials package: EDA track (WIP) 
* [Scatter Plot: Add lasso info into help and to community](https://community.datagrok.ai/t/extensions-to-the-scatter-plot-viewer/481) 
* Canceled previous commit fe4817c (commited by mistake) 
* JS API: textual descriptions of currently applied filters (RowList.filters) 
* Create Custom-Machine-Learning.md 
* (Bug) Grid: Sometimes the horizontal scroll bar does not show 
* Fix UI tests (WIP)
* Updated public token 
* (Bug) Send a stack from Datlas and Grok Connect exceptions  
* Update beta users list 
* (Bug) Line chart: Overview | Stacked Area Chart throws an exception 
* (Bug) Form: input highlighting is missing in the editing mode 
* Rename d4-input-editor class 
* SketchView: add missing tooltips 
* (Bug) Color coding: custom colors are not persisted 
* Update and rename Custom-Machine-Learning.md to custom-ml-models.md 
* (Bug) Scatter Plot: Marker Types in popup menu don't stay marked 
* (Bug) Bar chart: "Index out of range: index should be less than 1" when split by quarter 
* JS API: file importers 
* (Bug) Scatter plot: Row Source == All: non-filtered markers are rendered as "missing" 
* Add file handlers to function types 
* Fixed js-doc 
* Data on demand (WIP)
* resolves #80 Charts: WordCloud 
* OligoBatchCalculator: README.md 
* Scatter Plot: adjust current row / mouse over row appearance when the marker shape changes 


# 2021-08-20 Dev build 0.93.11

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.11`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Top menu and Popup menu. Add ".endgroup()" method for closing the group of items 
* DevTools package 
* (Bug) JS API: "combo-popup.js" example doesn't work 
* Update Tabs 2 column.js 
* Added local help docs to .gitignore (/xamgle/web/help/) 
* Adjust tooltip delay 
* JS API: ScatterPlot.onResetView 
* \#53 Tutorials package: Visualization 
* \#53 Tutorials package: code clean-up 
* Data on demand (WIP)
* Package content validation: export statement 
* JS API: optional parameter in ui.image() 
* \#53 Tutorials package: tutorial runner 
* JS API: update ScatterPlot event args types 
* (Bug) Bar chart: filter fails in a second bar chart  (WIP)
* (Bug) Datlas threads stuck unexpectedly 
* AddNewColumn: onMetadataChanged 


# 2021-08-18 Dev build 0.93.10

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.10`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* PowerPack: Described Power Widgets and Power Search 
* Document drawing events 
* Update README.md 
* (Bug) JSAPI: Can't pass null to the function because of validation 


# 2021-08-18 Dev build 0.93.9

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.9`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* datagrok-tools: update package template (WIP)
* \#53 Tutorials package 
* Functions: int (WIP)
* (Bug) Tools: False breakpoints in VS Code debugging (WIP)
* Pie Chart: improve default layout 
* (Bug) Grok connect: batch mode sometimes returns the first data frame instead of the last 
* Viewers: Legend usage refactoring 
* Wiki: Query View: Updated documentation 
* Tooling: Actualize datagrok-tools for the onboarding (WIP)
* Wiki: Switch to a relevant CVM before the new version 
* Core: Viewers: a mode with grayed-out column names for ColumnComboBox (WIP)
* Tutorials package (WIP)
* Grok connect: batch mode sometimes returns the first data frame instead of the last 
* (Bug) Can't select rows using match() with int  
* (Bug) "NullError: method not found: '_rowCount'" after query cancel 
* (Bug) "Throw of null." after query cancel 
* Grok connect: don't show canceled query as an error 
* Grok connect: debug queries using inspector and property panel 
* (Bug) Package Properties: default values are ignored 
* (Bug) Bar chart: filter fails in a second bar chart  
* (Bug) onColumnsAdded/onColumnsRemoved args don't have a list of columns 
* (Bug) Histogram: Show X Axis property hides X and Y axes 
* AddNewColumn: set the formula tag before adding a column to the dataframe 
* JS API: corrected layout of ui.switchInput 
* (Bug) Bar Chart: Loss of zoom upon selection 
* Datlas: WatchDog timer (WIP)
* Update tooltip-events.js 
* Bar Chart: Ability to drag chart after zooming 
* Update README.md 
* Fixed markup for some links to videos 
* Fixed a bug that created two views of the same object 


# 2021-08-13 Dev build 0.93.8

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.8`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Widgets update 
* Pie Chart: Auto determine whether to display labels on small segments 
* Cheminformatics: prepare an intro for developers (WIP)
* Minor cleanup 
* JS API: Fixed logic for "editor-for" annotation 
* Update README.md 
* Removed unused file. 
* Frontend track: WIP 
* Pie Chart: Sort types for pies 
* added widget ordering 
* image thumbnails for learn widget 
* df.dialogs.addNewColumn fix 
* Wiki: Dataframe (WIP)
* update learn widget, add new css 
* Change endpoints to use new divided cvm images (WIP)


# 2021-08-12 Dev build 0.93.7

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.7`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* DevTools: annotate API examples 
* Wiki: Parameterized queries (WIP)
* BatchCall function (WIP)
* (Bug) Data sync doesn't work if table was refreshed 
* (Bug) Queries: No use of input parameters in simple choice lists 
* Function parameter editor (WIP)


# 2021-08-12 Dev build 0.93.6

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.6`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Removed the "shift-drag on axes" action description 
* (Bug) Scatter Plot: Points are not displayed if you change the markers type from hamburger menu 
* \#49, #51: Ketcher and OpenChemLib sketchers: first working version 
* JS API: grok.shell.view should look for a view instead of a table view 
* (Bug) Scatter Plot: Color scheme does not change until you change the viewer size 
* Exercises: Simplify Exercises (WIP)
* JS API: Special logic of passing input to the ui.dialog() 
* Wiki: Dataframe (WIP)
* (Bug) Package Properties: default values are ignored (WIP)
* Packages: JavaScript formatting 
* Grok connect: allow calls without mainCallId 
* CLI: documentation improvements 
* (Bug) Color Coding | Conditional: For Scatter Plot inside the Trellis, colors are not turned on instantly, but only if you "pull" the corresponding plot 
* Fixed the typings. 
* (Bug) Sequence Translator: additional modifications are on the wrong side of antisense strand 
* Tooling: Actualize datagrok-tools for the onboarding (WIP)
* Introduction to cheminformatics 
* Pie Chart: Auto width for legend 


# 2021-08-10 Dev build 0.93.5

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.5`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) Pie Chart: If you add Shift, then area when hovering over the segments is displayed incorrectly 
* (Bug) Pie Chart: Chart is cut-off at the left if you add shift 
* New PowerPack grid styling 
* Function parameter editor (WIP)
* (Bug) Form: Title background remains after disabling the design mode 
* DevTools package (WIP)
* (Bug) Line Chart: "Default Line Color" property does not change the chart color 
* Update scripting.md 
* JS API: Pass InputBase to Dialog() insted of InpuBase.root 
* Add "Switch" control to JS API 
* Tooling: Actualize datagrok-tools for the onboarding (WIP)
* Add NC: Open column selector when type $ in expression 
* update ui.icons 
* Create toolbar.js 
* Update icons.js 
* New widgets styling 
* Update power-pack.css 
* toolbox example 
* Update toolbox.js 
* Change endpoints to use new divided cvm images (WIP)
* JS API: Move function initFormulaAccelerators from top level ui to tools or misc (WIP)
* Implement function scheduler 


# 2021-08-09 Dev build 0.93.4

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.4`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Wiki: Simplify Develop and Getting Started (WIP)
* Bar Chart: with the `include nulls` property off, the legend should not mention (no value) 
* Custom ML: show all columns as features and labels  
* Ability to adjust socket timeouts 
* MultiPlot: add category selection to the mulit linechart 
* Tooling: Actualize datagrok-tools for the onboarding (WIP)
* datagrok-tools: update package template (WIP)
* Functions: Categories order 
* (Bug) Conditional Coloring: Empty cells are coloured if all values in this column are empty 
* Merge branch 'master' of https://bitbucket.org/skalkin/reddata 
* Simplify info panels example 
* Viewers: render boolean values as categories on axes 
* (Bug) Grid: empty cells are colored in gray and cell borders are not visible for categorically and linearly colored columns 
* Wiki: Dataframe (WIP)
* (Bug) Bar chart: If you select the area with only the upper bars (by Alt+drag), the scroll bar will not appear 
* (Bug) Bar Chart: error when a boolean column is selected as a "stack" column 
* MultiPlot: minor bugs fix and add advanced example to Readme.md 
* MultiPlot: added file extension 
* Function parameter editor (WIP)
* (Bug) Bar chart: null values ​​are not indicated 
* Bar Chart \- minor code cleanup. 
* Exercises: Simplify Exercises (WIP)
* JS API: Focus for opened dialogs 
* (Bug) Grid: Error after changing the column width 
* Oligo Batch Calculator: add compatibility with Axolabs sequences 
* (Bug) Grid: scroll bar should be below the header (does not currently work for tall headers) 
* Package content validation: function parameters 
* Usage widget and new for UsageAnalysis queries 
* (Bug) Bar Chart: vertical scroll bar is visible as soon as you open a bar chart 
* \#49: Ketcher \- initial update 
* Fixed type annotations 
* JS API: ability to set current object to SemanticValue 
* (Bug) Pie Chart: Legend does not react to changes in the table 
* (Bug) Correlation Plot: When zooming Scatter Plot on the PP table becomes the current object 
* \#49: Ketcher \- work in progress 
* Chem: ability to switch sketchers (OpenChemLib, Marvin, Ketcher, etc) (WIP)
* Update the docs 
* Update the links 
* \#51: OpenChemLib package: initial update 
* \#49, #51: Ketcher and OpenChemLib sketchers: API harmonization 
* \#49, #51: Ketcher and OpenChemLib sketchers: code cleanup and harmonization 
* Bar Chart: leave the column combobox on the X axis empty when the aggregation is `count` 
* (Bug) Pie Chart: "Show Inner Percent" property works the other way around 
* \#49 Ketcher: import fix, add emitted js files to .gitignore 
* \#51 OpenChemLib sketcher: fix external modules imports 
* API: Column's .toString should return a boolean value 
* Change endpoints to use new divided cvm images (WIP)
* Fixed a test 
* Grok Connect: cancelable queries (WIP)
* Update learning-widget.ts 


# 2021-08-04 Dev build 0.93.3

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.3`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) View.grid throws an exception if there is no grid in the layout 


# 2021-08-04 Dev build 0.93.2

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.2`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) Bar chart: null values ​​are not indicated 
* MultiPlot: change color of marker of the status chart according to value and min and max columns 
* (Bug) Scatter plot \- showSizeSelector  
* JS API: ui.image, redirects on the new page when the link doesn't provided 
* (Bug) AddNC: Column names in formula are case sensitive 
* Update Index 
* datagrok-tools: update package template (WIP)
* (Bug) Bar chart: Exception if reduce sliders when zooming 
* Ability to adjust socket timeouts 
* (Bug) View.grid throws an exception if there is no grid in the layout 
* Datagrok API and system samples, add cross linking to improve connectivity (WIP)
* Documented the Color class. 


# 2021-08-03 Dev build 0.93.1

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.1`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* JS API: Some improvements of info-bar 
* Wiki: File Browsing and Sharing (WIP)
* Wiki: Simplify Develop and Getting Started (WIP)
* Wiki: Fix headings, links, and types; remove duplicates 
* Wiki: Add missing links 
* DevTools: integration with Inspector (WIP)
* MultiPlot: switch to several tables 
* Tooling: Actualize datagrok-tools for the onboarding (WIP)
* JS API: Moved toJs from grok_api to ui.ts 
* JS API: Create ui.makeDraggable() 
* Data on demand (WIP)
* Sockets diagram 
* Bar Chart: stack column selector 
* Multi linechart with variable number of charts 
* Wiki: a docker compose possible issue. 
* Wiki: "Defining semantic type detectors" reformatting and links 
* Updated beta users list. 
* Wiki: video-contents.md adjustments (cases, links, Dev Session 10, 11) 
* Wiki: video-contents.md adjustments 
* Wiki: links adjustments 
* JS API Examples: Drag and Drop 
* MultiPlot: manual categories setting in multi linechart 
* JS API Examples update 
* Update ui.md 
* NLP readme.md update 
* (Bug) JS API: FuncCall.setParamValue() don't work 
* Set new AddNC Dialog as standart add/edit dialog 
* (Bug) DDT: Error in determining quarters 
* MultiPlot: refactor and documenting (WIP)
* MultiPlot: make status chart to display same data as muli linechart 
* Add new icons 
* [Scatter Plot: zoom sliders for axes](https://community.datagrok.ai/t/extensions-to-the-scatter-plot-viewer/481) 
* Change endpoints to use new divided cvm images (WIP)
* (Bug) Bar Chart: an exception when a column with missing values is selected as "stack" 
* Updated index 
* Create property-panel.js 

# 2021-07-29 build 0.93.0

## Major features and improvements

* Serialization info panel for columns
* Ability to share files
* Ability to skip DF reading on server after Grok Connect request
* Add New Column: ColumnGrid Widget, minor performance improvements
* Add Repository: Add attribute "name" for Source Type selector
* Add string name indexing for columns
* Adde pubspec to speed-up packages resolving
* Bar Chart: adaptive font size
* Bar Chart: automatically zoom in to a reasonable number of categories in case of many categories
* Bar Chart: make scrolling smoother
* Binning Functions: BinBySpecificLimits
* Box Plot: better border color; showCategorySelector and showValueSelector properties
* Conditional color-coding: bring options from the property panel to column context menu
* Connectors: Redshift: Schema browsing
* Connectors: default schema to providers and use as condition for schema browsing
* Correlation plot: on cell click, show the corresponding scatter plot in the property panel 
* Grid: ShowVisibleColumnsInTooltip property; automatically pick up cell type from dataframe
* Histogram and BarChart: improved washed-out default colors
* JS API: Added grok.shell.startUri
* JS API: ScatterPlot: add viewport
* JS API: Stats.histogramsByCategories
* JS API: TabPane: new properties: root, content, parent
* JS API: Viewer.onContextMenu event 
* JS API: add event onAccordionConstructed and Accordion.context getter 
* JS API: add grok.shell.views and grok.shell.tableViews
* JS API: dialog.clear()
* JS API: ui.icons for commonly used icons
* JS API: ui.image
* JS API: ui.tools.setHoverVisibility(host, elements)
* Join tables: Don't include key fields by default
* Line Chart: MultiAxis: double-click (choosing a primary series) should have no effect
* Line Chart: move hit-testing functionality to series renderers
* Line Chart: support for multiple Y axes
* Math Functions: Log(arg1, arg2), Log10(), Median(), RandBetween(), Round10(), etc.
* Multi-value filters
* Query Builder dialog: Add "name" attributes to checkboxes for tables
* ScatterPlot custom renderer: disable overlay dot rendering on categories hover
* ScatterPlot custom renderer: changed back marker color; set default marker size to 10
* Scripting: Ability to set custom name for parameter
* Scripting: Ability to set postfix for parameter
* Scripting: Show execution results in the property panel
* Table View: Columns pane: add search
* Text Functions: RegExpExtract(), RegExpReplace(), StrFind(), etc.
* Trellis Plot: ability to enlarge individual in-trellis viewers
* UI: Ability to set dialog background
* UTC support in the datetime parser
* Viewers: implement 'dashboard' style for all standard viewers; support for styles
* datagrok-tools: allow skipping questions in `grok config`

## Bugs

* Add New Column: Dragging functions opens a drop-area for the table
* Add new column: Strings in functions cannot be enclosed in single quotes
* Bar Chart: Unexpected bar color change after filtering
* Bar Chart: Viewer coloring settings should take precedence over the grid coloring settings (WIP)
* Bar Chart: coloring cannot be disabled; coloring gets only applied when editing colors in a grid column
* Box Plot: an exception when stdev(value) = 0
* Box Plot: improve initial choice of value column (stdev > 0 if possible)
* Column format changes are not persisted in layouts
* Connectors: Impala: int32max instead of real values
* Connectors: Oracle: NullPointerException in DB table content
* Core: Bitset.falseCount returns the number of set bits
* Core: default tooltip config is no longer saved with layout
* Custom ML: apply function from different ML engine
* Data | Unpivot does not work
* DataQuery with choices throws an exception
* Dataframe: Detect column max significant digits in CSV loading
* Events: onViewerAdded and onViewerClosed are sent twice
* Excel import: an empty column is created as part of the dataframe
* Extra space breaks function annotation
* Filters: Multi-value filters does not turn off when corresponding checkbox is off
* Filters: There are no icons for sorting and searching on hover for Multi-value filters
* Filters: adding returns a filter for the column added previously instead of the currently selected one
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
* JS API: OnDialogClosed event fires twice
* JS API: Properties cannot be changed in JsViewer (JsViewer.props and JsViewer.setOptions result in errors)
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
* Query-driven dashboards: query controls do not show up when a project is open
* Radiobutton throws exception when created
* Scatter Plot: "axes follow filter" feature does not work
* Scatter Plot: Regression line appear without activation
* Share button doesn't work for projects
* View Layouts: Error balloon after deleting saved layout
* Viewers: Inconsistent column selection inside a viewer and its properties panel
* Viewers: textColor property misspelling
* Viewers: the legend colors are not synchronized
* Viewers: the menu item `Viewer` is not visible in uploaded projects
* Viewers: Fixes regression in 2D layout alignment of unknown origin
* grok.dapi.projects.where('bad filter').first() returns a Project instance with d == null

# 2021-07-29 Stable version 0.93.0 (autogenerated)

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.93.0`
  *  `docker pull datagrok/datagrok:stable`
* CVM: 
  *  `docker pull datagrok/cvm:0.93.0`
  *  `docker pull datagrok/cvm:stable`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* fixes regression in 2D layout alignment of unknown origin 
* Chem package version 0.9.8, bump. 
* BioSignals: make initial chart look similar to content of accordion panes 
* DevTools: add the `Snippets` pane 
* JS API: add Grid.getVisibleCells 
* Datagrok API and system samples, add cross linking to improve connectivity (WIP)
* BioSignals: Add url routing for Physionet records 
* (Bug) JS Samples: to-hierarchy.js example does not work 
* NGL viewer (WIP)
* UI: Ability to set dialog background (WIP)
* Fixed TS warning 
* Scripting: Ability to set custom name for parameter 
* Scripting: Ability to set postfix for parameter 
* JS: Ability to hide table params from FuncCall editor 
* Running an ulmo sample in PythonScripts via pip+conda env. 
* BioSignals: publish connection to Physionet demo files 
* Fixed warnings 
* SequenceTranslator: MM12 Sequence Conversion 
* JS API: add grok.shell.views and grok.shell.tableViews 
* SequenceTranslator: OP100 Sequence Conversion 
* Power Pack: Compare Columns 
* Fixed logins in some connections 
* RepertoireBrowser: MSA tabs (WIP)
* Add files via upload 
* Update chem_curator.py 
* UsageAnalysis: Added 'users by date' query 
* datagrok-tools: code clean-up 
* Viewer.serialize(includeDefaults): an option to include default values 
* update "Common Actions" section 
* JS API: add square brackets for columns 
* JS API:  add a possibility to use custom ScatterPlot renderer 
* ClinicalCase: harmonize menu items 
* ClinicalCase: updated description 
* JS API: add Accordion.removePane() 
* Updated package description. 
* SequenceTranslator: ABI Sequence Conversion 
* JS API: ui.info for a yellow info bar 
* ClinicalCase \- a plugin for analyzing SDTM-based clinical studies data (WIP)
* Update ui.ts 
* (Bug) JS API: ColumnListProxy:  Cannot convert a Symbol value to a number 
* Dataframe: Detect column max significant digits in CSV loading 
* ui.inputs: made the 'options' parameter optional 
* JS API: ScatterPlot: add viewport 
* PowerPack: rewrite Compare Columns to TypeScript 
* (Bug) Grok connect: java.lang.NullPointerException with pattern parameter 
* JS API: Stats.histogramsByCategories 
* Made Widget.root non-nullable. WIP. 
* Dialog should have ui-panel class 
* PowerPack: Compare Columns fixes 
* SequenceTranslator: Conversion to Axolabs sense/antisense strands 
* PowerPack: Compare Columns CSS fixes 
* Fix UI tests (WIP)
* Selenium: Add ui-test for Chord viewer 
* SequenceTranslator: flexible Axolabs sequence conversion 
* JS API: add scatterPlot.worldToScreen() method 
* PowerPack: Compare Columns: Add Column button 
* JS API: add events to the ColumnComboBox 
* chem curation description added 
* Demo of structural curation added 
* Minor style improvements 
* Fixed paths to images 
* JS API: Func.appy(params): Promise<TResult> 
* Widget-based Welcome page (WIP)
* JS API: Func.find(package, name, tags, returnType) 
* Added documentation for FASTA reader 
* Added tag for sequence field 
* Update README.md 
* Wiki: How to manipulate viewers (WIP)
* fetchProxy: better typing. 
* JS API: ui.link(text, url | action, tooltip?) 
* Power Pack: Community Widget 
* JS API: ui.iframe(src) 
* Power Pack: WebWidget 
* PLOT3D init 
* JS API: toDart(): ability for JS classes to define custom conversion 
* JS API: ui.icons for commonly used icons 
* JS API: ui.tools.setHoverVisibility(host, elements) 
* JS API: ui.image 
* Power Pack: Learn Widget (WIP)
* Global Search (WIP)
* Power Pack: Search: FunctionSearchProvider 
* Power Pack: Search: WikiSearch 
* Charts: Tree Viewer (WIP)
* info page polishing 
* chem Cureate: info page polishing 
* chem Curate: jpg files changed to png 
* PLOT3D: import three.js and show scene 
* Added beta flag to the old Usage analysis 
* (Bug) JS API: JsViewerHostCore is returned instead of Viewer instance 
* Packages: Curation, Beta (WIP)
* Power Pack: ServicesWidget 
* JS API: ui.cards: an easy way to build cards (WIP)
* PLOT3D: add look() function 
* PLOT3D: add functions and some refactor 
* Power Pack: Power Search: integrate with functions (WIP)
* PLOT3D: before refactor 
* PLOT3D: fix problem with getGeometry() and other 
* PLOT3D: delete unused files 
* JS API: dialog.clear() 
* WIP 
* PLOT3D: refactor 
* JS API: View.createByName('functions', {search: '#math'}): a way to create standard Dart views 
* Power Search: integrate with default views ("functions", "databases", "files", etc) 
* Charts: Timelines (WIP)
* Update icons.js 
* Code cleanup 
* Update ui.md 
* PLOT3D: working original mesh and shaders to WegGL2 
* Widgets: remember widget state in user settings (WIP)
* PLOT3D: working filter and fixed inverted normals 
* Power Search: community-curated, template-based, widget-driven search engine (WIP)
* PLOT3D: auto rotation 
* Calculation of actvity cliffs 
* Update ui.md and add new ui samples 
* Wiki: Onboarding Content (WIP)
* Fixed indentation 
* UsageAnalysis harmonization (WIP)
* Reordering ui samples files 
* JS wrap for activity cliffs script 
* JS API: add ui.markdown 
* PLOT3D: transparecy, some props, fixed filter 
* (Bug) Sequence Translator | Main: If you enter the wrong input value, then many balloons with stacktraces will appear 
* (Bug) Sequence Translator | AXOLABS: Name of an existing pattern is not displayed 
* Sequence Translator | AXOLABS: Add tooltips for field headers which do not fit completely 
* .multi-value-separated: added dot in front 
* NotificationsWidget 
* Power Pack: Recent Projects Widget 
* PLOT3D: added columns size and column color props and some clean code 
* PLOT3D: fix animation 
* Timelines: toggle zoom sliders 
* Timelines: add a heading visualizing data from LB domain (WIP)
* Added sudo prompt in datagrok-tools installation 
* Create markdown.js 
* (Bug) Connectors: Oracle: NullPointerException in DB table content  
* Extended meta package: an example 
* Wiki: link & typos correction 
* Power Search: deep integration with functions (WIP)
* Update generate-html-help.side 
* Deploy: test connections migration (WIP)
* Activity cliffs visual modification 
* Moved to table view 
* user friendly options added 
* About-widget 
* Power Search: JS integration 
* JS API: fluent API for plotting (WIP)
* Sequence Translator: detect font color suitable for background color 
* SequenceTranslator: svg color adjustments 
* Oligo Batch Calculator: create initial version and test on files 
* (Bug) Connectors: Impala: int32max instead of real values 
* Binning Functions: BinBySpecificLimits 
* bug with table view fixed 
* JS API: ability to read metadata for a package or a group of packages 
* Renamed "UsageAnalysisUI" to "UsageAnalysis" 
* AutoGen code rebuild (pub get and build execution) 
* Renamed PowerPack to meet both npm and DG conventions 
* Simplified datagrok installation instructions 
* [(Bug) JS API: df.replace won't work on a column with the same name](https://community.datagrok.ai/t/annot-replace-a-column-with-a-column-with-the-same-name/528) 
* (Bug) Connectors: output table files for debug 
* JS API: Advanced Func API (WIP)
* Sequence Translator: add modifications to beginning of both strands 
* Cast exception to string 
* Connections: add system connections  
* (Bug) Filters: range slider filters out nulls 
* (Bug) OnDialogClosed event fires twice 
* (Bug) File Sharing doesn't work 
* Created index 
* *BREAKING\* Log function renamed to Ln 
* (Bug) JS Editor: exception when ApiSamples package is not there 
* Adjusted demo query name to 'System:CoffeeCompany:StoresInState' 
* JS API: ui.remove(element): a proper way to dispose of UI elements that might contain widgets 
* Power Search: integrate with projects 
* Math Functions: Log(arg1, arg2) 
* Math Functions: Log10() 
* Math Functions: Round10() 
* Correct mistake in dart installation instructions 
* JS API: add simpleMode and move presentationMode, hideTabsInPresentationMode to windows 
* SequenceTranslator: refactoring of svg creation 
* Grok connect: add CSV chars escaping  
* JS API: improve the HttpDataSource.include method to accept entity properties 
* Math Functions: RandBetween() 
* (Bug) Functions: function search won't work on packages 
* Updated public token 
* Fixed project deploy 
* (Bug) Usage analysis: persistent loader 
* JS API: format(x, fmt) 
* Power Pack: Widgets: Kpi 
* Power Search: integrate with widgets 
* Power Search (WIP)
* JS API: Added grok.shell.startUri 
* JS API: TabPane: new properties: root, content, parent 
* Power Pack: Search: routing 
* (Bug) An exception when closing the last view 
* Not opening a default view at startup in debug mode (in favor of PowerPack's home view) 
* Not opening a default view at startup when the URL starts with '/search' (in favor of PowerPack's search) 
* Math Functions: Avg() 
* JS API: PropertyBag.setAll(object) 
* DG.tools.createElementFromHtml 
* Power Pack: HtmlWidget 
* Introduced Property.editor (used to be Property.info.editor) 
* Fixed a typo. 
* JS API: Document function roles (WIP)
* Wiki: additional sorting examples 
* Wiki: how-to pages title consistency 
* Math Functions: Median() 
* Update ui examples 
* Oligo Batch Calculator: calculate optical density(OD), nmole and mass from Yield Amount & Units 
* Add color blindess-friendly color palettes  (WIP)
* Text Functions: StrFind() 
* Text Functions: StrLeft() 
* (Bug) JS API: code editor view differs from the original one (WIP)
* Text Functions: StrRight() 
* Text Functions: StrRepeat() 
* Text Functions: RegExpReplace() 
* Text Functions: RegExpExtract() 
* Create color-checker.js 
* Layout templates for JS API Examples 
* SequenceTranslator: adjustments to svg creation 
* Updated help, section transform \- text functions 
* (Bug) Project | Links: Project link is generated incorrectly 
* JS API: Viewer.onContextMenu event 
* Added reuseList function 
* Chem Activity cliffs \- visulisation polishing 
* Added solution for git push to public repo issue 
* Removed an annotation comment from samples opened in JS editor 
* (Bug) Grok Script: Syntax error: "abc".trim().toUpperCase()  
* Wiki: details on color-coding 
* (Bug) Deploy ignores friendlyName in queries 
* Fixed typo 
* Core: Add key missing calculation functions 
* Wiki: expand an article on adding new columns 
* PhyloTreeViewer: JS-based phylogenetic tree visualisation (WIP)
* (Bug) Share button doesn't work for projects 
* Line-height fix 
* (Bug) JS API: JsViewerLoader.instance sometimes returns null instead of jsViewer  
* (Bug) Notebook sharing doesn't work 
* Sequence Translator: set default translation examples 
* Oligo Batch Calculator: normalize UI 
* Usage analysis: performance improvements  (WIP)
* Update help index 
* JS API: Viewer.network() 
* Fixed scrollbars hover 
* Usage analysis: lists of errors (WIP)
* (Bug) Grok connect: can't browse scheme of CompoundLookup  
* pictures for learn widget 
* Widgets updates 
* Package Settings Editors (WIP)
* Remove 'Dev' info panel in favor of the package-defined one 
* DevTools package (WIP)
* Added function to test code performance 
* Update about-widget.ts 
* ApiSamples: example rearrangement 
* Operator: Xor() 
* DevTools: annotate API examples (WIP)
* Conversion functions: Boolean() 
* JS API: add ScatterPlot.onBeforeDrawScene 
* 'Serialization' info panel for columns 
* Added class 'link-external' 
* Projects: merge Log analysis with Usage analysis  
* JS API: ui.link options 
* Update recent-projects-widget.ts 
* Property Panel: add an 'Edit' button for adjusting the formula of a calculated column 
* Removed obsolte parameter 
* Add New Column: minor performance improvements 
* (Bug) Include method doesn't work from TS 
* (Bug) Oligo Batch Calculator: units and ratios don't match with cases of original app 
* (Bug) JS API: ui.table returns null on JnJ 
* Fixed "request a demo" button 
* Update submodule tokens 
* Fixed typing 
* Fixed strict checking 
* Change 'link-external' icon 
* Change editor icon 
* Oligo Batch Calculator: nearest neighbor method for extinction coefficient (EC) calculation​ 
* Ability to declare Script as App 
* Help indexing 
* Wiki: file exporter example 
* Formula engine: propagate semantic type of the top-level function to column 
* Nucleotide sequences for exercises: a-h1n1, sars-cov-2 (fixing non-GATC symbols) 
* isEmpty and isNotEmpty \- true if null or empty string 
* JS API: Column.editFormula() 
* Oligo Batch Calculator: remove empty rows and spaces 
* Programming exercises (WIP)
* Sequence Translator: possibility to add phosphate linkage before first nucleotide 
* Fix a typo 
* Prefix a CSS class 
* Prefix for the link-external class 
* Package Content Validation: allow refering to the common directory 
* Sequence Translator: enumerate translated sequences with unique IDs, linked to pattern name 
* RepertoireBrowser: tree dashboard 
* (Bug) Add new column: Strings in functions cannot be enclosed in single quotes 
* Documented a few methods, and marked some as obsolete 
* JS API: GridColumn.getVisibleCells() 
* Package property panel: renamed "Content" to "Functions" 
* Formula engine: Ability to reference columns 
* JS API: onBeforeRunAction, onAfterRunAction 
* Ability to return multiple outputs in JS 
* Sequence Translator: using existing pattern and saving it as own 
* Add statistical functions 
* Types.getCommonElementType \- provided a description (still needs a good implementation) 
* Minor changes to the environment setup instructions. Added a couple of the potentially useful scripts. 
* Harmonize file names (data_actions -> functions, data_action -> func) 
* Ddt: harmonize folder structure 
* (Bug) Functions: Unexpected behavior of some statistical functions 
* Removed duplicate import 
* Modified default config 
* Grok connect: deserialization fix 
* PepLens: table with sequences application initiated 
* Sequence Translator: overhang options 
* (Bug) Datlas skips row count check for package queries  
* (Bug) Exception on SetAutoComplete 
* (Bug) Deploy failed if public is not checked out 
* (Bug) Query ExtractParameters clears existing params 
* (Bug) View Layouts: When saving layout in an open project, new layout is not generated, but the project layout is overwritten 
* Ability to use other parameters from choices and suggestions statements 
* Better event firing 
* Removed labs submodule 
* Fixed analyzer warnings 
* Fixed submodule host 
* Oligo Batch Calculator: save as CSV 
* Ability to switch to advanced mode in Function View 
* Chem: activity cliffs ported to JS 
* (Bug) Exception when running script as app 
* view polishing 
* Fixed test 
* Grok connect: add batch mode  
* (Bug) Functions: Variable "row" don't work 
* (Bug) "JS functions with multiple output params must return JS object" exception when there are no output parameters 
* Got rid of DataFrameViewer.view in favor of ViewerBase.view 
* Dataframe: version getter for columns 
* Set default admin password for docker-compose  
* SequenceTranslator: fix of Shifts To Align Numbers Inside Circle 
* Sequence Translator: in modification sections enumerate nucleotides only 
* Updated EC2 instructions 
* Updated help index 
* MPLOT: init 
* Minor code harmonization 
* MultiPlot: add Close button "X" and functionality 
* SequenceTranslator: minor changes 
* JS API: Stats.uniqueCount 
* Sequence Translator: align non-overhang nucleotides on the right edge of both strands 
* JS: DataFrame.plot. methods 
* Grok connect: add helpers for dataframe substitution  
* Functions: Parameter categories 
* Box Plot: showCategorySelector and showValueSelector properties 
* (Bug) Box Plot: an exception when stdev(value) = 0 
* (Bug) Box Plot: improve initial choice of value column (stdev > 0 if possible) 
* Box Plot: better border color 
* Viewers: support for styles 
* Viewers: implement 'dashboard' style for all standard viewers 
*  MultiPlot: show/hide controls using ">" icon 
* (Bug) MultiPlot: contorls not removing when plot remove 
* Clinical Case: Timelines View (WIP)
* Clinical Case: Patient Profile View (WIP)
* Removed devops position 
* Timelines: remove the additional visulization in favor of MultiPlot 
* (Bug) Functions View: Click on the selected category does not remove the check mark 
* Updated copyright messages. 
* JS: grok.functions.call should return JS object for multiple output parameters 
* Oligo Batch Calculator: create separate functions for every indicator 
* SequenceTranslator: check for characters/words outside of allowed set, distinguish rows (WIP)
* MultiPlot: clean code 
* (Bug) MultiPlot: fix package.js 
* (Bug) Extra space breaks function annotation 
* (Bug) grok.dapi.projects.where('bad filter').first() returns a Project instance with d == null 
* JS API: Project.open()  // in workspace 
* JS API: Simplify registering functions that return futures/promises 
* Clinical Case: Adverse Events View (WIP)
* Histogram and BarChart: washed-out default colors 
* MultiPlot: remove hardcode and fix some remarks 
* Line Chart: move hit-testing functionality to series renderers 
* (Bug) Line Chart: "stacked bars" mode does not work 
* Introduced Color.textColor 
* Bar Chart: adaptive font size 
* MultiPlot: settings to properties 
* Bar Chart: in case of many categories, automatically zoom in to a reasonable number of categories 
* Bar Chart: make scrolling smoother 
* Edit Column Dialog: validators change 
* JS API: Column.applyFormula() 
* JS API: DataFrame.dialogs methods (WIP)
* Sequence Translator: final UI adjustments 
* Calculated columns examples 
* ScatterPlot: changed back marker color 
* ScatterPlot: set default marker size to 10 
* (Bug) Layouts: Grid looses event handlers after layout restore 
* DevTools: project snippets 
* (Bug) MultiPlot: controls overlap if last plot is hidden 
* (Bug) MultiPlot:  exception (empty grid) when only timeLine chart is shown 
* (Bug) Filters: Multi-value filters does not turn off when corresponding checkbox is off 
* (Bug) Filters: multi-value filters have no square indicator on top to toggle category selection 
* (Bug) Filters: There are no icons for sorting and searching on hover for Multi-value filters 
* (Bug) Filters: clicking on "search" should open search field AND focus on it 
* FuncCall.inputs.values() fix 
* Examples for onBeforeRunAction, onAfterRunAction 
* Remove DG.toDart wrapper from string keys 
* Add missing parameter types 
* OligoBatchCalculator: calculator function returns multiple values 
* (Bug) Chem: proper cache invalidation 
* Fixed analyzer warning 
* \- 2D layouts are now normalized wrt their orientation \- F, Cl, and S elements are more visible on a green selection background \- maximum size of molecules is clipped \- font/molecule size ratio is adjusted 
* fixes typo 
* DevTools: convert to TypeScript 
* Chem, RDKit-based (WIP)
* makes sure that categories are up-to-date after updating the Structure column 
* FunctionsWidget 
* implemented several validation rules, added validation view, added error summary to study summary view 
* MultiPlot: make tooltips as in DG 
* Update public token 
* (Bug) JsViewer root is not correctly wrapped 
* (Bug) JS API: properties cannot be changed in JsViewer (JsViewer.props and JsViewer.setOptions result in errors) 
* (Bug) JS API: JsViewer is no longer in the context of the onContextMenu event 
* (Bug) JS API: JsViewerHostCore is returned instead of Viewer instance 
* grok 9125 prepared for biowasm integration 
* added validatin for DM domain, implemented filtering of validatin result based on selected rule and ighlighting of cell with violated rule 
* added check of dm dataframe for null 
* removed extracting domain from sdtm-tag 
* implemented get_morgan_fp_as_uint8array() 
* JS API: correct viewer type names 
* DevTools: generate code for viewers customization 
* Oligo Batch Calculator: check for characters/words outside of allowed set, distinguish rows 
* fixed typo that was causing a regression (SMILES did not align to scaffold anymore) 
* Got rid of dapi.entityToJs 
* Fixed TS warnings 
* JS API: DataFrame.plot.fromType return type correction 
* (Bug) Scripts: "suggestions" and "choices" don't work 
* Revert "Normalize and prettify 2D layouts" 
* Same as #PR 39 but with a bug fix for failures in reduce 
* Adde pubspec to speed-up packages resolving 
* Custom ML functions 
* \- do not add molecules twice to substructure library \- no need to overwrite cells if SMILES 
* PowerPack: Removed Function widget from Welcome screen 
* (Bug) Grok connect: "The method execute() cannot take arguments" error with query parameters  
* (Bug) Grok connect: DriverManager returns wrong driver 
* (Bug) MS SQL: Schemas are empty 
* Added debug prints to  ExternalDataProvider 
* MultiPlot: add click support 
* MultiPlot: move timeLines render to external file 
* Added a polyfill for Object.entries (won't work in Dartium without it, and did not find a way to specify the relevant target to TypeScript) 
* (Bug) JS API: Custom object handlers should get inserted after the generic JsObjectMeta in order to have effect 
* JS API: TableView.syncCurrentObject to control whether current columns/rows should be shown in the property panel 
* Removed debug printouts. 
* (Bug) Matcher: matching on multiple criteria ignores the "and/or" option 
* Clinical Case: Adverse Events: show events preceding the adverse event 
* updated rules table, implemented opening validation view for specific domain by link 
* (Bug) Core: Bitset.falseCount returns the number of set bits 
* Untied dataFrame setter and onFrameAttached 
* Updated help index, minor changes 
* Line Chart: support for multiple Y axes 
* Line Chart: MultiAxis: double-click (choosing a primary series) should have no effect 
* (Bug) Line Chart: point hit-testing doesn't work for points on the right 
* Scripting: GetFunc: support maps 
* Removed unnecessary toString() overloads, all column subclasses now simply return column name. 
* (Bug) Line Chart: NullReferenceError when changing X axis column 
* UI: Text Input that can be decorated with an icons 
* CachedComputation: generic way to cache computations on mutable objects (such as columns) 
* (Bug) DataFrame.onValuesChanged fires twice when user edits a value in the grid 
* Line Chart: do not show aggregation options if all X points are unique 
* (Bug) Scatter Plot: "axes follow filter" feature does not work 
* (Bug) Line Chart: "axes follow filter" feature does not work 
* Removed custom DataFrame property 
* Updated public token; Added beta-user 
* Add New Column: ColumnGrid Widget 
* Modified dev documentation 
* JS API: Added Extended StringInput 
* Sequence Translator: add possibility to copy sequences from HTML table 
* Sequence Translator: Modify sequence length of saved pattern 
* UI: GridColumn, bug with double lines at the bottom of searchbox 
* OligoBatchCalculator: minor changes 
* JS API: Create ui.searchBox() 
* Add DB cache invalidation (WIP)
* (Bug) Viewers: Inconsistent column selection inside a viewer and its properties panel (WIP)
* Updated ui.md (ui.stringInput, ui.searchInput) 
* (Bug) Custom ML: apply function from different ML engine  
* (Bug) DataQuery with choices throws an exception 
* (Bug) JS API: Label breaks layout of TextInput with icon 
* Removed obsolete features. 
* Code cleanup. 
* (Bug) Query-driven dashboards: query controls do not show up when a project is open 
* JS API: DataFrame.mouseOverRow 
* MultiPlot: trim long categories 
* MultiPlot: remove underline from plot type combobox 
* JS API: tableInput 
* UI Test: Script View 
* UI Test: Welcome View (WIP)
* Add file browser description and instructions on sharing a file 
* Specify return types for dataframe.plot methods 
* (Bug) JS API: properties cannot be changed in JsViewer (JsViewer.props and JsViewer.setOptions result in errors) 
* Move File Browser documentation to access docs 
* Rename the screenshots for File Browsing and Sharing doc 
* MultiPlot: select on click 
* SketcherBase: a common interface for molecular sketchers 
* MultiPlot: pass series settings to viewer constructor 
* GROK \- 9242 UI Test: Repertoire Browser 
* Temporary disabled broken code 
* Updated index 
* TS: Add DayJs 
* (Bug) Typescript: Elements Block25 (50 and 75) are of incorrect type 
* Fixed typo in ui.ts 
* (Bug) JS Viewers: error: NullError: method not found: 'where$1' on null 
* Selenium: Share dialog: Add name attribute for privilege level selector 
* Query Builder dialog: Add "name" attributes to checkboxes for tables 
* (Bug) Predictive Models: Model in browser is not displayed on PP after clicking on card 
* ROK-8325 Data | Open Text: View does not switch to created table view after clicking on "Done" 
* Add Repository: Add attribute "name" for Source Type selector 
* (Bug) View Layouts: Error balloon after deleting saved layout 
* GridColumn: got rid of the Prop() annotations 
* Wiki: update developer guides to reflect JsViewer API changes 
* JS API: Add helpUrl support in ui.dialog 
* (Bug) Custom ML: feature names are not available in an apply script 
* (Bug) Custom ML: don't show dataframe from apply 
* Wiki: File Browsing and Sharing, tweaks 
* MultiPlot: rectangle select (WIP)
* dayjs: import statement fix 
* Sequence Translator: Drag&Drop file in Axolabs 
* datagrok-tools: update package template (WIP)
* Updated public 
* (Bug) JSViewer: TypeError: v.setOptions is not a function 
* Custom ML: add hyperparameters 
* Intersection types for props (ObjectPropertyBag & any) 
* Add sourcemap files to .gitignore 
* MultiPlot: restore status chart 
* (Bug) Grok connect: "No suitable driver found for..." 
* datagrok-tools: allow skipping questions in `grok config` 
* JS API: Added util function getUniqueName() 
* Ability to skip DF reading on server after Grok Connect request 
* MultiPlot: move scripts according to the structure convention 
* Sequence Translator: README documentation 
* layout templates (update) 
* add helpUrl example 
* (Bug) Grid: switching global coloring on / off removes linear color-coding 
* (Bug) Bar Chart: Unexpected bar color change after filtering 
* ObservableMap: granular change notification (action, value) 
* (Bug) Viewers: the legend colors are not synchronized 
* MultiPlot: fix types 
* JS API: Get color for cell / category 
* Disable NGINX logging (and other logging) in DG docker image 
* impoved validation view 
* (Bug) Viewers: textColor property misspelling 
* (Bug) Scatter Plot: Regression line appear without activation 
* MultiPlot: type casting 
* (Bug) Grok connect: NullPointerException with meta.cache option 
* JS API: Added ui.makeDroppable() and Dialog.initDefaultHistory() 
* UI: Corrected Column Grid Widget Mode view 
* Functions: Corrected some function tags 
* UI: Corrected Function Grid Widget Mode view 
* Power Pack: Add New Column (WIP)
* MultiPlot: run with pure script 
* Ability to turn on GrokConnect DataFrame recompression 
* (Bug) DataSync doesn't work with Grok Connect 
* Multi linechart (WIP)
* (Bug) Oligo Batch Calculator: Changes in one nucleotide chain affect the results for others 
* Fixed a test 
* Peptides: biowasm attachment for peptides sequence alignment (WIP)
* Oligo Batch Calculator: README documentation 
* (Bug) datagrok-tools: `grok add` is incompatible with TypeScript packages 
* Fixed cvm version 
* CVM: attempting to restore the last build state with pinning Ubuntu 16.04 version. 
* Reduce screenshots in file browsing and sharing doc 
* MultiPlot: change timeLInes columns order 
* Add New Column: History bug \- inputs are not filled with history 
* Add New Column: Dragging functions opens a drop-area for the table 
* Bar Chart: relative values 
* JS API: Rename MakeDroppable -> UI_MakeDroppable 
* MultiPlot: all categories trimmed 
* Updated public index 
* added laboratory view, added survival chart to summary view 
* Wiki: Bar Chart features 
* (Bug) Datagrok to Python skips blank lines 
* NPE fixed 
* (Bug) Exception when user applies layout without grid 
* (Bug) R-Group Analysis returns 500 
* (Bug) Chem | Substructure Search: "Pattern" field should open "Molecular sketcher" 
* Login error handling 
* (Bug) Chem: Molecular filter does not work (WIP)
* Wiki: File Browsing and Sharing (WIP)
* (Bug) JS API: Classification column binning to small decimal doesn't work as expected 
* MultiPlot: refactor and documenting (WIP)
* MultiPlot: filter in options (required in mulit linechart) (WIP)
* Remove unused code 


# 2021-07-27 Dev build 0.92.35

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.35`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Disable NGINX logging (and other logging) in DG docker image 
* Multi linechart (WIP)


# 2021-07-27 Dev build 0.92.34

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.34`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* DevTools package (WIP)
* (Bug) Scatter Plot: Regression line appear without activation (WIP)
* (Bug) Grok connect: NullPointerException with meta.cache option 
* JS API: Added ui.makeDroppable() and Dialog.initDefaultHistory() 
* UI: Corrected Column Grid Widget Mode view 
* Functions: Corrected some function tags 
* UI: Corrected Function Grid Widget Mode view 
* Power Pack: Add New Column (WIP)
* Updated public token 
* MultiPlot: run with pure script 
* Ability to turn on GrokConnect DataFrame recompression 
* (Bug) DataSync doesn't work with Grok Connect 


# 2021-07-26 Dev build 0.92.33

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.33`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* MultiPlot: fix types 
* JS API: Get color for cell / category 
* datagrok-tools: update package template (WIP)
* Disable NGINX logging (and other logging) in DG docker image 
* MultiPlot: restore status chart (WIP)
* Sequence Translator: README documentation (WIP)
* (Bug) Viewers: textColor property misspelling 
* (Bug) Scatter Plot: Regression line appear without activation (WIP)
* Ability to skip DF reading on server after Grok Connect request 
* MultiPlot: type casting 


# 2021-07-25 Dev build 0.92.32

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.32`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* layout templates (update) 
* add helpUrl example 
* (Bug) Grid: switching global coloring on / off removes linear color-coding 
* (Bug) Bar Chart: Unexpected bar color change after filtering 
* ObservableMap: granular change notification (action, value) 
* (Bug) Viewers: the legend colors are not synchronized 


# 2021-07-25 Dev build 0.92.31

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.31`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Ability to skip DF reading on server after Grok Connect request 


# 2021-07-24 Dev build 0.92.30

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.30`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* MultiPlot: restore status chart (WIP)
* DevTools package (WIP)
* Ability to skip DF reading on server after Grok Connect request 


# 2021-07-24 Dev build 0.92.29

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.29`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Ability to skip DF reading on server after Grok Connect request 


# 2021-07-23 Dev build 0.92.28

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.28`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Ability to skip DF reading on server after Grok Connect request 


# 2021-07-23 Dev build 0.92.27

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.27`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* DevTools package (WIP)
* Sequence Translator: README documentation (WIP)
* Ability to skip DF reading on server after Grok Connect request 


# 2021-07-23 Dev build 0.92.26

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.26`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Ability to skip DF reading on server after Grok Connect request 
* DevTools package (WIP)


# 2021-07-23 Dev build 0.92.25

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.25`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) Grok connect: "No suitable driver found for..." 
* Updated public token 
* datagrok-tools: allow skipping questions in `grok config` 
* JS API: Added util function getUniqueName() 
* Ability to skip DF reading on server after Grok Connect request 
* (Bug) JSViewer: TypeError: v.setOptions is not a function 
* MultiPlot: move scripts according to the structure convention (WIP)


# 2021-07-23 Dev build 0.92.24

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.24`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) JSViewer: TypeError: v.setOptions is not a function 
* JS API: tableInput 
* Custom ML: add hyperparameters 
* Intersection types for props (ObjectPropertyBag & any) 
* Add sourcemap files to .gitignore 
* MultiPlot: restore status chart (WIP)


# 2021-07-22 Dev build 0.92.23

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.23`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) Custom ML: feature names are not available in an apply script 
* (Bug) Custom ML: don't show dataframe from apply 
* Wiki: File Browsing and Sharing, tweaks 
* MultiPlot: rectangle select (WIP)
* JS API: Add helpUrl support in ui.dialog 
* MultiPlot: pass series settings to viewer constructor 
* dayjs: import statement fix 
* ClinicalCase \- a plugin for analyzing SDTM-based clinical studies data (WIP)
* Sequence Translator: Drag&Drop file in Axolabs 
* datagrok-tools: update package template (WIP)
* DevTools package (WIP)
* DevTools: annotate API examples (WIP)
* Updated public 


# 2021-07-21 Dev build 0.92.22

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.22`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* MultiPlot: select on click 
* SketcherBase: a common interface for molecular sketchers 
* Code cleanup 
* MultiPlot: pass series settings to viewer constructor (WIP)
* GROK \- 9242 UI Test: Repertoire Browser 
* Temporary disabled broken code 
* Updated index 
* TS: Add DayJs 
* (Bug) Typescript: Elements Block25 (50 and 75) are of incorrect type 
* Fixed typo in ui.ts 
* (Bug) JS Viewers: error: NullError: method not found: 'where$1' on null (WIP)
* Selenium: Share dialog: Add name attribute for privilege level selector 
* Query Builder dialog: Add "name" attributes to checkboxes for tables 
* (Bug) Predictive Models: Model in browser is not displayed on PP after clicking on card 
* ROK-8325 Data | Open Text: View does not switch to created table view after clicking on "Done" 
* Add Repository: Add attribute "name" for Source Type selector 
* (Bug) View Layouts: Error balloon after deleting saved layout 
* GridColumn: got rid of the Prop() annotations 
* Wiki: update developer guides to reflect JsViewer API changes (WIP)
* JS API: Add helpUrl support in ui.dialog 


# 2021-07-20 Dev build 0.92.21

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.21`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Add New Column: ColumnGrid Widget 
* Modified dev documentation 
* UI: Text Input that can be decorated with an icons 
* MultiPlot: move timeLines render to external file 
* Charts: Timelines (WIP)
* DevTools: generate code for viewers customization 
* JS API: Added Extended StringInput 
* Sequence Translator: add possibility to copy sequences from HTML table 
* DevTools package (WIP)
* DevTools: annotate API examples (WIP)
* Sequence Translator: Modify sequence length of saved pattern 
* UI: GridColumn, bug with double lines at the bottom of searchbox 
* OligoBatchCalculator: minor changes 
* JS API: Create ui.searchBox() 
* Add DB cache invalidation (WIP)
* (Bug) Line Chart: "axes follow filter" feature does not work 
* (Bug) Viewers: Inconsistent column selection inside a viewer and its properties panel (WIP)
* Updated ui.md (ui.stringInput, ui.searchInput) 
* (Bug) Custom ML: apply function from different ML engine  
* (Bug) DataQuery with choices throws an exception 
* (Bug) JS API: Label breaks layout of TextInput with icon 
* Removed obsolete features. 
* Code cleanup. 
* (Bug) Query-driven dashboards: query controls do not show up when a project is open 
* JS API: DataFrame.mouseOverRow (WIP)
* MultiPlot: trim long categories 
* MultiPlot: remove underline from plot type combobox 
* JS API: tableInput 
* UI Test: Script View 
* UI Test: Welcome View (WIP)
* Sequence Translator: final UI adjustments 
* Specify return types for dataframe.plot methods 
* (Bug) JS API: properties cannot be changed in JsViewer (JsViewer.props and JsViewer.setOptions result in errors) 


# 2021-07-18 Dev build 0.92.20

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.20`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Line Chart: support for multiple Y axes 
* Line Chart: MultiAxis: double-click (choosing a primary series) should have no effect 
* (Bug) Line Chart: point hit-testing doesn't work for points on the right 
* Scripting: GetFunc: support maps 
* Removed unnecessary toString() overloads, all column subclasses now simply return column name. 
* (Bug) Line Chart: NullReferenceError when changing X axis column 
* UI: Text Input that can be decorated with an icons 
* CachedComputation: generic way to cache computations on mutable objects (such as columns) 
* (Bug) DataFrame.onValuesChanged fires twice when user edits a value in the grid 
* Line Chart: do not show aggregation options if all X points are unique 
* (Bug) Scatter Plot: "axes follow filter" feature does not work 
* (Bug) Line Chart: "axes follow filter" feature does not work 
* (Bug) Filters: clicking on "search" should open search field AND focus on it 
* Removed custom DataFrame property 
* Updated public token; Added beta-user 


# 2021-07-16 Dev build 0.92.19

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.19`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Added debug prints to  ExternalDataProvider 
* MultiPlot: add click support 
* MultiPlot: move timeLines render to external file (WIP)
* Added a polyfill for Object.entries (won't work in Dartium without it, and did not find a way to specify the relevant target to TypeScript) 
* (Bug) JS API: Custom object handlers should get inserted after the generic JsObjectMeta in order to have effect 
* JS API: TableView.syncCurrentObject to control whether current columns/rows should be shown in the property panel 
* Removed debug printouts. 
* (Bug) Matcher: matching on multiple criteria ignores the "and/or" option 
* Clinical Case: Adverse Events: show events preceding the adverse event 
* Update ui.md 
* ClinicalCase \- a plugin for analyzing SDTM-based clinical studies data (WIP)
* updated rules table, implemented opening validation view for specific domain by link 
* (Bug) Core: Bitset.falseCount returns the number of set bits 
* Untied dataFrame setter and onFrameAttached 
* Updated help index, minor changes 


# 2021-07-15 Dev build 0.92.18

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.18`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) Grok connect: "The method execute() cannot take arguments" error with query parameters  (WIP)
* ClinicalCase \- a plugin for analyzing SDTM-based clinical studies data (WIP)
* (Bug) Grok connect: DriverManager returns wrong driver 
* (Bug) MS SQL: Schemas are empty (WIP)


# 2021-07-15 Dev build 0.92.17

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.17`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Adde pubspec to speed-up packages resolving 
* Custom ML functions (WIP)
* \- do not add molecules twice to substructure library \- no need to overwrite cells if SMILES 
* Oligo Batch Calculator: check for characters/words outside of allowed set, distinguish rows 
* ClinicalCase \- a plugin for analyzing SDTM-based clinical studies data (WIP)
* PowerPack: Removed Function widget from Welcome screen 


# 2021-07-15 Dev build 0.92.16

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.16`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Revert "Normalize and prettify 2D layouts" 
* Same as #PR 39 but with a bug fix for failures in reduce 
* Chem, RDKit-based (WIP)
* Updated public token 


# 2021-07-14 Dev build 0.92.15

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.15`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Sequence Translator: set default translation examples 
* Oligo Batch Calculator: normalize UI 
* DevTools: add the `Snippets` pane (WIP)
* Usage analysis: performance improvements  (WIP)
* Wiki: expand an article on adding new columns 
* Update help index 
* JS API: Viewer.network() 
* Fixed scrollbars hover 
* PhyloTreeViewer: JS-based phylogenetic tree visualisation (WIP)
* Core: Add key missing calculation functions 
* Usage analysis: lists of errors (WIP)
* (Bug) Grok connect: can't browse scheme of CompoundLookup  
* pictures for learn widget 
* Widgets updates 
* Package Settings Editors (WIP)
* Remove 'Dev' info panel in favor of the package-defined one 
* DevTools package (WIP)
* JS API: Document function roles (WIP)
* Added function to test code performance 
* Update about-widget.ts 
* ApiSamples: example rearrangement 
* Operator: Xor() 
* Datagrok API and system samples, add cross linking to improve connectivity (WIP)
* DevTools: annotate API examples (WIP)
* Conversion functions: Boolean() 
* JS API: add ScatterPlot.onBeforeDrawScene 
* 'Serialization' info panel for columns 
* Deploy: test connections migration (WIP)
* Added class 'link-external' 
* Projects: merge Log analysis with Usage analysis  
* JS API: ui.link options 
* Update recent-projects-widget.ts 
* Property Panel: add an 'Edit' button for adjusting the formula of a calculated column (WIP)
* Removed obsolte parameter 
* Add New Column: minor performance improvements 
* (Bug) Include method doesn't work from TS 
* JS API: improve the HttpDataSource.include method to accept entity properties 
* (Bug) Oligo Batch Calculator: units and ratios don't match with cases of original app 
* (Bug) JS API: ui.table returns null on JnJ 
* Fixed "request a demo" button 
* Update submodule tokens 
* Fixed typing 
* Fixed strict checking 
* Change 'link-external' icon 
* Change editor icon 
* Oligo Batch Calculator: nearest neighbor method for extinction coefficient (EC) calculation​ 
* Ability to declare Script as App 
* Help indexing 
* Wiki: file exporter example 
* Formula engine: propagate semantic type of the top-level function to column 
* Nucleotide sequences for exercises: a-h1n1, sars-cov-2 (fixing non-GATC symbols) 
* isEmpty and isNotEmpty \- true if null or empty string 
* JS API: Column.editFormula() 
* Oligo Batch Calculator: remove empty rows and spaces 
* Programming exercises (WIP)
* Sequence Translator: possibility to add phosphate linkage before first nucleotide 
* Fix a typo 
* Prefix a CSS class 
* Prefix for the link-external class 
* Package Content Validation: allow refering to the common directory 
* Sequence Translator: enumerate translated sequences with unique IDs, linked to pattern name 
* RepertoireBrowser: tree dashboard 
* (Bug) Add new column: Strings in functions cannot be enclosed in single quotes 
* Documented a few methods, and marked some as obsolete 
* JS API: GridColumn.getVisibleCells() 
* Package property panel: renamed "Content" to "Functions" 
* Formula engine: Ability to reference columns 
* JS API: onBeforeRunAction, onAfterRunAction 
* Ability to return multiple outputs in JS 
* Sequence Translator: using existing pattern and saving it as own 
* Add statistical functions 
* Types.getCommonElementType \- provided a description (still needs a good implementation) 
* Minor changes to the environment setup instructions. Added a couple of the potentially useful scripts. 
* Harmonize file names (data_actions -> functions, data_action -> func) 
* Ddt: harmonize folder structure 
* (Bug) Functions: Unexpected behavior of some statistical functions 
* Removed duplicate import 
* ClinicalCase \- a plugin for analyzing SDTM-based clinical studies data (WIP)
* Modified default config 
* Wiki: Onboarding Content (WIP)
* Grok connect: deserialization fix 
* PepLens: table with sequences application initiated 
* Sequence Translator: overhang options 
* (Bug) Datlas skips row count check for package queries  
* (Bug) Exception on SetAutoComplete 
* (Bug) Deploy failed if public is not checked out 
* (Bug) Query ExtractParameters clears existing params 
* (Bug) View Layouts: When saving layout in an open project, new layout is not generated, but the project layout is overwritten 
* Ability to use other parameters from choices and suggestions statements 
* Better event firing 
* Removed labs submodule 
* Fixed analyzer warnings 
* Fixed submodule host 
* Oligo Batch Calculator: save as CSV 
* Ability to switch to advanced mode in Function View 
* Chem: activity cliffs ported to JS 
* (Bug) Exception when running script as app 
* view polishing 
* Fixed test 
* Grok connect: add batch mode  
* (Bug) Functions: Variable "row" don't work 
* (Bug) "JS functions with multiple output params must return JS object" exception when there are no output parameters 
* Got rid of DataFrameViewer.view in favor of ViewerBase.view 
* Dataframe: version getter for columns 
* Set default admin password for docker-compose  
* SequenceTranslator: fix of Shifts To Align Numbers Inside Circle 
* Sequence Translator: in modification sections enumerate nucleotides only 
* Updated EC2 instructions 
* Updated help index 
* MPLOT: init 
* Minor code harmonization 
* MultiPlot: add Close button "X" and functionality 
* SequenceTranslator: minor changes 
* JS API: Stats.uniqueCount 
* Sequence Translator: align non-overhang nucleotides on the right edge of both strands 
* JS: DataFrame.plot. methods 
* Grok connect: add helpers for dataframe substitution  
* Functions: Parameter categories 
* Box Plot: showCategorySelector and showValueSelector properties 
* (Bug) Box Plot: an exception when stdev(value) = 0 
* (Bug) Box Plot: improve initial choice of value column (stdev > 0 if possible) 
* Box Plot: better border color 
* Viewers: support for styles 
* Viewers: implement 'dashboard' style for all standard viewers 
*  MultiPlot: show/hide controls using ">" icon 
* (Bug) MultiPlot: contorls not removing when plot remove 
* Clinical Case: Timelines View (WIP)
* Clinical Case: Patient Profile View (WIP)
* Removed devops position 
* Timelines: remove the additional visulization in favor of MultiPlot 
* (Bug) Functions View: Click on the selected category does not remove the check mark 
* Updated copyright messages. 
* JS: grok.functions.call should return JS object for multiple output parameters 
* Oligo Batch Calculator: create separate functions for every indicator 
* SequenceTranslator: check for characters/words outside of allowed set, distinguish rows (WIP)
* MultiPlot: clean code 
* (Bug) MultiPlot: fix package.js 
* (Bug) Extra space breaks function annotation 
* (Bug) grok.dapi.projects.where('bad filter').first() returns a Project instance with d == null 
* JS API: Project.open()  // in workspace 
* JS API: Simplify registering functions that return futures/promises 
* Clinical Case: Adverse Events View (WIP)
* Histogram and BarChart: washed-out default colors 
* MultiPlot: remove hardcode and fix some remarks 
* Line Chart: move hit-testing functionality to series renderers 
* (Bug) Line Chart: "stacked bars" mode does not work 
* Introduced Color.textColor 
* Bar Chart: adaptive font size 
* MultiPlot: settings to properties 
* Bar Chart: in case of many categories, automatically zoom in to a reasonable number of categories 
* Bar Chart: make scrolling smoother 
* Edit Column Dialog: validators change 
* JS API: Column.applyFormula() 
* JS API: DataFrame.dialogs methods (WIP)
* Sequence Translator: final UI adjustments (WIP)
* Calculated columns examples 
* ScatterPlot: changed back marker color 
* ScatterPlot: set default marker size to 10 
* (Bug) Layouts: Grid looses event handlers after layout restore 
* DevTools: project snippets 
* Fix UI tests (WIP)
* Updated public token 
* (Bug) MultiPlot: controls overlap if last plot is hidden 
* (Bug) MultiPlot:  exception (empty grid) when only timeLine chart is shown 
* (Bug) Filters: Multi-value filters does not turn off when corresponding checkbox is off 
* (Bug) Filters: multi-value filters have no square indicator on top to toggle category selection 
* (Bug) Filters: There are no icons for sorting and searching on hover for Multi-value filters 
* (Bug) Filters: clicking on "search" should open search field AND focus on it 
* FuncCall.inputs.values() fix 
* Examples for onBeforeRunAction, onAfterRunAction 
* Remove DG.toDart wrapper from string keys 
* Add missing parameter types 
* OligoBatchCalculator: calculator function returns multiple values 
* (Bug) Chem: proper cache invalidation (WIP)
* Fixed analyzer warning 
* fixes typo 
* DevTools: convert to TypeScript 
* FunctionsWidget 
* implemented several validation rules, added validation view, added error summary to study summary view 
* MultiPlot: make tooltips as in DG (WIP)
* Update public token 
* (Bug) JsViewer root is not correctly wrapped 
* (Bug) JS API: properties cannot be changed in JsViewer (JsViewer.props and JsViewer.setOptions result in errors) 
* (Bug) JS API: JsViewer is no longer in the context of the onContextMenu event 
* (Bug) JS API: JsViewerHostCore is returned instead of Viewer instance 
* grok 9125 prepared for biowasm integration 
* added validatin for DM domain, implemented filtering of validatin result based on selected rule and ighlighting of cell with violated rule 
* added check of dm dataframe for null 
* removed extracting domain from sdtm-tag 
* JS API: correct viewer type names 
* DevTools: generate code for viewers customization (WIP)
* Oligo Batch Calculator: check for characters/words outside of allowed set, distinguish rows (WIP)


# 2021-06-24 Dev build 0.92.14

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.14`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) Share button doesn't work for projects 


# 2021-06-24 Dev build 0.92.13

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.13`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) JS API: JsViewerLoader.instance sometimes returns null instead of jsViewer  
* (Bug) Notebook sharing doesn't work 


# 2021-06-24 Dev build 0.92.12

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.12`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Fixed project deploy 
* Connections: add system connections  
* (Bug) Functions: function search won't work on packages 
* (Bug) Usage analysis: persistent loader 
* DevTools: add the `Snippets` pane (WIP)
* JS API: format(x, fmt) 
* Power Pack: Widgets: Kpi 
* Math Functions: Round10() 
* Math Functions: RandBetween() 
* Power Search: integrate with widgets 
* Power Search (WIP)
* JS API: Added grok.shell.startUri 
* JS API: TabPane: new properties: root, content, parent 
* Power Pack: Search: routing 
* (Bug) An exception when closing the last view 
* Not opening a default view at startup in debug mode (in favor of PowerPack's home view) 
* Not opening a default view at startup when the URL starts with '/search' (in favor of PowerPack's search) 
* Math Functions: Avg() 
* JS API: PropertyBag.setAll(object) 
* Power Search: community-curated, template-based, widget-driven search engine (WIP)
* DG.tools.createElementFromHtml 
* Power Pack: HtmlWidget 
* Introduced Property.editor (used to be Property.info.editor) 
* Fixed a typo. 
* JS API: Document function roles 
* Wiki: additional sorting examples 
* Wiki: how-to pages title consistency 
* Math Functions: Median() 
* Update ui examples 
* Oligo Batch Calculator: calculate optical density(OD), nmole and mass from Yield Amount & Units 
* Add color blindess-friendly color palettes  (WIP)
* Text Functions: StrFind() 
* Text Functions: StrLeft() 
* (Bug) JS API: code editor view differs from the original one (WIP)
* Text Functions: StrRight() 
* Text Functions: StrRepeat() 
* Text Functions: RegExpReplace() 
* Text Functions: RegExpExtract() 
* Create color-checker.js 
* Layout templates for JS API Examples 
* SequenceTranslator: adjustments to svg creation 
* Updated help, section transform \- text functions 
* (Bug) Project | Links: Project link is generated incorrectly 
* JS API: Viewer.onContextMenu event 
* Added reuseList function 
* Chem Activity cliffs \- visulisation polishing 
* Added solution for git push to public repo issue 
* Removed an annotation comment from samples opened in JS editor 
* (Bug) Grok Script: Syntax error: "abc".trim().toUpperCase()  
* Wiki: details on color-coding 
* (Bug) Deploy ignores friendlyName in queries 
* Deploy: test connections migration (WIP)
* Fixed typo 
* Core: Add key missing calculation functions (WIP)
* Wiki: expand an article on adding new columns (WIP)
* JS API: ability to read metadata for a package or a group of packages 
* PhyloTreeViewer: JS-based phylogenetic tree visualisation (WIP)
* (Bug) Share button doesn't work for projects 
* Line-height fix 


# 2021-06-18 Dev build 0.92.11

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.11`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) Grid: empty space on the right 
* (Bug) Filters: filter component state is different across pages 
* Update ui.md 
* SequenceTranslator: svg color adjustments 
* Oligo Batch Calculator: create initial version and test on files 
* Create datagrok installation instructions (v.2) 
* Correct Markdown view of instructions 
* Deploy: test connections migration (WIP)
* Connections: add system connections  
* Correct instruction points 
* Binning Functions: BinBySpecificLimits 
* DateTime extractor functions (year, month, day, etc) 
* (Bug) Connectors: Impala: int32max instead of real values 
* bug with table view fixed 
* JS API: ability to read metadata for a package or a group of packages 
* Renamed "UsageAnalysisUI" to "UsageAnalysis" 
* AutoGen code rebuild (pub get and build execution) 
* Renamed PowerPack to meet both npm and DG conventions 
* Simplified datagrok installation instructions 
* [(Bug) JS API: df.replace won't work on a column with the same name](https://community.datagrok.ai/t/annot-replace-a-column-with-a-column-with-the-same-name/528) 
* (Bug) Connectors: output table files for debug 
* JS API: Advanced Func API (WIP)
* Sequence Translator: add modifications to beginning of both strands 
* Cast exception to string 
* (Bug) Filters: range slider filters out nulls 
* (Bug) OnDialogClosed event fires twice 
* (Bug) File Sharing doesn't work 
* Created index 
* *BREAKING\* Log function renamed to Ln 
* (Bug) JS Editor: exception when ApiSamples package is not there 
* Adjusted demo query name to 'System:CoffeeCompany:StoresInState' 
* Math Functions: Log(arg1, arg2) 
* Math Functions: Log10() 
* Math Functions: Round10() 
* Correct mistake in dart installation instructions 
* JS API: add simpleMode and move presentationMode, hideTabsInPresentationMode to windows 
* SequenceTranslator: refactoring of svg creation 
* Grok connect: add CSV chars escaping  
* JS API: improve the HttpDataSource.include method to accept entity properties 
* DevTools: add the `Snippets` pane (WIP)
* Math Functions: RandBetween() 
* (Bug) Functions: function search won't work on packages (WIP)
* Updated public token 


# 2021-06-15 Dev build 0.92.10

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.10`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Update generate-html-help.side 
* Handling the case when an existing column is being replaced 
* Deploy: test connections migration (WIP)
* ClinicalCase \- a plugin for analyzing SDTM-based clinical studies data (WIP)
* Activity cliffs visual modification 
* Update beta users table 
* Sequence Translator | AXOLABS: Add tooltips for field headers which do not fit completely 
* Moved to table view 
* Updated public token 
* Reduce DataQueryCall tag payload (WIP)
* user friendly options added 
* (Bug) Miniconda hangs CVM builds (WIP)
* About-widget 
* Power Search: deep integration with functions (WIP)
* Add .vscode to gitignore 
* Power Search: JS integration 
* JS API: fluent API for plotting (WIP)
* (Bug) Deploy: projects migration bugs (WIP)
* Connections: add system connections  
* SequenceTranslator: flexible Axolabs sequence conversion 
* (Bug) Impossible to add a tile viewer to the table  
* Sequence Translator: detect font color suitable for background color 
* (Bug) Bar Chart: stacked mode: labels are rendered behind the bars 
* (Bug) Log export hangs 


# 2021-06-14 Dev build 0.92.9

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.9`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Sequence Translator | AXOLABS: Add tooltips for field headers which do not fit completely (WIP)
* Fixed int to double type mismatch 
* Timelines: add a heading visualizing data from LB domain (WIP)
* Added sudo prompt in datagrok-tools installation 
* Update help index 
* (Bug) Miniconda hangs CVM builds (WIP)
* Update ui.md 
* Create markdown.js 
* (Bug) Connectors: Oracle: NullPointerException in DB table content  
* Packages: Curation, Beta (WIP)
* Meta for ServerSettings 
* Deploy: test connections migration (WIP)
* (Bug) IntColumn.fromList(values) does not work with values outside of the int32 range 
* (Bug) Query runs forever 
* (Bug) JS API: Viewer.close() sometimes fails (WIP)
* Updated unit tests for integers 
* Extended meta package: an example 
* Wiki: link & typos correction 
* Typo correction 
* Code cleanup 
* Power Search: deep integration with functions (WIP)
* Fixed a test 
* Updated public token 


# 2021-06-10 Dev build 0.92.8

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.8`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) Func: not handled exception 
* Grok connect: tests (WIP)
* UsageAnalysis harmonization (WIP)
* Reordering ui samples files 
* JS wrap for activity cliffs script 
* JS API: add ui.markdown 
* Deploy: test connections migration (WIP)
* SequenceTranslator: flexible Axolabs sequence conversion (WIP)
* Power Search: community-curated, template-based, widget-driven search engine (WIP)
* (Bug) Current user is system while deploying 
* (Bug) Scatter Plot: axes should not pick up automatically assigned precision-related column format (WIP)
* (Bug) Miniconda hangs CVM builds (WIP)
* Conditional color-coding: bring options from the property panel to column context menu 
* (Bug) Hide function deselects the selected columns from the tooltip  
* (Bug) Bar Chart: Viewer coloring settings should take precedence over the grid coloring settings (WIP)
* PLOT3D: transparecy, some props, fixed filter 
* (Bug) Sequence Translator | Main: If you enter the wrong input value, then many balloons with stacktraces will appear (WIP)
* (Bug) Sequence Translator | AXOLABS: Name of an existing pattern is not displayed (WIP)
* Packages: Curation, Beta (WIP)
* Ability to drop DataQueryCache (WIP)
* (Bug) Data | Unpivot does not work 
* (Bug) Column format changes are not persisted in layouts 
* Sequence Translator | AXOLABS: Add tooltips for field headers which do not fit completely (WIP)
* BoxPlot: removed "Reset View" menu option 
* .multi-value-separated: added dot in front 
* (Bug) Deploy: deploy will fail if a project doesn't contain "tags"  
* NotificationsWidget 
* Widget-based Welcome page (WIP)
* Power Pack: Recent Projects Widget 
* FileSystemWidget 
* (Bug) Applications: Table views open both in the "Tables" section on side bar and in the open application section 
* PLOT3D: added columns size and column color props and some clean code 
* PLOT3D: fix animation 
* IntRle encoder: refusing to encode when max-min > 2^32 
* DataQuery View: a way to easily share this query 
* Timelines: toggle zoom sliders 


# 2021-06-07 Dev build 0.92.7

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.7`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* SequenceTranslator: flexible Axolabs sequence conversion (WIP)
* Function View (WIP)
* PLOT3D: add functions and some refactor 
* Power Pack: Power Search: integrate with functions (WIP)
* UTC support in the datetime parser (WIP)
* Package content validation: up-to-date bundle check 
* Grok connect: tests (WIP)
* PLOT3D: before refactor 
* PLOT3D: fix problem with getGeometry() and other 
* PLOT3D: delete unused files 
* Fixed exception on deploy 
* (Bug) Current user is system while deploying 
* (Bug) UsageAnalysis hangs backend (WIP)
* JS API: dialog.clear() 
* Add string name indexing for columns 
* GrokML: integrate WASM calls with Dart (pearsonsR for the correlation plot) (WIP)
* WIP 
* PLOT3D: refactor 
* JS API: View.createByName('functions', {search: '#math'}): a way to create standard Dart views 
* Power Search: integrate with default views ("functions", "databases", "files", etc) 
* Fixed a typo 
* (Bug) Can't close QueryView 
* Charts: Timelines (WIP)
* Update icons.js 
* Package content validation: ensure referenced files are in place (WIP)
* Widget-based Welcome page (WIP)
* Use current datagrok-api version for remote package build 
* Updated public token 
* Code cleanup 
* Fixed type casting 
* Update ui.md 
* Connections: add system connections  (WIP)
* DB: improve and fix migration scripts  (WIP)
* (Bug) Grok Script Tests: show errors and check for errors (WIP)
* PLOT3D: working original mesh and shaders to WegGL2 
* Widgets: remember widget state in user settings (WIP)
* PLOT3D: working filter and fixed inverted normals 
* Power Search: community-curated, template-based, widget-driven search engine (WIP)
* PLOT3D: auto rotation 
* (Bug) Miniconda hangs CVM builds (WIP)
* Calculation of actvity cliffs 
* Update ui.md and add new ui samples 
* Wiki: Onboarding Content (WIP)
* Fixed indentation 


# 2021-06-01 Dev build 0.92.6

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.6`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* chem curation description added 
* SequenceTranslator: flexible Axolabs sequence conversion (WIP)
* (Bug) JS API: Viewer.props does not work 
* ScatterPlot custom renderer: disable overlay dot rendering on categories hover  
* (Bug) Pie Chart: unnecessary datetime aggregation 
* (Bug) WebViewer won't load as part of the project (WIP)
* Demo of structural curation added 
* (Bug) JS API: toJs won't work on GrokPackage 
* Minor style improvements 
* Fixed paths to images 
* JS API: Func.find(package, name, tags, returnType) 
* JS API: Func.appy(params): Promise<TResult> 
* Widget-based Welcome page (WIP)
* Added documentation for FASTA reader 
* Added tag for sequence field 
* Update README.md 
* Wiki: How to manipulate viewers (WIP)
* JsWidgetBase: getJsObject now returns jsWidget 
* JS API:  add a possibility to use custom ScatterPlot renderer (WIP)
* fetchProxy: better typing. 
* JS API: ui.link(text, url | action, tooltip?) 
* Power Pack: Community Widget 
* JS API: ui.iframe(src) 
* Power Pack: WebWidget 
* PLOT3D init 
* JS API: toDart(): ability for JS classes to define custom conversion 
* JS API: ui.icons for commonly used icons 
* JS API: ui.tools.setHoverVisibility(host, elements) 
* JS API: ui.image 
* Power Pack: Learn Widget (WIP)
* Scripting: interactivity: show execution results in the property panel 
* Run script: added F5 shortcut to the tooltip 
* Global Search (WIP)
* Power Pack: Search: FunctionSearchProvider 
* Power Pack: Search: WikiSearch 
* Charts: Tree Viewer (WIP)
* info page polishing 
* chem Cureate: info page polishing 
* chem Curate: jpg files changed to png 
* Help: Ordering in chapters 
* PLOT3D: import three.js and show scene 
* Added beta flag to the old Usage analysis 
* Function View (WIP)
* (Bug) JS API: JsViewerHostCore is returned instead of Viewer instance 
* Data: activity_cliffs, chem_standards 
* (Bug) Miniconda hangs CVM builds (WIP)
* Packages: Curation, Beta (WIP)
* Property Panel expanding indicator 
* UI: Application pictures and Views grouping 
* Use current datagrok-api version for remote package build 
* Power Pack: ServicesWidget 
* JS API: ui.cards: an easy way to build cards (WIP)
* PLOT3D: add look() function 
* CVM: Octave to Matlab compliance (WIP)


# 2021-05-27 Dev build 0.92.5

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.5`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Help: Ordering in chapters 
* DevTools: add the `Snippets` pane (WIP)
* Dialog should have ui-panel class 
* Javascript vectorization 
* SequenceTranslator: Conversion to Axolabs sense/antisense strands 
* (Bug) Radiobutton throws exception when created 
* DeleteConnectionCache harmonization 
* Code cleanup 
* Function View (WIP)
* PowerPack: Compare Columns CSS fixes 
* Fix UI tests (WIP)
* Selenium: Add ui-test for Chord viewer 
* SequenceTranslator: flexible Axolabs sequence conversion (WIP)
* JS API: add scatterPlot.worldToScreen() method 
* PowerPack: Compare Columns: Add Column button 
* JS API: add events to the ColumnComboBox 
* Fixed analyzer warning 


# 2021-05-24 Dev build 0.92.4

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.4`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* PivotViewer 
* (Bug) Miniconda hangs CVM builds (WIP)
* DevTools: add the `Snippets` pane (WIP)
* PowerPack: Compare Columns fixes 


# 2021-05-22 Dev build 0.92.3

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.3`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) Miniconda hangs CVM builds (WIP)
* Dataframe: Detect column max significant digits in CSV loading (WIP)
* PowerPack: rewrite Compare Columns to TypeScript 
* Grid: automatically pick up cell type from data frame 
* JS API: Stats.histogramsByCategories 
* Statistics: Multi-histograms (WIP)
* Made Widget.root non-nullable. WIP. 
* ColumnsCorrelationMeta 
* Correlation plot: on cell click, show the corresponding scatter plot in the property panel 
* Histogram: multiple histograms mode (WIP)
* Dialog should have ui-panel class 


# 2021-05-21 Dev build 0.92.2

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.2`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* PowerPack: rewrite Compare Columns to TypeScript 
* (Bug) Log export hands 
* Help: Ordering in chapters (WIP)
* Grid: default column width doesn't respect format 
* (Bug) Grok connect: java.lang.NullPointerException with pattern parameter 
* Dataframe: Detect column max significant digits in CSV loading (WIP)


# 2021-05-21 Dev build 0.92.1

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.1`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* fixes regression in 2D layout alignment of unknown origin 
* Chem package version 0.9.8, bump. 
* BioSignals: make initial chart look similar to content of accordion panes 
* DINQ: Support JSON in WhereSmart method 
* DevTools: add the `Snippets` pane (WIP)
* JS API: context is sometimes missing in onAccordionConstructed() 
* JS API: add Grid.getVisibleCells 
* Datagrok API and system samples, add cross linking to improve connectivity (WIP)
* BioSignals: Add url routing for Physionet records 
* (Bug) JS Samples: to-hierarchy.js example does not work (WIP)
* NGL viewer (WIP)
* Scripting: Ability to set custom name for parameter 
* Scripting: Ability to set postfix for parameter 
* Scripting: Ability to set postfix for parameter 
* UI: Ability to set dialog background (WIP)
* Fixed TS warning 
* JS: Ability to hide table params from FuncCall editor 
* Updated public token 
* (Bug) Core: deepEquals method on Prop class always returns false 
* Running an ulmo sample in PythonScripts via pip+conda env. 
* Function editors: open editor with DataFrame parameters even when no tables are open (WIP)
* Function editors: allow to load input data frames after the dialog is open (WIP)
* BioSignals: publish connection to Physionet demo files 
* Fixed warnings 
* SequenceTranslator: MM12 Sequence Conversion 
* JS API: add grok.shell.views and grok.shell.tableViews 
* SequenceTranslator: OP100 Sequence Conversion 
* Function View (WIP)
* Power Pack: Compare Columns 
* Fixed logins in some connections 
* RepertoireBrowser: MSA tabs (WIP)
* Add files via upload 
* Update chem_curator.py 
* Ability to drop DataConnectionCache (WIP)
* UsageAnalysis: Added 'users by date' query 
* (Bug) Molecules not rendered for columns returned by scripts 
* datagrok-tools: code clean-up 
* (Bug) DINQ: Permission check on every object type 
* Fixed analyzer warning 
* Context commands: "Add derived column..." for multiple columns 
* Viewer.serialize(includeDefaults): an option to include default values 
* update "Common Actions" section 
* (Bug) Events: onViewerAdded and onViewerClosed are sent twice 
* JS API: add square brackets for columns 
* Prevent DG disk space usage inside docker container 
* (Bug) Invalid session message when running second node of DG 
* (Bug) Exceptions are not handled during DG start 
* Increase Func socket timeout 
* Handle AWS "SlowDown" exception 
* JS API:  add a possibility to use custom ScatterPlot renderer (WIP)
* ClinicalCase: harmonize menu items 
* PackageView: added a margin 
* ClinicalCase: updated description 
* JS API: add Accordion.removePane() 
* (Bug) Filters: if multi-value filters are present in the panel, the reset button doesn't work 
* Updated package description. 
* SequenceTranslator: ABI Sequence Conversion 
* JS API: ui.info for a yellow info bar (WIP)
* ClinicalCase \- a plugin for analyzing SDTM-based clinical studies data (WIP)
* Update ui.ts 
* UI: Copyright 2019 -> 2021. 
* (Bug) Chem: Tile Viewer issue (WIP)
* GrokML: integrate WASM calls with Dart (pearsonsR for the correlation plot) (WIP)
* (Bug) JS API: ColumnListProxy:  Cannot convert a Symbol value to a number 
* Dialog should have ui-panel class 
* Dataframe: Detect column max significant digits in CSV loading (WIP)
* Core: Add GrokJsObjectConst 
* (Bug) Miniconda hangs CVM builds (WIP)
* JS API: ScatterPlot: add viewport 


# 2021-05-10 Stable version 0.92.0

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.92.0`
  *  `docker pull datagrok/datagrok:stable`
* CVM: 
  *  `docker pull datagrok/cvm:0.92.0`
  *  `docker pull datagrok/cvm:stable`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) Chem: R-Groups Analysis with client RDKit 
* Ability to share files 
* StringColumn: batch editing: copying existing values 
* Chem: renderer selection checks. 
* BioSignals: add possibility to change previously chosen filter script 
* sequenceTranslator: use Ternary Operator when it brings clarity 
* BioSignals: choose input signal from context(WIP) 
* BioSignals: choose input signal from context 
* Proper dependencies 
* (Bug) Incorrect QNum parsing 
* do not crash on an empty substructure library 
* Chem: Robust invalid molecules handling in parallel Substructure Search 
* PhyloTreeViewer: JS-based phylogenetic tree visualisation (WIP)
* Bump commons-io from 2.6 to 2.7 in /connectors/grok_connect 
* Biosignals: add opened tables to App context 
* Connectors: Redshift: Schema browsing (WIP)
* Datagrok API and system samples, add cross linking to improve connectivity (WIP)
* Trellis Plot: ability to enlarge individual in-trellis viewers 
* NGL viewer (WIP)
* (Bug) JS API: Columns.byTags does not work 
* Widgets: MultiValueFilter 
* removed to LABS 
* removed src 
* files moved to LAB 
* Join tables: Don't include key fields by default 
* Icon Tooll 
* BioSignals: added app description, and ui.buttonsInput 
* BioSignals: rename folders with scripts to 'filters', 'extractors', 'indicators' 
* (Bug) Queries: "Converting object to an encodable object failed" with dataframe as parameter 
* Multi-value filters 
* (Bug) Excel import: an empty column is created as part of the dataframe 
* BarChart categories error 
* (Bug) Layouts: Undocked viewer remained after "Close All" 
* update package icons 
* update icon 
* JS API: BitSet: getBuffer() 
* MSA tabs added 
* BioSignals | AnnotatorViewer: generate timestamps from sampling frequency 
* Fixed grok_connect imports 
* update icons for script packages 
* change package name 
* RepertoireBrowser: msa rendering changed 
* RepertoireBrowser: added debouncing 
* Updated public toke 
* BioSignals: add possibility to enter sampling frequency for local files 
* Update package.png 
* (Bug) Grid: html cells are not re-rendered after sorting a column 
* Connectors: add default schema to providers and use as condition for schema browsing  
* Prevent socket memory consuming 
* RepertoireBrowser: gff annotations 
* (Bug) BioSignals: rename input variables of scripts to camelCase to show them properly in the UI (split into words) 
* Adjusted a unit test to reflect recent changes in the 'join' implementation 
* Added demo files 
* RepertoireBrowser: gff annotations (slight changes) 
* BioSignals: option to delete accordion pane 
* (Bug) Grok connect: can't connect to MsSQL using TLS10 
* JBio;Repertoir: Added GFF annotations for protein translation, tuned debounce timeout for better perf 
* Repertoir: code cleanup 
* Table View: Columns pane: add search 
* BioSignals: Preserve "Filtering, Extraction and Indicator" parameters while selecting "Physionet record" or "Column" 
* BioSignals: fixes and optimizations 
* BioSignals: AnnotatorViewer: round values and change tooltip formatter 
* JS API: ui.info for a yellow info bar (WIP)
* Docs: how to run Chrome with CORS disabled. 
* Docs: info about Labs repo and importDemoPackages setting for dev.datagrok.ai. 
* (Bug) Chem: several molecule filters won't work together (WIP)
* (Bug) Core: default tooltip config is no longer saved with layout 
* (Bug) JS Viewer: Incorrect height calculation 
* (Bug) Bar Chart: coloring cannot be disabled 
* BioSignals: fixes and minor changes 
* sequenceTranslator: fix of line feed bug; notify if sequence length is too short 
* (Bug) Filters: adding returns a filter for the column added previously instead of the currently selected one 
* (Bug) Bar Chart: coloring gets only applied when editing colors in a grid column 
* Packages: beta flag 
* Updated public token 
* Removed debug line 
* Chord Diagram (WIP)
* (Bug) Chem: Tile Viewer issue (WIP)
* sequenceTranslator: use functions to display default translations 
* datagrok-tools version 4.1.2, bump 
* (Bug) sequenceTranslator: molecule SVGs append after adding/removing characters to input sequence 
* (Bug) UI: new view appears below the current 
* Fixed TS warnings 
* ClinicalCase: minor adjustments 
* SequenceTranslator: OP100 Sequence Conversion (WIP)
* SequenceTranslator: MM12 Sequence Conversion (WIP)
* RepertoireBrowser: MSA tabs (WIP)
* Wiki: Harmonization (WIP)
* Dataframe: Detect column max significant digits in CSV loading (WIP)
* (Bug) Api Samples: "db-master-details-manual" example does not work 
* RepertoireBrowser: code formatting 
* Data: fixed demog redundant precision 
* JS API: add Column.isVirtual 
* (Bug) Chord Viewer: the diagram gets rendered twice when the same column is used for source and target 
* PhyloTreeViewer: node/branch tooltips based on dataframe rows (WIP)
* JS API: DG.TAGS.FORMAT tag as a preferred way to set the column format 
* JS API: add event onAccordionConstructed and Accordion.context getter  
* Minor cleanup 
* (Bug) BioSignals: error after calling Python filter first 
* Menu.fromCommand: show command description in the tooltip 
* (Bug) Viewers: the menu item `Viewer` is not visible in uploaded projects 
* Menu capitalization harmonization 
* Crash when adding columns to an existing DF/Grid 
* Grid: ShowVisibleColumnsInTooltip property 
* (Bug) Grid: drag-and-drop column reordering: provide drop zones 
* JS API: tags must only be strings 
* BioSignals: SequenceTranslator: use ui.info() for yellow information bars in views 
* Biosignals: parallel zoom for lineCharts (WIP)


# 2021-05-05 Dev build 0.91.10

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.91.10`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Chord Diagram (WIP)
* datagrok-tools version 4.1.2, bump 
* (Bug) sequenceTranslator: molecule SVGs append after adding/removing characters to input sequence 
* (Bug) UI: new view appears below the current 
* Fixed TS warnings 
* ClinicalCase: minor adjustments 
* sequenceTranslator: OP100 Sequence Conversion (WIP)
* Updated public token 
* sequenceTranslator: MM12 Sequence Conversion (WIP)


# 2021-05-04 Dev build 0.91.9

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.91.9`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* BarChart categories error 


# 2021-05-04 Dev build 0.91.8

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.91.8`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* PhyloTreeViewer: JS-based phylogenetic tree visualisation (WIP)
* Packages: beta flag (WIP)
* Updated public token 
* Removed debug line 
* Chord Diagram (WIP)
* (Bug) Chem: Tile Viewer issue (WIP)
* sequenceTranslator: use functions to display default translations 


# 2021-05-04 Dev build 0.91.7

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.91.7`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) Core: default tooltip config is no longer saved with layout (WIP)
* (Bug) JS Viewer: Incorrect height calculation 
* (Bug) Bar Chart: coloring cannot be disabled 
* BioSignals: fixes and minor changes 
* PhyloTreeViewer: JS-based phylogenetic tree visualisation (WIP)
* Datagrok API and system samples, add cross linking to improve connectivity (WIP)
* sequenceTranslator: fix of line feed bug; notify if sequence length is too short 
* (Bug) Filters: adding returns a filter for the column added previously instead of the currently selected one 
* (Bug) Bar Chart: coloring gets only applied when editing colors in a grid column 


# 2021-05-03 Dev build 0.91.6

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.91.6`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* BioSignals: option to delete accordion pane 
* (Bug) Grok connect: can't connect to MsSQL using TLS10 (WIP)
* JBio;Repertoir: Added GFF annotations for protein translation, tuned debounce timeout for better perf 
* Repertoir: code cleanup 
* Table View: Columns pane: add search 
* BioSignals: Preserve "Filtering, Extraction and Indicator" parameters while selecting "Physionet record" or "Column" 
* PhyloTreeViewer: JS-based phylogenetic tree visualisation (WIP)
* BioSignals: fixes and optimizations 
* DataGrok API and system samples, add cross linking to improve connectivity (WIP)
* BioSignals: AnnotatorViewer: round values and change tooltip formatter 
* JS API: ui.info for a yellow info bar 
* Docs: how to run Chrome with CORS disabled. 
* Docs: info about Labs repo and importDemoPackages setting for dev.datagrok.ai. 
* (Bug) Chem: several molecule filters won't work together (WIP)


# 2021-04-30 Dev build 0.91.5

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.91.5`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* RepertoireBrowser: gff annotations 
* BioSignals: rename input variables of scripts to camelCase to show them properly in the UI (split into words) 
* Adjusted a unit test to reflect recent changes in the 'join' implementation 
* Added demo files 
* JS API: BitSet: getBuffer() 
* RepertoireBrowser: gff annotations (slight changes) 
* BioSignals: option to delete accordion pane (WIP)
* DataGrok API and system samples, add cross linking to improve connectivity (WIP)
* (Bug) Queries: "Converting object to an encodable object failed" with dataframe as parameter 


# 2021-04-29 Dev build 0.91.4

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.91.4`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Update package.png 
* (Bug) Grid: html cells are not re-rendered after sorting a column 
* Connectors: add default schema to providers and use as condition for schema browsing  (WIP)
* Prevent socket memory consuming (WIP)


# 2021-04-29 Dev build 0.91.3

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.91.3`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* BioSignals: choose input signal from context 
* do not crash on an empty substructure library 
* Chem: Robust invalid molecules handling in parallel Substructure Search 
* PhyloTreeViewer: JS-based phylogenetic tree visualisation (WIP)
* Biosignals: add opened tables to App context 
* Connectors: Redshift: Schema browsing (WIP)
* DataGrok API and system samples, add cross linking to improve connectivity (WIP)
* Trellis Plot: ability to enlarge individual in-trellis viewers 
* NGL viewer (WIP)
* (Bug) JS API: Columns.byTags does not work 
* Widgets: MultiValueFilter (WIP)
* removed to LABS 
* removed src 
* files moved to LAB 
* Join tables: Don't include key fields by default 
* Icon Tooll 
* BioSignals: added app description, and ui.buttonsInput 
* BioSignals: rename folders with scripts to 'filters', 'extractors', 'indicators' 
* (Bug) Queries: "Converting object to an encodable object failed" with dataframe as parameter 
* Multi-value filters 
* (Bug) Excel import: an empty column is created as part of the dataframe 
* BarChart categories error (WIP)
* (Bug) Layouts: Undocked viewer remained after "Close All" 
* update package icons 
* update icon 
* JS API: BitSet: getBuffer() 
* MSA tabs added 
* BioSignals | AnnotatorViewer: generate timestamps from sampling frequency 
* Fixed grok_connect imports 
* update icons for script packages 
* change package name 
* RepertoireBrowser: msa rendering changed 
* RepertoireBrowser: added debouncing 
* Updated public toke 
* BioSignals: add possibility to enter sampling frequency for local files 


# 2021-04-26 Dev build 0.91.2

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.91.2`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Chem: renderer selection checks. 
* BioSignals: add possibility to change previously chosen filter script 
* sequenceTranslator: use Ternary Operator when it brings clarity 
* BioSignals: choose input signal from context(WIP) 
* BioSignals: choose input signal from context (WIP)
* Proper dependencies 
* Ability to share files 
* (Bug) Incorrect QNum parsing 


# 2021-04-23 Dev build 0.91.1

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.91.1`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) Chem: R-Groups Analysis with client RDKit 
* Ability to share files 
* StringColumn: batch editing: copying existing values 


# 2021-04-23 Stable version 0.91.0

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.91.0`
  *  `docker pull datagrok/datagrok:stable`
* CVM: 
  *  `docker pull datagrok/cvm:0.91.0`
  *  `docker pull datagrok/cvm:stable`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues



# 2021-04-23 Stable version 0.90.0

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.90.0`
  *  `docker pull datagrok/datagrok:stable`
* CVM: 
  *  `docker pull datagrok/cvm:0.90.0`
  *  `docker pull datagrok/cvm:stable`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* JS snippets: Viewer.property and Viewer.getInfo() 
* JS API: bitset.handleClick: a helper function that selects rows by passing a row predicate and a mouse event 
* Globe Viewer 
* Doing the right thing when get_mol throws; more optimal LRU cache use 
* RDKit-based searches 
* Viewer.grid(table) now returns an instance of Grid instead of Viewer 
* Replaced `onFrameAttached` method 
* Word Cloud (WIP)
* UI tests updates 
* Add UI test for new SPE dialog 
* Create UI test for "Open text" feature 
* Wiki: How to Build Custom Viewers 
* Grid: RDKit-based structure rendering (WIP)
* Chem package version 0.4.4, bump 
* Chem RDKitDemoPackage -> ChemPackage 
* Chord Diagram (WIP)
* datagrok-tools: URL suggestion fix 
* Chem Benchmark (WIP)
* Wiki: Harmonization (WIP)
* JS API: add possibility to filter columns  
* Added filter-columns.js snippet 
* Chem, RDKit-based (WIP)
* Fixed a couple of links 
* Wiki: Data Access 
* Normalize formatting: WebStorm JS default with 2 spaces 
* ... 
* Viewers: Sankey Diagram (WIP)
*  Get groups after groupBy 
* Added GroupByBuilder.getGroups() snippet 
* Added CONTRIB.md 
* Fixed the Copyright date 
* JS API: add Pearson and Spearman correlations  
* JS API: Dialog improvements 
* Internalize timeAsync to DG.timeAsync in DG js-api 
* VDJ tools (WIP)
* JS API: files 
* File editors 
* ui.onSizeChanged (+Dartium polyfill) 
* onSizeChanged dmeo 
* File editors: NGL Viewer 
* File editors: PDF viewer 
* (Bug) Table view slows down(or even freezes) selection 
* JS API: Row Matcher: ability to select or filter rows using structured queries 
* Rearranged the samples. 
* JS API: ability to provide 'where' clause for aggregation 
* (Bug) DateTime aggregation bug 
* JS API: Detector: a helper class for the efficient sampling of column content 
* NGL Viewer: added "PDB Viewer" widget 
* Remove ProteinViewer from the core in favor of the NglViewer package 
* Wiki: custom file viewers 
* NglViewer: more documentation. 
* Wiki improvements 
* JS API: ability to enumerate viewer properties (WIP)
* Wiki: Semantic Type Detectors 
* Expand the OpenApi article (WIP)
* (Bug) ui.divV(null) throws an exception 
* (Bug) Widgets: RadioButtonFilter does not work 
* Improved a sample code snippet. 
* Update query-transformations.md 
* DataFrame: virtual columns: support built-in types 
* Fixed links 
* fetchProxy method 
* Updated an example 
* Wiki: Computations 
* Added missing semicolons 
* Charts package (WIP)
* Support renv 
* Updated onSizeChanged method 
* Versioned RDKit + wasm caching 
* JS API: DataFrame.getSortedOrder 
* Charts: Tree Viewer (WIP)
* (Bug) Dialog.showModal doesn't work 
* JS API: custom viewers: a fluent way to define property properties 
* (Bug) Charts: Tree Map (WIP)
* Harmonize connection parameters 
* Wiki: Add an article about routing (WIP)
* Wiki: Embed viewer 
* Wiki: Table export 
* (Bug) User storage: saves, but doesn't load 
* UI test for new table tooltip 
* Wiki adjustments. 
* Updated connection parameters naming 
* JS API: fetchProxy example 
* Minor fixes of API examples 
* Markdown fix 
* fix error in routing.md 
* Fixed typo, and simplified the language 
* JS API: add TablesClient.uploadDataFrame() and TablesClient.getTable() 
* Added TablesClient.uploadDataFrame() and getTable() snippet 
* Adjust RDKit default molecule look 
* SDTM metadata \- WIP 
* RDKit WASM load on first use or explicitly 
* Updated webpack.config 
* Dumbbell Viewer 
* JS API: add getDensity method for dataframe  
* Commonly used Datagrok-specific libraries (that are supposed to be built into the packages that use them). 
* Refine searches and scoring API; Server + Client RDKit 
* RadioButtonFilter: Limited number of categories to 20 
* Bringing Widget closer to Viewer \- WIP. Viewer.onTableAttached now gets JS DataFrame 
* WidgetMeta (properties, events, functions) (WIP)
* Chem Benchmark (WIP)
* Refine searches and scoring API; Server + Client RDKit (WIP)
* Wiki: shortened the main JS API article in favor of dedicated how-to sections 
* Wiki: Custom Views How-To 
* Parallelize searchSubstructure 
* Add .vscode settings to grok create (WIP)
* (Bug) Structure search results column: substructure column checkboxes can be checked or unchecked 
* UI test for Projects with Data Sync from files 
* Adjusting the Chem samples 
* findSimilarServer won't work (WIP)
* Harmonize ui framework (WIP)
* Ability to get statistics for column subset: Stats.fromColumn(column, mask) (WIP)
* Added TAGS.ID 
* Cache molecules renders 
* Clean coordinate data by a render-through-smiles tag 
* Fixed a comparison 
* Wiki: Grid Customizations 
* Converted to ui.onSizeChanged 
* Render an envelope instead of an empty cell for failed molecules 
* Programming exercises (WIP)
* Update exercises.md 
* Updated fetchProxy examples 
* Simplified package `test` function 
* Fixed templates in datagrok-tools 
* Updated exercises.md 
* User Data Storage overview 
* Update user-data-storage.md 
* Fixed links in user-data-storage.md 
* Wiki: applications \- general, building, debugging. wip 
* Wiki: replaced outdated links 
* Wiki: how-to build an app. wip 
* Wiki: debugging packages. 
* JS API: added accordion persistence key 
* Simplified phrasing and linked to the main article 
* Better phrasing 
* DataFrame: JSDocs for append and unpivot. 
* JS API: DataFrame.append: ability to specify columns and rows to append (WIP)
* Wiki: how-to build an app. WIP: dataframe, credentials 
* Added new package ChemblBrowser 
* Wiki: Project types description 
* Wiki: switch to recommended package naming 
* draft of parametrized query by molregno setted up 
* Wiki: URL fix 
* Wiki: how-to build an app. WIP: restructuring 
* Wiki: how-to build an app. WIP: Semantic annotation and metadata 
* Wiki: how-to build an app. WIP: Dataframe aggregations 
* JS API: Column.init(value | indexToValueFunction) 
* Wiki: how-to build an app. WIP: User data storage 
* datagrok-tools: add `suffix` key for package version hash 
* Wiki: how-to build an app. WIP: Subscribing to events 
* Wiki: how-to build an app. WIP: Scripting 
* Single scaffold alignment mode 'chem-scaffold': simplified 
* (Bug) Scaffold highlights have offsets 
* ClinicalCase \- a plugin for analyzing SDTM-based clinical studies data (WIP)
* Wiki: how-to build an app. WIP: REST endpoints, groups, sharing 
* Wiki: minor updates 
* Wiki: how-to build an app. WIP: Cleaning up 
* Biosignals 
* Ability to store data in a layout 
* Charts: Timelines (WIP)
* JS API: Ability to insert columns 
* Fixed release notes 
* Fixed URLs 
* Fixed release history 
* Trellis/filter screencast 
* Filter layout screencast 
* null 
* Wiki: how-to build an app. WIP: working with databases, storing dataframes 
* NLP package (WIP)
* Code snippet for using viewers outside of the TableView 
* Wiki: how-to build an app. WIP: adding to viewers 
* Wiki: how-to build an app. WIP: Authentication; cleanup 
* Wiki: how-to build an app. WIP: Functions, calling scripts 
* Wiki: how-to build an app. WIP: Application lifecycle 
* Fixed a typo 
* Wiki: how-to build an app. WIP: typos, simplifications 
* Show group tooltips 
* Wiki: how-to build an app. WIP: simplifications 
* Wiki: docker-compose instructions update 
* Fixed package template 
* Overview / Navigation (WIP)
* Implement row selection 
* Viewers: retired Dumbbell 
* Package: Chembl Browser (WIP)
* \* bumped RDKit_minimal version number to 2021.03_05 \* simplified _molIsInMolBlock function to detect if the passed molString is a molblock and renamed to _isMolBlock \* now generate_aligned-coords will return the matched substructure JSON if successful and empty string if not, so no need to to call get_substruct_match beforehand \* get_new_coords is used to generate new coordinates rather than roundtripping through SMILES. This is more efficient and ensures that CoordGen is used for coordinate generation 
* delete the WASM object before reassignment 
* (Bug) Chem: Keep mol rendering on a malformed scaffold 
* if generate_aligned_coords() returns "" and mol has no coordinates, fall back to get_new_coords() 
* (Bug) Chem: Single scaffold and column scaffolds interfere 
* Do not overwrite existing coordinates if substructure match fails 
* Disable zoom reset on click 
* Default dimensions for viewers 
* Wiki: UI (WIP) 
* Update ui.md 
* JS API samples: indentation fix 
* JS API: Ability to access user home connection and project 
* (Bug) JS-based column property panels are not shown 
* Chem: Column features 
* JS API: add renderer to ComboPopup 
* JS API: Ability to post audit record 
* MapProxy: Tags deletion; Falsish traps on read 
* Chem: Highlight a filtering scaffold 
* \- set the acceptFailure parameter to false when calling generate_aligned_coords \- set the dummyIsotopeLabels parameter to false to hide isotopic labels on R groups and other dummy atoms 
* (Bug) UI: ui.button on property panel looks like ui.bigButton 
* JS API: EntityMeta (WIP)
* Samples: iterating over a dataframe with .values(). 
* Wiki: how-to build an app. WIP: dataframe iteration 
* Chem: demo file for scaffold alignments, selected from Chembl by @ptosco 
* Chem: demo file for scaffold alignments 
* (Bug) DG throws error on startup 
* fixes RDKit GitHub issue #3852 
* Wiki: Layout management 
* Chem: react in panels to tags changes 
* Harmonized form 
* Create Physionet records fileViewer 
* CONTRIB.md tweaks. 
* A sample on combining several variables in one line chart. 
* Wiki: create descriptions for all packages 
* A sample for settings up conditional color coding. 
* Optimize population of a dataframe with calculated columns 
* Wiki: how-to build an app. Terminology cleanup 
* New example code 
* New example 
* JS API: onValueChanged for all controls 
* Wiki: BioSignals types update 
* Wiki: BioSignals method categories described in table 
* Add annotator functionality (WIP)
* Release notes fix 
* Add to selenium test for Filters check for deleting rows 
* JS API: Make tags keys / objects available 
* Typos and spaces fixes in samples. 
* Renamed files 
* Wiki: BioSignals goals and signal types described in table 
* Chord Diagram: `include nulls` property 
* (Bug) files method rename doesnt work properly  
* GROK: quality-of-life updates 
* GROK: repo cleanup 
* Clinical Case: validate values on column level (WIP)
* Chem: docs updates. 
* Chem package version 0.9.0, bump. 
* Adjusted the color coding example. 
* Download Physionet record through package 
* JS: Ability to create ViewLayout from ViewState 
* Timelines: Add `Reset View` item to context menu 
* datagrok-tools: add package `init` function template 
* Improving the choice-input.js example. 
* JS: Log API 
* Adjusted init-column.js example 
* GROK: new sequence viewer added 
* JSDoc fix 
* JS API: EntityMeta: dynamic resolution of object handlers (WIP)
* Added logging 
* Ability to set dialog position and size 
* JS API: width and height parameters for dialogs 
* GROK: paratopes added 
* GROK: RepertoireBrowser: save and load rows functionality 
* ClinicalCase: features 
* ClinicalCase readme update 
* GROK: RepertoireBrowser: improved load rows combobox  with async filelit updates 
* Show indicators box plot by taking folder with signals as input 
* images for ui.md 
* Selenium: Edit UI tests relative to new "id" and "name" for elements (WIP)
* Bio: Create SaveLoadRows Viewer for jbio project  
* Chem: Reverted full name to 'Chem', it looks better this way. 'Cheminformatics support' is in the description. 
* Adding utils library. 
* Timelines: refine column selection mechanism for SDTM variables 
* Adapter to export logs into CloudWatch (WIP)
* UI Doc small changes 
* (Bug) VirtualView hides scrollbar 
* Js Viewer: grow to all available space 
* Create BioSignals app with same functionality as in dialog 
* JS API: Matchers 
* (Bug) BioSignals: No input preset when adding new Indicator 
* Ability to create TypeScript packages (grok create PackageName -ts) 
* datagrok-tools: package templates compatibility 
* (Bug) BioSignals: executes extractor with alphabetically first name, instead of chosen 
* Wiki: link fix 
* JS API: Migrate to Typescript 
* Code formatting 
* Utils: table content validation engine (WIP)
* datagrok-tools: version bump 
* Biosignals: scripts should take column as input 
* Biosignals: round values of signal to 4 digits after comma 
* Trellis Plot: show R group columns in the dropdown list 
* (Bug) Trellis Plot: can't select a column with 100+ categories from the column menu 
* Made colorcoding-related tags hidden (prefixed with '.'). 
* Tags: make `color-coding-type` and `color-coding-conditional` private 
* update images format jpg to png 
* TextArea inside the ui.inputs: Fix the ui.textArea left-margin size 
* Fixed warnings 
* Biosignals: possibility to add custom script 
* Biosignals: get columns using tags instead of names 
* Biosignals: add possibility to change input Physionet record 
* Biosignals: speed up Physionet record loading 
* Fixed sampleCategories method 
* Biosignals: Update sampling frequency after resampling 
* Data Connectors: update PostgreSQL driver version  
* Sequence: WebLogoViewer (WIP)
* ApiSamples bump to test package auto deployment. 
* Fixed warning 
* Wiki: Onboarding Content (WIP)
* (Bug) Biosignals: TaskBarProgressIndicator continues animation after error 
* NLP package: UI adjustment 
* Tile Viewer and Form: integrate with dynamic renderers (WIP)
* Usage Analysis \- New UI 
* Updated public token 
* Shell: added table(name) and view(name) functions 
* Biosignals: use scripting capabilities dynamically 
* display molecules or fragments that fail to kekulized with delocalized bonds rather than showing a cross 
* Sequence package: work in progress 
* Grid: API: onColumnResized / onRowsResized events 
* Added missing react import 
* added missing WASM file 
* Chem package version 0.9.1, bump. 
* Create sequenceTranslator 
* Added new example on how to subscribe to resize events 
* Fixed unit-tests crash 
* Reflect the change in event naming onRowsResized 
* added TS API Grid::onRowsResized() 
* Internal docs: public API autogenerated files. 
* Cosmetics changes, addd new line between functions 
* GROK: CDR3 regions added to sequence view 
* cosmetics 
* (Bug) Box Plot \- marker size becomes too small, impossible to show tooltips on hi-res display 
* GROK: table connection added 
* Scatter plot: conditional color-coding 
* Box Plot: minor improvements 
* Box plot: conditional color coding 
* fixes the CoordGen regression that was causing 10x slower performance compared to the past 
* JS API: onPackageLoaded event 
* JS API: Object Handlers (WIP)
* Ability to discover and choose JS-based grid cell renderers based on a sem type 
* Optimization: faster BitSet construction based on logical function and two arguments 
* Removed long-running query 
* Wiki: debugging packages, debugging in WebStorm. 
* Chem package version 0.9.2, bump. 
* (Bug) Heat Map: columns do not resize 
* (Bug) Exception after adding filters on Table View 
* (Bug) Column format changes are not persisted in layouts 
* (Bug) S3 AES: Cant post a file 
* Newick file viewer 
* (Bug) Ability to disable routing 
* Fixed CSS 
* (Bug) Persist undocked viewers as part of layout (WIP)
* Empty value edge case in sampleCategories 
* (Bug) UI boolInput: onChange event executes code twice.  
* Fixed the bool input debouncing issue. 
* (Bug) Filters panel opens with an error message 
* Grid: make the color of selected rows lighter 
* Bar Chart: ability to select categories by clicking on the category label 
* Bar Chart: tooltips on category labels 
* Bar Chart: disable zoom on the X axis 
* NGL viewer (WIP)
* Prepare a proposal for WASM use and adoption in Datagrok core (WIP)
* [Filters: split categories containing a list of values by pipes](https://community.datagrok.ai/t/visualization-related-updates/521/11) 
* (Bug) Color Coding | Conditional: Ranges "20-30" and ">40" are always automatically determined for any numeric columns 
* Create sequenceTranslator app with instructions 
* (Bug) Bar chart: Error after opening project with bar chart 
* (Bug) Column Context Menu: `format` submenu is empty when the first value in a column is null 
* [(Bug) Scatter Plot: Properties: Value of the "Color" and "Size" properties change when you just hover over the columns from the list](https://community.datagrok.ai/t/changing-the-marker-size-in-the-scatterplot-closes-the-properties-panel/511) 
* (Bug) Scatter Plot: selected color should override custom conditional color 
* (Bug) Filters: Selected rows indication is not updated as selection changes 
* Release Notes: 0.89.36. 
* Release Notes: removing strange line-breaks. 
* Update release-history.md 
* Release Notes: 0.89.36 
* (Bug) Query View: Multiple vertical scrolls appear in query area 
* Edit tooltip: close the dialog when "Sketch form" is clicked 
* CVM: adopting the newer version of Miniconda 
* PhyloTreeViewer: JS-based phylogenetic tree visualisation (WIP)
* (Bug) Filters: categorical column filters still handle "selection.changed" event after being manually detached 
* (Bug) Filters: significant performance degradation when resetting the filter (ESC key) 
* (Bug) Filters: clicking on a category name does not always select it 
* Filters: show tooltip with a category name on mouseover 
* Dock type clarification 
* Wiki: how-to build an app. Namespacing connections 
* Wiki: how-to build an app. Pushing credentials (updates) 
* Wiki: Getting Started. Moving the item to the top 
* UI: ui.wait should act as ui.box 
* Bar Chart: enable selection on Shift + drag 
* Update chem.ts 
* (Bug) Chem | Info Panels: 3D panel displayed empty (WIP)
* Column: Custom value comparers 
* JS API: Column: Custom value comparers (WIP)
* ignore files 
* Wiki: best practices for working with data (WIP)
* Comment clean-up 
* enables including explicit Hs in substructure searches 
* (Bug) Qualified Numbers: Export to JSON loses the qualifier 
* Chem package version 0.9.3, bump. 
* (Bug) Vertical scrolling won't appear on vertical growth 
* (Bug) Density plot: impossible to show the X axis 
* (Bug) Filters: Saved filters do not filter rows in table after applying 
* Ability to share files 
* Sequence: better detection of nucleotide sequences 
* Detector.sampleCategories: performance fix: return false immediately if a column is not a string column 
* Power Pack: a set of commonly used Datagrok tools (WIP)
* added MSA sample file (FASTA alignment) 
* Prevent rings connected by rotatable bonds to flip when they belong to a constrained scaffold 
* Chem: Support RDKit-JS substructure filtering 
* sequenceTranslator: show molecule view 
* (Bug) Scatter Plot: Filtered-out current row still visible in plot (WIP)
* DataGrok API and system samples, add cross linking to improve connectivity (WIP)
* Removed the obsolete demo 'Weather' application 
* JS API: Grid: expose "invalidate" method 
* BioSignals: call functions from DSP package after switch to dynamic scripts calling 
* JS API: Ability to use file exporters defined as functions from packages 
* Chem: export .SD file 
* Sequence package for FASTA import of alignments and features (WIP)
* Fixed possible repository publish exception 
* Added webpack-cli dependency 
* Chem: categorical filtering option for a column with molecules 
* sequenceTranslator: expand the column to see whole sequence 
* description of "current rows" 
* update titles for "current rows" 
* update "current rows" 
* Wiki: headings harmonization 
* (Bug) Postgres type Numeric converts to int when using PostgresNet 
* Chem: don't attempt to .delete() null-molecules (no-render) 
* Core: custom cell renderers in filters (WIP)
* Updated beta users list. 
* (Bug) Viewers: Tooltip | Edit view tooltip... does not have any effect 
* Trellis Plot: do not create empty plots 
* changed text "To bring a row..." 
* Viewers: TreeViewer 
* Connectors: Redshift: Schema browsing (WIP)
* GrokConnect: move JdbcDataProvider methods isFloat(), isInteger() to child classes  
* (Bug) Closing a viewer leaves ghost panel 
* labs token 
* Enabled MT Wasmp with pthreads and workers 
* sequenceTranslator: use HTML tables for sequences and CMO codes 


# 2021-04-22 Dev build 0.89.40

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.40`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* changed text "To bring a row..." 
* sequenceTranslator: expand the column to see whole sequence 
* Viewers: TreeViewer 
* PhyloTreeViewer: JS-based phylogenetic tree visualisation (WIP)
* Connectors: Redshift: Schema browsing (WIP)
* GrokConnect: move JdbcDataProvider methods isFloat(), isInteger() to child classes  
* (Bug) Closing a viewer leaves ghost panel 
* labs token 


# 2021-04-22 Dev build 0.89.39

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.39`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* DataGrok API and system samples, add cross linking to improve connectivity (WIP)
* (Bug) Scatter Plot: Filtered-out current row still visible in plot (WIP)
* JS API: Ability to use file exporters defined as functions from packages 
* Chem: export .SD file 
* Sequence package for FASTA import of alignments and features (WIP)
* Fixed possible repository publish exception 
* Added webpack-cli dependency 
* PhyloTreeViewer: JS-based phylogenetic tree visualisation (WIP)
* Chem, RDKit-based (WIP)
* Chem: categorical filtering option for a column with molecules 
* sequenceTranslator: expand the column to see whole sequence 
* description of "current rows" 
* update titles for "current rows" 
* update "current rows" 
* JS API: Migrate to Typescript 
* Wiki: headings harmonization 
* (Bug) Postgres type Numeric converts to int when using PostgresNet (WIP)
* Chem: don't attempt to .delete() null-molecules (no-render) 
* Core: custom cell renderers in filters (WIP)
* Updated beta users list. 
* (Bug) Viewers: Tooltip | Edit view tooltip... does not have any effect 
* Trellis Plot: do not create empty plots 


# 2021-04-20 Dev build 0.89.38

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.38`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) Filters: Saved filters do not filter rows in table after applying 
* Ability to share files 
* Sequence: better detection of nucleotide sequences 
* Detector.sampleCategories: performance fix: return false immediately if a column is not a string column 
* Create sequenceTranslator app with instructions 
* Prepare a proposal for WASM use and adoption in Datagrok core (WIP)
* Wiki: Onboarding Content (WIP)
* Power Pack: a set of commonly used Datagrok tools (WIP)
* NwkViewer: JS-based phylogenetic tree visualisation (WIP)
* added MSA sample file (FASTA alignment) 
* Prevent rings connected by rotatable bonds to flip when they belong to a constrained scaffold 
* Chem: Support RDKit-JS substructure filtering 
* sequenceTranslator: show molecule view 
* (Bug) Scatter Plot: Filtered-out current row still visible in plot (WIP)
* DataGrok API and system samples, add cross linking to improve connectivity (WIP)
* Removed the obsolete demo 'Weather' application 
* JS API: Grid: expose "invalidate" method 
* BioSignals: call functions from DSP package after switch to dynamic scripts calling 
* (Bug) Persist undocked viewers as part of layout (WIP)
* NGL viewer (WIP)


# 2021-04-16 Dev build 0.89.37

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.37`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Create sequenceTranslator app with instructions (WIP)
* (Bug) Bar chart: Error after opening project with bar chart 
* (Bug) Column Context Menu: `format` submenu is empty when the first value in a column is null 
* JS API: Object Handlers (WIP)
* [(Bug) Scatter Plot: Properties: Value of the "Color" and "Size" properties change when you just hover over the columns from the list](https://community.datagrok.ai/t/changing-the-marker-size-in-the-scatterplot-closes-the-properties-panel/511) 
* (Bug) Scatter Plot: selected color should override custom conditional color 
* (Bug) Filters: Selected rows indication is not updated as selection changes 
* Release Notes: 0.89.36. 
* Release Notes: removing strange line-breaks. 
* Update release-history.md 
* Release Notes: 0.89.36 
* (Bug) Query View: Multiple vertical scrolls appear in query area 
* Edit tooltip: close the dialog when "Sketch form" is clicked 
* CVM: adopting the newer version of Miniconda 
* NwkViewer: JS-based phylogenetic tree visualisation (WIP)
* Prepare a proposal for WASM use and adoption in Datagrok core (WIP)
* (Bug) Filters: categorical column filters still handle "selection.changed" event after being manually detached 
* (Bug) Filters: significant performance degradation when resetting the filter (ESC key) 
* (Bug) Filters: clicking on a category name does not always select it 
* Filters: show tooltip with a category name on mouseover 
* Dock type clarification 
* Wiki: Harmonization (WIP)
* Wiki: how-to build an app. Namespacing connections 
* Wiki: how-to build an app. Pushing credentials (updates) 
* Wiki: Getting Started. Moving the item to the top 
* UI: ui.wait should act as ui.box 
* Wiki: Onboarding Content (WIP)
* NGL viewer (WIP)
* Bar Chart: enable selection on Shift + drag 
* Update chem.ts 
* (Bug) Chem | Info Panels: 3D panel displayed empty (WIP)
* Column: Custom value comparers 
* JS API: Column: Custom value comparers (WIP)
* ignore files 
* Wiki: best practices for working with data (WIP)
* Comment clean-up 
* enables including explicit Hs in substructure searches 
* ClinicalCase \- a plugin for analyzing SDTM-based clinical studies data (WIP)
* (Bug) Qualified Numbers: Export to JSON loses the qualifier 
* Chem, RDKit-based (WIP)
* Chem package version 0.9.3, bump. 
* (Bug) Vertical scrolling won't appear on vertical growth 


# 2021-04-14 Dev build 0.89.36

In this release, we've focused on enriching both the experience of the platform end-users and on-platform developers, as well as on connecting these two groups. For instance, custom-built cell renderers, useful for displaying molecules, nucleotide sequences or experiments results, now expand the platform in many places such as [tooltips and tile viewers](https://community.datagrok.ai/t/cheminformatics-updates/457/8). There are also [new viewers](https://community.datagrok.ai/t/visualization-related-updates/521/4), [several bar chart](https://community.datagrok.ai/t/bar-chart-color-by-category/516) [features](https://community.datagrok.ai/t/visualization-related-updates/521/10), [color coding features](https://dev.datagrok.ai/js/samples/grid/color-coding-conditional), and a few dozen of [JS API]((https://datagrok.ai/help/govern/audit#javascript-api)) [enhancements](https://community.datagrok.ai/t/javascript-api-updates/526/8) and [bug fixes](##bug-fixes).

## Highlights

* Rendering of chemical molecules with RDKit [in tooltips, forms, viewers' axis and tile viewer](https://community.datagrok.ai/t/cheminformatics-updates/457/8)
* Bar Chart: support for DateTime [columns categorization](https://community.datagrok.ai/t/visualization-related-updates/521/10), functional improvements
* [Conditional color coding](https://dev.datagrok.ai/js/samples/grid/color-coding-conditional) for the grid, scatter plot, box plot
* New [Timelines Viewer](https://community.datagrok.ai/t/visualization-related-updates/521/4)
* JS API improvements: Typescript support, new [events](https://community.datagrok.ai/t/javascript-api-updates/526/8) and [methods](https://community.datagrok.ai/t/javascript-api-updates/526/5)
* Better JS editing in the [JS fiddle](https://public.datagrok.ai/js) with IntelliSense and `async/await`

## Major features and Improvements

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

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.36`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

# 2021-04-08 Dev build 0.89.35

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.35`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Demo data: added Chembl-4000 csv. 
* JS API: Migrate to Typescript 
* Biosignals: use scripting capabilities dynamically (WIP)
* JS: Ability to get HTML editor for function 
* display molecules or fragments that fail to kekulized with delocalized bonds rather than showing a cross 
* (Bug) Heat Map: columns do not resize 
* Ability to create TypeScript packages (grok create PackageName -ts) 
* Wiki: Harmonization (WIP)
* ClinicalCase \- a plugin for analyzing SDTM-based clinical studies data (WIP)
* Updated public token 


# 2021-04-07 Dev build 0.89.34

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.34`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* JS API: Migrate to Typescript 
* JS: Ability to get HTML editor for function 
* Fixed analyzer warnings 
* Fixed NGL Viewer size 
* (Bug) Bar Chart: category axis labels overlap with the column selector 
* (Bug) Core: importDemoPackages won't work if there's already DemoPackages 
* (Bug) Property panel won't open when viewer properties are reached from the context menu 
* Shell: added table(name) and view(name) functions 
* (Bug) Heat Map: columns do not resize 
* Fixed a typo 


# 2021-04-06 Dev build 0.89.33

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.33`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) Bar Chart: cannot plot date columns on a categorical axis (WIP)
* (Bug) Trellis Plot: error rendering empty molecules (empty strings) 
* JS API: Migrate to Typescript 
* CVM: NLOpt temporarily off to investigate the resulting image first. 


# 2021-04-06 Dev build 0.89.32

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.32`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* JS API: Migrate to Typescript 
* Wiki: Onboarding Content (WIP)
* GROK:8349 Improved Bar Chart to allow date-time plotting as categorical 
* Logging: Ability to send nested JS Objects 
* (Bug) Bar Chart: cannot plot date columns on a categorical axis (WIP)
* Custom HTML-based cell renderers: investigate performance issues (WIP)
* Updated public token 
* CVM: addressing the March numpy issue. 
* Viewers: ability to work with value functions (WIP)
* CVM: Octave 6.2.0, compilation (GraphicsMagick issue). 
* Wiki: Harmonization (WIP)
* (Bug) Biosignals: TaskBarProgressIndicator continues animation after error 
* JS: Ability to get HTML editor for function 
* RangeSelector: fixed integer range mode (for trellis plot) 
* NLP package: UI adjustment 
* Trellis Plot: set upper limit for displayed categories 
* Tile Viewer and Form: integrate with dynamic renderers (WIP)
* CVM: NLOpt legacy package from a .deb file 
* Usage Analysis \- New UI 
* Added TypeScript build before publish to NPM 


# 2021-04-02 Dev build 0.89.31

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.31`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* (Bug) DockManager: the viewer size reduces on resize 
* (Bug) View | Tooltip: After turning off all columns, the tooltip is still displayed (WIP)
* Biosignals: possibility to add custom script 
* Custom HTML-based cell renderers: investigate performance issues (WIP)
* JS API: Migrate to Typescript (WIP)
* Aggregation: ability to specify value transformation function 
* Biosignals: get columns using tags instead of names 
* Biosignals: add possibility to change input Physionet record 
* Biosignals: speed up Physionet record loading 
* Fixed sampleCategories method 
* Biosignals: Update sampling frequency after resampling 
* Data Connectors: update PostgreSQL driver version  
* Deploy ApiSamples as part of demo projects. 
* Sequence: WebLogoViewer (WIP)
* ApiSamples bump to test package auto deployment. 
* Deleted junk error logs 
* Updated public token 
* Fixed warning 


# 2021-03-31 Dev build 0.89.30

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.30`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Fixed warnings 
* Updated public token 


# 2021-03-31 Dev build 0.89.29

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.29`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Add annotator functionality (WIP)
* GROK: new sequence viewer added 
* Wiki: Layout management 
* JS API: EntityMeta: dynamic resolution of object handlers (WIP)
* Added logging 
* Biosignals (WIP)
* Ability to set dialog position and size 
* JS API: width and height parameters for dialogs 
* (Bug) JS Editor: tail comment leads to exception 
* GROK: paratopes added 
* GROK: RepertoireBrowser: save and load rows functionality 
* ClinicalCase: features 
* Update ui.md 
* ClinicalCase readme update 
* (Bug) JS API: tabControl.clear() doesn't remove tabs content 
* GROK: RepertoireBrowser: improved load rows combobox  with async filelit updates 
* TabControl: Calling Widget.clear on panels that are about to get removed 
* "Selection to Filter" command 
* Show indicators box plot by taking folder with signals as input 
* images for ui.md 
* Selenium: Edit UI tests relative to new "id" and "name" for elements (WIP)
* Updated public token 
* Adapter to export logs into CloudWatch (WIP)
* Trellis: virtual scrolling 
* Trellis: updating the viewer type indicator when loaded from the project 
* Bio: Create SaveLoadRows Viewer for jbio project  
* (Bug) Box Plot does not work on demo dataset beer 
* Molecules rendering in tooltips (WIP)
* [Scatter Plot: ability to show only selected rows](https://community.datagrok.ai/t/allow-to-show-only-marked-selected-rows-in-the-scatterplot/512) 
* Chem: Reverted full name to 'Chem', it looks better this way. 'Cheminformatics support' is in the description. 
* Chord Diagram (WIP)
* (Bug) UX: Credentials management UI regressed 
* Adding utils library. 
* Timelines: refine column selection mechanism for SDTM variables 
* Charts: Timelines (WIP)
* Ability to add email to the group 
* Fixed a typo 
* UI Doc small changes 
* (Bug) VirtualView hides scrollbar 
* Js Viewer: grow to all available space 
* Property panel: fixed the appearance of the search icon 
* Grid: "settings" icon improvement 
* (Bug) Viewers: Tooltip | Edit view tooltip... does not have any effect 
* HtmlCellRenderer: Fixed function signatures. 
* Property: added tags 
* Harmonized viewer menus. 
* (Bug) setOptions({ viewerType: ... }); does not work as expected for Trellis plots  
* Always show dock panel header  
* Create BioSignals app with same functionality as in dialog 
* (Bug) Viewers: Exception when adding "Box Plot" and "Tree Map" if table has one row 
* JS API: Matchers 
* (Bug) Box Plot \- marker size becomes too small, impossible to show tooltips on hi-res display 
* ClinicalCase \- a plugin for analyzing SDTM-based clinical studies data (WIP)
* (Bug) BioSignals: No input preset when adding new Indicator 
* Ability to create TypeScript packages (grok create PackageName -ts) 
* datagrok-tools: package templates compatibility 
* (Bug) BioSignals: executes extractor with alphabetically first name, instead of chosen 
* Wiki: link fix 
* JS API: Migrate to Typescript (WIP)
* (Bug) Exception after adding Histogram to "TSLA" test dataset  
* Code formatting 
* Utils: table content validation engine (WIP)
* datagrok-tools: version bump 
* Biosignals: scripts should take column as input 
* Biosignals: round values of signal to 4 digits after comma 
* Trellis Plot: show R group columns in the dropdown list 
* (Bug) Trellis Plot: can't select a column with 100+ categories from the column menu 
* (Bug) Chem: RDKit renderer: tooltip does not work 
* "Set default tooltip" menu item 
* Made colorcoding-related tags hidden (prefixed with '.'). 
* Tags: make `color-coding-type` and `color-coding-conditional` private 
* update images format jpg to png 
* (Bug) Bar Chart: unexpected indent on plots with many categories 
* NumericalMatcher: ability to match empty values 
* TextArea inside the ui.inputs: Fix the ui.textArea left-margin size 
* Conditional color coding editor: input validation 
* Conditional color coding: improve performance 
* (Bug) Grid: empty cells are colored with conditional color-coding on 
* (Bug) Visual Query does not work for all supported providers 
* Custom HTML-based cell renderers: investigate performance issues (WIP)
* CVM: NLOpt + Octave 6.2.0 (WIP)
* (Bug) Query View: Incorrect display of line numbering 
* (Bug) UX: Flickering scroll bars (WIP)


# 2021-03-18 Dev build 0.89.28

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.28`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Add to selenium test for Filters check for deleting rows 
* (Bug) Detectors load on each table open 
* Viewers (Bar Chart, etc) tooltips improvements (N UX/UI) (WIP)
* more tests for group-by ops 
* Add annotator functionality (WIP)
* Auto-generated release-history tweaking 
* JS API: Make tags keys / objects available 
* Typos and spaces fixes in samples. 
* IntelliSense for the JS fiddle 
* Package: Chembl Browser (WIP)
* (Bug) Layout apply doesn't work on renamed columns 
* Harmonize ui framework (WIP)
* (Bug) dapi.log.events returns empty array 
* (Bug) Chem: Column width with structures is not correct when using OCL 
* wrap script in async in JsEditorView 
* Renamed files 
* Programming exercises (WIP)
* Chord Diagram (WIP)
* Grid: Usability: delete selected rows resets the current grid position to top (WIP)
* Clicking on icons should not change currently focused element 
* Wiki: BioSignals goals and signal types described in table 
* Preserve filter state when closing filter panel (WIP)
* Grid: Harmonize popup menu (WIP)
* (Bug) "Normalize" action does not work 
* Improved code comments (minor) 
* Chord Diagram: `include nulls` property 
* ClinicalCase \- a plugin for analyzing SDTM-based clinical studies data (WIP)
* Color coding (WIP)
* Categorical color coding editor: ability to bin values 
* Document the Datagrok release process protocol 
* NumericMatcher: allowing spaces in range definitions 
* /info/packages route 
* Improve error handling mechanism for visualizations 
* (Bug) files method rename doesnt work properly  
* (Bug) Saved function annotation doesn't delete old parameters from database 
* (Bug) Settings apply drops password parameters 
* GROK: quality-of-life updates 
* GROK: repo cleanup 
* (Bug) User can't login with OpenId if someone shared something to him by email 
* Improved comments on DataFrameViewer::onFrameAttached(..) 
* Clinical Case: validate values on column level (WIP)
* Chem: docs updates. 
* Chem package version 0.9.0, bump. 
* JS API: pass JsViewer in onContextMenu events 
* (Bug) tabControl with ui.wait will concatenate contents of tabs and show loader after second tab click 
* Adjusted the color coding example. 
* Download Physionet record through package (WIP)
* JS: Ability to create ViewLayout from ViewState 
* Timelines: Add `Reset View` item to context menu 
* Biosignals (WIP)
* Limit panel drag-and-drop so that the header always remains visible (WIP)
* [(Bug) Histogram: duplicated Y axis](https://community.datagrok.ai/t/cannot-easily-change-the-aggregation-of-the-histogram/509/4) 
* datagrok-tools: add package `init` function template 
* Improving the choice-input.js example. 
* JS: Log API 
* Adjusted init-column.js example 
* GROK: new sequence viewer added 
* JSDoc fix 


# 2021-03-10 Dev build 0.89.27

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.27`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Release notes fix 
* Add annotator functionality (WIP)
* Packages: Init function 
* Harmonize ui framework (WIP)
* /info/packages route 


# 2021-03-10 Dev build 0.89.26

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.26`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Wiki: BioSignals types update 
* Wiki: Harmonization (WIP)
* Wiki: Programming exercises (WIP)
* Package: Chembl Browser (WIP)
* Biosignals (WIP)
* Charts: Timelines (WIP)
* ClinicalCase - a plugin for analyzing SDTM-based clinical studies data (WIP)
* Harmonize ui framework (WIP)
* Wiki: BioSignals method categories described in table 
* (Bug) Package Properties ignored on repository publish 
* (Bug) Grid View: Selection display is drawing incorrect on high res.displays (MacOS) 
* (Bug) error: NullError: method not found: '_ddt$_name' on null 
* Improved bitset notification logic to respect current level of notifications (aka updateLevel) 
* [Bar Chart: add setting for excluding null values category on bar segments](https://community.datagrok.ai/t/bar-chart-color-by-category/516)
* (Bug) Filters: deletion of rows in grid/df results in incorrect filters 
* Viewers tooltips improvements (WIP)
* Add annotator functionality (WIP)
* JS API: preserve metadata while changing column type (WIP)
* Trellis Plot: show more granular X and Y axis ticks on line charts 


# 2021-03-07 Dev build 0.89.25

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.25`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* JS API: onValueChanged for all controls 
* Create Physionet records fileViewer (WIP)
* Harmonize ui framework (WIP)


# 2021-03-07 Dev build 0.89.24

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.24`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* fixes RDKit GitHub issue #3852 
* Functions: minor performance improvements 
* (Bug) TileViewer doesn't work from JS 
* Harmonize ui framework (WIP)
* JS API: EntityMeta (WIP)
* Wiki: Layout management 
* Wiki: Harmonization (WIP)
* Chem: react in panels to tags changes 
* Package: Chembl Browser (WIP)
* Optimize population of a dataframe with calculated columns (WIP)
* (Bug) DG throws error on startup 
* Chord Diagram (WIP)
* (Bug) PostFileToCvm doesn't load credentials 
* Harmonized form 
* Create Physionet records fileViewer (WIP)
* CONTRIB.md tweaks. 
* A sample on combining several variables in one line chart. 
* Wiki: create descriptions for all packages 
* A sample for settings up conditional color coding. 
* Wiki: how-to build an app. Terminology cleanup 
* New example code 
* New example 
* Biosignals (WIP)
* Customer OAUTH integration (WIP)
* Fixed analyzer warning 


# 2021-03-02 Dev build 0.89.23

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.23`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Harmonize ui framework (WIP)
* (Bug) DG throws error on startup 


# 2021-03-02 Dev build 0.89.22

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.22`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* JS API: files 
* Added new users to Developers team 
* Viewers: retired Dumbbell 
* Biosignals (WIP)
* (Bug) Query Transformations:  Actions cannot be edited in already saved queries 
* Harmonize ui framework (WIP)
* Package: Chembl Browser (WIP)
* Wiki: How to Build Custom Viewers 
* Programming exercises (WIP)
* Wiki: Harmonization (WIP)
* * bumped RDKit_minimal version number to 2021.03_05 * simplified _molIsInMolBlock function to detect if the passed molString is a molblock and renamed to _isMolBlock * now generate_aligned-coords will return the matched substructure JSON if successful and empty string if not, so no need to to call get_substruct_match beforehand * get_new_coords is used to generate new coordinates rather than roundtripping through SMILES. This is more efficient and ensures that CoordGen is used for coordinate generation 
* delete the WASM object before reassignment 
* Chem, RDKit-based (WIP)
* (Bug) Chem: Keep mol rendering on a malformed scaffold 
* if generate_aligned_coords() returns "" and mol has no coordinates, fall back to get_new_coords() 
* (Bug) Chem: Single scaffold and column scaffolds interfere 
* Do not overwrite existing coordinates if substructure match fails 
* Disable zoom reset on click 
* Overview / Navigation (WIP)
* NLP package (WIP)
* (Bug) Direct link to file leads to 404 error 
* Default dimensions for viewers 
* Wiki improvements 
* Nucleotide sequences for exercises: a-h1n1, sars-cov-2 
* GROK 8056 JS API: added and tested rename method 
* Check or uncheck the checkbox by clicking on their label 
* GROK 8056 JS API: added example 
* Wiki: UI (WIP) 
* Update ui.md 
* GROK 8057 JS API: Search pattern implementation on server-side 
* GROK 8057 JS API: typo 
* Updated public token 
* GROK 8056 JS API: samples grooming 
* JS API samples: indentation fix 
* JS API: Ability to access user home connection and project 
* (Bug) JS-based column property panels are not shown 
* S3: Correct URI depending on region 
* Chem: Column features 
* updated public token 
* JS API: add renderer to ComboPopup 
* (Bug) Grid | Tooltip: After sketching from, the tooltip loses its style 
* Charts: Timelines (WIP)
* (Bug) "Converting object to an encodable object failed" returning a dataframe in JS 
* JS API: Ability to post audit record 
* (Bug) Filter: category not available as a Filter to add (world_mortality) 
* (Bug) Docking sizes/coordinates recalculation are incorrect after View is removed from a split 
* Small MD edits to remove ending dots and add notes on MacOS dev.setup 
* added script to help run pub get for setup dev environment 
* MapProxy: Tags deletion; Falsish traps on read 
* Chem: Highlight a filtering scaffold 
* Datlas: give a meaningful error description when a required environment variable is not found 
* Updates to dev. setup documentation 
* Added script to build/deploy dart sub-projects (work-in-progress) 
* Added unix +x permission to service script 
* - set the acceptFailure parameter to false when calling generate_aligned_coords - set the dummyIsotopeLabels parameter to false to hide isotopic labels on R groups and other dummy atoms 
* Got rid of the deprecated GROK_WEB_ROOT. 
* GROK 8331 UI: name setter of AccordionPane fixed 
* (Bug) name setter of TabControl works in the consol but not in UI 
* GROK 7949 Rendered ui cards with svg molecules inside and added eventListeners 
* (Bug) UI: ui.button on property panel looks like ui.bigButton 
* (Bug) Train Model View: Button with style "disabled" available for clicking 
* Improvements to env.settings docs 
* Improved dev.env install documentation 
* JS API: EntityMeta (WIP)
* Better build helper scripts 
* Improved tech.docs on dev.env setup 
* Samples: iterating over a dataframe with .values(). 
* Wiki: how-to build an app. WIP: dataframe iteration 
* Chord Diagram (WIP)
* (Bug) Core: Cell selection triggers column selection with no visual result 
* Improved (minor) dev environment setup docs 
* Charts: Tree Viewer (WIP)
* Chem: demo file for scaffold alignments, selected from Chembl by @ptosco 
* Chem: demo file for scaffold alignments 
* GROK 7543 JS API: Entity fixed 
* public token updated 
* Established stricker dependencies on packages 
* Improved script to check for errors 
* Improved script error handling 


# 2021-02-15 Dev build 0.89.21

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.21`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Charts: Timelines (WIP)
* Fixed release history 
* Updated public token 
* Harmonize ui framework (WIP)
* Fixed typo 
* Bar Chart: includeNulls property (WIP)
* (Bug) Trellis Plot: The viewer data will not change during filtering
Trellis Plot: show only filtered data for linked dataframes (type: selection to filter) 
* Trellis/filter screencast 
* Filter layout screencast 
* Chem, RDKit-based (WIP)
* GROK 7836 JS API: added class files to dapi.js 
* GROK 7836 JS API: files functions registered at grok_api.dart 
* Wiki: how-to build an app. WIP: working with databases, storing dataframes 
* Biosignals (WIP)
* GROK 7836 JS API: changed FileInfo->Dynamic and connection small fix 
* NLP package (WIP)
* Code snippet for using viewers outside of the TableView 
* GROK 7836 JS API: added connection GUID/path handler 
* JS API: files 
* Wiki: how-to build an app. WIP: adding to viewers 
* Wiki: how-to build an app. WIP: Authentication; cleanup 
* Wiki: how-to build an app. WIP: Functions, calling scripts 
* Wiki: how-to build an app. WIP: Application lifecycle 
* Wiki: Harmonization (WIP)
* (Bug) Random error on sharing a package 
* Fixed a typo 
* Inspector: ignoring tooltip-related events 
* Wiki: how-to build an app. WIP: typos, simplifications 
* Show group tooltips 
* Wiki: how-to build an app. WIP: simplifications 
* GROK 7836 JS API: small fixes 
* Wiki: docker-compose instructions update 
* GROK 7836 JS API: made workable readaAsBytes(String), list, write, delete 
* GROK 7836 JS API: minor fix in connectors.dart (action move part) 
* Ability to run DG without visual root 
* Fixed package template 
* Overview / Navigation (WIP)
* (Bug) Files Browser: Text file previews do not align to view width 
* Chord Diagram (WIP)
* Implement row selection 
* GROK 7836 JS API: workable all ops + added tests 
* GROK 7836 JS API: writeAsText added 
* GROK 7836 JS API: Dapi_UserFiles_WriteAsText added 


# 2021-02-10 Dev build 0.89.20

## Latest Docker Images

* Datagrok: 
  *  `docker pull datagrok/datagrok:0.89.20`
  *  `docker pull datagrok/datagrok:latest`
  
* [Docker-Compose](admin/docker-compose.md)

## Addressed Issues

* Charts: Timelines (WIP)
* JS API: Ability to insert columns (WIP)
* Programming exercises (WIP)
* Updated exercises.md 
* User Data Storage overview 
* Update user-data-storage.md 
* Fixed links in user-data-storage.md 
* Js Packages: Replace GUID hash with build number
* Packages: delete debug version when user publishes release 
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
* (Bug) JS API: column name search box disappears from ColumnComboBox (WIP)
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
* (Bug) Notebooks: New notebooks are created empty 
* Add Chem package properties 
* Single scaffold alignment mode 'chem-scaffold': simplified 
* (Bug) DockManager: width rounding error 
* [(Bug) Bar Chart: null values break color-coding of stacked bars](https://community.datagrok.ai/t/bar-chart-color-by-category/516) 
* (Bug) Scaffold highlights have offsets 
* ClinicalCase - a plugin for analyzing SDTM-based clinical studies data (WIP)
* Chem, RDKit-based (WIP)
* [Histogram: add a dropdown for 'Color Aggr Type' property](https://community.datagrok.ai/t/cannot-easily-change-the-aggregation-of-the-histogram/509) 
* Filters: improve sizing
* Filters: add DatePicker 
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
* (Bug) Table | Tooltip: After sketching form, fill of the values ​​in tooltip becomes white and their font changes 
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
