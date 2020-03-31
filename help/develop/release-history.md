<!-- TITLE: Release History -->
<!-- SUBTITLE: -->

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