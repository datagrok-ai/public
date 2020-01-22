<!-- TITLE: Release History -->
<!-- SUBTITLE: -->

# 2020.01.20 Service update

## Latest docker images
* Datagrok (new): 766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok:1.0.X-da74169
* CVM (old): 766822877060.dkr.ecr.us-east-2.amazonaws.com/grok_cvm:1.0.31-6261445

## Addressed issues
* GrokCon* nect: Denodo integration
* GrokConnect: ODBC driver
* JS API examples: manual master-detail linking
* Chembl project: moved to public repository
* GrokConnect: MSSQL collision on usage @<> syntax with query parameters
* GrokConnect: MSSQL provider - "the value is not set for the parameter #2"
* Support schemas for PostgreSql provider
* Query: Improve parameter editor
* Web provider: Editing some demo connections throws an exception
* Added options to project entity
* DataQuery: Parameter values casting is missing on server and java
* Data connections: Sparql connections testing does not work correctly
* Optimized DB queries
* Predictive Modeling: Train: Table field should be before the column selection fields (H2O)
* upload.dart fails with self-signed SSL certificates
* Predictive Modeling: Chemprop models do not train
* TreeView: Item context menu exception
* Grok Connect: Db Table PP: Content: error
* GrokConnect: Schema table content errors
* grok.query(query, params): make params optional

# 2020.01.15 Stable version

## Latest docker images
* Datagrok (new): 7766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok:1.0.65-8dc0463
* CVM (old): 766822877060.dkr.ecr.us-east-2.amazonaws.com/grok_cvm:1.0.31-6261445

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
* Datagrok (new): 766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok:1.0.X-0dc9ba6
* CVM (old): 766822877060.dkr.ecr.us-east-2.amazonaws.com/grok_cvm:1.0.31-6261445

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
* Datagrok (new): 766822877060.dkr.ecr.us-east-2.amazonaws.com/datagrok:1.0.48-5791556
* CVM (new): 766822877060.dkr.ecr.us-east-2.amazonaws.com/grok_cvm:1.0.31-6261445

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