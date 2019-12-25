<!-- TITLE: Tests: File connector -->
<!-- SUBTITLE: -->

# Tests: File connector

File shares are often used for storing and exchanging data, and Datagrok provides first-class support for them. Once a file share is mounted as a network drive on a server and registered within the platform, its content gets automatically indexed and can be browsed.

## Data exploring

1. Open "Connection to data" view

1. Expand the "Files" connector
   * There is demo connection
   
1. Expand demo connection in Source tree
   * There is files and directories from ```.../data/demo```

1. Select "northwind" folder for the [Property Panle](../features/property-panel.md)
   * "Details" tab shows file location relative to connection, size, date of creation and modification
   * In "Actions" tab folder can be downloaded locally as zip-archive

1. Click "Download Zip" action
   * Folder downloaded locally as zip-archive

1. Expand "northwind" folder in Source tree
   * There is structured list of files and folders inside "northwind" folder

1. Select "Northwind.xslx" file for the [Property Panle](../features/property-panel.md)
   * "Details" tab shows file location relative to connection, size, date of creation and modification
   * In "Actions" tab file can be downloaded locally, downloaded as zip-archive or opened in [Workspace](../features/workspace.md)

1. Click "Download" action
   * File downloaded locally

1. Click "Download Zip" action
   * File downloaded locally as zip-archive
   
1. Click "Open" action
   * File opened in [Workspace](../features/workspace.md)
   * Each sheet is opened as separate table
   
1. Expand "Northwind.xslx" file in Source tree
   * There is list of sheets inside file
   * Each sheet presented as separate table
   
1. Select "Categories" sheet for the [Property Panle](../features/property-panel.md)
   * "Details" tab shows number of rows and columns
   * In "Actions" tab sheet can be opened as separate table in [Workspace](../features/workspace.md) 

1. Click "Open" action
   * "Categories" table opened in [Workspace](../features/workspace.md)
   
1. Expand "Categories" sheet in Source tree
   * There is list of columns in "Categories" table
   * You can find columns data type selecting it on [Property Panle](../features/property-panel.md)
   
## Connection creating

1. Open "Add connection" dialog from context menu of "File" connector

1. Name new connection as "Test"

1. Set "Index Files" as true

1. In "Dir" field specify next path: ```/home/www/master/servergrok/data/formats``` and click "OK"
   * New file connection "Test" was created
   * File indexing from specified path has begun
   
1. Expand "Test" connection in Source tree after end of indexation
   * There is files and directories from ```...data/formats```
   
## Query creating

1. Open "Query View" from demo connection context menu
   * "Query View" is open
   * View shows list of files and folders from connection in form of tree
   
1. Click on "TSLA.csv" file in tree
   * Path to file "TSLA.csv" is added to "Path" relative to connection
   
1. Run query in "Query view" (icon on toolbar)
   * Query completed   
   * Result preview is shown in view
   
1. Click on "Northwind.xlsx" from "northwind" folder in tree
   * Path to file "Northwind.xlsx" is added to "Path" relative to connection  
   * "Sheet" field appeared for entering sheet from .xlsx file
   
1. Enter "Orders" to "Sheet" field and run query in "Query view" (icon on toolbar)
   * Query completed   
   * Result preview is shown in view ("Orders" table)
  
1. Click on "Products" sheet under "Northwind.xlsx" file
   * "Sheet" field filled with "Products" value
   
1. Save query with name "Test Query" 
   * New query appeared in tree under section "Queries" in demo connection
  
1. Run "Test Query" from Source Tree (by double-clicking)
   * "Products" table is added to [Workspace](../features/workspace.md)

## Common scenarios

1. Test functionality from previous sections with different file types. Supported file types can be found in [importing-data.md](../features/importing-data.md) file

1. Check common functionality for all data sources ([data-query-test](../entities/data-query-test.md))

1. Test transformations in this provider ([transformations-editor-test](../tests/transformation-editor-test.md))

1. Check non-functional aspects (UI, UX, Help, etc.)

See also:
 * [Files](../entities/connect/files.md)  
 * [Data connection](../entities/data-connection.md)
 * [Data connection Test](../entities/data-connection-test.md)
 * [Data query](../entities/data-query.md)
 * [Data query-test](../entities/data-query-test.md)
