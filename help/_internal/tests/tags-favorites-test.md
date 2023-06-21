<!-- TITLE: Tests: Tags -->
<!-- SUBTITLE: -->

# Tests: Tags

Platform provides ability to add custom tags for entities.

## Testing scenario

1. Open table in platform (from query or local file)

1. Open [Upload](../../datagrok/upload-project-test.md) dialog

1. Click on "Tag editor" button

* Tag edit field is open

1. Enter "test" in tag edit field and click ```Enter```

* \#test tag added to project

1. Upload project with name "Tags_Test"

1. Open "Datasets..." page on "Welcome" view and open "Tags_Test"
   on [Context Panel](../../datagrok/navigation.md#context-panel)

* \#test tag is present in "tags" field in "Details" tab

1. Click on "Tag editor" button near "tags" field in "Details" tab and add \#test2 tag

* One more tag \#test2 added to project "Tags_Test"

1. Add "Tags_Test" project to "Favorites" by clicking star near it's name
   on [Context Panel](../../datagrok/navigation.md#context-panel)

* On toolbar there appeared drop-down menu containing added objects to "favorite"
* When adding object to "Favorites" "star" icon flashes

1. Open **Admin | Data Connections**

1. Open "northwind" (PostgreSQL) connection on [Context Panel](../../datagrok/navigation.md#context-panel)

* There is \#demo tag in "tags" field in "Details" tab

1. Add * \#test tag by using "Tag editor" in "tags" field in "Details" tab

* \#test tag added to "northwind" (PostgreSQL) connection

1. Add "northwind" (PostgreSQL) connection to "Favorites" by clicking star near it's name
   on [Context Panel](../../datagrok/navigation.md#context-panel)

* "northwind" (PostgreSQL) connection added to dropdown menu "Favorites" on toolbar
* When adding object to "Favorites" "star" icon flashes

1. Open **Admin | Data Queries**

1. Open "Products" (PostgreSQL, northwind) query on [Context Panel](../../datagrok/navigation.md#context-panel)

* There is \#demo tag in "tags" field in "Details" tab

1. Add * \#test tag by using "Tag editor" in "tags" field in "Details" tab

* \#test tag added to "Products" (PostgreSQL, northwind) query

1. "Products" (PostgreSQL, northwind) query to "Favorites" by clicking star near it's name
   on [Context Panel](../../datagrok/navigation.md#context-panel)

* "Products" (PostgreSQL, northwind) query added to dropdown menu "Favorites" on toolbar
* When adding object to "Favorites" "star" icon flashes

1. Open **Admin | Data Jobs**

1. Open "Save Datlas Health statusd" job on [Context Panel](../../datagrok/navigation.md#context-panel)

* There is \#demo tag in "tags" field in "Details" tab

1. Add * \#test tag by using "Tag editor" in "tags" field in "Details" tab

* \#test tag added to "Save Datlas Health status" Job

1. Add "Save Datlas Health status" job to "Favorites" by clicking star near it's name
   on [Context Panel](../../datagrok/navigation.md#context-panel)

* "Save Datlas Health status" job added to dropdown menu "Favorites" on toolbar
* When adding object to "Favorites" "star" icon flashes

1. Open **Tools | Scripting | Browse Scripts**

1. Open "ACF" script on [Context Panel](../../datagrok/navigation.md#context-panel)

* There is \#demo tag in "tags" field in "Details" tab

1. Open "ACF" script for editing through double-click on it

* Script editing view opened in new tab

1. Add "test" tag to "\#tags:" line, separated by comma and save script changes

* \#test tag added to "ACF" script

1. Add "ACF" script to "Favorites" by clicking star near it's name
   on [Context Panel](../../datagrok/navigation.md#context-panel)

* "ACF" script added to dropdown menu "Favorites" on toolbar
* When adding object to "Favorites" "star" icon flashes

1. Open **Tools | Predictive Modeling | Browse Models**

1. Open "Predict country" model on [Context Panel](../../datagrok/navigation.md#context-panel)

* There is \#demo tag in "tags" field in "Details" tab

1. Add * \#test tag by using "Tag editor" in "tags" field in "Details" tab

* \#test tag added to "Save Datlas Health status" Job

1. Add "Predict country" model to "Favorites" by clicking star near it's name
   on [Context Panel](../../datagrok/navigation.md#context-panel)

* "Predict country" model added to dropdown menu "Favorites" on toolbar
* When adding object to "Favorites" "star" icon flashes

1. Open profile view

* Counter near "Favorites" tab on Toolbox shows number of favorite entities

1. Go to Favorites" tab on Toolbox

* In open view all entities, added in previous steps, are present

1. Choose arbitrary entity from favorite on [Context Panel](../../datagrok/navigation.md#context-panel)

1. Click on \#demo tag in "tags" field in "Details" tab

* Tree with entities was opened on [Context Panel](../../datagrok/navigation.md#context-panel), to which
  # test tag was added in previous steps.
