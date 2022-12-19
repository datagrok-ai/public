<!-- TITLE: Tests: Data connection -->
<!-- SUBTITLE: -->

# Tests: Data connection

[Data connection](data-connection.md) defines a way to connect to a particular data source.

## Testing scenarios

### Source tree

1. Open *"Data Sources"* view

1. Open *Add new connection* dialog from context menu for "PostgresDart" source

* *Add new connection* dialog is open
* Required fields are highlighted in red
* In *Data source* field selected ```PostgreSQL``` and closed for edit
* ```OK``` not available

1. Enter "test_postgres" to "Name" field

1. Click on ```Test``` button

* Warning "Connection "test_postgres" is not available" on red balun

1. Fill other fields with following data: ```Server = localhost```, ```DB = northwind```
   , ```login = postgres```, ```password = ********```

1. Cli1. Click on ```Test``` button

* Message "Connection "test_postgres" is available" on green balun

1. Click ```OK```

* New connection "test_postgres" added to tree

1. Open *"Edit connection"* dialog from context menu of "test_postgres" connection

1. Change name to "new_test_postgres" on click ```OK```

* "test_postgres" connection name changed to "new_test_postgres"

1. Change other parameters of "new_test_postgres" connection with arbitrary data and save changes

* "new_test_postgres" now has changed parameters

1. Repeat 2-10 steps for DB sources (Oracle, MariaDB, PostresNet, MySQLm MS SQL)

1. Open *Add new connection* dialog from context menu for "Sparql" source

* *Add new connection* dialog is open
* Required fields are highlighted in red
* In *Data source* field selected ```Sparql``` and closed for edit
* ```OK``` not available

1. Enter "test_sparql" to "Name" field

1. Click on ```Test``` button

* Warning "Connection "test_sparql" is not available" on red balun

1. Fill other fields with following data: ```Endpoint = http://data.ontotext.com/repositories/data-last```
   , ```Requires  Server = true```, *Prefixes* field left empty

1. Cli1. Click on ```Test``` button

* Message "Connection "test_sparql" is available" on green balun

1. Click ```OK```

* New connection "test_sparql" added to tree

1. Open *"Edit connection"* dialog from context menu of "test_sparql" connection

1. Change name to "new_test_sparql" on click ```OK```

* "test_sparql" connection name changed to "new_test_sparql"

1. Change other parameters of "new_test_sparql" connection with arbitrary data and save changes

* "new_test_sparql" now has changed parameters

1. Drag ```openweathermap.yaml``` file from local storage to platform on “Connect to data” view

* New connection "OpenWeatherMap" was created in  "Web" source

1. Open *"Edit connection"* dialog from context menu of "OpenWeatherMap" connection

1. Fill "ApiKey" with correct key

1. Cli1. Click on ```Test``` button

* Message "Connection "OpenWeatherMap" is available" on green balun

1. Change name to "test_web" on click ```OK```

* "OpenWeatherMap" connection name changed to "test_web"

### Connection browser

1. Open "Connection browser" from **Admin | Data Connections**

1. Apply "Created by me" filter from "Filters" tab on [Toolbox](../datagrok/navigation.md#toolbox)

* Only connections created current user are displayed in browser, including those created in previous paragraph

1. Open *"Add new connection"* dialog from "Actions" tab tab on [Toolbox](../datagrok/navigation.md#toolbox)

1. Select alternately different source in "Data source" field

* when changing data source, fields required for corresponding source are added to dialog

1. Select ```PostgreSQL``` for "Data source" field

* Fields "Name", "Server", "Database", "Login", "Password" appeared on dialog

1. Fill fields with following data: ```Name = test_posgtes_sql" Server = localhost```
   , ```DB = northwind```, ```login = postgres```, ```password = ********```

1. Cli1. Click on ```Test``` button

* Message "Connection "test_postgres_sql" is available" on green balun

1. Click ```OK```

* New connection "test_postgres_sql" created and displayed in connections browser

1. Go to "Connect to data" view

* Created connection "test_postgres_sql" from previous step is displayed in tree under "
  PostgreSQL" source.

1. Call context menu to "northwind" demo connection and select "Browse queries"

* "Query browser" is open with in queries belonging to selected connection
* Patten "connection = {connection_id}" is entered in browser search bar

1. Return to connection browser

1. Select "northwind" (PostgreSQL) connection to display on [Property Panel](../datagrok/navigation.md#properties)

1. Open the "Details" tab in [Property Panel](../datagrok/navigation.md#properties)

* "Details" tab is open
* The correct information for all fields is displayed (Description, Server, Database, Created by, Created, Updated,
  Tags)

1. Open a list of [Queries](data-query.md) on [Property Panel](../datagrok/navigation.md#properties)

* "noethwind" connection queries are displayed

1. Open a list of [Projects](../datagrok/project.md)
   on [Property Panel](../datagrok/navigation.md#properties)

* Projects that include "northwind" connection are displayed

1. Delete all test [connections](data-connection.md) created in previous steps from their context menu

* When deleting connection, confirmation occurs in "Are you sure?" dialog
* Deleted connections are no longer displayed in browser and on [Property Panel](../datagrok/navigation.md#properties)

1. Go to "Connect to data" view

* There are no test connections in source tree that were deleted

See also:

* [Data source test](data-source-test.md)
* [Data job test](data-job-test.md)
* [Data queries test](data-query-test.md)
* [Data connection](data-connection.md)
