<!-- TITLE: Tests: Web connector -->
<!-- SUBTITLE: -->

# Tests: Web connector

Web connector allows ability to connect to API services that provide external sources.

## Testing scenario

1. Open "Connection to data" view

1. Drag-and-drop swagger-file *"openweathermap_test.yaml"* from local storage to platform

* Process bar shows importing process
* After swagger-file importing message "Swagger file "openweathermap_test.yaml" is succesfully added into "Web"
  connectors" on the balun
* New connection "OpenWeatherMap_Test" was created for "Web" source

1. Expand "OpenWeatherMap_Test" connection in Source tree

* There is "5 day/3 hour Forecast By City Name" query in connection described in swagger-file

1. Double click on "5 day/3 hour Forecast By City Name" query

* Dialog with choice of query parameters is open
* Dialog have one field "Q" described in swagger-file
* "Q" field have default value "München,DE" described in swagger-file

1. Run "5 day/3 hour Forecast By City Name" query

* Query failed
* Message about empty API-key in balun
* In Source tree, under query, created failed query run (marked red)

1. Open "Edit connection" dialog for "OpenWeatherMap_Test" connection

1. Enter correct ApiKey to "API Key" field and save changes

1. Run again "5 day/3 hour Forecast By City Name" query with default parameters

* Query completed
* Returned table is open in platform
* In Source tree, under query, created completed query run (green)

1. Run again "5 day/3 hour Forecast By City Name" query with parameters that don't make sense (
   e.g. "test")

* Query completed
* Returned empty table that is open in platform

See also:

* [Data connection](../../access/data-connection.md)
* [Data connection Test](../../access/data-connection-test.md)
* [Data query](../../access/data-query.md)
* [Data query-test](../../access/data-query-test.md)
