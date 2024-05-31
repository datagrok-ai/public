<!-- TITLE: Tests: Data query -->
<!-- SUBTITLE: -->

# Tests: Data query

[Data query](access.md#data-query) defines which data should be extracted from the
[data source](databases/connectors/connectors.md). For databases, that would typically be the SQL query; for Excel file, that
would be sheet name, etc. The system allows you to create parameterized and non-parameterized [queries](access.md#data-query).

## Testing scenarios

1. Create a new non-parameterized [query](access.md#data-query)  with name = "Test"
   (example query: select * from Products)

1. Open "Data Queries" from "Admin" menu

* "Data Queries" is open.
* You can change the view, use sort and search.

1. Use search in "Data Queries" to find the "Test" query

* Query "Test" will find.

1. Call context menu for [data query](access.md#data-query). Test the functionality of the context menu items.

* Context menu is open.
* All its elements work correctly

1. Open "Query view" for "Test" query, edit it and run without saving

* Query completed.

1. Enter an incorrect query and run it

* Warning about errors in the query is shown

1. Open "Data Query Runs" from "Admin" menu

* "Data Queries Runs" is open
* You can change the view, use sort and search
* The list shows the status of run and start time

1. Run [query](access.md#data-query) from different places (context
   menu, [Context Panel](../datagrok/navigation/panels/panels.md#context-panel), by double-clicking)

1. Find in the list [queries](access.md#data-query) from the previous paragraphs

* Queries are present in the list
* The runtime is correct
* Queries with errors are marked with color coding

1. Open the "Details" tab in [Context Panel](../datagrok/navigation/panels/panels.md#context-panel)

* "Details" tab is open
* The correct information for all fields is displayed (Author, Created, Updated, Connection, Tags, Params – for
  parameterized queries)

1. Open the "Query" tab in [Context Panel](../datagrok/navigation/panels/panels.md#context-panel)

* "Query" tab is open
* Display current text of query

1. Open the "Activity" tab in [Context Panel](../datagrok/navigation/panels/panels.md#context-panel)

* "Activity" tab is open
* Display information about tables that were created after the query was executed.

1. Run the "Test" [query](access.md#data-query)

* Look on "Activity" tab
* New record about actually query is added here

1. Open the "Usage" tab in [Context Panel](../datagrok/navigation/panels/panels.md#context-panel)

* "Usage" tab is open. Display query runs statistics.

1. Delete "Test" [query](access.md#data-query)

* "Test" [query](access.md#data-query) deleted successfully
* No any errors and exceptions

1. Create a new parameterized [query](access.md#data-query)  with name = "Test_param"

* New parameterized [query](access.md#data-query) created
* Check this in "Data Queries" in menu "Admin"

1. Open the "Run" tab in [Context Panel](../datagrok/navigation/panels/panels.md#context-panel) for "
   Test_param" [query](access.md#data-query)

* "Run" tab is open, there are fields for entering parameters

1. Run parameterized [query](access.md#data-query) with correct parameters

* Parameterized [query](access.md#data-query) is run
* Check this in "Data Query Runs" in menu "Admin"

1. Run parameterized [query](access.md#data-query) with not-correct parameters (nulls, incorrect types of input data,
   incomplete data for strings, negative values for numbers, error datetime, etc.

* Warning about incorrect input parameters is shown to user.

The scenarios listed in clauses 2-17 must also be performed for parameterized queries

_example parameterized [query](access.md#data-query):_

```
--input: int employeeId = 5
--input: string shipVia = = 3 {pattern: int}
--input: double freight = 10.0
--input: string shipCountry = "France" {choices: ["France", "Germany", "USA", "Finland"]}
--input: string shipCity = "starts with r" {pattern: string}
--input: bool freightLess1000 = true
--input: datetime requiredDate = "1/1/1995"
--input: string orderDate = "after 1/1/1995" {pattern: datetime}
SELECT * FROM Orders WHERE (employeeId = @employeeId)
   AND (freight >= @freight)
   AND @shipVia(shipVia)
   AND ((freight < 1000) OR NOT @freightLess1000)
   AND (shipCountry = @shipCountry)
   AND @shipCity(shipCity)
   AND @orderDate(orderDate)
   AND (requiredDate >= @requiredDate)*
```

See also:

* [Data connections test](data-connection-test.md)
* [Data job test](data-job-test.md)
* [Data source test](data-source-test.md)
* [Data query](access.md#data-query)
