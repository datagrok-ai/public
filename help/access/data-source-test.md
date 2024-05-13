<!-- TITLE: Tests: Data source -->
<!-- SUBTITLE: -->

# Tests: Data source

[User](../govern/user.md) in system can connect to data from different sources.

## Testing scenarios

1. Call the context menu of the data source

* Adding a new [connection](access.md#data-connection) to the data source from context menu

1. Call context menu and select "Browse queries"

* "Query browser" is open with in queries belonging to selected connection
* Patten "connection.dataSource = {data_source_name}" is entered in browser search bar

1. Call context menu and select "Browse connections"

* "Connections browser" is open with in queries belonging to selected connection
* Patten "dataSource = {data_source_name}" is entered in browser search bar

1. Open a list of [connections](access.md#data-connection)
   on [Context Panel](../datagrok/navigation/panels/panels.md#context-panel)

1. Open a list of Jobs on Context Panel

1. Open a list of Job Runs on [Context Panel](../datagrok/navigation/panels/panels.md#context-panel)

* Job Runs list is open. Display start time, status and name of Job

1. Open a list of [Queries](access.md#data-query) on [Context Panel](../datagrok/navigation/panels/panels.md#context-panel)

1. Open a list of Query Runs on [Context Panel](../datagrok/navigation/panels/panels.md#context-panel)

* Query Runs list is open. Display start time, status and name

See also:

* [Data connections test](data-connection-test.md)
* [Data job test](data-job-test.md)
* [Data queries test](data-query-test.md)
