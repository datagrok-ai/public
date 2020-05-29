<!-- TITLE: Tests: Data Source -->
<!-- SUBTITLE: -->

# Tests: Data Source

[User](../govern/user.md) in system can connect to data from different sources. 

## Testing scenarios

1. Call the context menu of the data source
   * Adding a new [connection](data-connection.md) to the data source from context menu
   
1. Call context menu and select "Browse queries"   
   * "Query browser" is open with in queries belonging to selected connection
   * Patten "connection.dataSource = {data_source_name}" is entered in browser search bar  
   
1. Call context menu and select "Browse connections"   
   * "Connections browser" is open with in queries belonging to selected connection
   * Patten "dataSource = {data_source_name}" is entered in browser search bar      

1. Open a list of [connections](data-connection.md) on [Property Panel](../overview/property-panel.md)
  
1. Open a list of [Jobs](data-job.md) on Property Panel
   
1. Open a list of Job Runs on [Property Panel](../overview/property-panel.md)
   * Job Runs list is open. Display start time, status and name of Job

1. Open a list of [Queries](data-query.md) on [Property Panel](../overview/property-panel.md)
 
1. Open a list of Query Runs on [Property Panel](../overview/property-panel.md)
   * Query Runs list is open. Display start time, status and name

See also:
 * [Data Connections test](../tests/data-connection-test.md)
 * [Data Job test](../tests/data-job-test.md)
 * [Data Queries test](../tests/data-query-test.md)
