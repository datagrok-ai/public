<!-- TITLE: Tests: Data Source -->
<!-- SUBTITLE: -->

# Tests: Data Source

[User](../entities/user.md) in system can connect to data from different sources. 

[Data source](../entities/data-source.md) can represent a database, an Excel file, a CSV file, a web service, 
or basically anything that is capable of providing the data.

## Testing scenarios

1. Call the context menu of the data source
   * Adding a new [connection](../entities/data-connection.md) to the data source from context menu
   
1. Call context menu and select "Browse queries"   
   * "Query browser" is open with in queries belonging to selected connection
   * Patten "connection.dataSource = {data_source_name}" is entered in browser search bar  
   
1. Call context menu and select "Browse connections"   
   * "Connections browser" is open with in queries belonging to selected connection
   * Patten "dataSource = {data_source_name}" is entered in browser search bar      

1. Open a list of [connections](../entities/data-connection.md) on [Property Panel](../features/property-panel.md)
  
1. Open a list of [Jobs](../entities/data-job.md) on Property Panel
   
1. Open a list of Job Runs on [Property Panel](../features/property-panel.md)
   * Job Runs list is open. Display start time, status and name of Job

1. Open a list of [Queries](../entities/data-query.md) on [Property Panel](../features/property-panel.md)
 
1. Open a list of Query Runs on [Property Panel](../features/property-panel.md)
   * Query Runs list is open. Display start time, status and name

See also:
 * [Data Connections test](../tests/data-connection-test.md)
 * [Data Job test](../tests/data-job-test.md)
 * [Data Queries test](../tests/data-query-test.md)
 * [Data Sourse](../entities/data-source.md)
