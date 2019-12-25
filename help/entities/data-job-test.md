<!-- TITLE: Tests: Data Job -->
<!-- SUBTITLE: -->

# Tests: Data Job

[Data job](../entities/data-job.md) defines all actions that are needed to produce a dashboard. 
Each [data job](../entities/data-job.md) consists of the following parts that are executed to produce a dataset.

1. Open a "Create Data Job" dialog from context menu of "test" [data query](../entities/data-query.md). 

1. Name a new [job](../entities/data-job.md) "Test Data Job" and create them
   * New Data Job created with the name "Test Data Job" 

1. Open "Data Jobs" from "Admin" menu
   * "Data Jobs" is open
   * You can change the view, use sort and search 

1. Use search to find the "Test Data Job" 

1. Run "Test Data Job" 
   * Data Job completed
   * New dataset was created 

1. Open the "Details" tab in [Property Panel](../features/property-panel.md)
   * "Details" tab is open
   * The correct information for all fields is displayed (Author, Created, Updated, Connection, Queries, Last Run, Tags)

1. Open the "Queries" tab in [Property Panel](../features/property-panel.md)
   * "Queries" tab is open
   * Display queries which are included in Data Job

1. Open the "Run" tab in [Property Panel](../features/property-panel.md)
   * "Run" tab is open
   * You can select [query](../entities/data-query.md) parameters here, run Data Job and see the parameters history

1. Open the "History" tab in [Property Panel](../features/property-panel.md)
   * "History" tab is open
   * Display information about running of [jobs](../entities/data-job.md)
   * Here you can see the status and start time 

1. Open the "Statistics" tab in [Property Panel](../features/property-panel.md)
   * "Statistics" tab is open
   * Display information about runs count, average time, first and last runs

1. Open the "Activity" tab in [Property Panel](../features/property-panel.md)
   * "Activity" tab is open
   * Display information about actual actions with result of job execution

1. Open the "Shared with" tab in [Property Panel](../features/property-panel.md)
   * "Shared with" tab is open
   * Display users and users groups which this job is available 

See also:
 * [Data Connections test](../tests/data-connection-test.md)
 * [Data Query test](data-query-test.md)
 * [Data Sourse test](../tests/data-source-test.md)
 * [Data Job](data-job.md)
 