<-- TEST SCENARIO NAME should match the tested functionality name. In all other cases create a meaningful name -->
<-- The name of the test scenario should describe which part of the platform it is testing -->

# Manual Test: Edit Query

## Target functionality

Data query defines which data should be extracted from the data source. For databases, that would typically be the SQL query; for Excel file, that would be sheet name, etc. 

## Scenario

Read the table from the top, do not skip any steps until it is intentional, and you know what you are doing.
This scenario needs to be reproduced for different DBs (Oracle, MariaDB, MySQL, MS SQL)

| Step No. | Action                                                                             | Expected result                                          | Description              |
|:--------:|:-----------------------------------------------------------------------------------|:---------------------------------------------------------|--------------------------|
|    1     | Open context menu for query created in [adding query test case](./adding-query.md) | "Databases" view is open                                 |                          |
|    2     | Click on "Edit…" from context menu | Query View is open                            | Query view is open                                       |                          |
|    3     | Change name to "new_test_query"                                                    |                                                          |                          |
|    4     | Change the query text                                                              |                                                          | to: select * from orders |
|    5     | Click on "Play" button on toolbar                                                  | Query was executed and new result preview was displayed  |                          |
|    6     | Click on "Run query…" from "Actions" tab on Toolbox                                | Query was executed and new result table open in Datagrok |                          |
|    7     | Click “Save” button on toolbar                                                     | Query saved with changes                                 |                          |