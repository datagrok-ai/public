<-- TEST SCENARIO NAME should match the tested functionality name. In all other cases create a meaningful name -->
<-- The name of the test scenario should describe which part of the platform it is testing -->

# Manual Test: Delete Query

## Target functionality

Delete a data query that was previously created there

## Scenario

Read the table from the top, do not skip any steps until it is intentional, and you know what you are doing.
This scenario needs to be reproduced for different DBs (Oracle, MariaDB, MySQL, MS SQL)

| Step No. | Action                                                                 | Expected result            | Description                               |
|:--------:|:-----------------------------------------------------------------------|:---------------------------|-------------------------------------------|
|    1     | Find queries you create in [adding query test case](./adding-query.md) |                            | in queries browser or on "Databases" view |
|    2     | Call context menu  for for each of them                                |                            |                                           |
|    3     | Click on "Delete" from context                                         | Confirmation dialog open   |                                           |
|    4     | Click on YES in confirmation dialog                                    | Query successfully deleted |                                           |
|    5     | Check that queries has been deleted and is no longer present           |                            |                                           |
