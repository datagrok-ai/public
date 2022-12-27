<-- TEST SCENARIO NAME should match the tested functionality name. In all other cases create a meaningful name -->
<-- The name of the test scenario should describe which part of the platform it is testing -->

# Manual Test: Adding connections

## Target functionality

Deleting a previously created Data connection

## Scenario

Read the table from the top, do not skip any steps until it is intentional, and you know what you are doing.
This scenario needs to be reproduced for different DBs (Oracle, MariaDB, MySQL, MS SQL)

| Step No. | Action                                                                                      | Expected result               | Description |
|:--------:|:--------------------------------------------------------------------------------------------|:-------------------------------|:------------|
|    1     | Find connections you create in [Add connection test-case](./adding-connections.md)          |                                |             | 
|    2     | Call context menu  for for each of them                                                     |                                |             |
|    3     | Click on “Delete”                                                                           | Confirmation dialog open       |             |
|    4     | Click on YES in confirmation dialog                                                         | Connction successfully deleted |             |
|    5     | Check that connections has been deleted and is no longer present in the connections browser |                                |             |
