<-- TEST SCENARIO NAME should match the tested functionality name. In all other cases create a meaningful name -->
<-- The name of the test scenario should describe which part of the platform it is testing -->

# Manual Test: Adding connections

## Target functionality

Data connection is used for accessing data in a particular data source. Connection parameters depend on the data source. Typically, you would need to provide server name and login credentials.

## Scenario

Read the table from the top, do not skip any steps until it is intentional, and you know what you are doing.
This scenario needs to be reproduced for different DBs (Oracle, MariaDB, MySQL, MS SQL)

| Step No. | Action                                                                                          | Expected result                                        | Description |
|:--------:|:------------------------------------------------------------------------------------------------|:-------------------------------------------------------|:------------|
|    1     | Open context menu for connection created in [Add connection test-case](./adding-connections.md) |                                                        |             | 
|    2     | Click on "Edit connection" from context menu                                                    | "Edit connection" dialog is open                       |             |
|    3     | Change name to "new_test_postgres" on click OK                                                  | Сonnection now has a new name                          |             |
|    4     | Click on "Test" button                                                                          |                                                        |             |
|    5     | Change other parameters of "new_test_postgres" connection with arbitrary data and save changes  | Changes to all parameters are saved for the connection |             |
