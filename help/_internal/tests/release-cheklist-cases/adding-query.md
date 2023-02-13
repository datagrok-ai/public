<-- TEST SCENARIO NAME should match the tested functionality name. In all other cases create a meaningful name -->
<-- The name of the test scenario should describe which part of the platform it is testing -->

# Manual Test: Adding Query

## Target functionality

Data query defines which data should be extracted from the data source. For databases, that would typically be the SQL query; for Excel file, that would be sheet name, etc. 

## Scenario

Read the table from the top, do not skip any steps until it is intentional, and you know what you are doing.
This scenario needs to be reproduced for different DBs (Oracle, MariaDB, MySQL, MS SQL)

| Step No. | Action                                                                  | Expected result                                      | Description            |
|:--------:|:------------------------------------------------------------------------|:-----------------------------------------------------|------------------------|
|    1     | Open "Data" section on sidebar and then click on "Databases" item       | "Databases" view is open                             |                        | 
|    2     | Open context menu for "PostgresDart\northwind” and click on "Add query" | Query View is open                                   |                        |
|    3     | Enter "test_query" to "Name" field                                      |                                                      |                        |
|    4     | Enter query text to corresponding field                                 |                                                      | select * from products |
|    5     | Click on "Play" button on toolbar                                       | Query was executed and result preview was displayed  |                        |
|    6     | Click on "Run query…" from "Actions" tab on Toolbox                     | Query was executed and result table open in Datagrok |                        |
|    7     | Click “Save” button on toolbar                                          | New query saved                                      |                        |