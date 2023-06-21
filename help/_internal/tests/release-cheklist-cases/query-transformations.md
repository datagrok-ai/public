<-- TEST SCENARIO NAME should match the tested functionality name. In all other cases create a meaningful name -->
<-- The name of the test scenario should describe which part of the platform it is testing -->

# Manual Test: Query Transformations

## Target functionality

After the data is retrieved from the data provider by data query, it can be transformed using Datagrok functions.

## Scenario

Read the table from the top, do not skip any steps until it is intentional, and you know what you are doing.
This scenario needs to be reproduced for different DBs (Oracle, MariaDB, MySQL, MS SQL)

| Step No. | Action                                                                | Expected result                                                                    | Description |
|:--------:|:----------------------------------------------------------------------|:-----------------------------------------------------------------------------------|-------------|
|    1     | Open "Data" section on sidebar and then click on "Databases" item     | "Datasets" view is open                                                            |             |
|    2     | Call context menu for "Products" query from "PostgreDart > northwind" |                                                                                    |             |
|    3     | Click on "Edit…" from context menu                                    | Query view is open                                                                 |             |
|    4     | Click on "Transformations" tab in Query View                          | Query view switched to "Transformations" page                                      |             |
|    5     | Click on "Add new column" from actions list and add new column        | Added step "add new column" to transformation script                               |             |
|    6     | Open Cluster dialog from "ML | Custer…" menu                          | "Cluster" dialog is open                                                           |             |
|    7     | Click OK in "Cluster dialog" with default parameters                  | Added step "cluster" to transformation script                                      |             |
|    8     | On the top toolbar, click the **Play** button                         | Query was completed and as a result there are transformations added earlier        |             |
|    9     | On the top toolbar, click the **Save** button                         | Query saved |                                                                      |             |
|    10    | Run the query saved in step 11                                        | Returned dataframe contains the result of executing functions from transformations |             |