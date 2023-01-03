<-- TEST SCENARIO NAME should match the tested functionality name. In all other cases create a meaningful name -->
<-- The name of the test scenario should describe which part of the platform it is testing -->

# Manual Test: Query Browser

## Target functionality

Query browser - a special platfrom view on which all existing data queries are presented

## Scenario

Read the table from the top, do not skip any steps until it is intentional, and you know what you are doing.

| Step No. | Action                                                         | Expected result                                              | Description                 |
|:--------:|:---------------------------------------------------------------|:-------------------------------------------------------------|-----------------------------|
|    1     | Click on "Data" section on sidebar and then click on "Queries" | "Queries" view is open                                       |  This view is query browser |
|    2     | Type "new_test" in the search field                            | Query from [adding query test case](./adding-query.md) found |                             |
|    3     | Check all tabs on Property Panel for query                     | Content of all tabs is correct and up-to-date                |                             |
|    4     | Expand “Filters” tab on Toolbox and check its content          | Content of all tabs is correct                               |                             |
