<-- TEST SCENARIO NAME should match the tested functionality name. In all other cases create a meaningful name -->
<-- The name of the test scenario should describe which part of the platform it is testing -->

# Manual Test: Adding connections

## Target functionality

Connections browser - a special platfrom view on which all existing data connections are presented

## Scenario

Read the table from the top, do not skip any steps until it is intentional, and you know what you are doing.

| Step No. | Action                                                               | Expected result                       | Description                                             |
|:--------:|:---------------------------------------------------------------------|:--------------------------------------|:--------------------------------------------------------|
|    1     | Click on “Manage” section on sidebar and then click on “Connections” | "Connections" view is open            | This view is onnections browser                         |
|    2     | Type "new_test" in the search field                                  | "new_test_postgres" connections found | It connection from [test-case](./adding-connections.md) |
|    3     | Check all tabs on Property Panel for connection                      |                         |             | Select connection as current object                     |
|    4     | Expand “Filter” tab on Toolbox and check its content                 | Filter templates are displayed        |                                                         |
|    5     | Expand “Actions” tab on Toolbox and check its content                | "Add new connection" action is shown  |                                                         |
