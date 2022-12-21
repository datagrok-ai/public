<-- TEST SCENARIO NAME should match the tested functionality name. In all other cases create a meaningful name -->
<-- The name of the test scenario should describe which part of the platform it is testing -->

# Manual Test: Project Browser

## Target functionality

Service view of the platform, such as an entity gallery. This view is needed to store and display projects uploaded to the server.

## Scenario

Read the table from the top, do not skip any steps until it is intentional, and you know what you are doing.

| Step No. | Action                                                             | Expected result                                         | Description               |
|:--------:|:-------------------------------------------------------------------|:--------------------------------------------------------|:--------------------------|
|    1     | Open a "Projects" view for "Data" section on Sidebar               | "Projects" view is open                                 |                           |
|    2     | Use search to find the "demog" project | The example popup appears | Project "demog" displayed in the project browser        |                           |
|    3     | Open "Share" dialog for the "demog" project from its context menu  | "Share" dialog is open"                                 |                           |
|    4     | Check "Details" for the "demog" project from its context menu      | "Details" view for "demog" projectc is open             |                           |
|    5     | Check all tabs on Property Panel for "demog" project               | Tabs content is correct and there are no errors in them |                           |
