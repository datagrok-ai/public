<-- TEST SCENARIO NAME should match the tested functionality name. In all other cases create a meaningful name -->
<-- The name of the test scenario should describe which part of the platform it is testing -->

# Manual Test: Delete Project

## Target functionality

Delete a project from a server that was previously uploaded there

## Scenario

Read the table from the top, do not skip any steps until it is intentional, and you know what you are doing.

| Step No. | Action                                                                         | Expected result                    | Description                                      |
|:--------:|:-------------------------------------------------------------------------------|:-----------------------------------|:-------------------------------------------------|
|    1     | Find the project you upload in [Upload project test-case](./upload-project.md) | Project selected as current object |                                                  | 
|    2     | Call context menu for the found project                                        | Context menu opened                |                                                  |
|    3     | Click on "Delete project"                                                      | Confirmation dialog open           |                                                  |
|    4     | Click YES in confirmation dialog                                               | Project removed from server        | Project is no longer present in projects gallery |
