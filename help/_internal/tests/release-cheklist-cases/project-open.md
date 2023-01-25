<-- TEST SCENARIO NAME should match the tested functionality name. In all other cases create a meaningful name -->
<-- The name of the test scenario should describe which part of the platform it is testing -->

# Manual Test: Open Project

## Target functionality

Opening a project previously uploaded to the server and deserializing its contents

## Scenario

Read the table from the top, do not skip any steps until it is intentional, and you know what you are doing.

| Step No. | Action                                                                          | Expected result                        | Description                              |
|:--------:|:--------------------------------------------------------------------------------|:---------------------------------------|:-----------------------------------------|
|    1     | Find the project you upload in [Upload project test-case](./upload-project.md)  | Project selected as current object     |                                          | 
|    2     | Check all attributes added to the project (sharing, description, name, picture) |                                        |                                          |
|    3     | Double-click on it project                                                      | Project opened                         | that it should open in presentation mode |
|    4     | Press F7 to exit from presentation mode                                         | Project is displayed in standard mode  |                                          |
|    5     | Open "Sharing" section on sidebar to see open project                           | Open project is present on the Toolbox |                                          |
|    6     | Call context menu on sidebar and Close All                                      |                                        | Nothing left open in the platform        |