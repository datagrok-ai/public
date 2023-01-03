<-- TEST SCENARIO NAME should match the tested functionality name. In all other cases create a meaningful name -->
<-- The name of the test scenario should describe which part of the platform it is testing -->

# Manual Test: Upload Project

## Target functionality

Uploading the project to the server in order to be able to share it and store data and entities inside it

## Scenario

Read the table from the top, do not skip any steps until it is intentional, and you know what you are doing.

| Step No. | Action                                          | Expected result                                    | Description                                        |
|:--------:|:------------------------------------------------|:---------------------------------------------------|:---------------------------------------------------|
|    1     | Add the table in Datagrok                       | Dataset opened in Datagrok                         | from anywhere: files, queries, local storage, etc. |
|    2     | Click on "Share" section on Sidebar             | Share section is shown on Toolbox                  |                                                    |
|    3     | Click "Upload" button near "Scratchpad" project | "Upload" dialog is opened                          |                                                    |
|    4     | Enter a name for the new project                |                                                    | e.g. test_project                                  |
|    5     | Enter description for the new project           |                                                    | e.g. test_description                              |
|    6     | Turn on the "“"Presentation mode" switch        |                                                    |                                                    |
|    7     | Click “OK”                                      | "Upload" dialog is closed. "Sharing" dialog opened |                                                    |
|    8     | Click “OK” in “Share dialog”                    | Project has been uploaded to the server            |                                                    |