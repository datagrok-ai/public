<-- TEST SCENARIO NAME should match the tested functionality name. In all other cases create a meaningful name -->
<-- The name of the test scenario should describe which part of the platform it is testing -->

# Manual Test: Import Swagger

## Target functionality

OpenAPI, also known as swagger, is a popular format that describes the structure of the server APIs so that machines can read the document and use the service.

Datagrok integrates with OpenAPI really well. Once a swagger file is imported (you can simply drag-and-drop a yaml or json file into the app), its content gets translated to data connections, queries, and functions. 

## Scenario

Read the table from the top, do not skip any steps until it is intentional, and you know what you are doing.


| Step No. | Action                                                                                  | Expected result               | Description |
|:--------:|:----------------------------------------------------------------------------------------|:---------------------------|:------------|
|    1     | Find the "openweathermap.yaml" file in the Swaggers package                             |                            | (https://github.com/datagrok-ai/public/tree/master/packages/Swaggers/swaggers) | 
|    2     | Drag-and-drop  "openweathermap.yaml" file from your to "Datagrok"                       | New Web connection created |             |
|    3     | Click on "Data" section on sidebar and then click on "Webservices"                      | "Webservices" view is open |             |
|    4     | Make sure the new connection is created and follows the description in the Swagger File |                            | compare the structure in swagger file and the created connection in the platform |
