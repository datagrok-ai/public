<!-- TITLE: Tests: Transformation Editor -->
<!-- SUBTITLE: -->

# Tests: Transformation Editor

[Recipe Editor](../features/recipe-editor) lets you define data transformations. [Transformations](../features/recipe-editor.md) can be used as a post-run step for data
queries and data jobs, or can be executed separately. You can also share [Transformations](../features/recipe-editor.md) with other users.

## Testing scenario

1. Open [Transformations Editor](../features/recipe-editor.md) for *stores in @state* query. (*PostgreSQL* provider, *starbucks* connection)
   * New view [Transformations Editor](../features/recipe-editor.md) is added for *stores in @state* query 
   * First step of  transformation is selection of query parameters.
   * Table preview is empty because no parameters are defined.
    
1. Add "RenameColumn" action to second step
   * "RenameColumn" action added
   * Warning about not completed previous step
     
1. Determine value of parameter ```State = "NY"``` in first step of transformation
   * Table preview shows data according to parameter entered. (```State = "NY"```)
   * Warning for incomplete step not shown
   * Second step ("RenameColumn") is active for editing
   * Warning is displayed stating that "RenameColumn" action fields are not defined.
   * First step is highlighted in green

1. Click on ```►``` icon for action "RenameColumn"
   * Action fields edit area is expanded
   
1. Select column *"city"* in "column" field and enter value ```town``` in "New Name" field 
   * "RenameColumn" action step is highlighted in green
   * *"city"* column name was changed to *"town"* in table preview
   
1. Switch to list of actions categories on [Toolbox](../../overview/toolbox.md)

1. Choose "Enrich" actions category

1. Switch back to list of actions on [Toolbox](../../overview/toolbox.md)
   * Only "Enrich" category actions are shown
   
1. Choose "AddNewColumn" action for next step in transformation
   
1. Enter to "Expression" field: ```${street_address}, ${city}, ${state}``` and ```auto``` to "Name" field. ("Type" field leave default)
   * "AddNewColumn" action step is highlighted in green
   * New column *"Address"* (concatenation of *"street_address"*, *"town*" and *"state"*  columns) is added to table preview
   
1. Enter "AddressToCoordinates" in the search box on the toolbox. 
   * Search shows matching results as you type each character
   * After full input, results show "AddressToCoordinates" action "test"
   
1. Add "AddressToCoordinates" action to [Transformation](../../transform/recipe-editor.md) 
   * "address" field in "AddressToCoordinates" action step automatically filled by *"Address*" column
   * "AddressToCoordinates" action step is highlighted in green
   * New columns *"latitude"* and *"longitude"* are added to table preview
   
1. Click on first step in [Transformation](../features/recipe-editor.md)
   * Table preview was transferred to original
   
1. Delete "RenameColumn" action step 
   * "AddNewColumn" action step not performed, highlighted in red and warning with error is shown
   * "AddressToCoordinates" action step not performed, highlighted in red and warning about previous step failed is shown
   
1. Edit "Expression" field: ```${street_address}, ${city}, ${state}``` for "AddNewColumn" action step
   * "AddNewColumn" action step is executed and highlighted in green
   * "AddressToCoordinates" action step is executed and highlighted in green (*)
   
1. Close [Transformations Editor](../features/recipe-editor).md view, go to "Connect to data" view

1. Select *stores in @state* query. (*PostgreSQL* provider, *starbucks* connection) in [Property Panel](../features/property-panel.md)
   * "transformation" tab of [Property Panel](../features/property-panel.md) shows the script created earlier (**)

1. Run *stores in @state* query with value ```NY``` for @state parameter
   * Table contains *"Address"*, *"latitude"* and *"longitude"* columns described in query transformation
   
(\*) While the algorithm is being executed (transformation step is being executed), step is highlighted in blue.   
(\**) Script: ```AddNewColumn("StoresInState", "${street_address}, ${city}, ${state}", "Address", "auto", false)
                 AddressToCoordinates("StoresInState", "Address")```

 
See also:
  * [Recipe Editor](../../transform/features/recipe-editor.md)
  * [Console](../../overview/console.md)
  * [Console Test](../../overview/console-test.md)
