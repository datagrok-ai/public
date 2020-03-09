<!-- TITLE: Tests: Audit -->
<!-- SUBTITLE: -->

# Tests: Audit

Audit system is intended to store all meaningful user activity for future analysis. 

All of the changes made to the different [entities](../entities/entities.md) are tracked and can be analyzed using activity sections.
All events are joined to user session. Changes made to the entities are connected with corresponding entities.
Audit records can be posted both from client and server.

The implementation is based on the [grid](../viewers/grid.md), so many of the grid's features apply.

## Data Connection

1. Open *"Connect to data"* view on *"Welcome"* tab

1. Create new connection empty *"audit-test-connection"* for **PostgreSQL** source 
   * New connection *"audit-test-connection"* created
   
1. Open "Edit connection" dialog for *"audit-test-connection"* connection

1. Enter "Test_Server" to "Server" field and save changed
 
1. Open *"audit-test-connection"* on [Property Panel](../features/property-panel.md)    

1. Expand "Activity" tab on [Property Panel](../features/property-panel.md)
   * There is history from previous steps
   * "created" and "edited" actions arranged in chronological order with indication of user
   
1. Delete *"audit-test-connection"*
   * (expected result is presented in "Results" section of this document)


## Data Query

1. Open *"Connect to data"* view on *"Welcome"* tab

1. Create new empty query *"audit-test-query"* for *northwind* connection in **PostgreSQL** source 
   * New query *"audit-test-query"* created
   
1. Open "Query View" for edit *"audit-test-query"* query

1. Add query text on "Query" tab (*)

1. Save changes

1. Open "Query View" for edit *"audit-test-query"* query

1. Go to "Transformations" tab

1. Add transformation for query (for example, delete *"productid"* column)

1. Save transformations

1. Run *"audit-test-query"* query
   * Query completed
   
1. Open *"audit-test-query"* on [Property Panel](../features/property-panel.md)    

1. Expand "Activity" tab on [Property Panel](../features/property-panel.md)
   * There is history from previous steps
   * "created", "edited", "transformations edited", "started" actions arranged in chronological order with indication of user
   
1. Delete *"audit-test-query"*
   * (expected result is presented in "Results" section of this document)
   
## Data Job

1. Open **Admin | Data Jobs**

1. Click "Add new job" action
   * "New Job" view was open
   
1. Add steps for job ```TestData(demog, 1000) -> DeleteColumns(\["subj", "study"\])```

1. Save jobs with name *"audit-test-job"*

1. Run *"audit-test-job"* [data job](../access/data-job.md)
   * Job completed
   
1. Open "Edit Job" view for *"audit-test-job"* 

1. Add step ```RenameColumn("age", "New_Age")``` after ```DeleteColumns(\["subj", "study"\])``` step

1. Save job changes

1. Open *"audit-test-job"* on [Property Panel](../features/property-panel.md)    

1. Expand "Activity" tab on [Property Panel](../features/property-panel.md)
   * There is history from previous steps
   * "created", "started", "edited" actions arranged in chronological order with indication of user
   
1. Delete *"audit-test-job"*
   * (expected result is presented in "Results" section of this document)
   
## Scripts

1. Open "New Script" from **Tools | Scripting | New Script | R**
   * "New Script" view was open
   * R-script template contained in script text

1. Change script name to *"audit-test-script"* (`````#name: audit-test-script`````)
   
1. Click "Save" button on Toolbar
   * *"audit-test-script"* saved
   
1. Open "Scrip Browser" from **Tools | Scripting | Browse Scripts**

1. Open *"audit-test-script"* for editing

1. Add script text (**)

1. Click "Save" button on Toolbar
   * *"audit-test-script"* saved
   
1. Run *"audit-test-script"* with default parameters   
   *  Script completed
   
1. Open *"audit-test-script"* on [Property Panel](../features/property-panel.md)    

1. Expand "Activity" tab on [Property Panel](../features/property-panel.md)
   * There is history from previous steps
   * "created", "edited" and "started" actions arranged in chronological order with indication of user
   
1. Delete *"audit-test-script"*
   * (expected result is presented in "Results" section of this document)
   
## Predictive Models

1. Open *demog* table

1. Use **Tools | Data science | Missing Values Imputation** for impute data to null rows

1. Open "Train Model" view from **Tools | Predictive Modeling | Train**
   * "Predictive model" view was open

1. Enter *"audit-test-model"* to "Name" field

1. Select *"SEX"* column for "Predict" field

1. Select *"AGE"*,*"HEIGHT"*,*"WEIGHT"* columns for "Features" field

1. Leave other parameters with default values and train new model

1. Open "Model Browser" from **Tools | Predictive Modeling | Browse Models**

1. Open "Edit model" dialog for *"audit-test-model"*

1. Enter "Test_description" to "Description" field and add #test tag

1. Save model changes

1. Apply *"audit-test-model"* model to "demog" table
   * Applying complete
   
1. Open *"audit-test-model"* on [Property Panel](../features/property-panel.md)    

1. Expand "Activity" tab on [Property Panel](../features/property-panel.md)
   * There is history from previous steps
   * "created", "edited" and "ran" actions arranged in chronological order with indication of user
   
1. Delete *"audit-test-model"*
   * (expected result is presented in "Results" section of this document)
   
   
## Notebooks

1. Open "demog" table

1. Open "demog" table in Jupyter Notebook
   * Notebook with "demog" table opened in new view tab
   
1. Save Jupyter Notebook with "demog" table

1. Open "Notebook Browser" from **Tools | Notebooks | Browse Notebooks**

1. Rename created in previous steps notebook to *"audit-test-notebook"*

1. Apply *"audit-test-notebook"* to "demog" table
   * Applying complete
   
1. Open *"audit-test-notebook"* on [Property Panel](../features/property-panel.md)    

1. Expand "Activity" tab on [Property Panel](../features/property-panel.md)
   * There is history from previous steps
   * "created", "edited" and "ran" actions arranged in chronological order with indication of user
   
1. Delete *"audit-test-notebook"*
   * (expected result is presented in "Results" section of this document)
   
## Projects

1. Open "demog" table
   
1. Upload current workspace as project with name *"audit-test-project"*

1. Close all by **File | Close all**

1. Open *"Datasets..."* on "Welcome" view

1. Share *"audit-test-project"* project to "All users" from it's context menu

1. Open *"audit-test-project"*

1. Delete *"AGE"* column from "demog" table

1. Save project changes (re-upload)

1. Close all by **File | Close all**

1. Open *"Datasets..."* on "Welcome" view

1. Open *"audit-test-project"* on [Property Panel](../features/property-panel.md)    

1. Expand "Activity" tab on [Property Panel](../features/property-panel.md)
   * There is history from previous steps
   * "created", "shared", "edited" and "open" actions arranged in chronological order with indication of user
   
1. Delete *"audit-test-project"*
   * (expected result is presented in "Results" section of this document)

## Main Menu, Dialogs and Functions

1. Open "demog" table from local file using **File | Open**

1. Use **Tools | Data science | Missing Values Imputation** for impute data to null rows

1. Add [Scatter-Plot](../viewers/scatter-plot.md) viewer from **Add** menu

1. Open [Find and Replace](../transform/find-and-replace.md) from menu **Edit**

1. Replace ```M``` value of *"SEX"* to ```"Men"```
   * ```M``` value of *"SEX"* ​​replaced by ```"Men"```
   
1. Open **View | Columns** 
   * *"Columns"* is open

1. Open **View | Columns** 
   * *"Columns"* is open 
   
1. Open **Tools| Console** 
   * *"Console"* is open 
   
1. Call function ```eq("test","test")``` in [Console](../features/console.md)
   * Functions returned ```true```
   * result added to Variables
   
1. Open **Tools | Data science | Multivariate Analysis (PLS)**

1. Select *"AGE"* column for "Predict" field and *"HEIGHT"* and *"WEIGHT"* columns for "Features" and execute dialog

1. Open **Help | Wiki**
   * "Help" view was open
   
1. Open **Help | Functions**
   * "Functions" view was open
   
1. Call *DeleteColumns" function
   * Parameter input dialog was open
   
1. Select *"HEIGHT"* and *"WEIGHT"* columns for deleting and click ```OK```
   * *"HEIGHT"* and *"WEIGHT"* columns removed from "demog" table

## Results

1. Open *"Connect to data"* view on *"Welcome"* tab

1. Run *"Log Actions Summary by Hours"* query from "datagrok" connection in **PostgreSQL** source for last 2 hours
   * Query completed
   * Audit table added to workspace

1. Use filter for "login" column with current user

1. Inspect data presented in the table
   * There are events:
       ```query-created, query-edited, query-start, query-transformations-edited, query-transformations-edited, connection-created, connection-edited, job-created, job-edited, job-transformations-edited, job-start, script-created, script-edited, script-start, predictive-model-created, predictive-model-edited, predictive-model-start, notebook-created, notebook-edited, notebook-start, project-created, project-edited, project-opened```
      in *"name"* column for actions with entities in previous sections of this document
   * There are events:
        ```query-deleted, connection-deleted, job-deleted, script-deleted, predictive-model-deleted, notebook-deleted, project-deletet```
        in *"name"* column corresponding to entity deletion actions in previous sections of this document
   * *"name"* column contains entries ```dialog-ok``` corresponding clicking ```OK``` in  platform dialogs in previous section 
   * Table contains functions calls records whose names in *"name"* field and ```function``` in *"source"*
   * Table contains events corresponding to clicks on Main Menu and its items

See also:
 * [Audit](audit.md)
 
 
 (*):
 ```SELECT * FROM Products```
 
 (**):
 ```#name: t-test
    #description: Welch's t-test
    #help-url: https://en.wikipedia.org/wiki/Welch%27s_t-test
    #language: r
    #tags: demo
    #sample: TSLA.csv
    #input: dataframe data [Input data table]
    #input: column x {type:numerical} [X axis column name]
    #input: column y {type:numerical} [Y axis column name]
    #output: double pValue [P-value of t-statistics]
    
    require(stats)
    
    ttest = t.test(data[[x]], data[[y]])
    pValue = ttest$p.value
```