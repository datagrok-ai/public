<!-- TITLE: Tests: Audit -->
<!-- SUBTITLE: -->

# Tests: audit

Audit system is intended to store all meaningful user activity for future analysis.

All of the changes made to the different [entities](../overview/objects.md) are tracked and can be analyzed using
activity sections. All events are joined to user session. Changes made to the entities are connected with corresponding
entities. Audit records can be posted both from client and server.

The implementation is based on the [grid](../visualize/viewers/grid.md), so many of the grid's features apply.

## Data connection

1. Open *"Connect to data"* view on *"Welcome"* tab

2. Create new connection empty *"audit-test-connection"* for **PostgreSQL** source

    * New connection *"audit-test-connection"* created

3. Open "Edit connection" dialog for *"audit-test-connection"* connection

4. Enter "Test_Server" to "Server" field and save changed

5. Open *"audit-test-connection"* on [Property Panel](../overview/navigation.md#properties)

6. Expand "Activity" tab on [Property Panel](../overview/navigation.md#properties)

    * There is history from previous steps
    * "created" and "edited" actions arranged in chronological order with indication of user

7. Delete *"audit-test-connection"*

    * (expected result is presented in "Results" section of this document)

## Data query

1. Open *"Connect to data"* view on *"Welcome"* tab

2. Create new empty query *"audit-test-query"* for *northwind* connection in **PostgreSQL** source

    * New query *"audit-test-query"* created

3. Open "Query View" for edit *"audit-test-query"* query

4. Add query text on "Query" tab (*)

5. Save changes

6. Open "Query View" for edit *"audit-test-query"* query

7. Go to "Transformations" tab

8. Add transformation for query (for example, delete *"productid"* column)

9. Save transformations

10. Run *"audit-test-query"* query

    * Query completed

11. Open *"audit-test-query"* on [Property Panel](../overview/navigation.md#properties)

12. Expand "Activity" tab on [Property Panel](../overview/navigation.md#properties)

    * There is history from previous steps
    * "created", "edited", "transformations edited", "started" actions arranged in chronological order with indication
      of user

13. Delete *"audit-test-query"*

    * (expected result is presented in "Results" section of this document)

## Data job

1. Open **Admin | Data Jobs**

2. Click "Add new job" action

    * "New Job" view was open

3. Add steps for job ```TestData(demog, 1000) -> DeleteColumns(\["subj", "study"\])```

4. Save jobs with name *"audit-test-job"*

5. Run *"audit-test-job"* [data job](../access/data-job.md)

    * Job completed

6. Open "Edit Job" view for *"audit-test-job"*

7. Add step ```RenameColumn("age", "New_Age")``` after ```DeleteColumns(\["subj", "study"\])``` step

8. Save job changes

9. Open *"audit-test-job"* on [Property Panel](../overview/navigation.md#properties)

10. Expand "Activity" tab on [Property Panel](../overview/navigation.md#properties)

    * There is history from previous steps
    * "created", "started", "edited" actions arranged in chronological order with indication of user

11. Delete *"audit-test-job"*

    * (expected result is presented in "Results" section of this document)

## Scripts

1. Open "New Script" from **Tools | Scripting | New Script | R**

    * "New Script" view was open
    * R-script template contained in script text

2. Change script name to *"audit-test-script"* (`````#name: audit-test-script`````)

3. Click "Save" button on Toolbar

    * *"audit-test-script"* saved

4. Open "Scrip Browser" from **Tools | Scripting | Browse Scripts**

5. Open *"audit-test-script"* for editing

6. Add script text (**)

7. Click "Save" button on Toolbar

    * *"audit-test-script"* saved

8. Run *"audit-test-script"* with default parameters

    * Script completed

9. Open *"audit-test-script"* on [Property Panel](../overview/navigation.md#properties)

10. Expand "Activity" tab on [Property Panel](../overview/navigation.md#properties)

    * There is history from previous steps
    * "created", "edited" and "started" actions arranged in chronological order with indication of user

11. Delete *"audit-test-script"*

    * (expected result is presented in "Results" section of this document)

## Predictive models

1. Open *demog* table

2. Use **Tools | Data science | Missing Values Imputation** for impute data to null rows

3. Open "Train Model" view from **Tools | Predictive Modeling | Train**

    * "Predictive model" view was open

4. Enter *"audit-test-model"* to "Name" field

5. Select *"SEX"* column for "Predict" field

6. Select *"AGE"*,*"HEIGHT"*,*"WEIGHT"* columns for "Features" field

7. Leave other parameters with default values and train new model

8. Open "Model Browser" from **Tools | Predictive Modeling | Browse Models**

9. Open "Edit model" dialog for *"audit-test-model"*

10. Enter "Test_description" to "Description" field and add #test tag

11. Save model changes

12. Apply *"audit-test-model"* model to "demog" table

    * Applying complete

13. Open *"audit-test-model"* on [Property Panel](../overview/navigation.md#properties)

14. Expand "Activity" tab on [Property Panel](../overview/navigation.md#properties)

    * There is history from previous steps
    * "created", "edited" and "ran" actions arranged in chronological order with indication of user

15. Delete *"audit-test-model"*

    * (expected result is presented in "Results" section of this document)

## Notebooks

1. Open "demog" table

2. Open "demog" table in Jupyter Notebook

    * Notebook with "demog" table opened in new view tab

3. Save Jupyter Notebook with "demog" table

4. Open "Notebook Browser" from **Tools | Notebooks | Browse Notebooks**

5. Rename created in previous steps notebook to *"audit-test-notebook"*

6. Apply *"audit-test-notebook"* to "demog" table

    * Applying complete

7. Open *"audit-test-notebook"* on [Property Panel](../overview/navigation.md#properties)

8. Expand "Activity" tab on [Property Panel](../overview/navigation.md#properties)

    * There is history from previous steps
    * "created", "edited" and "ran" actions arranged in chronological order with indication of user

9. Delete *"audit-test-notebook"*

    * (expected result is presented in "Results" section of this document)

## Projects

1. Open "demog" table

2. Upload current workspace as project with name *"audit-test-project"*

3. Close all by **File | Close all**

4. Open *"Datasets..."* on "Welcome" view

5. Share *"audit-test-project"* project to "All users" from it's context menu

6. Open *"audit-test-project"*

7. Delete *"AGE"* column from "demog" table

8. Save project changes (re-upload)

9. Close all by **File | Close all**

10. Open *"Datasets..."* on "Welcome" view

11. Open *"audit-test-project"* on [Property Panel](../overview/navigation.md#properties)

12. Expand "Activity" tab on [Property Panel](../overview/navigation.md#properties)

    * There is history from previous steps
    * "created", "shared", "edited" and "open" actions arranged in chronological order with indication of user

13. Delete *"audit-test-project"*

    * (expected result is presented in "Results" section of this document)

## Main menu, dialogs and functions

1. Open "demog" table from local file using **File | Open**

2. Use **Tools | Data science | Missing Values Imputation** for impute data to null rows

3. Add [Scatter-Plot](../visualize/viewers/scatter-plot.md) viewer from **Add** menu

4. Open [Find and Replace](../transform/find-and-replace.md) from menu **Edit**

5. Replace `M` value of *"SEX"* to `"Men"`

    * `M` value of *"SEX"* ​​replaced by `"Men"`

6. Open **View | Columns**

    * *"Columns"* is open

7. Open **View | Columns**

    * *"Columns"* is open

8. Open **Tools| Console**

    * *"Console"* is open

9. Call function ```eq("test","test")``` in [Console](../overview/navigation.md#console)

    * Functions returned ```true```
    * result added to Variables

10. Open **Tools | Data science | Multivariate Analysis (PLS)**

11. Select *"AGE"* column for "Predict" field and *"HEIGHT"* and *"WEIGHT"* columns for "Features"
    and execute dialog

12. Open **Help | Wiki**

    * "Help" view was open

13. Open **Help | Functions**

    * "Functions" view was open

14. Call *DeleteColumns" function

    * Parameter input dialog was open

15. Select *"HEIGHT"* and *"WEIGHT"* columns for deleting and click ```OK```

    * *"HEIGHT"* and *"WEIGHT"* columns removed from "demog" table

## Results

1. Open *"Connect to data"* view on *"Welcome"* tab

2. Run *"Log Actions Summary by Hours"* query from "datagrok" connection in **PostgreSQL** source for last 2 hours

    * Query completed
    * Audit table added to workspace

3. Use filter for "login" column with current user

4. Inspect data presented in the table

    * There are events:
      `query-created`, `query-edited`, `query-start`, `query-transformations-edited`, `query-transformations-edited`
      , `connection-created`, `connection-edited`, `job-created`, `job-edited`, `job-transformations-edited`
      , `job-start`
      , `script-created`, `script-edited`, `script-start`, `predictive-model-created`, `predictive-model-edited`
      , `predictive-model-start`, `notebook-created`, `notebook-edited`, `notebook-start`, `project-created`
      , `project-edited`, `project-opened`
      in *"name"* column for actions with entities in previous sections of this document
    * There are events:
      `query-deleted`, `connection-deleted`, `job-deleted`, `script-deleted`, `predictive-model-deleted`
      , `notebook-deleted`
      , `project-deleted`
      in *"name"* column corresponding to entity deletion actions in previous sections of this document
    * *"name"* column contains entries `dialog-ok` corresponding clicking `OK` in platform dialogs in previous section
    * Table contains functions calls records whose names in *"name"* field and `function` in *"
      source"*
    * Table contains events corresponding to clicks on Main Menu and its items

See also:

* [Audit](audit.md)

(*):

```sql
SELECT * FROM Products
```

(**):

```
#name: t-test
#description: Welch's t-test
#help-url: https://en.wikipedia.org/wiki/welch%27s_t-test
#language: r
#tags: demo
#sample: TSLA.csv
#input: dataframe data [Input data table]
#input: column x {type:numerical} [X axis column name]
#input: column y {type:numerical} [Y axis column name]
#output: double pvalue [P-value of t-statistics]

    require(stats)

    ttest = t.test(data[[x]], data[[y]])
    pValue = ttest$p.value
```
