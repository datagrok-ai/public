<-- TEST SCENARIO NAME should match the tested functionality name. In all other cases create a meaningful name -->
<-- The name of the test scenario should describe which part of the platform it is testing -->

# Manual Test: Adding connections

## Target functionality

Data connection is used for accessing data in a particular data source. Connection parameters depend on the data source. Typically, you would need to provide server name and login credentials.

## Scenario

Read the table from the top, do not skip any steps until it is intentional, and you know what you are doing.
This scenario needs to be reproduced for different DBs (Oracle, MariaDB, MySQL, MS SQL)

| Step No. | Action                                                            | Expected result                                      | Description                         |
|:--------:|:------------------------------------------------------------------|:-----------------------------------------------------|:------------------------------------|
|    1     | Open “Data” section on sidebar and then click on "Databases" item | Source tree view is open                             |                                     | 
|    2     | Open context menu for "Postgres” and click on “Add new connection | "Add connection" dialog is open                      |                                     |
|    3     | Enter "test_postgres" to "Name" field                             |                                                      |                                     |
|    4     | Click on "Test" button                                            | Error ballon is shown                                | correctly handled error is expected |
|    5     | Fill other fields with following data*                            |                                                      |                                     |
|    6     | Click on "Test" button                                            | Balloon with a message about a successful connection |                                     |
|    7     | Click OK button                                                   | Data connection has been created                     |                                     |

/* Server = db.datagrok.ai:54322
DB = northwind
login = datagrok
password = datagrok
