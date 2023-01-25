<-- TEST SCENARIO NAME should match the tested functionality name. In all other cases create a meaningful name -->
<-- The name of the test scenario should describe which part of the platform it is testing -->

# Manual Test: Connection Schema

## Target functionality

Schema browser visualizes all tables with all columns at once, giving you a high-level overview of the database. 

## Scenario

Read the table from the top, do not skip any steps until it is intentional, and you know what you are doing.
This scenario needs to be reproduced for different DBs (Oracle, MariaDB, MySQL, MS SQL)

| Step No. | Action                                                                         | Expected result                                      | Description |
|:--------:|:-------------------------------------------------------------------------------|:-----------------------------------------------------|-------------|
|    1     | Open "Data" section on sidebar and then click on “Databases” item              | "Databases" view is open                             |             | 
|    2     | Expand the "PostgresDart" and call the context menu for "northwond" connection |                                                      |             |
|    3     | Click on "Browse schema" from context menu                                     | View with schema is open                             |             |
|    4     | Expand "View" and "Tables" sections on Toolbox                                 | Tab content is displayed correctly                   |             |
|    5     | Call context menu for any table on schema                                      | Context actions correspond to working with DB tables |             |
