<!-- TITLE: Tests: Join DB Tables -->
<!-- SUBTITLE: -->

# Tests: Data query

Aggregation Editor allows you to get table with aggregated values by columns from DB table.

## Testing scenario

1. Open "Connect to data..." tab on "Welcome" view

1. Expand ```PosgteSQL -> northwind -> Tables```

1. Open *New Aggregation Query* for **products** DB table (from it's context menu)

* Added special "Query View" with control for selection of aggregations

1. Select "*productid"* column for "Columns" field, "*suplplierid"* column for "Rows" field and select *avg(unitprice)*
   for "Measures" field

* Preview shows result table with selected aggregations
* First column correspond to values of "*suplplierid"* column from **products** DB table
* Other columns corresponds to values of "*productid"* column from **products** DB table
* Table values corresponds to average value of "*unitprice"* column for intersection of values of "*suplplierid"*
  and "*productid"* columns

1. Click on "Run query..." action on Toolbox

* Table with result of aggregations from "Aggregation Editor" is added to workspace

1. Return to view with [Source Tree](../../access/access.md#data-sources)

1. Click on "Get All" from context menu of **products** DB table

* Table "*products*" with all values and columns added to workspace

1. Open **Data | Aggregate Rows**

1. Configure fields of [Aggregate rows](../../transform/aggregate-rows.md) dialog as in step 4

1. Execute [Aggregate rows](../../transform/aggregate-rows.md) dialog

* Result of [Aggregate rows](../../transform/aggregate-rows.md) same as result in Aggregation Editor from step 7

1. Repeat previous steps for **Postgres**, **MySql**, **MS SQL**, **MariaDB**, **ORACLE**
   providers

1. Click on "Add result to workspace" from dialog menu (v)

* Result table of query from [Query Builder](../../access/databases/databases.md#join-tables) has been added to workspace

1. Click on "Get All" from context menu of **employees** DB table

* Table "*employees*" with all values and columns added to workspace

1. Repeat previous steps for **Postgres**, **MySql**, **MS SQL**, **MariaDB**, **ORACLE**
   providers

See also:

* [Data Sourse Test](../../access/data-source-test.md)
* [Data Query](../../access/access.md#data-query)
* [Join DB Tables Test](../tests/build-query-test.md)
