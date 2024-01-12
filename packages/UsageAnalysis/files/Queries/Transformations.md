1. Open “Data” section on sidebar and then click on “Databases” item
2. Expand "PostgreSQL > northwind” and then open context menu for “Products” query
3. Click on “Edit…” from context menu
4. Click on “Transformation” tab in Query View (opened in previous step)
5. Click on “Add new column” from actions list and add column which is equal to the “productid” column from the table (${productid})
6. Open Cluster dialog from “ML | Custer…” menu
7. Click OK in Cluster dialog with default parameters
8. Make sure that the actions from steps 6 and 7 are added to the transformation script
9. Click on “Play” button from toolbar
10. Click on “Run query…” from “Actions” tab on Toolbox
11. Click on “Save” button on toolbar
12. Run the query saved in step 11
13. Make sure that the functions added to the transformation script are executed for the returned dataframe