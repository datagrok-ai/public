1. Create a project:
    1. with two local files (SPGI, SPGI-linked1)
    2. with two files from the storage without datasync (SPGI, demog from the Home dir)
    3. with two files from the storage with datasync (SPGI, demog from the Home dir)
    4. with two files from the storage, use Join Tables and datasync - doesn't work for now
    5. with two files from the storage, use Link Tables and datasync (SPGI, SPGI-linked1 from the Home dir)
    6. with two query results, use Link Tables and datasync
    7. with a query result and a file from the storage, use Link Tables and datasync
    8. an SDF file with datasync
2. Close all
3. Go to Browse and check the preview of all created projects
4. Share created projects
5. Open created projects
6. Edit projects, save, and reopen:
  * add any viewer and save the original project, close all.
  * open the origianl project, add any viewer and save a copy with the **link**, close all, reopen it. Close all
  * open the origianl project, add any viewer and save a copy with the **clone**, close all, reopen it. Close all
7. Share the new created projects. Check ability to open the shared projects.

***

1. Link tables (SPGI, spgi-linked1, spgi-linked2 from here ).
1. Open SPGI-linked1.csv,  SPGI-linked2.csv
Go to Tables > SPGI-linked1.
1. On the menu ribbon, go to Data > Link Tables.
1. In the dialog, select the SPGI-linked2 table.
2. Set Link type to  'selection to filter' , columns: Sample Name, link column 1, link column 2, link column 3.
> Note: Alternatively, for steps 1-5, run the script:
```
//name: Link Test
//language: javascript
//input: dataframe df1 [Data table 1]
//input: dataframe df2 [Data table 2]
grok.data.linkTables(df1, df2, ['Sample Name', 'link column 1', 'link column 2', 'link column 3'], ['Sample Name', 'link column 1', 'link column 2', 'link column 3'], [DG.SYNC_TYPE.SELECTION_TO_FILTER]);
grok.shell.addTableView(df1);
grok.shell.addTableView(df2); 
```
6. Open SPGI-linked1
2. Add a line chart and set properties as follows:
   * Data: set Filter to ${link column 3}=="v ii" && ${link column 1} <30, 
   * Table to SPGI-linked2, 
   * Row source to filtered, 
   * Overview to link column 1
   * X: Set X to Value1
   * Y: set Y Axis Type to logarithmic, Split to link column 2
6. Select/deselect some/all rows. Verify that data displayed on the line chart is as expected (only filtered rows displayed, data not missing unexpectedly)
9. Save the project.
---
{
  "order": 5
}