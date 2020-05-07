<!-- TITLE: Demo: Chika Mexicana -->
<!-- SUBTITLE: -->

# Demo: Chika Mexicana 

Goal: Based on the restaurant sales data, analyze how time of day affects sales.

1. Get data from zip-archive "chika-mexicana.zip" with .slsx files.
    * Drag-and-drop zip-archive "chika-mexicana.zip" to platform window (or use "Open file" dialog, ```Ctrl+O```)
     
1. Combine all tables with data into one
   * Open *"Tables"* window from *"View"* menu
   * Select all tables (hold down ```Shift``` key and select all 13 tables)
   * In the context menu for selected tables, select **(13 tables | Append)** (or click on **Append** in *"Commands"* tab on Property Panel )
   
1. Delete empty rows from "result" table
   * Click on empty cell in table
   * Press the combination ```Shift + Enter``` from the keyboard
   * Click on **Delete Rows** in *"Actions"* tab on Property Panel (or ```Shift + Delete```)
   
1. Extract hours from  *"Item Created Date"* column to new column 
   * Call the context menu for any cell in *"Item Created Date"* column and select **(Extract | hour)**

1. Add [Line Chart](../../visualize/viewers/line-chart.md) on which the sum of prices by time is displayed
   * Select [Line Chart](../../visualize/viewers/line-chart.md) from *"Add"* menu (or click on [Line Chart](../../visualize/viewers/line-chart.md) pic on Toolbox)
   * Chose column *"hour(Ticket Created Date)"* for X-axis 
   * Double click on the [Line Chart](../../visualize/viewers/line-chart.md), on which the values of *"Price"* column (right selector near Y axis)
   * Select "sum" value for aggregate (left selector near Y axis)
   
1. Delete hour values for which the price value is zero
   * Hold down ```Shift```, select the area on [Line Chart](../../visualize/viewers/line-chart.md) that corresponds to the price value equal to zero
   * Click on **Delete Rows** in *"Actions"* tab on Property Panel (or ```Shift + Delete```)

1. Now you can visually evaluate the peaks of sales in terms of time of day
   * On the [Line Chart](../../visualize/viewers/line-chart.md) you can see that the peaks correspond to the lunch time (12:00 - 13:00 p.m.) and  evening ( 8:00 p.m.)
   
1. Add splitting by days of the week into [Line Chart](../../visualize/viewers/line-chart.md)   
   * Call context menu for any cell in *"Item Created Date"* column and select **(Extract | day of weak)**
   * Chose column *"day of week(Ticket Created Date)"* in *"Split by"* field on [Line Chart](../../visualize/viewers/line-chart.md) 
     (selector on top of [Line Chart](../../visualize/viewers/line-chart.md))
   
1. Now you can visually evaluate the peaks of sales in terms of day of week
      * On [Line Chart](../../visualize/viewers/line-chart.md) you can see that most of the sales occur on Friday, Sunday 
        sales are more "even" throughout the day, and the smallest sales occur on Tuesday and Wednesday

See also:

 * [Line Chart](../../visualize/viewers/line-chart.md) 
 * [Viewers](../../visualize/viewers.md)
