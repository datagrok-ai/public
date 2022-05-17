# PowerGrid

PowerGrid contains popular spreadsheet extensions, such as [sparklines](#sparklines)
and support for [images](#images).

## Sparklines

Sparklines are small charts that are built based on columns with a countable data. They help users to analyze and compare objects in data sets quickly by presenting them in a visual way. In Datagrok sparklines are available in four different reprentations: sparklines, bar charts, pie charts, and radar charts. 

Sparklines come in handy when you want to have a brief look at your data set and form a first impression. Lets check the example:

Imagine that we have a table filled with the patients data.

![](https://i.imgur.com/Ng4iR3v.png)

We want to quickly check the records of patients with a high Body Mass Index. 

![](https://i.imgur.com/2wEOT2m.png)

We add sparklines summary column that includes AGE, HEIGHT, and WEIGHT data. Now we can visually distinguish the records of people with high BMI by checking the sparklines on the right. We are looking for charts where the second and third dots are almost on the same level.

**To add sparklines column, do the following:**

<!--- A place for a gif --->

1. Open the file in a table view
2. Right click on any cell to open a context menu
3. Choose Add->Summary columns->sparklineCellRenderer
___

**Note:** you can also choose pie charts, bar charts, or radar charts as a rendering option.
___

By default sparklines include all the columns with the countable data. **To edit the list of included columns, do the following:**

<!--- A place for a gif --->

1. Click on the header of sparklines column
2. Expand **Renderer** options on the **properties panel**
3. Click on the columns list to open the **Edit menu**
4. Check the columns you want to include and click **OK** to confirm.

___

**Note:** if properties panel is missing, click on **Windows** icon on the side bar, and click on the **Properties** radio button. You can also toggle a **Properties panel** by pressing **F4** hotkey.
___

### Normalize sparklines

By default sparklines are presented in a Normalize mode. In this case the graphs are based on the minimum and maximum values among the records. You may turn off  **Normalize** in a **Renderer** options to set 0 as a starting point for the charts.


## Images

See also:

* [Viewers](https://datagrok.ai/help/visualize/viewers)
* [Develop custom cell renderers](https://datagrok.ai/help/develop/how-to/custom-cell-renderers)
* [Packages](https://datagrok.ai/help/develop/develop#packages)
