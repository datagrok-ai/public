<!-- TITLE: Select Tooltip Columns -->
<!-- SUBTITLE: -->

# Select Tooltip Columns

Select columns that should be put on a tooltip when the mouse cursor is over a data point. By default,
all columns are shown, if the table contains less than 21 columns, and no columns are shown, if a table
contains 21 columns or more.

## Viewer Tooltips

One of the unique [exploratory data analysis](exploratory-data-analysis.md)-related 
 features of Datagrok is the ability to interactively visualize
_group of rows_ on an arbitrary viewer when the mouse is over a visual that represents multiple rows.
 
Here are some examples of such visuals:
* Bars in a [histogram](../visualize/viewers/histogram.md)
* Bars in a [bar chart](../visualize/viewers/bar-chart.md)
* Segments and legend items in a [pie chart](../visualize/viewers/pie-chart.md)
* Category labels in a [trellis plot](../visualize/viewers/trellis-plot.md)
* Rectangular areas in a [tree map](../visualize/viewers/tree-map.md)
* Nodes and edges in a [network diagram](../visualize/viewers/network-diagram.md)
* Words in a [word cloud](../visualize/viewers/word-cloud.md)
* Categories and histogram bins in a [filter](../visualize/viewers/filters.md)

You can choose any viewer as a tooltip for the so-called "mouse-over rows". To do that, add a
viewer of interest (such as scatter plot), customize the settings accordingly, and then
select **Tooltip | Use as group tooltip**. After that, once a mouse cursor gets over the segment
that represents multiple rows, a scatter plot containing only mouse-over rows will be shown in a
tooltip. 

"Reset group tooltip" removes group tooltip.

See also
* [Viewers](../visualize/viewers.md)
