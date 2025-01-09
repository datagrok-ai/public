### Filtering 

1. Open SPGI
1. Add scatterplot, histogram, line chart, bar chart, pie chart, trellis plot, box plot
1. Set a legend for each viewer to Stereo Category
1. Open the Filter Panel and use:
   * Structure filters
   * Numerical filters
   * Categorical filter - **check the legend and viewer state after each type of filtering**
1. Save the layout
1. Apply the layout - check the legend
1. Reset the filter
1. For each viewer set in-viewer Filter to `${Stereo Category} in ["R_ONE", "S_UNKN"]` - **check the legend and viewer state**
1. Go to the Filter Panel and apply some filters - **check the legend and viewer state**
1. Save the layout
1. Apply the layout - check the legend
1. Set different values of Row Source for the viewers - **check the legend and viewer state**

Close All

1. Open SPGI
1. Add a scatterplot, histogram, line chart, bar chart, pie chart, trellis plot, box plot
1. Set a legend for each viewer
1. For the scatterplot:

   * Apply zooming to filter out some data - **check legends on the other viewers**
   * Double-click the scatterplot to reset the filtering
1. Go to the bar chart:

   * Set Data > On Click to Filter
   * Click bars and check legends on the other viewers
   * To reset the filtering, click the white space on the bar chart
1. Go to the pie chart:

   * Set Misc > On Click to Filter
   * Click the sections on the pie chart and check legends on the other viewers
   * To reset the filtering, click the white space on the pie chart
1. Close All

#### Scatterplot

1. Open SPGI
1. Add a scatterplot: Chemical space X/Y on the X and Y axes, respectively
1. Color: Primary scaffold name
1. Marker: Stereo category
1. Open the Filter Panel, add 'Primary scaffold name', and deselect some values -  ** check if data on the scatterplot is filtered, too; filtered out categories are no longer shown in the scatter plot legend**
1. Click 'R_ONE' in the scatter plot legend - previously filtered out points and legend categories should not re-appear in the scatterplot. The scatterplot should be further filtered by the selected marker category.**

#### Bar chart

1. Open SPGI
1. Add a bar chart, set:

   * Value (count) to CAST Idea ID 
   * Category to  Stereo Category 
   * Stack to Primary scaffold name
1. Uncheck **Value > Include nulls**
1. Open the Filter Panel and add Primary scaffold name
1. Filter by Primary scaffold name - **all displayed categories should be shown in the bar chart legend**

Close All

1. Open SPGI
1. Add a scatterplot
1. Set Marker to Stereo category
1. Set in-viewer Filter to `${Stereo Category} in ["R_ONE", "S_UNKN"]` - **check the legend**
1. Add another scatterplot
1. Set in-viewer Filter to `${Stereo Category} in ["R_ONE", "S_UNKN"]`
1. Set Marker to Stereo category - **check the legend**
1. Save the layout
1. Apply the layout - check the legend
---
{
  "order": 7,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}