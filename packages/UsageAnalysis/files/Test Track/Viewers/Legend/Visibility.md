### Color match, visibility, adjustment, category selection

1. Open SPGI
1. Add scatterplot, histogram, line chart, bar chart, pie chart, trellis plot, box plot
1. Set a legend for each viewer
1. Check:

   * The legend is visible
   * The colors on the legend match the colors on the viewer
   * Changing the column for Color (Split, Stack) setting changes the legend, respectively
   * You can adjust the legend size
   * The legend position is changed / the legend is hidden when you resize the viewer
   * You can select multiple categories via CTRL+click
1. Save the layout
1. Apply the saved layout - changes should be saved
1. Save the project
1. Close All
1. Open the saved project - changes should be saved

Close All

1. Open SPGI
Add a  line chart
1. Set Split to Primary series names, Series - check that different categories have different colors (within ten categories, colors should be more diverse)
1. On the Context Panel > Data, select Multi Axis - check that each line has its own set of categories in the legend
1. Save the layout
1. Apply the saved layout - changes should be saved
1. Save the project
1. Close All
1. Open the saved project - changes should be saved

Close All 

1. Open SPGI
1. Add a scatterplot
1. Set Color and Marker to Series - the legend must be combined
1. Check the color picker visibility and change some colors
1. Save the layout
1. Apply the saved layout - changes should be saved
1. Save the project
1. Close All
1. Open the saved project - changes should be saved
1. Go to the plot and use ‘+’ to set a new 1. Color value - the Add NewColumn dialog opens
1. Set: 

    * `if(${Stereo Category}=='S_UNKN', null, ${Average Mass})` - the linear legend should appear
    * `if(${Stereo Category}=='S_UNKN', null, ${Series})` - the categorical legend should appear
---
{
  "order": 1,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}