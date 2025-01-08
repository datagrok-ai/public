### Legend and viewer consistency when axes are changed

#### Scatterplot

1. Open SPGI
1. Add two columns with formulas:

   * col1: `if(${Stereo Category}!='S_UNKN', null, ${Average Mass})`
   * col2: `if(${Stereo Category}=='S_UNKN', null, ${Average Mass})`
1. Add a scatter plot and set col1 as X axis
1. Set color by Stereo Category
1. Change X axis to col2
1. Check for all filter and zoom options - legend categories should be updated according to the new data on the plot

#### Line chart

1. Open SPGI
1. Add a line chart
1. Configure two Y columns and select Data > Multi Axis
1. Add a split column - a legend appears on the plot
1. Change one of the Y columns using the in-plot column selector - **the legend should be updated when the column changes**
---
{
  "order": 6,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}