// Set default viewer settings for the specific dataframe

let df = grok.data.demo.demog();
let defaultScatterPlotSettings = {
  xColumnName: 'height',
  yColumnName: 'weight',
  sizeColumnName: 'age',
  colorColumnName: 'race',
}

// now, all newly created viewers attached to df will inherit the specified settings:
df.tags['.Viewer Template: ' + DG.VIEWER.SCATTER_PLOT] = JSON.stringify(defaultScatterPlotSettings);
df.tags['.Viewer Template: ' + DG.VIEWER.HISTOGRAM] = JSON.stringify({valueColumnName: 'height'});

let view = grok.shell.addTableView(df);
view.scatterPlot();
view.histogram();