let t = grok.data.demo.demog();
t.col('weight').setTag(DG.Tags.DefaultAxisType, DG.AxisType.logarithmic);
grok.shell.addTableView(t).scatterPlot();
