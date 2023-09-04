let t = grok.data.demo.demog();
t.col('weight').tags['.default-axis-type'] = 'logarithmic';
grok.shell.addTableView(t).scatterPlot();
