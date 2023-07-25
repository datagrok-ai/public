// https://datagrok.ai/help/visualize/viewers/correlation-plot

let view = grok.shell.addTableView(grok.data.demo.demog());

view.corrPlot({
  xs: ['age', 'weight', 'height'],
  ys: ['age', 'weight', 'height'],
});
