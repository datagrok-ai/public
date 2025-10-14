// https://datagrok.ai/help/visualize/viewers/scatter-plot


let eventId = 'd4-scatterplot-point-click';
let view = grok.shell.addTableView(grok.data.demo.demog());

let plot = view.scatterPlot({
  x: 'height',
  y: 'weight',
  size: 'age',
  color: 'race',
});

plot.onPointClicked.subscribe((e) => grok.shell.info(`${e.args.rowId} click`));
plot.onPointDoubleClicked.subscribe((e) => grok.shell.info(`${e.args.rowId} double click`));
plot.onZoomed.subscribe((e) => grok.shell.info(`zoomed`));

