// You can iterate over the added to the TableView viewers

let parentView = grok.shell.addTableView(grok.data.demo.demog());

parentView.scatterPlot({x: 'height', y: 'weight'});
parentView.histogram({value: 'age'});
parentView.lineChart();

for (let view of parentView.viewers)
  grok.shell.info(JSON.parse(view.getOptions()).type);