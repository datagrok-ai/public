// https://datagrok.ai/help/visualize/viewers/bar-chart

let view = grok.shell.addTableView(grok.data.demo.demog());

view.barChart({
  split: 'race',
  value: 'age',
  valueAggrType: 'avg'
});
