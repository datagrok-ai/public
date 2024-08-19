// https://datagrok.ai/help/visualize/viewers/radar

let view = grok.shell.addTableView(grok.data.demo.demog());

view.addViewer('radar', {
  split: 'race',
  value: 'age',
  valueAggrType: 'avg'
});