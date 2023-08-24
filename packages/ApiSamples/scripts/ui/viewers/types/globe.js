// https://datagrok.ai/help/visualize/viewers/globe

grok.data.getDemoTable('geo/earthquakes.csv').then((t) => {
  grok.shell.addTableView(t).addViewer('Globe');
});
