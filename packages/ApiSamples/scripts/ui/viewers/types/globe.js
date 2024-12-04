// https://datagrok.ai/help/visualize/viewers/globe

let t = await grok.data.getDemoTable('geo/earthquakes.csv');
grok.shell.addTableView(t).addViewer('Globe');