// https://datagrok.ai/help/visualize/viewers/globe

grok.data.getDemoTable('geo/world_pop_1990.csv').then((t) => {
    grok.shell.addTableView(t).addViewer('Globe');
});
