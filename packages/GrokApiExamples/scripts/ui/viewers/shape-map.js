// https://datagrok.ai/help/viewers/shape-map

grok.data.loadTable('https://public.datagrok.ai/demo//earnings-by-state.csv').then((t) => {
    let view = grok.shell.addTableView(t);
    view.shapeMap();
})
