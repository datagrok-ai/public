// https://datagrok.ai/help/viewers/shape-map

grok.data.loadDataFrame('https://public.datagrok.ai/demo//earnings-by-state.csv').then((t) => {
    let view = grok.shell.addTableView(t);
    view.shapeMap();
})
