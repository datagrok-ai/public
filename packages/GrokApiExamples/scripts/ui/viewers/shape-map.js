// https://datagrok.ai/help/viewers/shape-map

grok.loadDataFrame('https://public.datagrok.ai/demo//earnings-by-state.csv').then((t) => {
    let view = grok.addTableView(t);
    view.shapeMap();
})
