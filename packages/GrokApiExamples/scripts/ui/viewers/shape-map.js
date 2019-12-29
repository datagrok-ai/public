// https://datagrok.ai/help/viewers/shape-map

grok.loadDataFrame("/demo/earnings-by-state.csv").then((t) => {
    let view = grok.addTableView(t);
    view.shapeMap();
})
