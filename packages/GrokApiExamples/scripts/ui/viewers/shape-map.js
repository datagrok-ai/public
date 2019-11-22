// https://datagrok.ai/help/viewers/shape-map

gr.loadDataFrame("/demo/earnings-by-state.csv").then((t) => {
    let view = gr.addTableView(t);
    view.shapeMap();
})
