// https://datagrok.ai/help/viewers/scatter-plot

let view = grok.shell.addTableView(grok.data.demo.demog());

let plot = view.scatterPlot({
    x: 'height',
    y: 'weight',
    size: 'age',
    color: 'race',
});

plot.setOptions({
    showRegressionLine: true,
    markerType: 'square'
});