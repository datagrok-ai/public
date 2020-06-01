// https://datagrok.ai/help/viewers/scatter-plot

let view = grok.shell.addTableView(grok.data.testData('demog', 5000));

view.scatterPlot({
    x: 'height',
    y: 'weight',
    size: 'age',
    color: 'race',
});
