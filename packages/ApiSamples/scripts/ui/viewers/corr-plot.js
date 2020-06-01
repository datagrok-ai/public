// https://datagrok.ai/help/viewers/correlation-plot

let view = grok.shell.addTableView(grok.data.testData('demog', 5000));

view.corrPlot({
    xs: ['age', 'weight', 'height'],
    ys: ['age', 'weight', 'height'],
});
