// https://datagrok.ai/help/viewers/correlation-plot

let view = grok.addTableView(grok.testData('demog', 5000));

view.corrPlot({
    xs: ['age', 'weight', 'height'],
    ys: ['age', 'weight', 'height'],
});
