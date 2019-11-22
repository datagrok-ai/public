// https://datagrok.ai/help/viewers/correlation-plot

let view = gr.addTableView(gr.testData('demog', 5000));

view.corrPlot({
    xs: ['age', 'weight', 'height'],
    ys: ['age', 'weight', 'height'],
});
