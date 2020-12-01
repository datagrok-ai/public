// https://datagrok.ai/help/visualize/viewers/scatter-plot

let view = grok.shell.addTableView(grok.data.demo.demog());

let plot = view.scatterPlot({
    x: 'height',
    y: 'weight',
    size: 'age',
    color: 'race',
});

plot.onEvent().subscribe((e) => grok.shell.info(e.name));
