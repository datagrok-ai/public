let df = grok.data.demo.demog();
df.tags[DG.TAGS.TOOLTIP] = 'weight\nheight\nsex\nage';  // all viewers inherit it
let view = grok.shell.addTableView(df);
view.scatterPlot({title: 'Custom tooltip', rowTooltip: 'age\nsex'}); // override
view.scatterPlot({title: 'Derived from table'});
