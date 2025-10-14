let df = grok.data.demo.demog();
df.tags[DG.TAGS.TOOLTIP] = 'weight\nheight\nsex\nage';  // all viewers inherit it
let view = grok.shell.addTableView(df);
view.scatterPlot({title: 'Custom tooltip', rowTooltip: 'age\nsex', showTooltip: 'show custom tooltip', dataValues: 'Do not add'}); // override
view.scatterPlot({title: 'Derived from table'});
