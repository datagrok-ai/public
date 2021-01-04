let view = grok.shell.addTableView(grok.data.demo.demog());
let sp = view.scatterPlot();

var properties = sp.props.getProperties()
  .map((p) => p.propertyType + ' ' + p.name + ': ' + p.description + ' ' + p.columnFilter)
  .join('<br>')

grok.shell.info(properties);