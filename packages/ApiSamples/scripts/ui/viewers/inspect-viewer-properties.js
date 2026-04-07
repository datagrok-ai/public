let view = grok.shell.addTableView(grok.data.demo.demog());
let sp = view.scatterPlot();

// Use viewer.props.<property> to get or set a property
sp.props.xColumnName = 'race';
grok.shell.info(sp.props.xColumnName);

let descriptions = sp.props.getProperties()
  .map((p) => p.propertyType + ' ' + p.name + ': ' + p.description + ' ' + p.columnFilter)
  .join('<br>');

grok.shell.info(descriptions);