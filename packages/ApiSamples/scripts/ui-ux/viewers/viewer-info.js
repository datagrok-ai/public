// Gets and uses properties and info of the viewer

let view = grok.shell.addTableView(grok.data.demo.demog());

let bc = view.barChart();

let aggColSelector = bc.getInfo().aggColSelector;
let splitColSelector = bc.getInfo().splitColSelector;
let canvas = bc.getInfo().canvas;

grok.shell.info('splitColumnName name: ' + bc.props.getProperty('splitColumnName').name);
grok.shell.info('stackColumnName propertyType: ' + bc.props.getProperty('stackColumnName').propertyType);
grok.shell.info('valueColumnName semType:' + bc.props.getProperty('valueColumnName').semType);