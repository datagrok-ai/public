// Inspecting view contents

let view = grok.shell.addTableView(grok.data.demo.demog());

view.scatterPlot({x: 'height', y: 'weight'});
view.histogram({value: 'age'});
view.lineChart();

for (let viewer of view.viewers)
  grok.shell.info(viewer.type);