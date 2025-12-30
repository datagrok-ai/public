// Container

let c = ui.div([
  ui.h1('Header'), 
  ui.p('Paragraph text'), 
  'just text', 
  DG.Viewer.scatterPlot(grok.data.demo.demog())
]);

c.style.border = '1px dashed red'; 

grok.shell.newView('Containers', [c]);