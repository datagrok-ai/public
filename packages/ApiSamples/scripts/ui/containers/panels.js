// Panel

let panel = ui.panel([
  ui.h1('Header'), 
  ui.p('Paragraph text'), 
  'just text', 
  DG.Viewer.scatterPlot(grok.data.demo.demog())
]);

panel.style.border = '1px dashed red'; 

grok.shell.newView('Panel', [panel]);