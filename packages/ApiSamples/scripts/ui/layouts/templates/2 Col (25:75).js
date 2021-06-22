//Layout: 2 Column (25/75)

let windows = grok.shell.windows;
windows.showToolbox = false; //*hide tool box
windows.showProperties = false; //*hide property panel

let view = grok.shell.newView('2 Column (25/75)', [
  ui.block25(['Column 1 (25%)']),
  ui.block75(['Column 2 (75%)'])
],'layout');

view.box = true;

$('.layout div').css('border','1px dashed var(--grey-2)');