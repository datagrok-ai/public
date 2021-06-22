//Layout: 2 Column (75/25)

let windows = grok.shell.windows;
windows.showToolbox = false; //*hide tool box
windows.showProperties = false; //*hide property panel

let view = grok.shell.newView('2 Column (75/25)', [
    ui.block75(['Column 1 (75%)']),
    ui.block25(['Column 2 (25%)'])
  ],'layout');
  
  view.box = true;
  
$('.layout div').css('border','1px dashed var(--grey-2)');