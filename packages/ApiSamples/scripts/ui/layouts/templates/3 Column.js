//Layout: 3 Column

let windows = grok.shell.windows;
windows.showToolbox = false; //*hide tool box
windows.showProperties = false; //*hide property panel

let view = grok.shell.newView('3 Column', [
    ui.block(['Column 1 (33%)']),
    ui.block(['Column 2 (33%)']),
    ui.block(['Column 2 (33%)'])
  ],'layout');
  
  view.box = true;
  
$('.layout div').css('border','1px dashed var(--grey-2)');