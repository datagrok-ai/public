//Layout: 3 Column

let windows = grok.shell.windows;
windows.showToolbox = false; //*hide tool box
windows.showProperties = false; //*hide property panel

let view = grok.shell.newView('3 Column', [
    ui.block50(['Column 1']),
    ui.block(['Column 2']),
    ui.block(['Column 3'])
  ],'layout');
  
  view.box = true;
  
$('.layout div').css('border','1px dashed var(--grey-2)');