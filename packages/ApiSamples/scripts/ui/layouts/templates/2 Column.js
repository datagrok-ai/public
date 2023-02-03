//Layout: 2 Column

let windows = grok.shell.windows;
windows.showToolbox = false; //*hide tool box
windows.showContextPanel = false; //*hide context panel

let view = grok.shell.newView('2 Column', [
    ui.block25(['Column 1 (25%)']),
    ui.block75(['Column 2 (75%)']),
  ],'layout');
  
  view.box = true;
  
$('.layout div').css('border','1px dashed var(--grey-2)');