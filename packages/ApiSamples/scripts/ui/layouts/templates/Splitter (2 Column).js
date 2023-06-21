//Layout: Splitter 2 Column

let windows = grok.shell.windows;
windows.showToolbox = false; //*hide tool box
windows.showContextPanel = false; //*hide context panel

let view = grok.shell.newView('2 Column', [
    ui.splitH([
      ui.panel('Column 1'),
      ui.panel('Column 2')
    ])
  ],'layout');
  
  view.box = true;
  
$('.layout .ui-div').css('border','1px dashed var(--grey-2)');