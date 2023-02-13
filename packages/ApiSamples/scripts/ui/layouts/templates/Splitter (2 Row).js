//Layout: Splitter 2 Rows

let windows = grok.shell.windows;
windows.showToolbox = false; //*hide tool box
windows.showContextPanel = false; //*hide context panel

let view = grok.shell.newView('2 Rows', [
    ui.splitH([
      ui.splitV([
        ui.panel('Row 1'),
        ui.panel('Row 2')
      ])
    ])
  ],'layout');
  
  view.box = true;
  
$('.layout .ui-div').css('border','1px dashed var(--grey-2)');