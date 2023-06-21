//Layout: Splitter 1 Column 2 Row

let windows = grok.shell.windows;
windows.showToolbox = false; //*hide tool box
windows.showContextPanel = false; //*hide context panel

let view = grok.shell.newView('1 Column - 2 Rows', [
    ui.splitH([
      ui.panel('Column 2'),
      ui.splitV([
        ui.panel('Row 1'),
        ui.panel('Row 2')
      ])
    ])
  ],'layout');
  
  view.box = true;
  
$('.layout .ui-div').css('border','1px dashed var(--grey-2)');