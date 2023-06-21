//Layout: Splitter 2 - 3 Column

let windows = grok.shell.windows;
windows.showToolbox = false; //*hide tool box
windows.showContextPanel = false; //*hide context panel

let view = grok.shell.newView('Splitter 2 - 3 Column', [
    ui.splitV([
      ui.splitH([
        ui.panel('Row 1 Column 1'),
        ui.panel('Row 1 Column 2'),
      ]),
      ui.splitH([
        ui.panel('Row 2 Column 1'),
        ui.panel('Row 2 Column 2'),
        ui.panel('Row 2 Column 3'),
      ])
    ])
  ],'layout');
  
  view.box = true;
  
$('.layout .ui-div').css('border','1px dashed var(--grey-2)');