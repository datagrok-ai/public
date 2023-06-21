//Layout: Tabs 2 column

let windows = grok.shell.windows;
windows.showToolbox = false; //*hide tool box
windows.showContextPanel = false; //*hide context panel

let view = grok.shell.newView('Tabs 2 Column', [
    ui.splitH([
      ui.tabControl({
      'First': ui.panel('Tab 1, column 1'),
      'Second': ui.panel('Tab 2, column 1'),
        }).root,
      ui.panel('Column 2')
    ])
  ],'layout');
  
  view.box = true;
  
$('.layout .ui-box').css('border','1px dashed var(--grey-2)');