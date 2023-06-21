//Layout: Tabs navbar

let windows = grok.shell.windows;
windows.showToolbox = false; //*hide tool box
windows.showContextPanel = false; //*hide context panel

let view = grok.shell.newView('Tabs Navbar', [
    ui.tabControl({
      'First': ui.panel('First tab'),
      'Second': ui.panel('Second tab')
    })
  ],'layout');
  
  view.box = true;
  
$('.layout .ui-div').css('border','1px dashed var(--grey-2)');