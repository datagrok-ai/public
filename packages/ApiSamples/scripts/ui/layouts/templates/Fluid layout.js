//Layout: Fluid layout

let windows = grok.shell.windows;
windows.showToolbox = false; //*hide tool box
windows.showContextPanel = false; //*hide context panel

let view = grok.shell.newView('Fluid layout', [
    ui.block25(['block 1']),
    ui.block75(['block 2']),
    ui.block50(['block 3']),
    ui.block50(['block 4']),
    ui.block25(['block 5']),
    ui.block25(['block 6']),
    ui.block25(['block 7']),
    ui.block25(['block 8']),
  ],'layout');
  
$('.layout div').css('border','1px dashed var(--grey-2)');
$('.layout div').css('padding','10px');