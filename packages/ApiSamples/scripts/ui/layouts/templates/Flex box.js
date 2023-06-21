//Layout: Flex box

let windows = grok.shell.windows;
windows.showToolbox = false; //*hide tool box
windows.showContextPanel = false; //*hide context panel

let view = grok.shell.newView('Flex box', [
    ui.divV([
      ui.divH([ui.divText('item 1'), ui.divText('item 2')]),
      ui.divH([ui.divText('item 3'), ui.divText('item 4'),ui.divText('item 5')]),
      ui.divH([ui.divText('item 6'), ui.divText('item 7'),ui.divText('item 8'), ui.divText('item 9')]),
    ]),
  ],'layout');
  
$('.layout div:not(.d4-flex-col)').css('border','1px dashed var(--grey-2)');
$('.layout div:not(.d4-flex-col)').css('padding','10px');