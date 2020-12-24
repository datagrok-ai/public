// Creating custom dialogs


var dlg = ui.dialog('Windows');
dlg
  .add(ui.span(['People of Earth, your attention, please… ']))
  .onCancel(() => grok.shell.info('Cancel'))
  .onOK(() => grok.shell.info('OK!'))
  .addButton('HIDE', () => $(dlg.getButton('CANCEL'), 0, null).hide())
  .addContextAction('My Action', () => grok.shell.info('Action!'))
  .show({x: 300, y: 300});