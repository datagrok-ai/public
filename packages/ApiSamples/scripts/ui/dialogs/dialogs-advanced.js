// Creating custom dialogs

let dlg = ui.dialog({title: 'Boo', helpUrl: 'https://wikipedia.org'});
dlg
  .add(ui.span(['People of Earth, your attention, please… ']))
  .onCancel(() => grok.shell.info('Cancel'))
  .onOK(() => grok.shell.info('OK!'))
  .addButton('HIDE', () => $(dlg.getButton('CANCEL'), 0, null).hide())
  .addContextAction('My Action', () => grok.shell.info('Action!'))
  .show({x: 300, y: 300});

// Set the dialog help url
dlg.helpUrl = '/help/transform/add-new-column.md';