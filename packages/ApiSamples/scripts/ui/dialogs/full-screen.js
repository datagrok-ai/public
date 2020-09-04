// Creating custom dialogs

ui.dialog('Windows')
  .add(ui.span(['People of Earth, your attention, please… ']))
  .onOK(() => { grok.shell.info('OK!'); })
  .showModal(true);
