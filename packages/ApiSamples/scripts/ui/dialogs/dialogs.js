// Creating custom dialogs, and add the help url

ui.dialog({title: 'Windows', helpUrl: '/help/transform/add-new-column.md'})
  .add(ui.span(['People of Earth, your attention, please… ']))
  .onOK(() => {
    grok.shell.info('OK!');
  })
  .show();
