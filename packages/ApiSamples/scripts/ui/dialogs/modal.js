// Custom dialogs

ui.dialog({showHeader: false, showFooter: false})
  .add(ui.span(['Full-screen mode with no header or footer; click on the background, or press ESC to close.']))
  .onOK(() => grok.shell.info('OK!'))
  .showModal(true);