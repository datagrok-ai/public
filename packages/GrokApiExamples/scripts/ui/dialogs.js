// Creating custom dialogs

ui.dialog('Vogon Announcement')
  .add(ui.h1(''))
  .add(ui.span(['People of Earth, your attention, please… ']))
  .onOK(() => { gr.balloon.info('OK!'); })
  .show();