// Icons

let view = ui.divV([
  ui.h2('Font Awesome'),
  ui.ribbonPanel([
    ui.iconFA('question'),
    ui.iconFA('info'),
    ui.iconFA('cogs'),
    ui.iconFA('folder-open'),
    ui.iconFA('save')
  ]),
  ui.h2('Special'),
  ui.ribbonPanel([
    ui.icons.settings(() => grok.shell.info('click'), 'Settings'),
    ui.icons.help(() => grok.shell.info('click'), 'Help'),
    ui.icons.close(() => grok.shell.info('click'), 'Close'),
    ui.icons.edit(() => grok.shell.info('click'), 'Edit'),
    ui.icons.save(() => grok.shell.info('click'), 'Save'),
    ui.icons.copy(() => grok.shell.info('click'), 'Copy'),
    ui.icons.add(() => grok.shell.info('click'), 'Add'),
    ui.icons.remove(() => grok.shell.info('click'), 'Remove'),
    ui.icons.delete(() => grok.shell.info('click'), 'Delete'),
    ui.icons.undo(() => grok.shell.info('click'), 'Undo'),
    ui.icons.sync(() => grok.shell.info('click'), 'Sync'),
    ui.icons.info(() => grok.shell.info('click'), 'Info'),
    ui.icons.search(() => grok.shell.info('click'), 'Search'),
    ui.icons.filter(() => grok.shell.info('click'), 'Filter')
  ])
]);

ui.dialog('Icons')
  .add(view)
  .show();