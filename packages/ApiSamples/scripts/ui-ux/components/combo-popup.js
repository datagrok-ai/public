// Drop-down button

grok.shell.newView('tables', [
  ui.div([ui.comboPopupItems('FooBar', {
    'Foo': () => grok.shell.info('Foo!'),
    'Bar': () => grok.shell.info('Bar!'),
  })]),
  ui.div([ui.comboPopup('SpamEggs', ['Spam', 'Eggs'], (item) => grok.shell.info(item), (item) => ui.span([`Item: ${item}`]))])
]);
