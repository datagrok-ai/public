// Drop-down button

grok.shell.newView('tables', [
    ui.comboPopupItems('FooBar', {
        'Foo': () => grok.shell.info('Foo!'),
        'Bar': () => grok.shell.info('Bar!'),
    })
]);