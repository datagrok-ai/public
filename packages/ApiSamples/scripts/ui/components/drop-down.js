grok.shell.newView('Drop-Down', [
  ui.ribbonPanel([
    ui.dropDown('Element', () => ui.div('Element')),
    ui.dropDown('List', () => ui.list(['Element 1', 'Element 2', 'Element 3'])),
    ui.dropDown('Menu', () => ui.popupMenu({
      'Foo': () => grok.shell.info('foo'),
      'Bar': () => {}
    }, {show: false}).root)
  ])
]);
