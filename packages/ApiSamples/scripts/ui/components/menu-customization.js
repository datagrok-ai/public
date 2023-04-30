// Menu item customization. For top menu, see menu.js

let host = ui.divText('Click for context menu');
let menu = DG.Menu.popup()
  .item('Foo', () => grok.shell.info('Foo'), null, {
    description: 'Explore Foo',
    shortcut: 'Ctrl+E' })
  .group('Single choice')
  .items(['Apple', 'Banana', 'Peach'], (s) => grok.shell.info(s), { toString: s => s, isChecked: (s) => s === 'Banana'})
  .endGroup()
  .item('Bar', () => grok.shell.info('Bar'), null, {order: 0, check: true})
  .separator()
  .items(['Apple', 'Banana', 'Peach'], (s) => grok.shell.info(s), {radioGroup: 'Fruits'})
  .show();
