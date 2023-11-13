// Menu item customization. For top menu, see menu.js

const host = ui.divText('Click for context menu');
DG.Menu.popup()
  .item('Explore', () => grok.shell.info('Foo'), null, {
    description: 'Explore Foo',
    shortcut: 'Ctrl+E'})
  .group('Single choice')
  .items(['Apple', 'Banana', 'Peach'], (s) => grok.shell.info(s),
    {radioGroup: 'Fruits', isChecked: (s) => s === 'Banana'})
  .endGroup()
  .group('Multi choice')
  .items(['Apple', 'Banana', 'Peach'], (s) => grok.shell.info(s), {isChecked: (s) => s === 'Banana'})
  .endGroup()
  .separator()
  .item('Touch', () => grok.shell.info('Foo'), null, {isEnabled: () => 'Can\'t touch this!'})
  .bind(host)
  .show();