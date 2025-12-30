let showMenu = () => {
  let showBalloon = (item) => grok.shell.info(item);

  DG.Menu.popup()
    .item('Show info', () => grok.shell.info('Info'))
    .separator()
    .items(['First', 'Second'], showBalloon)
    .group('Group')
    .item('Inner item', () => grok.shell.info('item click'))
    .endGroup()
    .show();
};

let text = ui.divText('Click me');
text.onclick = () => showMenu();
grok.shell.newView('Popup Menu', [text]);