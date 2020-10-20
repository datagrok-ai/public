// Creating custom views

let v = grok.shell.newView('accordion');

v.root.appendChild(ui.h1('UI Toolkit'));

// accordion
let acc = ui.accordion();

acc.addPane('buttons', () => ui.div([
    ui.button('REGULAR'),
    ui.bigButton('BIG')]));

acc.addPane('headers', () => ui.divV([
    ui.h1('Header 1'),
    ui.h2('Header 2'),
    ui.h3('Header 3')]));

acc.addPane('tables', () => ui.tableFromMap({
    user: grok.shell.user,
    project: grok.shell.project,
    time: new Date(),
}));

acc.addPane('rendering', () => ui.span([
    ui.h1('Rendering'),
    'Currently, ', grok.shell.user, ' has the following project open: ', grok.shell.project
]));

acc.addPane('dialogs', () => ui.button('OPEN', () => {
    ui.dialog('Vogon Announcement')
      .add(ui.h1(''))
      .add(ui.span(['People of Earth, your attention, please… ']))
      .onOK(() => { grok.shell.info('OK!'); })
      .show();
}));

acc.addPane('menus', () => ui.button('SHOW', () => {
    ui.popupMenu({
        'About': () => grok.shell.info('About'),
        'File': {
            'Open': () => grok.shell.info('Open'),
            'Close': () => grok.shell.info('Close'),
        }
    });
}));

acc.addPane('text area', () => ui.textArea('multi\nline\ntext'));

acc.addPane('range ', () => ui.rangeSlider(0, 10, 2, 5).root);

//v.root.appendChild(acc.root);
v.append(acc.root);