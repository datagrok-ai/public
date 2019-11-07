// Creating custom views

let v = gr.newView('accordion');

v.root.appendChild(ui.h1('UI Toolkit'));

// accordion
var acc = ui.accordion();

acc.addPane(('buttons'), () => ui.div([
    ui.button('REGULAR'),
    ui.bigButton('BIG')]));

acc.addPane(('headers'), () => ui.divV([
    ui.h1('Header 1'),
    ui.h2('Header 2'),
    ui.h3('Header 3')]));

acc.addPane(('tables'), () => ui.tableFromMap({
    user: gr.user,
    project: gr.project,
    time: new Date(),
}));

acc.addPane(('rendering'), () => ui.span([
    ui.h1('Rendering'),
    'Currently, ', gr.user, ' has the following project open: ', gr.project
]));

acc.addPane(('dialogs'), () => ui.button('OPEN', () => {
    ui.dialog('Vogon Announcement')
      .add(ui.h1(''))
      .add(ui.span(['People of Earth, your attention, please… ']))
      .onOK(() => { gr.balloon.info('OK!'); })
      .show();
}));

acc.addPane(('menus'), () => ui.button('SHOW', () => {
    ui.popupMenu({
        'About': () => gr.balloon.info('About'),
        'File': {
            'Open': () => gr.balloon.info('Open'),
            'Close': () => gr.balloon.info('Close'),
        }
    });
}));

v.root.appendChild(acc.root);
