let v = grok.shell.newView('demo: context menu');

let showMenu = () => {
    let showBalloon = (item) => grok.shell.balloon.info(item);

    ui.Menu.popup()
        .item('Show info', () => grok.shell.balloon.info('Info'))
        .separator()
        .items(['First', 'Second'], showBalloon)
        .show();
};

let text = ui.divText('Clickable');
v.append(text);
text.addEventListener("click", showMenu);
