let v = grok.newView('demo: context menu');

let showMenu = () => {
    let showBalloon = (item) => grok.balloon.info(item);

    Menu.popup()
        .item('Show info', () => grok.balloon.info('Info'))
        .separator()
        .items(['First', 'Second'], showBalloon)
        .show();
};

let text = ui.divText('Clickable');
text.addEventListener("click", showMenu);
