let v = gr.newView('demo: context menu');

let showMenu = () => {
    let showBalloon = (item) => gr.balloon.info(item);

    Menu.popup()
        .item('Show info', () => gr.balloon.info('Info'))
        .separator()
        .items(['First', 'Second'], showBalloon)
        .show();
};

let text = ui.divText('Clickable');
text.addEventListener("click", showMenu);
