// Creating custom ribbons and toolboxes

let v = grok.shell.newView('Custom View', [
    ui.divText('See custom, view-specific elements in the ribbon.')
]);

v.append(ui.button('Close', () => v.close()));

v.setRibbonPanels([
    [
        ui.iconFA('search', () => grok.shell.balloon.info("clicked")),
        ui.iconFA('plus', () => grok.shell.balloon.info("plus"))
    ],
    [ ui.divText('Custom panel') ]
]);

let acc = ui.accordion();
acc.addPane('header 1', () => ui.divText('Dynamic content'));
acc.addPane('header 2', () => ui.divText('More dynamic content'));
v.toolbox = acc.root;