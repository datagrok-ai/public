// Creating custom ribbons

let v = grok.newView('toolbox demo', [ui.divText('See custom, view-specific elements on top below the menu.')]);

v.setRibbonPanels([
    [
        ui.iconFA('search', () => grok.balloon.info("clicked")),
        ui.iconFA('plus', () => grok.balloon.info("plus"))
    ],
    [ui.divText('Custom panel')]
]);