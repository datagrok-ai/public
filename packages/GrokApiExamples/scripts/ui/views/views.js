// Creating custom ribbons

let v = gr.newView('toolbox demo', [ui.divText('See custom, view-specific elements on top below the menu.')]);

v.setRibbonPanels([
    [
        ui.iconFA('search', () => gr.balloon.info("clicked")),
        ui.iconFA('plus', () => gr.balloon.info("plus"))
    ],
    [ui.divText('Custom panel')]
]);