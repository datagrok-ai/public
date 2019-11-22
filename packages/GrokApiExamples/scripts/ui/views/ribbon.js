// Creating custom views

let v = View.create();
v.name = 'JS View';
v.root.appendChild(ui.divText('Custom body'));

// ribbons
v.setRibbonPanels([
    [
        ui.iconFA('search', () => gr.balloon.info("clicked")),
        ui.iconFA('plus', () => gr.balloon.info("plus"))
    ],
    [ui.divText('Custom panel')]
]);

// accordion
var acc = ui.accordion();
acc.addPane(('header 1'), () => ui.divText('Dynamic content'));
acc.addPane(('header 2'), () => ui.divText('More dynamic content'));
v.root.appendChild(acc.root);

gr.addView(v);
