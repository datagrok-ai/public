// Creating custom views

let v = DG.View.create();
v.name = 'JS View';
v.root.appendChild(ui.divText('Custom body'));

// ribbons
v.setRibbonPanels([
  [
    ui.iconFA('search', () => grok.shell.info('clicked')),
    ui.iconFA('plus', () => grok.shell.info('plus'))
  ],
  [ui.divText('Custom panel')]
]);

// ribbon menu
v.ribbonMenu = DG.Menu.create()
  .group('Custom')
  .item('Foo!', () => grok.shell.info('Foo clicked'));

// accordion
let acc = ui.accordion();
acc.addPane('header 1', () => ui.divText('Dynamic content'));
acc.addPane('header 2', () => ui.divText('More dynamic content'));
v.root.appendChild(acc.root);

grok.shell.addView(v);
