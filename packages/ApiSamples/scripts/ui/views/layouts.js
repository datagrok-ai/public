// Layouts save, apply and serialization into user settings storage

const STORAGE_NAME = 'layouts-demo';

let view = grok.shell.addTableView(grok.data.testData('demog', 1000));
view.name = 'layouts demo';

let nameInput = ui.input.string('', {value: 'Default'});
nameInput.captionLabel.style.width = '0';

let save = ui.button('Save', () => {
  grok.userSettings.add(STORAGE_NAME, nameInput.value, view.saveLayout().toJson());
});

let load = ui.iconFA('bars', () => {
  let layouts = grok.userSettings.get(STORAGE_NAME);
  if (layouts !== null && Object.keys(layouts).length === 0)
    grok.shell.info('Storage is empty. Save some layouts to the storage');
  else {
    let menu = DG.Menu.popup();
    for (let layout of Object.keys(layouts)) {
      menu.item(layout, () => {
        view.loadLayout(DG.ViewLayout.fromJson(layouts[layout]));
        nameInput.stringValue = layout;
      });
    }
    menu.show();
  }
}, 'Select layout from storage');

let clear = ui.iconFA('trash-alt', () => {
  grok.userSettings.delete(STORAGE_NAME, null);
}, 'Clear storage');

let acc = view.toolboxPage.accordion;
let saveDiv = ui.divH([save, clear]);
saveDiv.style.marginLeft = '12px';
acc.addPane('Layouts demo', () => {
  return ui.div([
    ui.divH([nameInput.root, load]),
    saveDiv
  ], 'pure-form pure-form-aligned');
}, true, acc.panes[0]);

view.setRibbonPanels([[
  ui.iconFA('broom', () => view.resetLayout()),
]]);
