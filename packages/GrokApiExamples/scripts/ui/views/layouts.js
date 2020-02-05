// Layouts save, apply and serialization into user data storage

var STORAGE_NAME = 'layouts-demo';

var view = grok.addTableView(grok.testData('demog', 1000));
view.name = 'layouts demo';

var nameInput = ui.stringInput('', 'Default');
nameInput.captionLabel.style.width = '0';

let save = ui.button('Save', () => {
    grok.dapi.userDataStorage.postValue(STORAGE_NAME, nameInput.value, view.saveLayout().toJson());
});

var load = ui.iconFA('bars', () => {
    grok.dapi.userDataStorage.get(STORAGE_NAME).then((layouts) => {
        if (layouts !== null && Object.keys(layouts).length === 0)
            grok.balloon.info('Storage is empty. Save some layouts to the storage');
        else {
            let menu = Menu.popup();
            for (let layout of Object.keys(layouts)) {
                menu.item(layout, () => {
                    view.loadLayout(ViewLayout.fromJson(layouts[layout]));
                    nameInput.stringValue = layout;
                });
            }
            menu.show();
        }
    });
}, 'Select layout from storage');

let clear = ui.iconFA('trash-alt', () => {
    grok.dapi.userDataStorage.remove(STORAGE_NAME, null);
}, 'Clear storage');

let acc = view.toolboxPage.accordion;
var saveDiv = ui.divH([save, clear]);
saveDiv.style.marginLeft = '12px';
acc.addPane('Layouts demo', () => {
    return ui.div([
        ui.divH([nameInput.root, load]),
        saveDiv
    ], 'pure-form,pure-form-aligned');
}, true, acc.panes[0]);
