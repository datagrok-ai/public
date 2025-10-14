// Automatic dynamic loading of the package containing the handler

let view = grok.shell.newView();
grok.events.onPackageLoaded.subscribe((p) => console.log(p));

let handlers = await DG.ObjectHandler.forSemType('plate');

console.log(handlers);
view.append(ui.list(handlers));

var plate = {barcode: 'plate'};
var handler = DG.ObjectHandler.forEntity(plate);
view.append(handler.name);