// Automatic dynamic loading of the package containing the handler

let view = grok.shell.newView();
grok.events.onPackageLoaded.subscribe((p) => grok.shell.info(p.name));

DG.ObjectHandler.forSemType('plate').then((handlers) => {
  console.log(handlers);
  view.append(ui.list(handlers));
});

var plate = {barcode: '123456'};
var handler = DG.ObjectHandler.forEntity(plate);
view.append(handler.name);