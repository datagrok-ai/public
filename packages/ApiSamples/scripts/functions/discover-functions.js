// Discovering functions dynamically

let widgetFunctions = DG.Func.find({returnType: 'widget'});

let dialog = ui.dialog().show();
for (let f of widgetFunctions)
  dialog.add(ui.render(f));