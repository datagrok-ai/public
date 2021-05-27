// Discovering functions dynamically

let widgetFunctions = DG.Func.find({returnType: 'widget'});

let dialog = ui.dialog().show();
for (let f of widgetFunctions)
  f.apply().then((w) => dialog.add(w));