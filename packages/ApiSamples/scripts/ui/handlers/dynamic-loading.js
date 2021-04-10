// Automatic dynamic loading of the package containing the handler

DG.ObjectHandler
  .forSemType('plate')
  .then((handlers) => { grok.shell.newView().append(ui.list(handlers)); });
