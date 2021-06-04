let host = ui.divV([]);

let selector = ui.choiceInput('Type', DG.View.APPS,  DG.View.ALL, (x) => {
  host.appendChild(DG.View.createByType(x, {}));
});

ui.divV([selector, host])