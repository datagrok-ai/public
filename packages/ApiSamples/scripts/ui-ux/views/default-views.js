//tags: View
let host = ui.divV(null, {style: { width: '500px', height: '500px'}});

let selector = ui.choiceInput('View', DG.View.APPS,  DG.View.ALL_VIEW_TYPES, (x) => {
  $(host).empty();
  host.appendChild(DG.View.createByType(x).root);
});

ui.dialog().add(ui.divV([selector, host])).show();