  //tags: View
  let host = ui.divV(null, {style: {width: '500px', height: '500px'}});

  let selector = ui.input.choice('View', {items: DG.View.ALL_VIEW_TYPES, value: DG.View.APPS, onValueChanged: (x) => {
    $(host).empty();
    host.appendChild(DG.View.createByType(x).root);
  }});

  ui.dialog().add(ui.divV([selector, host])).show();
  host.appendChild(DG.View.createByType(DG.View.ALL_VIEW_TYPES[0]).root);