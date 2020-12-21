grok.dapi.admin.getServiceInfos().then(result => {
  let v = grok.shell.newView('list');

  v.root.appendChild(ui.table(result, (item, idx) =>
    [`${item.key}:`, item.status]));
});