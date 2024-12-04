let result = await grok.dapi.admin.getServiceInfos();
let v = grok.shell.newView("list");

v.root.appendChild(
  ui.table(result, (item, idx) => [`${item.key}:`, item.status])
);
