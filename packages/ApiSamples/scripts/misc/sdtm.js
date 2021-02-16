let links = {
  ae: { key: 'USUBJID', start: 'AESTDY', end: 'AEENDY', event: 'AETERM'},
  cm: { key: 'USUBJID', start: 'VISITDY', event: 'CMTRT'},
  ex: { key: 'USUBJID', start: 'EXSTDY', end: 'EXENDY', event: 'EXTRT'},
  lb: { key: 'USUBJID', start: 'LBDY', event: 'LBTEST'}
};

let result = null;

let getTable = function(domain) {
  let info = links[domain];
  let t = grok.shell
    .tableByName(domain)
    .clone(null, Object.keys(info).map(e => info[e]));
  t.columns.addNew('domain', DG.TYPE.STRING).init(domain);
  for (let name in info)
    t.col(info[name]).name = name;
  return t;
}

for (let domain in links) {
  let t = getTable(domain);
  if (result == null)
    result = t;
  else
    result.append(t, true);
}

grok.shell.addTableView(result);