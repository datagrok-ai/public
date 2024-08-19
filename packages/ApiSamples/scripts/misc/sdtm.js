let links = {
  ae: {key: 'USUBJID', start: 'AESTDY', end: 'AEENDY', event: 'AETERM'},
  sv: {key: 'USUBJID', start: 'VISITDY', event: 'VISIT'},
  lb: {key: 'USUBJID', start: 'LBDY', event: 'LBTEST'}
};

let result = null;

let getTable = async function(domain) {
  let info = links[domain];
  let df = await grok.data.files.openTable(`System:AppData/ClinicalCase/${domain}.csv`);
  let t = df.clone(null, Object.keys(info).map(e => info[e]));
  t.columns.addNew('domain', DG.TYPE.STRING).init(domain);
  for (let name in info)
    t.col(info[name]).name = name;
  return t;
};

for (let domain in links) {
  let t = await getTable(domain);
  if (result == null)
    result = t;
  else
    result.append(t, true);
}

grok.shell.addTableView(result);