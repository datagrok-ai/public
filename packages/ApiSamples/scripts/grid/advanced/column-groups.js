let table = grok.data.demo.demog();
table.meta.setGroups(
  {
    'Demographics' : { color: '#861616', columns: ['sex', 'race', 'age'] },
    'Biometrics': { color: '#1544A9FF', columns: ['height', 'weight'] }
  })

let view = grok.shell.addTableView(table);
