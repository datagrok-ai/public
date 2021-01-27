// Proxy a request via Datagrok's server with the same interface as "fetch"

const url = 'https://jsonplaceholder.typicode.com/posts';
const data = { name: 'username', password: 'password' };

// POST
grok.dapi.fetchProxy(url, {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify(data)
}).then(response => grok.shell.info(response.ok));

// GET
grok.dapi.fetchProxy('https://public.datagrok.ai/demo/demog.csv')
  .then(response => response.text())
  .then(data => grok.shell.addTableView(DG.DataFrame.fromCsv(data)));
