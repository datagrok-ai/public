// Proxy a request via Datagrok's server with the same interface as "fetch"

const url = 'https://jsonplaceholder.typicode.com/posts';
let data = { name: 'username', password: 'password' };

// POST
let response = await grok.dapi.fetchProxy(url, {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify(data)
});
grok.shell.info(response.ok)

// GET
response = await grok.dapi.fetchProxy('https://dev.datagrok.ai/demo/demog.csv');
data = await response.text();
console.log(data);
grok.shell.addTableView(DG.DataFrame.fromCsv(data))