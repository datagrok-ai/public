//tags: DataQuery
//help-url: https://datagrok.ai/help/access/parameterized-queries
// An example of using parameterized query

grok.data.query('Demo:CoffeeCompany:StoresInState', {'state': 'NY'})
  .then(t => grok.shell.info('Stores: ' + t.rowCount));
