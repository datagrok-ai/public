//tags: DataQuery
//help-url: https://datagrok.ai/help/access/parameterized-queries
// An example of using parameterized query

grok.data.query('System:CoffeeCompany:StoresInState', {'state': 'NY'})
  .then(t => grok.shell.info('Stores: ' + t.rowCount));
