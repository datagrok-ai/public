// An example of using parameterized query

grok.data.query('Demo:CoffeeCompany:StoresInState', {'state': 'NY'})
    .then(t => grok.shell.info('Stores: ' + t.rowCount));
