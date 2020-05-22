// An example of using parameterized query

grok.data.query('Demo:CoffeeCompany:StoresInState', {'state': 'NY'})
    .then(t => grok.shell.balloon.info('Stores: ' + t.rowCount));
