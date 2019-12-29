// An example of using parameterized query

grok.query('Demo:CoffeeCompany:StoresInState', {'state': 'NY'})
    .then(t => grok.balloon.info('Stores: ' + t.rowCount));
