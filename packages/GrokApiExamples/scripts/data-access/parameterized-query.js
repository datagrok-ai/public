// An example of using parameterized query

gr.query('Demo:CoffeeCompany:StoresInState', {'state': 'NY'})
    .then(t => gr.balloon.info('Stores: ' + t.rowCount));
