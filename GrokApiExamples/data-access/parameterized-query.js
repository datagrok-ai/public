// An example of using parametrized query

gr.query('StoresInState', {'state': 'NY'})
    .then(t => gr.balloon.info('Stores: ' + t.rowCount));
