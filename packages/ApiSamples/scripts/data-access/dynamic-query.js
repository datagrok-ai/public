// Creating a dynamic query and getting data

// Get a connection. Also, you can use grok.dapi.connection to read it
let connection = await grok.functions.eval('System:Datagrok')

// Create a query
let q = connection.query('query name', 'select 1 as hello_world');

// Get data
let data = await q.apply();

// Add to workspace
grok.shell.addTableView(data);
