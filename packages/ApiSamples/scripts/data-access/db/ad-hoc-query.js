//tags: query, ad-hoc query
// All queries execution are saved by default for the history. If you don't want this behaviour use special ad-hoc flag

// Get a connection
let connection = await grok.functions.eval('System:Datagrok')

// Create a query
let q = connection.query('query name', 'select 1 as hello_world');

// Create FuncCall and set adHoc flag
let call = q.prepare();
call.adHoc = true;

// Execute call
await call.call();
grok.shell.addTableView(call.getOutputParamValue());
