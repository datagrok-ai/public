const connection = await grok.dapi.connections.first();

// Get the connection's parameters
grok.shell.info(`${connection.name}: ${JSON.stringify(connection.parameters)}`);

// Find the connection's queries
const queries = await grok.dapi.queries.filter(`connection.id = "${connection.id}"`).list();
grok.shell.info(`Found queries: ${queries.length}`);
