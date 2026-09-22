const connection = await grok.dapi.connections.first();

// Get the connection's parameters
grok.shell.info(`${connection.name}: ${JSON.stringify(connection.parameters)}`);

// A database, or a file share / secret store (S3, Files, ...)
grok.shell.info(`${connection.dataSource}, database: ${connection.isDatabase}`);

// Find the connection's queries
const queries = await grok.dapi.queries.filter(`connection.id = "${connection.id}"`).list();
grok.shell.info(`Found queries: ${queries.length}`);
