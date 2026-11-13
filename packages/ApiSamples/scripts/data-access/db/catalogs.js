// The catalogs (databases) of a database server, and the schemas of one of them.
// Postgres, MS SQL, Snowflake, Databricks and some others have catalogs; a data source
// tells whether it does. Replace the connection with your own.

const connection = await grok.functions.eval('System:Datagrok');

if (DG.ConnectionDataSource.byType(connection.dataSource).supportCatalogs) {
  // One request, names only
  const catalogs = await grok.dapi.connections.getCatalogs(connection);
  grok.shell.info(`${catalogs.length} catalog(s): ${catalogs.join(', ')}`);

  // The schemas of a catalog: here the one the connection itself points to
  const own = connection.parameters['db'];
  const schemas = await grok.dapi.connections.getSchemas(connection, own);
  grok.shell.info(`Schemas of ${own}: ${schemas.join(', ')}`);
}
else
  grok.shell.info(`${connection.dataSource} has no catalogs`);
