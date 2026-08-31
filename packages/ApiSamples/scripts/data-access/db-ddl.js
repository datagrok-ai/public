// Structured DDL against a database connection via grok.data.db.ddl(...).
// Creates a scratch table, loads rows into it, reads them back, and drops it —
// showing the dryRun plan and the confirmDestructive contract on the way.
// Requires a connection whose provider supports DDL (Postgres, MySQL/MariaDB, MSSQL,
// Oracle) and the DataConnection.CreateTable / DropTable privileges on it.
// Replace CONNECTION with your own.

const CONNECTION = 'Samples:PostgresNorthwind';
const table = `apisamples_ddl_${Date.now()}`;
const ddl = grok.data.db.ddl(CONNECTION).schema('public');

// Every operation can dry-run first: the exact SQL, nothing executed.
const plan = await ddl.createTable(table, [
  {name: 'id', type: DG.TYPE.INT, nullable: false},
  {name: 'name', type: DG.TYPE.STRING},
  {name: 'qty', type: DG.TYPE.INT},
], {primaryKey: ['id']}).dryRun();
grok.shell.info(`Will run: ${plan.statements.join('; ')}`);

await ddl.createTable(table, [
  {name: 'id', type: DG.TYPE.INT, nullable: false},
  {name: 'name', type: DG.TYPE.STRING},
  {name: 'qty', type: DG.TYPE.INT},
], {primaryKey: ['id']}).execute();

try {
  await grok.data.db.table(CONNECTION, `public.${table}`).insert([
    {id: 1, name: 'alpha', qty: 10},
    {id: 2, name: 'beta', qty: 20},
    {id: 3, name: 'gamma', qty: 30},
  ]);

  const readBack = await grok.data.db.query(CONNECTION, `select * from public.${table} order by id`);
  grok.shell.addTableView(readBack);
} finally {
  // Dropping a table with live rows is destructive: an unconfirmed execute() rejects with
  // DdlConfirmationRequiredError carrying the plan (statements + live row counts);
  // pass confirmDestructive to proceed.
  try {
    await ddl.dropTable(table).execute();
  } catch (e) {
    if (e instanceof DG.DdlConfirmationRequiredError)
      grok.shell.info(`Confirmation required: ${e.plan.destructive.map((d) => `${d.code} (${d.count} rows)`).join(', ')}`);
    await ddl.dropTable(table).execute({confirmDestructive: true});
  }
}
