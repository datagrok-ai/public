// Structured writes to a database table via grok.data.db.table(...).
// Inserts a small DataFrame into a scratch table and reads it back.
// Requires a Postgres connection whose provider supports writes and
// DataConnection.Write permission on it. Replace CONNECTION with your own.

const CONNECTION = 'Samples:PostgresNorthwind';
const table = `public.apisamples_writes_${Date.now()}`;

// Create a scratch table (raw SQL runs through the query path).
await grok.data.db.query(CONNECTION, `create table ${table} (id int primary key, name text)`);

try {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'id', [1, 2, 3]),
    DG.Column.fromStrings('name', ['alpha', 'beta', 'gamma']),
  ]);

  const result = await grok.data.db.table(CONNECTION, table).insert(df);
  grok.shell.info(`Inserted rows: ${result.affectedRows}`);

  const readBack = await grok.data.db.query(CONNECTION, `select * from ${table} order by id`);
  grok.shell.addTableView(readBack);
} finally {
  await grok.data.db.query(CONNECTION, `drop table if exists ${table}`);
}
