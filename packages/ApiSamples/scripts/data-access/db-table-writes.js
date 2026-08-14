// Structured writes to a database table via grok.data.db.table(...).
// Inserts a small DataFrame into a scratch table and reads it back.
// Requires a Postgres connection whose provider supports writes and the
// DataConnection.AddRows permission on it (insert). Replace CONNECTION with your own.

const CONNECTION = 'Samples:PostgresNorthwind';
const table = `public.apisamples_writes_${Date.now()}`;

// Create a scratch table (raw SQL runs through the query path).
await grok.data.db.query(CONNECTION, `create table ${table} (id int primary key, name text, created timestamp)`);

try {
  const dbTable = grok.data.db.table(CONNECTION, table);

  // Insert a typed DataFrame (bulk).
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'id', [1, 2, 3]),
    DG.Column.fromStrings('name', ['alpha', 'beta', 'gamma']),
    DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'created', [new Date(), new Date(), new Date()]),
  ]);
  const dfResult = await dbTable.insert(df);
  grok.shell.info(`Inserted from DataFrame: ${dfResult.affectedRows}`);

  // Insert an array of row objects. Each column is typed by scanning all of its values:
  // a Date becomes a real datetime column and a null is preserved as SQL NULL.
  const objResult = await dbTable.insert([
    {id: 4, name: 'delta', created: new Date()},
    {id: 5, name: null, created: null},
  ]);
  grok.shell.info(`Inserted from object[]: ${objResult.affectedRows}`);

  const readBack = await grok.data.db.query(CONNECTION, `select * from ${table} order by id`);
  grok.shell.addTableView(readBack);
} finally {
  await grok.data.db.query(CONNECTION, `drop table if exists ${table}`);
}
