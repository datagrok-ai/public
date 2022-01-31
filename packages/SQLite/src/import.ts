
import * as DG from 'datagrok-api/dg';

import * as sql from 'sql.js';

export function importSQLiteImpl(bytes: Uint8Array, SQL: sql.SqlJsStatic): DG.DataFrame[] {
  const db = new SQL.Database(bytes);
  const tableList = db.exec('SELECT `name` FROM `sqlite_master` WHERE type=\'table\';')[0].values;
  const result = [];

  for (const tableName of tableList) {
    const tableObjects = [];
    const statement = db.prepare(`SELECT * FROM ${tableName[0]};`);

    while (statement.step())
      tableObjects.push(statement.getAsObject());

    statement.reset();

    const currentDf = DG.DataFrame.fromObjects(tableObjects);
    currentDf!.name = tableName[0]!.toString();
    result.push(currentDf!);
  }
  db.close();

  return result;
}
