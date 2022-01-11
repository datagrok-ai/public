import * as DG from 'datagrok-api/dg';
//@ts-ignore: there are no types
import initSqlJs from './sql-wasm.js';
import * as sql from 'sql.js';

export const _package = new DG.Package();

//name: importSQLite
//description: Opens SQLite files
//tags: file-handler
//meta.ext: sqlite
//input: list bytes
//output: list tables
export async function importSQLite(bytes: number[]) {
  //TODO: use web worker?
  const SQL: sql.SqlJsStatic = await initSqlJs({locateFile: () => _package.webRoot + 'dist/sql-wasm.wasm'});
  const db = new SQL.Database(Uint8Array.from(bytes));
  const tableList = db.exec('SELECT `name` FROM `sqlite_master` WHERE type=\'table\';')[0].values;
  const result = [];

  for (const tableName of tableList) {
    const tableObjects = [];
    const statement = db.prepare(`SELECT * FROM ${tableName[0]};`);

    while (statement.step()) {
      tableObjects.push(statement.getAsObject());
    }
    statement.reset();

    const currentDf = DG.DataFrame.fromObjects(tableObjects);
    currentDf!.name = tableName[0]!.toString();
    result.push(currentDf);
  }
  db.close();

  return result;
}
