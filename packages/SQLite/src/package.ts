import * as DG from 'datagrok-api/dg';

import * as sql from 'sql.js';

//@ts-ignore: there are no types
import initSqlJs from './sql-wasm.js';
import {importSQLiteImpl} from './import';

export const _package = new DG.Package();
let SQL: sql.SqlJsStatic;

//tags: init
export async function sqlJsInit() {
  //TODO: use web worker?
  SQL = await initSqlJs({locateFile: () => _package.webRoot + 'dist/sql-wasm.wasm'});
}

//name: importSQLite
//description: Opens SQLite files
//tags: file-handler
//meta.ext: sqlite
//input: list bytes
//output: list tables
export function importSQLite(bytes: Uint8Array) {
  return importSQLiteImpl(bytes, SQL);
}
