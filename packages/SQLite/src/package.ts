import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import * as sql from 'sql.js';

//@ts-ignore: there are no types
import initSqlJs from './sql-wasm.js';
import {importSQLiteImpl} from './import';
export * from './package.g';
export const _package = new DG.Package();
let SQL: sql.SqlJsStatic;

export class PackageFunctions {
  @grok.decorators.init()
  static async sqlJsInit() {
    //TODO: use web worker?
    SQL = await initSqlJs({locateFile: () => _package.webRoot + 'dist/sql-wasm.wasm'});
  }

  @grok.decorators.fileHandler({
    'ext': 'sqlite',
    'description': 'Opens SQLite files',
  })
  static importSQLite(
    @grok.decorators.param({'type': 'list'}) bytes: Uint8Array): DG.DataFrame[] {
    return importSQLiteImpl(bytes, SQL);
  }
}
