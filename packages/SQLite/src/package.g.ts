import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: sqlJsInit
//tags: init
export async function sqlJsInit() {
  return PackageFunctions.sqlJsInit();
}

//name: importSQLite
//description: Opens SQLite files
//tags: file-handler
//input: list bytes 
//output: list<dataframe> result
//meta.ext: sqlite
export function importSQLite(bytes: Uint8Array) {
  return PackageFunctions.importSQLite(bytes);
}
