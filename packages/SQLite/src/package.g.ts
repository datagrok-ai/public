import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//tags: init
export async function sqlJsInit() : Promise<void> {
  await PackageFunctions.sqlJsInit();
}

//description: Opens SQLite files
//tags: file-handler
//input: list bytes 
//output: list<dataframe> result
//meta.ext: sqlite
export function importSQLite(bytes: Uint8Array) : any {
  return PackageFunctions.importSQLite(bytes);
}
