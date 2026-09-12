import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//meta.role: Init
export async function sqlJsInit() : Promise<void> {
  await PackageFunctions.sqlJsInit();
}

//description: Opens SQLite files
//input: list bytes 
//output: list<dataframe> result
//meta.role: FileHandler
//meta.ext: sqlite
export function importSQLite(bytes: Uint8Array) : any {
  return PackageFunctions.importSQLite(bytes);
}
