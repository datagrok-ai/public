import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Funcs {
  export async function sqlJsInit(): Promise<any> {
    return await grok.functions.call('Sqlite:SqlJsInit', {});
  }

  //Opens SQLite files
  export async function importSQLite(bytes: any): Promise<any> {
    return await grok.functions.call('Sqlite:ImportSQLite', { bytes });
  }
}
