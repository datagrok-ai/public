import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function sqlJsInit(): Promise<any> {
    return await grok.functions.call('@datagrok/sqlite:SqlJsInit', {});
  }

  //Opens SQLite files
  export async function importSQLite(bytes: any): Promise<any> {
    return await grok.functions.call('@datagrok/sqlite:ImportSQLite', { bytes });
  }
}
