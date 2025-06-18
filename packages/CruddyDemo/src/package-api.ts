import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function northwindDemo(): Promise<any> {
    return await grok.functions.call('cruddydemo:NorthwindDemo', {});
  }

  export async function chemblDemo(): Promise<any> {
    return await grok.functions.call('cruddydemo:ChemblDemo', {});
  }
}
