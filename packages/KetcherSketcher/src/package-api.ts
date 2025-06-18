import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function ketcherSketcher(): Promise<any> {
    return await grok.functions.call('@datagrok/ketcher-sketcher:KetcherSketcher', {});
  }
}
