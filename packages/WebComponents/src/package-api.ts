import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Funcs {
  export async function init(): Promise<any> {
    return await grok.functions.call('Webcomponents:Init', {});
  }
}
