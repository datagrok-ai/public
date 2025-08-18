import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function signalsApp(): Promise<DG.View> {
    return await grok.functions.call('SignalsLink:SignalsApp', {});
  }

  export async function signalsAppTreeBrowser(treeNode: any ): Promise<void> {
    return await grok.functions.call('SignalsLink:SignalsAppTreeBrowser', { treeNode });
  }

  export async function saveToSignalsEln(data: string ): Promise<void> {
    return await grok.functions.call('SignalsLink:SaveToSignalsEln', { data });
  }
}
