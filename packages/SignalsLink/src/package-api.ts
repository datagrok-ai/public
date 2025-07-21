import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function signalsApp(): Promise<any> {
    return await grok.functions.call('SignalsLink:SignalsApp', {});
  }

  export async function signalsAppTreeBrowser(treeNode: any, browseView: DG.View): Promise<any> {
    return await grok.functions.call('SignalsLink:SignalsAppTreeBrowser', { treeNode, browseView });
  }

  export async function saveToSignalsEln(data: string): Promise<any> {
    return await grok.functions.call('SignalsLink:SaveToSignalsEln', { data });
  }
}
