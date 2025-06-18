import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function signalsApp(): Promise<any> {
    return await grok.functions.call('@datagrok/signals-link:SignalsApp', {});
  }

  export async function signalsAppTreeBrowser(treeNode: any, browseView: DG.View): Promise<any> {
    return await grok.functions.call('@datagrok/signals-link:SignalsAppTreeBrowser', { treeNode, browseView });
  }

  export async function saveToSignalsEln(data: string): Promise<any> {
    return await grok.functions.call('@datagrok/signals-link:SaveToSignalsEln', { data });
  }
}
