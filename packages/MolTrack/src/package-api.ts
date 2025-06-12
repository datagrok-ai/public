import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Queries {
  export async function getCompounds(): Promise<DG.DataFrame> {
    return await grok.data.query('MolTrack:GetCompounds', {});
  }
}

export namespace Funcs {
  export async function molTrackApp(): Promise<any> {
    return await grok.functions.call('MolTrack:MolTrackApp', {});
  }

  export async function cddVaultAppTreeBrowser(treeNode: any, browseView: DG.View): Promise<any> {
    return await grok.functions.call('MolTrack:CddVaultAppTreeBrowser', { treeNode, browseView });
  }
}
