import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace queries {
  export async function getCompounds(): Promise<DG.DataFrame> {
    return await grok.data.query('molTrack:GetCompounds', {});
  }
}

export namespace funcs {
  export async function molTrackApp(): Promise<any> {
    return await grok.functions.call('molTrack:MolTrackApp', {});
  }

  export async function cddVaultAppTreeBrowser(treeNode: any, browseView: DG.View): Promise<any> {
    return await grok.functions.call('molTrack:CddVaultAppTreeBrowser', { treeNode, browseView });
  }
}
