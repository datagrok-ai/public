import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {

}

export namespace funcs {
  export async function clinicalCaseApp(): Promise<any> {
    return await grok.functions.call('ClinicalCase:ClinicalCaseApp', {});
  }

  export async function clinicalCaseAppTreeBrowser(treeNode: any, browseView: DG.View): Promise<any> {
    return await grok.functions.call('ClinicalCase:ClinicalCaseAppTreeBrowser', { treeNode, browseView });
  }

  export async function clinicalCaseFolderLauncher(folder: DG.FileInfo): Promise<any> {
    return await grok.functions.call('ClinicalCase:ClinicalCaseFolderLauncher', { folder });
  }
}
