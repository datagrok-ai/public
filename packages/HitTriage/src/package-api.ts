import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function hitTriageAppTreeBrowser(treeNode: any, browseView: DG.View): Promise<any> {
    return await grok.functions.call('@datagrok/hit-triage:HitTriageAppTreeBrowser', { treeNode, browseView });
  }

  export async function hitDesignAppTreeBrowser(treeNode: any, browseView: DG.View): Promise<any> {
    return await grok.functions.call('@datagrok/hit-triage:HitDesignAppTreeBrowser', { treeNode, browseView });
  }

  export async function peptiHitAppTreeBrowser(treeNode: any, browseView: DG.View): Promise<any> {
    return await grok.functions.call('@datagrok/hit-triage:PeptiHitAppTreeBrowser', { treeNode, browseView });
  }

  export async function hitTriageApp(): Promise<any> {
    return await grok.functions.call('@datagrok/hit-triage:HitTriageApp', {});
  }

  export async function hitDesignApp(): Promise<any> {
    return await grok.functions.call('@datagrok/hit-triage:HitDesignApp', {});
  }

  export async function peptiHitApp(): Promise<any> {
    return await grok.functions.call('@datagrok/hit-triage:PeptiHitApp', {});
  }

  export async function demoFileIngest(): Promise<any> {
    return await grok.functions.call('@datagrok/hit-triage:DemoFileIngest', {});
  }

  export async function demoFileIngest1(): Promise<any> {
    return await grok.functions.call('@datagrok/hit-triage:DemoFileIngest1', {});
  }

  export async function demoFileIngest2(numberOfMolecules: number): Promise<any> {
    return await grok.functions.call('@datagrok/hit-triage:DemoFileIngest2', { numberOfMolecules });
  }

  export async function demoFileSubmit(df: DG.DataFrame, molecules: string): Promise<any> {
    return await grok.functions.call('@datagrok/hit-triage:DemoFileSubmit', { df, molecules });
  }

  export async function gasteigerCellRenderer(): Promise<any> {
    return await grok.functions.call('@datagrok/hit-triage:GasteigerCellRenderer', {});
  }

  export async function htPackageSettingEditor(propList: any): Promise<any> {
    return await grok.functions.call('@datagrok/hit-triage:HtPackageSettingEditor', { propList });
  }
}
