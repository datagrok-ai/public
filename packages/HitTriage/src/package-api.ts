import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Funcs {
  export async function hitTriageAppTreeBrowser(treeNode: any, browseView: DG.View): Promise<any> {
    return await grok.functions.call('HitTriage:HitTriageAppTreeBrowser', { treeNode, browseView });
  }

  export async function hitDesignAppTreeBrowser(treeNode: any, browseView: DG.View): Promise<any> {
    return await grok.functions.call('HitTriage:HitDesignAppTreeBrowser', { treeNode, browseView });
  }

  export async function peptiHitAppTreeBrowser(treeNode: any, browseView: DG.View): Promise<any> {
    return await grok.functions.call('HitTriage:PeptiHitAppTreeBrowser', { treeNode, browseView });
  }

  export async function hitTriageApp(): Promise<any> {
    return await grok.functions.call('HitTriage:HitTriageApp', {});
  }

  export async function hitDesignApp(): Promise<any> {
    return await grok.functions.call('HitTriage:HitDesignApp', {});
  }

  export async function peptiHitApp(): Promise<any> {
    return await grok.functions.call('HitTriage:PeptiHitApp', {});
  }

  export async function demoFileIngest(): Promise<any> {
    return await grok.functions.call('HitTriage:DemoFileIngest', {});
  }

  export async function demoFileIngest1(): Promise<any> {
    return await grok.functions.call('HitTriage:DemoFileIngest1', {});
  }

  export async function demoFileIngest2(numberOfMolecules: number): Promise<any> {
    return await grok.functions.call('HitTriage:DemoFileIngest2', { numberOfMolecules });
  }

  export async function demoFileSubmit(df: DG.DataFrame, molecules: string): Promise<any> {
    return await grok.functions.call('HitTriage:DemoFileSubmit', { df, molecules });
  }

  export async function gasteigerCellRenderer(): Promise<any> {
    return await grok.functions.call('HitTriage:GasteigerCellRenderer', {});
  }

  export async function htPackageSettingEditor(propList: any): Promise<any> {
    return await grok.functions.call('HitTriage:HtPackageSettingEditor', { propList });
  }
}
