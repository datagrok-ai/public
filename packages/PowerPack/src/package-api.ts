import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function compareColumns(): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:CompareColumns', {});
  }

  export async function addNewColumnDialog(call: any): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:AddNewColumnDialog', { call });
  }

  export async function welcomeView(): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:WelcomeView', {});
  }

  export async function recentProjectsWidget(): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:RecentProjectsWidget', {});
  }

  export async function communityWidget(): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:CommunityWidget', {});
  }

  export async function webWidget(): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:WebWidget', {});
  }

  export async function htmlWidget(): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:HtmlWidget', {});
  }

  export async function learnWidget(): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:LearnWidget', {});
  }

  export async function kpiWidget(): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:KpiWidget', {});
  }

  export async function isFormulaColumn(col: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:IsFormulaColumn', { col });
  }

  export async function formulaWidget(col: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:FormulaWidget', { col });
  }

  //Functions
  export async function functionSearch(s: string): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:FunctionSearch', { s });
  }

  //Scripts
  export async function scriptsSearch(s: string): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:ScriptsSearch', { s });
  }

  //Users
  export async function usersSearch(s: string): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:UsersSearch', { s });
  }

  //Protein Data Bank
  export async function pdbSearch(s: string): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:PdbSearch', { s });
  }

  //PubChem
  export async function pubChemSearch(s: string): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:PubChemSearch', { s });
  }

  //PubChem
  export async function wikiSearch(s: string): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:WikiSearch', { s });
  }

  export async function newUsersSearchWidget(s: string): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:NewUsersSearchWidget', { s });
  }

  export async function formulaLinesDialog(src: DG.DataFrame): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:FormulaLinesDialog', { src });
  }

  export async function powerPackInit(): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:PowerPackInit', {});
  }

  //Windows Manager
  export async function windowsManager(): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:WindowsManager', {});
  }

  //Open "Viewer Gallery" dialog
  export async function viewerDialog(tv: any): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:ViewerDialog', { tv });
  }

  //ViewerGallery
  export async function viewerGallery(): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:ViewerGallery', {});
  }

  export async function markdownFileViewer(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:MarkdownFileViewer', { file });
  }

  //Opens Excel file
  export async function xlsxFileHandler(bytes: any, sheetName: string): Promise<any> {
    return await grok.functions.call('@datagrok/power-pack:XlsxFileHandler', { bytes, sheetName });
  }
}
