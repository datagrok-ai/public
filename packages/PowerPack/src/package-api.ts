import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function compareColumns(): Promise<void> {
    return await grok.functions.call('PowerPack:CompareColumns', {});
  }

  export async function addNewColumnDialog(call: any ): Promise<void> {
    return await grok.functions.call('PowerPack:AddNewColumnDialog', { call });
  }

  export async function welcomeView(): Promise<DG.View> {
    return await grok.functions.call('PowerPack:WelcomeView', {});
  }

  export async function recentProjectsWidget(): Promise<any> {
    return await grok.functions.call('PowerPack:RecentProjectsWidget', {});
  }

  export async function communityWidget(): Promise<any> {
    return await grok.functions.call('PowerPack:CommunityWidget', {});
  }

  export async function webWidget(): Promise<any> {
    return await grok.functions.call('PowerPack:WebWidget', {});
  }

  export async function htmlWidget(): Promise<any> {
    return await grok.functions.call('PowerPack:HtmlWidget', {});
  }

  export async function learnWidget(): Promise<any> {
    return await grok.functions.call('PowerPack:LearnWidget', {});
  }

  export async function kpiWidget(): Promise<any> {
    return await grok.functions.call('PowerPack:KpiWidget', {});
  }

  export async function isFormulaColumn(col: DG.Column ): Promise<boolean> {
    return await grok.functions.call('PowerPack:IsFormulaColumn', { col });
  }

  export async function formulaWidget(col: DG.Column ): Promise<any> {
    return await grok.functions.call('PowerPack:FormulaWidget', { col });
  }

  //Functions
  export async function functionSearch(s: string ): Promise<any> {
    return await grok.functions.call('PowerPack:FunctionSearch', { s });
  }

  //Scripts
  export async function scriptsSearch(s: string ): Promise<any> {
    return await grok.functions.call('PowerPack:ScriptsSearch', { s });
  }

  //Queries
  export async function queriesSearch(s: string ): Promise<any> {
    return await grok.functions.call('PowerPack:QueriesSearch', { s });
  }

  //Users
  export async function usersSearch(s: string ): Promise<any> {
    return await grok.functions.call('PowerPack:UsersSearch', { s });
  }

  //Groups
  export async function groupsSearch(s: string ): Promise<any> {
    return await grok.functions.call('PowerPack:GroupsSearch', { s });
  }

  //Dockers
  export async function dockerSearch(s: string ): Promise<any> {
    return await grok.functions.call('PowerPack:DockerSearch', { s });
  }

  //Help
  export async function helpSearch(s: string ): Promise<any> {
    return await grok.functions.call('PowerPack:HelpSearch', { s });
  }

  //Apps
  export async function appSearch(s: string ): Promise<any> {
    return await grok.functions.call('PowerPack:AppSearch', { s });
  }

  //Connections
  export async function connectionsSearch(s: string ): Promise<any> {
    return await grok.functions.call('PowerPack:ConnectionsSearch', { s });
  }

  //Protein Data Bank
  export async function pdbSearch(s: string ): Promise<any> {
    return await grok.functions.call('PowerPack:PdbSearch', { s });
  }

  //PubChem
  export async function pubChemSearch(s: string ): Promise<any> {
    return await grok.functions.call('PowerPack:PubChemSearch', { s });
  }

  //PubChem
  export async function wikiSearch(s: string ): Promise<any> {
    return await grok.functions.call('PowerPack:WikiSearch', { s });
  }

  export async function newUsersSearchWidget(s: string ): Promise<any> {
    return await grok.functions.call('PowerPack:NewUsersSearchWidget', { s });
  }

  //Denial Search
  export async function denialSearch(s: string ): Promise<any> {
    return await grok.functions.call('PowerPack:DenialSearch', { s });
  }

  export async function formulaLinesDialog(src: DG.DataFrame ): Promise<void> {
    return await grok.functions.call('PowerPack:FormulaLinesDialog', { src });
  }

  export async function powerPackInit(): Promise<void> {
    return await grok.functions.call('PowerPack:PowerPackInit', {});
  }

  //Windows Manager
  export async function windowsManager(): Promise<void> {
    return await grok.functions.call('PowerPack:WindowsManager', {});
  }

  //Open "Viewer Gallery" dialog
  export async function viewerDialog(tv: any ): Promise<void> {
    return await grok.functions.call('PowerPack:ViewerDialog', { tv });
  }

  //ViewerGallery
  export async function viewerGallery(): Promise<void> {
    return await grok.functions.call('PowerPack:ViewerGallery', {});
  }

  export async function markdownFileViewer(file: DG.FileInfo ): Promise<DG.View> {
    return await grok.functions.call('PowerPack:MarkdownFileViewer', { file });
  }

  //Opens Excel file
  export async function xlsxFileHandler(bytes: any , sheetName?: string ): Promise<any> {
    return await grok.functions.call('PowerPack:XlsxFileHandler', { bytes, sheetName });
  }
}
