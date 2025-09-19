import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//input: dynamic treeNode 
export async function hitTriageAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.hitTriageAppTreeBrowser(treeNode);
}

//input: dynamic treeNode 
export async function hitDesignAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.hitDesignAppTreeBrowser(treeNode);
}

//input: dynamic treeNode 
export async function peptiHitAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.peptiHitAppTreeBrowser(treeNode);
}

//name: Hit Triage
//tags: app
//output: view result
//meta.browsePath: Chem
export async function hitTriageApp() : Promise<any> {
  return await PackageFunctions.hitTriageApp();
}

//name: Hit Design
//tags: app
//output: view result
//meta.icon: images/icons/hit-design-icon.png
//meta.browsePath: Chem
export async function hitDesignApp() : Promise<any> {
  return await PackageFunctions.hitDesignApp();
}

//name: PeptiHit
//tags: app
//output: view result
//meta.icon: images/icons/pepti-hit-icon.png
//meta.browsePath: Peptides
export async function peptiHitApp() : Promise<any> {
  return await PackageFunctions.peptiHitApp();
}

//name: Demo Molecules 100
//tags: HitTriageDataSource
//output: dataframe result
export async function demoFileIngest() : Promise<any> {
  return await PackageFunctions.demoFileIngest();
}

//name: Demo Molecules 5000
//tags: HitTriageDataSource
//output: dataframe result
export async function demoFileIngest1() : Promise<any> {
  return await PackageFunctions.demoFileIngest1();
}

//name: Demo Molecules variable
//tags: HitTriageDataSource
//input: int numberOfMolecules { description: Molecules counts }
//output: dataframe result
export async function demoFileIngest2(numberOfMolecules: number) : Promise<any> {
  return await PackageFunctions.demoFileIngest2(numberOfMolecules);
}

//name: Demo File Submit
//tags: HitTriageSubmitFunction
//input: dataframe df { description: Dataframe }
//input: string molecules { description: Molecules column name }
export async function demoFileSubmit(df: DG.DataFrame, molecules: string) : Promise<void> {
  await PackageFunctions.demoFileSubmit(df, molecules);
}

//name: gasteigerRenderer
//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: customGasteigerPNG
//meta.columnTags: quality=customGasteigerPNG
export function gasteigerCellRenderer() : any {
  return PackageFunctions.gasteigerCellRenderer();
}

//name: Hit Triage package settings editor
//tags: packageSettingsEditor
//input: object propList 
//output: widget result
export async function htPackageSettingEditor(properties: any) : Promise<any> {
  return await PackageFunctions.htPackageSettingEditor(properties);
}
