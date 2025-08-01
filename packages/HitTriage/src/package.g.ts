import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: hitTriageAppTreeBrowser
//input: dynamic treeNode 
//input: view browseView 
//output: dynamic result
export async function hitTriageAppTreeBrowser(treeNode: any, browsePanel: any) {
  return PackageFunctions.hitTriageAppTreeBrowser(treeNode, browsePanel);
}

//name: hitDesignAppTreeBrowser
//input: dynamic treeNode 
//input: view browseView 
//output: dynamic result
export async function hitDesignAppTreeBrowser(treeNode: any, browsePanel: any) {
  return PackageFunctions.hitDesignAppTreeBrowser(treeNode, browsePanel);
}

//name: peptiHitAppTreeBrowser
//input: dynamic treeNode 
//input: view browseView 
//output: dynamic result
export async function peptiHitAppTreeBrowser(treeNode: any, browsePanel: any) {
  return PackageFunctions.peptiHitAppTreeBrowser(treeNode, browsePanel);
}

//name: Hit Triage
//tags: app
//output: view result
//meta.browsePath: Chem
export async function hitTriageApp() {
  return PackageFunctions.hitTriageApp();
}

//name: Hit Design
//tags: app
//output: view result
//meta.icon: images/icons/hit-design-icon.png
//meta.browsePath: Chem
export async function hitDesignApp() {
  return PackageFunctions.hitDesignApp();
}

//name: PeptiHit
//tags: app
//output: view result
//meta.icon: images/icons/pepti-hit-icon.png
//meta.browsePath: Peptides
export async function peptiHitApp() {
  return PackageFunctions.peptiHitApp();
}

//name: Demo Molecules 100
//tags: HitTriageDataSource
//output: dataframe result
export async function demoFileIngest() {
  return PackageFunctions.demoFileIngest();
}

//name: Demo Molecules 5000
//tags: HitTriageDataSource
//output: dataframe result
export async function demoFileIngest1() {
  return PackageFunctions.demoFileIngest1();
}

//name: Demo Molecules variable
//tags: HitTriageDataSource
//input: int numberOfMolecules 
//output: dataframe result
export async function demoFileIngest2(numberOfMolecules: number) {
  return PackageFunctions.demoFileIngest2(numberOfMolecules);
}

//name: Demo File Submit
//tags: HitTriageSubmitFunction
//input: dataframe df 
//input: string molecules 
export async function demoFileSubmit(df: DG.DataFrame, molecules: string) {
  return PackageFunctions.demoFileSubmit(df, molecules);
}

//name: gasteigerRenderer
//tags: cellRenderer
//output: dynamic result
//meta.cellType: customGasteigerPNG
//meta.columnTags: quality=customGasteigerPNG
export function gasteigerCellRenderer() {
  return PackageFunctions.gasteigerCellRenderer();
}

//name: Hit Triage package settings editor
//tags: packageSettingsEditor
//input: object propList 
//output: dynamic result
export async function htPackageSettingEditor(properties: any) {
  return PackageFunctions.htPackageSettingEditor(properties);
}
