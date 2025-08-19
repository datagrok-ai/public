import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: hitTriageAppTreeBrowser
//input: dynamic treeNode 
export async function hitTriageAppTreeBrowser(treeNode: any) {
  return PackageFunctions.hitTriageAppTreeBrowser(treeNode);
}

//name: hitDesignAppTreeBrowser
//input: dynamic treeNode 
export async function hitDesignAppTreeBrowser(treeNode: any) {
  return PackageFunctions.hitDesignAppTreeBrowser(treeNode);
}

//name: peptiHitAppTreeBrowser
//input: dynamic treeNode 
export async function peptiHitAppTreeBrowser(treeNode: any) {
  return PackageFunctions.peptiHitAppTreeBrowser(treeNode);
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
//input: int numberOfMolecules { description: Molecules counts }
//output: dataframe result
export async function demoFileIngest2(numberOfMolecules: number) {
  return PackageFunctions.demoFileIngest2(numberOfMolecules);
}

//name: Demo File Submit
//tags: HitTriageSubmitFunction
//input: dataframe df { description: Dataframe }
//input: string molecules { description: Molecules column name }
export async function demoFileSubmit(df: DG.DataFrame, molecules: string) {
  return PackageFunctions.demoFileSubmit(df, molecules);
}

//name: gasteigerRenderer
//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: customGasteigerPNG
//meta.columnTags: quality=customGasteigerPNG
export function gasteigerCellRenderer() {
  return PackageFunctions.gasteigerCellRenderer();
}

//name: Hit Triage package settings editor
//tags: packageSettingsEditor
//input: object propList 
//output: widget result
export async function htPackageSettingEditor(properties: any) {
  return PackageFunctions.htPackageSettingEditor(properties);
}
