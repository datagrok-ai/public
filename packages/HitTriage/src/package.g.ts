import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//input: dynamic treeNode 
//meta.app: Hit Triage
export async function hitTriageAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.hitTriageAppTreeBrowser(treeNode);
}

//input: dynamic treeNode 
//meta.app: Hit Design
export async function hitDesignAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.hitDesignAppTreeBrowser(treeNode);
}

//input: dynamic treeNode 
//meta.app: PeptiHit
export async function peptiHitAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.peptiHitAppTreeBrowser(treeNode);
}

//name: Hit Triage
//output: view result
//meta.role: app
//meta.browsePath: Chem
export async function hitTriageApp() : Promise<any> {
  return await PackageFunctions.hitTriageApp();
}

//name: Hit Design
//output: view result
//meta.role: app
//meta.icon: images/icons/hit-design-icon.png
//meta.browsePath: Chem
export async function hitDesignApp() : Promise<any> {
  return await PackageFunctions.hitDesignApp();
}

//name: PeptiHit
//output: view result
//meta.role: app
//meta.icon: images/icons/pepti-hit-icon.png
//meta.browsePath: Peptides
export async function peptiHitApp() : Promise<any> {
  return await PackageFunctions.peptiHitApp();
}

//name: Demo Molecules 100
//output: dataframe result
//meta.role: hitTriageDataSource
export async function demoFileIngest() : Promise<any> {
  return await PackageFunctions.demoFileIngest();
}

//name: Demo Molecules 5000
//output: dataframe result
//meta.role: hitTriageDataSource
export async function demoFileIngest1() : Promise<any> {
  return await PackageFunctions.demoFileIngest1();
}

//name: Demo Molecules variable
//input: int numberOfMolecules { description: Molecules counts }
//output: dataframe result
//meta.role: hitTriageDataSource
export async function demoFileIngest2(numberOfMolecules: number) : Promise<any> {
  return await PackageFunctions.demoFileIngest2(numberOfMolecules);
}

//name: Demo File Submit
//input: dataframe df { description: Dataframe }
//input: string molecules { description: Molecules column name }
//meta.role: hitTriageSubmitFunction
export async function demoFileSubmit(df: DG.DataFrame, molecules: string) : Promise<void> {
  await PackageFunctions.demoFileSubmit(df, molecules);
}

//name: registerMoleculesToViD
export async function registerMoleculesToViD() : Promise<void> {
  await PackageFunctions.registerMoleculesToViD();
}

//name: Hit Design V-iD
//input: semantic_value vid { semType: HIT_DESIGN_VID }
//output: widget result
//meta.role: panel
export function hitDesignVidPanel(vid: DG.SemanticValue) : any {
  return PackageFunctions.hitDesignVidPanel(vid);
}

//name: gasteigerRenderer
//output: grid_cell_renderer result
//meta.cellType: customGasteigerPNG
//meta.columnTags: quality=customGasteigerPNG
//meta.role: cellRenderer
export function gasteigerCellRenderer() : any {
  return PackageFunctions.gasteigerCellRenderer();
}

//name: Hit Triage package settings editor
//input: object propList 
//output: widget result
//meta.role: packageSettingsEditor
export async function htPackageSettingEditor(properties: any) : Promise<any> {
  return await PackageFunctions.htPackageSettingEditor(properties);
}
