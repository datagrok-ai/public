import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() : void {
  PackageFunctions.info();
}

//output: list<string> result
export async function getFolders() : Promise<string[]> {
  return await PackageFunctions.getFolders();
}

//name: ReinventEditor
//input: funccall call 
//meta.role: editor
export function reinventEditor(call: DG.FuncCall) : void {
  PackageFunctions.reinventEditor(call);
}

//input: string ligand { semType: Molecule }
//input: string optimize 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function runReinvent(ligand: string, optimize: string) : Promise<string> {
  return await PackageFunctions.runReinvent(ligand, optimize);
}

//name: Reinvent
//input: string ligand = 'OC(CN1CCCC1)NC(CCC1)CC1Cl' { semType: Molecule }
//input: string optimize { choices: Reinvent4:getFolders }
//output: dataframe result
//meta.role: hitDesignerFunction
//editor: Reinvent4:ReinventEditor
export async function reinvent(ligand: string, optimize: string) : Promise<any> {
  return await PackageFunctions.reinvent(ligand, optimize);
}

//input: string ligand = 'OC(CN1CCCC1)NC(CCC1)CC1Cl' { semType: Molecule }
//input: string optimize { choices: Reinvent4:getFolders }
//top-menu: Chem | Generate molecules...
//editor: Reinvent4:ReinventEditor
export async function reinventTopMenu(ligand: string, optimize: string) : Promise<void> {
  await PackageFunctions.reinventTopMenu(ligand, optimize);
}
