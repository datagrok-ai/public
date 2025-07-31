import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
//output: dynamic result
export function info() {
  return PackageFunctions.info();
}

//name: getFolders
//output: list result
export async function getFolders() {
  return PackageFunctions.getFolders();
}

//name: ReinventEditor
//tags: editor
//input: funccall call 
export function reinventEditor(call: DG.FuncCall) {
  return PackageFunctions.reinventEditor(call);
}

//name: runReinvent
//input: string ligand { semType: Molecule }
//input: string optimize 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function runReinvent(ligand: string, optimize: string) {
  return PackageFunctions.runReinvent(ligand, optimize);
}

//name: Reinvent
//tags: HitDesignerFunction
//input: string ligand { semType: Molecule; default: 'OC(CN1CCCC1)NC(CCC1)CC1Cl' }
//input: string optimize { choices: Reinvent4:getFolders }
//output: dataframe result
//editor: Reinvent4
export async function reinvent(ligand: string, optimize: string) {
  return PackageFunctions.reinvent(ligand, optimize);
}

//name: reinventTopMenu
//input: string ligand { semType: Molecule; default: 'OC(CN1CCCC1)NC(CCC1)CC1Cl' }
//input: string optimize { choices: Reinvent4:getFolders }
//output: dynamic result
//top-menu: Chem | Generate molecules...
//editor: Reinvent4
export async function reinventTopMenu(ligand: string, optimize: string) {
  return PackageFunctions.reinventTopMenu(ligand, optimize);
}
