import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Chemistry | Retrosynthesis
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
//meta.allowAddAsColumn: false
//condition: true
export function retroSynthesisPath(molecule: string) {
  return PackageFunctions.retroSynthesisPath(molecule);
}

//name: retrosynthesisTopMenu
export function retrosynthesisTopMenu() {
  return PackageFunctions.retrosynthesisTopMenu();
}

//name: Retrosynthesis Demo
//description: Generate retrosynthesis paths
//meta.demoPath: Cheminformatics | Retrosynthesis
export async function retrosynthesisDemo() {
  return PackageFunctions.retrosynthesisDemo();
}
