import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Chemistry | Retrosynthesis
//input: string smiles { semType: Molecule }
//output: widget result
//meta.allowAddAsColumn: false
//meta.role: widgets,panel
//meta.domain: chem
//condition: true
export function retroSynthesisPath(molecule: string) : any {
  return PackageFunctions.retroSynthesisPath(molecule);
}

//name: retrosynthesisTopMenu
export function retrosynthesisTopMenu() : void {
  PackageFunctions.retrosynthesisTopMenu();
}

//name: Retrosynthesis Demo
//description: Generate retrosynthesis paths
//meta.demoPath: Cheminformatics | Retrosynthesis
export async function retrosynthesisDemo() : Promise<void> {
  await PackageFunctions.retrosynthesisDemo();
}
