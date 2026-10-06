import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Chemistry | Retrosynthesis
//description: Predict retrosynthesis routes for a molecule using AiZynthFinder.
//input: string smiles { semType: Molecule; description: Target molecule to plan a synthesis route for }
//output: widget result
//meta.allowAddAsColumn: false
//meta.role: widgets,panel
//meta.domain: chem
//condition: true
export function retroSynthesisPath(molecule: string) : any {
  return PackageFunctions.retroSynthesisPath(molecule);
}

//name: Retrosynthesis Demo
//description: Generate retrosynthesis paths
//meta.demoPath: Cheminformatics | Retrosynthesis
export async function retrosynthesisDemo() : Promise<void> {
  await PackageFunctions.retrosynthesisDemo();
}
