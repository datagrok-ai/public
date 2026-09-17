import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//meta.role: init
export async function initDrugBank() : Promise<void> {
  await PackageFunctions.initDrugBank();
}

//name: Databases | DrugBank | Substructure Search
//input: string mol { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
//condition: true
export async function drugBankSubstructureSearchPanel(mol: string) : Promise<any> {
  return await PackageFunctions.drugBankSubstructureSearchPanel(mol);
}

//name: Databases | DrugBank | Similarity Search
//input: string mol { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
//condition: true
export async function drugBankSimilaritySearchPanel(mol: string) : Promise<any> {
  return await PackageFunctions.drugBankSimilaritySearchPanel(mol);
}

//name: Drug Name Molecule
//input: string id 
//output: string result { semType: Molecule }
//meta.role: converter
//meta.inputRegexp: (db\:.+)
//connection: DrugBank
export function drugNameMolecule(id: string) : string {
  return PackageFunctions.drugNameMolecule(id);
}
