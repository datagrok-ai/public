import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: initDrugBank
//tags: init
export async function initDrugBank() {
  return PackageFunctions.initDrugBank();
}

//name: Databases | DrugBank | Substructure Search
//tags: panel, widgets
//input: string mol { semType: Molecule }
//output: widget result
//condition: true
export async function drugBankSubstructureSearchPanel(mol: string) {
  return PackageFunctions.drugBankSubstructureSearchPanel(mol);
}

//name: Databases | DrugBank | Similarity Search
//tags: panel, widgets
//input: string mol { semType: Molecule }
//output: widget result
//condition: true
export async function drugBankSimilaritySearchPanel(mol: string) {
  return PackageFunctions.drugBankSimilaritySearchPanel(mol);
}

//name: Drug Name Molecule
//input: string id 
//output: string result { semType: Molecule }
//meta.role: converter
//meta.inputRegexp: (db\:.+)
export function drugNameMolecule(id: string) {
  return PackageFunctions.drugNameMolecule(id);
}
