import * as DG from 'datagrok-api/dg';

import {drugBankSearchWidget, drugNameMoleculeWidget} from './widgets';

export const _package = new DG.Package();

let dbdf: DG.DataFrame;

//tags: init
export async function initDrugBank(): Promise<void> {
  dbdf = DG.DataFrame.fromCsv(await _package.files.readAsText('db.csv'));
}

//name: DrugBank Substructure Search
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export async function drugBankSubstructureSearchPanel(mol: string): Promise<DG.Widget> {
  if (!dbdf)
    await initDrugBank();
  return drugBankSearchWidget(mol, 'substructure', dbdf);
}

//name: DrugBank Similarity Search
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export async function drugBankSimilaritySearchPanel(mol: string): Promise<DG.Widget> {
  if (!dbdf)
    await initDrugBank();
  return drugBankSearchWidget(mol, 'similarity', dbdf);
}

//name: Drug Name Molecule
//meta.role: converter
//meta.inputRegexp: (db\:.+)
//connection: DrugBank
//input: string id
//output: string smiles { semType: Molecule }
export async function drugNameMolecule(id: string): Promise<string> {
  if (!dbdf)
    await initDrugBank();
  return drugNameMoleculeWidget(id, dbdf) ?? '';
}
