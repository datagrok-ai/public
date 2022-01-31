import * as DG from 'datagrok-api/dg';

import {drugBankSearchWidget} from './widgets';

export const _package = new DG.Package();

let dbdf: DG.DataFrame;

//tags: init
export async function initDrugBank() {
  dbdf = DG.DataFrame.fromCsv(await _package.files.readAsText('db.csv'));
}

//name: DrugBank Substructure Search
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export async function drugBankSubstructureSearchPanel(mol: string) {
  return await drugBankSearchWidget(mol, 'substructure', dbdf);
}

//name: DrugBank Similarity Search
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export async function drugBankSimilaritySearchPanel(mol: string) {
  return await drugBankSearchWidget(mol, 'similarity', dbdf);
}
