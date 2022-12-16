import * as DG from 'datagrok-api/dg';

import {drugNameMoleculeConvert, searchWidget} from './widgets';

export const _package = new DG.Package();

let dbdf: DG.DataFrame;
let synonymsCol: DG.Column<string>;
let moleculeCol: DG.Column<string>;
let dbdfRowCount: number;

//tags: init
export async function initDrugBank(): Promise<void> {
  dbdf = (await _package.files.readBinaryDataFrames('drugbank-open-structures.d42'))[0];
  synonymsCol = dbdf.getCol('SYNONYMS');
  moleculeCol = dbdf.getCol('molecule');
  dbdfRowCount = dbdf.rowCount;
}

//name: DrugBank | Substructure Search
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export async function drugBankSubstructureSearchPanel(mol: string): Promise<DG.Widget> {
  return searchWidget(mol, 'substructure', dbdf);
}

//name: DrugBank | Similarity Search
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export async function drugBankSimilaritySearchPanel(mol: string): Promise<DG.Widget> {
  return searchWidget(mol, 'similarity', dbdf);
}

//name: Drug Name Molecule
//meta.role: converter
//meta.inputRegexp: (db\:.+)
//connection: DrugBank
//input: string id
//output: string smiles { semType: Molecule }
export function drugNameMolecule(id: string): string {
  return drugNameMoleculeConvert(id, dbdfRowCount, synonymsCol, moleculeCol);
}
