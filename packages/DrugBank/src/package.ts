import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {SEARCH_TYPE, drugNameMoleculeConvert, searchWidget} from './widgets';

export * from './package.g';
export const _package = new DG.Package();

let dbdf: DG.DataFrame;
let synonymsCol: DG.Column<string>;
let moleculeCol: DG.Column<string>;
let dbdfRowCount: number;


export * from './package.g';
export class PackageFunctions {
  @grok.decorators.init()
  static async initDrugBank(): Promise<void> {
    dbdf = (await _package.files.readBinaryDataFrames('drugbank-open-structures.d42'))[0];
    synonymsCol = dbdf.getCol('SYNONYMS');
    moleculeCol = dbdf.getCol('molecule');
    dbdfRowCount = dbdf.rowCount;
  }


  @grok.decorators.panel({
    'tags': [
      'widgets',
    ],
    'name': 'Databases | DrugBank | Substructure Search',
    'condition': 'true',
  })
  static async drugBankSubstructureSearchPanel(
    @grok.decorators.param({'options': {'semType': 'Molecule'}}) mol: string,
  ): Promise<DG.Widget> {
    return mol ? searchWidget(mol, SEARCH_TYPE.SUBSTRUCTURE, dbdf) : new DG.Widget(ui.divText('SMILES is empty'));
  }


  @grok.decorators.panel({
    'tags': [
      'widgets',
    ],
    'name': 'Databases | DrugBank | Similarity Search',
    'condition': 'true',
  })
  static async drugBankSimilaritySearchPanel(
    @grok.decorators.param({'options': {'semType': 'Molecule'}}) mol: string,
  ): Promise<DG.Widget> {
    return mol ? searchWidget(mol, SEARCH_TYPE.SIMILARITY, dbdf) : new DG.Widget(ui.divText('SMILES is empty'));
  }


  //Add after release 1.26.1
  //'connection': 'DrugBank',
  @grok.decorators.func({
    'meta': {
      'role': 'converter',
      'inputRegexp': '(db\\:.+)',
    },
    'name': 'Drug Name Molecule',
  })
  static drugNameMolecule(id: string): string {
    return drugNameMoleculeConvert(id, dbdfRowCount, synonymsCol, moleculeCol);
  }
}
