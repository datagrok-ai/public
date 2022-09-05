import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {category, expect, expectFloat, test} from '@datagrok-libraries/utils/src/test';
import {_testSearchSubstructure,
  _testSearchSubstructureAllParameters,
  _testSearchSubstructureSARSmall,
  loadFileAsText} from './utils';

import {_importSdf} from '../open-chem/sdf-importer';
import { SEMTYPE } from 'datagrok-api/dg';

category('top menu script based', () => {
  test('curate', async () => {
    let t = DG.DataFrame.fromCsv(`Name,smiles
    metal_non,CCC(=O)O[Na]
    metal_st,CCC(=O)[O-].[Na+]
    parent_non,[Na]OC(=O)c1ccccc1
    parent_st,O=C([O-])c1ccccc1
    norm_non,C[N+](C)=CC=C[O-]
    norm_st,CN(C)C=CC=O
    reion_non,C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O
    reion_st,O=S(O)c1ccc(S(=O)(=O)[O-])cc1
    charge_non,O=C([O-])c1ccccc1
    charge_st,O=C(O)c1ccccc1
    tau_non,C1(=CCCCC1)O
    tau_st,O=C1CCCCC1
    main_component_non,CCC1=C(C)C=CC(O)=N1.OC(=O)CCC(=O)O
    main_component_non_st,CCC1=C(C)C=CC(O)=N1`);
    const v = grok.shell.addTableView(t);

    t = await grok.functions.call('Chem:CurateChemStructures', {'data': t, 'molecules': 'smiles',
      'kekulization': false, 'normalization': true, 'reionization': true,
      'neutralization': true, 'tautomerization': true, 'mainFragment': true});
    v.close();
    expect(t.getCol('curated_molecule').get(0), 'CCC(=O)[O-]');
    expect(t.getCol('curated_molecule').get(1), 'CCC(=O)O');
    expect(t.getCol('curated_molecule').get(2), 'O=C([O-])c1ccccc1');
    expect(t.getCol('curated_molecule').get(3), 'O=C(O)c1ccccc1');
    expect(t.getCol('curated_molecule').get(4), 'CN(C)C=CC=O');
    expect(t.getCol('curated_molecule').get(5), 'CN(C)C=CC=O');
    expect(t.getCol('curated_molecule').get(6), 'O=S(O)c1ccc(S(=O)(=O)O)cc1');
    expect(t.getCol('curated_molecule').get(7), 'O=S(O)c1ccc(S(=O)(=O)O)cc1');
    expect(t.getCol('curated_molecule').get(8), 'O=C(O)c1ccccc1');
    expect(t.getCol('curated_molecule').get(9), 'O=C(O)c1ccccc1');
    expect(t.getCol('curated_molecule').get(10), 'O=C1CCCCC1');
    expect(t.getCol('curated_molecule').get(11), 'O=C1CCCCC1');
    expect(t.getCol('curated_molecule').get(12), 'CCc1[nH]c(=O)ccc1C');
    expect(t.getCol('curated_molecule').get(13), 'CCc1[nH]c(=O)ccc1C');
  });

  test('mutate', async () => {
    const mutations = 100;
    const t: DG.DataFrame = await grok.functions.call('Chem:Mutate', {
      'smiles': 'CN1C(CC(O)C1=O)C1=CN=CC=C1', 
      'steps': 1,
      'randomize': true, 
      'maxRandomResults': mutations
    });

    expect(t.rowCount, mutations);
    expect(t.getCol('mutations').semType, SEMTYPE.MOLECULE);
  });

});
