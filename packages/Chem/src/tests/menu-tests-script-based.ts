import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {category, expect, expectFloat, test} from '@datagrok-libraries/utils/src/test';
import {_testSearchSubstructure,
  _testSearchSubstructureAllParameters,
  _testSearchSubstructureSARSmall,
  loadFileAsText} from './utils';

import {_importSdf} from '../open-chem/sdf-importer';

category('top menu script based', () => {
  test('curate', async () => {
    const t = DG.DataFrame.fromCsv(`Name,smiles
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

    await grok.functions.call('Chem:CurateChemStructures', {'data': t, 'molecules': 'smiles',
      'kekulization': false, 'normalization': true, 'reionization': true,
      'neutralization': true, 'tautomerization': true, 'mainFragment': true});
    v.close();
  });
});
