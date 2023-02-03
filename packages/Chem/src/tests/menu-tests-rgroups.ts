import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, expect, expectFloat, test, delay, before} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {Fingerprint} from '../utils/chem-common';
import {createTableView, readDataframe} from './utils';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';

import {
  _testSearchSubstructure,
  _testSearchSubstructureAllParameters,
  _testSearchSubstructureSARSmall,
  loadFileAsText
} from './utils';
import {findSimilar, getSimilarities} from '../package';
import {chemDiversitySearch} from '../analysis/chem-diversity-viewer';

const t = DG.DataFrame.fromCsv(`smiles
O=C1CN=C(c2ccccc2N1)C3CCCCC3
CN1C(=O)CN=C(c2ccccc12)C3CCCCC3
CCCCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
CC(C)CCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
O=C1CN=C(c2ccccc2N1CC3CCCCC3)C4CCCCC4
O=C1CN=C(c2cc(Cl)ccc2N1)C3CCCCC3
CN1C(=O)CN=C(c2cc(Cl)ccc12)C3CCCCC3`);

category('top menu r-groups', () => {

  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('mcs', async () => {
    const mcs = await grok.functions.call('Chem:FindMCS', {'molecules': 'smiles', 'df': t, 'returnSmarts': false});
    expect(mcs, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1');
  });

  test('rgroups', async () => {
    const rgroups: DG.DataFrame = await grok.functions.call('Chem:FindRGroups', {
      'molecules': 'smiles', 
      'df': t, 
      'core': 'c1ccccc1',
      'prefix': 'R'
    });
    expect(rgroups.getCol('R1').get(0), '*C(=NCC(=O)N[1*])C1CCCCC1');
    expect(rgroups.getCol('R1').get(1), '*C(=NCC(=O)N([1*])C)C1CCCCC1');
    expect(rgroups.getCol('R1').get(2), '*C(=NCC(=O)N([1*])CCCC)C1CCCCC1');
    expect(rgroups.getCol('R1').get(3), '*C(=NCC(=O)N([1*])CCC(C)C)C1CCCCC1');
    expect(rgroups.getCol('R1').get(4), '*C(=NCC(=O)N([1*])CC1CCCCC1)C1CCCCC1');
    expect(rgroups.getCol('R1').get(5), '*C(=NCC(=O)N[1*])C1CCCCC1');
    expect(rgroups.getCol('R1').get(6), '*C(=NCC(=O)N([1*])C)C1CCCCC1');
    expect(rgroups.getCol('R2').get(5), '[2*]Cl');
    expect(rgroups.getCol('R2').get(6), '[2*]Cl');
  });

});

