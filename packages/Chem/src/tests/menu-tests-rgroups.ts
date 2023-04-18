import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, expect, test, before, after} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {readDataframe} from './utils';
import {findMCS, findRGroups} from '../scripts-api';


category('top menu r-groups', () => {
  let empty: DG.DataFrame;
  let malformed: DG.DataFrame;
  let coreEmpty: string;
  let coreMalformed: string;

  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    empty = await readDataframe('tests/sar-small_empty_vals.csv');
    await grok.data.detectSemanticTypes(empty);
    malformed = await readDataframe('tests/Test_smiles_malformed.csv');
    await grok.data.detectSemanticTypes(malformed);
    coreEmpty = await findMCS('smiles', empty);
    coreMalformed = await findMCS('canonical_smiles', malformed);
  });

  test('mcs', async () => {
    const mcs = await grok.functions.call('Chem:FindMCS', {molecules: 'smiles', df: t, returnSmarts: false});
    expect(mcs, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1');
  });

  test('rgroups.smiles', async () => {
    const rgroups: DG.DataFrame = await grok.functions.call('Chem:FindRGroups', {
      molecules: 'smiles',
      df: t,
      core: 'c1ccccc1',
      prefix: 'R',
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

  test('rgroups.molV2000', async () => {
    const df = DG.DataFrame.fromColumns([(await readDataframe('tests/spgi-100.csv')).getCol('Structure')]);
    await grok.functions.call('Chem:FindRGroups', {
      molecules: 'Structure',
      df: df,
      core: `
Actelion Java MolfileCreator 1.0

  5  4  0  0  0  0  0  0  0  0999 V2000
    11.6207  -11.5938    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    12.9198  -10.8438    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    14.2188  -11.5938    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    15.5177  -10.8438    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    16.8168  -11.5938    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
M  END
      `,
      prefix: 'R',
    });
  });

  test('rgroups.molV3000', async () => {
    const df = await readDataframe('tests/approved-drugs-100.csv');
    await grok.functions.call('Chem:FindRGroups', {
      molecules: 'molecule',
      df: df,
      core: `
Actelion Java MolfileCreator 1.0

  2  1  0  0  0  0  0  0  0  0999 V2000
    13.5692  -11.3437    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    14.8683  -12.0938    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
      `,
      prefix: 'R',
    });
  });

  test('rgroups.emptyValues', async () => {
    const res = await findRGroups('smiles', empty, coreEmpty, 'R');
    expect(res.getCol('R1').stats.valueCount, 13);
    expect(res.getCol('R2').stats.valueCount, 13);
  });

  test('rgroups.emptyInput', async () => {
    await findRGroups('smiles', empty, '', 'R');
  });

  test('rgroups.malformedData', async () => {
    const res = await findRGroups('canonical_smiles', malformed, coreMalformed, 'R');
    expect(res.getCol('R1').stats.valueCount, 32);
    expect(res.getCol('R2').stats.valueCount, 9);
    expect(res.getCol('R3').stats.valueCount, 11);
    expect(res.getCol('R4').stats.valueCount, 13);
    expect(res.getCol('R5').stats.valueCount, 3);
  });

  test('rgroups.malformedInput', async () => {
    await findRGroups('canonical_smiles', malformed, malformed.getCol('canonical_smiles').get(2), 'R');
  }, {skipReason: '#1491'});

  after(async () => {
    grok.shell.closeAll();
  });
});

const t = DG.DataFrame.fromCsv(`smiles
O=C1CN=C(c2ccccc2N1)C3CCCCC3
CN1C(=O)CN=C(c2ccccc12)C3CCCCC3
CCCCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
CC(C)CCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
O=C1CN=C(c2ccccc2N1CC3CCCCC3)C4CCCCC4
O=C1CN=C(c2cc(Cl)ccc2N1)C3CCCCC3
CN1C(=O)CN=C(c2cc(Cl)ccc12)C3CCCCC3`);
