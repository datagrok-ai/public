import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, expect, test, before, after} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {readDataframe} from './utils';
import {findRGroups} from '../scripts-api';
import {getMCS} from '../utils/most-common-subs';


category('top menu r-groups', () => {
  let empty: DG.DataFrame;
  let malformed: DG.DataFrame;
  let dfForMcs: DG.DataFrame;
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
    coreEmpty = await getMCS(empty.col('smiles')!, true, true)!;
    coreMalformed = await getMCS(malformed.col('canonical_smiles')!, true, true)!;
    dfForMcs = DG.Test.isInBenchmark ? await grok.data.files.openTable('System:AppData/Chem/tests/smiles_50K.csv') :
      await readDataframe('tests/spgi-100.csv');
  });

  test('mcs.exactAtomsExactBonds', async () => {
    const mcs = await getMCS(dfForMcs.col(DG.Test.isInBenchmark ? `smiles` : `Structure`)!, true, true);
    expect(mcs, DG.Test.isInBenchmark ? '[#17]' : '[#6]-[#6]-[#7]-[#6]');
  });

  test('mcs.anyAtomsExactBonds', async () => {
    const mcs = await getMCS(dfForMcs.col(DG.Test.isInBenchmark ? `smiles` : `Structure`)!, false, true);
    expect(mcs, DG.Test.isInBenchmark ? '[#17]' : '[#6,#7,#8,#9]-[#6,#7,#8]-[#7,#6](-[#6,#8,#9,#16])-[#6,#7,#9]');
  });

  test('mcs.exactAtomsAnyBonds', async () => {
    const mcs = await getMCS(dfForMcs.col(DG.Test.isInBenchmark ? `smiles` : `Structure`)!, true, false);
    expect(mcs, DG.Test.isInBenchmark ? '[#17]' : '[#6]-,:[#6]-,:[#7]-,:[#6]-,:[#6]');
  });

  test('mcs.anyAtomsAnyBonds', async () => {
    const mcs = await getMCS(dfForMcs.col(DG.Test.isInBenchmark ? `smiles` : `Structure`)!, false, false);
    expect(mcs, DG.Test.isInBenchmark ?
      `[#8,#6,#7,#9,#15,#16,#17,#35]-,:[#17,#5,#6,#7,#8,#14,#15,#16,#33,#34]` :
      `[#6,#8,#9]-,:[#6,#7,#8]-,:[#7,#6](-,:[#6,#7,#8,#16]-,:[#6,#7,#8]-,:[#6,#7]-,:[#7,#6,#8]-,:[#6,#7,#8,#16]-,:[#6,#7,#8,#16])-,:[#6,#7,#8]-,:[#6,#7]-,:[#6,#7,#8]-,:[#6,#7,#8,#9]`);
  });

  test('rgroups.smiles', async () => {
    if (DG.Test.isInBenchmark) {
      await grok.functions.call('Chem:FindRGroups', {
        molecules: 'canonical_smiles',
        df: await grok.data.files.openTable('System:AppData/Chem/tests/smiles_200K.zip'),
        core: 'c1ccccc1',
        prefix: 'R',
      });
      return;
    }
    const rgroups: DG.DataFrame = await grok.functions.call('Chem:FindRGroups', {
      molecules: 'smiles',
      df: t,
      core: 'c1ccccc1',
      prefix: 'R',
    });
    if (!DG.Test.isInBenchmark) {
      expect(rgroups.getCol('R1').get(0), '*C(=NCC(=O)N[1*])C1CCCCC1');
      expect(rgroups.getCol('R1').get(1), '*C(=NCC(=O)N([1*])C)C1CCCCC1');
      expect(rgroups.getCol('R1').get(2), '*C(=NCC(=O)N([1*])CCCC)C1CCCCC1');
      expect(rgroups.getCol('R1').get(3), '*C(=NCC(=O)N([1*])CCC(C)C)C1CCCCC1');
      expect(rgroups.getCol('R1').get(4), '*C(=NCC(=O)N([1*])CC1CCCCC1)C1CCCCC1');
      expect(rgroups.getCol('R1').get(5), '*C(=NCC(=O)N[1*])C1CCCCC1');
      expect(rgroups.getCol('R1').get(6), '*C(=NCC(=O)N([1*])C)C1CCCCC1');
      expect(rgroups.getCol('R2').get(5), '[2*]Cl');
      expect(rgroups.getCol('R2').get(6), '[2*]Cl');
    }
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
  }, {timeout: 60000});

  test('rgroups.emptyValues', async () => {
    const res = await findRGroups('smiles', empty, coreEmpty, 'R');
    expect(res.getCol('R1').stats.valueCount, 13);
    expect(res.getCol('R2').stats.valueCount, 13);
  }, {timeout: 60000});

  test('rgroups.emptyInput', async () => {
    await findRGroups('smiles', empty, '', 'R');
  });

  test('rgroups.malformedData', async () => {
    const res = await findRGroups('canonical_smiles', malformed, coreMalformed, 'R');
    expect(res.getCol('R1').stats.valueCount, 38);
    expect(res.getCol('R2').stats.valueCount, 1);
    expect(res.getCol('R3').stats.valueCount, 1);
    expect(res.getCol('R4').stats.valueCount, 1);
  });

  test('rgroups.malformedInput', async () => {
    const res = await findRGroups('canonical_smiles', malformed, malformed.getCol('canonical_smiles').get(2), 'R');
    expect(res.rowCount, 0);
  });

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
