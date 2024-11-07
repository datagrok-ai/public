import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import { category, expect, test, before, after } from '@datagrok-libraries/utils/src/test';
import { _package } from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import { readDataframe } from './utils';
import { getMCS } from '../utils/most-common-subs';
import { rGroupsMinilib } from '../analysis/r-group-analysis';


category('top menu r-groups', () => {
  let empty: DG.DataFrame;
  let malformed: DG.DataFrame;
  let dfForMcs: DG.DataFrame;
  let coreEmpty: string;
  let coreMalformed: string;
  const rGroupOpts = {
    matchingStrategy: 'Greedy',
    includeTargetMolInResults: true,
    onlyMatchAtRGroups: false,
  };

  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    empty = await readDataframe('tests/sar-small_empty_vals.csv');
    await grok.data.detectSemanticTypes(empty);
    malformed = await readDataframe('tests/sar-small_test_malformed.csv');
    await grok.data.detectSemanticTypes(malformed);
    coreEmpty = await getMCS(empty.col('smiles')!, true, true)!;
    coreMalformed = await getMCS(malformed.col('smiles')!, true, true)!;
    dfForMcs = DG.Test.isInBenchmark ? await readDataframe('smiles.csv') :
      await readDataframe('tests/spgi-100.csv');
  });

  test('mcs.exactAtomsExactBonds', async () => {
    const mcs = await getMCS(dfForMcs.col(DG.Test.isInBenchmark ? `canonical_smiles` : `Structure`)!, true, true);
    expect(mcs, DG.Test.isInBenchmark ? '[#6]' : '[#6]-[#6]-[#7]-[#6]');
  }, {benchmark: true});

  test('mcs.anyAtomsExactBonds', async () => {
    const mcs = await getMCS(dfForMcs.col(DG.Test.isInBenchmark ? `canonical_smiles` : `Structure`)!, false, true);
    expect(mcs, DG.Test.isInBenchmark ? '[#7,#6,#8,#9,#16,#17,#35,#53]-[#6,#7,#8,#16,#53]' : '[#6,#7,#8,#9]-[#6,#7,#8]-[#7,#6](-[#6,#8,#9,#16])-[#6,#7,#9]');
  }, {benchmark: true});

  test('mcs.exactAtomsAnyBonds', async () => {
    const mcs = await getMCS(dfForMcs.col(DG.Test.isInBenchmark ? `canonical_smiles` : `Structure`)!, true, false);
    expect(mcs, DG.Test.isInBenchmark ? '[#6]-,:[#6]' : '[#6]-,:[#6]-,:[#7]-,:[#6]-,:[#6]');
  }, {benchmark: true});

  test('mcs.anyAtomsAnyBonds', async () => {
    const mcs = await getMCS(dfForMcs.col(DG.Test.isInBenchmark ? `canonical_smiles` : `Structure`)!, false, false);
    expect(mcs, DG.Test.isInBenchmark ?
      `[#7,#6,#8,#9,#16,#17,#35,#53]-,:[#6,#7,#8,#16]-,:[#6,#7,#8,#16]:,-[#6,#7,#8,#16]:,-[#7,#6,#8,#16]:,-[#7,#6,#8,#9,#16,#17,#35]` :
      `[#6,#8,#9]-,:[#6,#7,#8]-,:[#7,#6](-,:[#6,#7,#8,#16]-,:[#6,#7,#8]-,:[#6,#7]-,:[#7,#6,#8]-,:[#6,#7,#8,#16]-,:[#6,#7,#8,#16])-,:[#6,#7,#8]-,:[#6,#7]-,:[#6,#7,#8]-,:[#6,#7,#8,#9]`);
  }, {benchmark: true});

  test('rgroups.smiles', async () => {
    if (DG.Test.isInBenchmark) {
      /*       await grok.functions.call('Chem:FindRGroups', {
              molecules: 'canonical_smiles',
              df: await grok.data.files.openTable('System:AppData/Chem/tests/smiles_200K.zip'),
              core: 'c1ccccc1',
              prefix: 'R',
            }); */
      const df = await readDataframe('smiles.csv');
      await rGroupsMinilib(df.col('canonical_smiles')!, 'c1ccccc1', false, 0, rGroupOpts);
      return;
    }
    /*     const rgroups: DG.DataFrame = await grok.functions.call('Chem:FindRGroups', {
          molecules: 'smiles',
          df: t,
          core: 'c1ccccc1',
          prefix: 'R',
        }); */
    const rgroups = DG.DataFrame.fromColumns((await rGroupsMinilib(t.col('smiles')!,
      'c1ccccc1', false, 0, rGroupOpts)).rGroups);
    if (!DG.Test.isInBenchmark) {
      expect(rgroups.getCol('R1').get(0), 'O=C(C/N=C(/C1CCCCC1)[*:4])N[*:1]');
      expect(rgroups.getCol('R1').get(1), '[H][*:1]');
      expect(rgroups.getCol('R1').get(2), 'CCCCN(C(=O)C/N=C(/C1CCCCC1)[*:1])[*:4]');
      expect(rgroups.getCol('R1').get(3), 'CC(C)CCN(C(=O)C/N=C(/C1CCCCC1)[*:1])[*:4]');
      expect(rgroups.getCol('R1').get(4), 'O=C(C/N=C(/C1CCCCC1)[*:1])N(CC1CCCCC1)[*:4]');
      expect(rgroups.getCol('R1').get(5), 'O=C(C/N=C(/C1CCCCC1)[*:4])N[*:1]');
      expect(rgroups.getCol('R1').get(6), 'CN(C(=O)C/N=C(/C1CCCCC1)[*:4])[*:1]');
      expect(rgroups.getCol('R2').get(5), 'Cl[*:2]');
      expect(rgroups.getCol('R2').get(6), 'Cl[*:2]');
    }
  }, {benchmark: true});

  test('rgroups.molV2000', async () => {
    const df = DG.DataFrame.fromColumns([(await readDataframe('tests/spgi-100.csv')).getCol('Structure')]);
    const core = `
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
`;
    /*     await grok.functions.call('Chem:FindRGroups', {
          molecules: 'Structure',
          df: df,
          core: core,
          prefix: 'R',
        }); */
    await rGroupsMinilib(df.col('Structure')!, core, false, 0, rGroupOpts);
  });

  test('rgroups.molV3000', async () => {
    const df = await readDataframe('tests/approved-drugs-100.csv');
    const core = `
Actelion Java MolfileCreator 1.0

  2  1  0  0  0  0  0  0  0  0999 V2000
   13.5692  -11.3437    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   14.8683  -12.0938    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
`;
    /*     await grok.functions.call('Chem:FindRGroups', {
          molecules: 'molecule',
          df: df,
          core: core,
          prefix: 'R',
        }); */
    await rGroupsMinilib(df.col('molecule')!, core, false, 0, rGroupOpts);
  }, { timeout: 60000 });

  test('rgroups.emptyValues', async () => {
    //const res = await findRGroups('smiles', empty, coreEmpty, 'R');
    const res = DG.DataFrame.fromColumns((await rGroupsMinilib(empty.col('smiles')!, coreEmpty,
      true, 0, rGroupOpts)).rGroups);
    expect(res.getCol('R1').stats.valueCount, 16);
    expect(res.getCol('R2').stats.valueCount, 16);
  }, { timeout: 60000 });

  test('rgroups.emptyInput', async () => {
    //await findRGroups('smiles', empty, '', 'R');
    let exception = false;
    try {
      await rGroupsMinilib(empty.col('smiles')!, '', false, 0, rGroupOpts);
    } catch (e: any) {
      exception = true;
      expect(e.message === 'No core was provided');
    }
    expect(exception, true, 'Exception has not been generated on empty core');
  });

  test('rgroups.malformedData', async () => {
    //const res = await findRGroups('canonical_smiles', malformed, coreMalformed, 'R');
    const res = DG.DataFrame.fromColumns((await rGroupsMinilib(malformed.col('smiles')!,
      coreMalformed, true, 0, rGroupOpts)).rGroups);
    expect(res.getCol('R1').stats.valueCount, 16);
    expect(res.getCol('R2').stats.valueCount, 16);
  });

  test('rgroups.malformedInput', async () => {
    //const res = await findRGroups('canonical_smiles', malformed, malformed.getCol('canonical_smiles').get(2), 'R');
    const malformedInput = malformed.getCol('smiles').get(0);
    let exception = false;
    try {
      const res = DG.DataFrame.fromColumns((await rGroupsMinilib(malformed.col('smiles')!,
        malformedInput, false, 0, rGroupOpts)).rGroups);
      expect(res.rowCount, 0);
    } catch (e: any) {
      exception = true;
      expect(e.message === 'Core is possibly malformed');
    }
    expect(exception, true, 'Exception has not been generated on malformed core');
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
