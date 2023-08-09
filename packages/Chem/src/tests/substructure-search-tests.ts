import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {readDataframe, _testSearchSubstructure, _testSearchSubstructureAllParameters,
  _testSearchSubstructureSARSmall} from './utils';
import {before, category, test} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {chemSubstructureSearchLibrary} from '../chem-searches';

export const testSubstructure = 'C1CC1';
//csv with C1CC1 as substructure in pos 0 and 2
export const testCsv = `smiles
CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C4CC4
COc1ccc2cc(ccc2c1)C(C)C(=O)Oc3ccc(C)cc3OC
Fc1cc2C(=O)C(=CN(C3CC3)c2cc1N4CCNCC4)c5oc(COc6ccccc6)nn5
CC(C(=O)NCCS)c1cccc(c1)C(=O)c2ccccc2
COc1ccc2c(c1)c(CC(=O)N3CCCC3C(=O)Oc4ccc(C)cc4OC)c(C)n2C(=O)c5ccc(Cl)cc5`;

category('substructure search', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('searchSubstructure.sar-small_awaitAll', async () => {
    await _testSearchSubstructureAllParameters(
      _testSearchSubstructureSARSmall);
  });

  test('searchSubstructureColWithoutDf_awaitAll', async () => {
    const col = DG.Column.fromStrings('smiles', [
      'CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C4CC4',
      'Fc1cc2C(=O)C(=CN(C3CC3)c2cc1N4CCNCC4)c5oc(COc6ccccc6)nn5',
    ]);
    await grok.chem.searchSubstructure(col, 'C1CC1');
  });

  //Number of molecules is smaller than a number of threads
  test('searchSubstructure.5_rows_awaitAll', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [0, 2]);
  });

  test('searchSubstructureWithMalformedMolString_awaitAll', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    df.columns.byName('smiles').set(4, 'qq');
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [0, 2]);
  });

  test('searchSubstructureWithNull_awaitAll', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    df.columns.byName('smiles').set(4, null);
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [0, 2]);
  });

  test('searchSubstructurePerformance_awaitAll', async () => {
    const df = grok.data.demo.molecules(50000);
    const startTime = performance.now();
    await grok.chem.searchSubstructure(df.columns.byName('smiles'), 'c1ccccc1')!;
    const midTime1 = performance.now();
    for (let i = 0; i < 20; i++)
      await grok.chem.searchSubstructure(df.columns.byName('smiles'), 'c1ccccc1')!;
    const midTime2 = performance.now();
    await grok.chem.searchSubstructure(df.columns.byName('smiles'), 'C1CC1')!;
    const endTime = performance.now();

    console.log(`first Call to searchSubstructure took ${midTime1 - startTime} milliseconds`);
    console.log(`loop Call to searchSubstructure took ${midTime2 - midTime1} milliseconds`);
    console.log(`last Call to searchSubstructure took ${endTime - midTime2} milliseconds`);
  }, {timeout: 90000});

  test('searchSubstructureAfterChanges_awaitAll', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [0, 2]);
    df.columns.byName('smiles').set(0, 'COc1ccc2cc(ccc2c1)C(C)C(=O)Oc3ccc(C)cc3OC');
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [2]);
  });

  test('searchSubstructureWith2Columns_awaitAll', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    const col = DG.Column.fromStrings('smiles1', [
      'COc1ccc2cc(ccc2c1)C(C)C(=O)Oc3ccc(C)cc3OC',
      'CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C4CC4',
      'COc1ccc2c(c1)c(CC(=O)N3CCCC3C(=O)Oc4ccc(C)cc4OC)c(C)n2C(=O)c5ccc(Cl)cc5',
      'COc1ccc2c(c1)c(CC(=O)N3CCCC3C(=O)Oc4ccc(C)cc4OC)c(C)n2C(=O)c5ccc(Cl)cc5',
      'Fc1cc2C(=O)C(=CN(C3CC3)c2cc1N4CCNCC4)c5oc(COc6ccccc6)nn5',
    ]);
    df.columns.add(col);
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [0, 2]);
    await _testSearchSubstructure(df, 'smiles1', testSubstructure, [1, 4]);
  });

  test('searchSubstructure2Dataframes_awaitAll', async () => {
    const df1 = grok.data.demo.molecules(50);
    const df2 = grok.data.demo.molecules(50);

    await _testSearchSubstructure(df1, 'smiles', 'c1ccncc1', [6, 26, 46]);
    await _testSearchSubstructure(df2, 'smiles', 'c1ccccc1',
      [0, 4, 5, 7, 8, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 24,
        25, 27, 28, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 44, 45, 47, 48]);
  });

  test('searchSubstructureExplicitHydrogen_awaitAll', async () => {
    const df = await readDataframe('tests/explicit_h_test.csv');
    await _testSearchSubstructure(df, 'smiles', `
  MJ201900                      

  5  5  0  0  0  0  0  0  0  0999 V2000
    -0.1562   -0.7750    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    -0.8707   -0.3625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -0.8707    0.4625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5582    0.4625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5582   -0.3625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  5  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
M  END`, [0, 1, 2, 3]);

    await _testSearchSubstructure(df, 'smiles', `
  MJ201900                      

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.5357    0.1401    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    -0.1787    0.5526    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -0.1787    1.3776    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2501    1.3776    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2501    0.5526    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5357   -0.6848    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  5  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  1  6  1  0  0  0  0
M  END
  `, [2, 3]);
  });


  test('searchSubstructureAromaticBond_awaitAll', async () => {
    const df = await readDataframe('tests/aromatic_bond_test.csv');
    await _testSearchSubstructure(df, 'smiles', `
  MJ201900                      

  6  6  0  0  0  0  0  0  0  0999 V2000
    -1.3616    0.4008    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -2.0761   -0.0116    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -2.0761   -0.8367    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -1.3616   -1.2491    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -0.6471   -0.8367    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    -0.6471   -0.0116    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  6  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  4  0  0  0  0
M  END
`, [2, 3]);
  });


  test('search_benzene_awaitAll', async () => {
    const df = DG.Test.isInBenchmark ? await grok.data.files.openTable('Demo:Files/chem/smiles_1M.zip') :
      grok.data.demo.molecules(500);
    await performanceTestWithConsoleLog(df.col('smiles')!, 'c1ccccc1');
  });

  test('search_2_benzene_rings_awaitAll', async () => {
    const df = DG.Test.isInBenchmark ? await grok.data.files.openTable('Demo:Files/chem/smiles_1M.zip') :
      grok.data.demo.molecules(500);
    await performanceTestWithConsoleLog(df.col('smiles')!, 'c1ccc2ccccc2c1');
  });

  test('search_complex_structure_awaitAll', async () => {
    const df = DG.Test.isInBenchmark ? await grok.data.files.openTable('Demo:Files/chem/smiles_1M.zip') :
      grok.data.demo.molecules(500);
    await performanceTestWithConsoleLog(df.col('smiles')!, 'c1nn2cnnc2s1');
  });

  test('search_non_ring_structure_awaitAll', async () => {
    const df = DG.Test.isInBenchmark ? await grok.data.files.openTable('Demo:Files/chem/smiles_1M.zip') :
      grok.data.demo.molecules(500);
    await performanceTestWithConsoleLog(df.col('smiles')!, 'CNC(C)=O');
  });

  async function performanceTestWithConsoleLog(molCol: DG.Column, query: string) {
    const startTime = performance.now();
    await chemSubstructureSearchLibrary(molCol, query, '');
    const midTime1 = performance.now();
    await chemSubstructureSearchLibrary(molCol, query, '');
    const midTime2 = performance.now();

    console.log(`first Call to substructure search took ${midTime1 - startTime} milliseconds`);
    console.log(`second Call to substructure search took ${midTime2 - midTime1} milliseconds`);
  }
});
