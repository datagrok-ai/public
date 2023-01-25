import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {category, expect, expectFloat, test} from '@datagrok-libraries/utils/src/test';

import {readDataframe, _testSearchSubstructure, _testSearchSubstructureAllParameters, _testSearchSubstructureSARSmall} from './utils';

export const testSubstructure = 'C1CC1';
//csv with C1CC1 as substructure in pos 0 and 2
export const testCsv = `smiles
CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C4CC4
COc1ccc2cc(ccc2c1)C(C)C(=O)Oc3ccc(C)cc3OC
Fc1cc2C(=O)C(=CN(C3CC3)c2cc1N4CCNCC4)c5oc(COc6ccccc6)nn5
CC(C(=O)NCCS)c1cccc(c1)C(=O)c2ccccc2
COc1ccc2c(c1)c(CC(=O)N3CCCC3C(=O)Oc4ccc(C)cc4OC)c(C)n2C(=O)c5ccc(Cl)cc5`;

category('substructure search', () => {
  test('searchSubstructure.sar-small', async () => {
    await _testSearchSubstructureAllParameters(
      _testSearchSubstructureSARSmall);
  });

  test('searchSubstructureColWithoutDf', async () => {
    const col = DG.Column.fromStrings('smiles', [
      'CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C4CC4',
      'Fc1cc2C(=O)C(=CN(C3CC3)c2cc1N4CCNCC4)c5oc(COc6ccccc6)nn5',
    ]);
    await grok.chem.searchSubstructure(col, 'C1CC1');
  });

  //Number of molecules is smaller than a number of threads
  test('searchSubstructure.5_rows', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [0, 2]);
  });

  test('searchSubstructureWithMalformedMolString', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    df.columns.byName('smiles').set(4, 'qq');
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [0, 2]);
  });

  test('searchSubstructureWithNull', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    df.columns.byName('smiles').set(4, null);
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [0, 2]);
  });

  test('searchSubstructurePerformance', async () => {
    const df = grok.data.demo.molecules(50000);

    const startTime = performance.now();
    await grok.chem.searchSubstructure(df.columns.byName('smiles'), 'c1ccccc1')!;
    const midTime1 = performance.now();
    for (let i = 0; i < 20; i++)
      await grok.chem.searchSubstructure(df.columns.byName('smiles'), 'c1ccccc1')!;
    const midTime2 = performance.now();
    await grok.chem.searchSubstructure(df.columns.byName('smiles'), 'C1CC1')!;
    const endTime = performance.now();

    console.log(`first Call to WithLibrary took ${midTime1 - startTime} milliseconds`);
    console.log(`loop Call to WithLibrary took ${midTime2 - midTime1} milliseconds`);
    console.log(`last Call to WithLibrary took ${endTime - midTime2} milliseconds`);
  });

  test('searchSubstructureAfterChanges', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [0, 2]);
    df.columns.byName('smiles').set(0, 'COc1ccc2cc(ccc2c1)C(C)C(=O)Oc3ccc(C)cc3OC');
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [2]);
  });

  test('searchSubstructureWith2Columns', async () => {
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

  test('searchSubstructureExplicitHydrogen', async () => {
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


  test('searchSubstructureAromaticBond', async () => {
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
});
