import {category, expect, expectFloat, test} from '@datagrok-libraries/utils/src/test';
import {_testSearchSubstructure,
  _testSearchSubstructureAllParameters,
  _testSearchSubstructureSARSmall, loadFileAsText} from './utils';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {_testDiversitySearchViewerOpen, _testSimilaritySearchFunctionality,
  _testSimilaritySearchViewerOpen} from './similarity-diversity-utils';

//csv with C1CC1 as substructure in pos 0 and 2
const testCsv = `smiles
CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C4CC4
COc1ccc2cc(ccc2c1)C(C)C(=O)Oc3ccc(C)cc3OC
Fc1cc2C(=O)C(=CN(C3CC3)c2cc1N4CCNCC4)c5oc(COc6ccccc6)nn5
CC(C(=O)NCCS)c1cccc(c1)C(=O)c2ccccc2
COc1ccc2c(c1)c(CC(=O)N3CCCC3C(=O)Oc4ccc(C)cc4OC)c(C)n2C(=O)c5ccc(Cl)cc5`;

const testSubstructure = 'C1CC1';

category('chem', () => {
  test('searchSubstructure.sar-small', async () => {
    await _testSearchSubstructureAllParameters(
      _testSearchSubstructureSARSmall);
  });

  test('searchSustructureColWithoutDf', async () => {
    const col = DG.Column.fromStrings('smiles', [
      'CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C4CC4',
      'Fc1cc2C(=O)C(=CN(C3CC3)c2cc1N4CCNCC4)c5oc(COc6ccccc6)nn5',
    ]);
    await grok.chem.searchSubstructure(col, 'C1CC1');
  });

  //Number of molecules is smaller than a number of threads
  test('searchSustructure.5_rows', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [0, 2]);
  });

  test('searchSustructureWithMalformedMolString', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    df.columns.byName('smiles').set(4, 'qq');
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [0, 2]);
  });

  test('searchSustructureWithNull', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    df.columns.byName('smiles').set(4, null);
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [0, 2]);
  });

  test('searchSustructurePerformance', async () => {
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

  test('searchSustructureAfterChanges', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [0, 2]);
    df.columns.byName('smiles').set(0, 'COc1ccc2cc(ccc2c1)C(C)C(=O)Oc3ccc(C)cc3OC');
    await _testSearchSubstructure(df, 'smiles', testSubstructure, [2]);
  });

  test('searchSustructureWith2Columns', async () => {
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

  test('findSimilar.sar-small', async () => {
    const dfInput = DG.DataFrame.fromCsv(await loadFileAsText('sar-small.csv'));
    const colInput = dfInput.columns.byIndex(0);
    const dfResult: DG.DataFrame = // shouldn't be null
      (await grok.chem.findSimilar(colInput, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1'))!;
    const numRowsOriginal = dfInput.rowCount;
    const numRows = dfResult.rowCount;
    const columnNames = [
      dfResult.columns.byIndex(0).name,
      dfResult.columns.byIndex(1).name,
      dfResult.columns.byIndex(2).name,
    ];
    const first5Rows: any[] = [];
    for (let i = 0; i < 5; ++i) {
      const molecule: string = dfResult.columns.byIndex(0).get(i);
      const score: number = dfResult.columns.byIndex(1).get(i);
      const index: number = dfResult.columns.byIndex(2).get(i);
      first5Rows[i] = {molecule, score, index};
    }
    expect(numRows, numRowsOriginal);
    expect(columnNames[0], 'molecule');
    expect(columnNames[1], 'score');
    expect(columnNames[2], 'index');
    const arr = first5Rows;
    expect(arr[0].molecule, 'O=C1CN=C(c2ccccc2N1)C3CCCCC3');
    expect(arr[1].molecule, 'O=C1CN=C(c2cc(I)ccc2N1)C3CCCCC3');
    expect(arr[2].molecule, 'O=C1CN=C(c2cc(Cl)ccc2N1)C3CCCCC3');
    expect(arr[3].molecule, 'O=C1CN=C(c2cc(F)ccc2N1)C3CCCCC3');
    expect(arr[4].molecule, 'O=C1CN=C(c2cc(Br)ccc2N1)C3CCCCC3');
    expectFloat(arr[0].score, 1.0000);
    expectFloat(arr[1].score, 0.6905);
    expectFloat(arr[2].score, 0.6744);
    expectFloat(arr[3].score, 0.6744);
    expectFloat(arr[4].score, 0.6744);
    expect(arr[0].index, 0);
    expect(arr[1].index, 30);
    expect(arr[2].index, 5);
    expect(arr[3].index, 20);
    expect(arr[4].index, 25);
  });

  test('getSimilarities.molecules', async () => {
    const df = grok.data.demo.molecules();
    const scores = (await grok.chem.getSimilarities(df.columns.byName('smiles'),
      'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1'))!;
    expectFloat(scores.get(0), 0.1034);
    expectFloat(scores.get(1), 0.07407);
    expectFloat(scores.get(2), 0.11111);
    expectFloat(scores.get(3), 0.11111);
    expectFloat(scores.get(4), 0.07042);
    expectFloat(scores.get(5), 0.06349);
  });

  test('testSubstructureSearch', async () => {
    const t = grok.data.demo.molecules();
    await grok.chem.searchSubstructure(t.col('smiles')!, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1');
  });

  test('testDescriptors', async () => {
    const t = grok.data.demo.molecules();
    await grok.chem.descriptors(t, 'smiles', ['MolWt', 'Lipinski']);
  });

  test('testDiversitySearch', async () => {
    const t = grok.data.demo.molecules();
    await grok.chem.diversitySearch(t.col('smiles')!);
  });

  test('testMcs', async () => {
    const t = DG.DataFrame.fromCsv(`smiles
O=C1CN=C(c2ccccc2N1)C3CCCCC3
CN1C(=O)CN=C(c2ccccc12)C3CCCCC3
CCCCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
CC(C)CCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
O=C1CN=C(c2ccccc2N1CC3CCCCC3)C4CCCCC4
O=C1CN=C(c2cc(Cl)ccc2N1)C3CCCCC3
CN1C(=O)CN=C(c2cc(Cl)ccc12)C3CCCCC3`);
    await grok.chem.mcs(t.col('smiles')!);
  });

  test('testMcsPanel', async () => {
    const t = DG.DataFrame.fromCsv(`smiles
      O=C1CN=C(c2ccccc2N1)C3CCCCC3
      CN1C(=O)CN=C(c2ccccc12)C3CCCCC3
      CCCCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
      CC(C)CCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
      O=C1CN=C(c2ccccc2N1CC3CCCCC3)C4CCCCC4
      O=C1CN=C(c2cc(Cl)ccc2N1)C3CCCCC3
      CN1C(=O)CN=C(c2cc(Cl)ccc12)C3CCCCC3`);
    const v = grok.shell.addTableView(t);
    await grok.functions.call('Chem:addMcsPanel', {col: t.columns.byName('smiles')});
    v.close();
  });

  test('testInchiPanel', async () => {
    const t = DG.DataFrame.fromCsv(`smiles
      COc1ccc(CN(Cc2ccccc2)Cc3ccc(Br)cc3)cc1O
      COc1ccc(CN(CCc2ccc(Br)cc2)Cc3ccccc3)cc1O
      CCCCCCCC(=O)NCCC1CC(CC)(CC)C(=O)O1
      CC1=C(C(C(=C(C)N1)C(=O)OCc2ccccc2)c3csc(n3)c4ccc(Cl)cc4)C(=O)
      CCCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCC)
      CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC(C)C)
      CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCC(C)C)`);
    const v = grok.shell.addTableView(t);
    await grok.functions.call('Chem:addInchisPanel', {col: t.columns.byName('smiles')});
    v.close();
  });

  test('testInchiKeysPanel', async () => {
    const t = DG.DataFrame.fromCsv(`smiles
      COc1ccc(CN(Cc2ccccc2)Cc3ccc(Br)cc3)cc1O
      COc1ccc(CN(CCc2ccc(Br)cc2)Cc3ccccc3)cc1O
      CCCCCCCC(=O)NCCC1CC(CC)(CC)C(=O)O1
      CC1=C(C(C(=C(C)N1)C(=O)OCc2ccccc2)c3csc(n3)c4ccc(Cl)cc4)C(=O)
      CCCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCC)
      CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC(C)C)
      CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCC(C)C)`);
    const v = grok.shell.addTableView(t);
    await grok.functions.call('Chem:addInchisKeysPanel', {col: t.columns.byName('smiles')});
    v.close();
  });

  test('testCurateTopMenu', async () => {
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
    await grok.functions.call('Chem:CurateChemStructures', {'data': t, 'smiles': 'smiles'});
    v.close();
  });

  test('testRendering', async () => {
    ui.dialog()
      .add(grok.chem.svgMol('O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1'))
      .show();
  });

  test('testRenderingCanvas', async () => {
    const root = ui.div();
    const width = 300;
    const height = 200;
    const canvas1 = ui.canvas(width, height);
    canvas1.id = 'canvas1';
    const canvas2 = ui.canvas(width, height);
    canvas2.id = 'canvas2';
    root.appendChild(canvas1);
    root.appendChild(canvas2);
    // TODO: should be syncronous
    await grok.chem.canvasMol(
      0, 0, width, height, canvas1,
      'COc1ccc2cc(ccc2c1)C(C)C(=O)OCCCc3cccnc3',
      'c1ccccc1');
    await grok.chem.canvasMol(
      0, 0, width, height, canvas2,
      'CN1CCC(O)(CC1)c2ccccc2');
    ui
      .dialog({title: 'Test molecule'})
      .add(root)
      .show();
  });

  test('similaritySearchViewerOpen', async () => {
    await _testSimilaritySearchViewerOpen();
  });

  test('similaritySearchFunctionality', async () => {
    await _testSimilaritySearchFunctionality('Tanimoto', 'Morgan');
    await _testSimilaritySearchFunctionality('Dice', 'Morgan');
    await _testSimilaritySearchFunctionality('Cosine', 'Morgan');
    await _testSimilaritySearchFunctionality('Euclidean', 'Morgan');
    await _testSimilaritySearchFunctionality('Hamming', 'Morgan');
  });

  test('diversitySearchViewerOpen', async () => {
    await _testDiversitySearchViewerOpen();
  });
});
