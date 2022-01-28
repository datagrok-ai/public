import {category, expect, expectFloat, test, testExpectFinish} from "@datagrok-libraries/utils/src/test";
import {_testSearchSubstructure, _testSearchSubstructureAllParameters, _testSearchSubstructureSARSmall, requireText} from "./utils";
import * as DG from "datagrok-api/dg";
import * as grok from "datagrok-api/grok";
import * as ui from "datagrok-api/ui";

category('chem', () => {

  testExpectFinish('searchSubstructure.sar_small', async () => {
    await _testSearchSubstructureAllParameters(
      _testSearchSubstructureSARSmall);
  });

// Number of molecules is smaller than a number of threads
  testExpectFinish('searchSustructure.5_rows', async () => {
    const targetSmiles = [
      'CCOC(=O)c1oc2cccc(OCCNCc3cccnc3)c2c1C4CC4',
      'Fc1cc2C(=O)C(=CN(C3CC3)c2cc1N4CCNCC4)c5oc(COc6ccccc6)nn5',
    ];
    const substructureSmiles = 'C1CC1';
    const csv =
      `smiles
${targetSmiles[0]}
COc1ccc2cc(ccc2c1)C(C)C(=O)Oc3ccc(C)cc3OC
${targetSmiles[1]}
CC(C(=O)NCCS)c1cccc(c1)C(=O)c2ccccc2
COc1ccc2c(c1)c(CC(=O)N3CCCC3C(=O)Oc4ccc(C)cc4OC)c(C)n2C(=O)c5ccc(Cl)cc5
`;
    await _testSearchSubstructureAllParameters(
      async (params: any) => await _testSearchSubstructure(
        csv, substructureSmiles, [0, 2], params,
      ),
    );
  });

  testExpectFinish('findSimilar.sar_small', async () => {
    const dfInput = DG.DataFrame.fromCsv(await requireText('sar_small.csv'));
    const colInput = dfInput.columns[0];
    const dfResult: DG.DataFrame = // shouldn't be null
      (await grok.chem.findSimilar(colInput, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1'))!;
    const numRowsOriginal = dfInput.rowCount;
    const numRows = dfResult.rowCount;
    const columnNames = [
      dfResult.columns[0].name,
      dfResult.columns[1].name,
      dfResult.columns[2].name,
    ];
    const first5Rows: any[] = [];
    for (let i = 0; i < 5; ++i) {
      const molecule: string = dfResult.columns[0].get(i);
      const score: number = dfResult.columns[1].get(i);
      const index: number = dfResult.columns[2].get(i);
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

  testExpectFinish('getSimilarities.molecules', async () => {
    let df = grok.data.demo.molecules();
    const scores = (await grok.chem.getSimilarities(df.columns['smiles'], 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1'))!;
    expectFloat(scores.get(0), 0.1034);
    expectFloat(scores.get(1), 0.07407);
    expectFloat(scores.get(2), 0.11111);
    expectFloat(scores.get(3), 0.11111);
    expectFloat(scores.get(4), 0.07042);
    expectFloat(scores.get(5), 0.06349);
  });

  testExpectFinish('testSubstructureSearch', async () => {
    let t = grok.data.demo.molecules();
    await grok.chem.searchSubstructure(t.col('smiles')!, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1');
  });

  testExpectFinish('testDescriptors', async () => {
    let t = grok.data.demo.molecules();
    let descriptors = await grok.chem.descriptors(t, 'smiles', ['MolWt', 'Lipinski']);
  });

  testExpectFinish('testDiversitySearch', async () => {
    let t = grok.data.demo.molecules();
    await grok.chem.diversitySearch(t.col('smiles')!);
  });

  testExpectFinish('testMcs', async () => {
    let t = DG.DataFrame.fromCsv(`smiles
O=C1CN=C(c2ccccc2N1)C3CCCCC3
CN1C(=O)CN=C(c2ccccc12)C3CCCCC3
CCCCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
CC(C)CCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
O=C1CN=C(c2ccccc2N1CC3CCCCC3)C4CCCCC4
O=C1CN=C(c2cc(Cl)ccc2N1)C3CCCCC3
CN1C(=O)CN=C(c2cc(Cl)ccc12)C3CCCCC3`);
    await grok.chem.mcs(t.col('smiles')!);
  });

  testExpectFinish('testMcsPanel', async () => {
    let t = DG.DataFrame.fromCsv(`smiles
      O=C1CN=C(c2ccccc2N1)C3CCCCC3
      CN1C(=O)CN=C(c2ccccc12)C3CCCCC3
      CCCCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
      CC(C)CCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
      O=C1CN=C(c2ccccc2N1CC3CCCCC3)C4CCCCC4
      O=C1CN=C(c2cc(Cl)ccc2N1)C3CCCCC3
      CN1C(=O)CN=C(c2cc(Cl)ccc12)C3CCCCC3`);
    let v = grok.shell.addTableView(t);
    await grok.functions.call('Chem:addMcsPanel', {col: t.columns['smiles']});
    v.close();
  });

  testExpectFinish('testInchiPanel', async () => {
    let t = DG.DataFrame.fromCsv(`smiles
      COc1ccc(CN(Cc2ccccc2)Cc3ccc(Br)cc3)cc1O
      COc1ccc(CN(CCc2ccc(Br)cc2)Cc3ccccc3)cc1O
      CCCCCCCC(=O)NCCC1CC(CC)(CC)C(=O)O1
      CC1=C(C(C(=C(C)N1)C(=O)OCc2ccccc2)c3csc(n3)c4ccc(Cl)cc4)C(=O)
      CCCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCC)
      CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC(C)C)
      CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCC(C)C)`);
    let v = grok.shell.addTableView(t);
    await grok.functions.call('Chem:addInchisPanel', {col: t.columns['smiles']});
    v.close();
  });

  testExpectFinish('testInchiKeysPanel', async () => {
    let t = DG.DataFrame.fromCsv(`smiles
      COc1ccc(CN(Cc2ccccc2)Cc3ccc(Br)cc3)cc1O
      COc1ccc(CN(CCc2ccc(Br)cc2)Cc3ccccc3)cc1O
      CCCCCCCC(=O)NCCC1CC(CC)(CC)C(=O)O1
      CC1=C(C(C(=C(C)N1)C(=O)OCc2ccccc2)c3csc(n3)c4ccc(Cl)cc4)C(=O)
      CCCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCC)
      CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC(C)C)
      CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCC(C)C)`);
    let v = grok.shell.addTableView(t);
    await grok.functions.call('Chem:addInchisKeysPanel', {col: t.columns['smiles']});
    v.close();
  });

  testExpectFinish('testCurateTopMenu', async () => {
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
    let v = grok.shell.addTableView(t);
    await grok.functions.call('Chem:CurateChemStructures', {'data': t, 'smiles': 'smiles'})
    v.close();
  });

  testExpectFinish('testRendering', async () => {
    ui.dialog()
      .add(grok.chem.svgMol('O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1'))
      .show();
  });

  testExpectFinish('testRenderingCanvas', async () => {
    let root = ui.div();
    const width = 300;
    const height = 200;
    let canvas1 = ui.canvas(width, height);
    canvas1.id = 'canvas1';
    let canvas2 = ui.canvas(width, height);
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

  testExpectFinish('testSimilaritySearchViewer', async () => {
    let table = grok.data.demo.molecules(1000);
    table.selection.set(0, true);
    let view = grok.shell.addTableView(table);
    const timer = (ms: number) => new Promise(res => setTimeout(res, ms));
    await timer(5000);
    let similaritySearchviewer = view.addViewer('SimilaritySearchViewer');
    table.currentRowIdx = 5;
    table.selection.set(2, true);
    table.selection.set(6, true);
    table.selection.set(5, true);
    table.rows.removeAt(2, 10);
    similaritySearchviewer.setOptions({
      distanceMetric: 'dice',
      limit: 100,
      minScore: 0.5,
      moleculeColumnName: 'smiles',
    });
    await timer(5000);
  });

  testExpectFinish('testDiversitySearchViewer', async () => {
    let table = grok.data.demo.molecules(1000);
    table.selection.set(0, true);
    let view = grok.shell.addTableView(table);
    let similaritySearchviewer = view.addViewer('DiversitySearchViewer');
    table.currentRowIdx = 5;
    table.selection.set(2, true);
    table.selection.set(6, true);
    table.selection.set(5, true);
    table.rows.removeAt(2, 100);
    similaritySearchviewer.setOptions({
      distanceMetric: 'dice',
      limit: 100,
      moleculeColumnName: 'smiles',
    });
  });
});
