import {category, expect, test} from "@datagrok-libraries/utils/src/test";
import {_testSearchSubstructure, _testSearchSubstructureAllParameters, _testSearchSubstructureSARSmall, requireText} from "./utils";
import * as DG from "datagrok-api/dg";
import * as grok from "datagrok-api/grok";
import * as ui from "datagrok-api/ui";

category('chem', () => {

  test('chem.searchSubstructure.sar_small', async () => {
    await _testSearchSubstructureAllParameters(
      _testSearchSubstructureSARSmall);
  });

// Number of molecules is smaller than a number of threads
  test('chem.searchSubstructure.5_rows', async () => {
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

  test('chem.findSimilar.sar_small', async () => {
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
    const areEqualFloat = (a: number, b: number) => Math.abs(a - b) < 0.001;
    const arr = first5Rows;
    expect(arr[0].molecule, 'O=C1CN=C(c2ccccc2N1)C3CCCCC3');
    expect(arr[1].molecule, 'O=C1CN=C(c2cc(I)ccc2N1)C3CCCCC3');
    expect(arr[2].molecule, 'O=C1CN=C(c2cc(Cl)ccc2N1)C3CCCCC3');
    expect(arr[3].molecule, 'O=C1CN=C(c2cc(F)ccc2N1)C3CCCCC3');
    expect(arr[4].molecule, 'O=C1CN=C(c2cc(Br)ccc2N1)C3CCCCC3');
    expect(areEqualFloat(arr[0].score, 1.0000), true);
    expect(areEqualFloat(arr[1].score, 0.6905), true);
    expect(areEqualFloat(arr[2].score, 0.6744), true);
    expect(areEqualFloat(arr[3].score, 0.6744), true);
    expect(areEqualFloat(arr[4].score, 0.6744), true);
    expect(arr[0].index, 0);
    expect(arr[1].index, 30);
    expect(arr[2].index, 5);
    expect(arr[3].index, 20);
    expect(arr[4].index, 25);
  });


  test('testSubstructureSearch', async () => {
    let t = grok.data.demo.molecules();
    await grok.chem.searchSubstructure(t.col('smiles')!, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1');
  });

  test('testDescriptors', async () => {
    let t = grok.data.demo.molecules();
    let descriptors = await grok.chem.descriptors(t, 'smiles', ['MolWt', 'Lipinski']);
  });

  test('testDiversitySearch', async () => {
    let t = grok.data.demo.molecules();
    await grok.chem.diversitySearch(t.col('smiles')!);
  });

  test('testSimilaritySearch', async () => {
    let t = grok.data.demo.molecules();
    let scores = await grok.chem.getSimilarities(t.col('smiles')!, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1');
  });

  test('testMcs', async () => {
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

  test('testRendering', async () => {
    ui.dialog()
      .add(grok.chem.svgMol('O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1'))
      .show();
  });

  test('testRenderingCanvas', async () => {
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


});
