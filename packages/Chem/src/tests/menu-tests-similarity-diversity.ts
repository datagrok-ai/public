import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import { category, expect, expectFloat, test, delay } from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {Fingerprint} from '../utils/chem-common';
import { createTableView, readDataframe } from './utils';

import { _testSearchSubstructure,
  _testSearchSubstructureAllParameters,
  _testSearchSubstructureSARSmall,
  loadFileAsText } from './utils';

const testSimilarityResults = {
  'Tanimoto/Morgan': [
    {smiles: 'O=C1CN=C(c2ccccc2N1)C3CCCCC3', indexes: 0, score: 1},
    {smiles: 'O=C1CN=C(c2cc(I)ccc2N1)C3CCCCC3', indexes: 30, score: 0.6904761791229248}],
  'Dice/Morgan': [
    {smiles: 'O=C1CN=C(c2ccccc2N1)C3CCCCC3', indexes: 0, score: 1},
    {smiles: 'O=C1CN=C(c2cc(I)ccc2N1)C3CCCCC3', indexes: 30, score: 0.8169013857841492}],
  'Cosine/Morgan': [
    {smiles: 'O=C1CN=C(c2ccccc2N1)C3CCCCC3', indexes: 0, score: 1},
    {smiles: 'O=C1CN=C(c2cc(I)ccc2N1)C3CCCCC3', indexes: 30, score: 0.8176316022872925}],
  'Euclidean/Morgan': [
    {smiles: 'O=C1CN=C(c2ccccc2N1)C3CCCCC3', indexes: 0, score: 1},
    {smiles: 'O=C1CN=C(c2cc(I)ccc2N1)C3CCCCC3', indexes: 30, score: 0.21712927520275116}],
  'Hamming/Morgan': [
    {smiles: 'O=C1CN=C(c2ccccc2N1)C3CCCCC3', indexes: 0, score: 1},
    {smiles: 'O=C1CN=C(c2cc(I)ccc2N1)C3CCCCC3', indexes: 30, score: 0.0714285746216774}],
};

category('top menu similarity/diversity', () => {

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

function getSearchViewer(viewer: DG.Viewer, name: string) {
  for (const v of viewer.view.viewers) {
    if (v.type === name)
      return v;
  }
}

async function _testSimilaritySearchViewerOpen() {
  const molecules = await createTableView('sar-small.csv');
  const viewer = molecules.addViewer('SimilaritySearchViewer');
  await delay(3000);
  const similaritySearchviewer = getSearchViewer(viewer, 'SimilaritySearchViewer');
  expect(similaritySearchviewer.fingerprint, 'Morgan');
  expect(similaritySearchviewer.distanceMetric, 'Tanimoto');
  expect(similaritySearchviewer.scores!.get(0), 1);
  expect(similaritySearchviewer.idxs!.get(0), 0);
  expect(similaritySearchviewer.molCol!.get(0), 'O=C1CN=C(c2ccccc2N1)C3CCCCC3');
  molecules.dataFrame.currentRowIdx = 1;
  await delay(2000);
  expect(similaritySearchviewer.curIdx, 1);
  expect(similaritySearchviewer.molCol!.get(0), 'CN1C(=O)CN=C(c2ccccc12)C3CCCCC3');
  similaritySearchviewer.close();
  molecules.close();
}

async function _testSimilaritySearchFunctionality(distanceMetric: string, fingerprint: string) {
  const molecules = await readDataframe('sar-small.csv');
  const moleculeColumn = molecules.col('smiles');
  grok.functions.call(
    `Chem:chemSimilaritySearch`, {
      'table': molecules,
      'smiles': moleculeColumn!,
      'molecule': moleculeColumn!.get(0),
      'metricName': distanceMetric,
      'limit': 10,
      'minSCore': 0.01,
      'fingerprint': fingerprint as Fingerprint,
    }).then((res) => {
    const similarityDf = res;
    //@ts-ignore
    const testResults = testSimilarityResults[`${distanceMetric}/${fingerprint}`];
    for (let i = 0; i < testResults.lenght; i++) {
      console.log(testResults);
      Object.keys(i).forEach((it) => expect(similarityDf.get(it, i), testResults[i][it]));
    }
  });
}


async function _testDiversitySearchViewerOpen() {
  const molecules = await createTableView('sar-small.csv');
  const viewer = molecules.addViewer('DiversitySearchViewer');
  await delay(3000);
  const diversitySearchviewer = getSearchViewer(viewer, 'DiversitySearchViewer');
  expect(diversitySearchviewer.fingerprint, 'Morgan');
  expect(diversitySearchviewer.distanceMetric, 'Tanimoto');
  expect(diversitySearchviewer.initialized, true);
  expect(diversitySearchviewer.renderMolIds.length > 0, true);
  diversitySearchviewer.close();
  molecules.close();
}
