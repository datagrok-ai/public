import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, expect, expectFloat, expectArray, test, delay, before} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {Fingerprint} from '../utils/chem-common';
import {createTableView, readDataframe, molV2000, molV3000} from './utils';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';

import {
  // _testSearchSubstructure,
  // _testSearchSubstructureAllParameters,
  // _testSearchSubstructureSARSmall,
  loadFileAsText,
} from './utils';
import {findSimilar, getSimilarities} from '../package';
import {chemDiversitySearch} from '../analysis/chem-diversity-viewer';
import {chemSimilaritySearch} from '../analysis/chem-similarity-viewer';
import {tanimotoSimilarity} from '@datagrok-libraries/ml/src/distance-metrics-methods';

const testSimilarityResults = {
  'Tanimoto/Morgan': [
    {smiles: 'O=C1CN=C(c2ccccc2N1)C3CCCCC3', indexes: 0, score: 1},
    {smiles: 'CN(C)c1ccc2NC(=O)CN=C(c2c1)C3CCCCC3', indexes: 14, score: 0.644444465637207}],
  'Dice/Morgan': [
    {smiles: 'O=C1CN=C(c2ccccc2N1)C3CCCCC3', indexes: 0, score: 1},
    {smiles: 'CN(C)c1ccc2NC(=O)CN=C(c2c1)C3CCCCC3', indexes: 14, score: 0.7837837934494019}],
  'Cosine/Morgan': [
    {smiles: 'O=C1CN=C(c2ccccc2N1)C3CCCCC3', indexes: 0, score: 1},
    {smiles: 'CN(C)c1ccc2NC(=O)CN=C(c2c1)C3CCCCC3', indexes: 14, score: 0.7863729000091553}],
  'Euclidean/Morgan': [
    {smiles: 'O=C1CN=C(c2ccccc2N1)C3CCCCC3', indexes: 0, score: 1},
    {smiles: 'CN(C)c1ccc2NC(=O)CN=C(c2c1)C3CCCCC3', indexes: 14, score: 0.20000000298023224}],
  'Hamming/Morgan': [
    {smiles: 'O=C1CN=C(c2ccccc2N1)C3CCCCC3', indexes: 0, score: 1},
    {smiles: 'CN(C)c1ccc2NC(=O)CN=C(c2c1)C3CCCCC3', indexes: 14, score: 0.05882352963089943}],
};

category('top menu similarity/diversity', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('findSimilar.chem.smiles', async () => {
    await _testFindSimilar(findSimilar);
  });

  // test('findSimilar.chem.molV2000', async () => {
  //   await _testFindSimilar(findSimilar, molV2000, 'V2000');
  // });

  test('findSimilar.chem.molV3000', async () => {
    await _testFindSimilar(findSimilar, molV3000, 'V3000');
  });

  test('getSimilarities.chem.molecules', async () => {
    await _testGetSimilarities(getSimilarities);
  });

  test('testDiversitySearch', async () => {
    const t = grok.data.demo.molecules();
    await chemDiversitySearch(t.col('smiles')!, tanimotoSimilarity, 10, 'Morgan' as Fingerprint);
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

export async function _testFindSimilar(findSimilarFunction: (...args: any) => Promise<DG.DataFrame | null>,
  molecule: string = 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1', notation?: string) {
  let dfInput: DG.DataFrame;
  if (notation)
    dfInput = DG.DataFrame.fromCsv(await loadFileAsText(notation === 'V2000' ? 'tests/spgi-100.csv' : 'tests/approved-drugs-100.csv'));
  else
    dfInput = DG.DataFrame.fromCsv(await loadFileAsText('sar-small.csv'));
  await grok.data.detectSemanticTypes(dfInput);
  const colInput = dfInput.columns.bySemType(DG.SEMTYPE.MOLECULE)!;
  const numRowsOriginal = dfInput.rowCount;
  const dfResult: DG.DataFrame = // shouldn't be null
    (await findSimilarFunction(colInput, molecule))!;
  const numRows = dfResult.rowCount;
  const columnNames = [
    dfResult.columns.byIndex(0).name,
    dfResult.columns.byIndex(1).name,
    dfResult.columns.byIndex(2).name,
  ];
  const arr: any[] = [];
  for (let i = 0; i < 5; ++i) {
    const molecule: string = dfResult.columns.byIndex(0).get(i);
    const score: number = dfResult.columns.byIndex(1).get(i);
    const index: number = dfResult.columns.byIndex(2).get(i);
    arr[i] = {molecule, score, index};
  }
  expectArray([numRows, columnNames[0], columnNames[1], columnNames[2]],
    [numRowsOriginal, 'molecule', 'score', 'index']);
  switch (notation) {
  case 'V2000':
    break;
  case 'V3000':
    expectArray([arr[0].molecule, arr[1].molecule, arr[2].molecule, arr[3].molecule, arr[4].molecule],
      [colInput.get(0), colInput.get(71), colInput.get(43), colInput.get(87), colInput.get(89)]);
    expectArray([arr[0].score, arr[1].score, arr[2].score, arr[3].score, arr[4].score],
      [1, 0.13725490868091583, 0.12903225421905518, 0.12820513546466827, 0.12820513546466827]);
    expectArray([arr[0].index, arr[1].index, arr[2].index, arr[3].index, arr[4].index],
      [0, 71, 43, 87, 89]);
    break;
  default:
    expectArray([arr[0].molecule, arr[1].molecule, arr[2].molecule, arr[3].molecule, arr[4].molecule],
      ['O=C1CN=C(c2ccccc2N1)C3CCCCC3', 'O=C1CN=C(c2cc(I)ccc2N1)C3CCCCC3', 'O=C1CN=C(c2cc(Cl)ccc2N1)C3CCCCC3',
        'O=C1CN=C(c2cc(F)ccc2N1)C3CCCCC3', 'O=C1CN=C(c2cc(Br)ccc2N1)C3CCCCC3']);
    expectArray([arr[0].score, arr[1].score, arr[2].score, arr[3].score, arr[4].score],
      [1, 0.6904761791229248, 0.6744186282157898, 0.6744186282157898, 0.6744186282157898]);
    expectArray([arr[0].index, arr[1].index, arr[2].index, arr[3].index, arr[4].index],
      [0, 30, 5, 20, 25]);  
    break;
  }
}

export async function _testGetSimilarities(getSimilaritiesFunction: (...args: any) => Promise<any>) {
  const df = grok.data.demo.molecules();
  let scores = (await getSimilaritiesFunction(df.columns.byName('smiles'),
    'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1'))! as any;
  if (scores instanceof DG.DataFrame)
    scores = scores.columns.byIndex(0);

  expectFloat(scores.get(0), 0.1034);
  expectFloat(scores.get(1), 0.07407);
  expectFloat(scores.get(2), 0.11111);
  expectFloat(scores.get(3), 0.11111);
  expectFloat(scores.get(4), 0.07042);
  expectFloat(scores.get(5), 0.06349);
}

async function _testSimilaritySearchViewerOpen() {
  const molecules = await createTableView('tests/sar-small_test.csv');
  const viewer = molecules.addViewer('Chem Similarity Search');
  await delay(1000);
  const similaritySearchviewer = getSearchViewer(viewer, 'Chem Similarity Search');
  expect(similaritySearchviewer.fingerprint, 'Morgan');
  expect(similaritySearchviewer.distanceMetric, 'Tanimoto');
  expect(similaritySearchviewer.scores!.get(0), 1);
  expect(similaritySearchviewer.idxs!.get(0), 0);
  expect(similaritySearchviewer.molCol!.get(0), 'O=C1CN=C(c2ccccc2N1)C3CCCCC3');
  molecules.dataFrame.currentRowIdx = 1;
  await delay(100);
  expect(similaritySearchviewer.curIdx, 1);
  expect(similaritySearchviewer.molCol!.get(0), 'CN1C(=O)CN=C(c2ccccc12)C3CCCCC3');
  similaritySearchviewer.close();
  molecules.close();
}

async function _testSimilaritySearchFunctionality(distanceMetric: string, fingerprint: string) {
  const molecules = await readDataframe('tests/sar-small_test.csv');
  const moleculeColumn = molecules.col('smiles');
  const similarityDf = await chemSimilaritySearch(molecules, moleculeColumn!, moleculeColumn!.get(0),
    distanceMetric, 10, 0.01, fingerprint as Fingerprint);
    //@ts-ignore
  const testResults = testSimilarityResults[`${distanceMetric}/${fingerprint}`];
  for (let i = 0; i < testResults.length; i++) {
    console.log(testResults);
    Object.keys(testResults[i]).forEach((it) => expect(similarityDf.get(it, i), testResults[i][it]));
  }
}

async function _testDiversitySearchViewerOpen() {
  const molecules = await createTableView('tests/sar-small_test.csv');
  const viewer = molecules.addViewer('Chem Diversity Search');
  await delay(500);
  const diversitySearchviewer = getSearchViewer(viewer, 'Chem Diversity Search');
  expect(diversitySearchviewer.fingerprint, 'Morgan');
  expect(diversitySearchviewer.distanceMetric, 'Tanimoto');
  expect(diversitySearchviewer.initialized, true);
  expect(diversitySearchviewer.renderMolIds.length > 0, true);
  diversitySearchviewer.close();
  molecules.close();
}
