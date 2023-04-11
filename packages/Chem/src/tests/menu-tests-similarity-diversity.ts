import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, expect, expectFloat, expectArray, test, delay, before, after, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {Fingerprint} from '../utils/chem-common';
import {createTableView, readDataframe, molV2000, molV3000} from './utils';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {findSimilar, getSimilarities} from '../package';
import {chemDiversitySearch, ChemDiversityViewer} from '../analysis/chem-diversity-viewer';
import {chemSimilaritySearch, ChemSimilarityViewer} from '../analysis/chem-similarity-viewer';
import {tanimotoSimilarity} from '@datagrok-libraries/ml/src/distance-metrics-methods';


category('top menu similarity/diversity', () => {
  let spgi100: DG.DataFrame;
  let approvedDrugs100: DG.DataFrame;
  let molecules: DG.DataFrame;
  let empty: DG.DataFrame;

  before(async () => {
    grok.shell.closeAll();
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    spgi100 = await readDataframe('tests/spgi-100.csv');
    approvedDrugs100 = await readDataframe('tests/approved-drugs-100.csv');
    molecules = grok.data.demo.molecules(100);
    empty = await readDataframe('tests/sar-small_empty_vals.csv');
    await grok.data.detectSemanticTypes(empty);
  });

  test('findSimilar.chem.smiles', async () => {
    await _testFindSimilar(findSimilar);
  });

  test('findSimilar.chem.molV2000', async () => {
    await _testFindSimilar(findSimilar, molV2000, 'V2000');
  });

  test('findSimilar.chem.molV3000', async () => {
    await _testFindSimilar(findSimilar, molV3000, 'V3000');
  });

  test('getSimilarities.chem.molecules', async () => {
    await _testGetSimilarities(getSimilarities, molecules);
  });

  test('similarity.emptyValues', async () => {
    empty.currentRowIdx = 1;
    const tv = grok.shell.addTableView(empty);
    tv.addViewer('Chem Similarity Search');
    await delay(500);
    const viewer = getSearchViewer(tv, 'Chem Similarity Search') as ChemSimilarityViewer;
    try {
      expectArray(viewer.scores!.toList(), [1, 0.6808510422706604, 0.6739130616188049, 0.6521739363670349, 0.6458333134651184,
        0.4821428656578064, 0.4736842215061188, 0.4736842215061188, 0.4655172526836395, 0.4576271176338196]);
    } finally {tv.close();}
  }, {skipReason: 'GROK-12227'});

  test('similarity.emptyInput', async () => {
    empty.currentRowIdx = 0;
    const tv = grok.shell.addTableView(empty);
    DG.Balloon.closeAll();
    tv.addViewer('Chem Similarity Search');
    try {
      await awaitCheck(() => document.querySelector('.d4-balloon-content')?.innerHTML ===
        'Empty molecule cannot be used for similarity search', 'cannot find error balloon', 2000);
    } finally {tv.close();}
  }, {skipReason: 'GROK-12227'});

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

  test('testDiversitySearch.molecules', async () => {
    await chemDiversitySearch(molecules.getCol('smiles'), tanimotoSimilarity, 10, 'Morgan' as Fingerprint);
  });

  test('testDiversitySearch.molV2000', async () => {
    await chemDiversitySearch(spgi100.getCol('Structure'), tanimotoSimilarity, 10, 'Morgan' as Fingerprint);
  });

  test('testDiversitySearch.molV3000', async () => {
    await chemDiversitySearch(approvedDrugs100.getCol('molecule'), tanimotoSimilarity, 10, 'Morgan' as Fingerprint);
  });

  test('diversity.emptyValues', async () => {
    const tv = grok.shell.addTableView(empty);
    tv.addViewer('Chem Diversity Search');
    await delay(500);
    const viewer = getSearchViewer(tv, 'Chem Diversity Search') as ChemDiversityViewer;
    try {
      expect(viewer.renderMolIds.length, 10);
    } finally {tv.close();}
  });

  test('diversitySearchViewerOpen', async () => {
    await _testDiversitySearchViewerOpen();
  });

  after(async () => {
    grok.shell.closeAll();
  });
});

function getSearchViewer(tv: DG.TableView, name: string): DG.Viewer {
  for (const v of tv.viewers) {
    if (v.type === name)
      return v;
  }
  throw 'Search viewer not found';
}

export async function _testFindSimilar(findSimilarFunction: (...args: any) => Promise<DG.DataFrame | null>,
  molecule: string = 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1', notation?: string) {
  let dfInput: DG.DataFrame;
  if (notation) dfInput = await readDataframe(notation === 'V2000' ? 'tests/spgi-100.csv' : 'tests/approved-drugs-100.csv');
  else dfInput = await readDataframe('sar-small.csv');
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
    expectArray([arr[0].molecule, arr[1].molecule, arr[2].molecule, arr[3].molecule, arr[4].molecule],
      [colInput.get(82), colInput.get(58), colInput.get(51), colInput.get(59), colInput.get(93)]);
    expectArray([arr[0].score, arr[1].score, arr[2].score, arr[3].score, arr[4].score],
      [0.16867469251155853, 0.14814814925193787, 0.14444445073604584, 0.14102564752101898, 0.13953489065170288]);
    expectArray([arr[0].index, arr[1].index, arr[2].index, arr[3].index, arr[4].index],
      [82, 58, 51, 59, 93]);
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

export async function _testGetSimilarities(getSimilaritiesFunction: (...args: any) => Promise<any>, df: DG.DataFrame) {
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
  molecules.addViewer('Chem Similarity Search');
  await delay(500);
  const similaritySearchviewer = getSearchViewer(molecules, 'Chem Similarity Search') as ChemSimilarityViewer;
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
  molecules.addViewer('Chem Diversity Search');
  await delay(500);
  const diversitySearchviewer = getSearchViewer(molecules, 'Chem Diversity Search') as ChemDiversityViewer;
  expect(diversitySearchviewer.fingerprint, 'Morgan');
  expect(diversitySearchviewer.distanceMetric, 'Tanimoto');
  expect(diversitySearchviewer.initialized, true);
  expect(diversitySearchviewer.renderMolIds.length > 0, true);
  diversitySearchviewer.close();
  molecules.close();
}

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
