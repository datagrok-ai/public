import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {delay, expect} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {Fingerprint} from '../utils/chem-common';
import { createTableView, readDataframe } from './utils';

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


function getSearchViewer(viewer: DG.Viewer, name: string) {
  for (const v of viewer.view.viewers) {
    if (v.type === name)
      return v;
  }
}

export async function _testSimilaritySearchViewerOpen() {
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

export async function _testSimilaritySearchFunctionality(distanceMetric: string, fingerprint: string) {
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


export async function _testDiversitySearchViewerOpen() {
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
