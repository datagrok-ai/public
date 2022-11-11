import {after, before, category, test, expect, delay} from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';
import {createTableView, readDataframe} from './utils';
import * as grok from 'datagrok-api/grok';
import {SequenceSimilarityViewer} from '../analysis/sequence-similarity-viewer';

let viewList: DG.ViewBase[];
let dfList: DG.DataFrame[];


category('similarity/diversity', async () => {

  before(async () => {
    viewList = [];
    dfList = [];
  });

  after(async () => {
    for (const view of viewList) view.close();
    for (const df of dfList) grok.shell.closeTable(df);
  });


  test('similaritySearchViewer', async () => {
    await _testSimilaritySearchViewer();
  });
  test('diversitySearchViewer', async () => {
    await _testDiversitySearchViewer();
  });
});

async function _testSimilaritySearchViewer() {
  const molecules = await createTableView('tests/sample_MSA_data.csv');
  const viewer = molecules.addViewer('SequenceSimilaritySearchViewer');
  await delay(100);
  const similaritySearchViewer = getSearchViewer(viewer, 'SequenceSimilaritySearchViewer');
  viewList.push(similaritySearchViewer);
  viewList.push(molecules);
  if (!similaritySearchViewer.molCol)
    await waitForCompute(similaritySearchViewer);
  expect(similaritySearchViewer.fingerprint, 'Morgan');
  expect(similaritySearchViewer.distanceMetric, 'Tanimoto');
  expect(similaritySearchViewer.scores!.get(0), DG.FLOAT_NULL);
  expect(similaritySearchViewer.idxs!.get(0), 0);
  expect(similaritySearchViewer.molCol!.get(0),
    'D-Tyr_Et/Tyr_ab-dehydroMe/dV/E/N/D-Orn/D-aThr//Phe_4Me');
  expect(similaritySearchViewer.scores!.get(1), 0.4722222089767456);
  expect(similaritySearchViewer.idxs!.get(1), 11);
  expect(similaritySearchViewer.molCol!.get(1),
    'meI/hHis//Aca/meM/Tyr_ab-dehydroMe/dV/E/N/D-Orn/D-aThr//Phe_4Me');
  const waiter = waitForCompute(similaritySearchViewer); /* subscribe for computeCompleted event before start compute */
  molecules.dataFrame.currentRowIdx = 1;
  await waiter;
  expect(similaritySearchViewer.targetMoleculeIdx, 1);
  expect(similaritySearchViewer.molCol!.get(0),
    'meI/hHis/Aca/Cys_SEt/T/dK/Thr_PO3H2/Aca/Tyr_PO3H2/D-Chg/dV/Phe_ab-dehydro/N/D-Orn/D-aThr//Phe_4Me');
}


async function _testDiversitySearchViewer() {
  const molecules = await createTableView('tests/sample_MSA_data.csv');
  const viewer = molecules.addViewer('SequenceDiversitySearchViewer');
  await delay(10);
  const diversitySearchviewer = getSearchViewer(viewer, 'SequenceDiversitySearchViewer');
  viewList.push(diversitySearchviewer);
  viewList.push(molecules);
  if (!diversitySearchviewer.renderMolIds)
    await waitForCompute(diversitySearchviewer);
  expect(diversitySearchviewer.fingerprint, 'Morgan');
  expect(diversitySearchviewer.distanceMetric, 'Tanimoto');
  expect(diversitySearchviewer.initialized, true);
  expect(diversitySearchviewer.renderMolIds.length > 0, true);
}

function getSearchViewer(viewer: DG.Viewer, name: string) {
  for (const v of viewer.view.viewers) {
    if (v.type === name)
      return v;
  }
}

async function waitForCompute(viewer: SequenceSimilarityViewer) {
  const t = new Promise((resolve, reject) => {
    viewer.computeCompleted.subscribe(async (_: any) => {
      try {
        resolve(true);
      } catch (error) {
        reject(error);
      }
    });
  });
  await t;
}
