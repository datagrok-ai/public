import {after, before, category, test, expect, delay} from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';
import {createTableView, readDataframe} from './utils';
import * as grok from 'datagrok-api/grok';
import {SequenceSimilarityViewer} from '../analysis/sequence-similarity-viewer';

category('similarity/diversity', async () => {
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
  const similaritySearchviewer = getSearchViewer(viewer, 'SequenceSimilaritySearchViewer');
  if (!similaritySearchviewer.molCol)
    await waitForCompute(similaritySearchviewer);
  expect(similaritySearchviewer.fingerprint, 'Morgan');
  expect(similaritySearchviewer.distanceMetric, 'Tanimoto');
  expect(similaritySearchviewer.scores!.get(0), DG.FLOAT_NULL);
  expect(similaritySearchviewer.idxs!.get(0), 0);
  expect(similaritySearchviewer.molCol!.get(0),
    'D-Tyr_Et/Tyr_ab-dehydroMe/dV/E/N/D-Orn/D-aThr//Phe_4Me');
  expect(similaritySearchviewer.scores!.get(1), 0.4722222089767456);
  expect(similaritySearchviewer.idxs!.get(1), 11);
  expect(similaritySearchviewer.molCol!.get(1),
    'meI/hHis//Aca/meM/Tyr_ab-dehydroMe/dV/E/N/D-Orn/D-aThr//Phe_4Me');
  molecules.dataFrame.currentRowIdx = 1;
  await delay(100);
  await waitForCompute(similaritySearchviewer);
  expect(similaritySearchviewer.targetMoleculeIdx, 1);
  expect(similaritySearchviewer.molCol!.get(0),
    'meI/hHis/Aca/Cys_SEt/T/dK/Thr_PO3H2/Aca/Tyr_PO3H2/D-Chg/dV/Phe_ab-dehydro/N/D-Orn/D-aThr//Phe_4Me');
  similaritySearchviewer.close();
  molecules.close();
}


async function _testDiversitySearchViewer() {
  const molecules = await createTableView('tests/sample_MSA_data.csv');
  const viewer = molecules.addViewer('SequenceDiversitySearchViewer');
  await delay(10);
  const diversitySearchviewer = getSearchViewer(viewer, 'SequenceDiversitySearchViewer');
  if (!diversitySearchviewer.renderMolIds)
    await waitForCompute(diversitySearchviewer);
  expect(diversitySearchviewer.fingerprint, 'Morgan');
  expect(diversitySearchviewer.distanceMetric, 'Tanimoto');
  expect(diversitySearchviewer.initialized, true);
  expect(diversitySearchviewer.renderMolIds.length > 0, true);
  diversitySearchviewer.close();
  molecules.close();
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
