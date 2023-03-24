import {after, before, category, test, expect, delay} from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';
import {createTableView, readDataframe} from './utils';
import * as grok from 'datagrok-api/grok';
import {SequenceSimilarityViewer} from '../analysis/sequence-similarity-viewer';
import {SequenceDiversityViewer} from '../analysis/sequence-diversity-viewer';
import {SequenceSearchBaseViewer} from '../analysis/sequence-search-base-viewer';

let viewList: DG.ViewBase[];
let viewerList: DG.Viewer[];
let dfList: DG.DataFrame[];


category('similarity/diversity', async () => {
  before(async () => {
    viewList = [];
    viewerList = [];
    dfList = [];
  });

  after(async () => {
    for (const viewer of viewerList) viewer.close();
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
  const viewer = molecules.addViewer('Sequence Similarity Search');
  await delay(100);
  const similaritySearchViewer: SequenceSimilarityViewer =
    getSearchViewer(viewer, 'Sequence Similarity Search')! as SequenceSimilarityViewer;
  viewerList.push(similaritySearchViewer);
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
  const viewer = molecules.addViewer('Sequence Diversity Search');
  await delay(10);
  const diversitySearchViewer: SequenceDiversityViewer =
    getSearchViewer(viewer, 'Sequence Diversity Search')! as SequenceDiversityViewer;
  viewerList.push(diversitySearchViewer);
  viewList.push(molecules);
  if (!diversitySearchViewer.renderMolIds)
    await waitForCompute(diversitySearchViewer);
  expect(diversitySearchViewer.fingerprint, 'Morgan');
  expect(diversitySearchViewer.distanceMetric, 'Tanimoto');
  expect(diversitySearchViewer.initialized, true);
  expect(!!diversitySearchViewer.renderMolIds, true);
  expect(diversitySearchViewer.renderMolIds!.length > 0, true);
}

function getSearchViewer(viewer: DG.Viewer, name: string): DG.Viewer | null {
  for (const v of (viewer.view as DG.TableView).viewers) {
    if (v.type === name)
      return v;
  }
  return null;
}

async function waitForCompute(viewer: SequenceSearchBaseViewer) {
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
