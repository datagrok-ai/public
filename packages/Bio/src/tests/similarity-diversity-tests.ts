import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, test, expect, awaitCheck, delay} from '@datagrok-libraries/utils/src/test';
import {SequenceSimilarityViewer} from '../analysis/sequence-similarity-viewer';
import {SequenceDiversityViewer} from '../analysis/sequence-diversity-viewer';

import {_package} from '../package-test';

category('similarity/diversity', async () => {
  test('similaritySearchViewer', async () => {
    await _testSimilaritySearchViewer();
  });

  test('diversitySearchViewer', async () => {
    await _testDiversitySearchViewer();
  });
});

async function _testSimilaritySearchViewer() {
  const csv = await _package.files.readAsText('tests/sample_MSA_data.csv');
  const df = DG.DataFrame.fromCsv(csv);
  await grok.data.detectSemanticTypes(df);
  const moleculesView = grok.shell.addTableView(df);
  const seqCol = moleculesView.dataFrame.getCol('MSA');
  expect(seqCol.semType, DG.SEMTYPE.MACROMOLECULE);

  const viewer: SequenceSimilarityViewer = (await moleculesView.dataFrame.plot
    .fromType('Sequence Similarity Search')) as SequenceSimilarityViewer;
  let computeCompleted: boolean = false;
  viewer.computeCompleted.subscribe((value: any) => {
    if (value) computeCompleted = true;
  });
  moleculesView.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Similarity');
  await viewer.renderPromise; // required to wait for computeCompleted

  await awaitCheck(() => getSearchViewer(moleculesView, 'Sequence Similarity Search') !== undefined,
    'Sequence Similarity Search viewer has not been created', 100);
  if (!viewer.initialized) throw new Error('The viewer is not initialized.');
  if (!viewer.targetColumn) throw new Error('The viewer has not molecule column (onTableAttached).');
  if (!viewer.beforeRender()) throw new Error('The viewer is not able to render.');
  if (!viewer.computeRequested) throw new Error('The viewer has not compute requested even.');
  if (!computeCompleted) throw new Error('The viewer has not compute completed.');

  const similaritySearchViewer = viewer;
  await awaitCheck(() => similaritySearchViewer.root.getElementsByClassName('d4-grid').length !== 0,
    'Sequence Similarity Search viewer grid has not been created', 100);

  /* eslint-disable max-len */
  const str0: string = 'D-Tyr_Et/Tyr_ab-dehydroMe/dV/E/N/D-Orn/D-aThr//Phe_4Me';
  const str1: string = 'Thr_PO3H2/Aca/D-Tyr_Et/Tyr_ab-dehydroMe/dV/Chg/N/D-Orn/D-aThr//Phe_4Me';
  const simStr1: string = 'meI/hHis/Aca/Cys_SEt/T/dK/Thr_PO3H2/Aca/Tyr_PO3H2/D-Chg/dV/Phe_ab-dehydro/N/D-Orn/D-aThr//Phe_4Me';
  /* eslint-enable max-len */

  expect(similaritySearchViewer.fingerprint, 'Morgan');
  expect(similaritySearchViewer.distanceMetric, 'Tanimoto');
  expect(similaritySearchViewer.scores!.get(0), DG.FLOAT_NULL);
  expect(similaritySearchViewer.idxs!.get(0), 0);
  expect(similaritySearchViewer.molCol!.get(0), str0);
  expect((similaritySearchViewer.scores!.get(1) as number).toFixed(2), '0.73');
  expect(similaritySearchViewer.idxs!.get(1), 4);
  expect(similaritySearchViewer.molCol!.get(1), str1);
  moleculesView.dataFrame.currentRowIdx = 1;
  await awaitCheck(() => similaritySearchViewer.targetMoleculeIdx === 1,
    'Target molecule has not been changed', 5000);
  await awaitCheck(() => similaritySearchViewer.molCol!.get(0) === simStr1,
    'Incorrect first similar molecule', 5000);
}

async function _testDiversitySearchViewer() {
  const csv = await _package.files.readAsText('tests/sample_MSA_data.csv');
  const df = DG.DataFrame.fromCsv(csv);
  const moleculesView = grok.shell.addTableView(df);
  await grok.data.detectSemanticTypes(df);
  const seqCol = moleculesView.dataFrame.getCol('MSA');
  expect(seqCol.semType, DG.SEMTYPE.MACROMOLECULE);

  const viewer: SequenceDiversityViewer = (await moleculesView.dataFrame.plot
    .fromType('Sequence Diversity Search')) as SequenceDiversityViewer;
  let computeCompleted: boolean = false;
  viewer.computeCompleted.subscribe((value: boolean) => {
    if (value) computeCompleted = true;
  });
  moleculesView.dockManager.dock(viewer, DG.DOCK_TYPE.DOWN, null, 'Diversity');
  await viewer.renderPromise;

  await awaitCheck(() => getSearchViewer(moleculesView, 'Sequence Diversity Search') !== undefined,
    'Sequence Diversity Search viewer has not been created', 100);
  if (!viewer.initialized) throw new Error('The viewer is not initialized.');
  if (!viewer.targetColumn) throw new Error('The viewer has not molecule column (onTableAttached).');
  if (!viewer.beforeRender()) throw new Error('The viewer is not able to render.');
  if (!viewer.computeRequested) throw new Error('The viewer has not compute requested even.');
  if (!computeCompleted) throw new Error('The viewer has not compute completed.');

  //const diversitySearchViewer = getSearchViewer(viewer, 'Sequence Diversity Search');
  const diversitySearchViewer = viewer;
  await awaitCheck(() => diversitySearchViewer.root.getElementsByClassName('d4-grid').length !== 0,
    'Sequence Diversity Search viewer grid has not been created', 100);

  expect(diversitySearchViewer.fingerprint, 'Morgan');
  expect(diversitySearchViewer.distanceMetric, 'Tanimoto');
  expect(diversitySearchViewer.initialized, true);
  expect(diversitySearchViewer.renderMolIds!.length > 0, true);
}

function getSearchViewer(tv: DG.TableView, name: string): DG.Viewer | null {
  let res: DG.Viewer | null = null;
  for (const v of tv.viewers) {
    if (v.type === name)
      res = v;
  }
  return res;
}
