import {after, before, category, test, expect, awaitCheck} from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';
import {createTableView} from './utils';
import * as grok from 'datagrok-api/grok';
import {SequenceSimilarityViewer} from '../analysis/sequence-similarity-viewer';

category('similarity/diversity', async () => {
  before(async () => {
    grok.shell.closeAll();
  });

  after(async () => {
    grok.shell.closeAll();
  });

  test('similaritySearchViewer', async () => {
    await _testSimilaritySearchViewer();
  });
  test('diversitySearchViewer', async () => {
    await _testDiversitySearchViewer();
  });
});

async function _testSimilaritySearchViewer() {
  try {
    const molecules = await createTableView('tests/sample_MSA_data.csv');
    const viewer = molecules.addViewer('Sequence Similarity Search');
    await awaitCheck(() => getSearchViewer(viewer, 'Sequence Similarity Search') !== undefined,
      'Sequence Similarity Search has not been created', 5000);
    const similaritySearchViewer: SequenceSimilarityViewer = getSearchViewer(viewer, 'Sequence Similarity Search');
    await awaitCheck(() => similaritySearchViewer.root.getElementsByClassName('d4-grid').length !== 0,
      'Sequence Similarity Search has not been created', 5000);
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
    molecules.dataFrame.currentRowIdx = 1;
    await awaitCheck(() => similaritySearchViewer.targetMoleculeIdx === 1, 'Target molecule has not been changed', 5000);
    await awaitCheck(() => similaritySearchViewer.molCol!.get(0) === 
    'meI/hHis/Aca/Cys_SEt/T/dK/Thr_PO3H2/Aca/Tyr_PO3H2/D-Chg/dV/Phe_ab-dehydro/N/D-Orn/D-aThr//Phe_4Me',
    'Incorrect first similar molecule', 5000);
    
  } finally {
    grok.shell.closeAll();
  }
}

async function _testDiversitySearchViewer() {
  try {
    const molecules = await createTableView('tests/sample_MSA_data.csv');
    const viewer = molecules.addViewer('Sequence Diversity Search');
    await awaitCheck(() => getSearchViewer(viewer, 'Sequence Diversity Search') !== undefined,
      'Sequence Diversity Search has not been created', 5000);
    const diversitySearchviewer = getSearchViewer(viewer, 'Sequence Diversity Search');
    await awaitCheck(() => diversitySearchviewer.root.getElementsByClassName('d4-grid').length !== 0,
      'Sequence Diversity Search has not been created', 5000);
    expect(diversitySearchviewer.fingerprint, 'Morgan');
    expect(diversitySearchviewer.distanceMetric, 'Tanimoto');
    expect(diversitySearchviewer.initialized, true);
    expect(diversitySearchviewer.renderMolIds.length > 0, true);
  } finally {
    grok.shell.closeAll();
  }
}

function getSearchViewer(viewer: DG.Viewer, name: string) {
  //@ts-ignore
  for (const v of viewer.view.viewers) {
    if (v.type === name)
      return v;
  }
}

