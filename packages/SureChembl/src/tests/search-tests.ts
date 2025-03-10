import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, test, expect} from '@datagrok-libraries/utils/src/test';
import {ensureContainerRunning} from '@datagrok-libraries/utils/src/test-container-utils';

export const CONTAINER_TIMEOUT = 900000;

category('Search tests', () => {
  test('Substructure search', async () => {
    await ensureContainerRunning('surechembl', CONTAINER_TIMEOUT);
    const df: DG.DataFrame | null = await grok.functions.call('Surechembl:sureChemblSubstructureSearch', {
      molecule: 'FC(F)(F)c1ccc(OC2CCNCC2)cc1',
      limit: 10,
    });
    expect(df?.rowCount, 17);
    const patentIdx = df!.col('doc_id')?.toList().findIndex((it) => it === 'US-11918575-B2');
    expect(patentIdx !== undefined && patentIdx !== -1, true);
    expect(df!.get('smiles', patentIdx!), 'FC(F)(F)c1ccc(OC2CCNCC2)cc1');
    expect(df!.get('doc_surechembl_id', patentIdx!), 81759);

    if (df != null)
      grok.shell.closeTable(df);
  }, {timeout: 30000 + CONTAINER_TIMEOUT});

  test('Similarity search', async () => {
    await ensureContainerRunning('surechembl', CONTAINER_TIMEOUT);
    const df: DG.DataFrame | null = await grok.functions.call('Surechembl:sureChemblSimilaritySearch', {
      molecule: 'FC(F)(F)c1ccc(OC2CCNCC2)cc1',
      limit: 10,
      similarityThreshold: 0.6,
    });
    expect(df?.rowCount, 19);
    const patentIdx = df!.col('doc_id')?.toList().findIndex((it) => it === 'US-11918575-B2');
    expect(patentIdx !== undefined && patentIdx !== -1, true);
    expect(df!.get('smiles', patentIdx!), 'FC(F)(F)c1ccc(OC2CCNCC2)cc1');
    expect(df!.get('doc_surechembl_id', patentIdx!), 81759);

    if (df != null)
      grok.shell.closeTable(df);
  }, {timeout: 30000 + CONTAINER_TIMEOUT});
});
