import * as grok from 'datagrok-api/grok';
import {category, test, expect} from '@datagrok-libraries/utils/src/test';
import {ensureContainerRunning} from '@datagrok-libraries/utils/src/test-container-utils';
import {ReactionData, Tree} from '../aizynth-api';

export const CONTAINER_TIMEOUT = 900000;

category('retrosynthesis', async () => {
  const molStr = 'CC(C(=O)OCCCc1cccnc1)c2cccc(c2)C(=O)c3ccccc3';

  test('retrosynthesis', async () => {
    await ensureContainerRunning('retrosynthesis-aizynthfinder', CONTAINER_TIMEOUT);
    const result = await grok.functions.call('Retrosynthesis:run_aizynthfind',
      {molecule: molStr, config: ''});
    const reactionData: ReactionData = JSON.parse(result);
    const paths: Tree[] = reactionData?.data?.[0]?.trees;
    expect(paths.length > 0, true);
    expect(paths[0].smiles, 'CC(C(=O)OCCCc1cccnc1)c1cccc(C(=O)c2ccccc2)c1');
    expect(paths[0].scores['state score'], 0.9940398539);
  }, {timeout: 60000 + CONTAINER_TIMEOUT});
});
