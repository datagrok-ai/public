import * as grok from 'datagrok-api/grok';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {ReactionData, Tree} from '../aizynth-api';
import {before, timeout} from '@datagrok-libraries/test/src/test';


category('retrosynthesis', async () => {
  const molStr = 'CC(C(=O)OCCCc1cccnc1)c2cccc(c2)C(=O)c3ccccc3';

  before(async () => {
    try {
      await timeout(
        () => grok.functions.call('Retrosynthesis:check_health', {}),
        180000,
        'Health check timed out',
      );
    } catch (err: any) {
      throw new Error(`Health check failed: ${err}`);
    }
  });

  test('retrosynthesis', async () => {
    const result = await grok.functions.call('Retrosynthesis:run_aizynthfind',
      {molecule: molStr, config: ''});
    const reactionData: ReactionData = JSON.parse(result);
    const paths: Tree[] = reactionData?.data?.[0]?.trees;
    expect(paths.length > 0, true);
    expect(paths[0].smiles, 'CC(C(=O)OCCCc1cccnc1)c1cccc(C(=O)c2ccccc2)c1');
    expect(paths[0].scores['state score'], 0.9940398539);
  }, {timeout: 60000});
});
