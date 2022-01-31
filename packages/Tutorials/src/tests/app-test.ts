import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import { category, test, assure } from '@datagrok-libraries/utils/src/test';


category('Tutorials', () => {
  test('Launch app', async () => {
    await grok.functions.call('Tutorials:trackOverview');
    assure.notNull(document.querySelector('div.panel-content > div.tutorials-root'));
  });
});
