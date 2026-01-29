import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, test, assure, delay, before} from '@datagrok-libraries/test/src/test';
// import {_package} from '../package-test';

category('Tutorials App', () => {
  test('Launch app', async () => {
    await grok.functions.call('Tutorials:trackOverview');
    assure.notNull(document.querySelector('div.panel-content > div.tutorials-root'));
  });
});
