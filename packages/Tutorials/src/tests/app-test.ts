import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, test, assure, delay, before} from '@datagrok-libraries/utils/src/test';
// import {_package} from '../package-test';

category('Tutorials App', () => {
  test('Launch app', async () => {
    await grok.functions.call('Tutorials:trackOverview');
    assure.notNull(document.querySelector('div.panel-content > div.tutorials-root'));
  });
});

category('Demo', () => {
  before(async () => {
    grok.shell.lastError = '';
  });

  const demos = DG.Func.find({package: 'Tutorials', meta: {'demoPath': null}});
  for (const demo of demos) {
    test(demo.friendlyName, async () => {
      await demo.apply();
      await delay(2000);
      if (grok.shell.lastError) {
        const err = grok.shell.lastError;
        grok.shell.lastError = '';
        throw new Error(err);
      }
    });
  }
});
