import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {before, after, category, test, expect, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {ViewHandler} from '../view-handler';

category('App', () => {
  const tabs = ['Overview', 'Packages', 'Functions', 'Events', 'Log'];
  const num = [4, 2, 2, 4, 3];

  before(async () => {
    ViewHandler.getInstance();
    if (!grok.shell.view(ViewHandler.UAname))
      await ViewHandler.getInstance().init();
    grok.shell.windows.showContextPanel = false;
  });

  test('open', async () => {
    expect(grok.shell.v.name == 'Usage Analysis', true, 'cannot find Usage Analysis view');
  });

  for (let i = 0; i < 5; i++) {
    test(tabs[i], async () => {
      ViewHandler.changeTab(tabs[i]);
      await awaitCheck(() => document.querySelectorAll('canvas').length === num[i],
        `expected ${num[i]}, got ${document.querySelectorAll('canvas').length}`, 60000);
    });
  }

  after(async () => {
    grok.shell.closeAll();
    DG.Balloon.closeAll();
    grok.shell.windows.showContextPanel = true;
  });
}, {clear: false, timeout: 60000});
