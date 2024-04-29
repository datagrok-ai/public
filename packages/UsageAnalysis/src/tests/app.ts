import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {before, after, category, test, expect, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {ViewHandler} from '../view-handler';

category('App', () => {
  const tabs = ['Overview', 'Packages', 'Functions', 'Events', 'Log', 'Tests', 'Errors', 'Reports'];
  const num = [4, 4, 5, 4, 3, 5, 3, 3];
  let initTime: number = 0;

  before(async () => {
    if (grok.shell.sidebar.panes.every((p) => p.name != ViewHandler.UAname)) {
      await ViewHandler.getInstance().init();
      grok.shell.addView(ViewHandler.UA);
      grok.shell.sidebar.currentPane = grok.shell.sidebar.getPane(ViewHandler.UAname);
    }
  });

  test('open', async () => {
    expect(grok.shell.v.name === 'Usage Analysis', true, 'cannot find Usage Analysis view');
  });

  for (let i = 0; i < 8; i++) {
    test(tabs[i], async () => {
      ViewHandler.changeTab(tabs[i]);
      const view = grok.shell.v.root;
      let err = null;
      const s = ['Log', 'Errors', 'Reports'].includes(tabs[i]) ? '.grok-wait + .d4-grid canvas' : 'canvas';
      try {
        await awaitCheck(() => view.querySelectorAll(s).length === num[i], '', 10000);
      } catch (e) {
        err = 'Tab failed to load in 10 seconds';
      }
      await awaitCheck(() => view.querySelectorAll(s).length === num[i],
        `expected ${num[i]}, got ${document.querySelectorAll(s).length}`, 45000);
      if (err)
        throw new Error(err);
    });
  }

  after(async () => {
    grok.shell.closeAll();
    DG.Balloon.closeAll();
    grok.shell.windows.showContextPanel = true;
  });
}, {clear: false, timeout: 60000});
