import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {before, after, category, test, expect, awaitCheck} from '@datagrok-libraries/test/src/test';
import {ViewHandler} from '../view-handler';


category('App', () => {
  let handler: ViewHandler;

  before(async () => {
    handler = new ViewHandler();
    grok.shell.addView(handler.view);
    await handler.init();
  });

  test('open', async () => {
    expect(handler.view != null, true, 'view not created');
    expect(handler.view.tabs != null, true, 'tabs not initialized');
  });

  const tabs = ['Overview', 'Packages', 'Functions', 'Events', 'Clicks', 'Log', 'Projects'];

  for (const tab of tabs) {
    test(tab, async () => {
      handler.changeTab(tab);
      const view = handler.getCurrentView();
      expect(view.name === tab, true, `expected tab "${tab}", got "${view.name}"`);
      await awaitCheck(() => view.root.children.length > 0, `"${tab}" failed to initialize`, 30000);
    });
  }

  after(async () => {
    grok.shell.closeAll();
    DG.Balloon.closeAll();
  });
}, {clear: false, timeout: 60000});
