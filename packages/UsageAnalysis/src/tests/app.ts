import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {before, after, category, test, expect, awaitCheck} from '@datagrok-libraries/test/src/test';
import {ViewHandler} from '../view-handler';
import {queries} from '../package-api';


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

  const tabs = ['Overview', 'Packages', 'Functions', 'Events', 'Clicks', 'Log', 'System Activity', 'Capture', 'Timeline',
    'Projects'];

  for (const tab of tabs) {
    test(tab, async () => {
      handler.changeTab(tab);
      const view = handler.getCurrentView();
      expect(view.name === tab, true, `expected tab "${tab}", got "${view.name}"`);
      await awaitCheck(() => view.root.children.length > 0, `"${tab}" failed to initialize`, 30000);
    });
  }
}, {clear: false, timeout: 60000});

category('Capture', () => {
  test('CaptureRules columns', async () => {
    const t = await queries.captureRules('this year');
    for (const name of ['rule', 'author', 'subject', 'scope', 'reason', 'active', 'events', 'status', 'id'])
      expect(t.col(name) != null, true, `column "${name}" is missing`);
  });

  test('Clicks has request ids', async () => {
    const allUsers = (await grok.dapi.groups.getGroupsLookup('All users'))[0].id;
    const t = await queries.clicks('this week', [allUsers]);
    for (const name of ['event_time', 'user', 'event_type', 'description', 'request_id'])
      expect(t.col(name) != null, true, `column "${name}" is missing`);
  });

  test('Timeline of an unknown action is empty', async () => {
    const t: DG.DataFrame = await grok.functions.call('Timeline', {spec: JSON.stringify({action: '01M3NKD0PV60BG2W3EDPSC6ZZZ'})});
    expect(t.rowCount, 0);
    expect(t.col('requestId') != null, true, 'column "requestId" is missing');
  });
}, {timeout: 60000});
