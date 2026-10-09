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

  const tabs = ['Overview', 'Packages', 'Functions', 'Events', 'Clicks', 'Log', 'System Activity', 'Errors', 'Projects'];

  for (const tab of tabs) {
    test(tab, async () => {
      handler.changeTab(tab);
      const view = handler.getCurrentView();
      expect(view.name === tab, true, `expected tab "${tab}", got "${view.name}"`);
      await awaitCheck(() => view.root.children.length > 0, `"${tab}" failed to initialize`, 30000);
    });
  }

  test('Errors: top errors and their incident', async () => {
    const errors = (await grok.dapi.admin.getMetrics({date: 'this month', limit: 5, errorsBy: 'signature'})).errors;
    expect(errors.now.count >= errors.top.length, true);
    for (const e of errors.top)
      expect(typeof e.signature, 'string');
    const signature = errors.top[0]?.signature ?? '00000000-0000-0000-0000-000000000000';
    expect((await queries.errorAlerts(signature)).col('alert') != null, true);
  });

  test('Log: events carry request ids for the timeline', async () => {
    const allUsers = (await grok.dapi.groups.getGroupsLookup('All users'))[0].id;
    const t = await queries.logTail('today', [allUsers], ['all']);
    expect(t.col('request_id') != null, true, 'LogTail has no request_id');
    const requestId = t.col('request_id')?.toList().find((id) => id);
    if (requestId)
      expect((await grok.dapi.log.getTimeline({action: requestId.split('.')[0], limit: 5})).length > 0, true);
  });
}, {clear: false, timeout: 60000});
