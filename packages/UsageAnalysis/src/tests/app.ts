import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {before, after, category, test, expect, waitUntil} from '@datagrok-libraries/test/src/test';
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

  const tabs = ['Overview', 'Packages', 'Functions', 'Events', 'Clicks', 'Log', 'System Activity', 'Errors', 'Projects'];

  for (const tab of tabs) {
    test(tab, async () => {
      handler.changeTab(tab);
      const view = handler.getCurrentView();
      expect(view.name === tab, true, `expected tab "${tab}", got "${view.name}"`);
      await waitUntil(() => view.root.children.length > 0, `"${tab}" failed to initialize`, 30000);
    });
  }

  test('Errors: top errors and their incident', async () => {
    const errors = (await grok.dapi.admin.getMetrics({date: 'this month', limit: 5, errorsBy: 'signature'})).errors;
    expect(errors.now.count >= errors.top.length, true);
    for (const e of errors.top)
      expect(typeof e.signature, 'string');
    const signature = errors.top[0]?.signature ?? '000000000000';
    const problem = await grok.dapi.admin.getProblem(`error-incident:${signature.substring(0, 12)}`);
    expect(problem === null || problem.kind === 'error-incident', true);
  });

  test('Log: the timeline of the current session', async () => {
    const session = (await grok.dapi.users.currentSession()).id;
    const rows = await grok.dapi.log.getTimeline({session, from: new Date(Date.now() - 3600000), to: new Date(), limit: 5});
    expect(Array.isArray(rows), true);
  });
  test('a project reopens the tab, its sub-view and the filters', async () => {
    // a tag filter is shown on the Functions tab only, so it is applied only once that tab is current
    const tag: string | undefined = (handler as any).toolbox?.tagsInput?.choices?._store?.choices?.[0]?.value;
    const reopened = new ViewHandler();
    reopened.view.loadStateMap({tab: 'Functions', rout: '/ExecutionTime', filters: {date: 'this month', ...(tag ? {tags: [tag]} : {})}});
    expect(reopened.saveState()?.rout, '/ExecutionTime', 'a save while the app loads keeps the state it restores');
    grok.shell.addView(reopened.view);
    try {
      await reopened.init();
      await waitUntil(() => reopened.getCurrentView().name === 'Functions' && reopened.getCurrentView().rout === '/ExecutionTime',
        'the tab and its sub-view were not restored', 30000);
      expect(reopened.view.path.includes('date=this%20month'), true, `the date filter was not applied: ${reopened.view.path}`);
      if (tag)
        expect(reopened.view.path.includes('tags='), true, `the tag filter of the Functions tab was not applied: ${reopened.view.path}`);
    } finally {
      reopened.view.close();
    }
  });
}, {clear: false, timeout: 60000});
