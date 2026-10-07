import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {before, after, category, test, expect, awaitCheck} from '@datagrok-libraries/test/src/test';
import {ViewHandler} from '../view-handler';
import {queries} from '../package-api';
import {ErrorsView} from '../tabs/errors';
import {CaptureView} from '../tabs/capture';
import {TimelineView} from '../tabs/timeline';


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

  const tabs = ['Overview', 'Packages', 'Functions', 'Events', 'Clicks', 'Log', 'System Activity', 'Errors', 'Capture',
    'Timeline', 'Projects'];

  for (const tab of tabs) {
    test(tab, async () => {
      handler.changeTab(tab);
      const view = handler.getCurrentView();
      expect(view.name === tab, true, `expected tab "${tab}", got "${view.name}"`);
      await awaitCheck(() => view.root.children.length > 0, `"${tab}" failed to initialize`, 30000);
    });
  }

  test('Timeline tab of a real click', async () => {
    const allUsers = (await grok.dapi.groups.getGroupsLookup('All users'))[0].id;
    const clicks = await queries.clicks('this month', [allUsers]);
    let key = 'action';
    let id = clicks.col('action id')!.toList().find((v) => v);
    if (!id) {
      const errors = await grok.dapi.log.getErrors({since: '30d', limit: 50});
      key = 'request';
      id = errors.map((e) => e.requestId).find((v) => v);
    }
    expect(id != null, true, 'no click or error with a request id this month');
    handler.getCurrentView().openTimeline(key, id);
    const view = handler.getCurrentView() as TimelineView;
    await awaitCheck(() => view.host.querySelector('.d4-grid, .d4-viewer-error, .ua-empty') != null,
      'Timeline did not load', 30000);
    expect(view.host.querySelector('.d4-viewer-error')?.textContent ?? '', '');
    expect(view.host.querySelector('.d4-grid') != null, true, `no records for ${key} ${id}`);
  });
}, {clear: false, timeout: 60000});

category('Capture', () => {
  test('CaptureRules columns', async () => {
    const t = await queries.captureRules('this year');
    for (const name of ['rule', 'author', 'subject', 'scope', 'reason', 'active', 'events', 'status', 'id', 'stopped_by',
      'stop_reason'])
      expect(t.col(name) != null, true, `column "${name}" is missing`);
  });

  test('Clicks has request ids', async () => {
    const allUsers = (await grok.dapi.groups.getGroupsLookup('All users'))[0].id;
    const t = await queries.clicks('this week', [allUsers]);
    for (const name of ['time', 'user', 'type', 'element', 'action id'])
      expect(t.col(name) != null, true, `column "${name}" is missing`);
  });

  test('Refusals use the dialog labels', async () => {
    expect(CaptureView.refusal('subject.type is one of user, group'), 'Subject is one of user, group');
    expect(CaptureView.refusal('debugFlags: x not one of db'), 'Debug flags: x not one of db');
    expect(CaptureView.refusal('maxEvents is a whole number from 1'), 'Max events is a whole number from 1');
  });

  test('Debug flags come from the server, without credentials', async () => {
    const flags = await CaptureView.loadDebugFlags();
    expect(flags.includes('query'), true, `no "query" in ${flags}`);
    expect(flags.includes('credentials'), false, 'credentials is offered');
  });

  test('Timeline of an unknown action is empty', async () => {
    expect((await grok.dapi.log.getTimeline({action: '01M3NKD0PV60BG2W3EDPSC6ZZZ'})).length, 0);
  });
}, {timeout: 60000});

category('Errors', () => {
  const stats = async (spec: {[key: string]: string | number}): Promise<DG.DataFrame> =>
    ErrorsView.frame(await grok.dapi.log.getErrors(spec));

  test('Errors occurrences', async () => {
    const t = await stats({since: '7d', limit: 5});
    expect(t.rowCount <= 5, true);
    if (t.rowCount > 0) {
      for (const name of ['time', 'user', 'service', 'signature', 'error', 'package', 'version', 'route', 'server', 'requestId'])
        expect(t.col(name) != null, true, `column "${name}" is missing`);
    }
  });

  test('Errors by signature and package', async () => {
    const t = await stats({since: '7d', by: 'signature,package', trend: 'day'});
    if (t.rowCount === 0)
      return;
    for (const name of ['signature', 'package', 'count', 'users', 'firstVersion', 'firstSeen', 'trend', 'state', 'newInRange'])
      expect(t.col(name) != null, true, `column "${name}" is missing`);
    expect(t.col('count')!.type, DG.TYPE.INT);
    expect(t.get('trend', 0).split(' ').length, 7);
  });

  test('Errors refuses a fourth dimension', async () => {
    let message = '';
    try {
      await stats({since: '7d', by: 'signature,package,user,server'});
    }
    catch (e: any) {
      message = `${e?.message ?? e}`;
    }
    expect(message.includes('up to three'), true, `unexpected: "${message}"`);
  });

  test('Drill-down queries', async () => {
    const none = ['00000000-0000-0000-0000-000000000000'];
    expect((await queries.errorSessions(none, ['admin'], '2026-01-01T00:00:00Z', '2100-01-01T00:00:00Z')).rowCount, 0);
    expect((await queries.errorReports(none)).col('number') != null, true);
    expect((await queries.errorAlerts(none)).col('status') != null, true);
    expect((await queries.errorSample(none[0])).rowCount, 0);
  });

  test('ClicksFollowedByError columns', async () => {
    const allUsers = (await grok.dapi.groups.getGroupsLookup('All users'))[0].id;
    const t = await queries.clicksFollowedByError('this month', [allUsers], '');
    for (const name of ['element', 'clicks', 'users', 'followed_by_error', 'followed_by_error_pct'])
      expect(t.col(name) != null, true, `column "${name}" is missing`);
    const errors = await queries.clickErrors('this month', [allUsers], 'no such element', '');
    for (const name of ['signature', 'error', 'clicks', 'action'])
      expect(errors.col(name) != null, true, `column "${name}" is missing`);
  });

  test('Trend and state helpers', async () => {
    expect(ErrorsView.buckets('0 1 2 4 8').join(' '), '0 1 2 4 8');
    expect(ErrorsView.buckets('').length, 0);
    const summed = ErrorsView.buckets(Array(120).fill('1').join(' '));
    expect(summed.length, 60);
    expect(summed[0], 2);
    expect(ErrorsView.shortSignature('a41f9c3e-0000-0000-0000-000000000000'), 'a41f9c');
    const logins = [...Array(12).keys()].map((k) => `u${k}`);
    expect(ErrorsView.usersText(logins), 'u0, u1, u2, u3, u4, u5, u6, u7, u8, u9 +2 more');
    expect(ErrorsView.usersText([null, '']), 'none');
    expect(ErrorsView.emptyHint({since: '7d', by: 'user', minUsers: 2, group: undefined}),
      'Min users is 2. Clear it, or choose a longer Since');
    expect(ErrorsView.emptyHint({since: '7d'}), 'Choose a longer Since');
    const t = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('package', ['Chem', 'core']),
      DG.Column.fromStrings('route', ['GET /queries/{id}', 'socket']),
    ]);
    expect(JSON.stringify(ErrorsView.rowFilters(t, {by: 'package,route'}, 0)), '{"package":"Chem","route":"GET /queries/{id}"}');
    expect(JSON.stringify(ErrorsView.rowFilters(t, {by: 'package,route'}, 1)), '{"package":"core"}');
  });

  test('Rows become typed columns', async () => {
    const t = ErrorsView.frame([{signature: 'a41f9c', count: 3, lastSeen: '2026-10-01T10:00:00.000Z', trend: [1, 2],
      regressed: true, mttrMinutes: 4}]);
    expect(t.col('count')!.type, DG.TYPE.INT);
    expect(t.col('lastSeen')!.type, DG.TYPE.DATE_TIME);
    expect(t.col('regressed')!.type, DG.TYPE.BOOL);
    expect(t.col('mttrMinutes')!.type, DG.TYPE.FLOAT);
    expect(t.get('trend', 0), '1 2');
    expect(ErrorsView.frame([]).rowCount, 0);
  });

  test('CSV export neutralises formulas', async () => {
    const t = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('error', ['=1+2', '+x', '-y', '@z', 'plain']),
      DG.Column.fromInt32Array('count', new Int32Array([1, 2, 3, 4, 5])),
    ]);
    const lines = ErrorsView.toCsv(t).trim().split('\n');
    expect(lines.slice(1).map((l) => l.split(',')[0]).join(' '), '\'=1+2 \'+x \'-y \'@z plain');
    expect(t.get('error', 0), '=1+2');
  });

  test('Export takes the grid\'s columns and headers', async () => {
    const grid = DG.Viewer.grid(DG.DataFrame.fromColumns([
      DG.Column.fromStrings('signature', ['a41f9c3e-0000-0000-0000-000000000000']),
      DG.Column.fromInt32Array('count', new Int32Array([7])),
      DG.Column.fromStrings('requestId', ['x']),
    ]));
    grid.col('count')!.name = 'occurrences';
    grid.col('requestId')!.visible = false;
    const t = ErrorsView.exportTable(grid);
    expect(t.columns.names().join(','), 'signature,occurrences');
    expect(t.get('signature', 0), 'a41f9c');
    expect(t.get('occurrences', 0), 7);
  });
}, {timeout: 60000});
