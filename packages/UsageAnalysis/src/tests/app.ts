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
    for (const name of ['event_time', 'user', 'event_type', 'description', 'request_id'])
      expect(t.col(name) != null, true, `column "${name}" is missing`);
  });

  test('Refusals use the dialog labels', async () => {
    expect(CaptureView.refusal('subject.type is one of user, group'), 'Subject is one of user, group');
    expect(CaptureView.refusal('debugFlags: x not one of db'), 'Debug flags: x not one of db');
    expect(CaptureView.refusal('maxEvents is a whole number from 1'), 'Max events is a whole number from 1');
  });

  test('Timeline splits the server out of the source', async () => {
    const t = DG.DataFrame.fromColumns([DG.Column.fromStrings('source', ['client', 'server', 'A'])]);
    TimelineView.splitSource(t);
    expect(t.col('source')!.toList().join(','), 'client,server,server');
    expect(t.col('server')!.toList().join(','), ',,A');
  });

  test('Timeline of an unknown action is empty', async () => {
    const t: DG.DataFrame = await grok.functions.call('Timeline', {spec: JSON.stringify({action: '01M3NKD0PV60BG2W3EDPSC6ZZZ'})});
    expect(t.rowCount, 0);
    expect(t.col('requestId') != null, true, 'column "requestId" is missing');
  });
}, {timeout: 60000});

category('Errors', () => {
  const stats = async (spec: object): Promise<DG.DataFrame> =>
    await grok.functions.call('ErrorStats', {spec: JSON.stringify(spec)});

  test('ErrorStats occurrences', async () => {
    const t = await stats({since: '7d', limit: 5});
    for (const name of ['time', 'user', 'service', 'signature', 'error', 'package', 'version', 'route', 'server', 'requestId'])
      expect(t.col(name) != null, true, `column "${name}" is missing`);
    expect(t.rowCount <= 5, true);
  });

  test('ErrorStats by signature and package', async () => {
    const t = await stats({since: '7d', by: 'signature,package', trend: 'day'});
    for (const name of ['signature', 'package', 'count', 'users', 'firstVersion', 'firstSeen', 'trend', 'state', 'newInRange'])
      expect(t.col(name) != null, true, `column "${name}" is missing`);
    expect(t.col('count')!.type, DG.TYPE.INT);
    if (t.rowCount > 0)
      expect(t.get('trend', 0).split(' ').length, 7);
  });

  test('ErrorStats refuses a fourth dimension', async () => {
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

  test('Trend, state and schedule helpers', async () => {
    expect(ErrorsView.buckets('0 1 2 4 8').join(' '), '0 1 2 4 8');
    expect(ErrorsView.buckets('').length, 0);
    const summed = ErrorsView.buckets(Array(120).fill('1').join(' '));
    expect(summed.length, 60);
    expect(summed[0], 2);
    expect(ErrorsView.shortSignature('a41f9c3e-0000-0000-0000-000000000000'), 'a41f9c');
    expect(ErrorsView.cron('MON 07:00'), '0 7 * * 1');
    expect(ErrorsView.cron('weekdays 7:30'), '30 7 * * 1-5');
    expect(ErrorsView.cron('0 7 * * *'), '0 7 * * *');
    expect(ErrorsView.cron('someday'), null);
    const t = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('package', ['Chem', 'core']),
      DG.Column.fromStrings('route', ['GET /queries/{id}', 'socket']),
    ]);
    expect(JSON.stringify(ErrorsView.rowFilters(t, {by: 'package,route'}, 0)), '{"package":"Chem","route":"GET /queries/{id}"}');
    expect(JSON.stringify(ErrorsView.rowFilters(t, {by: 'package,route'}, 1)), '{"package":"core"}');
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
}, {timeout: 60000});
