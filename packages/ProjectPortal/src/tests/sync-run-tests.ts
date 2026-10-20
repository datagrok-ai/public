import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {loadAllPortfolios, portfolioByCode, readPortfolioFiles} from '../portal';
import {validatePortfolio} from '../model/config';
import {
  describeChange, describeCounts, describeSchedule, isDue, nextRun, parseSchedule, standingIssues, summarizeRun,
  validateSchedule, worthNotifying,
} from '../model/sync-runs';
import {SyncPlan} from '../model/sync-plan';
import {dueState, notifyOwners, recordedSync, syncDue, syncRunsOf} from '../sync-runs';
import {syncHistoryPanel} from '../views/sync-history';
import {_package} from '../package-test';

const at = (iso: string) => Date.parse(iso);
const iso = (t: number) => new Date(t).toISOString().slice(0, 16);

const plan = (over: Partial<SyncPlan> = {}): SyncPlan => ({portfolio: 'P', system: 'S', inserts: [], updates: [], unchanged: [],
  stageChanges: [], unresolvedStages: [], missing: [], lineage: [], overlaps: [], rejected: [], warnings: [], ...over});

category('Portal: sync schedule rules', () => {
  test('a schedule is an interval or a five-field cron line; anything else is refused with the reason', async () => {
    expect(JSON.stringify(['every 6h', '30m', '1 day', 'every 90 minutes'].map((s) => (parseSchedule(s) as any).ms / 60000)),
      '[360,30,1440,90]');
    expect(parseSchedule('0 6 * * 1-5').kind, 'cron');
    expect(parseSchedule('@daily').kind, 'cron');
    const reason = (s: string) => {
      try {
        parseSchedule(s);
        return '';
      } catch (e: any) {
        return e.message;
      }
    };
    expect(reason('5m'), '"5m" is shorter than 15 minutes');
    expect(reason('0 6 * *'), '"0 6 * *" is neither an interval (every 6h) nor a five-field cron line (0 6 * * 1-5)');
    expect(reason('0 25 * * *'), 'hour "25" is outside 0-23');
    expect(reason('0 6 * * mon'), 'day of week "mon" is not a number, range, list or step');
    expect(describeSchedule(parseSchedule('6h')), 'every 6h');
    expect(describeSchedule(parseSchedule('0 6 * * 1-5')), 'cron 0 6 * * 1-5 (UTC)');
  });

  test('cron in UTC: the next minute it names, both day fields either-or', async () => {
    const next = (s: string, after: string) => iso(nextRun(parseSchedule(s), at(after), 0));
    expect(next('30 6 * * 1-5', '2026-10-02T07:00:00Z'), '2026-10-05T06:30', 'Friday after 06:30 → Monday');
    expect(next('0 0 1 * 1', '2026-10-02T00:00:00Z'), '2026-10-05T00:00', 'a Monday before the 1st');
    expect(next('@daily', '2026-10-02T23:59:30Z'), '2026-10-03T00:00');
    expect(next('*/15 * * * *', '2026-10-02T10:07:00Z'), '2026-10-02T10:15');
    expect(next('0 0 29 2 *', '2026-03-01T00:00:00Z'), '2028-02-29T00:00', 'a leap day');
  });

  test('due: never run, or the next slot after the last start has come', async () => {
    const daily = parseSchedule('0 6 * * *');
    expect(isDue(daily, null, at('2026-10-02T01:00:00Z')), true, 'never run');
    expect(isDue(daily, at('2026-10-02T06:00:00Z'), at('2026-10-03T05:59:00Z')), false);
    expect(isDue(daily, at('2026-10-02T06:00:00Z'), at('2026-10-03T06:00:00Z')), true);
    const six = parseSchedule('every 6h');
    expect(isDue(six, at('2026-10-02T06:00:00Z'), at('2026-10-02T11:59:00Z')), false);
    expect(iso(nextRun(six, at('2026-10-02T06:00:00Z'), 0)), '2026-10-02T12:00');
  });

  test('a run\'s counts and per-program changes', async () => {
    const r = summarizeRun({
      programs: [plan({inserts: [{code: 'P1', name: 'One'} as any], updates: [{id: 'i2', code: 'P2', changes: {target: 'T', name: 'N'}}],
        stageChanges: [{code: 'P2', programId: 'i2', fromStageId: 's1', toStageId: 's2'}], missing: [{id: 'i3', code: 'P3', status: 'missing_from_feed'}],
        unresolvedStages: [{code: 'P4', spelling: 'Lead Gen'}], overlaps: [{code: 'P5', otherPortfolio: 'Q'}]})],
      units: [{portfolio: 'P', kind: 'series', feed: 'f', inserts: [{program_id: 'i2', code: 'S1', name: 'S1', details: null}],
        updates: [{id: 'u', key: 'P2/S2', changes: {}}], unchanged: [], missing: [{id: 'm', key: 'P2/S3'}], returned: [], rejected: [], warnings: []}],
      compounds: {portfolio: 'P', system: 'S', compoundInserts: [{code: 'C1'} as any], compoundUpdates: [],
        linkInserts: [{programId: 'i2', programCode: 'P2', compoundCode: 'C1', role: 'hit'}],
        roleChanges: [{linkId: 'l', programId: 'i2', programCode: 'P2', compoundCode: 'C2', from: 'hit', to: 'lead'}],
        linkMissing: [], linkReturned: [], unchanged: 0, rejected: [], warnings: []},
      programCodes: new Map([['i2', 'P2']]),
    });
    expect(JSON.stringify(r.counts), JSON.stringify({created: 4, updated: 3, flagged: 2, unresolved: 1,
      unresolvedSpellings: ['P4="Lead Gen"'], overlaps: ['P5 also in Q']}));
    expect(r.changes.map((c) => c.program).join(), 'P1,P2,P3');
    expect(describeChange(r.changes[1]), 'P2: changed name, target; stage moved; new series S1; updated series S2; ' +
      'missing series S3; compounds added C1; C2 hit → lead');
    expect(describeCounts(r.counts), '4 created, 3 updated, 2 flagged, 1 unresolved spelling (P4="Lead Gen"), 1 overlap (P5 also in Q)');
    expect(describeCounts(r.counts, 'boom'), 'failed: boom');
  });

  test('owners hear about a change, a failure or a new standing issue; not the same run again', async () => {
    const quiet = {created: 0, updated: 0, flagged: 0, unresolved: 1, unresolvedSpellings: ['X="y"'], overlaps: []};
    expect(worthNotifying({counts: quiet}, standingIssues(quiet)), false);
    expect(worthNotifying({counts: quiet}, null), true, 'the first run tells');
    expect(worthNotifying({counts: quiet}, []), true, 'a spelling that was not there before');
    expect(worthNotifying({counts: {...quiet, updated: 1}}, standingIssues(quiet)), true);
    expect(worthNotifying({counts: quiet, error: 'down'}, standingIssues(quiet)), true);
  });

  test('schedule and owners in the portfolio file are checked; the examples declare both', async () => {
    expect(validateSchedule('every 2m', ['a|b']).join('; '),
      'schedule: "every 2m" is shorter than 15 minutes; owner "a|b" is not a user login');
    const all = await readPortfolioFiles(_package);
    for (const c of all.filter((x) => x.code.startsWith('EX-'))) {
      expect(validatePortfolio(c).join('; '), '', c.code);
      expect(!!c.schedule && (c.owners ?? []).length > 0, true, `${c.code} declares a schedule and owners`);
    }
    expect(all.find((c) => c.code === 'EX-SMD')!.schedule, 'every 6h');
  });
});

category('Portal: scheduled sync', () => {
  before(async () => {
    await loadAllPortfolios(_package);
  });

  test('a sync by hand is recorded: trigger, who, finished, counts, the report', async () => {
    const r = await recordedSync('EX-SMD', 'manual');
    expect(r.error, null);
    expect(r.historyError ?? '', '');
    const smd = await portfolioByCode('EX-SMD');
    expect(smd.schedule, 'every 6h', 'the schedule is loaded with the portfolio');
    const last = (await syncRunsOf(smd, 1))[0];
    expect(last.id, r.runId!);
    expect(`${last.trigger} ${last.run_by === DG.User.current().id} ${!!last.finished_on}`, 'manual true true');
    expect(last.summary, r.report);
    expect(r.report.startsWith('EX-SMD: '), true, r.report);
  });

  test('syncDuePortfolios syncs only what is due, by schedule, and says why it skipped the rest', async () => {
    await recordedSync('EX-SMD', 'manual');
    const smd = await portfolioByCode('EX-SMD');
    const soon = await syncDue(['EX-SMD', 'NOPE'], Date.now() + 3600 * 1000);
    expect(soon.ran.length, 0, 'an hour after a run, a 6-hour schedule is not due');
    expect(soon.skipped.some((s) => s.startsWith('EX-SMD: next due')), true, soon.skipped.join('; '));
    expect(soon.skipped.includes('NOPE: not loaded'), true, soon.skipped.join('; '));
    const later = await syncDue(['EX-SMD'], Date.now() + 7 * 3600 * 1000);
    expect(later.ran.map((r) => `${r.portfolio} ${r.trigger} ${r.error}`).join(), 'EX-SMD schedule null');
    expect(!!later.ran[0].notified, true);
    const last = (await syncRunsOf(smd, 1))[0];
    expect(`${last.id === later.ran[0].runId} ${last.trigger}`, 'true schedule');
    const state = await dueState(smd);
    expect(state.due, false, 'just ran');
  });

  test('the platform notification reaches an owner; an unknown login is reported', async () => {
    const me = DG.User.current();
    const text = `Project Portal test notification ${Date.now()}`;
    const r = await notifyOwners({owners: `${me.login}|no-such-login-xyz`} as any, text);
    expect(r, `${me.login}: notified; no-such-login-xyz: no such user`);
    // the platform's js-api on this stand has no notifications client: read them over REST
    const call = (path: string, method = 'GET') => fetch(`${grok.dapi.root}/users/notifications/current${path}`,
      {method, credentials: 'include', headers: {Authorization: grok.dapi.token}});
    const mine: {id: string, text: string}[] = await (await call('?limit=5000')).json();
    const hit = mine.find((n) => n.text === text);
    expect(!!hit, true, `the notification is in the current user's list (${mine.length} read)`);
    if (hit)
      await call(`/${hit.id}/read`, 'POST');
    expect(await notifyOwners({owners: null} as any, 'x'), 'no owners declared');
  });

  test('the function and the history panel show the runs, newest first', async () => {
    await grok.functions.call('ProjectPortal:syncPortfolio', {code: 'EX-BIOD', dryRun: false});
    const df: DG.DataFrame = await grok.functions.call('ProjectPortal:syncHistory', {portfolio: 'EX-BIOD'});
    expect(df.rowCount > 0, true);
    expect(df.get('trigger', 0), 'manual');
    const panel = await syncHistoryPanel(await portfolioByCode('EX-BIOD'));
    expect(panel.querySelector('[data-u2-name="sync-schedule"]')!.textContent!.startsWith('Scheduled: cron 0 6 * * * (UTC)'), true);
    expect(!!panel.querySelector('[data-u2-name="sync-run-0"]'), true);
  });
});
