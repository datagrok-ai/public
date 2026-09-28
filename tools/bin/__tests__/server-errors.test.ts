import {describe, it, expect, beforeEach, afterEach, vi} from 'vitest';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';
import {handleErrors, errorFilters, normalizeRoute, byArg, aggregateRow, occurrenceRow, classifyHosts, toParquet} from '../commands/server-errors';
import {cronFromSchedule, sparkline, shortRequestId, fmtMinutes} from '../utils/obs-format';
import {mockConnect, captureOutput, localIso} from './obs-helpers';

const NOW = new Date(2026, 8, 28, 11, 0);

describe('filters', () => {
  it('maps flags to query parameters, route without /api', () => {
    expect(errorFilters({since: '1h', route: 'POST /api/queries/{id}/run', 'min-users': 2, regressed: true, group: 'Chemists', service: 'client'}))
      .toEqual({since: '1h', route: 'POST /queries/{id}/run', minUsers: 2, regressed: true, group: 'Chemists', service: 'client'});
    expect(errorFilters({})).toEqual({since: '24h'});
    expect(errorFilters({regressed: false, package: 'Chem', version: '1.14.2'})).toEqual({since: '24h', package: 'Chem', version: '1.14.2'});
  });

  it('normalizes routes', () => {
    expect(normalizeRoute('POST /api/queries/{id}/run')).toBe('POST /queries/{id}/run');
    expect(normalizeRoute('/api/files/{path}')).toBe('/files/{path}');
    expect(normalizeRoute('GET /queries')).toBe('GET /queries');
    expect(normalizeRoute('/apix/y')).toBe('/apix/y');
  });

  it('takes --from/--to as absolute or relative times, never with --since', () => {
    expect(errorFilters({from: '-7d', to: '2026-09-27'}, {}, NOW))
      .toEqual({from: new Date(2026, 8, 21, 11, 0).toISOString(), to: new Date(2026, 8, 27).toISOString()});
    expect(() => errorFilters({since: '7d', from: '-1d'})).toThrow(/either --since/);
    expect(() => errorFilters({since: 'week'})).toThrow(/duration/);
    expect(() => errorFilters({service: 'db'})).toThrow(/server or client/);
  });

  it('refuses --from in a saved job, where --to is the destination', () => {
    expect(errorFilters({since: '7d', to: 'System:AppData/x/'}, {range: false})).toEqual({since: '7d'});
    expect(() => errorFilters({from: '-7d', to: 'System:AppData/x/'}, {range: false})).toThrow(/--since/);
  });

  it('validates --by', () => {
    expect(byArg('signature,package')).toEqual(['signature', 'package']);
    expect(byArg(undefined, 'signature')).toEqual(['signature']);
    expect(() => byArg('signature,team')).toThrow(/Unknown --by dimension 'team'/);
    expect(() => byArg('user,group,package,version')).toThrow(/at most three/);
  });
});

describe('formats', () => {
  it('draws one block per bucket, scaled to the busiest', () => {
    expect(sparkline([0, 1, 2, 3, 4, 5, 6, 7])).toBe('▁▂▃▄▅▆▇█');
    expect(sparkline([0, 0, 0])).toBe('▁▁▁');
    expect(sparkline(undefined)).toBe('');
  });

  it('shortens request ids to the action', () => {
    expect(shortRequestId('mfz3k2a1b9x8y7kq.3')).toBe('mfz3…kq');
    expect(shortRequestId('01J9ABCDEFGHJKMNPQRSTVWX7K')).toBe('01J9…7K');
    expect(shortRequestId(null)).toBe('');
  });

  it('prints a median time to resolve', () => {
    expect(fmtMinutes(190)).toBe('3 h 10 min');
    expect(fmtMinutes(12)).toBe('12 min');
    expect(fmtMinutes(null)).toBe('—');
  });

  it('converts schedules to cron and passes cron through', () => {
    expect(cronFromSchedule('MON 07:00')).toBe('0 7 * * 1');
    expect(cronFromSchedule('daily 18:30')).toBe('30 18 * * *');
    expect(cronFromSchedule('WEEKDAYS 06:05')).toBe('5 6 * * 1-5');
    expect(cronFromSchedule('0 7 * * 1')).toBe('0 7 * * 1');
    expect(() => cronFromSchedule('FUNDAY 07:00')).toThrow(/MON..SUN/);
    expect(() => cronFromSchedule('MON 25:00')).toThrow(/MON..SUN/);
    expect(() => cronFromSchedule('weekly')).toThrow(/five-field/);
  });

  it('renders signature rows with first seen, sparkline and state', () => {
    const r = {signature: 'a41f9c3e-1111-2222-3333-444455556666', package: 'Chem', users: 9, count: 57, firstVersion: '1.14.2',
      firstSeen: localIso(10, 2), lastSeen: localIso(10, 47), trend: [0, 0, 0, 0, 0, 0, 57], state: 'muted', stateVersion: '1.14.3',
      topError: 'TypeError: Cannot read properties of undefined (reading \'molfile\')'};
    const row = aggregateRow(r, ['package', 'signature'], 'day');
    expect(Object.keys(row)).toEqual(['SIG', 'PACKAGE', 'ERROR', 'USERS', 'COUNT', 'FIRST SEEN', 'LAST', 'TREND (daily)', 'STATE']);
    expect(row.SIG).toBe('a41f9c');
    expect(row['FIRST SEEN']).toMatch(/^1\.14\.2 · \d\d-\d\d$/);
    expect(row['TREND (daily)']).toBe('▁▁▁▁▁▁█');
    expect(row.STATE).toBe('muted → 1.14.3');
    expect(aggregateRow({...r, state: 'muted', stateVersion: null, stateUntil: localIso(14, 0)}, ['signature'], 'hour').STATE).toBe('muted until 14:00');
    expect(aggregateRow({...r, state: undefined}, ['signature'], 'day').STATE).toBe('');
  });

  it('renders rollup rows when signature is not a dimension', () => {
    const row = aggregateRow({connection: 'Snowflake:PROD', signatures: 1, users: 14, count: 61, newInRange: 1, firstVersion: '1.28.2',
      trend: [1, 2], topError: 'SQLTimeoutException: warehouse suspended'}, ['connection'], 'hour');
    expect(Object.keys(row)).toEqual(['CONNECTION', 'SIGNATURES', 'USERS', 'COUNT', 'NEW', 'FIRST SEEN IN', 'LAST', 'TREND (hourly)', 'TOP ERROR']);
    expect(row['TOP ERROR']).toBe('SQLTimeoutException: warehouse suspended');
  });

  it('renders occurrence rows', () => {
    const row = occurrenceRow({time: localIso(10, 14), user: 'alice', service: 'client', signature: 'a41f9c3e', error: 'boom',
      package: 'Chem', version: '1.14.2', route: 'POST /projects/{id}/save', server: 'datlas-1', requestId: 'mfz3k2a1b9x8y7kq.2'});
    expect(row).toEqual({TIME: '10:14', USER: 'alice', SOURCE: 'client', SIG: 'a41f9c', ERROR: 'boom', PACKAGE: 'Chem',
      VERSION: '1.14.2', ROUTE: 'POST /projects/{id}/save', SERVER: 'datlas-1', REQ: 'mfz3…kq'});
  });
});

describe('classifyHosts', () => {
  it('splits signatures into only-on-a, only-on-b and both, ranked by users', () => {
    const a = [{signature: 's1', users: 2}, {signature: 's2', users: 1}, {signature: 's9', users: 5}];
    const b = [{signature: 's2', users: 3}, {signature: 's3', users: 6}, {signature: 's4', users: 1}];
    const {onlyA, onlyB, both} = classifyHosts(a, b);
    expect(onlyA.map((r) => r.signature)).toEqual(['s9', 's1']);
    expect(onlyB.map((r) => r.signature)).toEqual(['s3', 's4']);
    expect(both).toEqual(['s2']);
  });
});

describe('parquet', () => {
  it('names the missing packages when they are not installed', () => {
    const missing = () => { throw new Error('Cannot find module'); };
    expect(() => toParquet([{a: 1}], missing)).toThrow(/apache-arrow and parquet-wasm/);
  });

  it('flattens nested values and writes through arrow IPC', () => {
    const seen: any = {};
    const arrow = {tableFromJSON: (rows: any[]) => { seen.rows = rows; return 'T'; }, tableToIPC: (t: any, f: string) => { seen.ipc = [t, f]; return 'IPC'; }};
    const parquet = {Table: {fromIPCStream: (b: any) => ({b})}, writeParquet: (t: any) => new Uint8Array([80, 65, 82, 49, t.b.length])};
    const bytes = toParquet([{sig: 'a', trend: [1, 2]}], (n) => n === 'apache-arrow' ? arrow : parquet);
    expect(seen.rows).toEqual([{sig: 'a', trend: '[1,2]'}]);
    expect(seen.ipc).toEqual(['T', 'stream']);
    expect([...bytes]).toEqual([80, 65, 82, 49, 3]);
  });
});

describe('handleErrors', () => {
  beforeEach(() => { vi.useFakeTimers({toFake: ['Date']}); vi.setSystemTime(NOW); });
  afterEach(() => vi.useRealTimers());

  it('asks top for the grouping, trend and limit and prints the table', async () => {
    const {connect, calls} = mockConnect(() => [{signature: '7c02e1aa', package: 'PowerGrid', users: 6, count: 212, firstVersion: '2.3.0',
      firstSeen: localIso(9, 0, 17), lastSeen: localIso(10, 0), trend: [3, 4, 3, 5, 4, 6, 5], state: 'open', topError: 'NPE'}]);
    const {out} = await captureOutput(() => handleErrors(connect, 'top', [], {since: '7d', by: 'signature,package', 'min-users': 2, limit: 5}, 'table'));
    expect(calls[0].path).toBe('/errors?since=7d&minUsers=2&by=signature%2Cpackage&trend=day&limit=5');
    expect(out[0].split(/\s{2,}/)).toEqual(['SIG', 'PACKAGE', 'ERROR', 'USERS', 'COUNT', 'FIRST SEEN', 'LAST', 'TREND (daily)', 'STATE']);
    expect(out[2]).toMatch(/^7c02e1\s+PowerGrid\s+NPE\s+6\s+212\s+2\.3\.0 · 09-11\s+10:00\s+▅▆▅▇▆█▇\s+open/);
  });

  it('exports csv from the server and json/parquet from the rows, to stdout or -O', async () => {
    const {connect, calls} = mockConnect((_m, p) => p.includes('format=csv') ? 'package,count\nChem,57\n' : [{package: 'Chem', count: 57}]);
    const csv = await captureOutput(() => handleErrors(connect, 'top', [], {since: '7d', by: 'package,group', format: 'csv'}, 'table'));
    expect(csv.out.join('')).toBe('package,count\nChem,57\n');
    const file = path.join(os.tmpdir(), `grok-errors-${process.pid}.json`);
    const json = await captureOutput(() => handleErrors(connect, 'export', [], {format: 'json', O: file}, 'table'));
    expect(JSON.parse(fs.readFileSync(file, 'utf8'))).toEqual([{package: 'Chem', count: 57}]);
    expect(json.out[0]).toMatch(/^Wrote \d+ bytes to /);
    fs.unlinkSync(file);
    expect(calls.map((c) => c.path)).toEqual([
      '/errors?since=7d&by=package%2Cgroup&trend=day&limit=20&format=csv',
      '/errors?since=24h&format=json',
    ]);
    await expect(handleErrors(connect, 'export', [], {format: 'xlsx'}, 'table')).rejects.toThrow(/csv, json or parquet/);
  });

  it('prints the show block', async () => {
    const {connect, calls} = mockConnect(() => ({signature: 'a41f9c3e', error: 'TypeError: x', package: 'Chem', firstVersion: '1.14.2',
      firstSeen: localIso(10, 2), lastSeen: localIso(10, 47), occurrences: 57, users: 9,
      groups: [{name: 'Chemists', users: 7}, {name: 'Biology', users: 2}], reports: [4815, 4817],
      alert: {kind: 'error-incident', key: 'a41f9c', status: 'open', openedAt: localIso(10, 5)},
      change: {type: 'package-published', package: 'Chem', version: '1.14.2', by: 'j.doe', at: localIso(9, 58)}}));
    const {out} = await captureOutput(() => handleErrors(connect, 'show', ['a41f9c'], {}, 'table'));
    expect(calls[0].path).toBe('/errors/a41f9c?since=24h');
    expect(out).toEqual([
      'signature    a41f9c  TypeError: x',
      'package      Chem    first seen in 1.14.2 at 10:02 · last seen 10:47',
      'occurrences  57      users 9      groups Chemists 7 · Biology 2',
      'reports      #4815 #4817          alert error-incident, open since 10:05',
      'change       package-published Chem 1.14.2 by j.doe at 09:58',
    ]);
  });

  it('prints the window diff', async () => {
    const top = (sig: string, extra: any = {}) => ({signature: sig, package: 'core', error: 'TableView.close', users: 4, ...extra});
    const {connect, calls} = mockConnect(() => ({
      categories: {new: {count: 12, top: top('3e91aa')}, gone: {count: 8, top: top('91bb0c')},
        risen: {count: 3, top: top('7c02e1', {package: 'PowerGrid', before: 64, after: 212})}, regressed: {count: 0}},
      incidents: {opened: 5, resolved: 5, medianTtrMinutes: 190}, reports: {total: 22, human: 9, auto: 13, linkedToNew: 6}, rows: []}));
    const {out} = await captureOutput(() => handleErrors(connect, 'diff', [],
      {before: '2026-09-14..2026-09-20', after: '2026-09-21..2026-09-27', package: 'core'}, 'table'));
    expect(calls[0].path).toBe('/errors/diff?package=core&before=2026-09-14..2026-09-20&after=2026-09-21..2026-09-27');
    expect(out[0]).toMatch(/^NEW\s+12\s+first seen in 2026-09-21..2026-09-27\s+top: 3e91aa core TableView.close · 4 users$/);
    expect(out[2]).toMatch(/^RISEN\s+3\s+≥ 2× occurrences\s+top: 7c02e1 PowerGrid · 212 vs 64$/);
    expect(out[3]).toMatch(/^REGRESSED\s+0\s+muted, back on a newer version$/);
    expect(out[4]).toBe('incidents 5 opened · 5 resolved · median time to resolve 3 h 10 min');
    expect(out[5]).toBe('reports   22 (human 9, auto 13) · 6 linked to a NEW signature');
    await expect(handleErrors(connect, 'diff', [], {before: '2026-09-14', after: 'x..y'}, 'table')).rejects.toThrow(/<from>..<to>/);
  });

  it('compares two deployments', async () => {
    const rows: Record<string, any[]> = {
      prod: [{signature: '91bb0c11', package: 'core', users: 2, topError: 'ScatterPlot.render'}, {signature: 'bbbb', users: 1}],
      val: [{signature: '3e91aa22', package: 'core', users: 6, topError: 'TableView.close'}, {signature: 'bbbb', users: 1}],
    };
    const {connect, calls} = mockConnect((_m, p, _b, host) => p === '/info/server' ? {Version: host === 'prod' ? '1.28.1' : '1.28.3'} : rows[host]);
    const {out} = await captureOutput(() => handleErrors(connect, 'diff', [], {since: '7d', host: ['prod', 'val']}, 'table'));
    expect(calls.filter((c) => c.path.startsWith('/errors')).map((c) => c.path))
      .toEqual(['/errors?since=7d&by=signature%2Cpackage&limit=10000', '/errors?since=7d&by=signature%2Cpackage&limit=10000']);
    expect(out[0]).toMatch(/^ONLY ON prod \(1\.28\.1\)\s+1\s+top: 91bb0c core ScatterPlot.render · 2 users$/);
    expect(out[1]).toMatch(/^ONLY ON val \(1\.28\.3\)\s+1\s+top: 3e91aa core TableView.close · 6 users$/);
    expect(out[2]).toMatch(/^ON BOTH\s+1$/);
    await expect(handleErrors(connect, 'diff', [], {host: ['a', 'b', 'c']}, 'table')).rejects.toThrow(/exactly two/);
  });

  it('saves a scheduled job into a folder', async () => {
    const {connect, calls} = mockConnect(() => ({id: 'J1', name: 'Errors by team, weekly'}));
    const {out} = await captureOutput(() => handleErrors(connect, 'save', ['Errors by team, weekly'],
      {since: '7d', by: 'group,package', schedule: 'MON 07:00', to: 'System:AppData/Ops/errors/'}, 'table'));
    expect(calls[0]).toMatchObject({method: 'POST', path: '/errors/jobs', body: {name: 'Errors by team, weekly',
      spec: {since: '7d', by: 'group,package'}, format: 'csv', path: 'System:AppData/Ops/errors/errors-by-team-weekly-{date}.csv', cron: '0 7 * * 1'}});
    expect(out).toEqual(['saved job "Errors by team, weekly" 0 7 * * 1 → System:AppData/Ops/errors/errors-by-team-weekly-{date}.csv']);
    await expect(handleErrors(connect, 'save', ['x'], {to: 'System:AppData/x.parquet', format: 'parquet'}, 'table')).rejects.toThrow(/csv or json/);
  });

  it('lists occurrences from several hosts', async () => {
    const {connect, calls} = mockConnect(() => [{time: localIso(10, 0), user: 'alice', signature: 'a41f9c3e', error: 'boom'}]);
    const {out} = await captureOutput(() => handleErrors(connect, 'list', [], {host: ['prod', 'val'], user: 'alice'}, 'table'));
    expect(calls.map((c) => c.path)).toEqual(['/errors?since=24h&user=alice&limit=50', '/errors?since=24h&user=alice&limit=50']);
    expect(out[0].split(/\s{2,}/).slice(0, 3)).toEqual(['HOST', 'TIME', 'USER']);
  });
});
