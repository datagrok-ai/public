import {describe, it, expect, beforeEach, afterEach, vi} from 'vitest';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';
import {handleErrors, errorFilters, normalizeRoute, byArg, aggregateRow, occurrenceRow, classifyHosts, toParquet} from '../commands/server-errors';
import {mockConnect, captureOutput, utcIso} from './obs-helpers';

const NOW = new Date(Date.UTC(2026, 8, 28, 11, 0));

describe('filters', () => {
  it('maps flags to query parameters, route without /api', () => {
    expect(errorFilters({since: '1h', route: 'POST /api/queries/{id}/run', 'min-users': 2, regressed: true, group: 'Chemists', service: 'client'}))
      .toEqual({since: '1h', route: 'POST /queries/{id}/run', minUsers: 2, regressed: true, group: 'Chemists', service: 'client'});
    expect(errorFilters({})).toEqual({since: '24h'});
    expect(errorFilters({regressed: false, package: 'Chem', version: '1.14.2'})).toEqual({since: '24h', package: 'Chem', version: '1.14.2'});
    expect(errorFilters({since: '-7d', signature: ''})).toEqual({since: '7d'});
  });

  it('normalizes routes', () => {
    expect(normalizeRoute('POST /api/queries/{id}/run')).toBe('POST /queries/{id}/run');
    expect(normalizeRoute('/api/files/{path}')).toBe('/files/{path}');
    expect(normalizeRoute('GET /queries')).toBe('GET /queries');
    expect(normalizeRoute('/apix/y')).toBe('/apix/y');
  });

  it('takes --from/--to as absolute or relative times, never with --since', () => {
    expect(errorFilters({from: '-7d', to: '2026-09-27'}, {}, NOW))
      .toEqual({from: new Date(Date.UTC(2026, 8, 21, 11, 0)).toISOString(), to: new Date(Date.UTC(2026, 8, 27)).toISOString()});
    expect(() => errorFilters({since: '7d', from: '-1d'})).toThrow(/either --since/);
    expect(() => errorFilters({since: 'week'})).toThrow(/duration/);
    expect(() => errorFilters({service: 'db'})).toThrow(/server or client/);
    for (const v of ['many', true, 0, 1.5])
      expect(() => errorFilters({'min-users': v})).toThrow(/--min-users is a whole number/);
    expect(() => errorFilters({'min-count': 'x'})).toThrow(/--min-count is a whole number/);
  });

  it('validates --by', () => {
    expect(byArg('signature,package')).toEqual(['signature', 'package']);
    expect(byArg(undefined, 'signature')).toEqual(['signature']);
    expect(() => byArg('signature,team')).toThrow(/Unknown --by dimension 'team'/);
    expect(() => byArg('user,group,package,version')).toThrow(/at most three/);
  });

});

describe('rows', () => {
  it('renders signature rows with first seen, sparkline and state', () => {
    const r = {signature: 'a41f9c3e-1111-2222-3333-444455556666', package: 'Chem', users: 9, count: 57, firstVersion: '1.14.2',
      firstSeen: utcIso(10, 2), lastSeen: utcIso(10, 47), trend: [0, 0, 0, 0, 0, 0, 57], state: 'muted', stateVersion: '1.14.3',
      topError: 'TypeError: Cannot read properties of undefined (reading \'molfile\')'};
    const row = aggregateRow(r, ['package', 'signature'], 'day');
    expect(Object.keys(row)).toEqual(['SIG', 'PACKAGE', 'ERROR', 'USERS', 'COUNT', 'FIRST SEEN', 'LAST', 'TREND (daily)', 'STATE']);
    expect(row.SIG).toBe('a41f9c');
    expect(row['FIRST SEEN']).toMatch(/^1\.14\.2 · \d\d-\d\d$/);
    expect(row['TREND (daily)']).toBe('▁▁▁▁▁▁█');
    expect(row.STATE).toBe('muted → 1.14.3');
    expect(aggregateRow({...r, state: 'muted', stateVersion: null, stateUntil: utcIso(14, 0)}, ['signature'], 'hour').STATE).toBe('muted until 14:00Z');
    expect(aggregateRow({...r, state: undefined}, ['signature'], 'day').STATE).toBe('');
  });

  it('renders rollup rows when signature is not a dimension', () => {
    const row = aggregateRow({connection: 'Snowflake:PROD', signatures: 1, users: 14, count: 61, newInRange: 1, firstVersion: '1.28.2',
      trend: [1, 2], topError: 'SQLTimeoutException: warehouse suspended'}, ['connection'], 'hour');
    expect(Object.keys(row)).toEqual(['CONNECTION', 'SIGNATURES', 'USERS', 'COUNT', 'NEW', 'FIRST SEEN IN', 'LAST', 'TREND (hourly)', 'TOP ERROR']);
    expect(row['TOP ERROR']).toBe('SQLTimeoutException: warehouse suspended');
  });

  it('renders occurrence rows', () => {
    const row = occurrenceRow({time: utcIso(10, 14), user: 'alice', service: 'client', signature: 'a41f9c3e', error: 'boom',
      package: 'Chem', version: '1.14.2', route: 'POST /projects/{id}/save', server: 'datlas-1', requestId: 'mfz3k2a1b9x8y7kq.2'});
    expect(row).toEqual({TIME: '10:14Z', USER: 'alice', SOURCE: 'client', SIG: 'a41f9c', ERROR: 'boom', PACKAGE: 'Chem',
      VERSION: '1.14.2', ROUTE: 'POST /projects/{id}/save', SERVER: 'datlas-1', REQ: '…x8y7kq.2'});
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
  it('writes rows the real libraries read back', () => {
    const bytes = toParquet([{sig: 'a41f9c', users: 9, trend: [0, 57]}, {sig: '7c02e1', users: 6, trend: [3]}]);
    expect(bytes.subarray(0, 4).toString()).toBe('PAR1');
    const arrow = require('apache-arrow');
    const parquet = require('parquet-wasm');
    const table = arrow.tableFromIPC(parquet.readParquet(bytes).intoIPCStream());
    expect(table.toArray().map((r: any) => r.toJSON())).toEqual([{sig: 'a41f9c', users: 9, trend: '[0,57]'}, {sig: '7c02e1', users: 6, trend: '[3]'}]);
  }, 60000);
});

describe('handleErrors', () => {
  beforeEach(() => { vi.useFakeTimers({toFake: ['Date']}); vi.setSystemTime(NOW); });
  afterEach(() => vi.useRealTimers());

  it('asks top for the grouping, trend and limit and prints the table', async () => {
    const {connect, calls} = mockConnect(() => [{signature: '7c02e1aa', package: 'PowerGrid', users: 6, count: 212, firstVersion: '2.3.0',
      firstSeen: utcIso(9, 0, 17), lastSeen: utcIso(10, 0), trend: [3, 4, 3, 5, 4, 6, 5], state: 'open', topError: 'NPE'}]);
    const {out} = await captureOutput(() => handleErrors(connect, 'top', [], {since: '7d', by: 'signature,package', 'min-users': 2, limit: 5}, 'table'));
    expect(calls[0].path).toBe('/errors?since=7d&minUsers=2&by=signature%2Cpackage&trend=day&limit=5');
    expect(out[0].split(/\s{2,}/)).toEqual(['SIG', 'PACKAGE', 'ERROR', 'USERS', 'COUNT', 'FIRST SEEN', 'LAST', 'TREND (daily)', 'STATE']);
    expect(out[2]).toMatch(/^7c02e1\s+PowerGrid\s+NPE\s+6\s+212\s+2\.3\.0 · 09-11\s+10:00Z\s+▅▆▅▇▆█▇\s+open/);
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
      firstSeen: utcIso(10, 2), lastSeen: utcIso(10, 47), occurrences: 57, users: 9,
      groups: [{name: 'Chemists', users: 7}, {name: 'Biology', users: 2}], reports: [4815, 4817],
      alert: {kind: 'error-incident', key: 'a41f9c', status: 'open', openedAt: utcIso(10, 5)},
      change: {type: 'package-published', package: 'Chem', version: '1.14.2', by: 'j.doe', at: utcIso(9, 58)}}));
    const {out} = await captureOutput(() => handleErrors(connect, 'show', ['a41f9c'], {}, 'table'));
    expect(calls[0].path).toBe('/errors/a41f9c?since=24h');
    expect(out).toEqual([
      'signature    a41f9c  TypeError: x',
      'package      Chem    first seen in 1.14.2 at 10:02Z · last seen 10:47Z',
      'occurrences  57      users 9      groups Chemists 7 · Biology 2',
      'reports      #4815 #4817          alert error-incident, open since 10:05Z',
      'change       package-published Chem 1.14.2 by j.doe at 09:58Z',
    ]);
  });

  it('widens the show block for a long package name; says when the breakdowns are sampled', async () => {
    const {connect} = mockConnect(() => ({signature: 'a41f9c3e', error: 'TypeError: x', package: 'UsageAnalysis',
      firstVersion: '2.6.1', firstSeen: utcIso(10, 2), lastSeen: utcIso(10, 47), occurrences: 1234567, users: 1234567890,
      groups: [], reports: [4815, 4816, 4817, 4818, 4819], alert: null, sampled: 5000}));
    const {out} = await captureOutput(() => handleErrors(connect, 'show', ['a41f9c'], {}, 'table'));
    expect(out).toEqual([
      'signature    a41f9c         TypeError: x',
      'package      UsageAnalysis  first seen in 2.6.1 at 10:02Z · last seen 10:47Z',
      'occurrences  1234567        users 1234567890  groups (none)',
      'reports      #4815 #4816 #4817 #4818 #4819',
      'sampled      breakdowns from the latest 5,000 occurrences',
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
    const csv = await captureOutput(() => handleErrors(connect, 'diff', [], {since: '7d', host: ['prod', 'val']}, 'csv'));
    expect(csv.out[0]).toBe('host,signature,package,users,topError');
    expect(csv.out.slice(1)).toEqual(['prod,91bb0c11,core,2,ScatterPlot.render', 'val,3e91aa22,core,6,TableView.close']);
    const quiet = await captureOutput(() => handleErrors(connect, 'diff', [], {since: '7d', host: ['prod', 'val']}, 'quiet'));
    expect(quiet.out).toEqual(['91bb0c11', '3e91aa22']);
    await expect(handleErrors(connect, 'diff', [], {host: ['a', 'b', 'c']}, 'table')).rejects.toThrow(/exactly two/);
  });

  it('lists occurrences from several hosts', async () => {
    const {connect, calls} = mockConnect(() => [{time: utcIso(10, 0), user: 'alice', signature: 'a41f9c3e', error: 'boom'}]);
    const {out} = await captureOutput(() => handleErrors(connect, 'list', [], {host: ['prod', 'val'], user: 'alice'}, 'table'));
    expect(calls.map((c) => c.path)).toEqual(['/errors?since=24h&user=alice&limit=50', '/errors?since=24h&user=alice&limit=50']);
    expect(out[0].split(/\s{2,}/).slice(0, 3)).toEqual(['HOST', 'TIME', 'USER']);
  });

  it('prints occurrences as csv in full and their signatures in quiet', async () => {
    const error = `TableView.close: ${'the view was already detached '.repeat(3).trim()}`;
    const {connect} = mockConnect(() => [{time: utcIso(10, 0, 3), user: 'alice', signature: 'a41f9c3e', error}]);
    const csv = await captureOutput(() => handleErrors(connect, 'list', [], {}, 'csv'));
    expect(csv.out).toEqual(['time,user,signature,error', `${utcIso(10, 0, 3)},alice,a41f9c3e,${error}`]);
    expect((await captureOutput(() => handleErrors(connect, 'list', [], {}, 'quiet'))).out).toEqual(['a41f9c3e']);
  });
});
