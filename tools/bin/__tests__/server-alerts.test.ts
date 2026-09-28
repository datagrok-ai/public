import {describe, it, expect, beforeEach, afterEach, vi} from 'vitest';
import {handleAlerts, muteBody, alertRow, detectionRows, recentServers, deploymentKey} from '../commands/server-alerts';
import {hostList, singleHost} from '../utils/server-client';
import {mockConnect, captureOutput, localIso, apiError} from './obs-helpers';

const ALERT = {id: 'a41f9c00-0000-4000-8000-000000000001', kind: 'error-incident', key: 'a41f9c', severity: 'warning',
  audience: 'platform', status: 'open', openedAt: localIso(10, 5), openedOnServerName: 'datlas-2',
  summary: 'TypeError: Cannot read properties of undefined (reading \'molfile\') in Chem, 9 users in 15 min'};

describe('muteBody', () => {
  const now = new Date(2026, 9, 1, 9, 30);

  it('takes exactly one of --for, --until, --until-version, --forever', () => {
    expect(() => muteBody({reason: 'r'}, now)).toThrow(/exactly one/);
    expect(() => muteBody({reason: 'r', for: '2h', forever: true}, now)).toThrow(/exactly one/);
    expect(() => muteBody({for: '2h'}, now)).toThrow(/--reason/);
    expect(() => muteBody({reason: 'r', 'until-version': ''}, now)).toThrow(/exactly one/);
  });

  it('resolves --for and --until to an absolute time', () => {
    expect(muteBody({reason: 'r', for: '2h'}, now).body.until).toBe(new Date(2026, 9, 1, 11, 30).toISOString());
    expect(muteBody({reason: 'r', until: '14:00'}, now)).toEqual({body: {reason: 'r', until: new Date(2026, 9, 1, 14, 0).toISOString()}, until: 'until 14:00'});
    expect(muteBody({reason: 'r', until: '2026-10-04T06:00'}, now).body.until).toBe(new Date(2026, 9, 4, 6, 0).toISOString());
    expect(muteBody({reason: 'r', until: '2026-10-04T06:00'}, now).until).toBe('until 10-04 06:00');
    expect(() => muteBody({reason: 'r', until: '08:00'}, now)).toThrow(/in the past/);
    expect(() => muteBody({reason: 'r', until: 'tomorrow'}, now)).toThrow(/ISO time/);
  });

  it('passes a version and forever through', () => {
    expect(muteBody({reason: 'fixed', 'until-version': '1.14.3'}, now)).toEqual({body: {reason: 'fixed', untilVersion: '1.14.3'}, until: 'until 1.14.3'});
    expect(muteBody({reason: 'noise', forever: true}, now)).toEqual({body: {reason: 'noise', forever: true}, until: 'forever'});
  });
});

describe('rows', () => {
  it('prints the fixed alert columns, summary cut at 60', () => {
    const row = alertRow(ALERT);
    expect(Object.keys(row)).toEqual(['KIND', 'KEY', 'SEV', 'AUDIENCE', 'STATUS', 'OPENED', 'BY', 'SUMMARY']);
    expect(row).toMatchObject({KIND: 'error-incident', KEY: 'a41f9c', OPENED: '10:05', BY: 'datlas-2'});
    expect(row.SUMMARY.length).toBe(60);
    expect(row.SUMMARY.endsWith('…')).toBe(true);
  });

  it('marks the lease holder', () => {
    const rows = detectionRows({holder: 'S2', servers: [
      {id: 'S1', name: 'datlas-1', host: 'h1', version: '1.28.3', lastSeen: localIso(10, 5), live: false, eligible: true},
      {id: 'S2', name: 'datlas-2', host: 'h2', version: '1.28.3', lastSeen: localIso(10, 8), live: true, eligible: true},
    ]});
    expect(rows.map((r) => [r.SERVER, r.LIVE, r.OWNER])).toEqual([['datlas-1', 'no', ''], ['datlas-2', 'yes', '*']]);
    expect(Object.keys(rows[0])).toEqual(['SERVER', 'HOST NAME', 'VERSION', 'LAST SEEN', 'LIVE', 'ELIGIBLE', 'OWNER']);
  });

  it('keeps live servers and those gone within the hour unless --all', () => {
    const now = Date.parse('2026-10-01T10:00:00Z');
    const lease = {servers: [
      {id: 'L', live: true, lastSeen: '2026-09-01T00:00:00Z'},
      {id: 'S', live: false, stoppedAt: '2026-10-01T09:30:00Z', lastSeen: '2026-10-01T09:30:00Z'},
      {id: 'C', live: false, lastSeen: '2026-10-01T09:10:00Z'},
      {id: 'OLD', live: false, stoppedAt: '2026-09-30T09:00:00Z', lastSeen: '2026-10-01T09:59:00Z'},
      {id: 'GONE', live: false, lastSeen: '2026-09-20T00:00:00Z'},
    ]};
    expect(recentServers(lease, false, now).map((s) => s.id)).toEqual(['L', 'S', 'C']);
    expect(recentServers(lease, true, now).length).toBe(5);
    expect(deploymentKey({servers: [{id: 'b'}, {id: 'a'}]})).toBe(deploymentKey({holder: 'x', servers: [{id: 'a'}, {id: 'b'}]}));
  });
});

describe('hosts', () => {
  it('normalizes --host as minimist leaves it', () => {
    expect(hostList(undefined)).toEqual([]);
    expect(hostList('prod')).toEqual(['prod']);
    expect(hostList(['prod', 'val'])).toEqual(['prod', 'val']);
    expect(singleHost({host: 'prod'}, 'x')).toBe('prod');
    expect(() => singleHost({host: ['a', 'b']}, 'alerts mute')).toThrow(/takes one --host/);
  });
});

describe('handleAlerts', () => {
  beforeEach(() => { vi.useFakeTimers({toFake: ['Date']}); vi.setSystemTime(new Date(2026, 9, 1, 9, 30)); });
  afterEach(() => vi.useRealTimers());

  it('lists with the status filter and merges several hosts under a HOST column', async () => {
    const {connect, calls} = mockConnect((_m, _p, _b, host) => [{...ALERT, key: `${host}-sig`}]);
    const {out} = await captureOutput(() => handleAlerts(connect, 'list', [],
      {status: 'open,acknowledged', host: ['prod', 'val']}, 'table'));
    expect(calls.map((c) => [c.host, c.path])).toEqual([['prod', '/alerts?status=open%2Cacknowledged'], ['val', '/alerts?status=open%2Cacknowledged']]);
    expect(out[0].split(/\s+/).slice(0, 4)).toEqual(['HOST', 'KIND', 'KEY', 'SEV']);
    expect(out.some((l) => l.startsWith('prod') && l.includes('prod-sig'))).toBe(true);
    expect(out.some((l) => l.startsWith('val') && l.includes('val-sig'))).toBe(true);
  });

  it('keeps answering when one host fails, and exits 1', async () => {
    const {connect} = mockConnect((_m, _p, _b, host) => { if (host === 'val') throw apiError(403, 'Forbidden'); return [ALERT]; });
    const {out, err, exitCode} = await captureOutput(() => handleAlerts(connect, 'list', [], {host: ['prod', 'val']}, 'json'));
    expect(err[0]).toMatch(/^val: Forbidden/);
    expect(JSON.parse(out.join('\n'))).toEqual([{host: 'prod', ...ALERT}]);
    expect(exitCode).toBe(1);
  });

  it('looks a kind:key id up among open alerts, then mutes it by id', async () => {
    const {connect, calls} = mockConnect((m) => m === 'GET'
      ? [{id: 'E1', kind: 'connection', key: 'ELN:Prod', status: 'open'}]
      : {id: 'E1', kind: 'connection', key: 'ELN:Prod', status: 'muted'});
    const {out} = await captureOutput(() => handleAlerts(connect, 'mute', ['connection:ELN:Prod'],
      {until: '2026-10-04T06:00', reason: 'monthly ELN maintenance'}, 'table'));
    expect(calls.map((c) => [c.method, c.path])).toEqual([
      ['GET', '/alerts?kind=connection&key=ELN%3AProd&status=open%2Cacknowledged%2Cmuted'],
      ['POST', '/alerts/E1/mute'],
    ]);
    expect(calls[1].body).toEqual({reason: 'monthly ELN maintenance', until: new Date(2026, 9, 4, 6, 0).toISOString()});
    expect(out).toEqual(['muted connection:ELN:Prod until 10-04 06:00 — monthly ELN maintenance']);
  });

  it('refuses a kind:key that matches no open alert, or several', async () => {
    const none = mockConnect(() => []);
    await expect(handleAlerts(none.connect, 'ack', ['report:4820'], {}, 'table')).rejects.toThrow('No open alert report:4820');
    expect(none.calls.length).toBe(1);
    const two = mockConnect(() => [{id: 'A1', kind: 'error-incident', key: 'a41f9c01', status: 'open'},
      {id: 'A2', kind: 'error-incident', key: 'a41f9c02', status: 'muted'}]);
    await expect(handleAlerts(two.connect, 'get', ['error-incident:a41f9c'], {}, 'table'))
      .rejects.toThrow(/Several open alerts match error-incident:a41f9c; pass an id:\n {2}A1 .*\n {2}A2 /);
    expect(two.calls.length).toBe(1);
  });

  it('refuses a mute without a reason before calling the server', async () => {
    const {connect, calls} = mockConnect(() => ({}));
    await expect(handleAlerts(connect, 'mute', ['error-incident:a41f9c'], {'until-version': '1.14.3'}, 'table')).rejects.toThrow(/--reason/);
    expect(calls).toEqual([]);
  });

  it('acknowledges and resolves with an optional reason', async () => {
    const {connect, calls} = mockConnect((_m, path) => ({...ALERT, status: path.endsWith('/ack') ? 'acknowledged' : 'resolved'}));
    const ack = await captureOutput(() => handleAlerts(connect, 'ack', ['a41f9c00'], {}, 'table'));
    const res = await captureOutput(() => handleAlerts(connect, 'resolve', ['a41f9c00'], {reason: 'fixed'}, 'table'));
    expect(calls.map((c) => [c.path, c.body])).toEqual([['/alerts/a41f9c00/ack', {reason: undefined}], ['/alerts/a41f9c00/resolve', {reason: 'fixed'}]]);
    expect(ack.out).toEqual(['acknowledged error-incident:a41f9c']);
    expect(res.out).toEqual(['resolved error-incident:a41f9c — fixed']);
  });

  it('refuses several hosts for a transition', async () => {
    const {connect} = mockConnect(() => ({}));
    await expect(handleAlerts(connect, 'ack', ['x'], {host: ['a', 'b']}, 'table')).rejects.toThrow(/takes one --host/);
  });

  it('shows the detection lease per server', async () => {
    const {connect, calls} = mockConnect(() => ({holder: 'S1', holderName: 'datlas-1', epoch: 3,
      servers: [{id: 'S1', name: 'datlas-1', host: 'h1', version: '1.28.3', lastSeen: localIso(9, 29), live: true, eligible: true}]}));
    const {out} = await captureOutput(() => handleAlerts(connect, 'detection', [], {}, 'table'));
    expect(calls[0].path).toBe('/alerts/detection');
    expect(out[2].trimEnd()).toMatch(/^datlas-1\s+h1\s+1\.28\.3\s+09:29\s+yes\s+yes\s+\*$/);
  });

  it('names rows by --host alias and prints a shared database once', async () => {
    const server = (id: string, name: string) => ({id, name, host: 'pc-alex2', version: '1.28.3', lastSeen: localIso(9, 29), live: true, eligible: true});
    const shared = {holder: 'S1', servers: [server('S1', 'datlas-1'), server('S2', 'datlas-2')]};
    const {connect} = mockConnect((_m, _p, _b, host) => host === 'c' ? {holder: 'S9', servers: [server('S9', 'datlas-9')]} : shared);
    const {out} = await captureOutput(() => handleAlerts(connect, 'detection', [], {host: ['a', 'b', 'c']}, 'table'));
    expect(out[0].trimEnd().split(/\s{2,}/)).toEqual(['HOST', 'SERVER', 'HOST NAME', 'VERSION', 'LAST SEEN', 'LIVE', 'ELIGIBLE', 'OWNER']);
    expect(out.slice(2).map((l) => l.split(/\s{2,}/).slice(0, 2))).toEqual([['a, b', 'datlas-1'], ['a, b', 'datlas-2'], ['c', 'datlas-9']]);
    const json = await captureOutput(() => handleAlerts(connect, 'detection', [], {host: ['a', 'b', 'c']}, 'json'));
    expect(JSON.parse(json.out.join('\n')).map((d: any) => d.hosts)).toEqual([['a', 'b'], ['c']]);
  });

  it('answers an unknown verb with the usage', async () => {
    const {connect} = mockConnect(() => ({}));
    const {err, result} = await captureOutput(() => handleAlerts(connect, 'frobnicate', [], {}, 'table'));
    expect(result).toBe(false);
    expect(err.join('\n')).toMatch(/Usage: grok s alerts/);
  });
});
