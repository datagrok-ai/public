import {describe, it, expect, vi} from 'vitest';
import {NodeDapi} from '../utils/node-dapi';
import {handleObserve, parseDuration, parseTime, setPath} from '../commands/server-observe';

describe('parseDuration', () => {
  it('reads m as minutes, and h, d, w', () => {
    expect(parseDuration('30m', '--for')).toBe(30 * 60000);
    expect(parseDuration('2h', '--for')).toBe(2 * 3600000);
    expect(parseDuration('-7d', '--since')).toBe(7 * 86400000);
    expect(parseDuration('1w', '--for')).toBe(604800000);
  });

  it('refuses anything else', () => {
    for (const bad of ['', '0m', '5', '2y', 'soon'])
      expect(() => parseDuration(bad, '--for')).toThrow(/--for expects a duration/);
  });
});

describe('parseTime', () => {
  const now = new Date('2026-10-08T12:00:00Z');

  it('reads ISO as UTC unless it has an offset, and relative times', () => {
    expect(parseTime('2026-10-04T06:00', '--from', now).toISOString()).toBe('2026-10-04T06:00:00.000Z');
    expect(parseTime('2026-10-04T06:00+02:00', '--from', now).toISOString()).toBe('2026-10-04T04:00:00.000Z');
    expect(parseTime('2026-10-04', '--from', now).toISOString()).toBe('2026-10-04T00:00:00.000Z');
    expect(parseTime('-2h', '--from', now).toISOString()).toBe('2026-10-08T10:00:00.000Z');
  });

  it('refuses anything else', () => {
    for (const bad of ['yesterday', '10:00', '2026-13-45'])
      expect(() => parseTime(bad, '--to', now)).toThrow(/--to expects an ISO time/);
  });
});

describe('setPath', () => {
  it('sets a nested value from JSON, else a string, and returns the top key', () => {
    const target: any = {exportSettings: [{types: 'error'}]};
    expect(setPath(target, 'exportSettings.0.types=error,audit')).toBe('exportSettings');
    expect(setPath(target, 'debugFlags=["db"]')).toBe('debugFlags');
    setPath(target, 'a.b=3');
    expect(target).toEqual({exportSettings: [{types: 'error,audit'}], debugFlags: ['db'], a: {b: 3}});
  });

  it('refuses an assignment without a path', () => {
    expect(() => setPath({}, '=1')).toThrow(/<path>=<json>/);
  });
});

describe('logger set', () => {
  it('adds an expiring entry to userGroupSettings, keeping the others', async () => {
    const calls: {method: string; path: string; body: any}[] = [];
    const client: any = {
      async request(method: string, path: string, body?: any) {
        calls.push({method, path, body});
        if (path.startsWith('/public/v1/groups/lookup'))
          return [{id: 'g1', name: 'Chemists'}];
        return method === 'GET' ? {userGroupSettings: {g0: {debugFlags: []}}} : null;
      },
      get(path: string) { return this.request('GET', path); },
    };
    const log = vi.spyOn(console, 'log').mockImplementation(() => {});
    try {
      await handleObserve(new NodeDapi(client), 'logger', 'set', [],
        {group: 'Chemists', set: 'debugFlags=["db"]', for: '30m', reason: 'ticket 1'}, 'table');
    }
    finally {
      log.mockRestore();
    }
    const saved = calls.find((c) => c.method === 'POST')!;
    expect(saved.path).toBe('/admin/plugins/logger/settings');
    expect(Object.keys(saved.body.userGroupSettings)).toEqual(['g0', 'g1']);
    expect(saved.body.userGroupSettings.g1.debugFlags).toEqual(['db']);
    expect(saved.body.userGroupSettings.g1.reason).toBe('ticket 1');
    expect(Date.parse(saved.body.userGroupSettings.g1.expiresAt) - Date.now()).toBeGreaterThan(29 * 60000);
  });
});
