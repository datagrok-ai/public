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
        return method === 'GET' ? {'#key': 'logger', settings: {userGroupSettings: {g0: {debugFlags: []}}}} : null;
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

describe('timeline', () => {
  async function timelinePath(argv: any): Promise<URLSearchParams> {
    const paths: string[] = [];
    const client: any = {
      async request(method: string, path: string) {
        paths.push(path);
        return [];
      },
      get(path: string) { return this.request('GET', path); },
    };
    const log = vi.spyOn(console, 'log').mockImplementation(() => {});
    try {
      await handleObserve(new NodeDapi(client), 'timeline', undefined, [], argv, 'json');
    }
    finally {
      log.mockRestore();
    }
    expect(paths[0].startsWith('/log/timeline?')).toBe(true);
    return new URLSearchParams(paths[0].split('?')[1]);
  }

  it('a session without --from reads its last 15 minutes', async () => {
    const q = await timelinePath({session: 's1'});
    expect(q.get('session')).toBe('s1');
    expect(Date.parse(q.get('to')!) - Date.parse(q.get('from')!)).toBe(15 * 60000);
  });

  it('a session with --from leaves --to to the server unless given', async () => {
    const q = await timelinePath({session: 's1', from: '2026-10-09T10:00', limit: 5});
    expect([q.get('from'), q.has('to'), q.get('limit')]).toEqual(['2026-10-09T10:00:00.000Z', false, '5']);
  });

  it('needs a session', async () => {
    await expect(handleObserve(new NodeDapi({} as any), 'timeline', undefined, [], {}, 'json')).rejects.toThrow('--session');
  });
});

describe('logs', () => {
  async function run(verb: string, rest: string[], argv: any, answers: Record<string, any> = {}) {
    const paths: string[] = [];
    const client: any = {
      async request(method: string, path: string) {
        paths.push(path);
        return answers[path.split('?')[0]] ?? [];
      },
      get(path: string) { return this.request('GET', path); },
    };
    const lines: string[] = [];
    const log = vi.spyOn(console, 'log').mockImplementation((s) => { lines.push(String(s)); });
    try {
      await handleObserve(new NodeDapi(client), 'logs', verb, rest, argv, 'json');
    }
    finally {
      log.mockRestore();
    }
    return {queries: paths.map((p) => [p.split('?')[0], new URLSearchParams(p.split('?')[1])] as [string, URLSearchParams]),
      printed: lines.join('\n')};
  }

  it('cloud reads the first group of the instance for the last hour, as JSON rows', async () => {
    const {queries} = await run('cloud', [], {filter: 'ERROR'}, {'/log/cloud/groups': ['/datagrok/dev', '/eks/pods/x']});
    expect(queries[0][0]).toBe('/log/cloud/groups');
    const [path, q] = queries[1];
    expect(path).toBe('/log/cloud/events');
    expect([q.get('group'), q.get('filter'), q.get('format'), q.has('end')]).toEqual(['/datagrok/dev', 'ERROR', 'json', false]);
    expect(Date.now() - Date.parse(q.get('start')!)).toBeGreaterThanOrEqual(3600000 - 1000);
  });

  it('cloud takes a group and an absolute window', async () => {
    const {queries} = await run('cloud', [], {group: '/g', from: '2026-10-09T10:00', to: '2026-10-09T11:00', limit: 5});
    expect(queries.length).toBe(1);
    const q = queries[0][1];
    expect([q.get('group'), q.get('start'), q.get('end'), q.get('limit')])
      .toEqual(['/g', '2026-10-09T10:00:00.000Z', '2026-10-09T11:00:00.000Z', '5']);
    await expect(run('cloud', [], {group: '/g', since: '1h', from: '-2h'})).rejects.toThrow('--since');
  });

  it('archive list keeps the objects modified within --since', async () => {
    const recent = new Date(Date.now() - 3600000).toISOString();
    const {queries, printed} = await run('archive', ['list'], {prefix: 'cw/', since: '1d'}, {'/log/archive/objects': [
      {key: 'cw/new.gz', modified: recent, size: 10}, {key: 'cw/old.gz', modified: '2020-01-01T00:00:00.000Z', size: 10}]});
    expect([queries[0][1].get('prefix'), queries[0][1].get('format')]).toEqual(['cw/', 'json']);
    expect(JSON.parse(printed).map((o: any) => o.key)).toEqual(['cw/new.gz']);
  });

  it('archive read decodes one key', async () => {
    const {queries} = await run('archive', ['read', 'cw/tenant-a/app/2026/10/09/x.gz'], {connection: 'c1'});
    expect(queries[0][0]).toBe('/log/archive/events');
    expect([queries[0][1].get('key'), queries[0][1].get('connection'), queries[0][1].get('format')])
      .toEqual(['cw/tenant-a/app/2026/10/09/x.gz', 'c1', 'json']);
    await expect(run('archive', ['read'], {})).rejects.toThrow('Usage');
  });
});
