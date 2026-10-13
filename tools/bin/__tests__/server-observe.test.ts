import {describe, it, expect, vi} from 'vitest';
import {NodeDapi} from '../utils/node-dapi';
import {OutputFormat} from '../utils/server-output';
import {ALL_USERS} from '../utils/migrate/walker';
import {handleObserve, parseDuration, parseTime, setPath} from '../commands/server-observe';

type Call = {method: string; path: string; body: any};

/** Runs one observe command against a fake client that answers with [answer]; returns its calls and what it printed. */
async function observe(command: string, verb: string | undefined, rest: string[], argv: any,
                       answer: (method: string, path: string) => any = () => [], output: OutputFormat = 'json') {
  const calls: Call[] = [];
  const client: any = {
    async request(method: string, path: string, body?: any) {
      calls.push({method, path, body});
      return answer(method, path);
    },
    get(path: string) { return this.request('GET', path); },
  };
  const lines: string[] = [];
  const log = vi.spyOn(console, 'log').mockImplementation((s) => { lines.push(String(s)); });
  try {
    await handleObserve(new NodeDapi(client), command, verb, rest, argv, output);
  }
  finally {
    log.mockRestore();
  }
  return {calls, printed: lines.join('\n'), saved: calls.find((c) => c.method === 'POST')?.body?.userGroupSettings};
}

/** The server answers the settings as text, not application/json. */
const loggerSettings = (userGroupSettings: any) => (method: string, path: string) => {
  if (path.startsWith('/public/v1/groups/lookup'))
    return path.includes('alice') ? [{id: 'p1', name: 'alice', personal: true}] : [{id: 'g1', name: 'Chemists'}];
  if (path.startsWith('/public/v1/groups/'))
    return {id: 'p1', friendlyName: 'Alice'};
  return method === 'GET' ? JSON.stringify({'#key': 'logger', settings: {userGroupSettings}}) : null;
};

describe('parseDuration', () => {
  it('reads m as minutes, and h, d, w', () => {
    expect(parseDuration('30m', '--for')).toBe(30 * 60000);
    expect(parseDuration('2h', '--for')).toBe(2 * 3600000);
    expect(parseDuration('7d', '--since')).toBe(7 * 86400000);
    expect(parseDuration('1w', '--for')).toBe(604800000);
  });

  it('refuses anything else, negative values included', () => {
    for (const bad of ['', '0m', '5', '2y', 'soon', '-7d'])
      expect(() => parseDuration(bad, '--for')).toThrow(/--for expects a duration/);
  });
});

describe('parseTime', () => {
  const now = new Date('2026-10-08T12:00:00Z');

  it('reads ISO as UTC unless it has an offset, and relative times', () => {
    expect(parseTime('2026-10-04T06:00', '--from', now).toISOString()).toBe('2026-10-04T06:00:00.000Z');
    expect(parseTime('2026-10-04T06:00+02:00', '--from', now).toISOString()).toBe('2026-10-04T04:00:00.000Z');
    expect(parseTime('2026-10-04T06:00+02', '--from', now).toISOString()).toBe('2026-10-04T04:00:00.000Z');
    expect(parseTime('2026-10-04T06:00-0130', '--from', now).toISOString()).toBe('2026-10-04T07:30:00.000Z');
    expect(parseTime('2026-10-04', '--from', now).toISOString()).toBe('2026-10-04T00:00:00.000Z');
    expect(parseTime('-2h', '--from', now).toISOString()).toBe('2026-10-08T10:00:00.000Z');
  });

  it('refuses anything else', () => {
    for (const bad of ['yesterday', '10:00', '2026-13-45', '--2h'])
      expect(() => parseTime(bad, '--to', now)).toThrow(/--to expects/);
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

describe('problems', () => {
  it('reads an id without a colon by its route, encoded once', async () => {
    const {calls} = await observe('problems', 'get', ['a b'], {}, () => ({id: 'a b'}));
    expect(calls.map((c) => c.path)).toEqual(['/problems/a%20b']);
  });

  it('looks kind:key up by the list filters', async () => {
    const {calls} = await observe('problems', 'history', ['health:Jupyter@x'], {limit: 3},
      (_, path) => path.startsWith('/problems?') ? [{id: 'p 1'}] : []);
    expect(calls.map((c) => c.path)).toEqual(['/problems?kind=health&key=Jupyter%40x&limit=1', '/problems/p%201/history?limit=3']);
  });

  it('refuses an empty kind or key', async () => {
    for (const bad of [':x', 'health:', ':', ''])
      await expect(observe('problems', 'get', [bad], {})).rejects.toThrow('Usage');
  });

  it('get prints one field per line in a table, details by path', async () => {
    const {printed} = await observe('problems', 'get', ['p1'], {},
      () => ({id: 'p1', details: {error: 'x'.repeat(80), nested: {n: 1}}}), 'table');
    expect(printed).toContain('x'.repeat(80));
    expect(printed).toMatch(/details\.nested\.n\s+1/);
  });
});

describe('errors', () => {
  const top = [{message: 'a', source: 'client'}, {message: 'b', source: 'server'}, {message: 'c', source: 'server'}];

  it('top keeps the global top --limit of the per-source answer', async () => {
    const {printed} = await observe('errors', 'top', [], {limit: 2}, () => ({errors: {top}}));
    expect(JSON.parse(printed).top.map((e: any) => e.message)).toEqual(['a', 'b']);
  });

  it('get reads a problem\'s signature and searches both kinds of signature', async () => {
    const {calls, printed} = await observe('errors', 'get', ['error-incident:sig'], {limit: 5}, (_, path) =>
      path.startsWith('/problems?') ? [{id: 'p1'}]
        : path === '/problems/p1' ? {id: 'p1', details: {signatureId: 's-1'}}
          : path.includes('/users/') ? {login: 'system'}
            : [{id: 'e1', eventTime: 't', description: 'boom', parameters: [
              {parameter: {name: 'user'}, value: 'u1'}, {parameter: {name: 'function'}, value: 'Test'}]}]);
    const search = decodeURIComponent(calls.find((c) => c.path.startsWith('/log?'))!.path);
    expect(search).toContain('(eventType.id = "s-1" or eventType.errorStackTraceHash = "s-1") and eventTime > "');
    expect(search).toContain('limit=5');
    expect(JSON.parse(printed)[0]).toMatchObject({time: 't', user: 'system', function: 'Test', message: 'boom'});
  });

  it('get refuses a problem that is not about an error', async () => {
    await expect(observe('errors', 'get', ['p1'], {}, () => ({id: 'p1', details: {}}))).rejects.toThrow('no signatureId');
  });
});

describe('logs container', () => {
  const lines = ['2026-10-13T17:21:14.1+02:00 Kernel started', '  at frame', '2026-10-13T17:25:00.5+02:00 Kernel shutdown',
    'no stamp at the start'].join('\n');

  it('cuts the window in UTC across the container\'s offset, keeping continuation lines', async () => {
    const {calls, printed} = await observe('logs', 'container', ['c1'], {from: '2026-10-13T15:21:00', to: '2026-10-13T15:22:00'},
      (_, path) => path.startsWith('/docker/containers?') ? [{id: 'id-1', name: 'c1'}] : lines, 'table');
    expect(calls.map((c) => c.path)).toEqual(['/docker/containers?text=c1&limit=100', '/docker/containers/id-1/logs?limit=20000']);
    expect(printed).toBe('2026-10-13T17:21:14.1+02:00 Kernel started\n  at frame');
  });

  it('--grep filters the lines', async () => {
    const {printed} = await observe('logs', 'container', ['8b6f2a1e-0000-4000-8000-000000000001'], {grep: 'shutdown'},
      () => lines, 'table');
    expect(printed).toBe('2026-10-13T17:25:00.5+02:00 Kernel shutdown');
  });

  it('names similar containers when the name does not match', async () => {
    await expect(observe('logs', 'container', ['py'], {}, () => [{id: 'i', name: 'scripting-python_1'}]))
      .rejects.toThrow("No container 'py'; similar: scripting-python_1");
  });
});

describe('logger set', () => {
  const set = async (argv: any, userGroupSettings: any) =>
    (await observe('logger', 'set', [], {group: 'Chemists', ...argv}, loggerSettings(userGroupSettings), 'table')).saved;

  it('starts a new entry from All users, without its reason', async () => {
    const saved = await set({set: 'debugFlags=["db"]'}, {[ALL_USERS]: {saveLevels: ['error', 'audit'], reason: 'base'}});
    expect(saved.g1).toEqual({saveLevels: ['error', 'audit'], debugFlags: ['db']});
  });

  it('--remove deletes the entry, keeping the others', async () => {
    expect(await set({remove: true}, {g1: {debugFlags: ['db']}, g0: {}})).toEqual({g0: {}});
    await expect(set({remove: true}, {})).rejects.toThrow('no logger entry');
  });

  it('refuses an expiry over an entry that does not expire', async () => {
    await expect(set({set: 'debugFlags=["db"]', for: '1h'}, {g1: {debugFlags: []}})).rejects.toThrow('does not expire');
  });

  it('refuses --until in the past', async () => {
    await expect(set({set: 'debugFlags=["db"]', until: '2020-01-01T00:00'}, {})).rejects.toThrow('in the past');
  });

  it('adds an expiring entry to userGroupSettings, keeping the others', async () => {
    const {calls} = await observe('logger', 'set', [], {group: 'Chemists', set: 'debugFlags=["db"]', for: '30m',
      reason: 'ticket 1'}, (method, path) => path.startsWith('/public/v1/groups/lookup') ? [{id: 'g1', name: 'Chemists'}]
      : method === 'GET' ? {'#key': 'logger', settings: {userGroupSettings: {g0: {debugFlags: []}}}} : null, 'table');
    const saved = calls.find((c) => c.method === 'POST')!;
    expect(saved.path).toBe('/admin/plugins/logger/settings');
    expect(Object.keys(saved.body.userGroupSettings)).toEqual(['g0', 'g1']);
    expect(saved.body.userGroupSettings.g1.debugFlags).toEqual(['db']);
    expect(saved.body.userGroupSettings.g1.reason).toBe('ticket 1');
    expect(Date.parse(saved.body.userGroupSettings.g1.expiresAt) - Date.now()).toBeGreaterThan(29 * 60000);
  });
});

describe('logger history', () => {
  it('keeps a diff the server truncated as text', async () => {
    const event = (diff: string) => ({eventTime: '2026-10-09T10:00:00Z',
      parameters: [{parameter: {name: 'user'}, value: 'admin'}, {parameter: {name: 'diff'}, value: diff}]});
    const {printed} = await observe('logger', 'history', [], {},
      () => [event('[{"path": "saveHttpRequests"}]'), event('[{"path": "saveHt…')]);
    expect(JSON.parse(printed).map((r: any) => r.diff)).toEqual([[{path: 'saveHttpRequests'}], '[{"path": "saveHt…']);
  });
});

describe('logs capture', () => {
  const run = (verb: string | undefined, rest: string[], argv: any, userGroupSettings: any) =>
    observe('logs', 'capture', verb === undefined ? rest : [verb, ...rest], argv, loggerSettings(userGroupSettings));

  it('adds an expiring entry seeded from All Users that saves requests, keeping the others', async () => {
    const {saved} = await run('alice', [], {for: '2h'},
      {[ALL_USERS]: {saveLevels: ['error'], postLevels: ['error'], reason: 'baseline'}, g0: {debugFlags: ['db']}});
    expect(Object.keys(saved)).toEqual([ALL_USERS, 'g0', 'p1']);
    expect(saved.p1.saveHttpRequests).toBe(true);
    expect(saved.p1.saveLevels).toEqual(['error', 'usage']);
    expect(saved.p1.postLevels).toEqual(['error', 'usage']);
    expect(saved.p1.reason).toBeUndefined();
    expect(Date.parse(saved.p1.expiresAt) - Date.now()).toBeGreaterThan(119 * 60000);
  });

  it('waits out a save refused while the previous one is applied', async () => {
    vi.useFakeTimers();
    try {
      let posts = 0;
      const answer = loggerSettings({});
      const done = observe('logs', 'capture', ['alice'], {}, (method, path) => {
        if (method === 'POST' && posts++ === 0)
          throw Object.assign(new Error('wait'), {apiError: {error: 'wait', errorCode: 429}});
        return answer(method, path);
      });
      await vi.advanceTimersByTimeAsync(5000);
      const {calls, saved} = await done;
      expect(calls.filter((c) => c.method === 'POST').length).toBe(2);
      expect(saved.p1.saveHttpRequests).toBe(true);
    }
    finally {
      vi.useRealTimers();
    }
  });

  it('falls back to a group and defaults to 30 minutes', async () => {
    const {saved} = await run('Chemists', [], {reason: 'ticket 1'}, {});
    expect(saved.g1.reason).toBe('ticket 1');
    const left = Date.parse(saved.g1.expiresAt) - Date.now();
    expect(left).toBeGreaterThan(29 * 60000);
    expect(left).toBeLessThanOrEqual(30 * 60000);
  });

  it('keeps the later of the entry\'s expiry and --for', async () => {
    const later = '2099-01-01T00:00:00.000Z';
    expect((await run('alice', [], {for: '1h'}, {p1: {debugFlags: ['db'], expiresAt: later}})).saved.p1.expiresAt).toBe(later);
    const {saved} = await run('alice', [], {for: '1h'}, {p1: {expiresAt: new Date(Date.now() + 60000).toISOString()}});
    expect(Date.parse(saved.p1.expiresAt) - Date.now()).toBeGreaterThan(59 * 60000);
  });

  it('leaves an entry that does not expire alone', async () => {
    await expect(run('alice', [], {}, {p1: {debugFlags: ['db']}})).rejects.toThrow('does not expire');
  });

  it('stop removes the entry, keeping the others', async () => {
    const {saved} = await run('stop', ['alice'], {}, {[ALL_USERS]: {debugFlags: []},
      p1: {saveHttpRequests: true, saveLevels: ['usage'], debugFlags: [], expiresAt: '2099-01-01T00:00:00Z'}, g0: {debugFlags: ['db']}});
    expect(Object.keys(saved)).toEqual([ALL_USERS, 'g0']);
    await expect(run('stop', ['alice'], {}, {})).rejects.toThrow('No capture');
  });

  it('stop leaves an entry that does not expire, or holds other settings, to logger set', async () => {
    await expect(run('stop', ['alice'], {}, {p1: {saveHttpRequests: true}})).rejects.toThrow('does not expire');
    await expect(run('stop', ['alice'], {}, {p1: {saveHttpRequests: true, debugFlags: ['db'], expiresAt: '2099-01-01T00:00:00Z'}}))
      .rejects.toThrow('also holds debugFlags; remove it with grok s o logger set --user alice --remove');
  });

  it('lists the active captures', async () => {
    const {saved, printed} = await run(undefined, [], {}, {p1: {saveHttpRequests: true, expiresAt: '2099-01-01T00:00:00Z'},
      g0: {debugFlags: ['db'], expiresAt: '2099-01-01T00:00:00Z'}});
    expect(saved).toBeUndefined();
    expect(JSON.parse(printed)).toEqual([{id: 'p1', name: 'Alice', expiresAt: '2099-01-01T00:00:00Z', reason: ''}]);
  });
});

const query = (path: string): [string, URLSearchParams] => [path.split('?')[0], new URLSearchParams(path.split('?')[1])];

describe('timeline', () => {
  async function timelineQuery(argv: any): Promise<URLSearchParams> {
    const {calls} = await observe('timeline', undefined, [], argv);
    const [path, q] = query(calls[0].path);
    expect(path).toBe('/log/timeline');
    return q;
  }

  it('a session without --from reads its last 15 minutes', async () => {
    const q = await timelineQuery({session: 's1'});
    expect(q.get('session')).toBe('s1');
    expect(Date.parse(q.get('to')!) - Date.parse(q.get('from')!)).toBe(15 * 60000);
  });

  it('a session with --from leaves --to to the server unless given', async () => {
    const q = await timelineQuery({session: 's1', from: '2026-10-09T10:00', limit: 5});
    expect([q.get('from'), q.has('to'), q.get('limit')]).toEqual(['2026-10-09T10:00:00.000Z', false, '5']);
  });

  it('needs a session, and --from for --to', async () => {
    await expect(observe('timeline', undefined, [], {})).rejects.toThrow('--session');
    await expect(observe('timeline', undefined, [], {session: 's1', to: '2026-10-09T10:00'})).rejects.toThrow('--to needs --from');
  });
});

describe('errors top', () => {
  it('refuses --date with --from', async () => {
    await expect(observe('errors', 'top', [], {date: 'today', from: '-2h'})).rejects.toThrow('--date is not used with --from');
  });
});

describe('logs', () => {
  async function run(verb: string, rest: string[], argv: any, answers: Record<string, any> = {}) {
    const {calls, printed} = await observe('logs', verb, rest, argv, (_, path) => answers[path.split('?')[0]] ?? []);
    return {queries: calls.map((c) => query(c.path)), printed};
  }

  it('cloud reads the first group of the instance for the last hour, as JSON rows', async () => {
    const {queries} = await run('cloud', [], {filter: 'ERROR'}, {'/log/cloud/groups': ['/datagrok/dev', '/eks/pods/x']});
    expect(queries[0][0]).toBe('/log/cloud/groups');
    const [path, q] = queries[1];
    expect(path).toBe('/log/cloud/events');
    expect([q.get('group'), q.get('filter'), q.get('format'), q.has('end')]).toEqual(['/datagrok/dev', 'ERROR', 'json', false]);
    const age = Date.now() - Date.parse(q.get('start')!);
    expect(age).toBeGreaterThanOrEqual(3600000 - 1000);
    expect(age).toBeLessThan(3600000 + 60000);
  });

  it('cloud takes a group and an absolute window', async () => {
    const {queries} = await run('cloud', [], {group: '/g', from: '2026-10-09T10:00', to: '2026-10-09T11:00', limit: 5});
    expect(queries.length).toBe(1);
    const q = queries[0][1];
    expect([q.get('group'), q.get('start'), q.get('end'), q.get('limit')])
      .toEqual(['/g', '2026-10-09T10:00:00.000Z', '2026-10-09T11:00:00.000Z', '5']);
    await expect(run('cloud', [], {group: '/g', since: '1h', from: '-2h'})).rejects.toThrow('--since');
    await expect(run('cloud', [], {group: '/g', to: '2026-10-09T11:00'})).rejects.toThrow('--to needs --from');
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

describe('an unknown command', () => {
  it('fails with the usage', async () => {
    await expect(observe('problem', 'list', [], {})).rejects.toThrow(/Unknown command 'problem list'\nUsage: grok s observe/);
  });
});
