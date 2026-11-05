import {describe, it, expect, afterEach, vi} from 'vitest';
import {NodeApiClient, NodeDapi, NodeFuncsDataSource} from '../utils/node-dapi';
import {logSpecOf} from '../utils/server-client';
import {handleObserve, printRequestLog} from '../commands/server-observe';

const URL = 'http://stand/api';
const LOG_ID = '0b8a3a2e-6c1d-4f0e-9a51-3f1c2d4e5a6b';

function stub(body: any = {id: 'x'}) {
  const spy = vi.fn(async (_url: string, _opts: any) =>
    new Response(JSON.stringify(body), {headers: {'content-type': 'application/json'}}));
  vi.stubGlobal('fetch', spy);
  return spy;
}

// JSON lines, oldest first, as the server gathers them.
const rows = [
  {time: '2026-10-19T20:18:02.921Z', level: 'info', source: 'Server', message: 'Refreshing connection...', params: {}},
  {time: '2026-10-19T20:18:02.939Z', level: 'debug', flag: 'db', source: 'server', message: 'Executed select 1', params: {}},
];

function fakeDapi(answer: any) {
  const calls: string[] = [];
  const client: any = {
    async request(method: string, path: string) { calls.push(`${method} ${path}`); return answer; },
    get(path: string) { return this.request('GET', path); },
  };
  return {dapi: new NodeDapi(client), calls};
}

afterEach(() => vi.unstubAllGlobals());

describe('--log', () => {
  it('sends the logId and the spec on every request', async () => {
    const spy = stub();
    const client = new NodeApiClient(URL, 'token');
    client.logId = LOG_ID;
    client.logSpec = 'verbose,db';
    await client.get('/projects/x');
    await client.getBytes('/tables/x/data').catch(() => null);
    for (const call of spy.mock.calls) {
      expect(call[1].headers['datlas-log-id']).toBe(LOG_ID);
      expect(call[1].headers['datlas-log']).toBe('verbose,db');
    }
  });

  it('sends neither without --log', async () => {
    const spy = stub();
    await new NodeApiClient(URL, 'token').get('/projects/x');
    expect(spy.mock.calls[0][1].headers['datlas-log-id']).toBeUndefined();
    expect(spy.mock.calls[0][1].headers['datlas-log']).toBeUndefined();
  });

  it('is verbose when given without a value', () => {
    expect(logSpecOf({log: true})).toBe('verbose');
    expect(logSpecOf({log: 'debug,db'})).toBe('debug,db');
    expect(logSpecOf({})).toBeUndefined();
  });

  it('prints the collected log oldest first to stderr', async () => {
    const {dapi, calls} = fakeDapi(rows);
    const out: string[] = [];
    const write = vi.spyOn(process.stderr, 'write').mockImplementation((s: any) => (out.push(String(s)), true));
    try {
      await printRequestLog(dapi, LOG_ID, 'table');
    }
    finally {
      write.mockRestore();
    }
    expect(calls).toEqual([`GET /log/requests/${LOG_ID}`]);
    const text = out.join('');
    expect(text).toContain('2 records');
    expect(text.indexOf('Refreshing connection')).toBeLessThan(text.indexOf('Executed select 1'));
    expect(text).toContain('[db server] Executed select 1');
  });

  it('says so when the server cannot return the log', async () => {
    const client: any = {async request() { throw new Error('Not found'); }, get() { return this.request(); }};
    const out: string[] = [];
    const write = vi.spyOn(process.stderr, 'write').mockImplementation((s: any) => (out.push(String(s)), true));
    try {
      await printRequestLog(new NodeDapi(client), LOG_ID, 'table');
    }
    finally {
      write.mockRestore();
    }
    expect(out.join('')).toContain(`could not read the log of ${LOG_ID}: Not found`);
  });
});

describe('grok s o logs request', () => {
  it('reads one logId as the server returns it', async () => {
    const {dapi} = fakeDapi(rows);
    const lines: string[] = [];
    const log = vi.spyOn(console, 'log').mockImplementation((s) => { lines.push(String(s)); });
    try {
      await handleObserve(dapi, 'logs', 'request', [LOG_ID], {}, 'json');
    }
    finally {
      log.mockRestore();
    }
    const printed = JSON.parse(lines.join('\n'));
    expect(printed.map((r: any) => r.message)).toEqual(['Refreshing connection...', 'Executed select 1']);
    expect(printed[1].flag).toBe('db');
  });
});

describe('functions run --debug', () => {
  it('asks for a debug call', async () => {
    const spy = stub(42);
    expect(await new NodeFuncsDataSource(new NodeApiClient(URL, 'token')).run('Pkg:query', {}, true)).toBe(42);
    expect(spy.mock.calls[0][0]).toBe(`${URL}/public/v1/functions/Pkg.query/call?debug=true`);
  });
});
