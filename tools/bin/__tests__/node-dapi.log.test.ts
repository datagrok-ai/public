import {describe, it, expect, afterEach, vi} from 'vitest';
import {NodeApiClient, NodeFuncsDataSource, CapturedLog, LOG_HEADER} from '../utils/node-dapi';
import {logSpecOf} from '../utils/server-client';

const URL = 'http://stand/api';

function envelope(body: string | Buffer, contentType: string | null, extra: Record<string, string> = {}, status = 200) {
  const binary = Buffer.isBuffer(body);
  return new Response(JSON.stringify({
    contentType, encoding: binary ? 'base64' : null, body: binary ? body.toString('base64') : body,
    sessionId: 's1', dropped: 0, log: [{time: '2026-10-19T10:00:00.123Z', level: 'debug', flag: 'db', message: 'select 1'}],
  }), {status, headers: {'content-type': 'application/json', [LOG_HEADER]: 'envelope', ...extra}});
}

function stub(response: () => Response) {
  const spy = vi.fn(async (_url: string, _opts: any) => response());
  vi.stubGlobal('fetch', spy);
  return spy;
}

function client(logSpec?: string) {
  const c = new NodeApiClient(URL, 'token');
  const logs: {request: string, log: CapturedLog}[] = [];
  c.logSpec = logSpec;
  c.onLog = (request, log) => logs.push({request, log});
  return {c, logs};
}

afterEach(() => vi.unstubAllGlobals());

describe('--log', () => {
  it('sends the spec on every request and hands back the unwrapped body', async () => {
    const spy = stub(() => envelope('{"id":"x"}', 'application/json'));
    const {c, logs} = client('verbose,db');
    expect(await c.get('/projects/x')).toEqual({id: 'x'});
    expect(spy.mock.calls[0][1].headers[LOG_HEADER]).toBe('verbose,db');
    expect(logs).toHaveLength(1);
    expect(logs[0].request).toBe('GET /projects/x');
    expect(logs[0].log.log[0].message).toBe('select 1');
  });

  it('sends no header without --log, and passes a plain response through', async () => {
    const spy = stub(() => new Response('{"id":"x"}', {headers: {'content-type': 'application/json'}}));
    const {c, logs} = client();
    expect(await c.get('/projects/x')).toEqual({id: 'x'});
    expect(spy.mock.calls[0][1].headers[LOG_HEADER]).toBeUndefined();
    expect(logs).toHaveLength(0);
  });

  it('still reports an error the envelope carries, after printing its log', async () => {
    stub(() => envelope('{"#type":"ApiError","message":"No such function"}', 'application/json', {'api-error-code': '404'}));
    const {c, logs} = client('verbose');
    await expect(c.post('/public/v1/functions/X/call', {})).rejects.toThrow('No such function');
    expect(logs).toHaveLength(1);
  });

  it('restores a binary body from base64', async () => {
    stub(() => envelope(Buffer.from([0x44, 0x34, 0x32, 0x00, 0xff]), 'application/octet-stream'));
    const {c} = client('verbose');
    expect([...await c.getBytes('/tables/x/data')]).toEqual([0x44, 0x34, 0x32, 0x00, 0xff]);
  });

  it('reads an empty body as no result', async () => {
    stub(() => envelope('', null));
    const {c} = client('verbose');
    expect(await c.post('/users/x/block')).toBeNull();
  });

  it('is verbose when given without a value', () => {
    expect(logSpecOf({log: true})).toBe('verbose');
    expect(logSpecOf({log: 'debug,db'})).toBe('debug,db');
    expect(logSpecOf({})).toBeUndefined();
  });
});

describe('functions run --debug', () => {
  it('asks for a debug call and prints the log it brings back', async () => {
    const spy = stub(() => envelope('42', 'application/json'));
    const {c, logs} = client();
    expect(await new NodeFuncsDataSource(c).run('Pkg:query', {}, true)).toBe(42);
    expect(spy.mock.calls[0][0]).toBe(`${URL}/public/v1/functions/Pkg.query/call?debug=true`);
    expect(logs).toHaveLength(1);
  });

  it('runs a plain call without it', async () => {
    const spy = stub(() => new Response('42', {headers: {'content-type': 'application/json'}}));
    await new NodeFuncsDataSource(client().c).run('Pkg:query', {});
    expect(spy.mock.calls[0][0]).toBe(`${URL}/public/v1/functions/Pkg.query/call`);
  });
});
