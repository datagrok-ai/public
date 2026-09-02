import {describe, it, expect, afterEach, vi} from 'vitest';
import {InternalDataSource, NodeApiClient} from '../utils/node-dapi';

interface Call {method: string; path: string; body?: any}

function makeMock(responder: (method: string, path: string, body?: any) => any) {
  const calls: Call[] = [];
  const client: any = {
    async request(method: string, path: string, body?: any) {
      calls.push({method, path, body});
      return responder(method, path, body);
    },
    get(path: string) { return this.request('GET', path); },
    post(path: string, body?: any) { return this.request('POST', path, body); },
    del(path: string) { return this.request('DELETE', path); },
  };
  return {client, calls};
}

describe('InternalDataSource', () => {
  it('lists with the smart filter and an explicit page size', async () => {
    const {client, calls} = makeMock(() => []);
    await new InternalDataSource(client, '/connectors/queries').listAll({text: 'name = "NW"'});
    expect(calls[0].path).toBe('/connectors/queries?text=name%20%3D%20%22NW%22&limit=500&page=1');
  });

  it('pages from 1 — the server treats page 0 as page 1 — until a short page comes back', async () => {
    const page1 = Array.from({length: 2}, (_, i) => ({id: `a${i}`}));
    const {client, calls} = makeMock((_m, path) => path.includes('page=2') ? [{id: 'b'}] : page1);
    const all = await new InternalDataSource(client, '/scripts').listAll({}, 2);
    expect(all.map((e) => e.id)).toEqual(['a0', 'a1', 'b']);
    expect(calls.map((c) => c.path)).toEqual(['/scripts?limit=2&page=1', '/scripts?limit=2&page=2']);
  });

  it('finds by id with include', async () => {
    const {client, calls} = makeMock(() => ({id: 'x'}));
    await new InternalDataSource(client, '/projects').find('x', 'relations');
    expect(calls[0].path).toBe('/projects/x?include=relations');
  });

  it('returns null instead of throwing when the entity is gone', async () => {
    const {client} = makeMock(() => { throw Object.assign(new Error('Not Found'), {apiError: {error: 'Not Found', errorCode: 404}}); });
    expect(await new InternalDataSource(client, '/projects').find('x')).toBeNull();
  });

  it('rethrows anything that is not a 404', async () => {
    const {client} = makeMock(() => { throw Object.assign(new Error('boom'), {apiError: {error: 'boom', errorCode: 500}}); });
    await expect(new InternalDataSource(client, '/projects').find('x')).rejects.toThrow('boom');
  });

  it('saves with a query string and generates a missing id', async () => {
    const {client, calls} = makeMock((_m, _p, body) => body);
    await new InternalDataSource(client, '/projects').save({name: 'X'}, 'saveRelations=true');
    expect(calls[0].method).toBe('POST');
    expect(calls[0].path).toBe('/projects?saveRelations=true');
    expect(calls[0].body.id).toMatch(/^[0-9a-f-]{36}$/);
  });

  it('deletes by id', async () => {
    const {client, calls} = makeMock(() => 'ok');
    await new InternalDataSource(client, '/tables').delete('abc');
    expect(calls[0]).toMatchObject({method: 'DELETE', path: '/tables/abc'});
  });
});

describe('NodeApiClient.getBytes', () => {
  afterEach(() => vi.unstubAllGlobals());

  it('returns the response body as a Buffer', async () => {
    const payload = Buffer.from([0x29, 0xa0, 0x00, 0x01]);
    vi.stubGlobal('fetch', vi.fn(async (url: string, opts: any) => {
      expect(url).toBe('http://h/api/tables/t1/data');
      expect(opts.headers['Authorization']).toBe('tok');
      return new Response(payload, {status: 200});
    }));
    const bytes = await new NodeApiClient('http://h/api', 'tok').getBytes('/tables/t1/data');
    expect(Buffer.compare(bytes, payload)).toBe(0);
  });

  it('throws the API error shape on a failed response', async () => {
    vi.stubGlobal('fetch', vi.fn(async () =>
      new Response('{"message":"kaboom","errorCode":500}', {status: 500, headers: {'content-type': 'application/json'}})));
    await expect(new NodeApiClient('http://h/api', 'tok').getBytes('/tables/t1/data')).rejects.toThrow('kaboom');
  });
});
