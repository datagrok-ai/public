import {describe, it, expect} from 'vitest';
import {
  throwIfApiError, apiPath, mapPositionalParams, NodeUsersDataSource, NodeFuncsDataSource,
  NodeTablesDataSource, NodeFilesDataSource, NodeHttpDataSource,
} from '../utils/node-dapi';

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

const TABLE_ID = 'dddddddd-dddd-dddd-dddd-dddddddddddd';
const CONN_ID = 'eeeeeeee-eeee-eeee-eeee-eeeeeeeeeeee';

describe('throwIfApiError', () => {
  it('throws on an ApiError object and carries the message and code', () => {
    expect(() => throwIfApiError({'#type': 'ApiError', message: 'Not Found', errorCode: 404}))
      .toThrow('Not Found');
    try { throwIfApiError({'#type': 'ApiError', message: 'Not Found', errorCode: 404}); }
    catch (e: any) { expect(e.apiError.errorCode).toBe(404); }
  });

  it('throws on an ApiError serialized as a JSON string', () => {
    expect(() => throwIfApiError('{"#type":"ApiError","message":"Login already exists"}'))
      .toThrow('Login already exists');
  });

  it('passes everything else through untouched', () => {
    const obj = {'#type': 'User', login: 'a'};
    expect(throwIfApiError(obj)).toBe(obj);
    expect(throwIfApiError('a,b\n1,2')).toBe('a,b\n1,2');
    expect(throwIfApiError(null)).toBeNull();
    expect(throwIfApiError([1, 2])).toEqual([1, 2]);
  });
});

describe('apiPath', () => {
  it('drops a leading /api so the path works against /api bases and bare Datlas', () => {
    expect(apiPath('/api/users/current')).toBe('/users/current');
    expect(apiPath('/users/current')).toBe('/users/current');
    expect(apiPath('users/current')).toBe('/users/current');
    expect(apiPath('/api')).toBe('/');
  });

  it('does not touch paths that merely start with "api"', () => {
    expect(apiPath('/apikeys')).toBe('/apikeys');
  });
});

describe('mapPositionalParams', () => {
  const infos = {
    condition: {name: 'condition', isInput: null},
    ifTrue: {name: 'ifTrue'},
    ifFalse: {name: 'ifFalse'},
    result: {name: 'result', isInput: false},
  };

  it('maps positional values onto the inputs in declared order, skipping outputs', () => {
    expect(mapPositionalParams({0: true, 1: 1, 2: 2}, infos, 'If')).toEqual({condition: true, ifTrue: 1, ifFalse: 2});
  });

  it('keeps named keys and fills the rest by position', () => {
    expect(mapPositionalParams({0: true, ifFalse: 2}, infos, 'If')).toEqual({condition: true, ifFalse: 2});
  });

  it('leaves an all-named map untouched without needing parameterInfos', () => {
    const named = {x: 1};
    expect(mapPositionalParams(named, undefined, 'Sin')).toBe(named);
  });

  it('rejects more positional values than inputs', () => {
    expect(() => mapPositionalParams({0: 1, 1: 2}, {x: {name: 'x'}, result: {name: 'result', isInput: false}}, 'Sin'))
      .toThrow('Sin takes 1 input (x), got 2');
  });
});

describe('NodeHttpDataSource paging', () => {
  it('sends a 1-based page to the list route and the filter to /count', async () => {
    const {client, calls} = makeMock(() => []);
    const ds = new NodeHttpDataSource(client, 'users', '/users');
    await ds.filter('status = "active"').by(10).page(2).list();
    expect(calls[0].path).toBe('/users?text=status%20%3D%20%22active%22&limit=10&page=3');
    await ds.count();
    expect(calls[1].path).toBe('/users/count?text=status%20%3D%20%22active%22');
  });

  it('falls back to the public route without a list route', async () => {
    const {client, calls} = makeMock(() => []);
    await new NodeHttpDataSource(client, 'packages').by(5).list();
    expect(calls[0].path).toBe('/public/v1/packages?limit=5&page=1');
  });
});

describe('NodeUsersDataSource.delete', () => {
  it('refuses: the server has no user deletion', async () => {
    const {client, calls} = makeMock(() => ({}));
    await expect(new NodeUsersDataSource(client).delete('alice')).rejects.toThrow('grok s users block alice');
    expect(calls).toHaveLength(0);
  });
});

describe('NodeFuncsDataSource.delete', () => {
  it('routes a script to /scripts and a query to /connectors/queries', async () => {
    const {client, calls} = makeMock((method, path) => {
      if (path === '/public/v1/functions/Pkg.s') return {'#type': 'Script', id: 's1', name: 's'};
      if (path === '/public/v1/functions/Pkg.q') return {'#type': 'DataQuery', id: 'q1', name: 'q'};
      return null;
    });
    const ds = new NodeFuncsDataSource(client);
    await ds.delete('Pkg:s');
    await ds.delete('Pkg:q');
    expect(calls.filter((c) => c.method === 'DELETE').map((c) => c.path)).toEqual(['/scripts/s1', '/connectors/queries/q1']);
  });

  it('refuses a plain or package function', async () => {
    const {client} = makeMock(() => ({'#type': 'Func', id: 'f', name: 'If'}));
    await expect(new NodeFuncsDataSource(client).delete('If')).rejects.toThrow('Only scripts and queries');
  });
});

describe('NodeTablesDataSource', () => {
  const tables = [
    {'#type': 'TableInfo', id: TABLE_ID, name: 'Demo', namespace: 'Admin:Demo:'},
    {'#type': 'TableInfo', id: '11111111-1111-1111-1111-111111111111', name: 'DemoOld', namespace: 'Admin:DemoOld:'},
  ];

  it('finds by bare name, full name, or id', async () => {
    const {client, calls} = makeMock((_m, path) => path.startsWith('/tables/') ? tables[0] : tables);
    const ds = new NodeTablesDataSource(client);
    expect((await ds.find('Demo')).id).toBe(TABLE_ID);
    expect((await ds.find('Admin:Demo:Demo')).id).toBe(TABLE_ID);
    await ds.find(TABLE_ID);
    expect(calls[2].path).toBe(`/tables/${TABLE_ID}`);
  });

  it('reports a missing or ambiguous name', async () => {
    const {client} = makeMock(() => [tables[0], {...tables[0], id: '22222222-2222-2222-2222-222222222222', namespace: 'Bob:Demo:'}]);
    const ds = new NodeTablesDataSource(client);
    await expect(ds.find('Nope')).rejects.toThrow("No table named 'Nope'");
    await expect(ds.find('Demo')).rejects.toThrow('Multiple tables match');
  });

  it('downloads and deletes through the resolved id', async () => {
    const {client, calls} = makeMock((_m, path) => path.startsWith('/public/v1/tables/') ? 'a,b\n1,2' : tables);
    const ds = new NodeTablesDataSource(client);
    expect(await ds.download('Demo')).toBe('a,b\n1,2');
    await ds.delete('Demo');
    expect(calls.map((c) => `${c.method} ${c.path}`)).toContain(`GET /public/v1/tables/${TABLE_ID}`);
    expect(calls.map((c) => `${c.method} ${c.path}`)).toContain(`DELETE /tables/${TABLE_ID}`);
  });
});

describe('NodeFilesDataSource.list', () => {
  it('resolves the connector and lists through the internal files route', async () => {
    const {client, calls} = makeMock((_m, path) => {
      if (path === '/public/v1/connections/System.DemoFiles') return {id: CONN_ID};
      return [{'#type': 'FileInfo', path: 'geo/a.csv', isFile: true}];
    });
    const files = await new NodeFilesDataSource(client).list('System:DemoFiles/geo', true);
    expect(files[0].path).toBe('geo/a.csv');
    expect(calls[1].path).toBe(`/connectors/connections/${CONN_ID}/files/geo/?recursive=true`);
  });

  it('lists the connector root with a trailing slash', async () => {
    const {client, calls} = makeMock((_m, path) => path.includes('/public/v1/') ? {id: CONN_ID} : []);
    await new NodeFilesDataSource(client).list('System:AppData');
    expect(calls[1].path).toBe(`/connectors/connections/${CONN_ID}/files/`);
  });
});
