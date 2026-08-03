import {describe, it, expect} from 'vitest';
import {NodeConnectionsDataSource} from '../utils/node-dapi';

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

const CONN = {name: 'TestPg', dataSource: 'PostgreSQL', server: 'localhost', port: 5432, db: 'datagrok'};

describe('NodeConnectionsDataSource.save', () => {
  it('POSTs to /public/v1/connections without saveCredentials by default', async () => {
    const {client, calls} = makeMock((method, path, body) => {
      if (method === 'POST' && path === '/public/v1/connections') return {...body, id: 'new-id'};
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeConnectionsDataSource(client);
    const saved = await ds.save(CONN);
    expect(saved.id).toBe('new-id');
    expect(calls[0].path).toBe('/public/v1/connections');
    expect(calls[0].body).toEqual(CONN);
  });

  it('appends saveCredentials=true when requested', async () => {
    const {client, calls} = makeMock((method, path, body) => {
      if (method === 'POST' && path === '/public/v1/connections?saveCredentials=true') return {...body, id: 'new-id'};
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeConnectionsDataSource(client);
    await ds.save(CONN, true);
    expect(calls[0].path).toBe('/public/v1/connections?saveCredentials=true');
  });

  it('omits saveCredentials query param when false', async () => {
    const {client, calls} = makeMock((method, path) => {
      if (method === 'POST' && path === '/public/v1/connections') return {id: 'new-id'};
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeConnectionsDataSource(client);
    await ds.save(CONN, false);
    expect(calls[0].path).toBe('/public/v1/connections');
  });
});

describe('NodeConnectionsDataSource.test', () => {
  it('POSTs the connection body to /public/v1/connections/test and resolves on "ok"', async () => {
    const {client, calls} = makeMock((method, path) => {
      if (method === 'POST' && path === '/public/v1/connections/test') return 'ok';
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeConnectionsDataSource(client);
    await ds.test(CONN);
    expect(calls[0].body).toEqual(CONN);
  });

  it('strips surrounding quotes on a JSON-encoded "ok"', async () => {
    const {client} = makeMock(() => '"ok"');
    const ds = new NodeConnectionsDataSource(client);
    await expect(ds.test(CONN)).resolves.toBeUndefined();
  });

  it('throws the server error text when the response is not ok', async () => {
    const {client} = makeMock(() => 'connection refused: tcp://localhost:5432');
    const ds = new NodeConnectionsDataSource(client);
    await expect(ds.test(CONN)).rejects.toThrow(/connection refused/);
  });

  it('throws a generic message on empty response', async () => {
    const {client} = makeMock(() => '');
    const ds = new NodeConnectionsDataSource(client);
    await expect(ds.test(CONN)).rejects.toThrow(/failed/i);
  });
});
