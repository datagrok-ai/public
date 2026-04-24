import {describe, it, expect} from 'vitest';
import {NodeUsersDataSource, ensureBodyId} from '../utils/node-dapi';

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

const UUID_RE = /^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/i;

describe('ensureBodyId', () => {
  it('assigns a UUID when id is missing', () => {
    const body = {name: 'alice'};
    ensureBodyId(body);
    expect((body as any).id).toMatch(UUID_RE);
  });

  it('preserves an existing id', () => {
    const body = {id: 'preset-id', name: 'alice'};
    ensureBodyId(body);
    expect(body.id).toBe('preset-id');
  });
});

describe('NodeUsersDataSource.save', () => {
  it('POSTs to /public/v1/users and generates an id when missing', async () => {
    const {client, calls} = makeMock((method, path, body) => {
      if (method === 'POST' && path === '/public/v1/users') return {...body, name: 'alice'};
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeUsersDataSource(client);
    const result = await ds.save({firstName: 'Alice', lastName: 'A', email: 'a@x.com'});
    expect(calls[0].path).toBe('/public/v1/users');
    expect(calls[0].body.id).toMatch(UUID_RE);
    expect(result.name).toBe('alice');
  });

  it('preserves an existing id in the body', async () => {
    const {client, calls} = makeMock((method, path, body) => {
      if (method === 'POST' && path === '/public/v1/users') return body;
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeUsersDataSource(client);
    await ds.save({id: 'existing-uuid', firstName: 'Alice'});
    expect(calls[0].body.id).toBe('existing-uuid');
  });
});
