import {describe, it, expect} from 'vitest';
import {NodeSharesDataSource} from '../utils/node-dapi';

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

const ENTITY_UUID = 'aaaaaaaa-aaaa-aaaa-aaaa-aaaaaaaaaaaa';

describe('NodeSharesDataSource.share', () => {
  it('POSTs to /public/v1/entities/<name>/shares with groups and access in the query', async () => {
    const {client, calls} = makeMock((method, path) => {
      if (method === 'POST' && path.startsWith('/public/v1/entities/')) return {status: 'success'};
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeSharesDataSource(client);
    await ds.share('JohnDoe:MyConnection', 'Chemists,Admins', 'Edit');
    expect(calls[0].method).toBe('POST');
    expect(calls[0].path).toBe('/public/v1/entities/JohnDoe.MyConnection/shares?groups=Chemists%2CAdmins&access=Edit');
    expect(calls[0].body).toBeUndefined();
  });

  it('converts colon in grok name to dot for the URL path', async () => {
    const {client, calls} = makeMock(() => ({status: 'success'}));
    const ds = new NodeSharesDataSource(client);
    await ds.share('Pkg:Q', 'Chemists', 'View');
    expect(calls[0].path).toMatch(/\/public\/v1\/entities\/Pkg\.Q\/shares\?/);
  });

  it('defaults access to View when not specified', async () => {
    const {client, calls} = makeMock(() => ({status: 'success'}));
    const ds = new NodeSharesDataSource(client);
    await ds.share('Pkg:Q', 'Chemists');
    expect(calls[0].path).toContain('access=View');
  });

  it('passes a UUID entity id through unchanged in the path', async () => {
    const {client, calls} = makeMock(() => ({status: 'success'}));
    const ds = new NodeSharesDataSource(client);
    await ds.share(ENTITY_UUID, 'Chemists', 'View');
    expect(calls[0].path).toBe(`/public/v1/entities/${ENTITY_UUID}/shares?groups=Chemists&access=View`);
  });
});

describe('NodeSharesDataSource.list', () => {
  it('GETs /privileges/permissions?entityId=<id>', async () => {
    const {client, calls} = makeMock((method, path) => {
      if (method === 'GET' && path === `/privileges/permissions?entityId=${ENTITY_UUID}`)
        return [{id: 'p1', userGroup: {id: 'g1', friendlyName: 'Chemists'}, permission: {name: 'View'}}];
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeSharesDataSource(client);
    const perms = await ds.list(ENTITY_UUID);
    expect(perms).toHaveLength(1);
    expect(perms[0].userGroup.friendlyName).toBe('Chemists');
    expect(calls[0].path).toBe(`/privileges/permissions?entityId=${ENTITY_UUID}`);
  });
});
