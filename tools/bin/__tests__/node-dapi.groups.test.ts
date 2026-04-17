import {describe, it, expect, beforeEach} from 'vitest';
import {NodeGroupsDataSource} from '../utils/node-dapi';

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

const PARENT_ID = 'aaaaaaaa-aaaa-aaaa-aaaa-aaaaaaaaaaaa';
const ALICE_ID = 'bbbbbbbb-bbbb-bbbb-bbbb-bbbbbbbbbbbb';
const BOB_ID = 'cccccccc-cccc-cccc-cccc-cccccccccccc';

function makeParent(children: any[] = []) {
  return {id: PARENT_ID, friendlyName: 'Admins', children};
}

describe('NodeGroupsDataSource.resolve', () => {
  it('treats a UUID as a direct ID and calls find()', async () => {
    const {client, calls} = makeMock((_m, path) => {
      if (path.startsWith(`/public/v1/groups/${PARENT_ID}`)) return makeParent();
      throw new Error(`unexpected ${path}`);
    });
    const ds = new NodeGroupsDataSource(client);
    const g = await ds.resolve(PARENT_ID);
    expect(g.id).toBe(PARENT_ID);
    expect(calls[0].path).toBe(`/public/v1/groups/${PARENT_ID}`);
  });

  it('uses /lookup for non-UUID names and returns the unique match', async () => {
    const {client, calls} = makeMock((_m, path) => {
      if (path.startsWith('/public/v1/groups/lookup'))
        return [{id: ALICE_ID, name: 'alice', friendlyName: 'alice', personal: true}];
      throw new Error(`unexpected ${path}`);
    });
    const ds = new NodeGroupsDataSource(client);
    const g = await ds.resolve('alice');
    expect(g.id).toBe(ALICE_ID);
    expect(calls[0].path).toBe('/public/v1/groups/lookup?query=alice');
  });

  it('errors with a list when lookup returns multiple matches', async () => {
    const {client} = makeMock((_m, path) => {
      if (path.startsWith('/public/v1/groups/lookup'))
        return [
          {id: 'g1', friendlyName: 'alice', personal: false},
          {id: 'g2', friendlyName: 'alice', personal: true},
        ];
      throw new Error(`unexpected ${path}`);
    });
    const ds = new NodeGroupsDataSource(client);
    await expect(ds.resolve('alice')).rejects.toThrow(/Multiple groups match 'alice'/);
  });

  it('personalOnly filters lookup results to personal groups', async () => {
    const {client} = makeMock((_m, path) => {
      if (path.startsWith('/public/v1/groups/lookup'))
        return [
          {id: 'g1', friendlyName: 'alice', personal: false},
          {id: ALICE_ID, friendlyName: 'alice', personal: true},
        ];
      throw new Error(`unexpected ${path}`);
    });
    const ds = new NodeGroupsDataSource(client);
    const g = await ds.resolve('alice', {personalOnly: true});
    expect(g.id).toBe(ALICE_ID);
  });

  it('errors when personalOnly yields no matches', async () => {
    const {client} = makeMock((_m, path) => {
      if (path.startsWith('/public/v1/groups/lookup'))
        return [{id: 'g1', friendlyName: 'alice', personal: false}];
      throw new Error(`unexpected ${path}`);
    });
    const ds = new NodeGroupsDataSource(client);
    await expect(ds.resolve('alice', {personalOnly: true})).rejects.toThrow(/No group matching 'alice' \(personal\)/);
  });
});

describe('NodeGroupsDataSource.addMembers', () => {
  it('appends new relations and POSTs with saveRelations=true', async () => {
    const {client, calls} = makeMock((method, path, _body) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent();
      if (method === 'GET' && path.startsWith('/public/v1/groups/lookup'))
        return [{id: ALICE_ID, friendlyName: 'alice', personal: true}];
      if (method === 'POST' && path === '/public/v1/groups?saveRelations=true') return {id: PARENT_ID};
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeGroupsDataSource(client);
    const results = await ds.addMembers(PARENT_ID, ['alice']);
    expect(results).toEqual([{member: 'alice', status: 'added'}]);

    const post = calls.find((c) => c.method === 'POST')!;
    expect(post.path).toBe('/public/v1/groups?saveRelations=true');
    expect(post.body.children).toEqual([{parent: {id: PARENT_ID}, child: {id: ALICE_ID}, isAdmin: false}]);
  });

  it('supports --admin by setting isAdmin on the new relation', async () => {
    const {client, calls} = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent();
      if (method === 'GET' && path.startsWith('/public/v1/groups/lookup'))
        return [{id: ALICE_ID, friendlyName: 'alice', personal: true}];
      if (method === 'POST') return {id: PARENT_ID};
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeGroupsDataSource(client);
    await ds.addMembers(PARENT_ID, ['alice'], true);
    const post = calls.find((c) => c.method === 'POST')!;
    expect(post.body.children[0].isAdmin).toBe(true);
  });

  it('is idempotent: re-adding with the same role reports noop and skips POST', async () => {
    const existing = {parent: {id: PARENT_ID}, child: {id: ALICE_ID}, isAdmin: false};
    const {client, calls} = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent([existing]);
      if (method === 'GET' && path.startsWith('/public/v1/groups/lookup'))
        return [{id: ALICE_ID, friendlyName: 'alice', personal: true}];
      if (method === 'POST') return {id: PARENT_ID};
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeGroupsDataSource(client);
    const results = await ds.addMembers(PARENT_ID, ['alice']);
    expect(results).toEqual([{member: 'alice', status: 'noop'}]);
    expect(calls.find((c) => c.method === 'POST')).toBeUndefined();
  });

  it('flips isAdmin when re-adding an existing member with a different role', async () => {
    const existing = {parent: {id: PARENT_ID}, child: {id: ALICE_ID}, isAdmin: false};
    const {client, calls} = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent([existing]);
      if (method === 'GET' && path.startsWith('/public/v1/groups/lookup'))
        return [{id: ALICE_ID, friendlyName: 'alice', personal: true}];
      if (method === 'POST') return {id: PARENT_ID};
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeGroupsDataSource(client);
    const results = await ds.addMembers(PARENT_ID, ['alice'], true);
    expect(results).toEqual([{member: 'alice', status: 'updated'}]);
    const post = calls.find((c) => c.method === 'POST')!;
    expect(post.body.children).toHaveLength(1);
    expect(post.body.children[0].isAdmin).toBe(true);
  });

  it('processes a mixed batch: one valid, one unresolvable', async () => {
    const {client, calls} = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent();
      if (method === 'GET' && path === '/public/v1/groups/lookup?query=alice')
        return [{id: ALICE_ID, friendlyName: 'alice', personal: true}];
      if (method === 'GET' && path === '/public/v1/groups/lookup?query=nobody')
        return [];
      if (method === 'POST') return {id: PARENT_ID};
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeGroupsDataSource(client);
    const results = await ds.addMembers(PARENT_ID, ['alice', 'nobody']);
    expect(results[0]).toEqual({member: 'alice', status: 'added'});
    expect(results[1].member).toBe('nobody');
    expect(results[1].status).toBe('error');
    expect(results[1].error).toMatch(/No group matching 'nobody'/);
    // POST should still happen for the valid member
    const post = calls.find((c) => c.method === 'POST')!;
    expect(post.body.children).toHaveLength(1);
    expect(post.body.children[0].child.id).toBe(ALICE_ID);
  });

  it('handles multiple members in a single call', async () => {
    const {client, calls} = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent();
      if (method === 'GET' && path === '/public/v1/groups/lookup?query=alice')
        return [{id: ALICE_ID, friendlyName: 'alice', personal: true}];
      if (method === 'GET' && path === '/public/v1/groups/lookup?query=bob')
        return [{id: BOB_ID, friendlyName: 'bob', personal: true}];
      if (method === 'POST') return {id: PARENT_ID};
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeGroupsDataSource(client);
    const results = await ds.addMembers(PARENT_ID, ['alice', 'bob']);
    expect(results.map((r) => r.status)).toEqual(['added', 'added']);
    const post = calls.find((c) => c.method === 'POST')!;
    expect(post.body.children.map((c: any) => c.child.id)).toEqual([ALICE_ID, BOB_ID]);
    // exactly one POST regardless of batch size
    expect(calls.filter((c) => c.method === 'POST')).toHaveLength(1);
  });
});

describe('NodeGroupsDataSource.removeMembers', () => {
  it('filters the relation out and POSTs with saveRelations=true', async () => {
    const existing = {parent: {id: PARENT_ID}, child: {id: ALICE_ID}, isAdmin: false};
    const {client, calls} = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent([existing]);
      if (method === 'GET' && path.startsWith('/public/v1/groups/lookup'))
        return [{id: ALICE_ID, friendlyName: 'alice', personal: true}];
      if (method === 'POST') return {id: PARENT_ID};
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeGroupsDataSource(client);
    const results = await ds.removeMembers(PARENT_ID, ['alice']);
    expect(results).toEqual([{member: 'alice', status: 'removed'}]);
    const post = calls.find((c) => c.method === 'POST')!;
    expect(post.path).toBe('/public/v1/groups?saveRelations=true');
    expect(post.body.children).toEqual([]);
  });

  it('reports not-member and skips POST when the member isn\'t in the group', async () => {
    const {client, calls} = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent();
      if (method === 'GET' && path.startsWith('/public/v1/groups/lookup'))
        return [{id: ALICE_ID, friendlyName: 'alice', personal: true}];
      if (method === 'POST') return {id: PARENT_ID};
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeGroupsDataSource(client);
    const results = await ds.removeMembers(PARENT_ID, ['alice']);
    expect(results).toEqual([{member: 'alice', status: 'not-member'}]);
    expect(calls.find((c) => c.method === 'POST')).toBeUndefined();
  });
});

describe('NodeGroupsDataSource.getMembers', () => {
  it('calls /members without an admin query when admin is undefined', async () => {
    const {client, calls} = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent();
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}/members`) return [];
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeGroupsDataSource(client);
    await ds.getMembers(PARENT_ID);
    expect(calls.map((c) => c.path)).toContain(`/public/v1/groups/${PARENT_ID}/members`);
  });

  it('passes admin=true', async () => {
    const {client, calls} = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent();
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}/members?admin=true`) return [];
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeGroupsDataSource(client);
    await ds.getMembers(PARENT_ID, true);
    expect(calls.map((c) => c.path)).toContain(`/public/v1/groups/${PARENT_ID}/members?admin=true`);
  });

  it('passes admin=false', async () => {
    const {client, calls} = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent();
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}/members?admin=false`) return [];
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeGroupsDataSource(client);
    await ds.getMembers(PARENT_ID, false);
    expect(calls.map((c) => c.path)).toContain(`/public/v1/groups/${PARENT_ID}/members?admin=false`);
  });
});

describe('NodeGroupsDataSource.getMemberships', () => {
  it('calls /memberships with admin flag when provided', async () => {
    const {client, calls} = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent();
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}/memberships?admin=true`) return [];
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new NodeGroupsDataSource(client);
    await ds.getMemberships(PARENT_ID, true);
    expect(calls.map((c) => c.path)).toContain(`/public/v1/groups/${PARENT_ID}/memberships?admin=true`);
  });
});
