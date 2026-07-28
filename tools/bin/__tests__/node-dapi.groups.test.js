"use strict";

var _vitest = require("vitest");
var _nodeDapi = require("../utils/node-dapi");
function makeMock(responder) {
  const calls = [];
  const client = {
    async request(method, path, body) {
      calls.push({
        method,
        path,
        body
      });
      return responder(method, path, body);
    },
    get(path) {
      return this.request('GET', path);
    },
    post(path, body) {
      return this.request('POST', path, body);
    },
    del(path) {
      return this.request('DELETE', path);
    }
  };
  return {
    client,
    calls
  };
}
const PARENT_ID = 'aaaaaaaa-aaaa-aaaa-aaaa-aaaaaaaaaaaa';
const ALICE_ID = 'bbbbbbbb-bbbb-bbbb-bbbb-bbbbbbbbbbbb';
const BOB_ID = 'cccccccc-cccc-cccc-cccc-cccccccccccc';
function makeParent(children = []) {
  return {
    id: PARENT_ID,
    friendlyName: 'Admins',
    children
  };
}
(0, _vitest.describe)('NodeGroupsDataSource.resolve', () => {
  (0, _vitest.it)('treats a UUID as a direct ID and calls find()', async () => {
    const {
      client,
      calls
    } = makeMock((_m, path) => {
      if (path.startsWith(`/public/v1/groups/${PARENT_ID}`)) return makeParent();
      throw new Error(`unexpected ${path}`);
    });
    const ds = new _nodeDapi.NodeGroupsDataSource(client);
    const g = await ds.resolve(PARENT_ID);
    (0, _vitest.expect)(g.id).toBe(PARENT_ID);
    (0, _vitest.expect)(calls[0].path).toBe(`/public/v1/groups/${PARENT_ID}`);
  });
  (0, _vitest.it)('uses /lookup for non-UUID names and returns the unique match', async () => {
    const {
      client,
      calls
    } = makeMock((_m, path) => {
      if (path.startsWith('/public/v1/groups/lookup')) return [{
        id: ALICE_ID,
        name: 'alice',
        friendlyName: 'alice',
        personal: true
      }];
      throw new Error(`unexpected ${path}`);
    });
    const ds = new _nodeDapi.NodeGroupsDataSource(client);
    const g = await ds.resolve('alice');
    (0, _vitest.expect)(g.id).toBe(ALICE_ID);
    (0, _vitest.expect)(calls[0].path).toBe('/public/v1/groups/lookup?query=alice');
  });
  (0, _vitest.it)('errors with a list when lookup returns multiple matches', async () => {
    const {
      client
    } = makeMock((_m, path) => {
      if (path.startsWith('/public/v1/groups/lookup')) return [{
        id: 'g1',
        friendlyName: 'alice',
        personal: false
      }, {
        id: 'g2',
        friendlyName: 'alice',
        personal: true
      }];
      throw new Error(`unexpected ${path}`);
    });
    const ds = new _nodeDapi.NodeGroupsDataSource(client);
    await (0, _vitest.expect)(ds.resolve('alice')).rejects.toThrow(/Multiple groups match 'alice'/);
  });
  (0, _vitest.it)('personalOnly filters lookup results to personal groups', async () => {
    const {
      client
    } = makeMock((_m, path) => {
      if (path.startsWith('/public/v1/groups/lookup')) return [{
        id: 'g1',
        friendlyName: 'alice',
        personal: false
      }, {
        id: ALICE_ID,
        friendlyName: 'alice',
        personal: true
      }];
      throw new Error(`unexpected ${path}`);
    });
    const ds = new _nodeDapi.NodeGroupsDataSource(client);
    const g = await ds.resolve('alice', {
      personalOnly: true
    });
    (0, _vitest.expect)(g.id).toBe(ALICE_ID);
  });
  (0, _vitest.it)('errors when personalOnly yields no matches', async () => {
    const {
      client
    } = makeMock((_m, path) => {
      if (path.startsWith('/public/v1/groups/lookup')) return [{
        id: 'g1',
        friendlyName: 'alice',
        personal: false
      }];
      throw new Error(`unexpected ${path}`);
    });
    const ds = new _nodeDapi.NodeGroupsDataSource(client);
    await (0, _vitest.expect)(ds.resolve('alice', {
      personalOnly: true
    })).rejects.toThrow(/No group matching 'alice' \(personal\)/);
  });
});
(0, _vitest.describe)('NodeGroupsDataSource.save', () => {
  (0, _vitest.it)('POSTs the group body to /public/v1/groups and auto-generates an id when missing', async () => {
    const {
      client,
      calls
    } = makeMock((method, path, body) => {
      if (method === 'POST' && path === '/public/v1/groups') return {
        ...body,
        friendlyName: 'Chemists'
      };
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new _nodeDapi.NodeGroupsDataSource(client);
    await ds.save({
      friendlyName: 'Chemists'
    });
    (0, _vitest.expect)(calls[0].path).toBe('/public/v1/groups');
    (0, _vitest.expect)(calls[0].body.id).toMatch(/^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/i);
  });
  (0, _vitest.it)('appends saveRelations=true when requested', async () => {
    const {
      client,
      calls
    } = makeMock((method, path) => {
      if (method === 'POST' && path === '/public/v1/groups?saveRelations=true') return {
        id: 'g'
      };
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new _nodeDapi.NodeGroupsDataSource(client);
    await ds.save({
      id: 'g',
      friendlyName: 'Chemists'
    }, true);
    (0, _vitest.expect)(calls[0].path).toBe('/public/v1/groups?saveRelations=true');
  });
});
(0, _vitest.describe)('NodeGroupsDataSource.addMembers', () => {
  (0, _vitest.it)('appends new relations and POSTs with saveRelations=true', async () => {
    const {
      client,
      calls
    } = makeMock((method, path, _body) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent();
      if (method === 'GET' && path.startsWith('/public/v1/groups/lookup')) return [{
        id: ALICE_ID,
        friendlyName: 'alice',
        personal: true
      }];
      if (method === 'POST' && path === '/public/v1/groups?saveRelations=true') return {
        id: PARENT_ID
      };
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new _nodeDapi.NodeGroupsDataSource(client);
    const results = await ds.addMembers(PARENT_ID, ['alice']);
    (0, _vitest.expect)(results).toEqual([{
      member: 'alice',
      status: 'added'
    }]);
    const post = calls.find(c => c.method === 'POST');
    (0, _vitest.expect)(post.path).toBe('/public/v1/groups?saveRelations=true');
    (0, _vitest.expect)(post.body.children).toEqual([{
      parent: {
        id: PARENT_ID
      },
      child: {
        id: ALICE_ID
      },
      isAdmin: false
    }]);
  });
  (0, _vitest.it)('supports --admin by setting isAdmin on the new relation', async () => {
    const {
      client,
      calls
    } = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent();
      if (method === 'GET' && path.startsWith('/public/v1/groups/lookup')) return [{
        id: ALICE_ID,
        friendlyName: 'alice',
        personal: true
      }];
      if (method === 'POST') return {
        id: PARENT_ID
      };
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new _nodeDapi.NodeGroupsDataSource(client);
    await ds.addMembers(PARENT_ID, ['alice'], true);
    const post = calls.find(c => c.method === 'POST');
    (0, _vitest.expect)(post.body.children[0].isAdmin).toBe(true);
  });
  (0, _vitest.it)('is idempotent: re-adding with the same role reports noop and skips POST', async () => {
    const existing = {
      parent: {
        id: PARENT_ID
      },
      child: {
        id: ALICE_ID
      },
      isAdmin: false
    };
    const {
      client,
      calls
    } = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent([existing]);
      if (method === 'GET' && path.startsWith('/public/v1/groups/lookup')) return [{
        id: ALICE_ID,
        friendlyName: 'alice',
        personal: true
      }];
      if (method === 'POST') return {
        id: PARENT_ID
      };
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new _nodeDapi.NodeGroupsDataSource(client);
    const results = await ds.addMembers(PARENT_ID, ['alice']);
    (0, _vitest.expect)(results).toEqual([{
      member: 'alice',
      status: 'noop'
    }]);
    (0, _vitest.expect)(calls.find(c => c.method === 'POST')).toBeUndefined();
  });
  (0, _vitest.it)('flips isAdmin when re-adding an existing member with a different role', async () => {
    const existing = {
      parent: {
        id: PARENT_ID
      },
      child: {
        id: ALICE_ID
      },
      isAdmin: false
    };
    const {
      client,
      calls
    } = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent([existing]);
      if (method === 'GET' && path.startsWith('/public/v1/groups/lookup')) return [{
        id: ALICE_ID,
        friendlyName: 'alice',
        personal: true
      }];
      if (method === 'POST') return {
        id: PARENT_ID
      };
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new _nodeDapi.NodeGroupsDataSource(client);
    const results = await ds.addMembers(PARENT_ID, ['alice'], true);
    (0, _vitest.expect)(results).toEqual([{
      member: 'alice',
      status: 'updated'
    }]);
    const post = calls.find(c => c.method === 'POST');
    (0, _vitest.expect)(post.body.children).toHaveLength(1);
    (0, _vitest.expect)(post.body.children[0].isAdmin).toBe(true);
  });
  (0, _vitest.it)('processes a mixed batch: one valid, one unresolvable', async () => {
    const {
      client,
      calls
    } = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent();
      if (method === 'GET' && path === '/public/v1/groups/lookup?query=alice') return [{
        id: ALICE_ID,
        friendlyName: 'alice',
        personal: true
      }];
      if (method === 'GET' && path === '/public/v1/groups/lookup?query=nobody') return [];
      if (method === 'POST') return {
        id: PARENT_ID
      };
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new _nodeDapi.NodeGroupsDataSource(client);
    const results = await ds.addMembers(PARENT_ID, ['alice', 'nobody']);
    (0, _vitest.expect)(results[0]).toEqual({
      member: 'alice',
      status: 'added'
    });
    (0, _vitest.expect)(results[1].member).toBe('nobody');
    (0, _vitest.expect)(results[1].status).toBe('error');
    (0, _vitest.expect)(results[1].error).toMatch(/No group matching 'nobody'/);
    // POST should still happen for the valid member
    const post = calls.find(c => c.method === 'POST');
    (0, _vitest.expect)(post.body.children).toHaveLength(1);
    (0, _vitest.expect)(post.body.children[0].child.id).toBe(ALICE_ID);
  });
  (0, _vitest.it)('handles multiple members in a single call', async () => {
    const {
      client,
      calls
    } = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent();
      if (method === 'GET' && path === '/public/v1/groups/lookup?query=alice') return [{
        id: ALICE_ID,
        friendlyName: 'alice',
        personal: true
      }];
      if (method === 'GET' && path === '/public/v1/groups/lookup?query=bob') return [{
        id: BOB_ID,
        friendlyName: 'bob',
        personal: true
      }];
      if (method === 'POST') return {
        id: PARENT_ID
      };
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new _nodeDapi.NodeGroupsDataSource(client);
    const results = await ds.addMembers(PARENT_ID, ['alice', 'bob']);
    (0, _vitest.expect)(results.map(r => r.status)).toEqual(['added', 'added']);
    const post = calls.find(c => c.method === 'POST');
    (0, _vitest.expect)(post.body.children.map(c => c.child.id)).toEqual([ALICE_ID, BOB_ID]);
    // exactly one POST regardless of batch size
    (0, _vitest.expect)(calls.filter(c => c.method === 'POST')).toHaveLength(1);
  });
});
(0, _vitest.describe)('NodeGroupsDataSource.removeMembers', () => {
  (0, _vitest.it)('filters the relation out and POSTs with saveRelations=true', async () => {
    const existing = {
      parent: {
        id: PARENT_ID
      },
      child: {
        id: ALICE_ID
      },
      isAdmin: false
    };
    const {
      client,
      calls
    } = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent([existing]);
      if (method === 'GET' && path.startsWith('/public/v1/groups/lookup')) return [{
        id: ALICE_ID,
        friendlyName: 'alice',
        personal: true
      }];
      if (method === 'POST') return {
        id: PARENT_ID
      };
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new _nodeDapi.NodeGroupsDataSource(client);
    const results = await ds.removeMembers(PARENT_ID, ['alice']);
    (0, _vitest.expect)(results).toEqual([{
      member: 'alice',
      status: 'removed'
    }]);
    const post = calls.find(c => c.method === 'POST');
    (0, _vitest.expect)(post.path).toBe('/public/v1/groups?saveRelations=true');
    (0, _vitest.expect)(post.body.children).toEqual([]);
  });
  (0, _vitest.it)('reports not-member and skips POST when the member isn\'t in the group', async () => {
    const {
      client,
      calls
    } = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent();
      if (method === 'GET' && path.startsWith('/public/v1/groups/lookup')) return [{
        id: ALICE_ID,
        friendlyName: 'alice',
        personal: true
      }];
      if (method === 'POST') return {
        id: PARENT_ID
      };
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new _nodeDapi.NodeGroupsDataSource(client);
    const results = await ds.removeMembers(PARENT_ID, ['alice']);
    (0, _vitest.expect)(results).toEqual([{
      member: 'alice',
      status: 'not-member'
    }]);
    (0, _vitest.expect)(calls.find(c => c.method === 'POST')).toBeUndefined();
  });
});
(0, _vitest.describe)('NodeGroupsDataSource.getMembers', () => {
  (0, _vitest.it)('calls /members without an admin query when admin is undefined', async () => {
    const {
      client,
      calls
    } = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent();
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}/members`) return [];
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new _nodeDapi.NodeGroupsDataSource(client);
    await ds.getMembers(PARENT_ID);
    (0, _vitest.expect)(calls.map(c => c.path)).toContain(`/public/v1/groups/${PARENT_ID}/members`);
  });
  (0, _vitest.it)('passes admin=true', async () => {
    const {
      client,
      calls
    } = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent();
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}/members?admin=true`) return [];
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new _nodeDapi.NodeGroupsDataSource(client);
    await ds.getMembers(PARENT_ID, true);
    (0, _vitest.expect)(calls.map(c => c.path)).toContain(`/public/v1/groups/${PARENT_ID}/members?admin=true`);
  });
  (0, _vitest.it)('passes admin=false', async () => {
    const {
      client,
      calls
    } = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent();
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}/members?admin=false`) return [];
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new _nodeDapi.NodeGroupsDataSource(client);
    await ds.getMembers(PARENT_ID, false);
    (0, _vitest.expect)(calls.map(c => c.path)).toContain(`/public/v1/groups/${PARENT_ID}/members?admin=false`);
  });
});
(0, _vitest.describe)('NodeGroupsDataSource.getMemberships', () => {
  (0, _vitest.it)('calls /memberships with admin flag when provided', async () => {
    const {
      client,
      calls
    } = makeMock((method, path) => {
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}`) return makeParent();
      if (method === 'GET' && path === `/public/v1/groups/${PARENT_ID}/memberships?admin=true`) return [];
      throw new Error(`unexpected ${method} ${path}`);
    });
    const ds = new _nodeDapi.NodeGroupsDataSource(client);
    await ds.getMemberships(PARENT_ID, true);
    (0, _vitest.expect)(calls.map(c => c.path)).toContain(`/public/v1/groups/${PARENT_ID}/memberships?admin=true`);
  });
});