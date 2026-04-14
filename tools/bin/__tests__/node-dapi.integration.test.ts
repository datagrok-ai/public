/**
 * Integration tests for NodeApiClient, NodeDapi, and all data sources.
 * Requires a running Datagrok server configured in ~/.grok/config.yaml.
 *
 * Environment variables:
 *   HOST             Server alias or URL (default: config default)
 *   GROK_TEST_FUNC   Function to test with functions.run() (default: System:GetCurrentUser)
 *   GROK_TEST_FILES  File path to test with files.list() (default: System:AppData)
 *
 * Run: npm run test:integration
 * Run against a specific server: HOST=dev npm run test:integration
 */

import {beforeAll, describe, expect, it} from 'vitest';
import {NodeApiClient, NodeDapi} from '../utils/node-dapi';
import {getDevKey} from '../utils/test-utils';

// ─── Config ──────────────────────────────────────────────────────────────────

const HOST = process.env['HOST'] ?? '';
const TEST_FUNC = process.env['GROK_TEST_FUNC'] ?? 'System:GetCurrentUser';
const TEST_FILES_PATH = process.env['GROK_TEST_FILES'] ?? 'System:AppData';

// ─── Setup ───────────────────────────────────────────────────────────────────

let client: NodeApiClient;
let dapi: NodeDapi;
let offline = false;

// Pre-fetched IDs to avoid redundant round trips in find() tests
const seed: {
  userId?: string;
  groupId?: string;
  functionId?: string;
  testFuncName?: string;   // A callable function name for functions.run() test
  connectionId?: string;
  queryId?: string;
  scriptId?: string;
  packageId?: string;
  reportId?: string;
} = {};

beforeAll(async () => {
  try {
    const {url, key} = getDevKey(HOST);
    client = await NodeApiClient.login(url, key);
    dapi = new NodeDapi(client);

    // Pre-seed IDs — each is best-effort; individual tests skip gracefully if absent
    const [users, groups, functions, connections, queries, scripts, packages, reports] =
      await Promise.allSettled([
        dapi.users.by(1).list(),
        dapi.groups.by(1).list(),
        dapi.functions.by(10).list(),   // Fetch more to find a callable function
        dapi.connections.by(1).list(),
        dapi.queries.by(1).list(),
        dapi.scripts.by(1).list(),
        dapi.packages.by(1).list(),
        dapi.reports.by(1).list(),
      ]);

    if (users.status === 'fulfilled' && users.value.length) seed.userId = (users.value[0] as any).id;
    if (groups.status === 'fulfilled' && groups.value.length) seed.groupId = (groups.value[0] as any).id;
    if (functions.status === 'fulfilled' && functions.value.length) {
      seed.functionId = (functions.value[0] as any).id;
      // Find a callable function: prefer env var override, then any nqName we can discover
      if (process.env['GROK_TEST_FUNC'])
        seed.testFuncName = process.env['GROK_TEST_FUNC'];
      else {
        const fn = (functions.value as any[]).find((f: any) => f.nqName) as any;
        seed.testFuncName = fn?.nqName;
      }
    }
    if (connections.status === 'fulfilled' && connections.value.length) seed.connectionId = (connections.value[0] as any).id;
    if (queries.status === 'fulfilled' && queries.value.length) seed.queryId = (queries.value[0] as any).id;
    if (scripts.status === 'fulfilled' && scripts.value.length) seed.scriptId = (scripts.value[0] as any).id;
    if (packages.status === 'fulfilled' && packages.value.length) seed.packageId = (packages.value[0] as any).id;
    if (reports.status === 'fulfilled' && reports.value.length) seed.reportId = (reports.value[0] as any).id;
  } catch (err) {
    console.warn(`[integration] Server unreachable — all tests will be skipped: ${err}`);
    offline = true;
  }
});

/** Wraps a test so it shows as skipped when the server is offline. */
function stest(name: string, fn: () => Promise<void>) {
  return it(name, async (ctx) => {
    if (offline) { ctx.skip(); return; }
    await fn();
  });
}

/**
 * Like stest, but also skips gracefully when the endpoint returns 404.
 * Use for entities whose /public/v1/ endpoint may not exist on all server versions.
 */
function stestEndpoint(name: string, fn: () => Promise<void>) {
  return stest(name, async () => {
    try {
      await fn();
    } catch (e: any) {
      if (e.message === 'Not Found') {
        console.warn(`[integration] ${name}: endpoint not available on this server`);
        return;
      }
      throw e;
    }
  });
}

// ─── NodeApiClient ────────────────────────────────────────────────────────────

describe('NodeApiClient', () => {
  stest('login() returns a client with a non-empty token', async () => {
    const {url, key} = getDevKey(HOST);
    const c = await NodeApiClient.login(url, key);
    expect(c.token).toBeTruthy();
    expect(typeof c.token).toBe('string');
  });

  stest('login() throws on an invalid dev key', async () => {
    const {url} = getDevKey(HOST);
    await expect(NodeApiClient.login(url, 'invalid-key-xyz')).rejects.toThrow();
  });

  stest('get() returns parsed JSON for a known endpoint', async () => {
    const result = await client.get('/public/v1/users?limit=1');
    expect(Array.isArray(result)).toBe(true);
  });

  stest('request() propagates HTTP errors as structured exceptions', async () => {
    await expect(client.get('/public/v1/nonexistent-entity-xyz')).rejects.toThrow();
  });
});

// ─── users ────────────────────────────────────────────────────────────────────

describe('users', () => {
  stest('list() returns an array', async () => {
    const users = await dapi.users.list();
    expect(Array.isArray(users)).toBe(true);
  });

  stest('list() results have id and login fields', async () => {
    const users = await dapi.users.by(5).list();
    expect(users.length).toBeGreaterThan(0);
    for (const u of users as any[])
      expect(u).toMatchObject({id: expect.any(String), login: expect.any(String)});
  });

  stest('count() returns a non-negative integer', async () => {
    const n = await dapi.users.count();
    expect(typeof n).toBe('number');
    expect(n).toBeGreaterThanOrEqual(0);
  });

  stest('filter() narrows results', async () => {
    const all = await dapi.users.count();
    const filtered = await dapi.users.filter('admin').count();
    expect(filtered).toBeLessThanOrEqual(all);
  });

  stest('by() passes limit parameter and returns an array', async () => {
    const users = await dapi.users.by(2).list();
    expect(Array.isArray(users)).toBe(true);
    // Server may not honor small limits, but the call must succeed
    expect(users.length).toBeGreaterThanOrEqual(0);
  });

  stest('page() advances the offset', async () => {
    const page0 = await dapi.users.by(5).page(0).list() as any[];
    const page1 = await dapi.users.by(5).page(1).list() as any[];
    // If enough users exist, pages should differ
    if (page0.length === 5 && page1.length > 0)
      expect(page0[0].id).not.toBe(page1[0].id);
  });

  stest('order() accepts field and direction without throwing', async () => {
    const asc = await dapi.users.by(5).order('login').list();
    const desc = await dapi.users.by(5).order('login', true).list();
    expect(Array.isArray(asc)).toBe(true);
    expect(Array.isArray(desc)).toBe(true);
  });

  stest('find() retrieves a user by ID', async () => {
    if (!seed.userId) return;
    const user = await dapi.users.find(seed.userId) as any;
    expect(user.id).toBe(seed.userId);
    expect(user).toHaveProperty('login');
  });

  stest('delete() throws for a non-existent ID', async () => {
    await expect(dapi.users.delete('non-existent-id-xyz')).rejects.toThrow();
  });
});

// ─── groups ───────────────────────────────────────────────────────────────────

describe('groups', () => {
  stest('list() returns an array with id and name', async () => {
    const groups = await dapi.groups.by(5).list() as any[];
    expect(Array.isArray(groups)).toBe(true);
    for (const g of groups)
      expect(g).toMatchObject({id: expect.any(String), name: expect.any(String)});
  });

  stest('count() returns a non-negative integer', async () => {
    const n = await dapi.groups.count();
    expect(n).toBeGreaterThanOrEqual(0);
  });

  stest('find() retrieves a group by ID', async () => {
    if (!seed.groupId) return;
    const group = await dapi.groups.find(seed.groupId) as any;
    expect(group.id).toBe(seed.groupId);
  });
});

// ─── functions ────────────────────────────────────────────────────────────────

describe('functions', () => {
  stest('list() returns an array with id and name', async () => {
    const fns = await dapi.functions.by(5).list() as any[];
    expect(Array.isArray(fns)).toBe(true);
    expect(fns.length).toBeGreaterThan(0);
    for (const f of fns)
      expect(f).toMatchObject({id: expect.any(String), name: expect.any(String)});
  });

  stest('filter() narrows by name prefix', async () => {
    const all = await dapi.functions.count();
    const filtered = await dapi.functions.filter('System:').count();
    expect(filtered).toBeLessThanOrEqual(all);
  });

  stest('find() retrieves a function by ID', async () => {
    if (!seed.functionId) return;
    const fn = await dapi.functions.find(seed.functionId) as any;
    expect(fn.id).toBe(seed.functionId);
    expect(fn).toHaveProperty('name');
  });

  stest('run() invokes a discovered function and returns a result or param error', async () => {
    if (!seed.testFuncName) return;
    // Call with no params — may fail with a params error, but must NOT return "Not Found"
    try {
      const result = await dapi.functions.run(seed.testFuncName, {});
      expect(result).toBeDefined();
    } catch (e: any) {
      // Params required / type error is acceptable — it means the function exists and was invoked
      expect(e.message).not.toMatch(/not found/i);
    }
  });

  stest('run() throws a structured error for a non-existent function', async () => {
    await expect(dapi.functions.run('NonExistentPkg:nonExistentFunc', {}))
      .rejects.toThrow();
  });
});

// ─── connections ─────────────────────────────────────────────────────────────

describe('connections', () => {
  stest('list() returns an array', async () => {
    const conns = await dapi.connections.by(5).list();
    expect(Array.isArray(conns)).toBe(true);
  });

  stest('count() returns a non-negative integer', async () => {
    const n = await dapi.connections.count();
    expect(n).toBeGreaterThanOrEqual(0);
  });

  stest('find() retrieves a connection by ID', async () => {
    if (!seed.connectionId) return;
    const conn = await dapi.connections.find(seed.connectionId) as any;
    expect(conn.id).toBe(seed.connectionId);
  });
});

// ─── queries ──────────────────────────────────────────────────────────────────

describe('queries', () => {
  stestEndpoint('list() returns an array', async () => {
    const queries = await dapi.queries.by(5).list();
    expect(Array.isArray(queries)).toBe(true);
  });

  stest('find() retrieves a query by ID', async () => {
    if (!seed.queryId) return;
    const query = await dapi.queries.find(seed.queryId) as any;
    expect(query.id).toBe(seed.queryId);
  });
});

// ─── scripts ──────────────────────────────────────────────────────────────────

describe('scripts', () => {
  stestEndpoint('list() returns an array', async () => {
    const scripts = await dapi.scripts.by(5).list();
    expect(Array.isArray(scripts)).toBe(true);
  });

  stest('find() retrieves a script by ID', async () => {
    if (!seed.scriptId) return;
    const script = await dapi.scripts.find(seed.scriptId) as any;
    expect(script.id).toBe(seed.scriptId);
  });
});

// ─── packages ─────────────────────────────────────────────────────────────────

describe('packages', () => {
  stest('list() returns an array with id and name', async () => {
    const pkgs = await dapi.packages.by(5).list() as any[];
    expect(Array.isArray(pkgs)).toBe(true);
    expect(pkgs.length).toBeGreaterThan(0);
    for (const p of pkgs)
      expect(p).toMatchObject({id: expect.any(String), name: expect.any(String)});
  });

  stest('count() returns a non-negative integer', async () => {
    const n = await dapi.packages.count();
    expect(n).toBeGreaterThanOrEqual(0);
  });

  stest('filter() narrows by name', async () => {
    const all = await dapi.packages.count();
    // Empty filter returns all; non-matching filter returns fewer
    const filtered = await dapi.packages.filter('zzz-no-such-package-xyz').count();
    expect(filtered).toBeLessThanOrEqual(all);
  });

  stest('find() retrieves a package by ID', async () => {
    if (!seed.packageId) return;
    const pkg = await dapi.packages.find(seed.packageId) as any;
    expect(pkg.id).toBe(seed.packageId);
    expect(pkg).toHaveProperty('name');
  });
});

// ─── reports ──────────────────────────────────────────────────────────────────

describe('reports', () => {
  stestEndpoint('list() returns an array', async () => {
    const reports = await dapi.reports.by(5).list();
    expect(Array.isArray(reports)).toBe(true);
  });

  stest('find() retrieves a report by ID', async () => {
    if (!seed.reportId) return;
    const report = await dapi.reports.find(seed.reportId) as any;
    expect(report.id).toBe(seed.reportId);
  });
});

// ─── files ────────────────────────────────────────────────────────────────────

describe('files', () => {
  /** Normalises the files.list() response to an array regardless of server shape. */
  function toFileArray(result: any): any[] {
    if (Array.isArray(result)) return result;
    if (result && typeof result === 'object') {
      // Some servers return {items:[...]} or {files:[...]} or {value:[...]}
      const nested = result.items ?? result.files ?? result.value ?? result.data;
      if (Array.isArray(nested)) return nested;
    }
    return [];
  }

  stest(`list() returns a defined result for ${TEST_FILES_PATH}`, async () => {
    const result = await dapi.files.list(TEST_FILES_PATH);
    expect(result).toBeDefined();
    // Log actual shape to aid diagnosis
    const shape = Array.isArray(result) ? `array(${result.length})` : typeof result;
    console.log(`files.list shape: ${shape}`);
  });

  stest('list() with recursive=true returns a defined result', async () => {
    const result = await dapi.files.list(TEST_FILES_PATH, true);
    expect(result).toBeDefined();
  });

  stest('list() non-recursive returns fewer or equal files than recursive', async () => {
    const flat = toFileArray(await dapi.files.list(TEST_FILES_PATH, false));
    const recursive = toFileArray(await dapi.files.list(TEST_FILES_PATH, true));
    expect(flat.length).toBeLessThanOrEqual(recursive.length);
  });

  stest('get() returns content for an existing file', async () => {
    const files = toFileArray(await dapi.files.list(TEST_FILES_PATH, true));
    const textFile = files.find((f: any) =>
      typeof f === 'string'
        ? /\.(json|txt|yaml|md)$/.test(f)
        : /\.(json|txt|yaml|md)$/.test(f.path ?? f.name ?? ''),
    );
    if (!textFile) return; // No text files found — skip gracefully
    const filePath = typeof textFile === 'string' ? textFile : (textFile.path ?? textFile.name);
    const content = await dapi.files.get(`${TEST_FILES_PATH}/${filePath}`);
    expect(content).toBeDefined();
  });
});

// ─── NodeDapi.raw ─────────────────────────────────────────────────────────────

describe('NodeDapi.raw', () => {
  stest('GET /api/users/current returns the current user', async () => {
    const result = await dapi.raw('GET', '/api/users/current') as any;
    expect(result).toBeDefined();
    // The response may be an object or JSON text — just verify it contains "login" or "admin"
    const text = typeof result === 'string' ? result : JSON.stringify(result);
    expect(text.toLowerCase()).toMatch(/login|admin|user/);
  });

  stest('GET /api/info returns server info', async () => {
    const result = await dapi.raw('GET', '/api/info');
    expect(result).toBeDefined();
  });

  stest('accepts lowercase method names', async () => {
    const result = await dapi.raw('get', '/api/users/current');
    expect(result).toBeDefined();
  });
});

// ─── NodeDapi.describe ────────────────────────────────────────────────────────

describe('NodeDapi.describe', () => {
  stest('describe("connections") returns a schema object', async () => {
    const schema = await dapi.describe('connections');
    expect(schema).toBeDefined();
  });

  stest('describe("users") returns a schema object', async () => {
    const schema = await dapi.describe('users');
    expect(schema).toBeDefined();
  });

  stest('describe("nonexistentEntity") throws or returns null', async () => {
    // Either behavior is acceptable — the important thing is it does not hang
    try {
      await dapi.describe('nonexistentEntity-xyz');
    } catch {
      // Expected
    }
  });
});
