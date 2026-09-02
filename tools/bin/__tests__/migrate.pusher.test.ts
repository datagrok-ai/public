import {describe, it, expect, beforeEach, afterEach} from 'vitest';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';
import {NodeDapi} from '../utils/node-dapi';
import {Bundle, BundleEntity, normalize, read, write} from '../utils/migrate/bundle';
import {plan, push} from '../utils/migrate/pusher';
import {loadCreds} from '../commands/server-migrate';

interface Call {method: string; path: string; body?: any}

const CONN_ID = 'bbbbbbbb-1111-2222-3333-444444444444';
const TABLE_ID = 'f2444470-6fc0-11f1-b307-bdcf313f748a';
const PROJECT_ID = 'efe0b1f0-6fc0-11f1-b275-83ec2160b5e9';
const GROUP_ID = 'dddddddd-1111-2222-3333-444444444444';
const QUERY_ID = 'aaaaaaaa-1111-2222-3333-444444444444';
const TWIN_ID = '99999999-1111-2222-3333-444444444444';

let dir: string;
beforeEach(() => { dir = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-push-')); });
afterEach(() => fs.rmSync(dir, {recursive: true, force: true}));

function makeDapi(responder: (method: string, path: string, body?: any) => any,
                  overrides: Record<string, any> = {}) {
  const calls: Call[] = [];
  const canned: Record<string, any> = {
    '/info/server': {Version: '1.27.9', Commit: 'abc'},
    '/users/current': {project: {name: 'Admin'}},
    '/packages': [],
    '/privileges/permissions': [],
    ...overrides,
  };
  const client: any = {
    baseUrl: 'http://h/api',
    async request(method: string, path: string, body?: any) {
      calls.push({method, path, body});
      for (const prefix of Object.keys(canned))
        if (path.startsWith(prefix))
          return typeof canned[prefix] === 'function' ? canned[prefix](method, path, body) : canned[prefix];
      return responder(method, path, body);
    },
    get(path: string) { return this.request('GET', path); },
    post(path: string, body?: any) { return this.request('POST', path, body); },
    del(path: string) { return this.request('DELETE', path); },
    async putBytes(path: string, bytes: Buffer) { calls.push({method: 'POST', path, body: bytes}); return ''; },
  };
  return {dapi: new NodeDapi(client), calls};
}

function bundleOf(entities: [string, BundleEntity][], bytes = new Map<string, Buffer>()): Bundle {
  const map = new Map(entities);
  const manifest = write(dir, map, {
    source: {url: 'http://src/api', version: '1.27.9', userNamespace: 'Admin:'}, args: [], packages: [],
  }, {}, bytes);
  const normalized = new Map<string, BundleEntity>();
  for (const e of manifest.order)
    normalized.set(e.id, {type: e.type, json: normalize(e.type, map.get(e.id)!.json), file: e.file});
  return {dir, manifest, entities: normalized, idmap: {}};
}

const script = (id: string, name: string, body: string) =>
  [id, {type: 'Script', json: {'#type': 'Script', id, name, namespace: 'Admin:', script: body}}] as [string, BundleEntity];

const notFound = () => { throw Object.assign(new Error('Not Found'), {apiError: {error: 'Not Found', errorCode: 404}}); };

describe('plan', () => {
  it('reports identical when the target hash matches and update when it differs', async () => {
    const bundle = bundleOf([script('s1', 'A', 'print(1)'), script('s2', 'B', 'print(2)')]);
    const {dapi} = makeDapi((_m, path) => {
      if (path.startsWith('/scripts/s1')) return {'#type': 'Script', id: 's1', name: 'A', namespace: 'Admin:', script: 'print(1)', updatedOn: 'now'};
      if (path.startsWith('/scripts/s2')) return {'#type': 'Script', id: 's2', name: 'B', namespace: 'Admin:', script: 'CHANGED'};
      return notFound();
    });
    const {rows, ops} = await plan(dapi, bundle, {onConflict: 'fail'});
    expect(rows.map((r) => r.action)).toEqual(['identical', 'update']);
    expect(rows[1].detail).toBe('script');
    expect(ops.map((o) => o.id)).toEqual(['s2']);
  });

  it('ignores the sync_id stamp and the server-derived columns of a table', async () => {
    const json = {'#type': 'TableInfo', id: TABLE_ID, name: 'Cereal', namespace: 'Admin:', metaParams: {}};
    const bundle = bundleOf([[TABLE_ID, {type: 'TableInfo', json}]]);
    const {dapi} = makeDapi((_m, path) => path.startsWith(`/tables/${TABLE_ID}`) ? {
      ...json, metaParams: {sync_id: TABLE_ID}, updatedOn: 'now',
      columns: [{'#type': 'ColumnInfo', name: 'name', type: 'string'}],
    } : notFound());
    const {rows, ops} = await plan(dapi, bundle, {onConflict: 'fail'});
    expect(rows.map((r) => r.action)).toEqual(['identical']);
    expect(ops).toEqual([]);
  });

  it('refuses a connection a hand-edited bundle should never carry', async () => {
    const conn = (id: string, name: string, namespace: string, parameters: any = {}): [string, BundleEntity] =>
      [id, {type: 'DataConnection', json: {'#type': 'DataConnection', id, name, namespace, parameters}}];
    const bundle = bundleOf([conn('c1', 'AppData', 'System:'), conn('c2', 'Home', 'Admin:'),
      conn('c3', 'Files', 'Spc:', {isProject: true})]);
    const {dapi, calls} = makeDapi((_m, path) => path.startsWith('/entities?') ? [] : notFound());
    const {rows, ops} = await plan(dapi, bundle, {onConflict: 'fail'});
    expect(rows.map((r) => `${r.action}:${r.reason}`).sort())
      .toEqual(['skip:personal_storage', 'skip:platform_connection', 'skip:space_files_connection']);
    expect(ops).toEqual([]);
    expect(calls.some((c) => c.method === 'POST')).toBe(false);
  });

  it('updates a project the target links less than the bundle, and keeps the target-only link', async () => {
    const other = '77777777-1111-2222-3333-444444444444';
    const project = {'#type': 'Project', id: PROJECT_ID, name: 'Dash', namespace: 'Admin:', isDashboard: true,
      relations: [{id: 'r1', isLink: false, entity: {'#type': 'EntityRecord', id: TABLE_ID}}]};
    const onTarget = {...project, relations: [{id: 'rx', isLink: false, entity: {'#type': 'EntityRecord', id: other}}]};
    const bundle = bundleOf([[PROJECT_ID, {type: 'Project', json: project}]]);
    const {dapi} = makeDapi((_m, path) => path.startsWith(`/projects/${PROJECT_ID}`) ? onTarget : notFound());
    const {rows} = await plan(dapi, bundle, {onConflict: 'fail'});
    expect(rows.map((r) => `${r.action}:${r.reason}`).slice(0, 2))
      .toEqual(['update:relations missing on target', 'info:relation_not_removed']);
    expect(rows[0].detail).toBe('relations');
    expect(rows[1].detail).toBe(other);
  });

  it('reports create when neither the id nor the name is taken', async () => {
    const bundle = bundleOf([script('s1', 'A', 'print(1)')]);
    const {dapi} = makeDapi((_m, path) => path.startsWith('/entities?') ? [] : notFound());
    const {rows} = await plan(dapi, bundle, {onConflict: 'fail'});
    expect(rows.map((r) => r.action)).toEqual(['create']);
  });

  it('looks the twin up under the pusher namespace when the source namespace was personal', async () => {
    const bundle = bundleOf([script('s1', 'A', 'print(1)')]);
    const {dapi, calls} = makeDapi((_m, path) => path.startsWith('/entities?') ? [] : notFound());
    await plan(dapi, bundle, {onConflict: 'fail'});
    expect(calls.find((c) => c.path.startsWith('/entities?'))!.path).toBe('/entities?namespace=Admin%3A&name=A');
  });
});

describe('conflict policies', () => {
  const twin = (_m: string, path: string) =>
    path.startsWith('/entities?') ? [{'#type': 'Script', id: 'other', name: 'A', namespace: 'Admin:'}] : notFound();

  it('fail aborts before any write, listing the conflicts', async () => {
    const bundle = bundleOf([script('s1', 'A', 'print(1)')]);
    const {dapi, calls} = makeDapi(twin);
    await expect(push(dapi, bundle, {onConflict: 'fail'}, () => {})).rejects.toThrow(/Name conflicts/);
    expect(calls.some((c) => c.method === 'POST')).toBe(false);
  });

  it('skip writes nothing for the conflicted entity', async () => {
    const bundle = bundleOf([script('s1', 'A', 'print(1)')]);
    const {dapi, calls} = makeDapi(twin);
    const result = await push(dapi, bundle, {onConflict: 'skip'}, () => {});
    expect(result.items.map((r) => r.action)).toEqual(['skip']);
    expect(calls.some((c) => c.method === 'POST')).toBe(false);
  });
});

describe('push', () => {
  it('issues no POST at all on a dry run', async () => {
    const bundle = bundleOf([script('s1', 'A', 'print(1)')]);
    const {dapi, calls} = makeDapi((_m, path) => path.startsWith('/entities?') ? [] : notFound());
    const result = await push(dapi, bundle, {dryRun: true, onConflict: 'fail'}, () => {});
    expect(result.status).toBe('dry-run');
    expect(calls.some((c) => c.method === 'POST')).toBe(false);
  });

  it('pushes table bytes before the table metadata', async () => {
    const table: [string, BundleEntity] = [TABLE_ID, {type: 'TableInfo', json: {'#type': 'TableInfo', id: TABLE_ID, name: 'Cereal', namespace: 'Admin:'}}];
    const bundle = bundleOf([table], new Map([[TABLE_ID, Buffer.from([1, 2, 3])]]));
    const saved: any = {'#type': 'TableInfo', id: TABLE_ID, name: 'Cereal', namespace: 'Admin:'};
    let exists = false;
    const {dapi, calls} = makeDapi((method, path, body) => {
      if (path.startsWith('/entities?')) return [];
      if (method === 'POST') { exists = true; return body; }
      return exists ? saved : notFound();
    });
    await push(dapi, bundle, {onConflict: 'fail'}, () => {});
    const posts = calls.filter((c) => c.method === 'POST').map((c) => c.path);
    expect(posts).toEqual([`/tables/data?id=${TABLE_ID}`, '/tables']);
  });

  it('saves project relations after every entity, and retries a deadlock', async () => {
    const project: [string, BundleEntity] = [PROJECT_ID, {type: 'Project', json: {
      '#type': 'Project', id: PROJECT_ID, name: 'Dash', namespace: 'Admin:', isDashboard: true,
      relations: [{id: 'r1', isLink: false, entity: {'#type': 'TableInfo', id: TABLE_ID}}],
    }}];
    const table: [string, BundleEntity] = [TABLE_ID, {type: 'TableInfo', json: {'#type': 'TableInfo', id: TABLE_ID, name: 'Cereal', namespace: 'Admin:'}}];
    const bundle = bundleOf([project, table]);
    const rows: Record<string, any> = {};
    let deadlocks = 0;
    const {dapi, calls} = makeDapi((method, path, body) => {
      if (path.startsWith('/entities?')) return [];
      if (method === 'POST' && path.startsWith('/projects?saveRelations=true') && deadlocks++ < 1)
        throw new Error('PostgreSQL error 40P01: deadlock detected');
      if (method === 'POST') {
        const stored = {...body};
        if (!path.includes('saveRelations=true')) delete stored.relations;
        rows[body.id] = stored;
        return stored;
      }
      const id = path.split('?')[0].split('/').pop()!;
      return rows[id] ?? notFound();
    });
    const result = await push(dapi, bundle, {onConflict: 'fail'}, () => {});
    expect(result.counts).toMatchObject({create: 2});
    const posts = calls.filter((c) => c.method === 'POST').map((c) => c.path);
    expect(posts).toEqual(['/projects', '/tables', '/projects?saveRelations=true', '/projects?saveRelations=true']);
    expect(calls.find((c) => c.path === '/projects?saveRelations=true')!.body.relations)
      .toEqual([{id: 'r1', entity: {'#type': 'EntityRecord', id: TABLE_ID}, isLink: false}]);
  });

  it('keeps a relation the target has and the bundle does not', async () => {
    const other = '77777777-1111-2222-3333-444444444444';
    const project: [string, BundleEntity] = [PROJECT_ID, {type: 'Project', json: {
      '#type': 'Project', id: PROJECT_ID, name: 'Dash', namespace: 'Admin:', isDashboard: true, friendlyName: 'new',
      relations: [{id: 'r1', isLink: false, entity: {'#type': 'EntityRecord', id: TABLE_ID}}],
    }}];
    const table: [string, BundleEntity] = [TABLE_ID, {type: 'TableInfo', json: {'#type': 'TableInfo', id: TABLE_ID, name: 'Cereal', namespace: 'Admin:'}}];
    const onTarget = {'#type': 'Project', id: PROJECT_ID, name: 'Dash', namespace: 'Admin:', isDashboard: true, friendlyName: 'old',
      relations: [{id: 'rx', isLink: false, entity: {'#type': 'EntityRecord', id: other}}]};
    const {dapi, calls} = makeDapi((method, path, body) => {
      if (path.startsWith('/entities?')) return [];
      if (method === 'POST') return body;
      if (path.startsWith(`/projects/${PROJECT_ID}`)) return onTarget;
      return notFound();
    });
    await push(dapi, bundleOf([project, table]), {onConflict: 'fail'}, () => {});
    const saved = calls.find((c) => c.path === '/projects?saveRelations=true')!.body.relations;
    expect(saved.map((r: any) => r.entity.id).sort()).toEqual([other, TABLE_ID].sort());
    expect(saved.find((r: any) => r.entity.id === other).id).toBe('rx');
  });

  it('does not re-save relations of a project the target already matches', async () => {
    const relations = [{id: 'r1', isLink: false, entity: {'#type': 'EntityRecord', id: TABLE_ID}}];
    const project = {'#type': 'Project', id: PROJECT_ID, name: 'Dash', namespace: 'Admin:', isDashboard: true, relations};
    const bundle = bundleOf([[PROJECT_ID, {type: 'Project', json: project}]]);
    const {dapi, calls} = makeDapi((_m, path) => path.startsWith('/projects/') ? {...project, updatedOn: 'now'} : notFound());
    const result = await push(dapi, bundle, {onConflict: 'fail'}, () => {});
    expect(result.items[0].action).toBe('identical');
    expect(calls.some((c) => c.method === 'POST')).toBe(false);
  });

  it('stamps sync_id and reports a connection whose password was stripped', async () => {
    const conn: [string, BundleEntity] = [CONN_ID, {type: 'DataConnection', json: {
      '#type': 'DataConnection', id: CONN_ID, name: 'NW', namespace: 'Admin:', metaParams: {},
      parameters: {db: 'northwind', password: '_____________'},
    }}];
    const bundle = bundleOf([conn]);
    const {dapi, calls} = makeDapi((method, path, body) => {
      if (path.startsWith('/entities?')) return [];
      if (method === 'POST') return body;
      return calls.some((c) => c.method === 'POST') ? {'#type': 'DataConnection', id: CONN_ID, name: 'NW'} : notFound();
    });
    const result = await push(dapi, bundle, {onConflict: 'fail'}, () => {});
    const posted = calls.find((c) => c.method === 'POST')!.body;
    expect(posted.parameters).toEqual({db: 'northwind'});
    expect(posted.metaParams.sync_id).toBe(CONN_ID);
    expect(posted._credentials).toBeUndefined();
    expect(result.items.map((r) => r.action)).toEqual(['create', 'needs-credentials']);
    expect(result.items[1].detail).toBe('password');
  });

  it('announces the missing credentials of a connection already on the target', async () => {
    const json = {'#type': 'DataConnection', id: CONN_ID, name: 'NW', namespace: 'Admin:', parameters: {db: 'northwind'}};
    const bundle = bundleOf([[CONN_ID, {type: 'DataConnection', json}]]);
    const {dapi} = makeDapi((_m, path) =>
      path.startsWith(`/connectors/connections/${CONN_ID}`) ? {...json, parameters: {db: 'other'}} : notFound());
    const result = await push(dapi, bundle, {dryRun: true, onConflict: 'fail'}, () => {});
    expect(result.items.map((r) => r.action)).toEqual(['update', 'needs-credentials']);
    expect(result.items[1].detail).toBe('credentials are not migrated');
  });

  it('skips a connection whose every parameter was masked', async () => {
    const conn: [string, BundleEntity] = [CONN_ID, {type: 'DataConnection', json: {
      '#type': 'DataConnection', id: CONN_ID, name: 'NW', namespace: 'Admin:', parameters: {password: '_____________'},
    }}];
    const bundle = bundleOf([conn]);
    const {dapi, calls} = makeDapi((_m, path) => path.startsWith('/entities?') ? [] : notFound());
    const result = await push(dapi, bundle, {onConflict: 'fail'}, () => {});
    expect(result.items.map((r) => r.action)).toEqual(['skip']);
    expect(calls.some((c) => c.method === 'POST')).toBe(false);
  });

  it('warns on a major.minor version mismatch and fails a save the target did not keep', async () => {
    const bundle = bundleOf([script('s1', 'A', 'print(1)')]);
    bundle.manifest.source.version = '1.28.0';
    const {dapi} = makeDapi((method, path, body) => {
      if (path.startsWith('/entities?')) return [];
      if (method === 'POST') return body;
      return notFound();
    });
    const result = await push(dapi, bundle, {onConflict: 'fail'}, () => {});
    expect(result.items[0]).toMatchObject({action: 'warn', reason: 'version_mismatch'});
    expect(result.items[1]).toMatchObject({action: 'failed', reason: 'Save reported success but not on target'});
    expect(result.status).toBe('failed');
  });
});

/** Stores what it is POSTed and serves it back by id, like the server does. */
function storingDapi(twins: any[] = []) {
  const saved: Record<string, any> = {};
  return {
    ...makeDapi((method, path, body) => {
      if (path.startsWith('/entities?')) {
        const name = decodeURIComponent(path.split('name=')[1] ?? '');
        return twins.filter((t) => t.name === name);
      }
      if (method === 'POST') {
        // the server touches relations only when asked to
        const stored = {...body};
        if (!path.includes('saveRelations=true')) delete stored.relations;
        saved[body.id] = stored;
        return stored;
      }
      const id = path.split('?')[0].split('/').pop()!;
      return saved[id] ?? notFound();
    }),
    saved,
  };
}

describe('conflict policy adopt', () => {
  const conn: [string, BundleEntity] = [CONN_ID, {type: 'DataConnection', json: {
    '#type': 'DataConnection', id: CONN_ID, name: 'NW', namespace: 'Admin:', parameters: {db: 'northwind'},
  }}];
  const query: [string, BundleEntity] = [QUERY_ID, {type: 'DataQuery', json: {
    '#type': 'DataQuery', id: QUERY_ID, name: 'Q', namespace: 'Admin:', connection: {id: CONN_ID},
  }}];
  const twin = {'#type': 'DataConnection', id: TWIN_ID, name: 'NW', namespace: 'Admin:'};

  it('saves into the twin, rewrites later references, and records idmap.json', async () => {
    const bundle = bundleOf([conn, query]);
    const {dapi, calls} = storingDapi([twin]);
    const result = await push(dapi, bundle, {onConflict: 'adopt'}, () => {});

    expect(result.items[0].reason).toBe(`adopted ${TWIN_ID}`);
    const posts = calls.filter((c) => c.method === 'POST');
    expect(posts[0].body.id).toBe(TWIN_ID);
    expect(posts[1].body.connection).toEqual({id: TWIN_ID});
    expect(JSON.parse(fs.readFileSync(path.join(dir, 'idmap.json'), 'utf8'))[CONN_ID]).toBe(TWIN_ID);
  });

  it('re-pushes the adopted bundle without writing anything', async () => {
    bundleOf([conn, query]);
    const {dapi, calls} = storingDapi([twin]);
    await push(dapi, read(dir), {onConflict: 'adopt'}, () => {});
    const before = calls.length;
    const result = await push(dapi, read(dir), {onConflict: 'adopt'}, () => {});
    expect(result.counts).toEqual({identical: 2});
    expect(calls.slice(before).some((c) => c.method === 'POST')).toBe(false);
  });

  it('resolves every twin before rewriting, so a forward reference lands on the adopted id', async () => {
    const SUB_ID = '55555555-1111-2222-3333-444444444444';
    const SUB_TWIN = '66666666-1111-2222-3333-444444444444';
    const root: [string, BundleEntity] = [PROJECT_ID, {type: 'Project', json: {
      '#type': 'Project', id: PROJECT_ID, name: 'Root', namespace: '', isDashboard: false, isRoot: true,
      relations: [{id: 'r1', isLink: false, entity: {'#type': 'EntityRecord', id: SUB_ID}}],
    }}];
    const sub: [string, BundleEntity] = [SUB_ID, {type: 'Project', json: {
      '#type': 'Project', id: SUB_ID, name: 'Sub', namespace: 'Root:', isDashboard: false,
    }}];
    const bundle = bundleOf([root, sub]);
    const {dapi, calls} = storingDapi([
      {'#type': 'Project', id: TWIN_ID, name: 'Root', namespace: ''},
      {'#type': 'Project', id: SUB_TWIN, name: 'Sub', namespace: 'Root:'},
    ]);
    await push(dapi, bundle, {onConflict: 'adopt'}, () => {});
    const rootSave = calls.filter((c) => c.method === 'POST' && c.body?.id === TWIN_ID);
    expect(rootSave[0].body.relations[0].entity.id).toBe(SUB_TWIN);
    expect(calls.some((c) => c.path === '/projects?saveRelations=true' && c.body.id === TWIN_ID)).toBe(true);
  });

  it('mints fresh ids for the nested rows of an adopted entity', async () => {
    const relation = {id: 'cccccccc-1111-2222-3333-444444444444', isLink: false, entity: {'#type': 'EntityRecord', id: TABLE_ID}};
    const project: [string, BundleEntity] = [PROJECT_ID, {type: 'Project', json: {
      '#type': 'Project', id: PROJECT_ID, name: 'Dash', namespace: 'Admin:', isDashboard: true, relations: [relation],
    }}];
    const table: [string, BundleEntity] = [TABLE_ID, {type: 'TableInfo', json: {'#type': 'TableInfo', id: TABLE_ID, name: 'Cereal', namespace: 'Admin:'}}];
    const bundle = bundleOf([project, table]);
    const {dapi, calls} = storingDapi([{'#type': 'Project', id: TWIN_ID, name: 'Dash', namespace: 'Admin:'}]);
    await push(dapi, bundle, {onConflict: 'adopt'}, () => {});
    const relations = calls.find((c) => c.path === '/projects?saveRelations=true')!.body.relations;
    expect(relations[0].entity.id).toBe(TABLE_ID);
    expect(relations[0].id).not.toBe(relation.id);
  });
});

describe('conflict policy duplicate', () => {
  it('creates the bundle entity under its own id, next to the twin', async () => {
    const bundle = bundleOf([script('s1', 'A', 'print(1)')]);
    const {dapi, calls} = storingDapi([{'#type': 'Script', id: TWIN_ID, name: 'A', namespace: 'Admin:'}]);
    const result = await push(dapi, bundle, {onConflict: 'duplicate'}, () => {});
    expect(result.items.map((r) => r.action)).toEqual(['create']);
    expect(calls.find((c) => c.method === 'POST')!.body.id).toBe('s1');
  });
});

describe('tags', () => {
  const tagged: [string, BundleEntity] = [QUERY_ID, {type: 'Script', json: {
    '#type': 'Script', id: QUERY_ID, name: 'A', namespace: 'Admin:', script: 'print(1)',
    entityTags: [{id: 'tag-row', tag: 'migrate-test', entity: {id: QUERY_ID}}],
  }}];

  it('tags the target only when the tag is missing there', async () => {
    const bundle = bundleOf([tagged]);
    expect(bundle.entities.get(QUERY_ID)!.json._tags).toEqual(['migrate-test']);
    const {dapi, calls} = storingDapi();
    await push(dapi, bundle, {onConflict: 'fail'}, () => {});
    expect(calls.filter((c) => c.method === 'POST').map((c) => c.path))
      .toEqual(['/scripts', '/entities/tag?tag=migrate-test']);
    expect(calls.find((c) => c.path.startsWith('/entities/tag'))!.body).toEqual([QUERY_ID]);
  });

  it('writes no tag POST when the target already carries it', async () => {
    const bundle = bundleOf([tagged]);
    const {dapi, calls, saved} = storingDapi();
    saved[QUERY_ID] = {...bundle.entities.get(QUERY_ID)!.json, entityTags: [{id: 'other-row', tag: 'migrate-test'}]};
    const result = await push(dapi, bundle, {onConflict: 'fail'}, () => {});
    expect(result.items.map((r) => r.action)).toEqual(['identical']);
    expect(calls.some((c) => c.method === 'POST')).toBe(false);
  });
});

describe('save routes', () => {
  it('posts a group to the trailing-slash route, the only one the server exposes', async () => {
    const bundle = bundleOf([[GROUP_ID, {type: 'UserGroup', json: {
      '#type': 'UserGroup', id: GROUP_ID, name: 'Chemists', friendlyName: 'Chemists', parents: [], children: [],
    }}]]);
    const {dapi, calls} = storingDapi();
    const result = await push(dapi, bundle, {onConflict: 'fail'}, () => {});
    expect(result.items.map((r) => r.action)).toEqual(['create']);
    expect(calls.filter((c) => c.method === 'POST').map((c) => c.path)).toEqual(['/groups/']);
  });
});

describe('memberships and grants', () => {
  const group = (over: any = {}) => ({
    '#type': 'UserGroup', id: GROUP_ID, name: 'Chemists', friendlyName: 'Chemists',
    parents: [], children: [], ...over,
  });

  function groupBundle(members: any[], children: any[] = []) {
    const bundle = bundleOf([[GROUP_ID, {type: 'UserGroup', json: {...group(), _members: members}}]]);
    const saved = group({children});
    const {dapi, calls} = makeDapi((method, path, body) => {
      if (path.startsWith('/public/v1/groups/lookup')) {
        const query = decodeURIComponent(path.split('query=')[1]);
        return query === 'alice' ? [{id: 'alice-group', friendlyName: 'alice', personal: true}] : [];
      }
      if (path.startsWith('/public/v1/groups')) return method === 'POST' ? body : saved;
      if (path.startsWith('/entities?')) return [];
      if (path.startsWith(`/groups/${GROUP_ID}`)) return saved;
      if (method === 'POST') return body;
      return notFound();
    });
    return {bundle, dapi, calls};
  }

  it('replays a member by login and reports one the target does not have', async () => {
    const {bundle, dapi, calls} = groupBundle([
      {kind: 'user', login: 'alice', isAdmin: false},
      {kind: 'user', login: 'ghost', isAdmin: false},
    ]);
    const result = await push(dapi, bundle, {onConflict: 'fail'}, () => {});
    expect(calls.some((c) => c.path === '/public/v1/groups/lookup?query=alice')).toBe(true);
    expect(calls.filter((c) => c.method === 'POST').map((c) => c.path))
      .toEqual(['/public/v1/groups?saveRelations=true']);
    expect(result.items.find((r) => r.reason === 'member_not_found')!.detail).toMatch(/^ghost:/);
    expect(result.items[0].detail).toContain('members: 1 matched, 1 not on remote');
  });

  it('writes no membership POST when the target already has the member', async () => {
    const {bundle, dapi, calls} = groupBundle([{kind: 'user', login: 'alice', isAdmin: false}],
      [{id: 'rel1', parent: {id: GROUP_ID}, child: {id: 'alice-group'}, isAdmin: false}]);
    await push(dapi, bundle, {onConflict: 'fail'}, () => {});
    expect(calls.some((c) => c.method === 'POST')).toBe(false);
  });

  function grantBundle(grants: any[], existing: any[] = []) {
    const json = {'#type': 'Project', id: PROJECT_ID, name: 'Dash', namespace: 'Admin:', isDashboard: true};
    const bundle = bundleOf([[PROJECT_ID, {type: 'Project', json: {...json, _grants: grants}}]]);
    let saved = false;
    const {dapi, calls} = makeDapi((method, path, body) => {
      if (path.startsWith('/public/v1/groups/lookup'))
        return decodeURIComponent(path.split('query=')[1]) === 'Chemists' ? [{id: GROUP_ID, friendlyName: 'Chemists'}] : [];
      if (path.startsWith(`/groups/${GROUP_ID}`)) return group();
      if (path.startsWith('/entities?')) return [];
      if (method === 'POST') { saved = true; return body; }
      return saved ? json : notFound();
    }, {'/privileges/permissions': () => existing});
    return {bundle, dapi, calls};
  }

  it('shares only the pairs the target is missing', async () => {
    const {bundle, dapi, calls} = grantBundle([{group: 'Chemists', permission: 'View'}]);
    const result = await push(dapi, bundle, {onConflict: 'fail'}, () => {});
    expect(calls.map((c) => c.path)).toContain(`/public/v1/entities/${PROJECT_ID}/shares?groups=Chemists&access=View`);
    expect(result.items.some((r) => r.reason === 'not_visible')).toBe(false);
  });

  it('writes no share when the grant is already there', async () => {
    const existing = [{entityId: PROJECT_ID, userGroup: {id: GROUP_ID}, permission: {name: 'View'}}];
    const {bundle, dapi, calls} = grantBundle([{group: 'Chemists', permission: 'View'}], existing);
    await push(dapi, bundle, {onConflict: 'fail'}, () => {});
    expect(calls.some((c) => c.path.includes('/shares?'))).toBe(false);
  });

  it('reports an unsupported permission, an unknown group, and invisible content', async () => {
    const {bundle, dapi, calls} = grantBundle([
      {group: 'Chemists', permission: 'Delete'},
      {group: 'Ghosts', permission: 'View'},
    ]);
    const result = await push(dapi, bundle, {onConflict: 'fail'}, () => {});
    expect(calls.some((c) => c.path.includes('/shares?'))).toBe(false);
    expect(result.items.map((r) => r.reason)).toContain('unsupported_grant');
    expect(result.items.map((r) => r.reason)).toContain('group_not_found');
    expect(result.items.some((r) => r.reason === 'not_visible')).toBe(true);
  });
});

describe('a skipped entity takes its dependants with it', () => {
  const conn: [string, BundleEntity] = [CONN_ID, {type: 'DataConnection', json: {
    '#type': 'DataConnection', id: CONN_ID, name: 'NW', namespace: 'Admin:', parameters: {db: 'nw'},
  }}];
  const table: [string, BundleEntity] = [TABLE_ID, {type: 'TableInfo', json: {
    '#type': 'TableInfo', id: TABLE_ID, name: 'Cereal', namespace: 'Admin:',
  }}];
  const group: [string, BundleEntity] = [GROUP_ID, {type: 'UserGroup', json: {
    '#type': 'UserGroup', id: GROUP_ID, name: 'Chemists', friendlyName: 'Chemists', namespace: 'Admin:',
    parents: [], children: [],
  }}];

  /** Everything is new on the target except the entity named `taken`, whose name another id holds. */
  function pushWith(entities: [string, BundleEntity][], taken: string, type: string) {
    const bundle = bundleOf(entities);
    const {dapi, calls} = makeDapi((_m, path) => {
      if (path.startsWith('/entities?'))
        return decodeURIComponent(path.split('name=')[1] ?? '') === taken ? [{'#type': type, id: TWIN_ID, name: taken, namespace: 'Admin:'}] : [];
      return notFound();
    });
    return {bundle, dapi, calls};
  }

  const rowOf = (result: any, type: string) => result.items.find((r: any) => r.entityType === type);

  it('fails a query whose connection was skipped, without calling the server', async () => {
    const query: [string, BundleEntity] = [QUERY_ID, {type: 'DataQuery', json: {
      '#type': 'DataQuery', id: QUERY_ID, name: 'Q', namespace: 'Admin:', connection: {id: CONN_ID},
    }}];
    const {bundle, dapi, calls} = pushWith([conn, query], 'NW', 'DataConnection');
    const result = await push(dapi, bundle, {onConflict: 'skip'}, () => {});
    expect(rowOf(result, 'DataConnection').action).toBe('skip');
    expect(rowOf(result, 'DataQuery')).toMatchObject({action: 'failed', reason: 'dependency_skipped', detail: 'Admin:NW'});
    expect(calls.some((c) => c.method === 'POST')).toBe(false);
  });

  it('fails a view whose table was skipped', async () => {
    const view: [string, BundleEntity] = ['v1', {type: 'ViewInfo', json: {
      '#type': 'ViewInfo', id: 'v1', name: 'Cereal', namespace: 'Admin:', table: {id: TABLE_ID},
    }}];
    const {bundle, dapi} = pushWith([table, view], 'Cereal', 'TableInfo');
    const result = await push(dapi, bundle, {onConflict: 'skip'}, () => {});
    expect(rowOf(result, 'TableInfo').action).toBe('skip');
    expect(rowOf(result, 'ViewInfo')).toMatchObject({action: 'failed', reason: 'dependency_skipped'});
  });

  it('fails a project whose related child was skipped', async () => {
    const project: [string, BundleEntity] = [PROJECT_ID, {type: 'Project', json: {
      '#type': 'Project', id: PROJECT_ID, name: 'Dash', namespace: 'Admin:', isDashboard: true,
      relations: [{id: 'r1', isLink: false, entity: {'#type': 'EntityRecord', id: TABLE_ID}}],
    }}];
    const {bundle, dapi, calls} = pushWith([project, table], 'Cereal', 'TableInfo');
    const result = await push(dapi, bundle, {onConflict: 'skip'}, () => {});
    expect(rowOf(result, 'Project')).toMatchObject({action: 'failed', reason: 'dependency_skipped'});
    expect(calls.some((c) => c.method === 'POST')).toBe(false);
  });

  it('fails a project whose grant-holding group was skipped', async () => {
    const project: [string, BundleEntity] = [PROJECT_ID, {type: 'Project', json: {
      '#type': 'Project', id: PROJECT_ID, name: 'Dash', namespace: 'Admin:', isDashboard: true,
      _grants: [{group: 'Chemists', permission: 'View'}],
    }}];
    const {bundle, dapi} = pushWith([group, project], 'Chemists', 'UserGroup');
    const result = await push(dapi, bundle, {onConflict: 'skip'}, () => {});
    expect(rowOf(result, 'UserGroup').action).toBe('skip');
    expect(rowOf(result, 'Project')).toMatchObject({action: 'failed', reason: 'dependency_skipped', detail: 'Admin:Chemists'});
  });

  it('fails a file whose connection was skipped', async () => {
    const file: [string, BundleEntity] = ['f1', {type: 'FileInfo', json: {
      '#type': 'FileInfo', id: 'f1', name: 'a.csv', isFile: true, connection: {id: CONN_ID},
    }}];
    const {bundle, dapi} = pushWith([conn, file], 'NW', 'DataConnection');
    const result = await push(dapi, bundle, {onConflict: 'skip'}, () => {});
    expect(rowOf(result, 'FileInfo')).toMatchObject({action: 'failed', reason: 'dependency_skipped'});
  });
});

describe('files, spaces and models', () => {
  const file = (over: any = {}): [string, BundleEntity] => ['f1', {type: 'FileInfo', json: {
    '#type': 'FileInfo', id: 'f1', name: 'a.csv', friendlyName: 'a.csv', path: 'a.csv', isFile: true, ...over,
  }}];

  it('writes the metadata first and the bytes under the id the target answered with', async () => {
    const bundle = bundleOf([file()], new Map([['f1', Buffer.from([1, 2, 3])]]));
    const saved: any = {'#type': 'FileInfo', id: 'existing-id', name: 'a.csv'};
    const {dapi, calls} = makeDapi((method, path) => {
      if (path.startsWith('/entities?')) return [];
      if (method === 'POST' && path === '/files') return saved;
      return path.startsWith('/files/existing-id') ? saved : notFound();
    });
    const result = await push(dapi, bundle, {onConflict: 'fail'}, () => {});
    expect(result.items.map((r) => r.action)).toEqual(['create']);
    expect(calls.filter((c) => c.method === 'POST').map((c) => c.path)).toEqual(['/files', '/files/data/existing-id']);
    expect(JSON.parse(fs.readFileSync(path.join(dir, 'idmap.json'), 'utf8'))['f1']).toBe('existing-id');
  });

  it('reports that a model travels without its trained blob', async () => {
    const model: [string, BundleEntity] = ['m1', {type: 'PredictiveModelInfo', json: {
      '#type': 'PredictiveModelInfo', id: 'm1', name: 'Model', namespace: 'Admin:', trainedOn: {id: TABLE_ID},
    }}];
    const {dapi, calls} = storingDapi();
    const result = await push(dapi, bundleOf([model]), {onConflict: 'fail'}, () => {});
    expect(calls.find((c) => c.method === 'POST')!.body.trainedOn).toEqual({id: TABLE_ID});
    expect(result.items.map((r) => `${r.action}:${r.reason}`)).toContain('info:model_blob_skipped');
  });
});

describe('a masked connection that is already on the target blocks nothing', () => {
  it('keeps pushing the query of a connection skipped for having no parameters left', async () => {
    const conn = {'#type': 'DataConnection', id: CONN_ID, name: 'NW', namespace: 'Admin:',
      parameters: {password: '_____________'}};
    const query = {'#type': 'DataQuery', id: QUERY_ID, name: 'Q', namespace: 'Admin:', connection: {id: CONN_ID}};
    const bundle = bundleOf([
      [CONN_ID, {type: 'DataConnection', json: conn}],
      [QUERY_ID, {type: 'DataQuery', json: query}],
    ]);
    const {dapi, calls} = makeDapi((method, path, body) => {
      if (path.startsWith(`/connectors/connections/${CONN_ID}`))
        return {...conn, parameters: {password: '_____________', db: 'nw'}};
      if (path.startsWith('/entities?')) return [];
      if (method === 'POST') return body;
      return calls.some((c) => c.method === 'POST') ? query : notFound();
    });
    const result = await push(dapi, bundle, {onConflict: 'fail'}, () => {});
    expect(result.items.map((r) => r.action)).toEqual(['skip', 'create']);
    expect(calls.filter((c) => c.method === 'POST').map((c) => c.path)).toEqual(['/connectors/queries']);
  });
});

describe('--creds', () => {
  const masked: [string, BundleEntity] = [CONN_ID, {type: 'DataConnection', json: {
    '#type': 'DataConnection', id: CONN_ID, name: 'NW', namespace: 'Admin:',
    parameters: {db: 'northwind', password: '_____________'},
  }}];

  it('merges the target-side secret into the payload and drops needs-credentials', async () => {
    const bundle = bundleOf([masked]);
    const {dapi, calls} = storingDapi();
    const result = await push(dapi, bundle, {onConflict: 'fail', creds: {'Admin:NW': {password: 'secret'}}}, () => {});
    expect(calls.find((c) => c.method === 'POST')!.body.parameters).toEqual({db: 'northwind', password: 'secret'});
    expect(result.items.map((r) => r.action)).toEqual(['create']);
    expect(JSON.parse(fs.readFileSync(path.join(dir, 'DataConnection/Admin.NW.json'), 'utf8')).parameters.password).toBeUndefined();
  });

  it('rotates the secret of a connection the target already matches', async () => {
    const conn = {'#type': 'DataConnection', id: CONN_ID, name: 'NW', namespace: 'Admin:', parameters: {db: 'northwind'}};
    const {dapi, calls} = makeDapi((method, path, body) => {
      if (path.startsWith(`/connectors/connections/${CONN_ID}`)) return {...conn, updatedOn: 'now'};
      if (path.startsWith('/entities?')) return [];
      return method === 'POST' ? body : notFound();
    });
    const bundle = bundleOf([[CONN_ID, {type: 'DataConnection', json: conn}]]);
    const result = await push(dapi, bundle, {onConflict: 'fail', creds: {'Admin:NW': {password: 'secret'}}}, () => {});
    expect(result.items.map((r) => `${r.action}:${r.reason}`)).toEqual(['update:credentials']);
    expect(calls.find((c) => c.method === 'POST')!.body.parameters).toEqual({db: 'northwind', password: 'secret'});
  });

  it('reports needs-credentials for a connection the file does not cover', async () => {
    const {dapi} = storingDapi();
    const result = await push(dapi, bundleOf([masked]), {onConflict: 'fail', creds: {'Admin:Other': {password: 'secret'}}}, () => {});
    expect(result.items.map((r) => r.action)).toEqual(['create', 'needs-credentials']);
  });

  it('accepts the key under the pusher namespace when the source namespace was personal', async () => {
    const source: [string, BundleEntity] = [CONN_ID, {type: 'DataConnection', json: {
      '#type': 'DataConnection', id: CONN_ID, name: 'NW', namespace: 'Askalkin:', parameters: {db: 'northwind'},
    }}];
    const map = new Map([source]);
    const manifest = write(dir, map, {source: {url: 'http://src/api', version: '1.27.9', userNamespace: 'Askalkin:'}, args: [], packages: []}, {replace: true});
    const entities = new Map<string, BundleEntity>();
    for (const e of manifest.order)
      entities.set(e.id, {type: e.type, json: normalize(e.type, map.get(e.id)!.json), file: e.file});
    const {dapi, calls} = storingDapi();
    const result = await push(dapi, {dir, manifest, entities, idmap: {}}, {onConflict: 'fail', creds: {'Admin:NW': {password: 'secret'}}}, () => {});
    expect(result.items.map((r) => r.action)).toEqual(['create']);
    expect(calls.find((c) => c.method === 'POST')!.body.parameters.password).toBe('secret');
  });
});

describe('loadCreds', () => {
  const write = (text: string) => {
    const file = path.join(dir, 'creds.yaml');
    fs.writeFileSync(file, text);
    return file;
  };

  it('substitutes the environment into the values', () => {
    process.env['MIGTEST_PASSWORD'] = 's3cret$1';
    expect(loadCreds(write('Admin:NW:\n  password: ${MIGTEST_PASSWORD}\n')))
      .toEqual({'Admin:NW': {password: 's3cret$1'}});
    delete process.env['MIGTEST_PASSWORD'];
  });

  it('names every missing variable instead of writing anything', () => {
    delete process.env['MIGTEST_ABSENT'];
    expect(() => loadCreds(write('Admin:NW:\n  password: ${MIGTEST_ABSENT}\n')))
      .toThrow(/cannot find environment variable "MIGTEST_ABSENT"/);
  });

  it('rejects a top-level value that is not a map, naming the key', () => {
    expect(() => loadCreds(write('Admin:NW: s3cret\n'))).toThrow(/"Admin:NW" must be a map/);
    expect(() => loadCreds(write('Admin:NW: [a, b]\n'))).toThrow(/"Admin:NW" must be a map/);
  });
});
