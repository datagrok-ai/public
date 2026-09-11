import {describe, it, expect, beforeEach, afterEach} from 'vitest';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';
import {NodeDapi} from '../utils/node-dapi';
import {Bundle, BundleEntity, normalize, read, write, writeShares} from '../utils/migrate/bundle';
import {plan, push} from '../utils/migrate/pusher';
import {loadCreds} from '../commands/server-migrate';

interface Call {method: string; path: string; body?: any}

const CONN_ID = 'bbbbbbbb-1111-2222-3333-444444444444';
const TABLE_ID = 'f2444470-6fc0-11f1-b307-bdcf313f748a';
const PROJECT_ID = 'efe0b1f0-6fc0-11f1-b275-83ec2160b5e9';
const GROUP_ID = 'dddddddd-1111-2222-3333-444444444444';
const QUERY_ID = 'aaaaaaaa-1111-2222-3333-444444444444';
const TWIN_ID = '99999999-1111-2222-3333-444444444444';
const VIEW_ID = '77777777-1111-2222-3333-444444444444';

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
      // Placement reads the rows on their own; a fixture only has to describe the project.
      const forProject = /^\/projects\/relations\?.*projectId=([^&]+)/.exec(path);
      if (forProject) {
        const project = await this.request('GET', `/projects/${forProject[1]}`).catch(() => null);
        return project?.relations ?? [];
      }
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

  it('takes the intended name back once the entity is placed', async () => {
    const project: [string, BundleEntity] = [PROJECT_ID, {type: 'Project', json: {
      '#type': 'Project', id: PROJECT_ID, name: 'Dash', namespace: 'Admin:', isDashboard: true,
      relations: [{id: 'r1', isLink: false, entity: {'#type': 'EntityRecord', id: TABLE_ID}}],
    }}];
    const table: [string, BundleEntity] = [TABLE_ID, {type: 'TableInfo', json: {'#type': 'TableInfo', id: TABLE_ID, name: 'Cereal', namespace: 'Admin:'}}];
    const stored: Record<string, any> = {};
    let placed = false;
    const {dapi} = makeDapi((method, path, body) => {
      if (path.startsWith('/entities?')) return [];
      if (path === '/projects?saveRelations=true') { placed = true; stored[body.id] = body; return body; }
      // Until the project is placed the name is taken, exactly as the server sees it.
      if (method === 'POST') {
        stored[body.id] = {...body, name: body.name === 'Dash' && !placed ? 'Dash_1' : body.name, relations: []};
        return stored[body.id];
      }
      const id = path.split('?')[0].split('/').pop()!;
      return stored[id] ?? notFound();
    });
    const result = await push(dapi, bundleOf([project, table]), {onConflict: 'fail'}, () => {});
    expect(stored[PROJECT_ID].name).toBe('Dash');
    const row = result.items.find((r) => ['renamed', 'name_restored'].includes(r.reason))!;
    expect(row).toMatchObject({action: 'info', reason: 'name_restored'});
  });

  it('places the entity again after taking its name back', async () => {
    // Saving an entity re-homes it into the pusher's own root, and the namespace is derived from
    // containment — so restoring a name after placement quietly undoes the placement.
    const project: [string, BundleEntity] = [PROJECT_ID, {type: 'Project', json: {
      '#type': 'Project', id: PROJECT_ID, name: 'Dash', namespace: 'Skalkin:', isDashboard: true,
      relations: [{id: 'r1', isLink: false, entity: {'#type': 'EntityRecord', id: TABLE_ID}}],
    }}];
    const table: [string, BundleEntity] = [TABLE_ID, {type: 'TableInfo', json: {
      '#type': 'TableInfo', id: TABLE_ID, name: 'Cereal', namespace: 'Skalkin:'}}];
    const stored: Record<string, any> = {};
    let placed = false;
    const order: string[] = [];
    const {dapi} = makeDapi((method, path, body) => {
      if (path.startsWith('/entities?')) return [];
      if (path === '/projects?saveRelations=true') {
        placed = true; order.push('place'); stored[body.id] = body; return body;
      }
      if (method === 'POST') {
        if (body?.name === 'Dash') order.push('save-project');
        // The entity save drops it back out of the project it was just placed in.
        stored[body.id] = {...body, name: body.name === 'Dash' && !placed ? 'Dash_1' : body.name, relations: []};
        return stored[body.id];
      }
      const id = path.split('?')[0].split('/').pop()!;
      return stored[id] ?? notFound();
    });
    await push(dapi, bundleOf([project, table]), {onConflict: 'fail'}, () => {});
    expect(order.lastIndexOf('place')).toBeGreaterThan(order.lastIndexOf('save-project'));
    expect(stored[PROJECT_ID].relations.map((r: any) => r.entity.id)).toContain(TABLE_ID);
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
      if (path.startsWith(`/tables/${TABLE_ID}`)) return table[1].json;
      return notFound();
    });
    await push(dapi, bundleOf([project, table]), {onConflict: 'fail'}, () => {});
    const saved = calls.find((c) => c.path === '/projects?saveRelations=true')!.body.relations;
    expect(saved.map((r: any) => r.entity.id).sort()).toEqual([other, TABLE_ID].sort());
    expect(saved.find((r: any) => r.entity.id === other).id).toBe('rx');
  });

  it('drops the one relation the target refuses and writes the rest', async () => {
    const refused = '88888888-1111-2222-3333-444444444444';
    const project: [string, BundleEntity] = [PROJECT_ID, {type: 'Project', json: {
      '#type': 'Project', id: PROJECT_ID, name: 'Dash', namespace: 'Admin:', isDashboard: true,
      relations: [{id: 'r1', isLink: false, entity: {'#type': 'EntityRecord', id: TABLE_ID}},
        {id: 'r2', isLink: false, entity: {'#type': 'EntityRecord', id: refused}}],
    }}];
    const table: [string, BundleEntity] = [TABLE_ID, {type: 'TableInfo', json: {'#type': 'TableInfo', id: TABLE_ID, name: 'Cereal', namespace: 'Admin:'}}];
    const other: [string, BundleEntity] = [refused, {type: 'TableInfo', json: {'#type': 'TableInfo', id: refused, name: 'Broken', namespace: 'Admin:'}}];
    const stored: Record<string, any> = {};
    const {dapi, calls} = makeDapi((method, path, body) => {
      if (path.startsWith('/entities?')) return [];
      if (path === '/projects?saveRelations=true') {
        if (body.relations.some((r: any) => r.entity.id === refused))
          throw new Error(`Unable to add entity ${refused} to the project`);
        return body;
      }
      if (method === 'POST') { stored[body.id] = {...body, relations: []}; return body; }
      const id = path.split('?')[0].split('/').pop()!;
      return stored[id] ?? notFound();
    });
    const result = await push(dapi, bundleOf([project, table, other]), {onConflict: 'fail'}, () => {});
    const written = calls.filter((c) => c.path === '/projects?saveRelations=true');
    expect(written[written.length - 1].body.relations.map((r: any) => r.entity.id)).toEqual([TABLE_ID]);
    expect(result.items.find((r) => r.reason === 'relations_refused')!.detail).toContain(refused);
    expect(result.items.some((r) => r.action === 'failed')).toBe(false);
  });

  it('leaves a platform connection where the target keeps it', async () => {
    const demoFiles = '99999999-1111-2222-3333-444444444444';
    const project: [string, BundleEntity] = [PROJECT_ID, {type: 'Project', json: {
      '#type': 'Project', id: PROJECT_ID, name: 'Dash', namespace: 'Admin:', isDashboard: true,
      relations: [{id: 'r1', isLink: false, entity: {'#type': 'EntityRecord', id: TABLE_ID}},
        {id: 'r2', isLink: false, entity: {'#type': 'EntityRecord', id: demoFiles}}],
    }}];
    const table: [string, BundleEntity] = [TABLE_ID, {type: 'TableInfo', json: {'#type': 'TableInfo', id: TABLE_ID, name: 'Cereal', namespace: 'Admin:'}}];
    const stored: Record<string, any> = {};
    const {dapi, calls} = makeDapi((method, path, body) => {
      if (path.startsWith('/entities?')) return [];
      if (path.startsWith(`/entities/${demoFiles}`))
        return {'#type': 'DataConnection', id: demoFiles, name: 'DemoFiles', namespace: 'System:', dataSource: 'Files'};
      if (method === 'POST') { stored[body.id] = {...body, relations: []}; return body; }
      const id = path.split('?')[0].split('/').pop()!;
      return stored[id] ?? notFound();
    });
    await push(dapi, bundleOf([project, table]), {onConflict: 'fail'}, () => {});
    const written = calls.filter((c) => c.path === '/projects?saveRelations=true');
    expect(written[written.length - 1].body.relations.map((r: any) => r.entity.id)).toEqual([TABLE_ID]);
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

  // A platform group is on every instance under its own id, so a grant naming it always lands.
  it('does not fail a project granted to a platform group the bundle still carries', async () => {
    const builtin: [string, BundleEntity] = ['a4b45840-9a50-11e6-9cc9-8546b8bf62e6', {type: 'UserGroup', json: {
      '#type': 'UserGroup', id: 'a4b45840-9a50-11e6-9cc9-8546b8bf62e6', name: 'AllUsers', friendlyName: 'All users',
    }}];
    const project: [string, BundleEntity] = [PROJECT_ID, {type: 'Project', json: {
      '#type': 'Project', id: PROJECT_ID, name: 'Dash', namespace: 'Admin:', isDashboard: true,
      _grants: [{group: 'All users', permission: 'View'}],
    }}];
    const {dapi} = storingDapi();
    const result = await push(dapi, bundleOf([builtin, project]), {onConflict: 'fail'}, () => {});
    expect(rowOf(result, 'UserGroup')).toMatchObject({action: 'skip', reason: 'platform_group'});
    expect(rowOf(result, 'Project').action).toBe('create');
  });

  // A reference the bundle does not carry is only harmless when the target holds that id itself.
  it('tells a dangling reference apart from one the target already has', async () => {
    const ON_TARGET = '11111111-2222-3333-4444-555555555555';
    const DANGLING = '66666666-7777-8888-9999-000000000000';
    const query = (id: string, name: string, connId: string): [string, BundleEntity] =>
      [id, {type: 'DataQuery', json: {'#type': 'DataQuery', id, name, namespace: 'Admin:',
        query: 'select 1', connection: {id: connId}}}];
    const bundle = bundleOf([query(QUERY_ID, 'Good', ON_TARGET), query(TWIN_ID, 'Bad', DANGLING)]);
    const {dapi} = makeDapi((_m, path) => {
      if (path.startsWith(`/entities/${ON_TARGET}`)) return [{'#type': 'DataConnection', id: ON_TARGET, name: 'Shared'}];
      if (path.startsWith('/entities?')) return [];
      return notFound();
    });
    const {rows} = await plan(dapi, bundle, {onConflict: 'fail'});
    expect(rows.find((r) => r.name === ON_TARGET)).toMatchObject({action: 'info', reason: 'orphan_ref'});
    const missing = rows.find((r) => r.name === DANGLING)!;
    expect(missing).toMatchObject({action: 'warn', reason: 'dependency_missing'});
    expect(missing.detail).toContain('Admin:Bad');
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
    '#type': 'FileInfo', id: 'f1', name: 'a.csv', friendlyName: 'a.csv', path: 'a.csv', isFile: true,
    connection: {id: CONN_ID}, ...over,
  }}];

  it('writes a datasync table\'s share file back into the share of the same name', async () => {
    const table: [string, BundleEntity] = [TABLE_ID, {type: 'TableInfo', json: {
      '#type': 'TableInfo', id: TABLE_ID, name: 'Compounds', namespace: 'Admin:',
      metaParams: {'.data-sync': 'sync', '.script': 'Compounds = OpenFile("User:Home/Test projects/compounds.csv")'},
    }}];
    const b = bundleOf([table]);
    writeShares(b.dir, new Map([['User:Home/Test projects/compounds.csv', Buffer.from('smiles,id\nCCO,1\n')]]));
    const stored: Record<string, any> = {};
    const {dapi, calls} = makeDapi((method, path, body) => {
      if (path.startsWith('/entities?')) return [];
      if (method === 'POST') { stored[body.id] = body; return body; }
      const id = path.split('?')[0].split('/').pop()!;
      return stored[id] ?? notFound();
    });
    const result = await push(dapi, b, {onConflict: 'fail'}, () => {});
    const upload = calls.find((c) => c.path.startsWith('/public/v1/files/'));
    expect(upload!.path).toBe('/public/v1/files/User.Home/Test projects/compounds.csv');
    expect(upload!.body.toString()).toBe('smiles,id\nCCO,1\n');
    expect(result.items.find((r) => r.entityType === 'File')).toMatchObject({action: 'create', reason: 'share_file'});
  });

  it('names the users the target is missing before anything is written', async () => {
    const mine: [string, BundleEntity] = ['p-a', {type: 'Project', json: {
      '#type': 'Project', id: 'p-a', name: 'Alice', isRoot: true, _personalOf: 'alice', relations: [],
    }}];
    const theirs: [string, BundleEntity] = ['p-b', {type: 'Project', json: {
      '#type': 'Project', id: 'p-b', name: 'Bob', isRoot: true, _personalOf: 'bob', relations: [],
    }}];
    const {dapi} = makeDapi((_m, path) => {
      if (path.startsWith('/entities?')) return [];
      return notFound();
    }, {'/users?': [{login: 'alice', project: {id: 'target-alice'}}]});
    const {rows} = await plan(dapi, bundleOf([mine, theirs]), {onConflict: 'skip'});
    const warn = rows.find((r) => r.reason === 'user_missing')!;
    expect(warn).toMatchObject({entityType: 'User', action: 'warn', name: 'bob'});
  });

  it('survives an external whose type the registry does not carry', async () => {
    // A bundle references types the walker never migrates (UserReport among them); once a mapping
    // for one is in the idmap, re-verifying it must not assume a registry entry exists.
    const SOURCE = '11111111-eeee-0000-0000-000000000001';
    const query: [string, BundleEntity] = [QUERY_ID, {type: 'DataQuery', json: {
      '#type': 'DataQuery', id: QUERY_ID, name: 'Orders', namespace: 'Skalkin:', query: 'select 1',
    }}];
    const b = bundleOf([query]);
    b.manifest.externals = [{id: SOURCE, type: 'UserReport', nqName: 'Admin:Report1'}];
    b.idmap[SOURCE] = '22222222-eeee-0000-0000-000000000002';
    const stored: Record<string, any> = {};
    const {dapi} = makeDapi((method, path, body) => {
      if (path.startsWith('/entities?')) return [];
      if (method === 'POST') { stored[body.id] = body; return body; }
      const id = path.split('?')[0].split('/').pop()!;
      return stored[id] ?? notFound();
    });
    const result = await push(dapi, b, {onConflict: 'fail'}, () => {});
    expect(result.items.some((r) => r.action === 'failed')).toBe(false);
  });

  it('re-points a reference at the package entity the target made under its own id', async () => {
    const SOURCE_CONN = '11111111-2222-3333-4444-555555555555';
    const TARGET_CONN = '66666666-7777-8888-9999-000000000000';
    const query: [string, BundleEntity] = [QUERY_ID, {type: 'DataQuery', json: {
      '#type': 'DataQuery', id: QUERY_ID, name: 'Orders', namespace: 'Skalkin:',
      connection: {id: SOURCE_CONN}, query: 'select 1',
    }}];
    const b = bundleOf([query]);
    b.manifest.externals = [{id: SOURCE_CONN, type: 'DataConnection', nqName: 'Samples:PostgresNorthwind'}];
    const stored: Record<string, any> = {};
    const {dapi, calls} = makeDapi((method, path, body) => {
      if (path.startsWith('/entities?namespace=Samples%3A'))
        return [{'#type': 'DataConnection', id: TARGET_CONN, name: 'PostgresNorthwind', namespace: 'Samples:'}];
      if (path.startsWith('/entities?')) return [];
      if (method === 'POST') { stored[body.id] = body; return body; }
      const id = path.split('?')[0].split('/').pop()!;
      return stored[id] ?? notFound();
    });
    await push(dapi, b, {onConflict: 'fail'}, () => {});
    const saved = calls.find((c) => c.method === 'POST' && c.path === '/connectors/queries')!;
    expect(saved.body.connection.id).toBe(TARGET_CONN);
  });

  it('reports a view the source itself can no longer resolve as left behind, not as a failure', async () => {
    const DEAD = 'deadbeef-0000-0000-0000-000000000001';
    const view: [string, BundleEntity] = [VIEW_ID, {type: 'ViewInfo', json: {
      '#type': 'ViewInfo', id: VIEW_ID, name: 'Demog', namespace: 'Skalkin:', table: {id: DEAD},
    }}];
    const b = bundleOf([view]);
    b.manifest.dangling = [DEAD];
    const {dapi} = makeDapi((method, path) => {
      if (path.startsWith('/entities?')) return [];
      if (method === 'POST') throw new Error('Operation caused an exception');
      return notFound();
    });
    const result = await push(dapi, b, {onConflict: 'fail'}, () => {});
    const row = result.items.find((r) => r.entityType === 'ViewInfo' && r.name.includes('Demog'))!;
    expect(row).toMatchObject({action: 'skip', reason: 'dead_on_source'});
    expect(result.items.some((r) => r.action === 'failed')).toBe(false);
  });

  it('keeps a timeout a failure even when the entity references a dead id', async () => {
    const DEAD = 'deadbeef-0000-0000-0000-000000000003';
    const view: [string, BundleEntity] = [VIEW_ID, {type: 'ViewInfo', json: {
      '#type': 'ViewInfo', id: VIEW_ID, name: 'Demog', namespace: 'Skalkin:', table: {id: DEAD},
    }}];
    const b = bundleOf([view]);
    b.manifest.dangling = [DEAD];
    const {dapi} = makeDapi((method, path) => {
      if (path.startsWith('/entities?')) return [];
      if (method === 'POST') throw new Error('POST /views: no answer in 60000ms');
      return notFound();
    });
    const result = await push(dapi, b, {onConflict: 'fail'}, () => {});
    expect(result.items.find((r) => r.entityType === 'ViewInfo')!.action).toBe('failed');
  });

  it('still fails a refusal caused by a reference the source has and the target lacks', async () => {
    const ABSENT = 'deadbeef-0000-0000-0000-000000000002';
    const view: [string, BundleEntity] = [VIEW_ID, {type: 'ViewInfo', json: {
      '#type': 'ViewInfo', id: VIEW_ID, name: 'Demog', namespace: 'Skalkin:', table: {id: ABSENT},
    }}];
    const {dapi} = makeDapi((method, path) => {
      if (path.startsWith('/entities?')) return [];
      if (method === 'POST') throw new Error('Operation caused an exception');
      return notFound();
    });
    const result = await push(dapi, bundleOf([view]), {onConflict: 'fail'}, () => {});
    const row = result.items.find((r) => r.entityType === 'ViewInfo' && r.name.includes('Demog'))!;
    expect(row.action).toBe('failed');
    expect(row.detail).toContain('install the package first');
  });

  it('re-asserts the placement of a project whose own row is unchanged', async () => {
    const SPACE = '11111111-cccc-0000-0000-000000000001';
    const rel = {id: 'r1', isLink: false, entity: {'#type': 'EntityRecord', id: TABLE_ID}};
    const space: [string, BundleEntity] = [SPACE, {type: 'Project', json: {
      '#type': 'Project', id: SPACE, name: 'Space', namespace: 'Skalkin:', relations: [rel],
    }}];
    const table: [string, BundleEntity] = [TABLE_ID, {type: 'TableInfo', json: {
      '#type': 'TableInfo', id: TABLE_ID, name: 'Demog', namespace: 'Skalkin:',
    }}];
    // Identical at planning time, relation and all — then the relation is gone by the time
    // placement runs, exactly as the real owner's write would have taken it away. Skipping the
    // project here is what made a second push of the same bundle shed what the first placed.
    let reads = 0;
    const stored: Record<string, any> = {
      [TABLE_ID]: {'#type': 'TableInfo', id: TABLE_ID, name: 'Demog', namespace: 'Skalkin:'},
    };
    const {dapi, calls} = makeDapi((method, path, body) => {
      if (path.startsWith('/entities?')) return [];
      if (path === '/projects?saveRelations=true') { stored[body.id] = body; return body; }
      if (method === 'POST') { stored[body.id] = {...body, relations: []}; return stored[body.id]; }
      const id = path.split('?')[0].split('/').pop()!;
      if (id === SPACE)
        return {'#type': 'Project', id: SPACE, name: 'Space', namespace: 'Skalkin:',
          relations: reads++ === 0 ? [rel] : []};
      return stored[id] ?? notFound();
    });
    const result = await push(dapi, bundleOf([space, table]), {onConflict: 'skip'}, () => {});
    expect(result.items.find((r) => r.entityType === 'Project' && r.name.includes('Space'))!.action).toBe('identical');
    const written = calls.filter((c) => c.path === '/projects?saveRelations=true');
    expect(written.map((c) => c.body.id)).toContain(SPACE);
    expect(written[0].body.relations.map((r: any) => r.entity.id)).toContain(TABLE_ID);
  });

  it('lets the owner claim an entity instead of writing a release the server ignores', async () => {
    const OUTER = '11111111-dddd-0000-0000-000000000001';
    const INNER = '11111111-dddd-0000-0000-000000000002';
    // Both claim the table; the deeper one owns it, so the outer one has to hold it as a link or
    // the owner's write deletes the outer row and nothing ever puts it back.
    const outerRels = [
      {id: 'r1', isLink: false, entity: {'#type': 'EntityRecord', id: INNER}},
      {id: 'r2', isLink: false, entity: {'#type': 'EntityRecord', id: TABLE_ID}},
    ];
    const outer: [string, BundleEntity] = [OUTER, {type: 'Project', json: {
      '#type': 'Project', id: OUTER, name: 'Outer', namespace: 'Skalkin:', relations: outerRels,
    }}];
    const inner: [string, BundleEntity] = [INNER, {type: 'Project', json: {
      '#type': 'Project', id: INNER, name: 'Inner', namespace: 'Skalkin:',
      relations: [{id: 'r3', isLink: false, entity: {'#type': 'EntityRecord', id: TABLE_ID}}],
    }}];
    const table: [string, BundleEntity] = [TABLE_ID, {type: 'TableInfo', json: {
      '#type': 'TableInfo', id: TABLE_ID, name: 'Demog', namespace: 'Skalkin:',
    }}];
    // The outer project is identical, so nothing rewrites its rows — and it holds the table as a
    // containment claim the bundle gives to the inner one.
    const stored: Record<string, any> = {
      [OUTER]: {'#type': 'Project', id: OUTER, name: 'Outer', namespace: 'Skalkin:',
        relations: JSON.parse(JSON.stringify(outerRels))},
    };
    const {dapi, calls} = makeDapi((method, path, body) => {
      if (path.startsWith('/entities?')) return [];
      if (path === '/projects?saveRelations=true') { stored[body.id] = body; return body; }
      if (method === 'POST') { stored[body.id] = {...body, relations: []}; return stored[body.id]; }
      const id = path.split('?')[0].split('/').pop()!;
      return stored[id] ?? notFound();
    });
    const result = await push(dapi, bundleOf([outer, inner, table]), {onConflict: 'skip'}, () => {});
    expect(result.items.find((r) => r.name.includes('Outer'))!.action).toBe('identical');
    const wrote = (id: string) => calls.filter((c) => c.path === '/projects?saveRelations=true' && c.body.id === id).pop();
    // The owner claims it, and the server drops the other holder's containment row itself. Writing
    // the outer project to say "link" is ignored and recomputed, so it is not worth tens of seconds.
    expect(wrote(INNER)!.body.relations.find((r: any) => r.entity.id === TABLE_ID).isLink).toBe(false);
    expect(wrote(OUTER)).toBeUndefined();
  });

  it('links a dashboard inside a personal space, not just the space itself', async () => {
    const ROOT = '11111111-aaaa-bbbb-cccc-000000000001';
    const DASH = '22222222-aaaa-bbbb-cccc-000000000002';
    const TARGET_ROOT = '99999999-aaaa-bbbb-cccc-000000000009';
    const root: [string, BundleEntity] = [ROOT, {type: 'Project', json: {
      '#type': 'Project', id: ROOT, name: 'Skalkin', isRoot: true, _personalOf: 'skalkin',
      relations: [{id: 'r0', isLink: false, entity: {'#type': 'EntityRecord', id: DASH}}],
    }}];
    const dash: [string, BundleEntity] = [DASH, {type: 'Project', json: {
      '#type': 'Project', id: DASH, name: 'Accelerometer', namespace: 'Skalkin:', isDashboard: true,
      relations: [{id: 'r1', isLink: false, entity: {'#type': 'EntityRecord', id: TABLE_ID}}],
    }}];
    const table: [string, BundleEntity] = [TABLE_ID, {type: 'TableInfo', json: {
      '#type': 'TableInfo', id: TABLE_ID, name: 'Accelerometer', namespace: 'Skalkin:',
    }}];
    const stored: Record<string, any> = {[TARGET_ROOT]: {'#type': 'Project', id: TARGET_ROOT, name: 'Skalkin', relations: []}};
    const {dapi, calls} = makeDapi((method, path, body) => {
      if (path.startsWith('/entities?')) return [];
      if (path === '/projects?saveRelations=true') { stored[body.id] = body; return body; }
      if (method === 'POST') { stored[body.id] = {...body, relations: []}; return stored[body.id]; }
      const id = path.split('?')[0].split('/').pop()!;
      return stored[id] ?? notFound();
    }, {'/users?': [{login: 'skalkin', project: {id: TARGET_ROOT}}]});
    await push(dapi, bundleOf([root, dash, table]), {onConflict: 'skip'}, () => {});
    const written = calls.filter((c) => c.path === '/projects?saveRelations=true').map((c) => c.body.id);
    expect(written).toContain(TARGET_ROOT);
    expect(written).toContain(DASH);
  });

  it('writes a container\'s relations before those of the projects inside it', async () => {
    const ROOT = '11111111-cccc-bbbb-aaaa-000000000001';
    const DASH = '22222222-cccc-bbbb-aaaa-000000000002';
    // Manifest order deliberately puts the dashboard first — writing the root after it would
    // rebuild the containment and drop what the dashboard had just been given.
    const dash: [string, BundleEntity] = [DASH, {type: 'Project', json: {
      '#type': 'Project', id: DASH, name: 'Dash', namespace: 'Admin:', isDashboard: true,
      relations: [{id: 'r1', isLink: false, entity: {'#type': 'EntityRecord', id: TABLE_ID}}],
    }}];
    const root: [string, BundleEntity] = [ROOT, {type: 'Project', json: {
      '#type': 'Project', id: ROOT, name: 'Space', isRoot: true,
      relations: [{id: 'r0', isLink: false, entity: {'#type': 'EntityRecord', id: DASH}}],
    }}];
    const table: [string, BundleEntity] = [TABLE_ID, {type: 'TableInfo', json: {
      '#type': 'TableInfo', id: TABLE_ID, name: 'Cereal', namespace: 'Admin:',
    }}];
    const stored: Record<string, any> = {};
    const {dapi, calls} = makeDapi((method, path, body) => {
      if (path.startsWith('/entities?')) return [];
      if (path === '/projects?saveRelations=true') { stored[body.id] = body; return body; }
      if (method === 'POST') { stored[body.id] = {...body, relations: []}; return stored[body.id]; }
      const id = path.split('?')[0].split('/').pop()!;
      return stored[id] ?? notFound();
    });
    await push(dapi, bundleOf([dash, root, table]), {onConflict: 'fail'}, () => {});
    const written = calls.filter((c) => c.path === '/projects?saveRelations=true').map((c) => c.body.id);
    expect(written).toEqual([ROOT, DASH]);
  });

  it('pushes a file that has no share of its own when its bytes travelled', async () => {
    // `files_service.dart` saves a connection-less FileInfo as a GUID-addressed blob
    // (`addToUserProject: f.connection == null`), so the blob is the entity.
    const orphan = file({connection: undefined});
    const bundle = bundleOf([orphan], new Map([['f1', Buffer.from([1, 2, 3])]]));
    const stored: Record<string, any> = {};
    const {dapi, calls} = makeDapi((method, path, body) => {
      if (path.startsWith('/entities?')) return [];
      if (method === 'POST') { stored[body.id] = body; return body; }
      const id = path.split('?')[0].split('/').pop()!;
      return stored[id] ?? notFound();
    });
    const result = await push(dapi, bundle, {onConflict: 'fail'}, () => {});
    expect(result.items.find((r) => r.entityType === 'FileInfo')!.action).toBe('create');
    expect(calls.some((c) => c.method === 'POST' && c.path.startsWith('/files'))).toBe(true);
  });

  it('skips a file that has neither a share nor its bytes', async () => {
    const orphan = file({connection: undefined});
    const {dapi, calls} = makeDapi((_m, path) => path.startsWith('/entities?') ? [] : notFound());
    const result = await push(dapi, bundleOf([orphan]), {onConflict: 'fail'}, () => {});
    expect(result.items.find((r) => r.entityType === 'FileInfo')).toMatchObject({
      action: 'skip', reason: 'file_bytes_missing'});
    expect(calls.some((c) => c.method === 'POST' && c.path.startsWith('/files'))).toBe(false);
  });

  it('writes the metadata first and the bytes under the id the target answered with', async () => {
    const bundle = bundleOf([file()], new Map([['f1', Buffer.from([1, 2, 3])]]));
    const saved: any = {'#type': 'FileInfo', id: 'existing-id', name: 'a.csv'};
    const {dapi, calls} = makeDapi((method, path) => {
      if (path.startsWith('/entities?')) return [];
      if (method === 'POST' && path === '/files') return saved;
      return path.startsWith('/files/existing-id') ? saved : notFound();
    });
    const result = await push(dapi, bundle, {onConflict: 'fail'}, () => {});
    expect(result.items.find((r) => r.entityType === 'FileInfo')).toMatchObject({action: 'create'});
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
