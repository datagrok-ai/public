/**
 * Integration tests for `grok s pull | push | bundle ls`.
 * Requires a running Datagrok server configured in ~/.grok/config.yaml.
 *
 * Environment variables:
 *   HOST   Server alias or URL (default: config default — never point this at a shared server)
 *
 * Run: npm run test:integration
 *
 * Everything it creates is named MigTest* and tagged `migrate-test`, and is deleted in afterAll.
 */

import {afterAll, beforeAll, describe, expect, it} from 'vitest';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';
import {randomUUID} from 'crypto';
import {NodeApiClient, NodeDapi} from '../utils/node-dapi';
import {getDevKey} from '../utils/test-utils';
import * as bundle from '../utils/migrate/bundle';
import {BytesKind, TYPES, typeOptionsOf} from '../utils/migrate/registry';
import {Selection, expand, pullBytes, select} from '../utils/migrate/walker';
import {ConflictPolicy, push} from '../utils/migrate/pusher';
import {reId} from './migrate-helpers';

const HOST = process.env['HOST'] ?? '';
const TMP = fs.mkdtempSync(path.join(os.tmpdir(), 'migtest-'));
const BUNDLE = path.join(TMP, 'bundle');
const COPY = path.join(TMP, 'bundle-copy');
const ADOPT = path.join(TMP, 'bundle-adopt');
const SPACE = path.join(TMP, 'space');
const SPACE_COPY = path.join(TMP, 'space-copy');
const FILE_BYTES = Buffer.from('MigTest file bytes\n');

let dapi: NodeDapi;
let offline = false;

const ids = {conn: randomUUID(), query: randomUUID(), dash: randomUUID(), group: randomUUID(),
  space: randomUUID(), spaceConn: randomUUID(), file: randomUUID(), notebook: randomUUID()};
const copyId = randomUUID();
const adoptTwin = randomUUID();
const copies: Record<string, string> = {};
const adopted: Record<string, string> = {};
const spaceCopies: Record<string, string> = {};

beforeAll(async () => {
  try {
    const {url, key} = getDevKey(HOST);
    dapi = new NodeDapi(await NodeApiClient.login(url, key));
  } catch (err) {
    console.warn(`[integration] Server unreachable — all tests will be skipped: ${err}`);
    offline = true;
    return;
  }
  await dapi.internal('/connectors/connections').save({
    '#type': 'DataConnection', id: ids.conn, name: 'MigTestConn', friendlyName: 'MigTest Conn',
    dataSource: 'PostgresDart', parameters: {server: 'localhost:5432', db: 'northwind', ssl: false, password: 'migtest'},
  });
  await dapi.internal('/connectors/queries').save({
    '#type': 'DataQuery', id: ids.query, name: 'MigTestQuery', friendlyName: 'MigTest Query',
    query: 'select 1 as x', connection: {id: ids.conn},
  });
  await dapi.internal('/projects').save({
    '#type': 'Project', id: ids.dash, name: 'MigTestDash', friendlyName: 'MigTest Dash', isDashboard: true,
  });
  await dapi.internal('/groups/').save({
    '#type': 'UserGroup', id: ids.group, name: 'MigTestGroup', friendlyName: 'MigTestGroup',
  });
  await dapi.internal('/connectors/connections').save({
    '#type': 'DataConnection', id: ids.spaceConn, name: 'MigTestSpaceConn', friendlyName: 'MigTest Space Conn',
    dataSource: 'PostgresDart', parameters: {server: 'localhost:5432', db: 'northwind', ssl: false},
  });
  // A file entity is a blob addressed by its own id; a file inside a share has no entity to migrate.
  await dapi.internal('/files').save({
    '#type': 'FileInfo', id: ids.file, isFile: true, size: FILE_BYTES.length,
    name: 'MigTestFile.txt', friendlyName: 'MigTestFile.txt', path: 'MigTestFile.txt',
  });
  await dapi.client.putBytes('/files/data/' + ids.file, FILE_BYTES);
  await dapi.internal('/projects').save({
    '#type': 'Project', id: ids.space, name: 'MigTestSpace', friendlyName: 'MigTestSpace',
    isDashboard: false, isRoot: true,
  });
  const space = await dapi.internal('/projects').find(ids.space);
  delete space.storage;
  space.relations = [
    ...(space.relations ?? []).map((r: any) => ({id: r.id, entity: {'#type': 'EntityRecord', id: r.entity.id}, isLink: false})),
    {id: randomUUID(), entity: {'#type': 'EntityRecord', id: ids.spaceConn}, isLink: false},
    {id: randomUUID(), entity: {'#type': 'EntityRecord', id: ids.file}, isLink: false},
  ];
  await dapi.internal('/projects').save(space, 'saveRelations=true');
  await dapi.internal('/notebooks').save({
    '#type': 'Notebook', id: ids.notebook, name: 'MigTestNotebook', friendlyName: 'MigTest Notebook',
    description: 'migrate test',
    notebook: {cells: [], metadata: {kernelspec: {name: 'python3', display_name: 'Python 3'}}, nbformat: 4, nbformat_minor: 2},
  });
  await dapi.groups.addMembers(ids.group, ['admin'], false, true);
  await dapi.shares.share(ids.dash, 'MigTestGroup', 'View');
  await dapi.client.post('/entities/tag?tag=migrate-test', [ids.conn, ids.query, ids.dash]);
});

afterAll(async () => {
  try {
    if (offline) return;
    // A group holding permissions cannot be deleted, and a permission whose entity is
    // already gone can no longer be revoked — so the grants go first, entities second.
    for (const group of [ids.group, copies[ids.group], adopted[ids.group]].filter(Boolean))
      for (const perm of await dapi.client.get(`/privileges/permissions?groupId=${group}&all=true`).catch(() => []))
        if (perm.userGroup?.id === group)
          await dapi.internal('/privileges/permissions').delete(perm.id).catch(() => {});
    const routes: [string, string][] = [
      ['/notebooks', ids.notebook], ['/notebooks', spaceCopies[ids.notebook]],
      ['/files', ids.file], ['/files', spaceCopies[ids.file]],
      ['/projects', ids.space], ['/projects', spaceCopies[ids.space]],
      ['/connectors/connections', ids.spaceConn], ['/connectors/connections', spaceCopies[ids.spaceConn]],
      ['/projects', ids.dash], ['/projects', copyId], ['/projects', copies[ids.dash]], ['/projects', adopted[ids.dash]],
      ['/connectors/queries', ids.query], ['/connectors/queries', copies[ids.query]], ['/connectors/queries', adopted[ids.query]],
      ['/connectors/connections', ids.conn], ['/connectors/connections', copies[ids.conn]],
      ['/connectors/connections', adopted[ids.conn]], ['/connectors/connections', adoptTwin],
      ['/groups', ids.group], ['/groups', copies[ids.group]], ['/groups', adopted[ids.group]],
    ];
    for (const [route, id] of routes)
      if (id) await dapi.internal(route).delete(id).catch(() => {});
  } finally {
    fs.rmSync(TMP, {recursive: true, force: true});
  }
});

function stest(name: string, fn: () => Promise<void>) {
  return it(name, async (ctx) => {
    if (offline) { ctx.skip(); return; }
    await fn();
  });
}

function selection(names: string[]): Selection {
  return {types: [], names};
}

async function pullInto(dir: string, sel: Selection, opts: {replace?: boolean} = {}, kinds: BytesKind[] = ['tables']) {
  const notes: any[] = [];
  const note = (row: any) => notes.push(row);
  const entities = await expand(dapi, await select(dapi, sel, note), note);
  const info = await dapi.serverInfo();
  const user = await dapi.client.get('/users/current');
  const manifest = bundle.write(dir, entities, {
    source: {url: dapi.client.baseUrl, version: info.version, commit: info.commit, userNamespace: `${user.project.name}:`},
    args: ['pull', ...sel.names],
    packages: [],
  }, opts, await pullBytes(dapi, entities, note, kinds));
  return {entities, manifest, notes};
}

const pull = (names: string[], opts: {replace?: boolean} = {}) => pullInto(BUNDLE, selection(names), opts);

/** Fresh ids + a MigTest* → MigTest<Suffix>* rename, so the copy can be pushed to the same server. */
function writeCopy(from: string, dir: string, suffix: string, into: Record<string, string>): void {
  const source = bundle.read(from);
  const {entities, idmap} = reId(source.entities, (name) => name.replace(/^MigTest/, `MigTest${suffix}`));
  Object.assign(into, idmap);
  const bytes = new Map<string, Buffer>();
  for (const [id, e] of source.entities) {
    const kind = TYPES[e.type].bytes?.kind;
    const file = kind ? bundle.bytesPath(from, kind, id) : '';
    if (file && fs.existsSync(file))
      bytes.set(idmap[id], fs.readFileSync(file));
  }
  bundle.write(dir, entities, {source: source.manifest.source, args: ['copy'], packages: []}, {replace: true}, bytes);
}

async function countingPush(dir: string, onConflict: ConflictPolicy = 'fail') {
  const client: any = dapi.client;
  const original = client.request.bind(client);
  let posts = 0;
  client.request = (method: string, p: string, body?: any) => {
    if (method === 'POST') posts++;
    return original(method, p, body);
  };
  try {
    return {result: await push(dapi, bundle.read(dir), {onConflict}, () => {}), posts};
  } finally {
    client.request = original;
  }
}

describe('grok s pull / push', () => {
  stest('pulls a dashboard and a query with its connection', async () => {
    const {entities, manifest} = await pull(['Admin:MigTestDash', ids.query]);
    expect([...entities.keys()].sort()).toEqual([ids.conn, ids.query, ids.dash, ids.group].sort());

    for (const file of ['Project/Admin.MigTestDash.json', 'DataQuery/Admin.MigTestQuery.json',
      'DataConnection/Admin.MigTestConn.json', 'UserGroup/MigTestGroup.json'])
      expect(fs.existsSync(path.join(BUNDLE, file)), file).toBe(true);

    expect(manifest.source.version).toMatch(/^\d+\.\d+/);
    expect(manifest.order.map((e) => e.type)).toEqual(['UserGroup', 'DataConnection', 'Project', 'DataQuery']);
  });

  stest('never writes a credential into the bundle', async () => {
    const conn = JSON.parse(fs.readFileSync(path.join(BUNDLE, 'DataConnection/Admin.MigTestConn.json'), 'utf8'));
    expect(conn.parameters.password).toBeUndefined();
    expect(conn._credentials).toEqual(['password']);
    expect(JSON.stringify(conn)).not.toContain('migtest');
    expect(JSON.stringify(conn)).not.toContain('_____________');
  });

  stest('reduces the nested query connection to an id reference', async () => {
    const query = JSON.parse(fs.readFileSync(path.join(BUNDLE, 'DataQuery/Admin.MigTestQuery.json'), 'utf8'));
    expect(query.connection).toEqual({id: ids.conn});
  });

  stest('records the tags of every tagged entity, and no tag rows', async () => {
    const dash = JSON.parse(fs.readFileSync(path.join(BUNDLE, 'Project/Admin.MigTestDash.json'), 'utf8'));
    expect(dash._tags).toEqual(['migrate-test']);
    expect(dash.entityTags).toBeUndefined();
  });

  stest('records the grant on the dashboard and the group members by login', async () => {
    const dash = JSON.parse(fs.readFileSync(path.join(BUNDLE, 'Project/Admin.MigTestDash.json'), 'utf8'));
    expect(dash._grants).toContainEqual({group: 'MigTestGroup', permission: 'View'});

    const group = JSON.parse(fs.readFileSync(path.join(BUNDLE, 'UserGroup/MigTestGroup.json'), 'utf8'));
    expect(group.parents).toEqual([]);
    expect(group.children).toEqual([]);
    expect(group._members.map((m: any) => m.login)).toContain('admin');
  });

  stest('plans an unchanged bundle as entirely identical, writing nothing', async () => {
    const result = await push(dapi, bundle.read(BUNDLE), {dryRun: true, onConflict: 'fail'}, () => {});
    expect(result.status).toBe('dry-run');
    expect(result.items.map((r) => `${r.entityType}:${r.action}`).sort())
      .toEqual(['DataConnection:identical', 'DataQuery:identical', 'Project:identical', 'UserGroup:identical']);
  });

  stest('merges a second pull into the same directory without losing anything', async () => {
    const before = fs.readdirSync(path.join(BUNDLE, 'DataConnection'));
    const {manifest} = await pull([ids.query]);
    expect(manifest.pulls.length).toBe(2);
    expect(manifest.order.length).toBe(4);
    expect(fs.readdirSync(path.join(BUNDLE, 'DataConnection'))).toEqual(before);
    expect(fs.existsSync(path.join(BUNDLE, 'Project/Admin.MigTestDash.json'))).toBe(true);
    expect(manifest.pulls[1].ids.sort()).toEqual([ids.conn, ids.query].sort());
  });

  stest('refuses to merge a bundle pulled from another server', async () => {
    expect(() => bundle.write(BUNDLE, new Map(), {
      source: {url: 'http://elsewhere/api', version: '1.27.9', userNamespace: 'Admin:'}, args: [], packages: [],
    }, {})).toThrow(/use --replace/);
  });

  stest('creates a dashboard that is not on the target yet', async () => {
    const copyDir = path.join(TMP, 'bundle-single');
    const source = bundle.read(BUNDLE);
    const project = JSON.parse(JSON.stringify(source.entities.get(ids.dash)!.json));
    project.id = copyId;
    project.name = 'MigTestDashCopy';
    project.friendlyName = 'MigTest Dash Copy';
    delete project._grants;
    bundle.write(copyDir, new Map([[copyId, {type: 'Project', json: project}]]),
      {source: source.manifest.source, args: ['copy'], packages: []}, {replace: true});

    const result = await push(dapi, bundle.read(copyDir), {onConflict: 'fail'}, () => {});
    expect(result.items.map((r) => r.action)).toEqual(['create', 'warn']);
    expect(result.items[1].reason).toBe('not_visible');
    expect(result.status).toBe('ok');
    const created = await dapi.internal('/projects').find(copyId);
    expect(created.name).toBe('MigTestDashCopy');
    expect(created.metaParams.sync_id).toBe(copyId);
  });

  stest('pushes a re-id\'d copy with its group, members and grant', async () => {
    writeCopy(BUNDLE, COPY, 'Copy', copies);
    const {result} = await countingPush(COPY);
    expect(result.items.filter((r) => r.action === 'failed')).toEqual([]);
    expect(result.counts['create']).toBe(4);

    const group = await dapi.internal('/groups').find(copies[ids.group]);
    expect(group.name).toBe('MigTestCopyGroup');
    expect(result.items.find((r) => r.entityType === 'UserGroup')!.detail).toContain('members: 1 matched');
    const [personal] = (await dapi.groups.lookup('admin')).filter((g: any) => g.personal);
    expect((group.children ?? []).map((c: any) => c.child?.id)).toContain(personal.id);

    const perms = await dapi.internal('/privileges/permissions').list({entityId: copies[ids.dash], all: 'true'});
    expect(perms.some((p: any) => p.entityId === copies[ids.dash] && p.userGroup?.id === copies[ids.group])).toBe(true);

    const query = await dapi.internal('/connectors/queries').find(copies[ids.query]);
    expect(query.connection.id).toBe(copies[ids.conn]);
  });

  stest('re-pushes the same copy without writing anything', async () => {
    const {result, posts} = await countingPush(COPY);
    expect(result.items.map((r) => `${r.entityType}:${r.action}`).sort())
      .toEqual(['DataConnection:identical', 'DataQuery:identical', 'Project:identical', 'UserGroup:identical']);
    expect(posts).toBe(0);
  });

  stest('replays the tags onto the copy', async () => {
    const copy = await dapi.internal('/projects').find(copies[ids.dash]);
    expect((copy.entityTags ?? []).map((t: any) => t.tag)).toEqual(['migrate-test']);
  });

  stest('adopts a same-name twin and points the bundle references at it', async () => {
    writeCopy(BUNDLE, ADOPT, 'Adopt', adopted);
    await dapi.internal('/connectors/connections').save({
      '#type': 'DataConnection', id: adoptTwin, name: 'MigTestAdoptConn', friendlyName: 'MigTestAdopt Conn',
      dataSource: 'PostgresDart', parameters: {server: 'localhost:5432', db: 'northwind', ssl: false},
    });

    const {result} = await countingPush(ADOPT, 'adopt');
    expect(result.items.filter((r) => r.action === 'failed')).toEqual([]);
    expect(result.items.find((r) => r.entityType === 'DataConnection')!.reason).toBe(`adopted ${adoptTwin}`);
    expect(JSON.parse(fs.readFileSync(path.join(ADOPT, 'idmap.json'), 'utf8'))[adopted[ids.conn]]).toBe(adoptTwin);

    expect(await dapi.internal('/connectors/connections').find(adopted[ids.conn])).toBeNull();
    const query = await dapi.internal('/connectors/queries').find(adopted[ids.query]);
    expect(query.connection.id).toBe(adoptTwin);
  });

  stest('re-pushes the adopted bundle as entirely identical', async () => {
    const {result, posts} = await countingPush(ADOPT, 'adopt');
    expect(result.counts).toEqual({identical: 4});
    expect(posts).toBe(0);
  });
});

describe('spaces, files and notebooks', () => {
  const spaceSelection = (): Selection => ({
    types: ['Project', 'DataConnection', 'FileInfo'], names: [], space: 'MigTestSpace',
    typeOptions: {Project: typeOptionsOf('space')!},
  });

  stest('pulls a space with the connection and the file under it', async () => {
    const {entities, notes} = await pullInto(SPACE, spaceSelection(), {replace: true}, ['tables', 'files']);
    expect([...entities.keys()].sort()).toEqual([ids.space, ids.spaceConn, ids.file].sort());
    // the space's own Files connection is re-created by the target, so it never travels
    expect(notes.some((n) => n.reason === 'space_files_connection')).toBe(true);

    const file = JSON.parse(fs.readFileSync(path.join(SPACE, 'FileInfo/MigTestFile.txt.json'), 'utf8'));
    expect(file.isFile).toBe(true);
    expect(fs.readFileSync(bundle.bytesPath(SPACE, 'files', ids.file))).toEqual(FILE_BYTES);
  });

  stest('pushes the renamed space with its members, and the file bytes travel', async () => {
    writeCopy(SPACE, SPACE_COPY, 'Spc', spaceCopies);
    const {result} = await countingPush(SPACE_COPY);
    expect(result.items.filter((r) => r.action === 'failed')).toEqual([]);
    expect(result.counts['create']).toBe(3);

    const under = await dapi.internal('/entities').list({namespace: 'MigTestSpcSpace:'});
    expect(under.map((e: any) => `${e['#type']}:${e.name}`).sort())
      .toEqual(['DataConnection:Files', 'DataConnection:MigTestSpcSpaceConn', 'FileInfo:MigTestSpcFile.txt'].sort());
    expect(await dapi.client.getBytes(`/files/data/${spaceCopies[ids.file]}`)).toEqual(FILE_BYTES);
  });

  stest('re-pushes the space copy without writing anything', async () => {
    const {result, posts} = await countingPush(SPACE_COPY);
    expect(result.items.filter((r) => ['create', 'update'].includes(r.action))).toEqual([]);
    expect(posts).toBe(0);
  });

  stest('round-trips a notebook', async () => {
    const dir = path.join(TMP, 'notebook');
    await pullInto(dir, {types: [], names: ['Admin:MigTestNotebook']}, {replace: true});
    const pulled = JSON.parse(fs.readFileSync(path.join(dir, 'Notebook/Admin.MigTestNotebook.json'), 'utf8'));
    expect(pulled.notebook.nbformat).toBe(4);

    writeCopy(dir, path.join(TMP, 'notebook-copy'), 'Spc', spaceCopies);
    const {result} = await countingPush(path.join(TMP, 'notebook-copy'));
    expect(result.items.filter((r) => r.action === 'failed')).toEqual([]);
    const copy = await dapi.internal('/notebooks').find(spaceCopies[ids.notebook]);
    expect(copy.name).toBe('MigTestSpcNotebook');
    expect(copy.notebook).toEqual(pulled.notebook);
  });
});
