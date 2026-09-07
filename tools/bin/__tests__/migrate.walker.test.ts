import {describe, it, expect} from 'vitest';
import {NodeDapi} from '../utils/node-dapi';
import {Selection, compileFilter, expand, matchesName, normalizeSince, pullBytes, select} from '../utils/migrate/walker';

interface Call {method: string; path: string; body?: any}

function makeDapi(responder: (method: string, path: string, body?: any) => any) {
  const calls: Call[] = [];
  const client: any = {
    baseUrl: 'http://h/api',
    async request(method: string, path: string, body?: any) {
      calls.push({method, path, body});
      return responder(method, path, body);
    },
    get(path: string) { return this.request('GET', path); },
    post(path: string, body?: any) { return this.request('POST', path, body); },
    del(path: string) { return this.request('DELETE', path); },
  };
  return {dapi: new NodeDapi(client), calls};
}

function sel(over: Partial<Selection> = {}): Selection {
  return {types: [], names: [], ...over};
}

const PROJECT_ID = 'efe0b1f0-6fc0-11f1-b275-83ec2160b5e9';
const TABLE_ID = 'f2444470-6fc0-11f1-b307-bdcf313f748a';
const VIEW_ID = 'f24a10d0-6fc0-11f1-e396-fd8fa45dfcf4';
const QUERY_ID = 'aaaaaaaa-1111-2222-3333-444444444444';
const CONN_ID = 'bbbbbbbb-1111-2222-3333-444444444444';
const SYSTEM_CONN_ID = 'cccccccc-1111-2222-3333-444444444444';
const GROUP_ID = 'dddddddd-1111-2222-3333-444444444444';
const PERSONAL_ID = 'eeeeeeee-1111-2222-3333-444444444444';

describe('compileFilter', () => {
  it('compiles each flag', () => {
    expect(compileFilter(sel({namespace: 'Chem'}))).toBe('namespace starts "Chem:"');
    expect(compileFilter(sel({namespace: 'Chem:'}))).toBe('namespace starts "Chem:"');
    expect(compileFilter(sel({author: 'alice'}))).toBe('author = "alice"');
    expect(compileFilter(sel({since: '-2w'}))).toBe('updatedOn > -2w');
    expect(compileFilter(sel({filter: 'dataSource = "Postgres"'}))).toBe('(dataSource = "Postgres")');
    expect(compileFilter(sel())).toBe('');
  });

  it('escapes a quote in a user-supplied value', () => {
    expect(compileFilter(sel({author: 'a"b'}))).toBe('author = "a\\"b"');
  });

  it('sends --name as free text, but only when it is the whole filter', () => {
    expect(compileFilter(sel({name: 'Cereal'}))).toBe('Cereal');
    expect(compileFilter(sel({name: 'Cereal*'}))).toBe('Cereal');
    expect(compileFilter(sel({name: 'Cereal', namespace: 'Chem'}))).toBe('namespace starts "Chem:"');
  });

  it('ANDs the clauses in flag order', () => {
    expect(compileFilter(sel({namespace: 'Chem', author: 'alice', since: '-1d', filter: 'a = 1'})))
      .toBe('namespace starts "Chem:" and author = "alice" and updatedOn > -1d and (a = 1)');
  });
});

describe('matchesName', () => {
  it('matches either name, case-insensitively, and honours the glob', () => {
    expect(matchesName('CerealDemog', {name: 'CerealDemog', friendlyName: 'cereal, demog'})).toBe(true);
    expect(matchesName('cerealdemog', {name: 'CerealDemog'})).toBe(true);
    expect(matchesName('cereal, demog', {name: 'CerealDemog', friendlyName: 'cereal, demog'})).toBe(true);
    expect(matchesName('Cereal', {name: 'CerealDemog'})).toBe(false);
    expect(matchesName('Cereal*', {name: 'CerealDemog'})).toBe(true);
  });
});

describe('normalizeSince', () => {
  it('reads the sugar form and rejects the value minimist ate', () => {
    expect(normalizeSince(undefined)).toBeUndefined();
    expect(normalizeSince('2w')).toBe('-2w');
    expect(normalizeSince('-2w')).toBe('-2w');
    expect(normalizeSince('2026-01-01')).toBe('"2026-01-01"');
    expect(() => normalizeSince(true)).toThrow(/needs a value/);
    expect(() => normalizeSince('last tuesday')).toThrow(/timespan/);
  });
});

describe('select', () => {
  it('resolves a positional UUID through /entities and unwraps the array', async () => {
    const {dapi, calls} = makeDapi((_m, path) => {
      if (path === `/entities/${TABLE_ID}`) return [{'#type': 'TableInfo', id: TABLE_ID, name: 'Cereal'}];
      if (path.startsWith(`/tables/${TABLE_ID}`)) return {'#type': 'TableInfo', id: TABLE_ID, name: 'Cereal', rowCount: 77};
      throw new Error(`unexpected ${path}`);
    });
    const picked = await select(dapi, sel({names: [TABLE_ID]}), () => {});
    expect([...picked.keys()]).toEqual([TABLE_ID]);
    expect(picked.get(TABLE_ID)!.json.rowCount).toBe(77);
    expect(calls[0].path).toBe(`/entities/${TABLE_ID}`);
  });

  it('resolves a positional nqName through namespace + name', async () => {
    const {dapi, calls} = makeDapi((_m, path) => {
      if (path.startsWith('/entities?')) return [{'#type': 'Project', id: PROJECT_ID, name: 'CerealDemog', namespace: 'Admin:'}];
      return {'#type': 'Project', id: PROJECT_ID, name: 'CerealDemog', namespace: 'Admin:'};
    });
    await select(dapi, sel({names: ['Admin:CerealDemog']}), () => {});
    expect(calls[0].path).toBe('/entities?namespace=Admin%3A&name=CerealDemog');
  });

  it('sends the empty namespace of a root entity and ignores namespaced homonyms', async () => {
    const {dapi, calls} = makeDapi((_m, path) => {
      if (path.startsWith('/entities?')) return [
        {'#type': 'Project', id: PROJECT_ID, name: 'Shared', namespace: 'Admin:'},
        {'#type': 'Project', id: 'root-id', name: 'Shared'},
      ];
      return {'#type': 'Project', id: 'root-id', name: 'Shared'};
    });
    const picked = await select(dapi, sel({names: ['Shared']}), () => {});
    expect(calls[0].path).toBe('/entities?namespace=&name=Shared');
    expect([...picked.keys()]).toEqual(['root-id']);
  });

  it('lists each requested type with the compiled filter and tag param', async () => {
    const {dapi, calls} = makeDapi(() => []);
    await select(dapi, sel({types: ['Script', 'ViewLayout'], tag: 'demo', author: 'alice'}), () => {});
    expect(calls.map((c) => c.path)).toEqual(['/scripts?text=author%20%3D%20%22alice%22&tags=demo&limit=500&page=1']);
  });

  it('applies the --name glob client-side to what the free-text search returned', async () => {
    const {dapi} = makeDapi((_m, path) => {
      if (path.startsWith('/scripts?text=')) return [
        {'#type': 'Script', id: 's1', name: 'Descriptors'},
        {'#type': 'Script', id: 's2', name: 'DescriptorsDemo'},
      ];
      return {'#type': 'Script', id: 's1', name: 'Descriptors'};
    });
    const picked = await select(dapi, sel({types: ['Script'], name: 'Descriptors'}), () => {});
    expect([...picked.keys()]).toEqual(['s1']);
  });

  it('accepts a match that carries no namespace — the server-side filter already applied it', async () => {
    const job = {'#type': 'DataJob', id: 'j1', name: 'Nightly'};
    const {dapi} = makeDapi((_m, path) => path.startsWith('/entities?') ? [job] : {...job, transformations: 'x'});
    const picked = await select(dapi, sel({names: ['Admin:Nightly']}), () => {});
    expect([...picked.keys()]).toEqual(['j1']);
  });

  it('warns instead of listing a type whose router has no tags param', async () => {
    const notes: any[] = [];
    const {dapi} = makeDapi(() => []);
    await select(dapi, sel({types: ['ViewLayout'], tag: 'demo'}), (r) => notes.push(r));
    expect(notes).toEqual([{name: '', entityType: 'ViewLayout', action: 'warn', reason: 'tag_unsupported'}]);
  });

  it('lists spaces with their own clause, root spaces included', async () => {
    const {dapi, calls} = makeDapi(() => []);
    await select(dapi, sel({types: ['Project'], namespace: 'MigTestSpace',
      typeOptions: {Project: {clause: 'isDashboard = false and isEntity = false', params: {includeRoot: 'true'}}}}), () => {});
    expect(calls.map((c) => c.path)).toEqual([
      '/projects?text=namespace%20starts%20%22MigTestSpace%3A%22%20and%20isDashboard%20%3D%20false%20and%20isEntity%20%3D%20false&includeRoot=true&limit=500&page=1',
    ]);
  });

  it('lists files by type id — there is no type-wide /files listing', async () => {
    const {dapi, calls} = makeDapi(() => []);
    await select(dapi, sel({types: ['FileInfo'], namespace: 'MigTestSpace'}), () => {});
    expect(calls[0].path).toBe(
      '/entities?text=namespace%20starts%20%22MigTestSpace%3A%22&typeId=34d75630-e870-11e6-bfe1-590ff6f10d14&limit=500&page=1');
  });

  it('drops what never travels even when the walk is skipped (--no-deps)', async () => {
    const conns: any[] = [
      {'#type': 'DataConnection', id: 'c1', name: 'AppData', namespace: 'System:'},
      {'#type': 'DataConnection', id: 'c2', name: 'Home', namespace: 'Admin:'},
      {'#type': 'DataConnection', id: 'c3', name: 'Files', namespace: 'Spc:', parameters: {isProject: true}},
      {'#type': 'DataConnection', id: 'c4', name: 'NW', namespace: 'Admin:'},
    ];
    const notes: any[] = [];
    const {dapi} = makeDapi((_m, path) => path.startsWith('/connectors/connections?')
      ? conns : conns.find((c) => path.startsWith(`/connectors/connections/${c.id}`)));
    const picked = await select(dapi, sel({types: ['DataConnection']}), (r) => notes.push(r));
    expect([...picked.keys()]).toEqual(['c4']);
    expect(notes.map((n) => `${n.action}:${n.reason}`))
      .toEqual(['info:platform_connection', 'warn:personal_storage', 'info:space_files_connection']);
  });

  it('drops a package-owned entity and a file that lives in a share', async () => {
    const notes: any[] = [];
    const {dapi} = makeDapi((_m, path) => {
      if (path.startsWith('/packages/published/p1')) return {name: 'Chem'};
      if (path.startsWith('/scripts?')) return [{'#type': 'Script', id: 's1'}];
      if (path.startsWith('/scripts/s1')) return {'#type': 'Script', id: 's1', name: 'D', namespace: 'Chem:', package: {id: 'p1'}};
      if (path.startsWith('/entities?')) return [{'#type': 'FileInfo', id: 'f1'}];
      return {'#type': 'FileInfo', id: 'f1', name: 'a.csv', connection: {id: CONN_ID}};
    });
    const picked = await select(dapi, sel({types: ['Script', 'FileInfo']}), (r) => notes.push(r));
    expect(picked.size).toBe(0);
    expect(notes.map((n) => n.reason)).toEqual(['package_entity', 'file_in_share_not_migratable']);
  });

  it('does not list anything when only positionals are given', async () => {
    const {dapi, calls} = makeDapi((_m, path) => {
      if (path.startsWith('/entities?')) return [{'#type': 'Script', id: 's1', name: 'S'}];
      return {'#type': 'Script', id: 's1', name: 'S'};
    });
    await select(dapi, sel({types: ['Script'], names: ['S']}), () => {});
    expect(calls.some((c) => c.path.startsWith('/scripts?'))).toBe(false);
  });
});

function dashboardResponder(): (method: string, path: string) => any {
  return (_m, path) => {
    if (path.startsWith('/privileges/permissions?')) return [];
    if (path.startsWith('/projects/relations?')) return [
      {entity: {'#type': 'TableInfo', id: TABLE_ID}},
      {entity: {'#type': 'ViewInfo', id: VIEW_ID}},
    ];
    if (path.startsWith('/views?')) return [{'#type': 'ViewInfo', id: VIEW_ID}];
    if (path.startsWith('/layouts?')) return [];
    if (path.startsWith(`/tables/${TABLE_ID}`)) return {'#type': 'TableInfo', id: TABLE_ID, name: 'Cereal'};
    if (path.startsWith(`/views/${VIEW_ID}`)) return {'#type': 'ViewInfo', id: VIEW_ID, name: 'Cereal', table: {id: TABLE_ID}};
    if (path.startsWith(`/connectors/connections/${CONN_ID}`)) return {'#type': 'DataConnection', id: CONN_ID, name: 'NW', namespace: 'Admin:'};
    if (path.startsWith(`/connectors/connections/${SYSTEM_CONN_ID}`)) return {'#type': 'DataConnection', id: SYSTEM_CONN_ID, name: 'Datagrok', namespace: 'System:'};
    throw new Error(`unexpected ${path}`);
  };
}

describe('expand', () => {
  it('pulls a dashboard\'s tables and views, de-duped by id', async () => {
    const {dapi, calls} = makeDapi(dashboardResponder());
    const project = {'#type': 'Project', id: PROJECT_ID, name: 'CerealDemog', isDashboard: true};
    const out = await expand(dapi, new Map([[PROJECT_ID, {type: 'Project', json: project}]]), () => {});
    expect([...out.keys()].sort()).toEqual([PROJECT_ID, TABLE_ID, VIEW_ID].sort());
    expect(calls.filter((c) => c.path.startsWith(`/views/${VIEW_ID}`)).length).toBe(1);
  });

  it('pulls the connection a query references and drops the System one', async () => {
    const notes: any[] = [];
    const {dapi} = makeDapi(dashboardResponder());
    const query = {'#type': 'DataQuery', id: QUERY_ID, name: 'Q', connection: {id: CONN_ID}};
    const systemQuery = {'#type': 'DataQuery', id: 'q2', name: 'Q2', connection: {id: SYSTEM_CONN_ID}};
    const out = await expand(dapi, new Map<string, any>([
      [QUERY_ID, {type: 'DataQuery', json: query}],
      ['q2', {type: 'DataQuery', json: systemQuery}],
    ]), (r) => notes.push(r));
    expect(out.has(CONN_ID)).toBe(true);
    expect(out.has(SYSTEM_CONN_ID)).toBe(false);
    expect(notes).toEqual([{name: 'System:Datagrok', entityType: 'DataConnection', action: 'info', reason: 'platform_connection', detail: undefined}]);
  });

  it('drops a package-owned script, naming the package behind the published id', async () => {
    const notes: any[] = [];
    const {dapi, calls} = makeDapi((_m, path) => {
      if (path.startsWith('/packages/published/p1')) return {name: 'Chem', package: {id: 'pkg1'}};
      throw new Error(`unexpected ${path}`);
    });
    const script = {'#type': 'Script', id: 's1', name: 'Descriptors', namespace: 'Chem:', package: {id: 'p1'}};
    const out = await expand(dapi, new Map([['s1', {type: 'Script', json: script}]]), (r) => notes.push(r));
    expect(out.size).toBe(0);
    expect(calls[0].path).toBe('/packages/published/p1');
    expect(notes[0]).toMatchObject({entityType: 'Script', action: 'warn', reason: 'package_entity', detail: 'Chem'});
  });

  it('falls back to the project relations of the project itself when the route errors', async () => {
    const notes: any[] = [];
    const project = {'#type': 'Project', id: PROJECT_ID, name: 'CerealDemog', isDashboard: true,
      relations: [{id: 'r1', isLink: false, entity: {'#type': 'TableInfo', id: TABLE_ID}}]};
    const {dapi} = makeDapi((_m, path) => {
      if (path.startsWith('/projects/relations'))
        return {'#type': 'ApiError', message: 'NoSuchMethodError: createNew on null'};
      if (path.startsWith(`/tables/${TABLE_ID}`)) return {'#type': 'TableInfo', id: TABLE_ID, name: 'Cereal'};
      return [];
    });
    const out = await expand(dapi, new Map([[PROJECT_ID, {type: 'Project', json: project}]]), (r) => notes.push(r));
    expect([...out.keys()].sort()).toEqual([PROJECT_ID, TABLE_ID].sort());
    expect(notes[0]).toMatchObject({entityType: 'Project', action: 'warn', reason: 'relations_degraded'});
  });

  it('drops a file that lives inside a share — only a blob of its own migrates', async () => {
    const notes: any[] = [];
    const file = {'#type': 'FileInfo', id: 'f1', name: 'a.csv', connection: {id: CONN_ID}};
    const {dapi} = makeDapi(() => { throw new Error('nothing else is fetched'); });
    const out = await expand(dapi, new Map([['f1', {type: 'FileInfo', json: file}]]), (r) => notes.push(r));
    expect(out.size).toBe(0);
    expect(notes[0]).toMatchObject({entityType: 'FileInfo', action: 'info', reason: 'file_in_share_not_migratable'});
  });

  it('is the only step that walks relations, so --no-deps issues none', async () => {
    const {dapi, calls} = makeDapi((_m, path) => {
      if (path === `/entities/${PROJECT_ID}`) return [{'#type': 'Project', id: PROJECT_ID, name: 'CerealDemog'}];
      if (path.startsWith(`/projects/${PROJECT_ID}`)) return {'#type': 'Project', id: PROJECT_ID, name: 'CerealDemog'};
      throw new Error(`unexpected ${path}`);
    });
    const picked = await select(dapi, sel({names: [PROJECT_ID]}), () => {});
    expect(picked.size).toBe(1);
    expect(calls.some((c) => c.path.startsWith('/projects/relations'))).toBe(false);
  });
});

describe('grants and groups', () => {
  const granted = (over: any = {}) => (_m: string, path: string) => {
    if (path.startsWith('/privileges/permissions?')) return [
      {entityId: PROJECT_ID, userGroup: {id: GROUP_ID}, permission: {name: 'View'}},
      {entityId: PROJECT_ID, userGroup: {id: PERSONAL_ID}, permission: {name: 'Edit'}},
      {entityId: 'other', userGroup: {id: GROUP_ID}, permission: {name: 'Edit'}},
    ];
    if (path === `/groups/${GROUP_ID}`) return {
      '#type': 'UserGroup', id: GROUP_ID, name: 'Chemists', friendlyName: 'Chemists',
      children: [{child: {id: PERSONAL_ID}, isAdmin: true}, {child: {id: 'analysts'}}],
      parents: [], ...over,
    };
    if (path === `/groups/${PERSONAL_ID}`) return {'#type': 'UserGroup', id: PERSONAL_ID, name: 'alice', personal: true};
    if (path === `/groups/${PERSONAL_ID}/user`) return {login: 'alice', email: 'alice@x.io'};
    if (path === '/groups/analysts') return {'#type': 'UserGroup', id: 'analysts', name: 'Analysts', friendlyName: 'Analysts'};
    if (path === '/groups/a4b45840-9a50-11e6-9cc9-8546b8bf62e6') return {'#type': 'UserGroup', id: 'a4b45840-9a50-11e6-9cc9-8546b8bf62e6', name: 'AllUsers', friendlyName: 'All users'};
    if (path.startsWith('/projects/relations?') || path.startsWith('/views?') || path.startsWith('/layouts?')) return [];
    throw new Error(`unexpected ${path}`);
  };

  it('records the grant, pulls the holding group bare, and resolves its members', async () => {
    const {dapi} = makeDapi(granted());
    const project = {'#type': 'Project', id: PROJECT_ID, name: 'CerealDemog', isDashboard: true};
    const out = await expand(dapi, new Map([[PROJECT_ID, {type: 'Project', json: project}]]), () => {});

    expect(out.get(PROJECT_ID)!.json._grants).toEqual([{group: 'Chemists', permission: 'View'}]);
    expect(out.get(GROUP_ID)!.type).toBe('UserGroup');
    expect(out.get(GROUP_ID)!.json._members).toEqual([
      {kind: 'user', login: 'alice', email: 'alice@x.io', isAdmin: true},
      {kind: 'group', name: 'Analysts', isAdmin: false},
    ]);
    expect(out.has(PERSONAL_ID)).toBe(false);
  });

  it('bundles a member group too, so the membership can be replayed on a fresh server', async () => {
    const {dapi} = makeDapi(granted());
    const project = {'#type': 'Project', id: PROJECT_ID, name: 'CerealDemog', isDashboard: true};
    const out = await expand(dapi, new Map([[PROJECT_ID, {type: 'Project', json: project}]]), () => {});
    expect(out.get('analysts')).toMatchObject({type: 'UserGroup', json: {friendlyName: 'Analysts'}});
  });

  it('pulls a visible parent group but never All users', async () => {
    const parents = [
      {parent: {id: 'a4b45840-9a50-11e6-9cc9-8546b8bf62e6'}},
      {parent: {id: 'analysts'}},
    ];
    const {dapi} = makeDapi(granted({parents, children: []}));
    const project = {'#type': 'Project', id: PROJECT_ID, name: 'CerealDemog', isDashboard: true};
    const out = await expand(dapi, new Map([[PROJECT_ID, {type: 'Project', json: project}]]), () => {});
    expect([...out.keys()].sort()).toEqual([GROUP_ID, PROJECT_ID, 'analysts'].sort());
    expect(out.has('a4b45840-9a50-11e6-9cc9-8546b8bf62e6')).toBe(false);
  });
});

describe('pullBytes', () => {
  const entities = new Map<string, any>([
    [TABLE_ID, {type: 'TableInfo', json: {id: TABLE_ID, name: 'Cereal'}}],
    ['f1', {type: 'FileInfo', json: {id: 'f1', name: 'blob.txt'}}],
    ['s1', {type: 'Script', json: {id: 's1', name: 'S'}}],
  ]);

  function bytesDapi() {
    const paths: string[] = [];
    const {dapi} = makeDapi(() => null);
    (dapi.client as any).getBytes = async (path: string) => { paths.push(path); return Buffer.from('x'); };
    return {dapi, paths};
  }

  it('fetches only the kinds the flags asked for', async () => {
    const {dapi, paths} = bytesDapi();
    expect([...(await pullBytes(dapi, entities, () => {}, ['tables'])).keys()]).toEqual([TABLE_ID]);
    expect(paths).toEqual([`/tables/${TABLE_ID}/data`]);
  });

  it('fetches file bytes only with --include-files', async () => {
    const {dapi, paths} = bytesDapi();
    expect([...(await pullBytes(dapi, entities, () => {}, ['tables', 'files'])).keys()]).toEqual([TABLE_ID, 'f1']);
    expect(paths).toContain('/files/data/f1');
    const none = bytesDapi();
    expect((await pullBytes(none.dapi, entities, () => {}, [])).size).toBe(0);
    expect(none.paths).toEqual([]);
  });
});
