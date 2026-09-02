import {describe, it, expect} from 'vitest';
import {
  PASSWORD_PLACEHOLDER, SAME_PASSWORD, TYPES, datasyncConnectionRefs, fileNameFor, nqNameOf, rankOf, resolveType,
  resolveTypes, stripCredentials, typeOptionsOf,
} from '../utils/migrate/registry';

describe('stripCredentials', () => {
  it('removes only the masked values and records what was lost', () => {
    const conn = {parameters: {server: 'localhost', password: PASSWORD_PLACEHOLDER, token: SAME_PASSWORD, ssl: false}, credentials: {x: 1}};
    stripCredentials(conn);
    expect(conn.parameters).toEqual({server: 'localhost', ssl: false});
    expect((conn as any)._credentials).toEqual(['password', 'token']);
    expect((conn as any).credentials).toBeUndefined();
  });

  it('leaves a connection with no masked parameters untouched', () => {
    const conn: any = {parameters: {db: 'northwind'}};
    stripCredentials(conn);
    expect(conn.parameters).toEqual({db: 'northwind'});
    expect(conn._credentials).toBeUndefined();
  });
});

describe('type strips', () => {
  it('reduces a nested DataQuery connection to its id', () => {
    const q: any = {connection: {id: 'c1', name: 'NW', parameters: {password: 'secret'}}};
    TYPES['DataQuery'].strip!(q);
    expect(q.connection).toEqual({id: 'c1'});
  });

  it('reduces a ViewInfo table to its id', () => {
    const v: any = {table: {id: 't1', metaParams: {}}};
    TYPES['ViewInfo'].strip!(v);
    expect(v.table).toEqual({id: 't1'});
  });

  it('drops project storage and compacts relations deterministically', () => {
    const p: any = {
      storage: {id: 's1'},
      relations: [
        {id: 'r2', isLink: false, entity: {'#type': 'ViewInfo', id: 'b', viewStateMap: {}}},
        {id: 'r1', entity: {'#type': 'TableInfo', id: 'a', rowCount: 5}},
      ],
    };
    TYPES['Project'].strip!(p);
    expect(p.storage).toBeUndefined();
    expect(p.relations).toEqual([
      {id: 'r1', entity: {'#type': 'EntityRecord', id: 'a'}, isLink: false},
      {id: 'r2', entity: {'#type': 'EntityRecord', id: 'b'}, isLink: false},
    ]);
  });

  it('blanks group parents and children', () => {
    const g: any = {parents: [{parent: {id: 'p'}}], children: [{child: {id: 'c'}}]};
    TYPES['UserGroup'].strip!(g);
    expect(g).toEqual({parents: [], children: []});
  });
});

describe('files, jobs, notebooks and models', () => {
  it('reduces a file connection to its id', () => {
    const f: any = {connection: {id: 'c1', namespace: 'MigTestSpace:', name: 'Files', parameters: {isProject: true}}, tables: [{id: 't1'}]};
    TYPES['FileInfo'].strip!(f);
    expect(f.connection).toEqual({id: 'c1'});
    expect(f.tables).toBeUndefined();
    expect(TYPES['FileInfo'].deps!(f)).toEqual([{type: 'DataConnection', id: 'c1'}]);
  });

  it('lists files polymorphically — there is no type-wide /files listing', () => {
    expect(TYPES['FileInfo'].listVia).toBe('entities');
    expect(TYPES['FileInfo'].typeId).toBe('34d75630-e870-11e6-bfe1-590ff6f10d14');
  });

  it('writes file bytes after the save and table bytes before it', () => {
    expect(TYPES['FileInfo'].bytes!.afterSave).toBe(true);
    expect(TYPES['FileInfo'].bytes!.put('f1')).toBe('/files/data/f1');
    expect(TYPES['TableInfo'].bytes!.afterSave).toBeUndefined();
  });

  it('reduces notebook tables and a model trainedOn to ids, and depends on them', () => {
    const n: any = {tables: [{id: 't1', rowCount: 5}, {id: 't2'}]};
    TYPES['Notebook'].strip!(n);
    expect(n.tables).toEqual([{id: 't1'}, {id: 't2'}]);
    expect(TYPES['Notebook'].deps!(n)).toEqual([{type: 'TableInfo', id: 't1'}, {type: 'TableInfo', id: 't2'}]);

    const m: any = {trainedOn: {id: 't1', namespace: 'Admin:', rowCount: 5}};
    TYPES['PredictiveModelInfo'].strip!(m);
    expect(m.trainedOn).toEqual({id: 't1'});
    expect(TYPES['PredictiveModelInfo'].deps!(m)).toEqual([{type: 'TableInfo', id: 't1'}]);
  });

  it('has no reliable dependency for a job', () => {
    expect(TYPES['DataJob'].deps).toBeUndefined();
    expect(TYPES['DataJob'].route).toBe('/connectors/jobs');
  });

  it('turns the space and dashboard aliases into listing rules', () => {
    expect(typeOptionsOf('dashboard')).toEqual({clause: 'isDashboard = true'});
    expect(typeOptionsOf('space')).toEqual({
      clause: 'isDashboard = false and isEntity = false', params: {includeRoot: 'true'}});
    expect(typeOptionsOf('project')).toBeUndefined();
  });
});

describe('datasyncConnectionRefs', () => {
  it('finds the connection of every OpenFile-family call in the creation script', () => {
    const t = {metaParams: {'.data-sync': 'sync', '.script': 'OpenServerFile("Chem:Files/a.csv")\nOpenFolder(\'Bio:Data/x\')'}};
    expect(datasyncConnectionRefs(t)).toEqual([
      {type: 'DataConnection', nqName: 'Chem:Files'},
      {type: 'DataConnection', nqName: 'Bio:Data'},
    ]);
  });

  it('ignores tables that are not datasync', () => {
    expect(datasyncConnectionRefs({metaParams: {'.script': 'OpenFile("Chem:Files/a.csv")'}})).toEqual([]);
  });
});

describe('rankOf / resolveType', () => {
  it('ranks groups before connections before projects before queries before tables before views', () => {
    expect(rankOf('UserGroup')).toBeLessThan(rankOf('DataConnection'));
    expect(rankOf('DataConnection')).toBeLessThan(rankOf('Project'));
    expect(rankOf('Project')).toBeLessThan(rankOf('DataQuery'));
    expect(rankOf('DataQuery')).toBe(rankOf('Script'));
    expect(rankOf('Script')).toBeLessThan(rankOf('TableInfo'));
    expect(rankOf('TableInfo')).toBeLessThan(rankOf('ViewInfo'));
    expect(rankOf('ViewInfo')).toBe(rankOf('ViewLayout'));
    expect(rankOf('ViewLayout')).toBeLessThan(rankOf('FileInfo'));
    expect(rankOf('FileInfo')).toBeLessThan(rankOf('DataJob'));
    expect(rankOf('DataJob')).toBe(rankOf('Notebook'));
    expect(rankOf('Notebook')).toBe(rankOf('PredictiveModelInfo'));
  });

  it('maps the CLI aliases', () => {
    expect(resolveType('conn')).toBe('DataConnection');
    expect(resolveType('dashboard')).toBe('Project');
    expect(resolveType('space')).toBe('Project');
    expect(resolveType('Script')).toBe('Script');
    expect(resolveType('file')).toBe('FileInfo');
    expect(resolveType('job')).toBe('DataJob');
    expect(resolveType('notebook')).toBe('Notebook');
    expect(resolveType('model')).toBe('PredictiveModelInfo');
    expect(() => resolveType('nope')).toThrow(/not supported/);
  });
});

describe('nqNameOf', () => {
  it('joins the namespace with the name, falling back to the id', () => {
    expect(nqNameOf({namespace: 'Admin:', name: 'CerealDemog'})).toBe('Admin:CerealDemog');
    expect(nqNameOf({name: 'blob.txt'})).toBe('blob.txt');
    expect(nqNameOf({id: 'f1'})).toBe('f1');
  });
});

describe('fileNameFor', () => {
  it('flattens an nqName', () => {
    expect(fileNameFor({namespace: 'Chem:', name: 'Chembl'})).toBe('Chem.Chembl');
  });

  it('flattens a file path', () => {
    expect(fileNameFor({namespace: 'Admin:', name: 'Data/a/b.csv'})).toBe('Admin.Data.a.b.csv');
  });

  it('falls back to the id for unnamed entities', () => {
    expect(fileNameFor({id: 'f24a10d0'})).toBe('f24a10d0');
  });

  it('keeps a Windows device name out of the bundle', () => {
    expect(fileNameFor({name: 'NUL', id: 'abc'})).toBe('NUL-abc');
    expect(fileNameFor({name: 'com1.csv', path: 'com1.csv', id: 'abc'})).toBe('com1-abc.csv');
    expect(fileNameFor({namespace: 'Admin:', name: 'NUL'})).toBe('Admin.NUL');
  });

  it('names a namespace-less file by its full path, so two basenames do not collide', () => {
    expect(fileNameFor({name: 'b.csv', path: 'dir1/b.csv'})).toBe('dir1.b.csv');
    expect(fileNameFor({name: 'b.csv', path: 'dir2/b.csv'})).toBe('dir2.b.csv');
  });
});

describe('resolveTypes', () => {
  it('dedupes aliases of the same type', () => {
    expect(resolveTypes(['conn', 'connection']).types).toEqual(['DataConnection']);
    expect(resolveTypes(['dashboard', 'dashboard'])).toEqual({types: ['Project'], typeOptions: {Project: typeOptionsOf('dashboard')}});
  });

  it('refuses two aliases that would list the same type differently', () => {
    expect(() => resolveTypes(['space', 'dashboard'])).toThrow(/both select Project but list it differently/);
    expect(() => resolveTypes(['project', 'space'])).toThrow(/both select Project but list it differently/);
  });
});
