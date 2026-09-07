import {describe, it, expect, beforeEach, afterEach} from 'vitest';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';
import {BundleEntity, hashOf, list, normalize, read, write} from '../utils/migrate/bundle';

let dir: string;

beforeEach(() => { dir = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-bundle-')); });
afterEach(() => fs.rmSync(dir, {recursive: true, force: true}));

const SOURCE = {url: 'http://h/api', version: '1.27.9', userNamespace: 'Admin:'};

function conn(id: string, name: string, extra: any = {}): [string, BundleEntity] {
  return [id, {type: 'DataConnection', json: {'#type': 'DataConnection', id, name, namespace: 'Admin:', ...extra}}];
}

function put(entities: [string, BundleEntity][], opts: {replace?: boolean} = {}, bytes = new Map<string, Buffer>()) {
  return write(dir, new Map(entities), {source: SOURCE, args: ['pull'], packages: []}, opts, bytes);
}

describe('bundle.write', () => {
  it('merges consecutive pulls: the newer copy wins, nothing is removed', () => {
    put([conn('c1', 'NW'), conn('c2', 'PG')]);
    const manifest = put([conn('c1', 'NW', {dataSource: 'PostgresDart'}), conn('c3', 'Files')]);

    expect(manifest.pulls.length).toBe(2);
    expect(manifest.order.map((e) => e.id).sort()).toEqual(['c1', 'c2', 'c3']);
    expect(JSON.parse(fs.readFileSync(path.join(dir, 'DataConnection/Admin.NW.json'), 'utf8')).dataSource).toBe('PostgresDart');
    expect(fs.existsSync(path.join(dir, 'DataConnection/Admin.PG.json'))).toBe(true);
  });

  it('gives two same-named files their own bundle file, paired with the right id', () => {
    const file = (id: string, dir: string): [string, BundleEntity] => [id, {type: 'FileInfo', json:
      {'#type': 'FileInfo', id, name: 'b.csv', path: `${dir}/b.csv`, isFile: true}}];
    const manifest = put([file('f1', 'one'), file('f2', 'one'), file('f3', 'two')]);
    const files = manifest.order.map((e) => e.file);
    expect(new Set(files).size).toBe(3);
    expect(files).toContain('FileInfo/one.b.csv.json');
    expect(files).toContain('FileInfo/two.b.csv.json');

    const {entities} = read(dir);
    for (const [id, e] of entities)
      expect(e.json.id).toBe(id);
  });

  it('orders by rank, then by the pull that first brought the entity in', () => {
    put([['p1', {type: 'Project', json: {id: 'p1', name: 'Dash', namespace: 'Admin:'}}]]);
    const manifest = put([conn('c1', 'NW'), ['t1', {type: 'TableInfo', json: {id: 't1', name: 'T', namespace: 'Admin:'}}]]);
    expect(manifest.order.map((e) => e.type)).toEqual(['DataConnection', 'Project', 'TableInfo']);
  });

  it('moves a renamed entity to its new file and deletes the stale one', () => {
    put([conn('c1', 'NW')]);
    put([conn('c1', 'Northwind')]);
    expect(fs.existsSync(path.join(dir, 'DataConnection/Admin.NW.json'))).toBe(false);
    expect(fs.existsSync(path.join(dir, 'DataConnection/Admin.Northwind.json'))).toBe(true);
    expect(read(dir).entities.size).toBe(1);
  });

  it('--replace wipes the directory first', () => {
    put([conn('c1', 'NW')]);
    const manifest = put([conn('c2', 'PG')], {replace: true});
    expect(manifest.pulls.length).toBe(1);
    expect(fs.existsSync(path.join(dir, 'DataConnection/Admin.NW.json'))).toBe(false);
  });

  it('writes table bytes next to the metadata', () => {
    put([['t1', {type: 'TableInfo', json: {id: 't1', name: 'T', namespace: 'Admin:'}}]], {},
      new Map([['t1', Buffer.from([1, 2, 3])]]));
    expect(fs.readFileSync(path.join(dir, 'tables/t1.d42'))).toEqual(Buffer.from([1, 2, 3]));
  });
});

describe('bundle.read / list', () => {
  it('round-trips the entities and reports one row per manifest entry', () => {
    put([conn('c1', 'NW')]);
    const bundle = read(dir);
    expect(bundle.entities.get('c1')!.json.name).toBe('NW');
    expect(bundle.idmap).toEqual({});
    expect(list(dir)).toEqual([{
      type: 'DataConnection', nqName: 'Admin:NW', file: 'DataConnection/Admin.NW.json',
      pulledOn: bundle.manifest.pulls[0].at, pull: 1,
    }]);
  });

  it('refuses a directory that is not a bundle', () => {
    expect(() => read(path.join(dir, 'nope'))).toThrow(/no manifest.json/);
  });
});

describe('normalize / hashOf', () => {
  const base = {'#type': 'Script', id: 's1', name: 'S', language: 'python', script: 'print(1)'};

  it('is insensitive to key order', () => {
    expect(hashOf('Script', {...base})).toBe(hashOf('Script', {script: 'print(1)', language: 'python', name: 'S', id: 's1', '#type': 'Script'}));
  });

  it('ignores volatile fields', () => {
    expect(hashOf('Script', {...base, createdOn: 'a', updatedOn: 'b', author: {id: 'u1'}, pictureId: 'p'})).toBe(hashOf('Script', base));
  });

  it('ignores bundle-only keys but keeps them in the written file', () => {
    const masked = {'#type': 'DataConnection', id: 'c1', name: 'NW', parameters: {password: '_____________', db: 'nw'}};
    expect(normalize('DataConnection', masked)._credentials).toEqual(['password']);
    expect(hashOf('DataConnection', masked)).toBe(hashOf('DataConnection', {'#type': 'DataConnection', id: 'c1', name: 'NW', parameters: {db: 'nw'}}));
  });

  it('changes when real content changes', () => {
    expect(hashOf('Script', {...base, script: 'print(2)'})).not.toBe(hashOf('Script', base));
  });

  it('ignores the namespace the target computes for itself, but writes it to the file', () => {
    expect(hashOf('Script', {...base, namespace: 'Askalkin:'})).toBe(hashOf('Script', {...base, namespace: 'Admin:'}));
    put([['s1', {type: 'Script', json: {...base, namespace: 'Askalkin:'}}]]);
    expect(JSON.parse(fs.readFileSync(path.join(dir, 'Script/Askalkin.S.json'), 'utf8')).namespace).toBe('Askalkin:');
  });
});
