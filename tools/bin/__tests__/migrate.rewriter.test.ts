import {describe, it, expect} from 'vitest';
import {BundleEntity} from '../utils/migrate/bundle';
import {nestedIds, rewrite} from '../utils/migrate/rewriter';
import {reId} from './migrate-helpers';

const CONN = 'bbbbbbbb-1111-2222-3333-444444444444';
const CONN2 = 'bbbbbbbb-9999-2222-3333-444444444444';
const TABLE = 'f2444470-6fc0-11f1-b307-bdcf313f748a';
const PROJECT = 'efe0b1f0-6fc0-11f1-b275-83ec2160b5e9';
const QUERY = 'aaaaaaaa-1111-2222-3333-444444444444';
const VIEW = 'f24a10d0-6fc0-11f1-e396-fd8fa45dfcf4';
const REL = 'cccccccc-1111-2222-3333-444444444444';
const GROUP = 'dddddddd-1111-2222-3333-444444444444';
const UNKNOWN = '11111111-2222-3333-4444-555555555555';

describe('rewrite', () => {
  it('replaces mapped ids in values, arrays, typed strings and keys', () => {
    const json = {
      id: CONN,
      list: [{connection: {id: CONN}}, `DataConnection:${CONN}`],
      byId: {[`DataConnection:${CONN}`]: 'x', [CONN]: 'y'},
      untouched: 'DataConnection:not-a-uuid',
    };
    const out: any = rewrite(json, {[CONN]: CONN2});
    expect(out.id).toBe(CONN2);
    expect(out.list[0].connection.id).toBe(CONN2);
    expect(out.list[1]).toBe(`DataConnection:${CONN2}`);
    expect(Object.keys(out.byId)).toEqual([`DataConnection:${CONN2}`, CONN2]);
    expect(out.untouched).toBe('DataConnection:not-a-uuid');
    expect(json.id).toBe(CONN);
  });

  it('leaves an unmapped id in place and reports it as an orphan', () => {
    const orphans = new Set<string>();
    const out: any = rewrite({id: CONN, table: {id: UNKNOWN}}, {[CONN]: CONN2}, orphans);
    expect(out.table.id).toBe(UNKNOWN);
    expect([...orphans]).toEqual([UNKNOWN]);
  });
});

describe('nestedIds', () => {
  it('collects rows that own a primary key and skips the entity itself and plain references', () => {
    const project = {
      id: PROJECT,
      connection: {id: CONN},
      relations: [{id: REL, entity: {'#type': 'EntityRecord', id: TABLE}, isLink: false}],
      params: [{id: QUERY, name: 'x'}],
    };
    expect(nestedIds(project).sort()).toEqual([QUERY, REL].sort());
  });
});

describe('reId', () => {
  const bundle = (): Map<string, BundleEntity> => new Map<string, BundleEntity>([
    [CONN, {type: 'DataConnection', json: {'#type': 'DataConnection', id: CONN, name: 'MigTestConn'}}],
    [PROJECT, {type: 'Project', json: {
      '#type': 'Project', id: PROJECT, name: 'MigTestDash', _grants: [{group: 'MigTestGroup', permission: 'View'}],
      relations: [{id: REL, entity: {'#type': 'EntityRecord', id: TABLE}, isLink: false}],
    }}],
    [TABLE, {type: 'TableInfo', json: {'#type': 'TableInfo', id: TABLE, name: 'MigTestTable'}}],
    [QUERY, {type: 'DataQuery', json: {'#type': 'DataQuery', id: QUERY, name: 'MigTestQuery', connection: {id: CONN}}}],
    [VIEW, {type: 'ViewInfo', json: {'#type': 'ViewInfo', id: VIEW, name: 'view', table: {id: TABLE}}}],
    [GROUP, {type: 'UserGroup', json: {
      '#type': 'UserGroup', id: GROUP, name: 'MigTestGroup', friendlyName: 'MigTestGroup',
      _members: [{kind: 'user', login: 'admin'}, {kind: 'group', name: 'MigTestChild'}],
    }}],
  ]);

  it('mints a fresh identity for every entity and keeps the references consistent', () => {
    const {entities, idmap} = reId(bundle());
    expect([...entities.keys()].sort()).toEqual([CONN, PROJECT, TABLE, QUERY, VIEW, GROUP].map((id) => idmap[id]).sort());

    const query = entities.get(idmap[QUERY])!.json;
    expect(query.id).toBe(idmap[QUERY]);
    expect(query.connection.id).toBe(idmap[CONN]);
    expect(entities.get(idmap[VIEW])!.json.table.id).toBe(idmap[TABLE]);
    const project = entities.get(idmap[PROJECT])!.json;
    expect(project.relations[0].entity.id).toBe(idmap[TABLE]);
  });

  it('mints fresh ids for nested rows so the source rows are never re-pointed', () => {
    const {entities, idmap} = reId(bundle());
    const relation = entities.get(idmap[PROJECT])!.json.relations[0];
    expect(relation.id).toBe(idmap[REL]);
    expect(relation.id).not.toBe(REL);
  });

  it('renames the entity and the group names grants and memberships refer to', () => {
    const {entities, idmap} = reId(bundle(), (name) => name.replace('MigTest', 'MigTestCopy'));
    expect(entities.get(idmap[CONN])!.json.name).toBe('MigTestCopyConn');
    const group = entities.get(idmap[GROUP])!.json;
    expect(group.friendlyName).toBe('MigTestCopyGroup');
    expect(group._members).toEqual([{kind: 'user', login: 'admin'}, {kind: 'group', name: 'MigTestCopyChild'}]);
    expect(entities.get(idmap[PROJECT])!.json._grants).toEqual([{group: 'MigTestCopyGroup', permission: 'View'}]);
  });

  it('renames the namespace segments, so a dependency follows its renamed owner', () => {
    const source = bundle();
    source.get(TABLE)!.json.namespace = 'Admin:MigTestDash:';
    source.set('f1', {type: 'FileInfo', json: {'#type': 'FileInfo', id: 'f1', name: 'a.csv', path: 'dir/MigTestA.csv'}});
    const {entities, idmap} = reId(source, (name) => name.replace('MigTest', 'MigTestCopy'));
    expect(entities.get(idmap[TABLE])!.json.namespace).toBe('Admin:MigTestCopyDash:');
    expect(entities.get(idmap['f1'])!.json.path).toBe('dir/MigTestCopyA.csv');
  });
});
