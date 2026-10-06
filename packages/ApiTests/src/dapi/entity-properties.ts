import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {thrown} from './domain-lifecycle';

// grok.meta — the curated property catalogs of platform entity types and of the
// read-only `Core` domain tables that serve their rows (entity-properties §10). The
// keystone rule under test is that ONE name serves both: what propertiesOf discovers is
// what the Core query route accepts (the drift probe) and what the platform's own grids
// and filter panels are built from.
category('Dapi: entity properties', () => {
  const SECRET = /pwd|password|token|salt|secret|hash/i;

  // The catalog of a type, failing here rather than at the first null dereference.
  async function catalog(type: string, options?: {filterable?: boolean}): Promise<_DG.EntityPropertyInfo[]> {
    const props = await grok.meta.propertiesOf(type, options);
    if (props == null)
      throw new Error(`${type} must have a catalog`);
    return props;
  }

  function prop(props: _DG.EntityPropertyInfo[], name: string): _DG.EntityPropertyInfo {
    const p = props.find((x) => x.name === name);
    if (p == null)
      throw new Error(`no ${name} among ${props.map((x) => x.name).join(', ')}`);
    return p;
  }

  async function location(type: string): Promise<_DG.CoreLocation> {
    const core = await grok.meta.coreLocationOf(type);
    if (core == null)
      throw new Error(`${type} must be served as a Core table`);
    return core;
  }

  const names = (props: _DG.EntityPropertyInfo[]) => props.map((p) => p.name);

  test('catalogs: present, ordered, labeled', async () => {
    expect(names(await catalog('User')).join(','),
      'login,firstName,lastName,status,createdOn,picture');
    expect(names(await catalog('UserSession')).join(','), 'user,started,ended,type');
    expect(names(await catalog('DataQuery')).join(','),
      'friendlyName,source,createdOn,updatedOn,description');

    // The label survived the joined -> createdOn harmonization.
    const joined = prop(await catalog('User'), 'createdOn');
    expect(joined.friendlyName, 'Joined');
    expect(joined.type, DG.TYPE.DATE_TIME);

    // Filterable is the database-backed subset: the picture is a rendering detail.
    const filterable = await catalog('User', {filterable: true});
    expect(names(filterable).join(','), 'login,firstName,lastName,status,createdOn');
    for (const p of filterable)
      expect(p.relationKind != null, true, `${p.name} must declare a relation kind`);

    // An uncurated type and an unknown one are the same answer, and neither is an error.
    expect(await grok.meta.propertiesOf('TableInfo'), null);
    expect(await grok.meta.propertiesOf('NoSuchTypeAtAll'), null);
  });

  test('references carry the type they point at and how they relate', async () => {
    const session = await catalog('UserSession');
    const user = prop(session, 'user');
    expect(user.refType, 'User');
    expect(user.relationKind, DG.RELATION_KIND.BELONGS_TO);
    expect(prop(session, 'started').refType, null);
    expect(prop(session, 'started').relationKind, DG.RELATION_KIND.FIELD);
  });

  test('a plugin domain table answers the same shape', async () => {
    const props = await catalog('apitests.item');
    expect(props.some((p) => p.name === 'sku'), true, 'the package schema must serve its columns');
    for (const p of props) {
      expect(typeof p.name, 'string');
      expect(typeof p.type, 'string');
      // A domain column is not an entity reference and carries no db annotation.
      expect(p.relationKind, null, p.name);
    }
    // Every domain column is filterable, so the flag narrows nothing here.
    expect(names(await catalog('apitests.item', {filterable: true})).join(','), names(props).join(','));
    // A table that does not exist misses exactly the way an uncurated type does.
    expect(await grok.meta.propertiesOf('apitests.nosuchtable'), null);
  });

  test('no catalog leaks a secret', async () => {
    // The all-provider sweep lives Dart-side (grok_shared/test/entity_catalog_test.dart,
    // 'no catalog exposes a secret'); this pins the same rule on the JS surface.
    const types = ['User', 'UserGroup', 'UserSession', 'DataConnection', 'DataQuery', 'Script',
      'Project', 'FileInfo', 'DataJob', 'FuncCall', 'Notebook', 'ViewLayout',
      'Core.users', 'Core.groups', 'Core.users_sessions'];
    for (const type of types)
      for (const p of (await grok.meta.propertiesOf(type)) ?? [])
        expect(SECRET.test(p.name), false, `${type}.${p.name} looks like a secret`);
    // Explicit, so a rename of the fields cannot quietly pass the pattern above.
    for (const n of ['pwdHash', 'pwdSalt', 'tokenHash', 'externalToken', 'token'])
      expect(names(await catalog('User')).includes(n), false, n);
  });

  test('coreLocationOf points at the registered table, null for the rest', async () => {
    expect(JSON.stringify(await location('User')), JSON.stringify({schema: 'Core', table: 'users'}));
    expect(JSON.stringify(await location('UserSession')),
      JSON.stringify({schema: 'Core', table: 'users_sessions'}));
    // Deferred registrations and uncurated types answer null, never a guess
    // (projects: ARCHITECTURE §7.1 ruling).
    expect(await grok.meta.coreLocationOf('Project'), null);
    expect(await grok.meta.coreLocationOf('NoSuchTypeAtAll'), null);
    // A registered provider that is not database-backed at all.
    expect(await grok.meta.coreLocationOf('Property'), null);
  });

  test('descriptor names match the Core manifest', async () => {
    for (const type of ['User', 'UserSession']) {
      const core = await location(type);
      const served = names(await catalog(`${core.schema}.${core.table}`));
      for (const p of await catalog(type))
        expect(served.includes(p.name), true,
          `${type}.${p.name} is not served as ${core.schema}.${core.table}.${p.name}`);
    }
  });

  test('drift probe: every filterable User property filters through Core', async () => {
    const users = grok.dapi.domains.table('Core.users');
    const total = await users.count();
    expect(total > 0, true, 'the Core users table must serve rows');
    for (const p of await catalog('User', {filterable: true})) {
      const count = await users.count({property: p.name, operator: 'is not', value: null});
      expect(count >= 0 && count <= total, true, `${p.name} filtered to ${count} of ${total}`);
      // login is NOT NULL for every user, so a route that accepted the discovered name
      // without actually applying it would still pass the range check above.
      if (p.name === 'login')
        expect(count, total, 'every user has a login');
    }
  });

  test('an undiscoverable column is refused with no oracle', async () => {
    const users = grok.dapi.domains.table('Core.users');
    const secret = await thrown(() => users.count({property: 'pwdHash', operator: 'is not', value: null}));
    const absent = await thrown(() => users.count({property: 'nosuchcolumn', operator: 'is not', value: null}));
    expect(secret instanceof DG.DomainError, true, `${secret?.constructor?.name}: ${secret?.message}`);
    expect(absent instanceof DG.DomainError, true, `${absent?.constructor?.name}: ${absent?.message}`);
    expect(secret.code, absent.code);
    expect(secret.message.replace('pwdHash', '#'), absent.message.replace('nosuchcolumn', '#'),
      'a secret column must be indistinguishable from a column that does not exist');
  });

  test('session -> user round trip agrees with dinq', async () => {
    const me = await grok.dapi.users.current();
    const sessions = grok.dapi.domains.table('Core.users_sessions');
    const hop = (await catalog('UserSession')).find((p) => p.refType === 'User');
    if (hop == null)
      throw new Error('UserSession must carry a User reference');
    const key = prop(await catalog('User', {filterable: true}), 'login').name;

    // The same discovered name serves dinq's `text` filter and the domain compiler.
    const viaDinq = await grok.dapi.users.filter(`${key} = "${me.login}"`).list();
    expect(viaDinq.length, 1, `dinq must resolve exactly one user by ${key}`);
    expect(await grok.dapi.domains.table('Core.users').count(
      {property: key, operator: '=', value: me.login}), 1);

    // ...and through the session's ref hop, where the Core registration is the
    // only thing that knows `user` points at Core.users.
    const count = async (filter: _DG.DomainFilter) =>
      (await sessions.facets({filter: filter, facets: [{id: 'n', kind: 'count'}]})).facets.n.count;
    const byLogin = await count([{property: `${hop.name}.${key}`, operator: '=', value: me.login}]);
    expect(byLogin, await count([{property: `${hop.name}.id`, operator: '=', value: me.id}]),
      'the discovered login hop and the id hop must select the same sessions');
    expect(byLogin > 0, true, 'the caller has at least the current session');
  });

  test('Core writes are refused with a structured 403', async () => {
    const login = `core-write-probe-${Date.now()}`;
    try {
      const err = await thrown(() => grok.dapi.domains.table('Core.users').insert({login: login}));
      expect(err instanceof DG.DomainError, true, `${err?.constructor?.name}: ${err?.message}`);
      expect(err.status, 403);
      expect(err.code, 'read-only');
    } finally {
      for (const u of await grok.dapi.users.filter(`login = "${login}"`).list())
        await grok.dapi.users.delete(u);
    }
  });
});
