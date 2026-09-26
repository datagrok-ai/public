import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {thrown, withRestrictedUser} from './domain-lifecycle';

// Editing a live external binding from a script: the access snapshot, access deltas riding the
// apply (a dry run, then the commit), the two concurrency tokens, and the schema-scoped draft an
// editor without CreateDomainSchema gets. Rides a throwaway binding cloned from the `extwlive`
// fixture's manifest — the scratch tables in test_db, bound read-only here, so nothing below
// writes them; Northwind is never touched. Skips cleanly without the fixture or the privilege.
category('Dapi: domain editing', () => {
  const fixture = 'extwlive';
  const name = `zzed${`${Date.now()}`.slice(-8)}`;
  const handle = grok.dapi.domains.schema(name);
  const allUsers = () => grok.dapi.groups.filter('friendlyName = "All users"').first();
  let manifest: _DG.DomainManifest;
  let skip: string | null = null;

  const tokens = async (): Promise<{ifVersion: string; ifIncarnation: string}> => {
    const m = await handle.manifest();
    return {ifVersion: m.version, ifIncarnation: m.incarnation};
  };

  before(async () => {
    if (!(await grok.dapi.domains.schemas.list()).some((s) => s.name === fixture)) {
      skip = `the ${fixture} fixture is not registered`;
      return;
    }
    const m: any = await grok.dapi.domains.schema(fixture).manifest();
    delete m['incarnation'];
    delete m['storage']['writable'];
    for (const t of Object.values(m['tables'] as {[name: string]: any}))
      delete t['writable'];
    manifest = m;
    try {
      await grok.dapi.domains.createSchema(name, {friendlyName: 'Editing probe', manifest});
    } catch (e: any) {
      if (e instanceof DG.DomainError && (e.code === 'forbidden' || e.status === 403))
        skip = 'no CreateDomainSchema privilege';
      else
        throw e;
    }
  });

  after(async () => {
    if (skip == null)
      await handle.delete().catch((e) => console.log(`cleanup of ${name} failed: ${e}`));
  });

  const skipped = (): boolean => {
    if (skip != null)
      console.log(`skipped: ${skip}`);
    return skip != null;
  };

  test('the access snapshot: tables with their direct grants, unrestricted columns, the manifest tokens', async () => {
    if (skipped())
      return;
    const [m, access] = await Promise.all([handle.manifest(), handle.access()]);
    expect(access.version, m.version);
    expect(access.incarnation, m.incarnation);
    expect(Object.keys(access.tables).sort().join(','), Object.keys(m.tables).sort().join(','));
    const thing = access.tables.thing;
    expect(thing.remote, 'extw_live_orders');
    expect(thing.entityId != null, true, JSON.stringify(thing));
    expect(thing.canShare, true);
    // The creator's bootstrap: the full set on the table, in one complete permission list.
    expect(thing.grants!.some((g) => ['View', 'Edit', 'Delete', 'Share'].every((p) => g.permissions.includes(p))),
      true, JSON.stringify(thing.grants));
    expect(thing.coreSchema.canShare, true);
    expect(thing.coreSchema.grants != null, true, 'the core schema grants are readable with Share on it');
    const secret = access.columns['thing.secret'];
    expect(secret.state, 'unrestricted');
    expect(secret.canShare, true);
    expect(secret.schemaId == null && secret.view == null, true, JSON.stringify(secret));
    expect(Object.keys(access.columns).filter((c) => c.startsWith('thing.')).length,
      Object.keys(m.tables.thing.columns).length);
  });

  test('access in the apply: a dry run lists every op with its effect and moves nothing; a commit lands', async () => {
    if (skipped())
      return;
    const group = (await allUsers()).id;
    const delta: _DG.DomainAccessDelta = {
      grant: [{table: 'thing', group, permission: 'View'}],
      revoke: [{table: 'thing', group, permission: 'Edit'}],
      restrict: [{table: 'thing', column: 'secret', grant: [{group, permission: 'View'}]}],
    };
    const initial = await handle.access();
    const plan = await handle.apply({...await tokens(), access: delta}, {dryRun: true});
    expect(plan.registrationOnly, true, JSON.stringify(plan));
    expect(plan.applied == null, true, 'a dry run is not applied');
    expect(plan.access!.grant[0].effect, 'grant');
    expect(plan.access!.grant[0].group.id, group);
    expect(plan.access!.revoke[0].effect, 'none', 'nothing to revoke: the group held no Edit');
    expect(plan.access!.restrict[0].effect, 'restrict');
    expect(plan.access!.restrict[0].grants[0].effect, 'grant');
    const unmoved = await handle.access();
    expect(unmoved.version, initial.version, 'a dry run moves no token');
    expect(unmoved.columns['thing.secret'].state, 'unrestricted', 'a dry run restricts nothing');

    const res = await handle.apply({...await tokens(), access: delta});
    expect(res.applied, true);
    expect(res.access!.grant[0].effect, 'grant');
    expect(res.access!.restrict[0].effect, 'restrict');
    const applied = await handle.access();
    expect(Number(applied.version), Number(initial.version) + 1);
    expect(applied.tables.thing.grants!.find((g) => g.group.id === group)?.permissions.join(','), 'View');
    const secret = applied.columns['thing.secret'];
    expect(secret.state, 'restricted');
    expect(secret.view!.map((g) => g.id).includes(group), true, JSON.stringify(secret));
    expect(secret.edit!.length, 0);

    // The same delta again: every op a no-op, so the apply writes nothing.
    const again = await handle.apply({...await tokens(), access: delta});
    expect(again.applied, false);
    expect(again.noop, true);
    expect(again.access!.grant[0].effect, 'none');
    expect(again.access!.restrict[0].effect, 'none');
    expect(again.access!.restrict[0].grants[0].effect, 'none');

    const back = await handle.apply({...await tokens(), access: {unrestrict: [{table: 'thing', column: 'secret'}]}});
    expect(back.access!.unrestrict[0].effect, 'unrestrict');
    expect((await handle.access()).columns['thing.secret'].state, 'unrestricted');
  });

  test('a stale version token is a typed conflict naming both versions', async () => {
    if (skipped())
      return;
    const t = await tokens();
    await handle.apply({...t, description: 'moved'});
    const stale = await thrown(() => handle.apply({...t, description: 'never'}));
    expect(stale instanceof DG.DomainVersionConflictError, true,
      `expected DomainVersionConflictError, got ${stale?.constructor?.name}: ${stale?.message}`);
    expect(stale.code, 'version-conflict');
    expect(stale.expectedVersion, t.ifVersion);
    expect(stale.currentVersion, `${Number(t.ifVersion) + 1}`);
    expect((await grok.dapi.domains.schemas.list()).find((s) => s.name === name)!.description, 'moved');
  });

  test('the incarnation travels with the version: after a delete and a re-create the old token conflicts', async () => {
    if (skipped())
      return;
    const old = await handle.manifest();
    await handle.delete();
    const created = await grok.dapi.domains.createSchema(name, {friendlyName: 'Editing probe', manifest});
    expect(created.version, '1');
    expect(created.incarnation !== old.incarnation, true, 'a re-created name has a new incarnation');
    const stale = await thrown(() =>
      handle.apply({ifVersion: '1', ifIncarnation: old.incarnation, description: 'stale'}));
    expect(stale instanceof DG.DomainVersionConflictError, true,
      `expected DomainVersionConflictError, got ${stale?.constructor?.name}: ${stale?.message}`);
    expect(stale.expectedIncarnation, old.incarnation);
    expect(stale.currentIncarnation, created.incarnation);
    const fresh = await handle.apply({ifVersion: '1', ifIncarnation: created.incarnation, description: 'fresh'});
    expect(fresh.applied, true);
    expect((await handle.manifest()).incarnation, created.incarnation);
  });

  test('a delete pinned to an incarnation is a conflict once the name was re-created, and keeps the schema', async () => {
    if (skipped())
      return;
    const old = (await handle.manifest()).incarnation;
    await handle.delete({ifIncarnation: old});
    const created = await grok.dapi.domains.createSchema(name, {friendlyName: 'Editing probe', manifest});
    const stale = await thrown(() => handle.delete({ifIncarnation: old}));
    expect(stale instanceof DG.DomainVersionConflictError, true,
      `expected DomainVersionConflictError, got ${stale?.constructor?.name}: ${stale?.message}`);
    expect(stale.status, 409);
    expect(stale.expectedIncarnation, old);
    expect(stale.currentIncarnation, created.incarnation);
    expect((await handle.manifest()).incarnation, created.incarnation, 'the schema re-created since is kept');
  });

  test('an unrestriction carries the ACL revision it was made from: a revoke since is a conflict', async () => {
    if (skipped())
      return;
    const group = (await allUsers()).id;
    await handle.apply({...await tokens(), access: {restrict: [{table: 'thing', column: 'secret',
      grant: [{group, permission: 'View'}]}]}});
    const loaded = (await handle.access()).columns['thing.secret'];
    expect(loaded.state, 'restricted');
    expect(typeof loaded.revision, 'string', JSON.stringify(loaded));
    await handle.apply({...await tokens(), access: {restrict: [{table: 'thing', column: 'secret',
      revoke: [{group, permission: 'View'}]}]}});
    const current = (await handle.access()).columns['thing.secret'];
    expect(current.revision !== loaded.revision, true, 'a revoke moves the revision');
    const t = await tokens();
    const stale = await thrown(() => handle.apply({...t, access: {unrestrict: [{table: 'thing', column: 'secret',
      from: 'restricted', revision: loaded.revision}]}}));
    expect(stale instanceof DG.DomainError, true, `expected DomainError, got ${stale?.constructor?.name}: ${stale?.message}`);
    expect(stale.status, 409);
    expect(stale.code, 'access-conflict');
    expect((await handle.access()).columns['thing.secret'].state, 'restricted', 'the stale unrestriction wrote nothing');
    const back = await handle.apply({...t, access: {unrestrict: [{table: 'thing', column: 'secret',
      from: 'restricted', revision: current.revision}]}});
    expect(back.access!.unrestrict[0].effect, 'unrestrict');
  });

  test('an editor without Share on a table reads its grants as unknown and cannot change them', async () => {
    if (skipped())
      return;
    const group = (await allUsers()).id;
    await withRestrictedUser('edit', async (probe) => {
      await handle.grant(probe.group, 'Edit');
      const snapshot = await probe.asUser(() => handle.access());
      expect(snapshot.tables.thing.canShare, false);
      expect(snapshot.tables.thing.grants, null);
      expect(snapshot.columns['thing.note'].canShare, false);
      const t = await tokens();
      const err = await thrown(() => probe.asUser(() => handle.apply({...t,
        access: {grant: [{table: 'thing', group, permission: 'View'}]}})));
      expect(err instanceof DG.DomainError, true,
        `expected DomainError, got ${err?.constructor?.name}: ${err?.message}`);
      expect(err.status, 403);
      expect(err.code, 'access-forbidden');
      expect(err.body['targets'].join(','), 'thing');
      // Edit alone carries a metadata change — and the apply grants the editor nothing.
      const ok = await probe.asUser(() => handle.apply({...t, description: 'by the editor'}));
      expect(ok.applied, true);
      expect((await handle.access()).tables.thing.grants!.some((g) => g.group.id === probe.group), false);
    });
  });

  test('the schema-scoped draft: schema Edit plus the connection privileges, no CreateDomainSchema', async () => {
    if (skipped())
      return;
    // The smart filter's `name` does not match a plugin connection; free text over the short
    // name does, and the nqName picks the one.
    const nqName = manifest.storage!.connection!;
    const connection = (await grok.dapi.connections.filter(nqName.split(':')[1]).list())
      .find((c) => c.nqName === nqName);
    expect(connection != null, true, `connection ${nqName} not found`);
    await withRestrictedUser('draft', async (probe) => {
      expect(await probe.asUser(() => grok.dapi.permissions.checkGlobal(DG.Permission.CREATE_DOMAIN_SCHEMA)), false);
      await handle.grant(probe.group, 'Edit');
      // A stand may hand every user the connection through "All users"; the invisible leg
      // is asserted only where the probe really cannot see it, and the share is made only
      // where a privilege is missing.
      const needed = ['View', DG.Permission.DATA_CONNECTION_GET_SCHEMA, DG.Permission.DATA_CONNECTION_QUERY];
      const held: string[] = [];
      for (const p of needed) {
        if (await probe.asUser(() => grok.dapi.permissions.check(connection!, p)))
          held.push(p);
      }
      if (held.length === 0) {
        // The connection is invisible to the editor until it is shared: no oracle.
        const invisible = await thrown(() => probe.asUser(() => handle.draft()));
        expect(invisible instanceof DG.DomainError, true, `${invisible?.constructor?.name}: ${invisible?.message}`);
        expect(invisible.code, 'unknown-connection');
      } else
        console.log(`the probe already holds ${held.join(', ')} on ${nqName}: the unknown-connection leg is skipped`);

      const group = await grok.dapi.groups.find(probe.group);
      const share = held.length < needed.length;
      if (share)
        await grok.dapi.permissions.grant(connection!, group, false);
      try {
        const draft = await probe.asUser(() => handle.draft({tables: ['extw_live_orders']}));
        expect(Object.keys(draft.manifest.tables).join(','), 'extw_live_orders', JSON.stringify(draft.manifest));
        expect(draft.manifest.storage!.connection, nqName);
        expect(draft.inventory.tables.some((t) => t.remote === 'extw_live_parent'), true,
          'the inventory covers the whole remote schema');
        await handle.revoke(probe.group, 'Edit');
        const forbidden = await thrown(() => probe.asUser(() => handle.draft()));
        expect(forbidden instanceof DG.DomainForbiddenError, true,
          `${forbidden?.constructor?.name}: ${forbidden?.message}`);
      } finally {
        if (share)
          await grok.dapi.permissions.revoke(connection!, group);
      }
    });
  });
});
