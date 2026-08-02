import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';

// WO-4a (GROK-20601): schema lifecycle handle (grok.dapi.domains.schema(name):
// manifest/apply/audit/delete + createSchema) and table-scoped admin surface
// (grants/grant/revoke, shareColumn/restrictColumn/restoreColumnVisibility).
// The lifecycle round-trip needs the CreateDomainSchema privilege and skips
// cleanly without it; grants/column tests ride the package's own 'apitests'
// schema and grant only to 'All users' with finally-cleanup.
category('Dapi: domain lifecycle', () => {
  async function thrown(action: () => Promise<any>): Promise<any> {
    try {
      await action();
      return null;
    } catch (e) {
      return e;
    }
  }

  const allUsers = async () => (await grok.dapi.groups.filter('friendlyName = "All users"').first());

  test('schema lifecycle: create → apply (dryRun first) → data → audit → typed errors → delete', async () => {
    const name = `zzlc${`${Date.now()}`.slice(-8)}`;
    try {
      await grok.dapi.domains.createSchema(name, {friendlyName: 'Lifecycle probe'});
    } catch (e: any) {
      if (e instanceof DG.DomainError && (e.code === 'forbidden' || e.status === 403)) {
        console.log('skipped: no CreateDomainSchema privilege');
        return;
      }
      throw e;
    }
    const handle = grok.dapi.domains.schema(name);
    try {
      const tables = {things: {columns: {name: {type: 'string', required: true}, qty: {type: 'int'}}}};
      const plan = await handle.apply({tables}, {dryRun: true});
      expect(plan != null, true, 'dryRun must return the change plan');
      await handle.apply({tables});

      const manifest = await handle.manifest();
      expect(manifest['tables']?.['things'] != null, true, JSON.stringify(manifest));

      // Data flows through the standard table client immediately after apply.
      const things = grok.dapi.domains.table(`${name}.things`);
      const [ins] = await things.insert({name: 'w1', qty: 2});
      expect(ins.created, true);
      expect((await things.query({filter: 'name = "w1"'})).length, 1);

      // Whole-schema audit: ddl registry events + the row insert, newest first.
      const audit = await handle.audit({limit: 100});
      expect(audit.some((a) => a.op === 'ddl'), true, 'ddl events expected');
      expect(audit.some((a) => a.op === 'insert'), true, 'row events expected');

      // Stale ifVersion → typed version conflict.
      const stale = await thrown(() => handle.apply({tables, ifVersion: '999'}));
      expect(stale instanceof DG.DomainVersionConflictError, true,
        `expected DomainVersionConflictError, got ${stale?.constructor?.name}: ${stale?.message}`);

      // Destructive apply without confirmation → typed error carrying the plan.
      const destructive = await thrown(() => handle.apply({dropTables: ['things']}));
      expect(destructive instanceof DG.DomainError, true,
        `expected DomainError, got ${destructive?.constructor?.name}: ${destructive?.message}`);
      expect(destructive.code, 'destructive-confirmation-required');
      expect(destructive.body['plan'] != null, true, JSON.stringify(destructive?.body));
    } finally {
      await handle.delete();
    }
    const schemas = await grok.dapi.domains.schemas.list();
    expect(schemas.some((s) => s.name === name), false, 'deleted schema must not be listed');
  });

  test('table grants: list, grant, revoke round-trip', async () => {
    const items = grok.dapi.domains.table('apitests.item');
    const group = await allUsers();
    const before = (await items.grants()).length;
    await items.grant(group.id, 'View');
    try {
      const listed = await items.grants();
      expect(listed.some((g) => g.group.id === group.id && g.permission === 'View'), true,
        JSON.stringify(listed));
    } finally {
      await items.revoke(group.id, 'View');
    }
    const after = await items.grants();
    expect(after.some((g) => g.group.id === group.id && g.permission === 'View'), false,
      'revoked grant must disappear from a second read');
    expect(after.length, before);
  });

  test('column security: shareColumn, idempotent restrictColumn, restoreColumnVisibility', async () => {
    const items = grok.dapi.domains.table('apitests.item');
    const group = await allUsers();
    // Sharing to All users first means a mid-test failure never hides the column.
    const shared = await items.shareColumn('name', group.id);
    try {
      expect(shared.id != null, true, JSON.stringify(shared));
      expect(shared.name, 'apitests.item.name');
      // The grant half really landed. Per-column property schemas have NO
      // entities rows (registry-entities rule), so grok.dapi.getEntities /
      // permissions cannot read them and no typed JS accessor takes arbitrary
      // registry ids (the WO-4a security boundary) — read the domains-native
      // grants route directly, the ApiTests authenticated-fetch pattern
      // (cf. user-settings-storage.ts).
      const res = await fetch(`${grok.dapi.root}/domains/grants/${shared.id}`,
        {credentials: 'include'});
      const granted = await res.json();
      expect(Array.isArray(granted), true, JSON.stringify(granted));
      expect(granted.some((g: any) => g.group?.id === group.id && g.permission === 'View'), true,
        'shareColumn must grant View on the per-column schema');
      const again = await items.restrictColumn('name');
      expect(again.id, shared.id, 'restrict is idempotent — same per-column schema');
    } finally {
      await items.restoreColumnVisibility('name');
    }
    // Restored: the core schema serves the column again — a query NAMING it
    // succeeding is the assertion (a restricted column would 400, no-oracle).
    const rows = await items.query({limit: 1, columns: ['name']});
    expect(Array.isArray(rows), true);
  });

  test('schema grants: list, grant, revoke round-trip', async () => {
    const schema = grok.dapi.domains.schema('apitests');
    const group = await allUsers();
    const before = (await schema.grants()).length;
    await schema.grant(group.id, 'View');
    try {
      const listed = await schema.grants();
      expect(listed.some((g) => g.group.id === group.id && g.permission === 'View'), true,
        JSON.stringify(listed));
    } finally {
      await schema.revoke(group.id, 'View');
    }
    const after = await schema.grants();
    expect(after.some((g) => g.group.id === group.id && g.permission === 'View'), false,
      'revoked schema grant must disappear from a second read');
    expect(after.length, before);
  });
});
