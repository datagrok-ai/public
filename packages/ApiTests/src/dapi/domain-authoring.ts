import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {thrown, withRestrictedUser} from './domain-lifecycle';

// Authoring an external binding from the product: draft a manifest over a database the
// caller may query, dry-run it, create the schema, validate the binding, delete it. Rides
// the stand's `NorthwindBinding:PostgresNorthwind` connection (the local demo container)
// and registers a throwaway schema of its own — Northwind itself is never touched, and the
// `northwind` schema of the ApiSamples fixture is neither read nor written. Skips cleanly
// without the connection or the CreateDomainSchema privilege.
category('Dapi: domain authoring', () => {
  const connection = 'NorthwindBinding:PostgresNorthwind';
  const request = {connection, schema: 'public', tables: ['orders', 'order_details']};
  let draft: _DG.DomainDraft | null = null;
  let skip: string | null = null;

  before(async () => {
    try {
      draft = await grok.dapi.domains.draft(request);
    } catch (e: any) {
      if (e instanceof DG.DomainError && e.code === 'unknown-connection')
        skip = `${e.code}: ${e.message}`;
      else
        throw e;
    }
  });

  const skipped = (): boolean => {
    if (skip != null)
      console.log(`skipped: ${skip}`);
    return skip != null;
  };

  test('draft: remote names, the primary key as the business key, one foreign key as a ref', async () => {
    if (skipped())
      return;
    const d = draft!;
    expect(Object.keys(d.manifest.tables).sort().join(','), 'order_details,orders', JSON.stringify(d.manifest));
    expect(d.manifest.storage.kind, 'external');
    expect(d.manifest.storage.connection, connection);
    expect(d.manifest.storage.schema, 'public');
    expect(d.manifest.tables.orders.businessKey.join(','), 'orderid');
    expect(d.manifest.tables.order_details.businessKey.join(','), 'orderid,productid');
    expect(d.manifest.tables.order_details.columns.orderid.type, 'ref');
    expect(d.manifest.tables.order_details.columns.orderid.ref, 'orders');
    expect(d.manifest.tables.orders.columns.orderid.required, true, 'a key column is required');

    // The inventory covers the whole remote schema, drafted or not — what else could be added.
    const byRemote = (remote: string) => d.inventory.tables.find((t) => t.remote === remote);
    expect(`${byRemote('orders')?.logical}:${byRemote('orders')?.key?.join('+')}`, 'orders:orderid',
      JSON.stringify(d.inventory.tables));
    expect(`${byRemote('order_details')?.logical}:${byRemote('order_details')?.key?.join('+')}`,
      'order_details:orderid+productid');
    expect(byRemote('customers')?.bindable, true, JSON.stringify(byRemote('customers')));
    expect('customers' in d.manifest.tables, false, 'an undrafted table stays out of the manifest');
    const refs = d.inventory.relations.filter((r) => r.status === 'ref');
    expect(refs.length, 1, JSON.stringify(d.inventory.relations));
    expect(`${refs[0].table}.${refs[0].column}->${refs[0].targetTable}.${refs[0].targetColumn}`,
      'order_details.orderid->orders.orderid');
    const out = d.inventory.relations.find((r) => r.table === 'orders' && r.column === 'customerid');
    expect(out?.status, 'plain', JSON.stringify(out));
    expect(out?.code, 'external-ref-out', 'the target is not in the draft');
    expect(d.diagnostics.length, 0, JSON.stringify(d.diagnostics));
  });

  test('dry run → create → manifest → validate → delete', async () => {
    if (skipped())
      return;
    const name = `auth_${Math.random().toString(36).slice(2, 8)}`;
    const manifest = draft!.manifest;

    const plan = await grok.dapi.domains.createSchema(name, {manifest, dryRun: true});
    expect(plan.status, 'ok', JSON.stringify(plan));
    expect(plan.issues.length, 0, JSON.stringify(plan));
    expect((await grok.dapi.domains.schemas.list()).some((s) => s.name === name), false,
      'a dry run registers nothing');

    const created = await grok.dapi.domains.createSchema(name, {friendlyName: 'Authoring probe', manifest});
    const handle = grok.dapi.domains.schema(name);
    try {
      expect(created.name, name, JSON.stringify(created));
      expect(created.pgSchema, `ext_${name}`);
      expect(created.binding?.status, 'ok', JSON.stringify(created));

      const stored = await handle.manifest();
      expect(stored.name, name);
      expect(stored.storage.kind, 'external');
      expect(Object.keys(stored.tables).sort().join(','), 'order_details,orders', JSON.stringify(stored));
      expect(stored.tables.order_details.columns.orderid.ref, 'orders');

      const validation = await handle.validate();
      expect(validation.schema, name);
      expect(validation.status, 'ok', JSON.stringify(validation));
      expect(validation.issues.length, 0);

      // The binding is live: the registered table reads the remote rows.
      expect(await grok.dapi.domains.table(`${name}.orders`).count(), 830);

      const again = await thrown(() => grok.dapi.domains.createSchema(name, {manifest, dryRun: true}));
      expect(again instanceof DG.DomainError, true, `${again?.constructor?.name}: ${again?.message}`);
      expect(again.code, 'schema-name-taken');
    } finally {
      await handle.delete();
    }
    expect((await grok.dapi.domains.schemas.list()).some((s) => s.name === name), false,
      'deleted schema must not be listed');
  });

  test('refusals by name, addressed by manifest path', async () => {
    if (skipped())
      return;
    const bad = await thrown(() => grok.dapi.domains.draft({connection, schema: 'pub.lic'}));
    expect(bad instanceof DG.DomainManifestValidationError, true, `${bad?.constructor?.name}: ${bad?.message}`);
    expect(bad.errors.map((e: _DG.DomainManifestIssue) => `${e.path}:${e.code}`).join(','),
      'schema:external-identifier', JSON.stringify(bad.body));

    const platform = await thrown(() => grok.dapi.domains.createSchema('auth_platform', {
      manifest: {tables: {things: {columns: {name: {type: 'string'}}}}}, dryRun: true}));
    expect(platform instanceof DG.DomainError, true, `${platform?.constructor?.name}: ${platform?.message}`);
    expect(platform.code, 'invalid-storage', JSON.stringify(platform?.body));
  });

  test('a user without CreateDomainSchema cannot draft', async () => {
    if (skipped())
      return;
    await withRestrictedUser('auth', async (probe) => {
      const err = await thrown(() => probe.asUser(() => grok.dapi.domains.draft(request)));
      expect(err instanceof DG.DomainForbiddenError, true, `${err?.constructor?.name}: ${err?.message}`);
      expect(err.status, 403);
    });
  });
});
