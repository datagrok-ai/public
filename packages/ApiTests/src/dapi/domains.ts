import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Tests for grok.dapi.domains against the 'apitests' domain schema that this
// package declares in databases/apitests/schema.json (deployed on publish).
category('Dapi: domains', () => {
  const items = () => grok.dapi.domains.table('apitests.item');
  const sku = () => `SKU-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;

  test('schemas listing', async () => {
    const schemas = await grok.dapi.domains.schemas.list();
    const s = schemas.find((x) => x.name === 'apitests');
    expect(s != null, true, 'apitests schema not registered');
    expect(s!.pgSchema, 'apitests');
    expect(s!.tables.some((t) => t.name === 'item'), true, 'item table not registered');
    expect(s!.tables.find((t) => t.name === 'item')!.securityMode, 'row');
  });

  test('insert, query, get, update, audit, delete', async () => {
    const key = sku();
    const [ins] = await items().insert({sku: key, name: 'Widget', quantity: 5, note: 'first'});
    expect(ins.created, true);
    const id: string = ins.id;
    try {
      const rows = await items().query({filter: `sku = "${key}"`});
      expect(rows.length, 1);
      expect(rows[0].name, 'Widget');
      expect(rows[0].note, 'first', 'jsonb property schema key not returned');

      const row = await items().get(id);
      expect(row.sku, key);
      expect(row.version, 1);

      const upd = await items().update(id, {name: 'Widget 2'}, {version: 1});
      expect(upd.version, 2);

      const audit = await items().audit(id);
      expect(audit.map((a) => a.op).join(','), 'insert,update');
      expect(audit[1].before.name, 'Widget');
      expect(audit[1].after.name, 'Widget 2');
    } finally {
      await items().delete(id);
    }
    expect(await items().get(id), null, 'soft-deleted row is still visible');
  });

  test('optimistic concurrency conflict', async () => {
    const [ins] = await items().insert({sku: sku(), name: 'A'});
    try {
      await items().update(ins.id, {name: 'B'}, {version: 1});
      let conflict = '';
      try {
        await items().update(ins.id, {name: 'C'}, {version: 1});
      } catch (e: any) {
        conflict = e.message ?? `${e}`;
      }
      expect(conflict.includes('Version conflict'), true, `unexpected error: '${conflict}'`);
    } finally {
      await items().delete(ins.id);
    }
  });

  test('duplicate business key reported per row', async () => {
    const key = sku();
    const [first] = await items().insert({sku: key});
    try {
      const [dup] = await items().insert({sku: key});
      expect(dup.status, 'duplicate');
      expect(dup.existingId, first.id);
    } finally {
      await items().delete(first.id);
    }
  });

  test('promote makes a row an entity', async () => {
    const [ins] = await items().insert({sku: sku()});
    try {
      const res = await items().promote(ins.id);
      expect(res.promoted, true);
      const audit = await items().audit(ins.id);
      expect(audit[audit.length - 1].op, 'promote');
    } finally {
      await items().delete(ins.id);
    }
  });

  test('dapi2 generated client: domains removed at parity, init export intact', async () => {
    // WO-4b: every dapi2 function belonged to the domains namespace — with it
    // removed at typed-surface parity the generated client has no value members
    // left (the chats namespace never emitted functions), so `dapi2` is
    // type-only now; dapi2Init survives as a value and the OpenAPI yaml keeps
    // every /domains/ route.
    expect(typeof grok.dapi2Init, 'function');
    // NB: compare against true — utils expect() treats a passed undefined
    // `expected` as its default (true), so expect(x, undefined) never passes.
    expect((grok as any).dapi2?.domains === undefined, true, 'dapi2.domains must be gone');
  });

  test('table name validation', async () => {
    let error = '';
    try {
      grok.dapi.domains.table('noseparator');
    } catch (e: any) {
      error = e.message;
    }
    expect(error.includes('<schema>.<table>'), true);
  });
}, {owner: 'askalkin@datagrok.ai'});

// Registry reflection (ui-js-api WO-2): grok.dapi.domains.registry — the runtime
// Property metadata, table info with FK-inverted child tables, and batched
// display-name resolution. Assertions run against this package's own 'apitests'
// schema; grit.issue assertions (the dogfood schema with choices/min/nameColumn)
// skip cleanly where Grit is not deployed.
category('Dapi: domain registry', () => {
  const registry = () => grok.dapi.domains.registry;
  const sku = () => `SKU-REG-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;

  async function gritDeployed(): Promise<boolean> {
    const schemas = await grok.dapi.domains.schemas.list();
    return schemas.some((s) => s.name === 'grit');
  }

  test('rowProperties: schema.json constraints round-trip', async () => {
    const props = await registry().rowProperties('apitests.item');
    const by = (n: string) => props.find((p) => p.name === n);
    expect(by('sku') != null, true, `sku property missing: ${props.map((p) => p.name).join(',')}`);
    expect(by('sku')!.nullable, false, 'required column must not be nullable');
    expect(by('name')!.nullable, true, 'optional column must be nullable');
    expect(by('quantity')!.min, 0, 'min constraint from schema.json');
    expect(by('quantity')!.propertyType, 'int');
  });

  test('rowProperties: grit.issue matches schema.json (choices/min/nullable)', async () => {
    if (!await gritDeployed()) {
      console.log('skipped: Grit is not deployed');
      return;
    }
    const props = await registry().rowProperties('grit.issue');
    const by = (n: string) => props.find((p) => p.name === n)!;
    expect(JSON.stringify(by('status').choices), JSON.stringify(['open', 'in progress', 'resolved', 'closed']));
    expect(JSON.stringify(by('priority').choices), JSON.stringify(['low', 'medium', 'high', 'critical']));
    expect(by('number').min, 1, 'min from schema.json');
    expect(by('number').nullable, false, 'required column must not be nullable');
    expect(by('title').nullable, false);
    expect(by('description').nullable, true);
    expect(by('project_id').semType, 'grit.project', 'ref column must carry the target row semType');
    expect(by('project_id').friendlyName, 'Project',
      'ref columns must carry the label the platform renders, not the wire name');
  });

  test('rowProperties: unknown table rejects with a typed validation error', async () => {
    let e: any = null;
    try {
      await registry().rowProperties('apitests.nosuch');
    } catch (x) {
      e = x;
    }
    expect(e instanceof DG.DomainValidationError, true,
      `expected DomainValidationError, got ${e?.constructor?.name}: ${e?.message}`);
  });

  test('tableInfo: identity, security, and FK-inverted childTables', async () => {
    const info = await registry().tableInfo('apitests.item');
    expect(JSON.stringify(info.businessKey), JSON.stringify(['sku']));
    // NB: nameColumn is not asserted here — deployed registries may carry an
    // isName drift on apitests.item (same manifest version, no reapply); the
    // grit.issue test pins nameColumn against its committed schema.json.
    expect(info.securityMode, 'row');
    expect(info.audit, true);
    expect(info.singularName, 'item', 'effective singular derived from the table name');
    expect(info.pluralName, 'items');
    const child = info.childTables.find((c) => c.table === 'item_event');
    expect(child != null, true, `item_event missing from childTables: ${JSON.stringify(info.childTables)}`);
    expect(child!.schema, 'apitests');
    expect(child!.fkColumn, 'item_id');
    expect(child!.label, 'Item', 'friendly FK label (the _id suffix dropped)');
  });

  test('tableInfo: grit.issue childTables list comment', async () => {
    if (!await gritDeployed()) {
      console.log('skipped: Grit is not deployed');
      return;
    }
    const info = await registry().tableInfo('grit.issue');
    expect(info.nameColumn, 'title');
    expect(info.childTables.some((c) => c.table === 'comment' && c.fkColumn === 'issue_id'), true,
      `comment missing from childTables: ${JSON.stringify(info.childTables)}`);
  });

  test('resolveNames: display-identity chain, null for unresolvable ids', async () => {
    const items = grok.dapi.domains.table('apitests.item');
    const info = await registry().tableInfo('apitests.item');
    const named = sku();
    const bare = sku();
    const ghost = '00000000-0000-0000-0000-000000000000';
    let inserted: {id: string}[] = [];
    try {
      // One row with a name value, one without: the second proves the
      // business-key fallback regardless of whether the registry declares a
      // name column for apitests.item (deployed registries drift on isName).
      // A single insert inside the try so no partial pair can leak.
      inserted = await items.insert([{sku: named, name: 'Resolve probe'}, {sku: bare}]);
      const [insNamed, insBare] = inserted;
      const names = await registry().resolveNames('apitests.item', [insNamed.id, insBare.id, ghost]);
      expect(names[insNamed.id], info.nameColumn != null ? 'Resolve probe' : named,
        `display identity must follow the declared name column (${info.nameColumn})`);
      expect(names[insBare.id], bare, 'empty name — the business key is the display identity');
      expect(Object.keys(names).includes(ghost), true, 'every requested id must be a key');
      expect(names[ghost] == null, true, 'unresolvable ids must map to null');
    } finally {
      for (const r of inserted)
        await items.delete(r.id);
    }
  });
}, {owner: 'askalkin@datagrok.ai'});
