import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {items, uniqueKey} from './domains-fixture';

// Tests for grok.dapi.domains against the 'apitests' domain schema that this
// package declares in databases/apitests/schema.json (deployed on publish).
category('Dapi: domains', () => {
  const sku = () => uniqueKey('SKU');

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

  test('dapi2 generated client smoke', async () => {
    grok.dapi2Init(grok.dapi.root, grok.dapi.token);
    const key = sku();
    const [ins] = await items().insert({sku: key, name: 'Gen'});
    try {
      const rows = await grok.dapi2.domains.queryRows('apitests', 'item', {filter: `sku = "${key}"`});
      expect(rows.length, 1);
      expect(rows[0].name, 'Gen');
    } finally {
      await items().delete(ins.id);
    }
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
