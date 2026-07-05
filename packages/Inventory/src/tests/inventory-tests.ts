import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {inventoryDb} from '../generated/db';
import {adjustStock} from '../inventory-app';

// CI tests for the 'inventory' domain schema (databases/inventory/schema.json),
// exercised through the typed clients that `grok api` generates (src/generated/db.ts).
category('Inventory', () => {
  const items = () => inventoryDb.items;
  const movements = () => inventoryDb.stockMovements;
  const prefix = () => `INV-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;

  async function purge(p: string): Promise<void> {
    const rows = await items().query({filter: `sku starts "${p}"`, limit: 1000});
    for (const r of rows)
      await items().delete(r.id); // stock_movements cascade with their item
  }

  test('schema registration', async () => {
    const schemas = await grok.dapi.domains.schemas.list();
    const s = schemas.find((x) => x.name === 'inventory');
    expect(s != null, true, 'inventory schema not registered');
    const itemsTable = s!.tables.find((t) => t.name === 'items');
    const movementsTable = s!.tables.find((t) => t.name === 'stock_movements');
    expect(itemsTable!.securityMode, 'table');
    expect(itemsTable!.businessKey.join(','), 'sku');
    expect(movementsTable!.securityMode, 'master');
    expect(movementsTable!.audit, false, 'stock_movements must be the audit-off example');
  });

  // Phase-2 milestone: "Inventory stock sync upserts by business key".
  test('stock sync upserts by business key (CSV)', async () => {
    const p = prefix();
    try {
      const first = await items().batch(
        `sku,name,quantity,location\n${p}-1,Acetone,10,A\n${p}-2,Ethanol,20,A\n${p}-3,Methanol,30,B`,
        {mode: 'upsert'});
      expect(first.inserted, 3);
      expect(first.updated, 0);
      expect(first.errorCount, 0);

      const second = await items().batch(
        `sku,name,quantity,location\n${p}-1,Acetone,15,A\n${p}-3,Methanol,5,B\n${p}-4,Toluene,40,B`,
        {mode: 'upsert'});
      expect(second.inserted, 1, 'only the new SKU must be inserted');
      expect(second.updated, 2, 'existing SKUs must be updated in place');
      expect(second.errorCount, 0);

      const rows = await items().query({filter: `sku starts "${p}"`, sort: 'sku'});
      expect(rows.length, 4);
      expect(rows[0].quantity, 15, 'upsert must merge by sku');
      expect(rows[2].quantity, 5);
    } finally {
      await purge(p);
    }
  });

  test('parquet import via Arrow', async () => {
    const p = prefix();
    const arrow = DG.Func.find({package: 'Arrow', name: 'fromParquet'});
    if (arrow.length === 0) {
      // Skip-with-reason: Arrow is not deployed — verify the clear-error fallback instead.
      let error = '';
      try {
        await items().batch(new Uint8Array([1, 2, 3]), {format: 'parquet'});
      } catch (e: any) {
        error = e.message ?? `${e}`;
      }
      expect(error.includes('Arrow'), true, `expected an error naming the Arrow package, got '${error}'`);
      return;
    }
    try {
      const df = DG.DataFrame.fromCsv(
        `sku,name,quantity,location\n${p}-1,Benzene,7,C\n${p}-2,Xylene,9,C`);
      const bytes: Uint8Array = await grok.functions.call('Arrow:toParquet', {table: df});
      const report = await items().batch(bytes, {format: 'parquet', mode: 'upsert'});
      expect(report.inserted, 2);
      expect(report.errorCount, 0);
      const rows = await items().query({filter: `sku starts "${p}"`, sort: 'sku'});
      expect(rows.length, 2);
      expect(rows[1].quantity, 9);
    } finally {
      await purge(p);
    }
  });

  test('optimistic concurrency: 409 surfaces, retry succeeds', async () => {
    const p = prefix();
    const [ins] = await items().insert({sku: `${p}-1`, name: 'Contended', quantity: 10, location: 'A'});
    try {
      await items().update(ins.id, {name: 'Contended 2'}, {version: 1});
      let conflict = '';
      try {
        await items().update(ins.id, {name: 'Contended 3'}, {version: 1});
      } catch (e: any) {
        conflict = e.message ?? `${e}`;
      }
      expect(conflict.includes('Version conflict'), true, `unexpected error: '${conflict}'`);

      // Two concurrent adjustments: the loser retries with the fresh version and wins.
      await Promise.all([
        adjustStock(ins.id, 5, 'received'),
        adjustStock(ins.id, -3, 'shipped'),
      ]);
      const row = await items().get(ins.id);
      expect(row.quantity, 12, 'both adjustments must land exactly once');

      const moves = await movements().query({filter: `item_id = "${ins.id}"`});
      expect(moves.length, 2);
      expect((await movements().audit(moves[0].id)).length, 0,
        'audit-off stock_movements must write no audit rows');
    } finally {
      await purge(p);
    }
  });

  test('stock on hand by location aggregation', async () => {
    const p = prefix();
    try {
      await items().batch([
        {sku: `${p}-1`, name: 'A1', quantity: 10, location: 'Warehouse A'},
        {sku: `${p}-2`, name: 'A2', quantity: 15, location: 'Warehouse A'},
        {sku: `${p}-3`, name: 'B1', quantity: 7, location: 'Warehouse B'},
      ]);
      const res = await items().aggregate({
        groupBy: ['location'],
        measures: [{fn: 'sum', column: 'quantity', as: 'stock_on_hand'}, {fn: 'count', as: 'items'}],
        filter: `sku starts "${p}"`,
        sort: 'location',
      });
      expect(res.length, 2);
      expect(res[0].location, 'Warehouse A');
      expect(res[0].stock_on_hand, 25);
      expect(res[0].items, 2);
      expect(res[1].stock_on_hand, 7);
    } finally {
      await purge(p);
    }
  });

  // Per-group column visibility itself cannot be asserted from this session (the CI
  // user is the publisher/admin and sees every column; granting a property schema to
  // another group has no exposed API yet — see the package README): this test pins the
  // registration and the per-schema column attribution the visibility check gates on.
  test('per-schema column attribution', async () => {
    const schemas = await grok.dapi.stickyMeta.getSchemas();
    for (const name of ['chemistry', 'procurement', 'quality'])
      expect(schemas.some((s) => s.name === name), true, `property schema '${name}' not registered`);

    const p = prefix();
    try {
      await items().insert({sku: `${p}-1`, name: 'Tagged', quantity: 1, location: 'A',
        cas_number: '67-64-1', hazard_class: 'flammable', unit_cost: 9.5, inspection_status: 'passed'});
      const df = await items().queryDf({filter: `sku starts "${p}"`});
      expect(df.col('cas_number')!.getTag('dbPropertySchema'), 'chemistry');
      expect(df.col('unit_cost')!.getTag('dbPropertySchema'), 'procurement');
      expect(df.col('inspection_status')!.getTag('dbPropertySchema'), 'quality');
      const core = df.col('sku')!.getTag('dbPropertySchema');
      expect(core !== 'chemistry' && core !== 'procurement' && core !== 'quality', true,
        'relational columns must belong to the core schema, not a department schema');
    } finally {
      await purge(p);
    }
  });

  test('demo setup creates the department groups', async () => {
    const names = ['Chemists', 'Procurement'];
    const existedBefore = new Set<string>();
    for (const n of names)
      if (await grok.dapi.groups.filter(`shortName = "${n}"`).first() != null)
        existedBefore.add(n);
    try {
      const summary: string = await grok.functions.call('Inventory:SetupInventoryDemo');
      expect(summary.length > 0, true);
      for (const n of names)
        expect(await grok.dapi.groups.filter(`shortName = "${n}"`).first() != null, true,
          `group '${n}' was not created`);
    } finally {
      for (const n of names) {
        if (existedBefore.has(n))
          continue;
        const g = await grok.dapi.groups.filter(`shortName = "${n}"`).first();
        if (g != null)
          await grok.dapi.groups.delete(g);
      }
    }
  });
}, {owner: 'askalkin@datagrok.ai'});
