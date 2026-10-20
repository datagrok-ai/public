import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, expectObject, test} from '@datagrok-libraries/test/src/test';

// Phase-2 surface of grok.dapi.domains (batch/upsert, transactions, aggregate,
// queryDf) against the 'apitests' schema from databases/apitests/schema.json.
category('Dapi: domains batch', () => {
  const items = () => grok.dapi.domains.table('apitests.item');
  const prefix = () => `BT-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;

  async function purge(p: string): Promise<void> {
    // One filtered bulk delete instead of the old query + per-row N+1 loop.
    for (let guard = 0; guard < 100; guard++)
      if (!(await items().deleteWhere(`sku starts "${p}"`)).hasMore)
        return;
  }

  test('batch CSV string', async () => {
    const p = prefix();
    try {
      const res = await items().batch(`sku,name,quantity\n${p}-1,Widget,1\n${p}-2,Widget,2`);
      expect(res.inserted, 2);
      expect(res.errorCount, 0);
      const rows = await items().query({filter: `sku starts "${p}"`});
      expect(rows.length, 2);
    } finally {
      await purge(p);
    }
  });

  test('batch DataFrame and d42 bytes', async () => {
    const p = prefix();
    try {
      const df = DG.DataFrame.fromCsv(`sku,name,quantity\n${p}-1,W,1\n${p}-2,W,2`);
      const res = await items().batch(df);
      expect(res.inserted, 2);
      const df2 = DG.DataFrame.fromCsv(`sku,name,quantity\n${p}-3,W,3`);
      const res2 = await items().batch(new Uint8Array(df2.toByteArray()));
      expect(res2.inserted, 1);
      const rows = await items().query({filter: `sku starts "${p}"`});
      expect(rows.length, 3);
    } finally {
      await purge(p);
    }
  });

  test('batch parquet', async () => {
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
      const df = DG.DataFrame.fromCsv(`sku,name,quantity\n${p}-1,P,1\n${p}-2,P,2`);
      const bytes: Uint8Array = await grok.functions.call('Arrow:toParquet', {table: df});
      const res = await items().batch(bytes, {format: 'parquet'});
      expect(res.inserted, 2);
      const rows = await items().query({filter: `sku starts "${p}"`, sort: 'sku'});
      expect(rows.length, 2);
      expect(rows[1].quantity, 2);
    } finally {
      await purge(p);
    }
  });

  test('batch upsert counts', async () => {
    const p = prefix();
    const mk = (sku: string, quantity: number) => ({sku, name: 'U', quantity});
    try {
      const first = await items().batch([mk(`${p}-1`, 1), mk(`${p}-2`, 2)]);
      expect(first.inserted, 2);
      const second = await items().batch([mk(`${p}-1`, 10), mk(`${p}-2`, 20), mk(`${p}-3`, 3)], {mode: 'upsert'});
      expect(second.inserted, 1);
      expect(second.updated, 2);
      expect(second.errorCount, 0);
      const rows = await items().query({filter: `sku = "${p}-1"`});
      expect(rows[0].quantity, 10);
    } finally {
      await purge(p);
    }
  });

  test('batch partial-success report', async () => {
    const p = prefix();
    try {
      const res = await items().batch(
        [{sku: `${p}-ok`, quantity: 1}, {sku: `${p}-bad`, quantity: -5}], {allOrNothing: false});
      expect(res.inserted, 1);
      expect(res.errorCount, 1);
      const bad = res.rows.find((r) => r.status === 'error');
      expect(bad != null, true, 'no error row in the report');
      expect(bad!.index, 1);
      expect(bad!.errors![0].column, 'quantity');
    } finally {
      await purge(p);
    }
  });

  test('batch allOrNothing abort', async () => {
    const p = prefix();
    const res = await items().batch([{sku: `${p}-ok`, quantity: 1}, {sku: `${p}-bad`, quantity: -5}]);
    expect(res.error, 'validation');
    expect(res.errorCount, 1);
    const rows = await items().query({filter: `sku starts "${p}"`});
    expect(rows.length, 0, 'allOrNothing batch must not apply any rows');
  });

  test('validateOnly predicts exactly what the commit then does', async () => {
    const p = prefix();
    try {
      // A payload with one of everything: a row that exists (so upsert updates it), two new
      // rows, and one the fixture refuses.
      await items().insert({sku: `${p}-1`, name: 'Seed', quantity: 1});
      const payload = [
        {sku: `${p}-1`, name: 'Existing', quantity: 9},
        {sku: `${p}-2`, name: 'New A', quantity: 2},
        {sku: `${p}-3`, name: 'New B', quantity: 3},
        {sku: `${p}-4`, name: 'Refused', quantity: -5},
      ];
      const before = (await items().query({filter: `sku starts "${p}"`})).length;
      const preview = await items().batch(payload, {mode: 'upsert', allOrNothing: false, validateOnly: true});
      expect(preview.validateOnly, true, `the preview must answer its own shape: ${JSON.stringify(preview)}`);
      expect(preview.rowCount, 4, JSON.stringify(preview));
      expect(preview.willInsert, 2, `two new business keys: ${JSON.stringify(preview)}`);
      expect(preview.willUpdate, 1, `one matched business key: ${JSON.stringify(preview)}`);
      expect(preview.errorCount, 1, `one refused row: ${JSON.stringify(preview)}`);
      // A predicted insert has no id — the id of a rolled-back insert does not exist.
      for (const row of preview.rows)
        expect('id' in row, false, `a preview row must carry no id: ${JSON.stringify(row)}`);
      const bad = preview.rows.find((r) => r.predicted === 'error');
      expect(bad != null, true, `no predicted error row: ${JSON.stringify(preview.rows)}`);
      expect(bad!.index, 3, 'the preview must name the payload row it refuses');
      expect(bad!.errors![0].column, 'quantity', JSON.stringify(bad));
      const updated = preview.rows.find((r) => r.predicted === 'update');
      expect(updated!.existingId != null, true,
        `a predicted update must name the row it would merge into: ${JSON.stringify(updated)}`);
      // Nothing was written.
      expect((await items().query({filter: `sku starts "${p}"`})).length, before,
        'a validate-only batch must write no rows');
      expect((await items().query({filter: `sku = "${p}-1"`}))[0].quantity, 1,
        'a validate-only batch must not update the row it predicted an update for');

      // Now commit the same payload: every verdict must match, one to one.
      const report = await items().batch(payload, {mode: 'upsert', allOrNothing: false});
      expect(report.inserted, preview.willInsert, 'willInsert disagrees with the commit');
      expect(report.updated, preview.willUpdate, 'willUpdate disagrees with the commit');
      expect(report.skipped, preview.willSkip, 'willSkip disagrees with the commit');
      expect(report.errorCount, preview.errorCount, 'errorCount disagrees with the commit');
      const predicted = new Map(preview.rows.map((r) => [r.index, r.predicted]));
      const committed: {[status: string]: string} =
        {inserted: 'insert', updated: 'update', duplicate: 'skip', error: 'error'};
      for (const row of report.rows)
        expect(predicted.get(row.index), committed[row.status],
          `row ${row.index}: predicted ${predicted.get(row.index)}, committed ${row.status}`);
    } finally {
      await purge(p);
    }
  });

  test('batch throwOnError rejects', async () => {
    const p = prefix();
    const payload = [{sku: `${p}-ok`, quantity: 1}, {sku: `${p}-bad`, quantity: -5}];
    try {
      const report = await items().batch(payload);
      let err: any = null;
      try {
        await items().batch(payload, {throwOnError: true});
      } catch (e) {
        err = e;
      }
      expect(err instanceof DG.DomainValidationError, true,
        `expected DomainValidationError, got ${err?.constructor?.name}: ${err?.message}`);
      expectObject(err.body, report);

      const partial = await items().batch(payload, {throwOnError: true, allOrNothing: false});
      expect(partial.errorCount, 1, 'partial mode is not an abort and must resolve');
    } finally {
      await purge(p);
    }
  });

  test('validateOnly still refuses a malformed request', async () => {
    const p = prefix();
    let error: any = null;
    try {
      await items().batch([{sku: `${p}-1`, nosuchcolumn: 1}], {validateOnly: true});
    } catch (e: any) {
      error = e;
    }
    // A request-level refusal precedes the transaction, so it is a refusal either way.
    expect(error != null, true, 'an unknown column must reject a preview as it rejects a commit');
    expect(error instanceof DG.DomainError, true,
      `expected a typed domain error, got ${error?.constructor?.name}: ${error?.message}`);
  });

  test('transaction $ref across tables', async () => {
    const key = `${prefix()}-tx`;
    const res = await grok.dapi.domains.transaction('apitests', [
      {op: 'insert', table: 'item', ref: 'i', values: {sku: key, name: 'TxWidget', quantity: 1}},
      {op: 'insert', table: 'item_event', values: {item_id: '$i', kind: 'created', amount: 5}},
      {op: 'update', table: 'item', id: '$i', values: {quantity: 2}, expectedVersion: 1},
    ]);
    try {
      expect(res.length, 3);
      // the mapped-tuple transaction types each result by its op — no casts
      expect(res[0].created, true);
      expect(res[1].created, true);
      expect(res[2].version, 2);
      const rows = await items().query({filter: `sku = "${key}"`, expand: ['details:item_event']});
      expect(rows[0].item_event.length, 1);
      expect(rows[0].item_event[0].kind, 'created');
    } finally {
      await items().delete(res[0].id); // cascades item_event
    }
  });

  test('transaction rollback carries the typed opIndex', async () => {
    const key = `${prefix()}-rb`;
    let error: any = null;
    try {
      await grok.dapi.domains.transaction('apitests', [
        {op: 'insert', table: 'item', ref: 'i', values: {sku: key}},
        {op: 'insert', table: 'item_event', values: {item_id: '$i'}}, // missing required 'kind'
      ]);
    } catch (e: any) {
      error = e;
    }
    expect(error instanceof DG.DomainValidationError, true,
      `expected DomainValidationError, got ${error?.constructor?.name}: ${error?.message}`);
    expect(error.opIndex, 1, 'the error must point at the failing op');
    const rows = await items().query({filter: `sku = "${key}"`});
    expect(rows.length, 0, 'first op was not rolled back');
  });

  test('aggregate', async () => {
    const p = prefix();
    try {
      await items().batch([
        {sku: `${p}-1`, name: 'AgA', quantity: 1}, {sku: `${p}-2`, name: 'AgA', quantity: 2},
        {sku: `${p}-3`, name: 'AgB', quantity: 3}, {sku: `${p}-4`, name: 'AgB', quantity: 4}]);
      const res = await items().aggregate<string, string>({
        groupBy: ['name'],
        measures: [{fn: 'count'}, {fn: 'sum', column: 'quantity'}, {fn: 'max', column: 'quantity', as: 'top'}],
        filter: `sku starts "${p}"`,
        sort: 'name',
      });
      expect(res.length, 2);
      expect(res[0].name, 'AgA');
      expect(res[0].count, 2);
      expect(res[0].sum_quantity, 3);
      expect(res[1].top, 4);
    } finally {
      await purge(p);
    }
  });

  test('queryDf values and tags', async () => {
    const p = prefix();
    try {
      await items().insert([
        {sku: `${p}-1`, name: 'DfA', quantity: 1, note: 'n1'},
        {sku: `${p}-2`, name: 'DfB', quantity: 2, note: 'n2'}]);
      const df = await items().queryDf({filter: `sku starts "${p}"`, sort: 'sku'});
      expect(df.rowCount, 2);
      expect(df.col('sku')!.get(0), `${p}-1`);
      expect(df.col('quantity')!.type, DG.COLUMN_TYPE.INT);
      expect(df.col('quantity')!.get(1), 2);
      expect(df.col('sku')!.getTag('dbPropertyName'), 'sku');
      expect(df.col('note')!.getTag('dbPropertySchema'), 'item_meta');
      expect(!df.col('id')!.getTag('dbPropertyName'), true, 'system columns must stay untagged');
    } finally {
      await purge(p);
    }
  });
}, {owner: 'askalkin@datagrok.ai'});
