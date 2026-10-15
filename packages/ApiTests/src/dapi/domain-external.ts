import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {thrown, withRestrictedUser} from './domain-lifecycle';

// The typed client over a READ-ONLY external binding: ApiSamples' `northwind`
// (databases/northwind/schema.json over ApiSamples:PostgresNorthwind). Rows are
// addressed by their business key; every write and every platform-only read
// feature answers `unsupported` by name. Northwind is never written. The category
// skips cleanly while the schema is not registered (a development convenience
// only: the release evidence runs it positively).
category('Dapi: domain external', () => {
  const orders = () => grok.dapi.domains.table('northwind.order');
  const lines = () => grok.dapi.domains.table('northwind.order_detail');
  let skip: string | null = null;
  let extwlive = false;

  before(async () => {
    const schemas = await grok.dapi.domains.schemas.list();
    if (!schemas.some((s) => s.name === 'northwind'))
      skip = 'the northwind schema is not registered';
    extwlive = schemas.some((s) => s.name === 'extwlive');
  });

  const skipped = (): boolean => {
    if (skip != null)
      console.log(`skipped: ${skip}`);
    return skip != null;
  };

  const refused = async (op: string, action: () => Promise<any>) => {
    const err = await thrown(action);
    expect(err instanceof DG.DomainUnsupportedError, true,
      `${op}: expected DomainUnsupportedError, got ${err?.constructor?.name}: ${err?.message}`);
    expect(err.op, op);
  };

  test('access: support declares what a read-only binding cannot do; keys stay readonly', async () => {
    if (skipped())
      return;
    const access = await orders().access();
    const support = access.support;
    expect(Object.keys(support).sort().join(','),
      'ancestors,audit,batch,captions,concurrency,deleted,filters,probe,restore,systemColumns,transaction,' +
      'updateWhere,version,watch,writes', `unexpected support keys: ${JSON.stringify(support)}`);
    for (const k of ['writes', 'transaction', 'probe', 'deleted', 'restore', 'audit', 'watch', 'captions',
      'updateWhere', 'version', 'ancestors'])
      expect((support as any)[k], false, `${k} must be false on a read-only binding: ${JSON.stringify(support)}`);
    expect(support.concurrency, 'none');
    expect(support.filters, 'basic');
    for (const option of ['upsert', 'partial', 'validate', 'skipDuplicates'])
      expect((support.batch as any)[option], false, `batch.${option} must be false: ${JSON.stringify(support.batch)}`);
    expect(support.systemColumns.join(','), 'id', 'an external table carries only id');
    // Without writes nothing is 'editable', so the key stays 'readonly' — the
    // 'immutable' rewrite applies to a writable binding only.
    expect(access.fields['order_id'], 'readonly', JSON.stringify(access.fields));
    expect(access.can.insert || access.can.edit || access.can.delete, false, JSON.stringify(access.can));
    expect((await grok.dapi.domains.registry.tableInfo('northwind.order')).rowAddress, 'id');
    expect((await grok.dapi.domains.registry.tableInfo('apitests.item')).rowAddress, 'businessKey');
  });

  test('identity: id is the key as a string, composite ids carry the order id', async () => {
    if (skipped())
      return;
    const rows = await orders().query({limit: 3, sort: 'order_id'});
    expect(rows.length, 3);
    expect(rows[0].id, '10248');
    expect(rows[0].order_id, 10248, 'the key column keeps its declared int type');
    expect((await orders().get('10248')).ship_country, 'France');
    expect(await orders().count(), 830);
    expect(await lines().count(), 2155);
    const line = await lines().get('10248,11');
    expect(line.id, '10248,11');
    expect(line.order_id, '10248', 'a ref value on the wire is the target id');
    expect(line.product_id, 11);
    const ofOrder = await lines().query({filter: 'order_id = "10248"'});
    expect(ofOrder.length, 3, JSON.stringify(ofOrder));
    for (const l of ofOrder)
      expect(l.id.startsWith('10248,'), true, l.id);
  });

  test('projection: a columns list returns id and those columns only', async () => {
    if (skipped())
      return;
    const [row] = await orders().query({limit: 1, columns: ['ship_country']});
    expect(Object.keys(row).sort().join(','), 'id,ship_country', JSON.stringify(row));
  });

  test('refusals: every platform-only feature answers unsupported by name', async () => {
    if (skipped())
      return;
    await refused('captions', () => lines().query({captions: ['order_id']}));
    await refused('deleted', () => orders().query({deleted: 'only'}));
    await refused('filter', () => lines().query({filter: 'id in ("10248,11", "10248,42")'}));
    await refused('filter', () => orders().query({filter: 'ship_country = "France" or freight > 10 and customer_id = "VINET"'}));
    await refused('literal', () => orders().query({filter: 'ship_country like "%a_b%"'}));
    await refused('ancestors', () => orders().pathTo('10248'));
    await refused('insert', () => orders().insert({order_id: 1}));
    await refused('updateWhere', () => orders().updateWhere('order_id = 1', {freight: 1}));
    await refused('restore', () => orders().restore('10248'));
    await refused('audit', () => orders().audit('10248'));
    await refused('watch', () => orders().watch('10248'));
  });

  test('facets and aggregate: bounded, exact', async () => {
    if (skipped())
      return;
    const res = await orders().facets({facets: [{id: 'c', kind: 'categories', column: 'ship_country'}]});
    const categories = (res.facets['c'] as _DG.DomainFacetCategoriesResult).categories;
    expect(categories.length > 0, true, JSON.stringify(res));
    const france = categories.find((c) => c.value === 'France');
    expect(france?.filtered, 77, `France: ${JSON.stringify(france)}`);
    expect(france?.total, null, 'an external facet has no unfiltered total');
    const [total] = await orders().aggregate({measures: [{fn: 'count', as: 'n'}]});
    expect(total['n'], 830, JSON.stringify(total));
  });

  // Runs in the browser, i.e. through nginx, which turns a `%3A` in the path into a raw colon
  // that the row route does not match: the client sends the id as one more-encoded segment.
  test('a datetime key reaches its row through the proxy', async () => {
    if (skipped())
      return;
    if (!extwlive) {
      console.log('skipped: the extwlive schema is not registered');
      return;
    }
    const id = '2026-09-18T07%3A08%3A09.123Z';
    const row = await grok.dapi.domains.table('extwlive.dtkeyed').get(id);
    expect(row?.id, id, `the canonical id is what the row answers: ${JSON.stringify(row)}`);
  });

  test('a restricted user: View on the table reads rows, no View reads nothing', async () => {
    if (skipped())
      return;
    await withRestrictedUser('ext', async (probe) => {
      await orders().grant(probe.group, 'View');
      try {
        expect((await probe.asUser(() => orders().query({limit: 1}))).length, 1,
          'a table View grant must read rows without any right on the connection');
        expect((await probe.asUser(() => lines().query({limit: 1}))).length, 0,
          'no View on order_detail answers no rows, not 403');
      } finally {
        await orders().revoke(probe.group, 'View');
        grok.dapi.domains.invalidateUiCaches();
      }
    });
  });
});
