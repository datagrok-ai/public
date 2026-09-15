import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {thrown} from './domain-lifecycle';

// DomainQuerySpec.search: a case-insensitive substring over the table's
// searchable columns — apitests.item declares `name` searchable (schema.json);
// apitests.item_event declares nothing and has no name column.
category('Dapi: domain search', () => {
  const items = () => grok.dapi.domains.table('apitests.item');
  const events = () => grok.dapi.domains.table('apitests.item_event');
  const stamp = () => `${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
  const skuLike = (prefix: string): any => ({property: 'sku', operator: 'like', value: `${prefix}%`});

  async function seed(prefix: string): Promise<void> {
    await items().insert([
      {sku: `${prefix}-0`, name: `Alpha ${prefix}`},
      {sku: `${prefix}-1`, name: `Beta ${prefix}`},
      {sku: `${prefix}-2`, name: `Gamma ${prefix}`},
    ]);
  }

  async function cleanup(prefix: string): Promise<void> {
    try {
      await items().deleteWhere(skuLike(prefix));
    } catch (e) {
      console.error(`search fixture ${prefix} not cleaned up: ${e}`);
    }
  }

  test('query search: a case-insensitive substring over the searchable columns, ANDed with filter', async () => {
    const prefix = `se-q-${stamp()}`;
    await seed(prefix);
    try {
      const rows = await items().query({filter: skuLike(prefix), search: 'ALPH'});
      expect(rows.length, 1, `search did not narrow: ${JSON.stringify(rows.map((r) => r.name))}`);
      expect(rows[0].name, `Alpha ${prefix}`, 'search matched the wrong row');
      // sku is not searchable: a substring only the sku holds matches nothing.
      expect((await items().query({filter: skuLike(prefix), search: `${prefix}-1`})).length, 0,
        'search matched a column that is not searchable');
      const df = await items().queryDf({filter: skuLike(prefix), search: 'beta'});
      expect(df.rowCount, 1, 'queryDf ignored search');
    } finally {
      await cleanup(prefix);
    }
  });

  test('builder .search() and count agree', async () => {
    const prefix = `se-b-${stamp()}`;
    await seed(prefix);
    try {
      const rows = await items().query().where('sku', 'like', `${prefix}%`).search('gamma');
      expect(rows.length, 1, 'the builder search did not narrow');
      expect(await items().query().where('sku', 'like', `${prefix}%`).search('gamma').count(), 1,
        'the builder count disagrees with its rows under search');
      expect(await items().count(skuLike(prefix), {search: 'a'}), 3, 'count() ignored search');
      expect(await items().count(skuLike(prefix), {search: 'beta'}), 1, 'count() did not narrow by search');
    } finally {
      await cleanup(prefix);
    }
  });

  test('tableInfo.searchableColumns lists the declared column', async () => {
    const info = await grok.dapi.domains.registry.tableInfo('apitests.item');
    expect(JSON.stringify(info.searchableColumns), JSON.stringify(['name']),
      `expected ['name']: ${JSON.stringify(info.searchableColumns)}`);
    const tag = await grok.dapi.domains.registry.tableInfo('apitests.tag');
    expect(JSON.stringify(tag.searchableColumns), JSON.stringify(['name']),
      `a table with a name column and no searchable one defaults to it: ${JSON.stringify(tag.searchableColumns)}`);
    const event = await grok.dapi.domains.registry.tableInfo('apitests.item_event');
    expect(event.searchableColumns.length, 0, `a table with neither lists none: ${JSON.stringify(event.searchableColumns)}`);
  });

  test('toQuery() refuses a search instead of dropping it', async () => {
    const builder = () => items().query().where('sku', 'like', 'se-%');
    let refusal: any = null;
    try {
      builder().search('alpha').toQuery();
    } catch (e) {
      refusal = e;
    }
    expect(refusal != null, true, 'toQuery() silently dropped the search');
    expect(`${refusal.message}`.includes('search'), true,
      `the refusal does not name the search: ${refusal.message}`);
    expect(builder().toQuery().filters!.length, 1, 'toQuery() without a search stopped working');
  });

  test('a table without a searchable column rejects with DomainFilterError', async () => {
    const e = await thrown(() => events().query({search: 'x'}));
    expect(e instanceof DG.DomainFilterError, true,
      `expected DomainFilterError, got ${e?.constructor?.name}: ${e?.message}`);
  });
}, {owner: 'askalkin@datagrok.ai'});
