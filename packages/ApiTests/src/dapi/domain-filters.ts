import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';

// Entity-filters (GROK-20591): the batched facets endpoint and saved filter
// presets through the JS API — grok.dapi.domains.table(...).facets(spec) and
// .filters list/save/delete. The fixture is this package's own 'apitests'
// domain schema (databases/apitests/schema.json, deployed on publish); every
// test creates its own uniquely-marked rows and removes them in cleanup, so
// leftover rows from other runs never perturb the counts (facet assertions
// isolate via a filter/search on the marker).
category('Dapi: domain filters', () => {
  const items = () => grok.dapi.domains.table('apitests.item');
  const events = () => grok.dapi.domains.table('apitests.item_event');
  const marker = `FCT${Date.now()}${Math.floor(Math.random() * 1e4)}`;
  const sku = () => `${marker}-${Math.floor(Math.random() * 1e9)}`;
  const itemIds: string[] = [];
  const eventIds: string[] = [];

  before(async () => {
    // Three quantified items (two sharing a name) + one extra parent pair for
    // the FK-path test: two events on P1, one on P2.
    const ins = await items().insert([
      {sku: sku(), name: `${marker} red`, quantity: 11},
      {sku: sku(), name: `${marker} red`, quantity: 19},
      {sku: sku(), name: `${marker} blue`, quantity: 30},
      {sku: sku(), name: `${marker} P1`},
      {sku: sku(), name: `${marker} P2`},
    ]);
    for (const r of ins)
      itemIds.push(r.id);
    const evs = await events().insert([
      {item_id: itemIds[3], kind: 'e1', amount: 1},
      {item_id: itemIds[3], kind: 'e1', amount: 2},
      {item_id: itemIds[4], kind: 'e2', amount: 3},
    ]);
    for (const r of evs)
      eventIds.push(r.id);
  });

  after(async () => {
    for (const id of eventIds)
      await events().delete(id).catch(() => {});
    for (const id of itemIds)
      await items().delete(id).catch(() => {});
  });

  test('facets round trip: categories, histogram, minMax, count in one batch', async () => {
    const filter: _DG.DomainConditionTree = [{property: 'name', operator: 'like', value: `${marker}%`}];
    const res = await items().facets({filter, facets: [
      {id: 'names', kind: 'categories', column: 'name', search: marker},
      {id: 'hist', kind: 'histogram', column: 'quantity', bins: 2, min: 10, max: 30},
      {id: 'mm', kind: 'minMax', column: 'quantity'},
      {id: 'n', kind: 'count'},
    ]});

    const names = res.facets['names'] as _DG.DomainFacetCategoriesResult;
    const byValue: {[v: string]: any} = {};
    for (const c of names.categories)
      byValue[c.value] = c;
    expect(byValue[`${marker} red`]?.total, 2, `red total: ${JSON.stringify(names)}`);
    expect(byValue[`${marker} blue`]?.total, 1, `blue total: ${JSON.stringify(names)}`);

    // Histogram BUCKETS and count are filter-aware (the name filter is not
    // their own column), so only the marked rows count. minMax — like plan
    // and histogram BOUNDS — is computed under the row predicate only (the
    // stable-axis rule ignores the filter tree), so leftover rows from other
    // suites may widen it: assert containment, not equality.
    const hist = res.facets['hist'] as _DG.DomainFacetHistogramResult;
    expect(hist.buckets.reduce((a: number, b: number) => a + b, 0), 3,
      `histogram buckets: ${JSON.stringify(hist)}`);
    const mm = res.facets['mm'] as _DG.DomainFacetMinMaxResult;
    expect((mm.min as number) <= 11, true, `minMax must contain the marked rows: ${JSON.stringify(mm)}`);
    expect((mm.max as number) >= 30, true, `minMax must contain the marked rows: ${JSON.stringify(mm)}`);
    expect((res.facets['n'] as _DG.DomainFacetCountResult).count, 5); // three quantified + two path parents
  });

  test('FK-path filter: dotted tree in query, path categories facet', async () => {
    // Query item_event through the FK path (compiles to EXISTS server-side).
    const rows = await events().query(
      {filter: [{property: 'item_id.name', operator: '=', value: [`${marker} P1`]}]});
    expect(rows.length, 2, `expected the two P1 events, got ${rows.length}`);
    expect(rows.every((r: any) => r.kind === 'e1'), true);

    // Path categories facet: groups by the referenced row's name.
    const res = await events().facets({facets: [
      {id: 'parents', kind: 'categories', column: 'item_id.name', search: marker},
    ]});
    const byValue: {[v: string]: any} = {};
    for (const c of res.facets['parents'].categories)
      byValue[c.value] = c;
    expect(byValue[`${marker} P1`]?.total, 2, JSON.stringify(res.facets['parents']));
    expect(byValue[`${marker} P2`]?.total, 1);
  });

  test('saved filters: save, list, update in place, share, delete', async () => {
    const name = `zz apitests preset ${marker}`;
    const states = {'Name': {type: 'categorical', column: 'name', selected: [`${marker} red`]}};
    const saved = await items().filters.save(name, states);
    expect(saved.id != null, true, `no id on the saved preset: ${JSON.stringify(saved)}`);
    try {
      const listed = (await items().filters.list()).find((f) => f.id === saved.id);
      expect(listed != null, true, 'saved preset is not listed');
      expect((listed!.states['Name'] as any).selected[0], `${marker} red`);

      // Same id → in-place update: one entity, new states.
      const updated = await items().filters.save(name, {
        'Name': {type: 'categorical', column: 'name', selected: [`${marker} blue`]}}, {id: saved.id});
      expect(updated.id, saved.id);
      const relisted = (await items().filters.list()).filter((f) => f.friendlyName === name);
      expect(relisted.length, 1, 'same-name save must not create a second preset');
      expect((relisted[0].states['Name'] as any).selected[0], `${marker} blue`);

      // Presets are ordinary entities — the standard sharing machinery applies.
      const [entity] = await grok.dapi.getEntities([saved.id]);
      const allUsers = await grok.dapi.groups.filter('friendlyName = "All users"').first();
      await grok.dapi.permissions.grant(entity, allUsers, false);
      const perms = await grok.dapi.permissions.get(entity);
      expect((perms as any).view.some((g: _DG.Group) => g.friendlyName === 'All users'), true,
        'the granted group must appear in the view permissions');
    } finally {
      await items().filters.delete(saved.id);
    }
    const gone = (await items().filters.list()).some((f) => f.id === saved.id);
    expect(gone, false, 'deleted preset is still listed');
  });
}, {owner: 'askalkin@datagrok.ai'});
