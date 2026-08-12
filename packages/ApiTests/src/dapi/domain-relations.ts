import type * as _DG from 'datagrok-api/dg';
declare let DG: typeof _DG;

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {apitestsDb, ItemColumn, ItemUpdate} from '../generated/db';

// Many-to-many relations (GROK-20298): the 'tags' relation this package declares on
// apitests.item (via the master-mode junction apitests.item_tag, target apitests.tag)
// through the JS API — expand, relation filter paths, inline set-semantics writes, the
// create-and-link transaction, and the batch refusal. Driven through the `grok api`
// generated client (src/generated/db.ts), so the generated relation types compile here.
// The read fixture is built once in before(); every write test links tags of its own, so
// the read assertions hold whatever order the tests run in.
category('Dapi: domain relations', () => {
  // The two typed clients, resolved once in before() — the schema handle builds a
  // fresh one on every property access.
  let items: typeof apitestsDb.items;
  let tags: typeof apitestsDb.tags;
  const marker = `REL${Date.now()}${Math.floor(Math.random() * 1e4)}`;
  const sku = () => `${marker}-${Math.floor(Math.random() * 1e9)}`;
  const itemIds: string[] = [];
  const tagIds: string[] = [];
  // red/green/blue belong to the read fixture and are never linked anywhere else;
  // the write tests link wtag1/wtag2 so they cannot perturb the counts above.
  let red: string;
  let blue: string;
  let green: string;
  let wtag1: string;
  let wtag2: string;
  let itemA: string; // linked to red and blue
  let itemG: string; // linked to green
  let itemBare: string; // no links at all

  // Display names of one item's links, as the expand returns them.
  async function linked(id: string): Promise<string[]> {
    const rows = await items.query().where('id', '=', id).expand('tags');
    return (rows[0].tags ?? []).map((l) => l.name).sort();
  }

  before(async () => {
    items = apitestsDb.items;
    tags = apitestsDb.tags;
    const ins = await tags.insert([{name: `${marker} red`}, {name: `${marker} blue`},
      {name: `${marker} green`}, {name: `${marker} wtag1`}, {name: `${marker} wtag2`}]);
    for (const r of ins)
      tagIds.push(r.id);
    [red, blue, green, wtag1, wtag2] = tagIds;

    // The read fixture: A carries a duplicate id (deduped on write), G one link, bare none.
    const rows = await items.insert([
      {sku: sku(), name: `${marker} A`, tags: [red, blue, red]},
      {sku: sku(), name: `${marker} G`, tags: [green]},
      {sku: sku(), name: `${marker} bare`},
    ]);
    for (const r of rows)
      itemIds.push(r.id);
    [itemA, itemG, itemBare] = itemIds;
  });

  after(async () => {
    for (const id of itemIds)
      await items.delete(id).catch(() => {}); // cascades its junction rows
    for (const id of tagIds)
      await tags.delete(id).catch(() => {});
  });

  test('insert with links, expand, and the d42 chips + id companion', async () => {
    // JSON expand: {id, name} links, deduped, ordered by display name, never null.
    const rows = await items.query().where('id', '=', itemA).expand('tags');
    const links = rows[0].tags!;
    expect(links.length, 2, `expected two deduped links, got ${JSON.stringify(links)}`);
    expect(links.map((l) => l.name).join(','), `${marker} blue,${marker} red`);
    expect(links.every((l) => l.id === red || l.id === blue), true, JSON.stringify(links));

    // An item with no links is [], not null.
    expect((await linked(itemBare)).length, 0);

    // d42: the display names joined by ', ' plus the '~tags.id' companion in the SAME
    // order — the ids are the source of truth for in-grid editing.
    const df = await items.queryDf({filter: [{property: 'id', operator: '=', value: itemA}],
      expand: ['tags']});
    const names = df.columns.byName('tags');
    const ids = df.columns.byName('~tags.id');
    expect(names != null && ids != null, true, `columns: ${df.columns.names().join(', ')}`);
    expect(names.get(0), `${marker} blue, ${marker} red`);
    expect(ids.get(0).split(', ').length, 2, `id companion: ${ids.get(0)}`);
    expect(ids.get(0).split(', ').includes(red), true);

    // No links in d42 is EMPTY, in both columns — a DataFrame string column has no
    // null slot of its own, so this is what "null" looks like there. Consumers must
    // test emptiness (a chips renderer must not paint one empty chip).
    const bare = await items.queryDf({filter: [{property: 'id', operator: '=', value: itemBare}],
      expand: ['tags']});
    expect(bare.columns.byName('tags').isNone(0), true, `bare tags: "${bare.columns.byName('tags').get(0)}"`);
    expect(bare.columns.byName('~tags.id').get(0), '');
  });

  test('filter: a relation segment and an id ANY-of list', async () => {
    const marked: _DG.DomainConditionTree<ItemColumn> =
      [{property: 'name', operator: 'like', value: `${marker}%`}];

    // 'tags.name' matches the owners that HAVE such a link (per-hop EXISTS).
    const byName = await items.query({filter: [...marked, 'and',
      {property: 'tags.name', operator: '=', value: `${marker} red`}]});
    expect(byName.length, 1, `expected the one red-tagged item, got ${byName.length}`);
    expect(byName[0].id, itemA);

    // The smart-filter STRING form takes the same paths (values quoted, as everywhere).
    const byString = await items.query({filter: `name like "${marker}" and tags.name = "${marker} red"`});
    expect(byString.length, 1, `string-form relation filter: got ${byString.length}`);
    expect(byString[0].id, itemA);

    // 'tags.id' with a list is the ANY-of shape the filter panel's checkbox list emits.
    const byIds = await items.query({filter: [...marked, 'and',
      {property: 'tags.id', operator: '=', value: [red, green]}], sort: 'name'});
    expect(byIds.map((r) => r.id).join(','), `${itemA},${itemG}`);
    // An empty selection matches nothing (never "everything").
    expect((await items.query({filter: [...marked, 'and',
      {property: 'tags.id', operator: '=', value: []}]})).length, 0);

    // '<relation>.id = null' is the facet's "(no value)" bucket: owners with NO
    // visible link (a NOT EXISTS — the one question per-hop EXISTS cannot ask).
    // Membership, not a count: the write tests leave unlinked items of their own.
    const unlinked = (await items.query({filter: [...marked, 'and',
      {property: 'tags.id', operator: '=', value: null}]})).map((r) => r.id);
    const anyLink = (await items.query({filter: [...marked, 'and',
      {property: 'tags.id', operator: '!=', value: null}]})).map((r) => r.id);
    expect(unlinked.includes(itemBare), true, 'the unlinked item belongs to the null bucket');
    expect(unlinked.includes(itemA) || unlinked.includes(itemG), false, JSON.stringify(unlinked));
    expect(anyLink.includes(itemA) && anyLink.includes(itemG), true, JSON.stringify(anyLink));
    expect(anyLink.includes(itemBare), false, 'the two sides never overlap');

    // '!=' with an id list is the facet's UNCHECK gesture: it EXCLUDES the owners
    // linked to any of them. The per-hop reading ("has some other tag") would keep
    // itemA — it carries blue as well — right after the user unchecked red.
    const notRed = (await items.query({filter: [...marked, 'and',
      {property: 'tags.id', operator: '!=', value: [red]}]})).map((r) => r.id);
    expect(notRed.includes(itemA), false, `unchecking red must drop itemA: ${JSON.stringify(notRed)}`);
    expect(notRed.includes(itemG) && notRed.includes(itemBare), true, JSON.stringify(notRed));
    // Excluding nothing constrains nothing.
    expect((await items.query({filter: [...marked, 'and',
      {property: 'tags.id', operator: '!=', value: []}]})).length >= 3, true);

    // Unchecking red AND the "(no value)" bucket — what the filter panel emits
    // for that gesture (an exclusion list minus null, AND '!= null'). The extra
    // term is redundant on a plain column but load-bearing here: the id leaf's
    // '!=' is a NOT EXISTS, true for an owner with no links at all.
    const notRedNorBare = (await items.query({filter: [...marked, 'and',
      {property: 'tags.id', operator: '!=', value: [red]}, 'and',
      {property: 'tags.id', operator: '!=', value: null}]})).map((r) => r.id);
    expect(notRedNorBare.includes(itemBare), false,
      `the null bucket must be excluded too: ${JSON.stringify(notRedNorBare)}`);
    expect(notRedNorBare.includes(itemA), false, JSON.stringify(notRedNorBare));
    expect(notRedNorBare.includes(itemG), true, JSON.stringify(notRedNorBare));

    // Facets count DISTINCT OWNERS per target, with the target's display name.
    const res = await items.facets({filter: marked, facets: [
      {id: 't', kind: 'categories', column: 'tags.id', search: marker}]});
    const byValue: {[v: string]: any} = {};
    for (const c of (res.facets['t'] as _DG.DomainFacetCategoriesResult).categories)
      byValue[c.value] = c;
    expect(byValue[red]?.total, 1, JSON.stringify(res.facets['t']));
    expect(byValue[red]?.display, `${marker} red`, JSON.stringify(res.facets['t']));
    expect(byValue[green]?.total, 1);
  });

  test('update: the list replaces the link set, [] clears it, absent leaves it alone', async () => {
    const [ins] = await items.insert({sku: sku(), name: `${marker} P`, tags: [wtag1]});
    itemIds.push(ins.id);

    // Set-replace: wtag1 out, wtag2 in — the generated <Table>Update names the relation.
    const values: ItemUpdate = {tags: [wtag2]};
    const upd = await items.update(ins.id, values);
    expect((await linked(ins.id)).join(','), `${marker} wtag2`);
    // A relations-only change bumps the row version, so expectedVersion covers the links.
    expect(upd.version, 2, 'a relations-only update must bump the version');

    // Re-sending the same set is a no-op, an ordinary column rides along, and an
    // ABSENT relation key leaves the links untouched.
    await items.update(ins.id, {tags: [wtag2]});
    await items.update(ins.id, {quantity: 7});
    expect((await linked(ins.id)).length, 1, 'an absent relation key must not clear the links');

    // [] clears the visible link set.
    await items.update(ins.id, {tags: []});
    expect((await linked(ins.id)).length, 0);

    // A target that is missing (or invisible — one code, no existence oracle)
    // fails the whole write, and leaves the link set as it was.
    let err: any = null;
    try {
      await items.update(ins.id, {tags: ['3f2504e0-4f89-11d3-9a0c-0305e82c3301']});
    } catch (e) {
      err = e;
    }
    expect(err instanceof DG.DomainValidationError, true,
      `expected DomainValidationError, got ${err?.constructor?.name}: ${err?.message}`);
    expect(err.rows[0].errors[0].code, 'not-visible-or-missing', JSON.stringify(err.body));
    expect(err.rows[0].errors[0].column, 'tags');
    expect((await linked(ins.id)).length, 0, 'a failed link write must roll back');
  });

  test('transaction: create the target and link it in one atomic operation', async () => {
    const name = `${marker} fresh`;
    const [tag, item] = await apitestsDb.transaction([
      {op: 'insert', table: 'tag', ref: 't1', values: {name: name}},
      {op: 'insert', table: 'item', values: {sku: sku(), name: `${marker} TX`, tags: ['$t1', wtag1]}},
    ]);
    tagIds.push(tag.id);
    itemIds.push(item.id);
    expect((await linked(item.id)).join(','), `${marker} fresh,${marker} wtag1`);
  });

  test('batch refuses relation keys and points at insert/update', async () => {
    let message = '';
    try {
      const report = await items.batch([{sku: sku(), name: `${marker} B`, tags: [wtag1]}]);
      message = report.error ?? '';
    } catch (e: any) {
      message = e.message ?? `${e}`;
    }
    expect(message.includes('not supported in batch mode'), true, `unexpected outcome: '${message}'`);
    expect(message.includes('insert or patch'), true, `unexpected outcome: '${message}'`);
  });
}, {owner: 'askalkin@datagrok.ai'});
