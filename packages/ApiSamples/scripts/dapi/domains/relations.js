// Many-to-many relations on a domain table: declared once in the manifest
// ("relations": {"tags": {"via": "item_tag", "target": "tag"}}) and then available
// as an expand key, a filter path, and an inline write value.
// In a package, `grok api` types them: <Table>Expand gets the link array,
// <Table>Insert / <Table>Update get the `tags?: string[]` link set.

if (!(await grok.dapi.domains.schemas.list()).some((s) => s.name === 'apitests'))
  return grok.shell.info('Deploy the ApiTests package first (it declares the apitests domain schema)');

const items = grok.dapi.domains.table('apitests.item');
const tags = grok.dapi.domains.table('apitests.tag');
const stamp = `REL-${Date.now()}`;

const [red, blue] = await tags.insert([{name: `${stamp} red`}, {name: `${stamp} blue`}]);

// WRITE: a relation travels as the full list of target ids, in the same
// transaction as the row itself.
const [item] = await items.insert({sku: `${stamp}-1`, name: 'Widget', tags: [red.id, blue.id]});

// READ: expand returns {id, name} links, capped at 100 and ordered by display name.
const [row] = await items.query({filter: [{property: 'id', operator: '=', value: item.id}],
  expand: ['tags']});
grok.shell.info(`tags: ${row.tags.map((t) => t.name).join(', ')}`);

// In a DataFrame the same links become a chips column plus a hidden id companion —
// '~tags.id' carries the ids in the same order and is the source of truth.
const df = await items.queryDf({filter: `sku starts "${stamp}"`, expand: ['tags']});
grok.shell.info(`d42 columns: ${df.columns.names().join(', ')}`);

// FILTER: a relation is one hop, so 'tags.name' and 'tags.id' are filter paths
// (an id list compiles to ANY-of — what a checkbox list emits).
const tagged = await items.query({filter: [{property: 'tags.id', operator: '=', value: [red.id]}]});
grok.shell.info(`${tagged.length} item(s) tagged red`);

// The smart-filter string takes the same paths (values quoted, as everywhere).
const byName = await items.query({filter: `tags.name = "${stamp} red"`});
grok.shell.info(`${byName.length} item(s) matched by the string filter`);

// Facet counts are per OWNER: an item with two tags counts once under each,
// and 'tags.id' also brings the target's display name back.
const facets = await items.facets({filter: `sku starts "${stamp}"`,
  facets: [{id: 't', kind: 'categories', column: 'tags.id'}]});
grok.shell.info(facets.facets['t'].categories.map((c) => `${c.display}: ${c.total}`).join(', '));

// UPDATE has SET semantics over the links you can SEE: the list becomes the link
// set, [] clears it, and an absent key leaves the relation untouched. Links whose
// junction or target row you cannot view are never diffed — a replace can only
// remove links its author knows about.
await items.update(item.id, {tags: [blue.id]});
grok.shell.info(`after replace: ${(await items.query({expand: ['tags'],
  filter: [{property: 'id', operator: '=', value: item.id}]}))[0].tags.length} link(s)`);

// Creating the target and linking it atomically is a transaction — '$ref'
// placeholders resolve element-wise inside the list.
const [fresh, linked] = await grok.dapi.domains.transaction('apitests', [
  {op: 'insert', table: 'tag', ref: 't', values: {name: `${stamp} fresh`}},
  {op: 'insert', table: 'item', values: {sku: `${stamp}-2`, name: 'Gadget', tags: ['$t', red.id]}},
]);

// Cleanup: deleting an owner takes its junction rows with it (cascade).
await items.delete(item.id);
await items.delete(linked.id);
for (const t of [red.id, blue.id, fresh.id])
  await tags.delete(t);
