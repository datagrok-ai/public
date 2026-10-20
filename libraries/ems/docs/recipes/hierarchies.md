# Recipe: a self-referencing table as a tree

Locations, categories, an org chart, a folder structure: one table whose rows point at rows of the
same table. Declare it in `schema.json` and three things follow — a row's ancestors, a subtree
filter term, and a `VirtualTree` that drives any other collection. Reference implementation:
**Stockroom**'s `Locations` app (`packages/Stockroom`), the location tree beside the containers
stored anywhere under the selected node.

## Declare it

```jsonc
"location": {
  "securityMode": "table",
  "businessKey": ["site", "name"],
  "friendlyName": "Locations",
  "hierarchy": true,                                          // ← the whole declaration
  "columns": {
    "name":      {"type": "string", "required": true, "isName": true, "searchable": true},
    "kind":      {"type": "string", "required": true, "choices": ["site", "building", "room", "cabinet"]},
    "parent_id": {"type": "ref", "ref": "location", "onDelete": "setnull", "friendlyName": "Parent"},
    "site":      {"type": "string", "required": true}
  }
}
```

`hierarchy: true` requires **exactly one** ref column targeting the table itself — that column
becomes the registry's `parentColumn`. Two self-refs, or none, and the schema is refused with
`invalid-hierarchy` at deploy time. A root row is one whose parent column is empty; nothing
enforces a single root, and a cycle does not hang anything (see below).

Two things the registry now reports, which everything else reads:

```ts
const locations = await domains.table('stockroom.location');
locations.info.hierarchy;      // true
locations.info.parentColumn;   // 'parent_id'
```

## The two server surfaces

**Ancestors** — a row's chain of parents, **root first**, and **without the row itself** (it is the
breadcrumb *in front of* the row, so a root answers `[]`):

```ts
const path = await grok.dapi.domains.table('stockroom.location').pathTo(shelf.id);
path.map((p) => p.name).join(' > ');     // 'Main campus > Building A > Lab 101'
```

**`under`** — a filter term matching a node's whole subtree, the node included. It reads through
the tree table's own `id`, or through any ref column pointing into the tree:

```ts
await locationClient.count(`id under "${root.id}"`);              // the nodes of the subtree
await containerClient.count(`location_id under "${root.id}"`);    // everything stored under it
```

Both walk the caller's View predicate at every level, which is what makes them safe to expose: an
ancestor the caller cannot see truncates the chain (no oracle on what is above it), a subtree they
cannot see is simply not matched, a cycle terminates, and the walk stops at depth 64. `under` is
server-only — it has no client-side mask, so a `dataframe`-scoped filter reports it as
`not-expressible` rather than silently answering the wrong rows.

The filter panel offers **is under** only where it means something: on a ref column whose target
table declares `hierarchy`, and on `id` of a hierarchy table.

## The tree, and the collection it drives

```ts
const tree = domains.tree(locations, {expandTo: row.id});
tree.selected;        // ReadonlySignal<RowView | null>
```

Roots are the rows whose parent column is null (`{property: parentColumn, operator: '=', value:
null}`), a branch loads its children when it is opened,
and both go through the plain `query` seam — no frame, no writer and no editor per node, which is
what makes a tree over a large table affordable. `expandTo` walks `ancestors(id)`, opens every one
of them and then selects the row. Nodes are labelled by the table's renderer caption, carry the
table's own row actions plus Open, and a table that is not a hierarchy is refused by name.

What a tree is FOR is driving a collection — the whole Stockroom `Locations` app:

```ts
//name: Locations
//tags: app
//output: view result
export async function stockroomLocations(): Promise<DG.ViewBase> {
  const [locations, containers] = await Promise.all([
    domains.table('stockroom.location'), domains.table('stockroom.container')]);
  const session = new SharedSession();
  const source = SharedSession.runWith(session, () =>
    containers.source({session, live: true, pageSize: 100}));
  const tree = domains.tree(locations);
  const table = domains.dataTable(source);
  tree.effect(() => {
    const node = tree.selected.value;
    source.query.value = node === null ? '' : `location_id under "${node.id}"`;
  });
  return appView({
    name: 'Locations',
    content: new Splitter([tree, table], {direction: 'horizontal', sizes: [0.35, 0.65]}),
    own: [source],
    status: source.summary,
    path: '/apps/Stockroom/Locations',
  });
}
```

Three things to copy from it:

- **Nothing selected is the whole table.** Clearing the query, not leaving the last node's filter
  in force, is what a user expects from deselecting.
- **The right pane is a `domains.dataTable`** — virtualized HTML rows over the source, with the
  editor's dirty and invalid colours; the platform grid (`domains.grid`) is the one that edits
  in place, and this pane does not.
- **`live: true`** polls one aggregate (`count` + `max(updated_on)`) under the source's filter and
  search every 30 s, and refreshes while the page is clean — a container another session moves
  shows up here without a reload. See [crud-app](crud-app.md) for the stale/Refresh half.

`domains.tree` and `domains.dataTable` are spec-registered too (`u2-domain-tree`,
`u2-domain-data-table`), so the same page can be laid out as a `dg-ui/1` spec — see
[spec-app](spec-app.md).

## Anti-patterns

- A `parent_id` column without `hierarchy: true` — the ref works, but `pathTo`, `under`, the tree
  and the filter panel's "is under" all refuse the table, and every consumer re-implements the
  walk.
- Fetching the whole table and building the tree in the browser — the tree loads one level at a
  time on purpose; a `where parent_id = …` per open branch is what the server is for.
- A recursive client-side ancestor walk (`while (row.parent_id) await client.get(…)`) — that is
  one round trip per level and it leaks which rows exist above the caller's reach; `pathTo` is one
  request and truncates instead.
- Filtering the collection by the selected node alone (`location_id = "<id>"`) — the point of a
  tree is that picking a site answers everything below it, which is what `under` says.
