// A table that declares `"hierarchy": true` in schema.json has exactly one ref column
// pointing at itself, and two surfaces over that tree: `pathTo(id)` — the row's
// ancestors, root first and the row itself excluded (the breadcrumb in front of it) —
// and the `under` filter term, which matches a node's whole subtree (through the tree
// table's own id, or through any ref column into it). Both walk the caller's View
// predicate at every level: an ancestor you cannot see truncates the chain, and a
// subtree you cannot see is simply not matched.
// The demo tree here is stockroom.location (site › room › shelf).

const deployed = (await grok.dapi.domains.schemas.list()).some((s) => s.name === 'stockroom');
const info = deployed ? await grok.dapi.domains.registry.tableInfo('stockroom.location') : null;

if (info == null || !info.hierarchy)
  grok.shell.info('stockroom.location is not deployed here, or does not declare "hierarchy": true');
else {
  const locations = grok.dapi.domains.table('stockroom.location');
  grok.shell.info(`${info.pluralName}: a tree over '${info.parentColumn}'`);

  // The breadcrumb of one nested node: its ancestors, root first ([] for a root row).
  const leaf = await locations.first({filter: `${info.parentColumn} != null`});
  const path = await locations.pathTo(leaf.id);
  grok.shell.info(`${path.map((p) => p.name).join(' > ')} > ${leaf.name}`);

  // The subtree of its root — that node and everything below it.
  const root = path.length === 0 ? leaf : path[0];
  const subtree = await locations.query({filter: `id under "${root.id}"`});
  grok.shell.info(`${subtree.length} locations under ${root.name}`);

  // The same term through a ref column INTO the tree: every container stored anywhere
  // below that node — what a tree node click filters the collection by.
  const hasContainers = (await grok.dapi.domains.registry.tableInfo('stockroom.container')
    .catch(() => null)) != null;
  if (!hasContainers)
    grok.shell.info('stockroom.container is not deployed here — skipping the ref-column half');
  else {
    const containers = grok.dapi.domains.table('stockroom.container');
    grok.shell.info(`${await containers.count(`location_id under "${root.id}"`)} containers under ${root.name}`);
  }
}
