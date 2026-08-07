// Domain registry reflection + permission-driven capabilities: the runtime metadata
// reflective UI components consume to work on ANY domain table without codegen.

const table = 'grit.issue'; // any '<schema>.<table>' registered on this server

// Property objects for the table's declared columns — the same metadata (type,
// semType, choices, min/max, nullable) that drives the built-in row editor and
// the server-side validation, so a form built from these matches the server.
const props = await grok.dapi.domains.registry.rowProperties(table);
grok.shell.info(props.map((p) => `${p.name}: ${p.propertyType}` +
  (p.choices?.length ? ` ${JSON.stringify(p.choices)}` : '')).join('\n'));

// Display identity, security mode, and the FK-inverted child tables that drive
// detail-table links ('comment' and 'issue_label' reference grit.issue).
const info = await grok.dapi.domains.registry.tableInfo(table);
grok.shell.info(`name column: ${info.nameColumn}, business key: ${info.businessKey}\n` +
  `children: ${info.childTables.map((c) => `${c.schema}.${c.table} via ${c.fkColumn} ("${c.label}")`).join(', ')}`);

// Effective capabilities of the CURRENT user — server-truth permission probes on
// the securing entity plus writable columns under column security. Gate UI
// affordances on this object instead of hand-wiring permission checks; it is
// cached until a grant change (grant/revoke through the client invalidate it,
// call grok.dapi.domains.invalidateUiCaches() after out-of-band changes).
const caps = await grok.dapi.domains.table(table).capabilities();
grok.shell.info(`insert: ${caps.canInsert}, edit: ${caps.canEdit}, delete: ${caps.canDelete}\n` +
  `writable columns: ${caps.writableColumns.join(', ')}`);

// Batched display-name resolution: ids → 'name column value → business key → id'.
// Requests coalesce into one fetch per table — cheap to call per rendered cell.
const [row] = await grok.dapi.domains.table(table).query({limit: 1});
if (row != null) {
  const names = await grok.dapi.domains.registry.resolveNames(table, [row.id]);
  grok.shell.info(`${row.id} → ${names[row.id]}`);
}
