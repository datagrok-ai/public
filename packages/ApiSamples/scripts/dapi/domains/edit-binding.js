//api: DG.DomainSchemaClient.manifest, DG.DomainSchemaClient.apply, DG.DomainSchemaClient.delete
// Editing a live domain schema: an apply on a user-managed schema names the version it was
// edited against, and may change the schema's friendlyName and description — and, on an
// external binding, whether it is writable. Everything else about a binding's storage
// (connection, remote schema, catalog, kind) is immutable: rebinding is a new schema.
// Runs over a throwaway platform-stored schema; the storage refusal is shown on it too.

const name = `ed_${Math.random().toString(36).slice(2, 8)}`;
await grok.dapi.domains.createSchema(name, {friendlyName: 'Edit probe'});
const handle = grok.dapi.domains.schema(name);
try {
  // 1. Every apply carries ifVersion: the manifest's version at the time of the edit.
  //    Without it the server answers 'version-required'; a stale one is a
  //    DomainVersionConflictError — of two concurrent applies exactly one commits.
  let {version} = await handle.manifest();
  await handle.apply({tables: {note: {columns: {text: {type: 'string'}}}}, ifVersion: version});

  // 2. Metadata travels in the same apply; a retained table keeps its grants exactly
  //    (bootstrap grants go only to what an apply creates).
  ({version} = await handle.manifest());
  await handle.apply({friendlyName: 'Edit probe, renamed', description: 'edited from a script', ifVersion: version});
  grok.shell.info(`renamed: ${(await grok.dapi.domains.schemas.list()).find((s) => s.name === name).friendlyName}`);

  // 3. `storage.writable` is the one storage key an apply may change, and only on an
  //    external binding; here the schema is platform-stored, so it is refused by name.
  ({version} = await handle.manifest());
  try {
    await handle.apply({storage: {writable: true}, ifVersion: version});
  } catch (e) {
    grok.shell.info(`refused: ${e.code}`);   // 'storage-immutable'
  }

  // 4. A dry run of a drop says what goes with the table — counts of direct grants,
  //    column restrictions, promoted rows and saved filters — before you confirm it.
  const plan = await handle.apply({dropTables: ['note'], ifVersion: version}, {dryRun: true});
  grok.shell.info(`dropping note loses ${JSON.stringify(plan.lost.tables.note)}`);
} finally {
  await handle.delete();
}
