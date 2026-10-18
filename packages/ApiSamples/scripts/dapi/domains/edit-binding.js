//api: DG.DomainSchemaClient.manifest, DG.DomainSchemaClient.access, DG.DomainSchemaClient.apply, DG.DomainSchemaClient.draft, DG.DomainSchemaClient.delete, DG.DomainVersionConflictError
// Editing a live domain schema: an apply names the version AND the incarnation it was edited
// against, and may change the schema's friendlyName, description — and, on an external binding,
// whether it is writable. Access changes ride the SAME apply as permission-triple deltas, built
// from the access snapshot, and land in one transaction with the manifest change. Everything else
// about a binding's storage (connection, remote schema, catalog, kind) is immutable: rebinding is
// a new schema. Runs over a throwaway platform-stored schema; the binding-only parts
// (`storage.writable`, `draft()`) are shown as the refusals they answer on it.

const name = `ed_${Math.random().toString(36).slice(2, 8)}`;
await grok.dapi.domains.createSchema(name, {friendlyName: 'Edit probe'});
const handle = grok.dapi.domains.schema(name);
const tokens = async () => {
  const {version, incarnation} = await handle.manifest();
  return {ifVersion: version, ifIncarnation: incarnation};
};
try {
  // 1. Every apply carries ifVersion: the manifest's version at the time of the edit; without it
  //    the server answers 'version-required'. ifIncarnation guards against the name having been
  //    deleted and re-created meanwhile. A stale token is a DomainVersionConflictError — of two
  //    concurrent applies exactly one commits.
  await handle.apply({...await tokens(), tables: {note: {columns: {text: {type: 'string'}, secret: {type: 'string'}}}}});

  // 2. Metadata travels in the same apply; a retained table keeps its grants exactly
  //    (bootstrap grants go only to what an apply creates).
  await handle.apply({...await tokens(), friendlyName: 'Edit probe, renamed', description: 'edited from a script'});
  grok.shell.info(`renamed: ${(await grok.dapi.domains.schemas.list()).find((s) => s.name === name).friendlyName}`);

  // 3. The access snapshot: per table its direct grants (complete permission sets), per column
  //    whether it is restricted and to whom; null where you may not Share the target.
  const access = await handle.access();
  grok.shell.info(`note: ${access.tables.note.grants.length} grantee(s), secret ${access.columns['note.secret'].state}`);

  // 4. A delta built from the snapshot, dry-run first: every op with the effect it would have
  //    ('none' where the triple is already there), nothing moved. Then the commit.
  const allUsers = await grok.dapi.groups.filter('friendlyName = "All users"').first();
  const delta = {
    grant: [{table: 'note', group: allUsers.id, permission: 'View'}],
    restrict: [{table: 'note', column: 'secret', grant: [{group: allUsers.id, permission: 'View'}]}],
  };
  const plan = await handle.apply({...await tokens(), access: delta}, {dryRun: true});
  grok.shell.info(`planned: grant ${plan.access.grant[0].effect}, restrict ${plan.access.restrict[0].effect}`);
  const applied = await handle.apply({...await tokens(), access: delta});
  grok.shell.info(`applied: ${applied.applied}, secret is now ${(await handle.access()).columns['note.secret'].state}`);

  //    Making a column visible to everyone again is no delta: the op names the state and the ACL
  //    revision the snapshot answered, so a grant or revoke made since is an 'access-conflict'.
  const secret = (await handle.access()).columns['note.secret'];
  await handle.apply({...await tokens(), access: {unrestrict: [{table: 'note', column: 'secret',
    from: 'restricted', revision: secret.revision}]}});

  // 5. A stale token: the conflict names both versions and both incarnations (present whenever
  //    the apply sent ifIncarnation; they differ only after a delete and a re-create of the name).
  const stale = await tokens();
  await handle.apply({...stale, description: 'moved on'});
  try {
    await handle.apply({...stale, description: 'never lands'});
  } catch (e) {
    grok.shell.info(`conflict: expected v${e.expectedVersion}, current v${e.currentVersion}`);
  }

  // 6. `storage.writable` is the one storage key an apply may change, and only on an external
  //    binding; `draft()` reads a binding's live catalog for its editor. Here the schema is
  //    platform-stored, so both are refused by name.
  try {
    await handle.apply({...await tokens(), storage: {writable: true}});
  } catch (e) {
    grok.shell.info(`refused: ${e.code}`);   // 'storage-immutable'
  }
  try {
    await handle.draft();
  } catch (e) {
    grok.shell.info(`refused: ${e.code}`);   // 'unsupported'
  }

  // 7. A dry run of a drop says what goes with the table — counts of direct grants, column
  //    restrictions, promoted rows and saved filters — before you confirm it.
  const drop = await handle.apply({...await tokens(), dropTables: ['note']}, {dryRun: true});
  grok.shell.info(`dropping note loses ${JSON.stringify(drop.lost.tables.note)}`);
} finally {
  // A delete pinned to the incarnation you looked at never removes a schema re-created since.
  await handle.delete({ifIncarnation: (await handle.manifest()).incarnation});
}
