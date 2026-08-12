// The platform's row dialogs and panes as one-line openers: create/edit/clone,
// delete confirmation, the lookup picker, the 409 conflict dialog, the audit
// history and the grants pane. Each one enters the SAME flow the built-in
// Domain View uses — registry-driven inputs, server-parity validation, column
// security, promote-on-share — so an app never rebuilds any of them.

const issues = grok.dapi.domains.table('grit.issue');
const handler = new DG.DomainObjectHandler('grit.issue');

// Create: resolves to whether a row was saved.
const created = await DG.DomainObjectHandler.createRow('grit.issue');
grok.shell.info(created ? 'Issue created' : 'Cancelled');

// The lookup picker: the target table's Domain View in a dialog (null on cancel).
const project = await DG.DomainObjectHandler.pickRow('grit.project');

const [first] = await issues.query().orderBy('created_on', true).top(1);
if (first == null)
  throw new Error('grit.issue has no rows — create one first (/domains/grit/issue).');
const row = await handler.getById(first.id);                 // DG.DomainRow

// Edit / clone / delete — the instance members delegate to the same statics,
// so a subclass replaces a flow by overriding one method.
await handler.editRow(row);                                  // true when saved
// await handler.cloneRow(row); await handler.deleteRow(row);

// The standard optimistic-concurrency dialog, for code that saves by itself:
// resolves 'reload' | 'overwrite' | null, and the caller applies the decision.
const stale = {id: row.id, version: row.version, title: row.values.title};
await issues.save({...stale});                 // bumps the version...
try {
  await issues.save(stale);                    // ...so this one is stale: 409
} catch (e) {
  if (!(e instanceof DG.DomainVersionConflictError))
    throw e;
  const decision = await DG.DomainObjectHandler.showConflictDialog(row.semValue);
  if (decision === 'overwrite')
    await issues.save({...stale, version: e.currentVersion});
  else if (decision === 'reload')
    grok.shell.info(`Server version is ${e.currentVersion}`);
}

// Embeddable panes: the row's audit trail, and the grants of the registry
// entity that governs the table (rows themselves share via handler.shareRow).
const table = (await grok.dapi.domains.schemas.include('tables').filter('name = "grit"').first())
  .tables.find((t) => t.name === 'issue');
grok.shell.newView('grit.issue dialogs', [
  ui.h2(handler.getCaption(row)),
  ui.divText(`Picked project: ${project?.semValue ?? '(none)'}`),
  ui.h3('History'), DG.DomainObjectHandler.auditPane(row),
  ui.h3('Access'), DG.DomainObjectHandler.grantsPane(table.id, 'grit.issue'),
]);
