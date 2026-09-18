// Validate an import before committing it: `batch(data, {validateOnly: true})` runs the WHOLE
// commit path on the server — coercions, per-row validation, required and auto-number checks,
// intra-batch and live business-key duplicates, FK existence, immutable columns, and the merge
// itself — then rolls the transaction back. So the verdicts are the ones a real commit would
// produce, not a second implementation that can drift from it.
//
// Nothing is written by a validation: no row, no audit entry, no auto-number, no notification,
// and no change to the table's version token. And no id on any row — the id of a rolled-back
// insert does not exist. The verdicts are a prediction against the table AS IT IS NOW; a
// concurrent write can change them before the commit.

const items = grok.dapi.domains.table('apitests.item');
const stamp = `IV-${Date.now()}`;

// One row that already exists (upsert will merge into it), two new ones, one the table refuses.
await items.insert({sku: `${stamp}-1`, name: 'Seed', quantity: 1});
const payload = [
  {sku: `${stamp}-1`, name: 'Existing', quantity: 9},
  {sku: `${stamp}-2`, name: 'New A', quantity: 2},
  {sku: `${stamp}-3`, name: 'New B', quantity: 3},
  {sku: `${stamp}-4`, name: 'Refused', quantity: -5},
];

// `allOrNothing` is forced false for a preview — it judges every row, not just the first bad one.
const preview = await items.batch(payload, {mode: 'upsert', validateOnly: true});
grok.shell.info(`${preview.rowCount} rows: ${preview.willInsert} insert, ${preview.willUpdate} update, ` +
  `${preview.willSkip} skip, ${preview.errorCount} error`);

// The per-row report, capped server-side at 1000 and ordered errors -> skips -> updates -> inserts.
// `predicted` is 'insert' | 'update' | 'skip' | 'error'; `existingId` names the row an update or
// a skip would land on.
const report = DG.DataFrame.fromObjects(preview.rows.map((r) => ({
  row: r.index,
  predicted: r.predicted,
  existing: r.existingId ?? null,
  problem: (r.errors ?? []).map((e) => `${e.column}: ${e.message}`).join('; '),
})));
grok.shell.addTableView(report);

// Cleanup: one filtered bulk delete (bounded; loop while hasMore for larger sets).
const cleanup = async () => {
  for (let guard = 0; guard < 100; guard++)
    if (!(await items.deleteWhere(`sku starts "${stamp}"`)).hasMore)
      return;
};

// Let the user decide with the verdicts in front of them, then commit the SAME payload.
// `allOrNothing: false` applies the good rows and reports the bad ones — the counts come back
// equal to the preview's unless someone else wrote in between.
ui.dialog('Import preview')
  .add(ui.divText(`${preview.willInsert + preview.willUpdate} row(s) will be applied, ` +
    `${preview.errorCount} will fail. Import?`))
  .onOK(async () => {
    const committed = await items.batch(payload, {mode: 'upsert', allOrNothing: false});
    grok.shell.info(`committed: ${committed.inserted} inserted, ${committed.updated} updated, ` +
      `${committed.errorCount} failed`);
    await cleanup();
  })
  .onCancel(cleanup)
  .show();
