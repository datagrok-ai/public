import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {DomainFrameEditor, DomainGrid, SERVICE_COLUMNS, STATE_COLUMN, CHANGES_COLUMN,
  ERRORS_COLUMN, validateCellValue} from '@datagrok-libraries/domain-ui';
import {withRestrictedUser} from './domain-lifecycle';

// ui-js-api WO-7: @datagrok-libraries/domain-ui — DomainFrameEditor (THE single
// writer of the '~state'/'~changes'/'~errors' service columns) and DomainGrid.
// The state-transition and op-builder tests run alongside the live loop here, as
// planned: they import the library, which only ApiTests bundles. NONE of them is
// pure — an editor probes the registry and the caller's capabilities on create,
// so every test below needs a live server with the apitests schema registered.
//
// Fixture: the package's own 'apitests.item' (sku required+unique, name,
// quantity int min 0). Every test inserts inside its try and deletes its rows in
// the finally with ONE filtered deleteWhere.
category('Dapi: domain frame editor', () => {
  const items = () => grok.dapi.domains.table('apitests.item');
  const stamp = () => `${Date.now()}-${Math.floor(Math.random() * 1e6)}`;

  /** Inserts [n] fixture rows sharing a sku prefix; resolves to their ids. */
  async function seed(prefix: string, n: number): Promise<string[]> {
    const rows = [];
    for (let i = 0; i < n; i++)
      rows.push({sku: `${prefix}-${i}`, name: `Item ${i}`, quantity: i});
    return (await items().insert(rows)).map((r) => r.id);
  }

  /** Drops every row of the fixture prefix (rows the test inserted AND rows a
   * save inserted through the editor). */
  async function cleanup(prefix: string): Promise<void> {
    try {
      await items().deleteWhere({property: 'sku', operator: 'like', value: `${prefix}%`});
    } catch (e) {
      console.error(`frame-editor fixture ${prefix} not cleaned up: ${e}`);
    }
  }

  /** An editor over exactly the rows of [prefix], oldest first. */
  async function editorFor(prefix: string): Promise<DomainFrameEditor> {
    const query = {filter: {property: 'sku', operator: 'like', value: `${prefix}%`} as any, sort: 'sku'};
    return await DomainFrameEditor.create(items() as any, {query: query});
  }

  test('service columns: attached, tagged, and absent from csv/binary exports', async () => {
    const prefix = `fe-cols-${stamp()}`;
    await seed(prefix, 2);
    try {
      const editor = await editorFor(prefix);
      const df = editor.dataFrame;
      for (const name of SERVICE_COLUMNS) {
        const col = df.columns.byName(name);
        expect(col != null, true, `${name} was not attached`);
        expect(col!.meta.includeInBinaryExport, false, `${name} is not excluded from binary export`);
        expect(col!.meta.includeInCsvExport, false, `${name} is not excluded from csv export`);
      }

      // A DIRTY frame is the interesting one: state present, and still invisible
      // to every writer that could persist or upload it.
      editor.setValue(0, 'name', 'Renamed');
      editor.addRow({sku: `${prefix}-new`, name: 'Added', quantity: 1});
      editor.markDeleted(1);
      expect(editor.isDirty, true, 'the editor did not go dirty');

      const csv = df.toCsv();
      for (const name of SERVICE_COLUMNS)
        expect(csv.includes(name), false, `${name} leaked into toCsv()`);

      const roundTrip = DG.DataFrame.fromByteArray(df.toByteArray());
      for (const name of SERVICE_COLUMNS)
        expect(roundTrip.columns.contains(name), false, `${name} leaked into toByteArray()`);
      // The data columns DID survive — otherwise the assertion above is vacuous.
      expect(roundTrip.columns.contains('sku'), true, 'the binary round trip lost the data columns');

      editor.detach();
    } finally {
      await cleanup(prefix);
    }
  });

  test('state transitions: edit, revert, re-edit, add, delete, undelete', async () => {
    const prefix = `fe-state-${stamp()}`;
    await seed(prefix, 2);
    try {
      const editor = await editorFor(prefix);
      const df = editor.dataFrame;
      const original = df.get('name', 0);

      editor.setValue(0, 'name', 'Edited once');
      expect(editor.stateOf(0), 'modified', 'the row did not become modified');
      expect(editor.isChanged(0, 'name'), true, 'the cell is not marked changed');
      expect(editor.changesOf(0)['name'], original, '~changes does not hold the ORIGINAL value');
      expect(editor.changeCount, 1, 'wrong change count after one edit');

      // Back to the original: the cell is clean again, and so is the row.
      editor.setValue(0, 'name', original);
      expect(editor.stateOf(0), '', 'reverting by hand left the row modified');
      expect(editor.isChanged(0, 'name'), false, 'reverting by hand left the cell marked');
      expect(editor.isDirty, false, 'the editor stayed dirty after a manual revert');

      // Re-edit, then revert through the editor.
      editor.setValue(0, 'name', 'Edited twice');
      expect(editor.changesOf(0)['name'], original, 'the second edit re-captured the wrong original');
      editor.revertCell(0, 'name');
      expect(df.get('name', 0), original, 'revertCell did not restore the value');
      expect(editor.stateOf(0), '', 'revertCell left the row modified');

      // An invalid value STAYS pending even at the original value, so its marker
      // survives (the platform's own in-grid semantics).
      editor.setValue(0, 'quantity', -5);
      expect(editor.errorOf(0, 'quantity')?.kind, 'error', 'a min violation was not marked');
      expect(editor.errorCount, 1, 'wrong blocking-error count');
      editor.revertCell(0, 'quantity');
      expect(editor.errorOf(0, 'quantity'), null, 'revertCell left the error marker');

      // New row.
      const added = editor.addRow({sku: `${prefix}-added`, name: 'Added'});
      expect(editor.stateOf(added), 'new', 'addRow did not mark the row new');
      expect(editor.isChanged(added, 'name'), true, 'every cell of a new row counts as changed');

      // Deleted rows stay in the frame and leave the filter.
      const rowsBefore = df.rowCount;
      editor.markDeleted(0);
      expect(df.rowCount, rowsBefore, 'markDeleted removed the row from the frame');
      expect(editor.stateOf(0), 'deleted', 'markDeleted did not set the state');
      expect(df.filter.get(0), false, 'a deleted row is still in the filter');
      editor.unmarkDeleted(0);
      expect(editor.stateOf(0), '', 'unmarkDeleted did not restore the state');
      expect(df.filter.get(0), true, 'unmarkDeleted did not bring the row back');

      // discard() takes everything back to the loaded state.
      editor.discard();
      expect(editor.isDirty, false, 'discard left the editor dirty');
      expect(df.rowCount, rowsBefore - 1, 'discard did not remove the new row');
      expect(df.get('name', 0), original, 'discard did not restore the cell');
      editor.detach();
    } finally {
      await cleanup(prefix);
    }
  });

  test('transaction ops: changed columns only, expectedVersion, insert/delete shapes', async () => {
    const prefix = `fe-ops-${stamp()}`;
    await seed(prefix, 3);
    try {
      const editor = await editorFor(prefix);
      const df = editor.dataFrame;
      const version = df.get('version', 0);

      editor.setValue(0, 'name', 'Only this column');
      editor.addRow({sku: `${prefix}-ins`, name: 'Inserted', quantity: 7});
      editor.markDeleted(2);

      const ops = editor.buildOps().map((p) => p.op);
      expect(ops.length, 3, 'wrong number of ops');

      const update = ops.find((o) => o.op === 'update')!;
      expect(update.table, 'item', 'the op does not name the table');
      expect(Object.keys(update.values!).join(','), 'name',
        'the update carries columns that did not change');
      expect(update.expectedVersion, version, 'the update lost its optimistic-concurrency guard');
      expect(update.id, df.get('id', 0), 'the update addresses the wrong row');

      const insert = ops.find((o) => o.op === 'insert')!;
      expect((insert.values as any).sku, `${prefix}-ins`, 'the insert lost its values');
      expect('id' in (insert.values as any), false, 'a new row must not send an id');

      const del = ops.find((o) => o.op === 'delete')!;
      expect(del.id, df.get('id', 2), 'the delete addresses the wrong row');

      // A row that was added and then deleted never reached the server: no op.
      const ghost = editor.addRow({sku: `${prefix}-ghost`, name: 'Ghost', quantity: 0});
      editor.markDeleted(ghost);
      expect(editor.buildOps().length, 3, 'an add-then-delete row produced an op');

      // Undeleting it puts it back where it was: a 'new' row (it has no id —
      // buildOps' own predicate — and nothing but the frame knows about it).
      editor.unmarkDeleted(ghost);
      expect(editor.stateOf(ghost), 'new', 'undelete did not restore the new row');
      expect(editor.buildOps().length, 4, 'the restored row is invisible to save');
      expect(await editor.save(), true, 'the batch did not save');
      expect(await items().count({property: 'sku', operator: '=',
        value: `${prefix}-ghost`} as any), 1, 'the restored row never reached the server');
      editor.detach();
    } finally {
      await cleanup(prefix);
    }
  });

  test('save: one transaction (shared tx_id), version bumps, server readback', async () => {
    const prefix = `fe-save-${stamp()}`;
    const ids = await seed(prefix, 3);
    try {
      const editor = await editorFor(prefix);
      const df = editor.dataFrame;
      const versionBefore = df.get('version', 0);

      editor.setValue(0, 'name', 'Saved name');
      editor.setValue(1, 'quantity', 42);
      editor.addRow({sku: `${prefix}-added`, name: 'Added row', quantity: 5});
      editor.markDeleted(2);

      expect(await editor.save(), true, 'save() failed');
      expect(editor.isDirty, false, 'the editor stayed dirty after a successful save');

      // In-frame: versions bumped, the new row got its id, the deleted row is gone.
      expect(df.get('version', 0) > versionBefore, true, 'the version was not bumped in the frame');
      expect(editor.stateOf(0), '', 'the row state was not cleared');
      expect(df.rowCount, 3, 'the deleted row was not removed after the save');
      const addedId = df.get('id', 2);
      expect(addedId != null && `${addedId}` !== '', true, 'the inserted row did not get its id');

      // Server readback: value for value.
      const server = await items().query({filter: {property: 'sku', operator: 'like',
        value: `${prefix}%`} as any, sort: 'sku'});
      expect(server.length, 3, 'the server does not show the expected rows');
      expect(server.find((r: any) => r.id === ids[0])!.name, 'Saved name', 'the edit did not land');
      expect(server.find((r: any) => r.id === ids[1])!.quantity, 42, 'the second edit did not land');
      expect(server.some((r: any) => r.id === ids[2]), false, 'the deleted row is still there');
      expect(server.some((r: any) => r.sku === `${prefix}-added`), true, 'the inserted row is missing');

      // ONE transaction: every write of the batch shares a tx_id.
      const log = await items().auditLog({limit: 100});
      const mine = log.filter((e) => ids.includes(`${e.row_id}`) || `${e.row_id}` === `${addedId}`);
      expect(mine.length >= 4, true, `expected the seed + batch audit rows, got ${mine.length}`);
      // The seed insert is its own transaction; the batch is exactly one more.
      // Assert on the NON-NULL ids and demand that every batch row carries one —
      // a set of all-null tx_ids would also have size 1 and prove nothing.
      const batch = mine.filter((e) => e.op !== 'insert' || `${e.row_id}` === `${addedId}`);
      const batchTx = new Set(batch.filter((e) => e.tx_id != null).map((e) => `${e.tx_id}`));
      expect(batch.every((e) => e.tx_id != null), true, 'a batch audit row carries no tx_id');
      expect(batchTx.size, 1, `the batch spanned ${batchTx.size} transactions, not one`);
      editor.detach();
    } finally {
      await cleanup(prefix);
    }
  });

  test('save: a server validation error lands on the offending cell', async () => {
    const prefix = `fe-invalid-${stamp()}`;
    await seed(prefix, 2);
    try {
      const editor = await editorFor(prefix);
      // Duplicate the other row's sku: unique server-side, and the CLIENT has no
      // uniqueness rule — so this exercises the server → ~errors mapping, not the
      // inline validator.
      const duplicate = editor.dataFrame.get('sku', 1);
      editor.setValue(0, 'sku', duplicate);
      expect(editor.errorCount, 0, 'the client rejected the value before the server saw it');

      expect(await editor.save(), false, 'a duplicate sku was accepted');
      expect(editor.errorOf(0, 'sku')?.kind, 'error', 'the server error did not reach the cell');
      expect(editor.isDirty, true, 'the rejected edit was dropped instead of staying pending');
      // Nothing landed: the transaction rolled back.
      const server = await items().first({filter: {property: 'sku', operator: '=', value: duplicate} as any});
      expect(server != null, true, 'the readback row disappeared');
      editor.detach();
    } finally {
      await cleanup(prefix);
    }
  });

  test('version conflict: reload drops that row edits, overwrite lands the value', async () => {
    const prefix = `fe-409-${stamp()}`;
    const ids = await seed(prefix, 1);
    const original = (DG.DomainObjectHandler as any).showConflictDialog;
    try {
      const id = ids[0];
      let decision: 'reload' | 'overwrite' | null = 'reload';
      let asked = 0;
      // The dialog itself is covered by the WO-5 opener tests; here the OUTCOME
      // is what matters, so the standard dialog is answered programmatically.
      (DG.DomainObjectHandler as any).showConflictDialog = async () => {
        asked++;
        return decision;
      };

      // RELOAD: an out-of-band write makes the pending edit stale.
      let editor = await editorFor(prefix);
      // The grid snapshots the row when it becomes current, BEFORE the edit — so
      // the reload has a pre-reload snapshot to invalidate.
      editor.beginEdit(0);
      editor.setValue(0, 'name', 'My name');
      await items().update(id, {name: 'Their name'});
      expect(await editor.save(), true, 'the reload path did not finish the save');
      expect(asked, 1, 'the conflict dialog was not consulted');
      expect(editor.isDirty, false, 'the reload path left edits pending');
      expect((await items().get(id)).name, 'Their name', 'reload did not keep the server value');
      expect(editor.dataFrame.get('name', 0), 'Their name', 'the frame kept the discarded edit');

      // The pre-reload snapshot is gone: the NEXT in-grid edit records the value
      // the cell actually held after the reload, so a revert restores THAT.
      editor.beginEdit(0);
      editor.dataFrame.set('name', 0, 'Typed after the reload');
      editor.commitEdit(0, 'name');
      expect(editor.changesOf(0)['name'], 'Their name',
        'the edit after a reload captured a pre-reload original');
      editor.revertCell(0, 'name');
      expect(editor.dataFrame.get('name', 0), 'Their name', 'revertCell restored a stale value');
      editor.detach();

      // OVERWRITE: the same setup, the other answer.
      decision = 'overwrite';
      asked = 0;
      editor = await editorFor(prefix);
      editor.setValue(0, 'name', 'Mine wins');
      await items().update(id, {name: 'Theirs again'});
      expect(await editor.save(), true, 'the overwrite path did not finish the save');
      expect(asked, 1, 'the conflict dialog was not consulted on overwrite');
      expect((await items().get(id)).name, 'Mine wins', 'overwrite did not land the value');
      editor.detach();

      // DISMISSED: nothing is written and the row says why.
      decision = null;
      editor = await editorFor(prefix);
      editor.setValue(0, 'name', 'Never saved');
      await items().update(id, {name: 'Server again'});
      expect(await editor.save(), false, 'a dismissed conflict reported success');
      expect(editor.errorOf(0, 'name')?.kind, 'conflict', 'the dismissed conflict left no marker');
      expect(editor.isDirty, true, 'the dismissed conflict dropped the edit');
      expect((await items().get(id)).name, 'Server again', 'a dismissed conflict wrote anyway');
      editor.detach();
    } finally {
      (DG.DomainObjectHandler as any).showConflictDialog = original;
      await cleanup(prefix);
    }
  });

  test('refresh rebuilds and discards pending edits by design', async () => {
    const prefix = `fe-refresh-${stamp()}`;
    await seed(prefix, 2);
    try {
      const editor = await editorFor(prefix);
      const before = editor.dataFrame;
      editor.setValue(0, 'name', 'Doomed edit');
      editor.addRow({sku: `${prefix}-doomed`});
      editor.markDeleted(1);
      // Readable as dirty IMMEDIATELY before the call — this is the affordance
      // the contract tells callers to gate their refresh policy on.
      expect(editor.isDirty, true, 'the editor was not dirty before refresh()');

      const after = await editor.refresh();
      expect(after === before, false, 'refresh() did not rebuild the frame');
      expect(editor.dataFrame === after, true, 'the editor kept the stale frame');
      expect(editor.isDirty, false, 'refresh() carried the dirty state over');
      expect(editor.changeCount, 0, 'refresh() carried changes over');
      expect(after.rowCount, 2, 'the rebuilt frame does not match a fresh query');
      expect(after.get('name', 0), 'Item 0', 'the discarded edit survived the rebuild');
      for (let i = 0; i < after.rowCount; i++) {
        expect(editor.stateOf(i), '', `row ${i} kept service-column residue`);
        expect(after.getCol(CHANGES_COLUMN).get(i) ?? '', '', `row ${i} kept ~changes residue`);
        expect(after.getCol(ERRORS_COLUMN).get(i) ?? '', '', `row ${i} kept ~errors residue`);
      }
      // The server never saw any of it.
      expect((await items().count({property: 'sku', operator: 'like', value: `${prefix}%`} as any)), 2,
        'refresh() wrote something');
      editor.detach();
    } finally {
      await cleanup(prefix);
    }
  });

  test('onDirtyChanged: every transition, refresh() included', async () => {
    const prefix = `fe-dirty-${stamp()}`;
    await seed(prefix, 2);
    try {
      const editor = await editorFor(prefix);
      const seen: boolean[] = [];
      const sub = editor.onDirtyChanged.subscribe((dirty) => seen.push(dirty));
      try {
        editor.setValue(0, 'name', 'Edited');
        await editor.refresh();
        editor.setValue(0, 'name', 'Edited again');
        editor.discard();
        // The refresh transition is the one a save/discard prompt binds to: a
        // subscriber that saw `true` and never sees `false` prompts forever.
        expect(seen.join(','), 'true,false,true,false',
          `onDirtyChanged reported [${seen.join(', ')}]`);
        expect(editor.isDirty, false, 'the editor stayed dirty after discard()');
      } finally {
        sub.unsubscribe();
      }
      editor.detach();
    } finally {
      await cleanup(prefix);
    }
  });

  test('in-flight save: writes, discard and refresh are refused', async () => {
    const prefix = `fe-busy-${stamp()}`;
    await seed(prefix, 2);
    // Seam: pause the transaction so the whole save window is observable. The
    // data source is re-created on every `grok.dapi.domains` access, so the stub
    // goes on its prototype.
    const proto = Object.getPrototypeOf(grok.dapi.domains) as any;
    const realTransaction = proto.transaction;
    let release: () => void = () => {};
    const paused = new Promise<void>((resolve) => release = resolve);
    proto.transaction = async function(...args: any[]): Promise<any> {
      await paused;
      return await realTransaction.apply(this, args);
    };
    try {
      const editor = await editorFor(prefix);
      const grid = DomainGrid.forEditor(editor);
      const rowsBefore = editor.dataFrame.rowCount;
      editor.setValue(0, 'name', 'Saved under the lock');

      const saving = editor.save();
      expect(editor.isSaving, true, 'the editor did not close for the save');
      expect(grid.grid.props.allowEdit, false, 'the grid stayed editable during the save');

      // Every writer is refused — loudly, and without touching the batch the
      // in-flight transaction addresses.
      editor.setValue(1, 'name', 'Would be wiped');
      expect(editor.isChanged(1, 'name'), false, 'a write landed during the save');
      editor.addRow({sku: `${prefix}-nope`});
      expect(editor.dataFrame.rowCount, rowsBefore, 'addRow ran during the save');
      editor.discard();
      expect(editor.isDirty, true, 'discard() ran during the save');
      expect(await editor.refresh() === editor.dataFrame, true,
        'refresh() rebuilt the frame during the save');
      expect(editor.dataFrame.rowCount, rowsBefore, 'the frame moved under the save');

      release();
      expect(await saving, true, 'the save did not finish');
      expect(editor.isSaving, false, 'the editor stayed closed after the save');
      expect(grid.grid.props.allowEdit, grid.editable, 'the grid stayed locked after the save');
      expect(editor.isDirty, false, 'the save left the batch pending');
      const server = await items().query({filter: {property: 'sku', operator: 'like',
        value: `${prefix}%`} as any, sort: 'sku'});
      expect(server[0].name, 'Saved under the lock', 'the edit did not reach the server');
      grid.detach();
    } finally {
      proto.transaction = realTransaction;
      await cleanup(prefix);
    }
  });

  test('deleted rows stay hidden across a filter recompute', async () => {
    const prefix = `fe-filter-${stamp()}`;
    await seed(prefix, 3);
    try {
      const editor = await editorFor(prefix);
      const df = editor.dataFrame;
      editor.markDeleted(1);
      expect(df.filter.trueCount, 2, 'the deleted row was not excluded');

      // Another participant recomputes the filter (what a platform filter does
      // on every change): the cooperative AND must re-apply.
      df.filter.setAll(true, false);
      df.rows.requestFilter();
      expect(df.filter.get(1), false, 'a filter recompute un-hid the deleted row');
      expect(df.filter.trueCount, 2, 'the deleted row came back after a recompute');
      editor.detach();
    } finally {
      await cleanup(prefix);
    }
  });

  test('validateCellValue: server-parity messages', async () => {
    const properties = await grok.dapi.domains.registry.rowProperties('apitests.item');
    const quantity = properties.find((p) => p.name === 'quantity')!;
    const sku = properties.find((p) => p.name === 'sku')!;
    expect(validateCellValue(quantity, 3), null, 'a valid int was rejected');
    expect(validateCellValue(quantity, -1), 'Value should not be less than 0.0, passed: -1',
      'the min message does not match the server');
    expect(validateCellValue(quantity, 1.5), 'Types differ. Expected: int, passed: double',
      'a fractional int was not rejected');
    expect(validateCellValue(quantity, 'x'), 'Numerical value expected, passed: string',
      'a non-numeric int was not rejected');
    expect(validateCellValue(sku, ''), "Value can't be empty", 'a required empty value was accepted');
    expect(validateCellValue(sku, null), "Value can't be empty", 'a required null was accepted');
  });

  test('DomainGrid: platform decoration, hidden service columns, capability gating', async () => {
    const prefix = `fe-grid-${stamp()}`;
    await seed(prefix, 2);
    let grid: DomainGrid | null = null;
    try {
      grid = await DomainGrid.create(items() as any, {
        query: {filter: {property: 'sku', operator: 'like', value: `${prefix}%`} as any, sort: 'sku'},
      });
      const caps = grid.editor.capabilities;
      expect(grid.editable, caps.canEdit, 'the grid ignored the table capability');

      // The editing state is never visible or editable, whatever a handler does.
      for (const name of SERVICE_COLUMNS) {
        const gc = grid.grid.col(name);
        expect(gc != null, true, `${name} is missing from the grid`);
        expect(gc!.visible, false, `${name} is visible in the grid`);
        expect(gc!.editable, false, `${name} is editable in the grid`);
      }
      // Platform decoration ran: system columns hidden (the renderGrid contract).
      expect(grid.grid.col('id')?.visible, false, 'renderGrid did not hide the system columns');

      // Column security: only writable, non-reference columns take in-grid edits.
      if (caps.canEdit) {
        expect(grid.grid.props.allowEdit, true, 'an editable table produced a read-only grid');
        for (const p of grid.editor.properties) {
          const gc = grid.grid.col(p.name);
          if (gc != null && !caps.writableColumns.includes(p.name))
            expect(gc.editable, false, `${p.name} is editable without write access`);
        }
        expect(grid.grid.col('sku')?.editable, true, 'a writable column is not editable');
      }

      // The grid is the editor's mouth, never a second writer.
      grid.editor.setValue(0, 'name', 'Through the editor');
      expect(grid.editor.isChanged(0, 'name'), true, 'the edit was not tracked');
      expect(grid.dataFrame === grid.editor.dataFrame, true, 'the grid and editor drifted apart');
    } finally {
      grid?.detach();
      await cleanup(prefix);
    }
  });

  test('DomainGrid.decorate: an overriding handler of another type still wins', async () => {
    const seen: any[] = [];
    // A plugin handler that claims apitests.item rows through isApplicable while
    // declaring a type of its own — the shape the collapse rule must NOT drop.
    class ForeignHandler extends DG.ObjectHandler {
      retired = false;
      get type() { return this.retired ? 'apitests.foreign-retired' : 'apitests.foreign'; }
      isApplicable(x: any) {
        return !this.retired && x instanceof DG.DomainRow && x.typeName === 'apitests.item';
      }
      renderGrid(grid: _DG.Grid, options?: {items?: _DG.DataFrame}) {
        seen.push(options?.items ?? null);
        grid.columns.byName('quantity')!.visible = false;
      }
    }
    const handler = new ForeignHandler();
    DG.ObjectHandler.register(handler);
    const prefix = `fe-sentinel-${stamp()}`;
    await seed(prefix, 1);
    try {
      const df = await items().queryDf({filter: {property: 'sku', operator: 'like',
        value: `${prefix}%`} as any});
      const grid = DG.Grid.create(df);
      DomainGrid.decorate(grid, 'apitests.item', df);
      expect(seen.length, 1, 'the overriding handler was collapsed away');
      expect(grid.columns.byName('quantity')!.visible, false, 'its decoration did not land');

      // Retired: nothing claims the table, and the platform decoration is back
      // (the id column is a system column renderGrid hides).
      handler.retired = true;
      const plain = DG.Grid.create(df);
      DomainGrid.decorate(plain, 'apitests.item', df);
      expect(seen.length, 1, 'a retired handler still decorated');
      expect(plain.col('id')?.visible, false, 'the platform decoration did not run');
    } finally {
      handler.retired = true;
      await cleanup(prefix);
    }
  });

  test('DomainGrid: read-only degradation under a restricted user', async () => {
    const prefix = `fe-ro-${stamp()}`;
    await seed(prefix, 1);
    try {
      const outcome = await withRestrictedUser('wo7ro', async (probe) => {
        const table = items();
        // Give the throwaway user View only: rows are readable, nothing is writable.
        await table.grant(probe.group, 'View');
        try {
          await probe.asUser(async () => {
            grok.dapi.domains.invalidateUiCaches();
            const caps = await grok.dapi.domains.table('apitests.item').capabilities();
            expect(caps.canEdit, false, 'a View-only user reports canEdit');
            const grid = await DomainGrid.create(grok.dapi.domains.table('apitests.item') as any,
              {query: {filter: {property: 'sku', operator: 'like', value: `${prefix}%`} as any}});
            try {
              expect(grid.editable, false, 'the grid stayed editable for a View-only user');
              expect(grid.grid.props.allowEdit, false, 'allowEdit survived read-only degradation');
              for (const p of grid.editor.properties) {
                const gc = grid.grid.col(p.name);
                if (gc != null)
                  expect(gc.editable, false, `${p.name} is editable for a View-only user`);
              }
            } finally {
              grid.detach();
            }
          });
        } finally {
          // The grant is ours whatever the body did — never leave the probe
          // user holding View on the fixture table.
          await table.revoke(probe.group, 'View');
        }
        return 'checked';
      });
      // withRestrictedUser resolves null when it could not build a restricted
      // session — say so instead of reporting a green run that checked nothing.
      expect(outcome, 'checked',
        'read-only degradation was NOT verified: no restricted session (see the console for why)');
    } finally {
      grok.dapi.domains.invalidateUiCaches();
      await cleanup(prefix);
    }
  });
});
