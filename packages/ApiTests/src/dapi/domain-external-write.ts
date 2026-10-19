import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {thrown} from './domain-lifecycle';

// The typed client and DG.DomainFrameEditor over a WRITABLE external binding: the `extwlive`
// fixture (core/docs/features/ems/external-bindings/fixtures/extwlive), which the tester
// publishes before a run. Rows are addressed by their business key (`tid`), and an edit is
// guarded by `expected` — the columns' values as last read — never by `version`. The category
// skips cleanly while the fixture is not registered (a development convenience only: the
// release evidence runs it positively).
category('Dapi: domain external write', () => {
  const DomainFrameEditor = DG.DomainFrameEditor;
  const things = () => grok.dapi.domains.table('extwlive.thing');
  let skip: string | null = null;

  before(async () => {
    const schemas = await grok.dapi.domains.schemas.list();
    if (!schemas.some((s) => s.name === 'extwlive'))
      skip = 'the extwlive fixture is not registered';
  });

  const skipped = (): boolean => {
    if (skip != null)
      console.log(`skipped: ${skip}`);
    return skip != null;
  };

  const isA = (err: any, ctor: any) =>
    expect(err instanceof ctor, true, `expected ${ctor.name}, got ${err?.constructor?.name}: ${err?.message}`);

  test('update with expected: a matching guard writes, a stale one names what moved', async () => {
    if (skipped())
      return;
    const tid = 900000 + (Date.now() % 90000);
    const [ins] = await things().insert({tid, note: `guard ${tid}`, n: 1});
    try {
      expect(ins.id, `${tid}`, 'an external row is addressed by the key the caller supplied');
      await things().update(ins.id, {n: 2}, {expected: {n: 1}});
      expect((await things().get(ins.id)).n, 2);

      let err = await thrown(() => things().update(ins.id, {n: 3}, {expected: {n: 1}}));
      isA(err, DG.DomainVersionConflictError);
      expect(err.code, 'version-conflict');
      expect(err.id, ins.id);
      expect(err.body.expected.n, 1, JSON.stringify(err.body));
      expect(err.body.current.n, 2, JSON.stringify(err.body));

      err = await thrown(() => things().update(ins.id, {n: 4}, {expected: {price: 1.5}}));
      isA(err, DG.DomainUnsupportedError);
      expect(err.op, 'expected', 'a float cannot guard');

      err = await thrown(() => things().update(ins.id, {n: 4}, {version: 1}));
      isA(err, DG.DomainValidationError);
      expect(err.message.includes('applies to platform-stored tables'), true, err.message);
      expect((await things().get(ins.id)).n, 2, 'a refused guard writes nothing');

      await things().delete(ins.id);
      err = await thrown(() => things().update(ins.id, {n: 5}, {expected: {n: 2}}));
      isA(err, DG.DomainNotFoundError);
    } finally {
      try {
        await things().delete(ins.id);
      } catch (e) {
        if (!(e instanceof DG.DomainNotFoundError))
          throw e;
      }
    }
  });

  test('transaction: an update op with expected writes, a stale one rolls the batch back naming the op', async () => {
    if (skipped())
      return;
    const tid = 800000 + (Date.now() % 90000);
    const [ins] = await things().insert({tid, note: `tx guard ${tid}`, n: 1});
    try {
      const results = await grok.dapi.domains.transaction('extwlive',
        [{op: 'update', table: 'thing', id: ins.id, values: {n: 2}, expected: {n: 1}}]);
      expect(results[0].id, ins.id);
      expect((await things().get(ins.id)).n, 2);

      const err = await thrown(() => grok.dapi.domains.transaction('extwlive', [
        {op: 'update', table: 'thing', id: ins.id, values: {descr: 'first'}},
        {op: 'update', table: 'thing', id: ins.id, values: {n: 3}, expected: {n: 1}},
      ]));
      isA(err, DG.DomainVersionConflictError);
      expect(err.opIndex, 1, JSON.stringify(err.body));
      expect(err.expected.n, 1);
      expect(err.current.n, 2);
      const row = await things().get(ins.id);
      expect(row.n, 2, 'the stale op wrote nothing');
      // the read side collapses a NULL text column to '' (F-1)
      expect(row.descr == null || row.descr === '', true, 'the op before it was rolled back');
    } finally {
      await things().delete(ins.id);
    }
  });

  test('batch upsert: rows keyed by the caller land as merged, and again on a change', async () => {
    if (skipped())
      return;
    const base = 600000 + (Date.now() % 90000);
    const rows = [{tid: base, note: `merge ${base}`, n: 1}, {tid: base + 1, note: `merge ${base + 1}`, n: 1}];
    try {
      let report = await things().batch(rows, {mode: 'upsert'});
      expect(report.merged, 2, JSON.stringify(report));
      expect(report.inserted, 0);
      expect(report.rows.every((r) => r.status === 'merged'), true, JSON.stringify(report.rows));

      report = await things().batch([{...rows[0], n: 2}], {mode: 'upsert'});
      expect(report.merged, 1, JSON.stringify(report));
      expect(report.rows[0].status, 'merged');
      expect((await things().get(`${base}`)).n, 2, 'the second upsert did not update');
    } finally {
      for (const row of rows) {
        try {
          await things().delete(`${row.tid}`);
        } catch (e) {
          if (!(e instanceof DG.DomainNotFoundError))
            throw e;
        }
      }
    }
  });

  const conflictDialog = (): HTMLElement | undefined => Array.from(document.querySelectorAll<HTMLElement>('.d4-dialog'))
    .find((d) => d.textContent?.includes('changed since you read it'));
  const dialogButton = (dialog: HTMLElement | undefined, text: string): HTMLButtonElement | undefined =>
    Array.from(dialog?.querySelectorAll('button') ?? []).find((b) => b.textContent?.trim() === text);

  /** Waits for the editor's `expected` conflict dialog and clicks its [button]. */
  async function answerConflict(button: string): Promise<void> {
    let dialog: HTMLElement | undefined;
    for (let i = 0; i < 100 && dialog == null; i++) {
      await DG.delay(50);
      dialog = conflictDialog();
    }
    expect(dialog != null, true, 'the conflict dialog did not open');
    const target = dialogButton(dialog, button);
    expect(target != null, true, `the conflict dialog has no ${button} button`);
    target!.click();
  }

  /** What a failed wait must not leave behind: an open conflict dialog, a save in flight, a
   * bound editor. */
  async function release(editor: _DG.DomainFrameEditor | undefined, saving: Promise<boolean> | undefined):
      Promise<void> {
    dialogButton(conflictDialog(), 'CANCEL')?.click();
    await saving?.catch(() => false);
    editor?.detach();
  }

  test('editor: buildOps guards the changed columns by declared type, with the $$ escape', async () => {
    if (skipped())
      return;
    const tid = 700000 + (Date.now() % 90000);
    const [ins] = await things().insert({tid, note: `$guard ${tid}`, n: 1, price: 1.5});
    let editor: _DG.DomainFrameEditor | undefined;
    try {
      editor = await DomainFrameEditor.create(things() as any,
        {query: {filter: {property: 'tid', operator: '=', value: tid}}});
      expect(editor.access.support.concurrency, 'expected', 'the fixture declares the old-value guard');
      editor.setValue(0, 'n', 2);
      let ops = editor.buildOps().map((p) => p.op);
      expect(ops.length, 1);
      expect(ops[0].op, 'update');
      expect(JSON.stringify(ops[0].values), '{"n":2}');
      expect(JSON.stringify(ops[0].expected), '{"n":1}', 'the original of the changed int guards');
      expect('expectedVersion' in ops[0], false, 'no version guard on an external row');
      editor.revertCell(0, 'n');

      editor.setValue(0, 'price', 2.5);
      ops = editor.buildOps().map((p) => p.op);
      expect('expected' in ops[0], false, 'a float original cannot guard');
      editor.revertCell(0, 'price');

      editor.setValue(0, 'note', `changed ${tid}`);
      ops = editor.buildOps().map((p) => p.op);
      expect(JSON.stringify(ops[0].expected), JSON.stringify({note: `$$guard ${tid}`}),
        'the guard is escaped');
      editor.revertCell(0, 'note');

      // A guarded save lands, and the post-save re-read goes per id: a column moved behind
      // the editor comes back with it.
      editor.setValue(0, 'n', 2);
      await things().update(ins.id, {descr: 'behind'});
      expect(await editor.save(), true);
      expect(editor.isDirty, false);
      expect((await things().get(ins.id)).n, 2);
      expect(editor.dataFrame.get('descr', 0), 'behind', 'the per-id re-read did not land');
    } finally {
      editor?.detach();
      await things().delete(ins.id);
    }
  });

  test('editor: a guard is a literal — a draft-looking original goes out as it is, a uuid-shaped one does not guard', async () => {
    if (skipped())
      return;
    const tid = 730000 + (Date.now() % 90000);
    const uuid = '0f1e2d3c-4b5a-4978-8f6e-5d4c3b2a1900';
    const [ins] = await things().insert({tid, note: '~new:plain', descr: uuid, n: 1});
    let editor: _DG.DomainFrameEditor | undefined;
    try {
      editor = await DomainFrameEditor.create(things() as any,
        {query: {filter: {property: 'tid', operator: '=', value: tid}}});
      editor.setValue(0, 'note', `plain ${tid}`);
      let ops = editor.buildOps().map((p) => p.op);
      expect(JSON.stringify(ops[0].expected), '{"note":"~new:plain"}', 'the original is not a reference');
      editor.revertCell(0, 'note');

      editor.setValue(0, 'descr', `text ${tid}`);
      ops = editor.buildOps().map((p) => p.op);
      expect('expected' in ops[0], false, 'a uuid-shaped original cannot guard');
    } finally {
      editor?.detach();
      await things().delete(ins.id);
    }
  });

  test('editor: a stale guard conflicts — Overwrite retries the row unguarded, Reload takes the row', async () => {
    if (skipped())
      return;
    const tid = 710000 + (Date.now() % 90000);
    const [ins] = await things().insert({tid, note: `conflict ${tid}`, n: 1});
    let editor: _DG.DomainFrameEditor | undefined;
    let saving: Promise<boolean> | undefined;
    try {
      editor = await DomainFrameEditor.create(things() as any,
        {query: {filter: {property: 'tid', operator: '=', value: tid}}});
      const conflicts: _DG.DomainVersionConflictError[] = [];
      editor.onConflict.subscribe((e) => conflicts.push(e));

      // OVERWRITE: the stale guard is refused, the row's changed values go out again unguarded.
      editor.setValue(0, 'n', 5);
      await things().update(ins.id, {n: 99});
      saving = editor.save();
      await answerConflict('OVERWRITE');
      expect(await saving, true, 'the overwrite path did not finish the save');
      expect(conflicts.length, 1, 'the conflict did not reach onConflict');
      expect(conflicts[0].expected!.n, 1, JSON.stringify(conflicts[0].body));
      expect(conflicts[0].current!.n, 99, JSON.stringify(conflicts[0].body));
      expect(conflicts[0].opIndex, 0);
      expect((await things().get(ins.id)).n, 5, 'overwrite did not land the value');
      expect(editor.dataFrame.get('n', 0), 5);
      // Scoped to that save: the next edit of the same row is guarded again.
      editor.setValue(0, 'n', 6);
      expect(JSON.stringify(editor.buildOps()[0].op.expected), '{"n":5}', 'the overwrite outlived its save');
      editor.revertCell(0, 'n');

      // RELOAD: the edit is dropped and the row shows what the server holds.
      editor.setValue(0, 'n', 7);
      await things().update(ins.id, {n: 99});
      saving = editor.save();
      await answerConflict('RELOAD');
      expect(await saving, true, 'the reload path did not finish the save');
      expect(conflicts.length, 2);
      expect(editor.isDirty, false, 'the reload path left edits pending');
      expect(editor.dataFrame.get('n', 0), 99, 'the frame kept the discarded edit');
      expect((await things().get(ins.id)).n, 99, 'reload wrote anyway');

      // CANCEL: nothing is written and the cell says why.
      editor.setValue(0, 'n', 8);
      await things().update(ins.id, {n: 100});
      saving = editor.save();
      await answerConflict('CANCEL');
      expect(await saving, false, 'a dismissed conflict reported success');
      expect(editor.errorOf(0, 'n')?.kind, 'conflict', 'the dismissed conflict left no marker');
      expect((await things().get(ins.id)).n, 100, 'a dismissed conflict wrote anyway');
    } finally {
      await release(editor, saving);
      await things().delete(ins.id);
    }
  });
});
