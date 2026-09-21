/* The write UI over a WRITABLE external binding (EMS external bindings, phase C WO-C5): the same
   `domains.table(...).app()` over `extwlive.thing` — a new row typed with its key, an edit
   guarded by the changed columns' old values, the conflict dialog when they moved behind the app
   (Overwrite and Reload), a delete. Fixture: the `extwlive` scratch package
   (core/docs/features/ems/external-bindings/fixtures/extwlive), which the tester publishes before
   a run; the category skips itself while it is not registered. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {Control} from '@datagrok-libraries/u2';
import type {DomainSource} from '@datagrok-libraries/u2';
import {domains, DomainApp, DomainTable} from '@datagrok-libraries/u2/src/dg/index.js';

const SCHEMA = 'extwlive';
const THINGS = `${SCHEMA}.thing`;
const BASE = `/domains/${SCHEMA}/thing`;

category('U2: domain external write', () => {
  const things = () => grok.dapi.domains.table(THINGS);
  const tid = 720000 + (Date.now() % 90000);
  const id = `${tid}`;
  let skip: string | null = null;
  let seeded = false;
  let table: DomainTable;
  let app: DomainApp;

  async function ready(src: DomainSource): Promise<void> {
    for (let i = 0; i < 200 && (src.state.value === 'idle' || src.state.value === 'loading'); i++)
      await delay(50);
    if (src.state.value !== 'ready')
      throw new Error(`source is ${src.state.value}: ${String(src.error.value)}`);
  }

  async function until(check: () => boolean, what: string): Promise<void> {
    for (let i = 0; i < 100 && !check(); i++)
      await delay(50);
    if (!check())
      throw new Error(`timed out waiting for ${what}`);
  }

  const skipped = (): boolean => {
    if (skip != null)
      console.log(`skipped: ${skip}`);
    return skip != null;
  };

  const conflictDialog = (): HTMLElement | undefined =>
    [...document.querySelectorAll<HTMLElement>('.d4-dialog')]
      .find((d) => d.textContent?.includes('changed since you read it'));
  const dialogButton = (dialog: HTMLElement, text: string): HTMLButtonElement | undefined =>
    [...dialog.querySelectorAll('button')].find((b) => b.textContent?.trim() === text);

  /** Opens the row, moves `n` behind the app, edits it here, saves, answers the dialog with [button]. */
  async function conflict(behind: number, mine: number, button: string): Promise<boolean> {
    expect(await app.open(`${BASE}/${id}`), true);
    const source = app.entitySource.value!;
    await ready(source);
    await until(() => app.form.value?.input('n') !== undefined, 'the form');
    await things().update(id, {n: behind});
    app.form.value!.input('n')!.value.value = mine;
    const saving = app.session.save();
    await until(() => conflictDialog() !== undefined, 'the conflict dialog');
    const dialog = conflictDialog()!;
    expect(dialog.textContent?.includes(`${behind}`), true, 'the dialog shows what the row holds now');
    dialogButton(dialog, button)!.click();
    return await saving;
  }

  // The row the edit, conflict and delete tests work on; the New test inserts its own.
  before(async () => {
    const schemas = await grok.dapi.domains.schemas.list();
    if (!schemas.some((s) => s.name === SCHEMA)) {
      skip = `the ${SCHEMA} fixture is not registered`;
      return;
    }
    await things().insert({tid, note: `u2 write ${tid}`, n: 1});
    seeded = true;
    table = await domains.table(THINGS);
    app = domains.app({table, base: BASE, children: false});
    document.body.append(app.root);
    await ready(app.listSource);
  });

  after(async () => {
    app?.session.discard();
    app?.dispose();
    app?.root.remove();
    if (!seeded)
      return;
    try {
      await things().delete(id);
    } catch (e) {
      console.warn(`${THINGS} cleanup: ${e}`);
    }
  });

  test('New: the form takes the key; Save inserts the row under it', async () => {
    if (skipped())
      return;
    const newTid = tid + 3;
    try {
      expect(table.table.support.concurrency, 'expected');
      expect(await app.goTo('entity', DomainApp.NEW), true);
      await until(() => app.form.value?.input('tid') !== undefined, 'the draft form');
      app.form.value!.input('tid')!.value.value = newTid;
      app.form.value!.input('note')!.value.value = `u2 new ${newTid}`;
      expect(await app.session.save(), true);
      const row = await things().get(`${newTid}`);
      expect(row.note, `u2 new ${newTid}`);
    } finally {
      try {
        await things().delete(`${newTid}`);
      } catch (e) {
        if (!(e instanceof DG.DomainNotFoundError))
          throw e;
      }
    }
  });

  test('an edit saves guarded, and the row reads back', async () => {
    if (skipped())
      return;
    expect(await app.open(`${BASE}/${id}`), true);
    const source = app.entitySource.value!;
    await ready(source);
    await until(() => app.form.value?.input('n') !== undefined, 'the form');
    app.form.value!.input('n')!.value.value = 3;
    expect(await app.session.save(), true);
    expect((await things().get(id)).n, 3);
    expect(app.session.isDirty.value, false);
  });

  test('a change behind the app: the conflict dialog; Overwrite saves mine', async () => {
    if (skipped())
      return;
    expect(await conflict(99, 5, 'OVERWRITE'), true, 'the overwrite path did not finish the save');
    expect(conflictDialog() === undefined, true, 'the dialog closed');
    expect((await things().get(id)).n, 5, 'overwrite did not land the value');
    expect(app.form.value!.input('n')!.value.value, 5);
  });

  test('a change behind the app: Reload shows the server value and drops mine', async () => {
    if (skipped())
      return;
    expect(await conflict(99, 7, 'RELOAD'), true, 'the reload path did not finish the save');
    expect((await things().get(id)).n, 99, 'reload wrote anyway');
    await until(() => app.form.value?.input('n')?.value.value === 99, 'the reloaded value');
    expect(app.session.isDirty.value, false);
  });

  test('Delete removes the row', async () => {
    if (skipped())
      return;
    expect(await app.goTo('list'), true);
    await ready(app.listSource);
    app.listSource.query.value = `tid = ${tid}`;
    await delay(100);
    await ready(app.listSource);
    const row = app.listSource.rows.byKey(id)!;
    expect(row !== undefined, true, 'the row is listed');
    app.list.actionsFor(row).find((a) => a.name === 'Delete')!.run();
    expect(app.session.isDirty.value, true);
    expect(await app.session.save(), true);
    // get(id) resolves null for a gone row (the client's contract), it does not throw
    expect(await things().get(id), null, 'the row is still there');
  });

  test('Import: the wizard offers what the binding declares; an upsert lands and reports merged rows', async () => {
    if (skipped())
      return;
    const ids = [tid + 1, tid + 2].map(String);
    const source = DG.DataFrame.fromCsv(`tid,note\n${ids[0]},u2 import a\n${ids[1]},u2 import b`);
    const buttonNamed = (text: string) => [...document.querySelectorAll<HTMLButtonElement>('.u2-dialog button')]
      .find((b) => b.textContent === text)!;
    const named = (name: string) => document.querySelector(`[data-u2-name="${name}"]`);
    const report = () => document.querySelector('.u2-domain-import-report')?.textContent ?? '';
    let running: Promise<unknown> | undefined;
    try {
      running = domains.import(table, {source});
      await until(() => named('mode') !== null, 'the wizard');
      expect(named('allOrNothing') === null, true, 'no partial import on an external binding');
      expect(named('errorOnDuplicate') === null, true, 'a duplicate is always an error there');
      (Control.forElement(named('mode')!) as unknown as {value: {value: string}}).value.value = 'upsert';
      await delay(100);
      buttonNamed('NEXT').click();
      await delay(200);
      buttonNamed('NEXT').click();
      await delay(200);
      const preview = document.querySelector('.u2-domain-import-preview')?.textContent ?? '';
      expect(preview.includes('Rows are checked when imported'), true, preview);
      buttonNamed('NEXT').click();
      await until(() => /merged|failed|aborted/.test(report()), 'the report');
      expect(report().includes('0 inserted, 2 merged'), true, report());
      buttonNamed('CLOSE').click();
      const result = await running as {merged?: number, rows: {status: string}[]};
      expect(result.merged, 2);
      expect(result.rows.every((r) => r.status === 'merged'), true);
      for (const [i, rowId] of ids.entries())
        expect((await things().get(rowId)).note, `u2 import ${'ab'[i]}`);
    } finally {
      // a failed expect must not leave the wizard modal over the next category
      if (document.querySelector('.u2-domain-import') !== null)
        [buttonNamed('CANCEL'), buttonNamed('CLOSE')]
          .find((b) => b !== undefined && b.style.display !== 'none')?.click();
      await running;
      for (const rowId of ids) {
        try {
          await things().delete(rowId);
        } catch (e) {
          if (!(e instanceof DG.DomainNotFoundError))
            throw e;
        }
      }
    }
  });
});
