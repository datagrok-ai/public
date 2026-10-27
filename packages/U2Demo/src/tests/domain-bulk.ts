/* Bulk edit against a real domain table: the ⋯ menu's items, the dialog's field picker,
   "N selected" writing `id in (…)` and "all M matching" writing the list's filter —
   all of it through the server's `POST …/{table}/update`. Fixture: `apitests.item`, rows
   prefixed per run and deleted afterwards. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {Control, Tooltip} from '@datagrok-libraries/u2';
import type {DomainSource} from '@datagrok-libraries/ems';
import {domains, DomainApp, DomainTable} from '@datagrok-libraries/ems/src/dg';

const ITEMS = 'apitests.item';
const BASE = '/apps/U2Demo/bulk';

category('U2: domain bulk', () => {
  const items = () => grok.dapi.domains.table(ITEMS);
  const prefix = `U2B-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
  const query = `sku starts "${prefix}"`;
  let table: DomainTable;
  let view: DG.View;
  let app: DomainApp;
  let ids: string[];

  async function ready(src: DomainSource): Promise<void> {
    await waitUntil(() => src.state.value !== 'idle' && src.state.value !== 'loading',
      () => `source is ${src.state.value}: ${String(src.error.value)}`);
    if (src.state.value !== 'ready')
      throw new Error(`source is ${src.state.value}: ${String(src.error.value)}`);
  }

  async function reloaded(): Promise<void> {
    await ready(app.listSource);
  }

  type Field = {value: {value: unknown}, root: HTMLElement};
  /** A control in the dialog by name: the field picker (`fields`), "Apply to" (`target`), a column's editor. */
  const control = (name: string): Field | undefined => {
    const el = document.querySelector<HTMLElement>(`.u2-dialog [data-u2-name="${name}"]`);
    return el === null ? undefined : Control.forElement(el) as unknown as Field;
  };
  const shown = (name: string) =>
    document.querySelector(`.u2-dialog [data-u2-name="${name}"]`)?.closest('[hidden]') === null;
  const buttonNamed = (text: string) => [...document.querySelectorAll<HTMLButtonElement>('.u2-dialog button')]
    .find((b) => b.textContent === text)!;
  /** The primary button, "Update N …". */
  const ok = () => [...document.querySelectorAll<HTMLButtonElement>('.u2-dialog [data-u2-part="footer"] button')]
    .pop()!;
  /** A refusal, said inside the dialog. */
  const refusal = () => document.querySelector('.u2-dialog [data-u2-part="error-row"]')?.textContent ?? '';
  const quantities = async () => (await items().query({filter: query, sort: 'sku'}))
    .map((r) => r.quantity).join();

  /** What holds the button: its tooltip. */
  async function okReason(): Promise<string> {
    ok().dispatchEvent(new Event('pointerenter'));
    await waitUntil(() => Tooltip.root.style.display === 'block' && Tooltip.root.textContent !== '',
      'the button\'s tooltip shows');
    const text = Tooltip.root.textContent ?? '';
    ok().dispatchEvent(new Event('pointerleave'));
    return text;
  }

  /** The dialog up and rendered. The promise is handed back WRAPPED: awaiting an async function
   * that returns it would chain on it and only come back once the dialog is answered. */
  async function open(): Promise<{running: Promise<number | null>}> {
    const running = domains.bulkEdit(app.listSource);
    await waitUntil(() => control('fields') !== undefined, 'the bulk dialog is up', 5000);
    expect(control('fields') !== undefined, true, 'the bulk dialog is up');
    return {running};
  }

  /** Picks the fields to change, then fills each of them. */
  async function write(values: Record<string, unknown>): Promise<void> {
    control('fields')!.value.value = Object.keys(values);
    await waitUntil(() => Object.keys(values).every(shown), 'the chosen fields show');
    for (const [name, value] of Object.entries(values))
      control(name)!.value.value = value;
    await waitUntil(() => !ok().disabled, 'the button is on');
  }

  async function selectRows(keep: (id: string) => boolean): Promise<void> {
    const df = app.listSource.df.value as unknown as DG.DataFrame;
    const idColumn = df.columns.byName('id')!;
    let kept = 0;
    for (let i = 0; i < df.rowCount; i++) {
      const on = keep(idColumn.get(i));
      df.selection.set(i, on);
      kept += on ? 1 : 0;
    }
    const selected = () => app.listSource.selection.value;
    await waitUntil(() => selected().length === kept && selected().every((r) => keep(r.id)),
      'the selection to reach the source');
  }

  before(async () => {
    const inserted = await items().insert([
      {sku: `${prefix}-1`, name: 'Alpha', quantity: 1, origin: 'fixture'},
      {sku: `${prefix}-2`, name: 'Beta', quantity: 2, origin: 'fixture'},
      {sku: `${prefix}-3`, name: 'Gamma', quantity: 3, origin: 'fixture'},
    ]);
    ids = inserted.map((x) => x.id);
    table = await domains.table(ITEMS);
    view = table.app({name: 'Bulk demo', path: BASE, query}) as DG.View;
    app = DomainApp.of(view)!;
    grok.shell.addView(view);
    view.temp['ignoreCloseAll'] = true;
    await ready(app.listSource);
  });

  after(async () => {
    app.session.discard();
    view.close();
    await items().deleteWhere(query);
  });

  test('the ⋯ menu offers Import…, Bulk edit… and Show deleted', async () => {
    expect(app.menuActions().map((a) => a.name).join(), 'Import…,Bulk edit…,Show deleted');
  });

  test('the field picker decides what is written; the button says why it is off, never a balloon', async () => {
    await selectRows(() => false);
    const {running} = await open();
    expect(shown('quantity'), false, 'a field not chosen is not in the form');
    expect(ok().disabled, true, 'nothing chosen yet');
    expect(await okReason(), 'Choose the fields to change');
    expect((control('target')?.root.textContent ?? '').includes('All 3'), true,
      'no target the server would refuse: nothing is selected, so the filter is the one target');
    control('fields')!.value.value = ['quantity'];
    await waitUntil(() => shown('quantity'), 'the chosen field shows');
    expect(shown('quantity'), true, 'a chosen field shows its editor');
    expect(control('quantity')!.root.querySelector('.u2-input-error')?.textContent ?? '', '',
      'nothing here is required: a chosen field left empty CLEARS the column');
    expect(ok().disabled, false);
    buttonNamed('CANCEL').click();
    expect(await running, null);
  });

  test('the chosen column is written into the selection only', async () => {
    await selectRows((id) => id !== ids[1]);
    expect(app.listSource.selection.value.length, 2, 'two of the three rows are selected');
    const {running} = await open();
    await write({quantity: 9});
    ok().click();
    expect(await running, 2);
    await reloaded();
    expect(await quantities(), '9,2,9', 'the unselected row kept its quantity');
  });

  test('"all matching" writes everything the list\'s filter matches', async () => {
    await selectRows(() => false);
    const {running} = await open();
    await write({quantity: 4});
    ok().click();
    expect(await running, 3);
    await reloaded();
    expect(await quantities(), '4,4,4');
  });

  test('the live probe is ONE aggregate: count + max(updated_on), narrowed by search and deleted', async () => {
    const handle = table.table as unknown as
      {probe(spec: unknown): Promise<{count: number, last: string | null}>};
    const all = await handle.probe({filter: query});
    expect(all.count, 3);
    expect(typeof all.last, 'string', 'the aggregate compiler takes the system column updated_on as a measure');
    expect((await handle.probe({filter: query, search: 'Alpha'})).count, 1, 'search reaches the aggregate');
    const trashed = await handle.probe({filter: query, deleted: 'only'});
    expect(trashed.count, 0, 'deleted reaches the aggregate');
    expect(trashed.last, null, 'and an empty match dates to null');
  });

  test('an immutable column is refused by the server in the dialog, and nothing is written', async () => {
    const {running} = await open();
    await write({origin: 'rewritten', quantity: 7});
    ok().click();
    await waitUntil(() => refusal() !== '', 'the refusal is the dialog\'s error line');
    expect(shown('origin'), true, 'the dialog is still up, the fields still chosen');
    expect(await quantities(), '4,4,4', 'the whole call rolled back');
    buttonNamed('CANCEL').click();
    expect(await running, null);
  });
});
