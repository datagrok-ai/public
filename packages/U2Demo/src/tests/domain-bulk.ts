/* Bulk edit against a real domain table (WO 3-8): the ⋯ menu's items, the dialog's include
   checkboxes, "N selected" writing `id in (…)` and "all M matching" writing the list's filter —
   all of it through the server's `POST …/{table}/update` (WO 3-1). Fixture: `apitests.item`, rows
   prefixed per run and deleted afterwards. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {Control} from '@datagrok-libraries/u2';
import type {DomainSource} from '@datagrok-libraries/u2';
import {domains, DomainApp, DomainTable} from '@datagrok-libraries/u2/src/dg/index.js';

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
    for (let i = 0; i < 200 && (src.state.value === 'idle' || src.state.value === 'loading'); i++)
      await delay(50);
    if (src.state.value !== 'ready')
      throw new Error(`source is ${src.state.value}: ${String(src.error.value)}`);
  }

  async function reloaded(): Promise<void> {
    await delay(100);
    await ready(app.listSource);
  }

  const includeBox = (name: string) =>
    document.querySelector<HTMLInputElement>(`[data-u2-include="${name}"]`);
  const fieldRoot = (name: string) => includeBox(name)!.parentElement!.querySelector('.u2-input-root');
  const buttonNamed = (text: string) => [...document.querySelectorAll<HTMLButtonElement>('.u2-dialog button')]
    .find((b) => b.textContent === text)!;
  const quantities = async () => (await items().query({filter: query, sort: 'sku'}))
    .map((r) => r.quantity).join();

  /** The dialog up and rendered. The promise is handed back WRAPPED: awaiting an async function
   * that returns it would chain on it and only come back once the dialog is answered. */
  async function open(): Promise<{running: Promise<number | null>}> {
    const running = domains.bulkEdit(app.listSource);
    for (let i = 0; i < 100 && includeBox('quantity') === null; i++)
      await delay(50);
    expect(includeBox('quantity') !== null, true, 'the bulk dialog is up');
    return {running};
  }

  const fieldOf = (name: string) =>
    Control.forElement(fieldRoot(name)) as unknown as
      {value: {value: unknown}, enabled: boolean, root: HTMLElement};

  /** A real click, as a user's is — the box must not be inside the disabled input it enables. */
  function write(name: string, value: unknown): void {
    includeBox(name)!.click();
    fieldOf(name).value.value = value;
  }

  async function selectRows(keep: (id: string) => boolean): Promise<void> {
    const df = app.listSource.df.value as unknown as DG.DataFrame;
    const idColumn = df.columns.byName('id')!;
    for (let i = 0; i < df.rowCount; i++)
      df.selection.set(i, keep(idColumn.get(i)));
    await delay(100);
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

  test('the ⋯ menu offers Import…, Bulk edit… and Trash', async () => {
    expect(app.menuActions().map((a) => a.name).join(), 'Import…,Bulk edit…,Trash');
  });

  test('a real click on the include box enables its editor; OK is gated inline, not by balloon', async () => {
    await selectRows(() => false);
    const {running} = await open();
    expect(fieldOf('quantity').enabled, false, 'an unchecked column has no editor');
    expect(includeBox('quantity')!.parentElement!.classList.contains('u2-domain-bulk-row'), true,
      'the box sits OUTSIDE the input it enables — inside a disabled one it would be dead');
    expect(buttonNamed('OK').disabled, true, 'nothing checked yet');
    expect((document.querySelector('.u2-domain-bulk-reason')?.textContent ?? '')
      .includes('Check the fields to write'), true, 'and the dialog says so, in the dialog');
    expect(fieldOf('quantity').root.querySelector('.u2-input-error')?.textContent ?? '', '',
      'nothing here is required: a checked field left empty CLEARS the column');
    includeBox('quantity')!.click();
    await delay(100);
    expect(fieldOf('quantity').enabled, true, 'the click reached it');
    expect(buttonNamed('OK').disabled, false);
    const options = [...document.querySelectorAll<HTMLOptionElement>('[data-u2-name="target"] option')]
      .map((o) => o.textContent);
    expect(options.join('|'), 'all 3 matching', 'no empty option, and no target the server would refuse');
    buttonNamed('CANCEL').click();
    expect(await running, null);
  });

  test('the checked column is written into the selection only', async () => {
    await selectRows((id) => id !== ids[1]);
    expect(app.listSource.selection.value.length, 2, 'two of the three rows are selected');
    const {running} = await open();
    write('quantity', 9);
    buttonNamed('OK').click();
    expect(await running, 2);
    await reloaded();
    expect(await quantities(), '9,2,9', 'the unselected row kept its quantity');
  });

  test('"all matching" writes everything the list\'s filter matches', async () => {
    await selectRows(() => false);
    const {running} = await open();
    write('quantity', 4);
    buttonNamed('OK').click();
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

  test('an immutable column is refused by the server, and nothing is written', async () => {
    const {running} = await open();
    write('origin', 'rewritten');
    write('quantity', 7);
    buttonNamed('OK').click();
    expect(await running, null, 'the refusal reaches the caller');
    await reloaded();
    expect(await quantities(), '4,4,4', 'the whole call rolled back');
  });
});
