/* Trash and restore against a real domain table (3-5, 3-8): a row deleted in the app leaves the
   list, the ⋯ menu's Trash puts the list on `deleted: 'only'` — the deleted rows, read-only, with
   Restore as the only row action and `?trash=1` in the view path — and Restore brings the row
   back into the live list. Needs the server's `deleted` query flag and `POST …/{id}/restore`
   (WO 3-1) and the js-api `DomainTableClient.restore` (WO 3-4). Fixture: `apitests.item`, rows
   prefixed per run and deleted afterwards. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {Rows} from '@datagrok-libraries/u2';
import type {DomainSource} from '@datagrok-libraries/u2';
import {domains, DomainApp, DomainTable} from '@datagrok-libraries/u2/src/dg/index.js';

const ITEMS = 'apitests.item';
const BASE = '/apps/U2Demo/trash';

category('U2: domain trash', () => {
  const items = () => grok.dapi.domains.table(ITEMS);
  const prefix = `U2T-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
  const query = `sku starts "${prefix}"`;
  let table: DomainTable;
  let view: DG.View;
  let app: DomainApp;
  let goneId: string;

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

  const names = () => app.listSource.rows.items.value.map((r) => r.name).sort().join();

  before(async () => {
    const inserted = await items().insert([
      {sku: `${prefix}-1`, name: 'Alpha', quantity: 1},
      {sku: `${prefix}-2`, name: 'Beta', quantity: 2},
    ]);
    goneId = inserted[0].id;
    table = await domains.table(ITEMS);
    view = table.app({name: 'Trash demo', path: BASE, query}) as DG.View;
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

  test('a deleted row leaves the list; the ⋯ menu offers Trash', async () => {
    expect(names(), 'Alpha,Beta');
    expect(app.menuActions().map((a) => a.name).join(), 'Trash');
    const row = app.listSource.rows.byKey(goneId)!;
    const del = app.list.actionsFor(row).find((a) => a.name === 'Delete')!;
    expect(del !== undefined, true, 'Delete is the row action on a live row');
    del.run();
    expect(app.session.isDirty.value, true);
    expect(await app.session.save(), true);
    await reloaded();
    expect(names(), 'Beta');
  });

  test('Trash: the deleted rows, read-only, with Restore — and ?trash=1 in the path', async () => {
    expect(await app.setTrash(true), true);
    await reloaded();
    expect(app.listSource.deleted.value, 'only');
    expect(view.path, `${BASE}?q=${encodeURIComponent(query)}&trash=1`);
    expect(names(), 'Alpha');
    expect(app.summary.value, '1 deleted row');
    const row = app.listSource.rows.byKey(goneId)!;
    expect(Rows.isDeleted(row), true, 'the row carries ~is_deleted');
    expect(app.listSource.access.value.can('edit'), false, 'the trash is read-only');
    expect(app.listSource.access.value.row(row).can('edit'), false, 'and no row can lift that');
    expect(app.list.actionsFor(row).map((a) => a.name).join(), 'Restore');
    expect(app.list.root.querySelector(`[data-u2-row="${goneId}"]`)!
      .classList.contains('u2-domain-list-deleted'), true);
    expect(app.menuActions().map((a) => a.name).join(), 'Exit trash');
  });

  test('Restore brings the row back into the live list', async () => {
    const row = app.listSource.rows.byKey(goneId)!;
    app.list.actionsFor(row)[0].run();
    await reloaded();
    expect(names(), '', 'the restored row left the trash');
    expect(await app.setTrash(false), true);
    await reloaded();
    expect(names(), 'Alpha,Beta');
    expect(view.path, `${BASE}?q=${encodeURIComponent(query)}`);
    expect((await items().get(goneId)).name, 'Alpha', 'and the row is readable again');
  });
});
