/* `DomainTable.app()` against real domain tables (2-9): the view `grok.shell.addView` docks — the
   list renders, search narrows, the filter query and the entity show in `view.path`, an edit on the
   entity page and a new event under it save as one transaction, Back through the gate keeps the
   page and the changes, the history pane lists the update, `open()` restores a deep link,
   `DomainApp.activate` answers with the same view, and closing with pending changes asks first.
   Fixture: `apitests.item` and `apitests.item_event` (item_id → item), rows prefixed per run and
   deleted afterwards. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {Rows} from '@datagrok-libraries/u2';
import type {DomainSource} from '@datagrok-libraries/u2';
import {domains, DomainApp, DomainTable} from '@datagrok-libraries/u2/src/dg/index.js';

const ITEMS = 'apitests.item';
const EVENTS = 'apitests.item_event';
const BASE = '/apps/U2Demo/items';

category('U2: domain app', () => {
  const items = () => grok.dapi.domains.table(ITEMS);
  const events = () => grok.dapi.domains.table(EVENTS);
  const prefix = `U2A-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
  const query = `sku starts "${prefix}"`;
  let itemTable: DomainTable;
  let eventTable: DomainTable;
  let view: DG.View;
  let app: DomainApp;
  let alphaId: string;

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

  const dialog = () => document.querySelector('.u2-dialog');
  const dialogButton = (text: string) =>
    [...document.querySelectorAll<HTMLButtonElement>('.u2-dialog button')].find((b) => b.textContent === text);

  before(async () => {
    const inserted = await items().insert([
      {sku: `${prefix}-1`, name: 'Alpha', quantity: 1},
      {sku: `${prefix}-2`, name: 'Beta', quantity: 2},
    ]);
    alphaId = inserted[0].id;
    itemTable = await domains.table(ITEMS);
    eventTable = await domains.table(EVENTS);
    view = itemTable.app({name: 'Items', path: BASE, query}) as DG.View;
    app = DomainApp.of(view)!;
    grok.shell.addView(view);
    // the harness calls grok.shell.closeAll() after every test, and closing the view kills the app
    // with it (the `data-kill-on-close` contract); this category is one story over one view
    view.temp['ignoreCloseAll'] = true;
    await ready(app.listSource);
  });

  after(async () => {
    app.session.discard();
    view.close();
    const mine = await items().query({filter: query});
    for (const item of mine)
      await events().deleteWhere(`item_id = "${item.id}"`);
    await items().deleteWhere(query);
  });

  test('the list renders; search narrows; the filter query shows in the path', async () => {
    expect(view.name, 'Items');
    expect(app.page.value, 'list');
    expect(app.listSource.rows.items.value.length >= 2, true);
    expect(app.list.root.querySelectorAll('.u2-list-row').length >= 2, true, 'rows rendered');
    expect(view.path, `${BASE}?q=${encodeURIComponent(query)}`);
    app.listSource.search.value = 'Beta';
    await delay(100);
    await ready(app.listSource);
    expect(app.listSource.rows.items.value.map((r) => r.name).join(), 'Beta');
    app.listSource.search.value = '';
    await delay(100);
    await ready(app.listSource);
    app.listSource.query.value = `${query} and quantity > 1`;
    expect(view.path.includes('?q='), true);
    expect(decodeURIComponent(view.path).endsWith('quantity > 1'), true);
    await delay(100);
    await ready(app.listSource);
    expect(app.listSource.rows.items.value.map((r) => r.name).join(), 'Beta');
    app.listSource.query.value = query;
    await delay(100);
    await ready(app.listSource);
  });

  test('opening a row switches to the entity page; an edit and a new child event save as one transaction', async () => {
    expect(await app.goTo('entity', alphaId), true);
    expect(app.page.value, 'entity');
    expect(view.path, `${BASE}?entity=${alphaId}`);
    const source = app.entitySource.value!;
    await ready(source);
    await until(() => app.form.value?.input('name') !== undefined, 'the form');
    expect(source.session === app.session, true, 'the entity source shares the session');
    app.form.value!.input('name')!.value.value = 'Alpha A';
    const kids = eventTable.source({query: `item_id = "${alphaId}"`, defaults: {item_id: alphaId}, session: app.session});
    try {
      await ready(kids);
      const event = kids.newRow({kind: 'made', amount: 5});
      expect(Rows.isDraft(event), true);
      expect(app.session.changeCount.value, 2);
      expect(app.summary.value, '2 unsaved changes in 2 tables');
      expect(await app.session.save(), true);
      expect(app.session.isDirty.value, false);
      expect(Rows.isDraft(kids.currentRow.value!), false, 'the exact post-save re-point');
      const [itemAudit, eventAudit] = await Promise.all([
        items().audit(alphaId), events().audit(kids.currentRow.value!.id)]);
      expect(itemAudit.at(-1)!.tx_id, eventAudit.at(-1)!.tx_id, 'one transaction in the audit');
      expect((await items().get(alphaId)).name, 'Alpha A');
    } finally {
      kids.dispose();
    }
  });

  test('Back through the gate: cancel keeps the page and the dirty state; the history pane lists the update', async () => {
    expect(app.page.value, 'entity');
    await until(() => app.form.value?.input('name') !== undefined, 'the form');
    app.form.value!.input('name')!.value.value = 'Alpha B';
    expect(app.session.isDirty.value, true);
    const back = app.goTo('list');
    await until(() => dialog() !== null, 'the unsaved-changes dialog');
    dialogButton('CANCEL')!.click();
    expect(await back, false);
    expect(app.page.value, 'entity');
    expect(app.session.isDirty.value, true);
    expect(app.form.value!.input('name')!.value.value, 'Alpha B');
    app.session.discard();
    expect(app.session.isDirty.value, false);
    const history = app.panes.querySelector('[data-u2="domain-history"]')!;
    await until(() => history.querySelectorAll('.u2-domain-history-line').length > 0, 'the history lines');
    const lines = [...history.querySelectorAll('.u2-domain-history-line')].map((l) => l.textContent ?? '');
    expect(lines.some((l) => l.includes('updated') && l.includes('Alpha A')), true, 'the update is listed');
  });

  test('open() restores a deep link; DomainApp.activate answers with the same view', async () => {
    expect(await app.open(`${BASE}?q=`), true);
    expect(app.page.value, 'list');
    expect(await app.open(`?entity=${alphaId}`), true);
    expect(app.page.value, 'entity');
    expect(app.entity.value, alphaId);
    expect(view.path, `${BASE}?entity=${alphaId}`);
    await app.goTo('list');
    expect(DomainApp.activate(BASE, alphaId), true);
    await delay(50);
    expect(grok.shell.v.dart === view.dart, true, 'the view comes to the front');
    expect(app.entity.value, alphaId);
    expect([...grok.shell.views].some((v) => v.dart === view.dart), true);
  });

  test('closing with pending changes asks first; dismiss keeps the view', async () => {
    await until(() => app.form.value?.input('name') !== undefined, 'the form');
    app.form.value!.input('name')!.value.value = 'Alpha C';
    expect(app.session.isDirty.value, true);
    view.close();
    await until(() => dialog() !== null, 'the unsaved-changes dialog on close');
    dialogButton('CANCEL')!.click();
    await delay(50);
    expect([...grok.shell.views].some((v) => v.dart === view.dart), true, 'the view stays');
    expect(app.session.isDirty.value, true);
    app.session.discard();
  });
});
