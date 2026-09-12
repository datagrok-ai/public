/* The u2 domain stack against a real domain table (WO-8..11): `domains.table()` over
   `grok.dapi.domains`, a `DomainSource` over the frame the backend fetched with the js-api editor
   attached, edit → save as one transaction, discard, a pristine draft inserted, paging into the
   same frame, and the form, list and picker over it. Fixture: `apitests.item` (sku, name,
   quantity), rows prefixed per run and deleted afterwards. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {Rows, Scope} from '@datagrok-libraries/u2';
import type {DomainSource} from '@datagrok-libraries/u2';
import {domains, domainForm, domainList, DomainPick, DomainTable, saveButton, discardButton}
  from '@datagrok-libraries/u2/src/dg/index.js';

const TABLE = 'apitests.item';

category('U2: domain source', () => {
  const items = () => grok.dapi.domains.table(TABLE);
  const prefix = `U2-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
  const query = `sku starts "${prefix}"`;
  const ids: string[] = [];
  let table: DomainTable;

  /** The first load settled, or the reason it did not. */
  async function ready(src: DomainSource): Promise<void> {
    for (let i = 0; i < 200 && (src.state.value === 'idle' || src.state.value === 'loading'); i++)
      await delay(50);
    if (src.state.value !== 'ready')
      throw new Error(`source is ${src.state.value}: ${String(src.error.value)}`);
  }

  async function source(options: {pageSize?: number} = {}): Promise<DomainSource> {
    const src = table.source({query, pageSize: options.pageSize ?? 10});
    await ready(src);
    return src;
  }

  before(async () => {
    const inserted = await items().insert([
      {sku: `${prefix}-1`, name: 'Alpha', quantity: 1},
      {sku: `${prefix}-2`, name: 'Beta', quantity: 2},
    ]);
    ids.push(...inserted.map((r) => r.id));
    table = await domains.table(TABLE);
  });

  after(async () => {
    await items().deleteWhere(query);
  });

  test('handle: schema, info and access in one await', async () => {
    expect(table.address, TABLE);
    expect(table.info.nameColumn, 'name');
    expect(table.properties.some((p) => p.name === 'quantity'), true);
    expect(table.access.can('edit'), true);
    expect(table.access.can('insert'), true);
    expect(table.access.field('sku'), 'editable');
    expect(table.access.field('id'), 'readonly', 'a system column is listed, read-only');
    expect(table.access.field('nosuch'), 'hidden');
    expect(table.handler instanceof DG.DomainObjectHandler, true);
  });

  test('load: rows over a frame with the editor attached; the access columns ride the rows', async () => {
    const live = Scope.liveCount;
    const src = await source();
    try {
      expect(src.df.value instanceof DG.DataFrame, true, 'the frame path');
      expect(src.rows.items.value.length, 2);
      expect(src.total.value, 2);
      const row = src.rows.byKey(ids[0])!;
      expect(row.name, 'Alpha');
      expect(row.quantity, 1);
      expect(row['~can_edit'], true);
      expect(row['~can_delete'], true);
      expect(Object.keys(row).some((k) => k.startsWith('~')), false, 'service columns are read by name only');
      expect(src.edit.value!.isChanged(ids[0], 'name'), false);
      expect(src.isDirty.value, false);
      expect(src.summary.value, '2 items');
      const df = src.df.value as DG.DataFrame;
      for (const name of DG.DOMAIN_ACCESS_COLUMNS)
        expect(df.col(name)?.meta.includeInCsvExport, false, `${name} is tagged out of export`);
    } finally {
      src.dispose();
    }
    expect(Scope.liveCount, live, 'no scope left behind');
  });

  test('edit then save writes one transaction; discard restores the cell', async () => {
    const src = await source();
    try {
      const row = src.rows.byKey(ids[0])!;
      row.name = 'Alpha 2';
      expect(src.isDirty.value, true);
      expect(src.changeCount.value, 1);
      expect(src.edit.value!.isChanged(ids[0], 'name'), true);
      expect(await src.save(), true);
      expect(src.isDirty.value, false);
      const fresh = await items().get(ids[0]);
      expect(fresh.name, 'Alpha 2');
      expect(src.rows.byKey(ids[0])!.version, fresh.version, 'the version came back into the frame');
      row.name = 'Zzz';
      expect(src.isDirty.value, true);
      src.discard();
      expect(src.isDirty.value, false);
      expect(src.rows.byKey(ids[0])!.name, 'Alpha 2');
    } finally {
      src.dispose();
      await items().update(ids[0], {name: 'Alpha'});
    }
  });

  test('draft: a pristine draft, inserted by save', async () => {
    const draft = table.draft({sku: `${prefix}-3`, name: 'Gamma'});
    try {
      await ready(draft);
      const row = draft.currentRow.value;
      expect(row !== null, true, 'the draft is the current row');
      expect(Rows.isDraft(row!), true);
      expect(draft.isDraft, true);
      expect(draft.isDirty.value, false, 'pristine until touched');
      row!.quantity = 3;
      expect(draft.isDirty.value, true);
      expect(await draft.save(), true);
      const [saved] = await items().query({filter: `sku = "${prefix}-3"`});
      expect(saved?.quantity, 3);
      expect(saved?.name, 'Gamma');
      ids.push(saved.id);
    } finally {
      draft.dispose();
    }
  });

  test('loadMore appends the next page into the same frame', async () => {
    const src = await source({pageSize: 1});
    try {
      const df = src.df.value;
      expect(src.rows.items.value.length, 1);
      await src.loadMore();
      expect(src.df.value === df, true, 'the same frame');
      expect(src.rows.items.value.length, 2);
      expect(src.state.value, 'ready');
    } finally {
      src.dispose();
    }
  });

  test('form: an input per editable column, text for the rest, writes through the source', async () => {
    const src = await source();
    src.currentRow.value = src.rows.byKey(ids[0]);
    const form = domainForm(src);
    try {
      expect(form.form !== null, true);
      expect(form.input('name') !== undefined, true);
      expect(form.input('id') === undefined, true, 'a system column is not a field');
      const status = form.form!.getWidgetStatus().inputs;
      expect(status.find((f) => f.name === 'id') === undefined, true);
      expect(status.find((f) => f.name === 'name')?.access, 'editable');
      form.input('name')!.value.value = 'Alpha 3';
      await delay(10);
      expect(src.isDirty.value, true);
      expect(src.rows.byKey(ids[0])!.name, 'Alpha 3');
      src.discard();
      await delay(10);
      expect(form.input('name')!.value.value, 'Alpha', 'the form re-read the row');
      expect(src.isDirty.value, false);
    } finally {
      form.dispose();
      src.dispose();
    }
  });

  test('list: rows through the handler with Open and Delete; selection is the current row', async () => {
    const src = await source();
    const list = domainList(src);
    try {
      const names = list.actionsFor(src.rows.byKey(ids[0])!).map((a) => a.name);
      expect(names.includes('Open') && names.includes('Delete'), true, names.join(', '));
      list.list.selectedIndex.value = 1;
      await delay(10);
      expect(src.currentRow.value?.id, src.rows.items.value[1].id);
      src.currentRow.value = null;
      await delay(10);
      expect(list.list.selectedIndex.value, -1);
      list.actionsFor(src.rows.byKey(ids[1])!).find((a) => a.name === 'Delete')!.run();
      await delay(10);
      expect(src.rows.byKey(ids[1])!['~state'], 'deleted', 'marked deleted: kept in the rows, struck through');
      expect(src.summary.value, '1 deletion pending');
      expect(list.actionsFor(src.rows.byKey(ids[1])!).map((a) => a.name).join(), 'Restore');
      src.discard();
      await delay(10);
      expect(src.rows.byKey(ids[1])!['~state'] ?? '', '', 'discard restores it (an empty state cell reads null)');
    } finally {
      list.dispose();
      src.dispose();
    }
  });

  test('pick: candidates by the name column, an id resolved to its name', async () => {
    const found = await DomainPick.search(TABLE, 'Alp', {filter: query});
    expect(found.length, 1, JSON.stringify(found));
    expect(found[0].id, ids[0]);
    expect(found[0].name, 'Alpha');
    const resolved = await DomainPick.resolve(TABLE, ids[1]);
    expect(resolved.name, 'Beta');
  });

  test('Save and Discard buttons follow the session', async () => {
    const src = await source();
    const save = saveButton(src);
    const discard = discardButton(src);
    try {
      expect(save.session === src.session, true, 'the source\'s own session');
      expect(save.button.disabled, true);
      src.rows.byKey(ids[1])!.quantity = 5;
      await delay(10);
      expect(save.button.disabled, false);
      expect(discard.button.disabled, false);
      discard.button.click();
      await delay(10);
      expect(src.isDirty.value, false);
      expect(save.button.disabled, true);
    } finally {
      save.dispose();
      discard.dispose();
      src.dispose();
    }
  });
});
