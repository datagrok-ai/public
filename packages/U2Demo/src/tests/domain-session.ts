/* The u2 session against real domain tables (2-4b): two sources under one explicit `SharedSession`
   save a draft parent and a draft child referencing it as ONE transaction; `domains.grid` hosts the
   platform grid with the source's js-api editor attached; `search` narrows the rows and the total;
   `confirmDiscard` over a clean session answers without a dialog. Fixture: `apitests.item` and
   `apitests.item_event` (item_id → item), rows prefixed per run and deleted afterwards. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {Rows, SharedSession, confirmDiscard} from '@datagrok-libraries/u2';
import type {DomainSource} from '@datagrok-libraries/u2';
import {domains, DomainGrid, DomainTable} from '@datagrok-libraries/u2/src/dg/index.js';

const ITEMS = 'apitests.item';
const EVENTS = 'apitests.item_event';

category('U2: domain session', () => {
  const items = () => grok.dapi.domains.table(ITEMS);
  const events = () => grok.dapi.domains.table(EVENTS);
  const prefix = `U2S-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
  const query = `sku starts "${prefix}"`;
  let itemTable: DomainTable;
  let eventTable: DomainTable;

  async function ready(src: DomainSource): Promise<void> {
    for (let i = 0; i < 200 && (src.state.value === 'idle' || src.state.value === 'loading'); i++)
      await delay(50);
    if (src.state.value !== 'ready')
      throw new Error(`source is ${src.state.value}: ${String(src.error.value)}`);
  }

  before(async () => {
    await items().insert([
      {sku: `${prefix}-1`, name: 'Alpha', quantity: 1},
      {sku: `${prefix}-2`, name: 'Beta', quantity: 2},
    ]);
    itemTable = await domains.table(ITEMS);
    eventTable = await domains.table(EVENTS);
  });

  after(async () => {
    const mine = await items().query({filter: query});
    for (const item of mine)
      await events().deleteWhere(`item_id = "${item.id}"`);
    await items().deleteWhere(query);
  });

  test('one session, two tables: a draft item and a draft event referencing it land as one transaction', async () => {
    const session = new SharedSession();
    const parent = itemTable.source({query, pageSize: 10, session});
    let kids: DomainSource | undefined;
    try {
      await ready(parent);
      const draft = parent.newRow({sku: `${prefix}-3`, name: 'Gamma', quantity: 3});
      expect(Rows.isDraft(draft), true);
      // the master–detail shape: the child collection is queried by the parent's id, draft or real
      kids = eventTable.source({query: `item_id = "${draft.id}"`, defaults: {item_id: draft.id},
        pageSize: 10, session});
      await ready(kids);
      expect(kids.rows.items.value.length, 0, 'a query naming a draft id is not sent to the server');
      expect(session.sources.value.length, 2);
      const event = kids.newRow({kind: 'made', amount: 5});
      expect(event.item_id, draft.id, 'the child holds the parent\'s draft id');
      expect(session.isDirty.value, true);
      expect(session.changeCount.value, 2);
      expect(session.summary.value, '2 unsaved changes in 2 tables');
      expect(await session.save(), true);
      expect(session.isDirty.value, false);
      const savedItem = parent.currentRow.value!;
      expect(Rows.isDraft(savedItem), false, 'the current draft re-pointed to its row');
      expect(kids.query.value, `item_id = "${savedItem.id}"`, 'the child query names the saved parent');
      const savedEvent = kids.currentRow.value!;
      expect(Rows.isDraft(savedEvent), false);
      const fresh = await events().get(savedEvent.id);
      expect(fresh.item_id, savedItem.id, 'the event\'s item_id is the item\'s real id');
      expect(fresh.kind, 'made');
      expect((await items().get(savedItem.id)).name, 'Gamma');
      expect(kids.newRow({kind: 'used', amount: 1}).item_id, savedItem.id, 'a later draft takes the real id');
    } finally {
      parent.dispose();
      kids?.dispose();
    }
  });

  test('domains.grid: the source\'s editor attached, service columns hidden, a cell edit tracked and saved', async () => {
    const src = itemTable.source({query, pageSize: 10});
    let grid: DomainGrid | undefined;
    try {
      await ready(src);
      grid = domains.grid(src);
      const g = grid.grid;
      expect(g.dataFrame.dart === (src.df.value as DG.DataFrame).dart, true, 'over the source\'s frame');
      expect(g.editor !== null, true, 'an editor is attached');
      expect(g.editor === (src.edit.value as any).editor, true, 'the source\'s own');
      for (const col of g.dataFrame.columns) {
        if (col.name.startsWith('~'))
          expect(g.col(col.name)?.visible ?? false, false, `${col.name} hidden`);
      }
      expect(g.columns.byName('name')!.idx < g.columns.byName('quantity')!.idx, true, 'the name column first');
      const df = g.dataFrame;
      const at = df.col('name')!.toList().indexOf('Alpha');
      expect(at >= 0, true);
      g.cell('name', at).setValue('Alpha G');
      await delay(50);
      expect(src.isDirty.value, true, 'the in-cell edit reached the source');
      const id = df.get('id', at) as string;
      expect(src.edit.value!.isChanged(id, 'name'), true);
      expect(await src.session.save(), true);
      expect((await items().get(id)).name, 'Alpha G');
    } finally {
      grid?.dispose();
      src.dispose();
      const alpha = await items().query({filter: `sku = "${prefix}-1"`});
      if (alpha[0] !== undefined)
        await items().update(alpha[0].id, {name: 'Alpha'});
    }
  });

  test('search narrows the rows and the total', async () => {
    const src = itemTable.source({query, pageSize: 10});
    try {
      await ready(src);
      const all = src.rows.items.value.length;
      expect(all >= 2, true);
      expect(src.total.value, all);
      src.search.value = 'Beta';
      await delay(100);
      await ready(src);
      expect(src.rows.items.value.map((r) => r.name).join(), 'Beta');
      expect(src.total.value, 1, 'the count agrees with the rows');
      src.search.value = '';
      await delay(100);
      await ready(src);
      expect(src.rows.items.value.length, all);
    } finally {
      src.dispose();
    }
  });

  test('confirmDiscard over a clean session answers true without a dialog', async () => {
    const src = itemTable.source({query, pageSize: 10});
    try {
      await ready(src);
      const dialogs = document.querySelectorAll('.u2-dialog').length;
      expect(await confirmDiscard(src.session), true);
      expect(document.querySelectorAll('.u2-dialog').length, dialogs, 'no dialog');
    } finally {
      src.dispose();
    }
  });
});
