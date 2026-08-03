import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';

// ObjectHandler.renderGrid (ui-js-api WO-1): the grid-customization member, wired
// through both interop proxies (JS handler <- Dart dispatch, Dart meta <- JS call).
// Contract: presentation derives from column tags/semTypes only — the queryDf
// frames used here carry NO hidden '~item' object column. The legacy-handler
// guard (a JS handler without the member falls through to the Dart meta) is
// covered by the Dart dispatch suite (xamgle tests/domain_row_dispatch_test.dart).
category('ObjectHandler: renderGrid', () => {
  const items = () => grok.dapi.domains.table('apitests.item');
  const events = () => grok.dapi.domains.table('apitests.item_event');
  const typeOf = (h: any) => { try { return h.type; } catch { return null; } };

  test('base class defines a no-op default', async () => {
    expect(typeof (DG.ObjectHandler.prototype as any).renderGrid, 'function',
      'ObjectHandler.renderGrid missing from the base class');
  });

  test('registered JS handler wins dispatch, decorates a queryDf frame', async () => {
    let active = true;
    const seen: (_DG.DataFrame | null)[] = [];
    class ItemHandler extends DG.ObjectHandler {
      get type() { return 'apitests.item'; }
      isApplicable(x: any) { return active && x instanceof DG.DomainRow && x.typeName === 'apitests.item'; }
      renderGrid(grid: _DG.Grid, options?: {items?: _DG.DataFrame}) {
        seen.push(options?.items ?? null);
        grid.columns.byName('quantity')!.visible = false;
      }
    }
    const handler = new ItemHandler();
    DG.ObjectHandler.register(handler);
    const sku = `SKU-RG-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
    const [ins] = await items().insert({sku, name: 'RenderGrid probe', quantity: 3});
    try {
      const df = await items().queryDf({filter: `sku = "${sku}"`});
      const grid = DG.Grid.create(df);
      // The last-registered JS handler must lead the dispatch list for its type
      // (Dart-side demotion of the per-table meta) and round-trip identically.
      const winner = DG.ObjectHandler.list().find((h) => typeOf(h) === 'apitests.item');
      expect(winner === handler, true, 'registered handler does not win the apitests.item dispatch');
      winner!.renderGrid(grid, {items: df});
      expect(seen.length, 1, 'renderGrid not invoked');
      expect(seen[0] != null, true, 'options.items not delivered');
      expect(grid.columns.byName('quantity')!.visible, false, 'grid customization did not land');
    } finally {
      active = false; // neutralize the handler for the rest of the session
      await items().delete(ins.id);
    }
  });

  test('Dart meta renderGrid decorates a JS-created grid', async () => {
    const sku = `SKU-RG-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
    const kind = `rg-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
    const [item] = await items().insert({sku, name: 'RenderGrid parent'});
    try {
      await events().insert({item_id: item.id, kind, amount: 1});
      const df = await events().queryDf({filter: `kind = "${kind}"`});
      const grid = DG.Grid.create(df);
      // No JS handler is registered for item_event — the first handler of this
      // type is the per-table Dart DomainRowMeta behind an EntityMetaDartProxy.
      const meta = DG.ObjectHandler.list().find((h) => typeOf(h) === 'apitests.item_event');
      expect(meta != null, true, 'per-table Dart meta for apitests.item_event not registered');
      meta!.renderGrid(grid, {items: df});
      expect(grid.columns.byName('id')!.visible, false, 'system column not hidden');
      expect(df.col('item_id')!.getTag(DG.TAGS.FRIENDLY_NAME), 'Item',
        'ref-column caption not stamped');
    } finally {
      await items().delete(item.id); // cascades to item_event (onDelete: cascade)
    }
  });
}, {owner: 'askalkin@datagrok.ai'});
