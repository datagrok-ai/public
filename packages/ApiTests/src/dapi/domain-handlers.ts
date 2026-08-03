import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test, before} from '@datagrok-libraries/test/src/test';

// DG.DomainRow JS wrapper + ObjectHandler dispatch for domain-table rows (§8.1),
// including ObjectHandler.renderGrid (ui-js-api WO-1) — the grid-customization
// member wired through both interop proxies (JS handler <- Dart dispatch, Dart
// meta <- JS call). Contract: presentation derives from column tags/semTypes
// only — the queryDf frames used here carry NO hidden '~item' object column.
// The inherited-no-op guard (a class-based handler without a renderGrid
// override falls through to the Dart meta) is covered by the Dart dispatch
// suite (xamgle tests/domain_row_dispatch_test.dart).
category('JS: domain handlers', () => {
  const items = () => grok.dapi.domains.table('apitests.item');
  const events = () => grok.dapi.domains.table('apitests.item_event');
  const typeOf = (h: any) => { try { return h.type; } catch { return null; } };

  before(async () => {
    // grok test's fresh client profile boots with enableDomainDatabases off, so
    // the per-table Dart metas never registered — set the flag and force the
    // registration so dispatch sees the same session state as a flagged client.
    (grok.shell.settings as any).enableDomainDatabases = true;
    await (window as any).grok_DomainRowMeta_RegisterPerTableMetas();
  });

  test('DG.DomainRow wrapper exported', async () => {
    expect(typeof DG.DomainRow, 'function', 'DG.DomainRow is not exported');
    const props = Object.getOwnPropertyNames((DG.DomainRow as any).prototype);
    for (const p of ['schemaName', 'tableName', 'typeName', 'semValue', 'values', 'id', 'permissions'])
      expect(props.includes(p), true, `DG.DomainRow.${p} member missing`);
  });

  test('grit.issue ObjectHandler resolves', async () => {
    // Requires the Grit package (which registers the grit.issue handler) to be
    // deployed; skip cleanly where it is not so the suite stays green everywhere.
    const handlers = await DG.ObjectHandler.forSemType('grit.issue');
    if (handlers.length === 0)
      return;
    expect(handlers.some((h) => h.type === 'grit.issue'), true,
      'grit.issue handler not resolved by forSemType');
  });

  test('renderGrid: base class defines a sentinel-marked no-op', async () => {
    const base = (DG.ObjectHandler.prototype as any).renderGrid;
    expect(typeof base, 'function', 'ObjectHandler.renderGrid missing from the base class');
    // Without the sentinel, EVERY class-based handler inherits a member that
    // wins the Dart dispatch and silently swallows the default decoration.
    expect(base.isPlatformDefault, true, 'base renderGrid is not marked isPlatformDefault');
  });

  test('renderGrid: registered JS handler wins dispatch, decorates a queryDf frame', async () => {
    const seen: (_DG.DataFrame | null)[] = [];
    class ItemHandler extends DG.ObjectHandler {
      retired = false;
      get type() { return this.retired ? 'apitests.item-retired' : 'apitests.item'; }
      isApplicable(x: any) { return !this.retired && x instanceof DG.DomainRow && x.typeName === 'apitests.item'; }
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
      // There is no unregister API — retire the handler's type so later
      // type-keyed resolution in this session can never pick it up again.
      handler.retired = true;
      await items().delete(ins.id);
    }
  });

  test('renderGrid: Dart meta renderGrid decorates a JS-created grid', async () => {
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
