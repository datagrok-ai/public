import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test, before, after} from '@datagrok-libraries/test/src/test';

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
  const stamp = () => `${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
  /** The PLATFORM meta for [type] — a JS handler is never an EntityMetaDartProxy,
   * so this picks the Dart one whatever else is registered. */
  const dartMetaFor = (type: string): any => DG.ObjectHandler.list()
    .find((h) => h instanceof DG.EntityMetaDartProxy && typeOf(h) === type);
  /** Resolves to the thrown error, or null when the action succeeded. */
  const thrown = async (action: () => Promise<any>): Promise<any> => {
    try {
      await action();
      return null;
    } catch (e) {
      return e;
    }
  };

  // enableDomainDatabases is a client-profile setting shared with the rest of the
  // session — capture and restore it (precedent: domain_filters_test.dart).
  let flagBefore: any;

  before(async () => {
    // grok test's fresh client profile boots with enableDomainDatabases off, so
    // the per-table Dart metas never registered — set the flag and force the
    // registration so dispatch sees the same session state as a flagged client.
    flagBefore = (grok.shell.settings as any).enableDomainDatabases;
    (grok.shell.settings as any).enableDomainDatabases = true;
    await (window as any).grok_DomainRowMeta_RegisterPerTableMetas();
  });

  after(async () => {
    (grok.shell.settings as any).enableDomainDatabases = flagBefore;
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
      // There is no unregister API — retire the handler's type and applicability
      // so later type-keyed resolution and forEntity skip it. NB retirement does
      // not undo registration side effects (ElementRenderer.renderers keys on the
      // type captured AT registration time).
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

  // ─────────── DomainObjectHandler (ui-js-api WO-4) ───────────
  // The reflective per-table handler: registry-driven metadata, capability-gated
  // actions, and render defaults that DELEGATE to the platform's Dart meta — so
  // an override-nothing subclass is never a regression.

  test('DomainObjectHandler: overriding nothing renders like the platform', async () => {
    const handler = new DG.DomainObjectHandler('apitests.item');
    const sku = `SKU-DOH-${stamp()}`;
    const [ins] = await items().insert({sku, name: 'Handler probe', quantity: 7});
    try {
      // getById is the JS acquisition path for a DomainRow (there was none before).
      const row = await handler.getById(ins.id) as _DG.DomainRow;
      expect(row instanceof DG.DomainRow, true, 'getById must resolve a DG.DomainRow');
      expect(row.id, ins.id);
      expect(row.typeName, 'apitests.item');
      const dart = dartMetaFor('apitests.item');
      expect(dart != null, true, 'per-table Dart meta for apitests.item not registered');
      expect(handler.getCaption(row), dart.getCaption(row), 'caption differs from the platform default');
      for (const member of ['renderCard', 'renderMarkup', 'renderTooltip'])
        expect((handler as any)[member](row).outerHTML, dart[member](row).outerHTML,
          `${member} differs from the platform default`);
      expect(handler.isApplicable(row), true, 'own row not claimed');
      expect(handler.isApplicable(DG.SemanticValue.fromValueType(row, 'apitests.item')), true,
        'SemanticValue-wrapped row not claimed');
      expect(handler.isApplicable('apitests.item'), false, 'a bare string must not be claimed');
      expect(handler.deepLink(row)!.startsWith('/domains/apitests/item/'), true,
        `unexpected deep link: ${handler.deepLink(row)}`);
    } finally {
      await items().delete(ins.id);
    }
  });

  test('DomainObjectHandler: registering it keeps the platform meta (CRUD commands)', async () => {
    class ItemDoh extends DG.DomainObjectHandler {
      retired = false;
      constructor() { super('apitests.item'); }
      get type(): string { return this.retired ? 'apitests.item-doh-retired' : this.table; }
      isApplicable(x: any): boolean { return !this.retired && super.isApplicable(x); }
      renderCard(x: any): HTMLElement {
        const d = document.createElement('div');
        d.textContent = `custom ${this.rowOf(x)!.values.sku}`;
        return d;
      }
    }
    const handler = new ItemDoh();
    const sku = `SKU-DOH-${stamp()}`;
    const [ins] = await items().insert({sku, name: 'Dispatch probe'});
    DG.ObjectHandler.register(handler);
    try {
      const row = await handler.getById(ins.id) as _DG.DomainRow;
      expect(DG.ObjectHandler.list().find((h) => typeOf(h) === 'apitests.item') === handler, true,
        'the registered DomainObjectHandler must win the dispatch');
      // Platform CRUD commands live on the Dart meta and are registration-based:
      // a JS handler winning the RENDERING dispatch must not displace it.
      expect(dartMetaFor('apitests.item') != null, true,
        'the Dart meta (owner of the Edit/Delete/Share commands) was displaced');
      expect(DG.ObjectHandler.forEntity(row)!.renderCard(row).textContent, `custom ${sku}`,
        'the override did not win forEntity');
    } finally {
      // No unregister API — retire the type and applicability (see the note above).
      handler.retired = true;
      await items().delete(ins.id);
    }
  });

  test('DomainObjectHandler: reflective properties, detail tabs, and editor', async () => {
    const handler = new DG.DomainObjectHandler('apitests.item');
    const props = await handler.getProperties();
    for (const c of ['sku', 'name', 'quantity'])
      expect(props.some((p) => p.name === c), true, `rowProperties must describe ${c}`);
    const caps = await handler.capabilities();
    expect(caps.securityMode, 'row');
    const sku = `SKU-DOH-${stamp()}`;
    const kind = `doh-${stamp()}`;
    const [ins] = await items().insert({sku, name: 'Detail probe'});
    try {
      await events().insert({item_id: ins.id, kind, amount: 2});
      const row = await handler.getById(ins.id) as _DG.DomainRow;
      const tabs = await handler.getDetailTabs(row);
      const tab = tabs.find((t) => t.table === 'apitests.item_event');
      expect(tab != null, true, `child-table tab missing: ${JSON.stringify(tabs.map((t) => t.table))}`);
      expect(tab!.fkColumn, 'item_id');
      const children = await events().query({filter: tab!.filter});
      expect(children.length, 1, 'the detail-tab filter must scope the child table to the parent row');
      // Writable columns only — the form mirrors column security, and system
      // columns are never editable.
      const editor = await handler.renderEditor(row);
      const inputs = editor.querySelectorAll('.ui-input-root');
      expect(inputs.length, caps.writableColumns.length,
        `editor inputs must match writableColumns (${JSON.stringify(caps.writableColumns)})`);
      expect(editor.textContent!.includes('version'), false, 'system columns must not be editable');
      // A brand-new unsaved row is editable the same way (create form).
      expect((await handler.renderEditor()).querySelectorAll('.ui-input-root').length,
        caps.writableColumns.length, 'create form must offer the same writable columns');
    } finally {
      await items().delete(ins.id); // cascades to item_event
    }
  });

  test('DomainObjectHandler: row permissions and ribbon actions follow a real grant', async () => {
    const handler = new DG.DomainObjectHandler('apitests.item');
    const sku = `SKU-DOH-${stamp()}`;
    const [ins] = await items().insert({sku, name: 'Permission probe'});
    const current = await grok.dapi.users.current();
    const me = await grok.dapi.users.include('group').filter(`login = "${current.login}"`).first();
    expect(me?.group?.id != null, true, 'current user personal group not resolved');
    const group = me!.group.id;
    // The fixture's BASELINE grants must survive: revoking a Delete grant the
    // environment already had would 403 every later row delete in the suite.
    const hadDelete = (await items().grants())
      .some((g) => g.group.id === group && g.permission === 'Delete');
    try {
      const row = await handler.getById(ins.id) as _DG.DomainRow;
      const names = async () => (await handler.getRibbonActions(row)).map((a) => a.name);
      // Both directions of one real grant round-trip (grant/revoke drop the caches).
      const gate = async (granted: boolean) => {
        granted ? await items().grant(group, 'Delete') : await items().revoke(group, 'Delete');
        return {perms: await row.permissions(), actions: await names()};
      };
      const denied = await gate(false);
      const allowed = await gate(true);
      expect(typeof denied.perms.edit, 'boolean',
        `permissions() must resolve flags: ${JSON.stringify(denied.perms)}`);
      expect(allowed.perms.delete, true,
        `a Delete grant must reach row permissions: ${JSON.stringify(allowed.perms)}`);
      expect(allowed.actions.includes('Delete'), true, 'the Delete action must appear with the grant');
      // An admin SESSION short-circuits the probe to all-true (Auth.adminMode) —
      // the denial side is only observable where the probe actually ran.
      if (!denied.perms.delete)
        expect(denied.actions.includes('Delete'), false, 'Delete must be hidden without the grant');
      // Always offered; Share only for row-mode tables (apitests.item is one).
      for (const n of ['History', 'Copy link', 'Share...'])
        expect(allowed.actions.includes(n), true, `${n} action missing: ${JSON.stringify(allowed.actions)}`);
      // Unsaved rows have no securing entity, hence no permissions — for admins too.
      const p = await handler.newRow().permissions();
      expect(p.edit || p.delete || p.share, false,
        `an unsaved row must hold no permissions: ${JSON.stringify(p)}`);
      expect((await handler.getRibbonActions(handler.newRow() as any)).some((a) => a.name === 'Delete'),
        false, 'an unsaved row must offer no Delete');
    } finally {
      // Delete needs the grant; restore the baseline only afterwards.
      try {
        await items().grant(group, 'Delete');
      } catch (_) { /* the probe above may have failed before the grant existed */ }
      await items().delete(ins.id);
      if (!hadDelete)
        await items().revoke(group, 'Delete');
      grok.dapi.domains.invalidateUiCaches();
    }
  });

  test('DomainObjectHandler: malformed and unknown tables fail loudly', async () => {
    let ctor: any = null;
    try {
      new DG.DomainObjectHandler('nodot');
    } catch (e) {
      ctor = e;
    }
    expect(ctor instanceof Error, true, 'a malformed address must throw at construction');
    const unknown = new DG.DomainObjectHandler('apitests.nosuch');
    for (const [member, action] of [['getProperties', () => unknown.getProperties()],
      ['capabilities', () => unknown.capabilities()]] as [string, () => Promise<any>][]) {
      const e = await thrown(action);
      expect(e instanceof DG.DomainValidationError, true,
        `${member} on an unknown table must reject with DomainValidationError, got ${e?.constructor?.name}: ${e?.message}`);
    }
  });
}, {owner: 'askalkin@datagrok.ai'});
