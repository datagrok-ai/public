import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test, before, after, awaitCheck} from '@datagrok-libraries/test/src/test';

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

  test('renderGrid: a plain handler falls through to the platform from JS too', async () => {
    // The base ObjectHandler.renderGrid delegates to the platform meta for its
    // type (Meta_DartForType), so a handler that overrides NOTHING decorates a
    // grid exactly like the built-in view — the JS-side mirror of the Dart
    // fall-through. Not registered: this asserts the member, not the dispatch.
    class PlainHandler extends DG.ObjectHandler {
      get type() { return 'apitests.item_event'; }
      isApplicable(_x: any) { return false; }
    }
    const sku = `SKU-RG-${stamp()}`;
    const kind = `rg-plain-${stamp()}`;
    const [item] = await items().insert({sku, name: 'Root fall-through parent'});
    try {
      await events().insert({item_id: item.id, kind, amount: 1});
      const decorate = async (render: (g: _DG.Grid, df: _DG.DataFrame) => void) => {
        const df = await events().queryDf({filter: `kind = "${kind}"`});
        const grid = DG.Grid.create(df);
        render(grid, df);
        const visible: string[] = [];
        for (let i = 0; i < grid.columns.length; i++) {
          const gc = grid.columns.byIndex(i);
          if (gc?.visible)
            visible.push(gc.name);
        }
        return {visible: visible, caption: df.col('item_id')!.getTag(DG.TAGS.FRIENDLY_NAME)};
      };
      const dart = dartMetaFor('apitests.item_event');
      expect(dart != null, true, 'per-table Dart meta for apitests.item_event not registered');
      const viaPlatform = await decorate((g, df) => dart.renderGrid(g, {items: df}));
      const viaHandler = await decorate((g, df) => new PlainHandler().renderGrid(g, {items: df}));
      expect(viaHandler.visible.join(','), viaPlatform.visible.join(','),
        `a plain handler must decorate like the platform meta (${viaHandler.visible.join(',')})`);
      expect(viaHandler.visible.includes('id'), false, 'system columns must be hidden');
      expect(viaHandler.caption, viaPlatform.caption, 'ref-column caption not stamped');
    } finally {
      await items().delete(item.id); // cascades to item_event
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
    // Everything after the insert runs inside the try — a throw in any of the
    // lookups below would otherwise leak the row.
    let group: string | undefined;
    let hadDelete = false;
    try {
      const current = await grok.dapi.users.current();
      const me = await grok.dapi.users.include('group').filter(`login = "${current.login}"`).first();
      expect(me?.group?.id != null, true, 'current user personal group not resolved');
      group = me!.group.id;
      // The fixture's BASELINE grants must survive: revoking a Delete grant the
      // environment already had would 403 every later row delete in the suite.
      hadDelete = (await items().grants())
        .some((g) => g.group.id === group && g.permission === 'Delete');
      const row = await handler.getById(ins.id) as _DG.DomainRow;
      const names = async () => (await handler.getRibbonActions(row)).map((a) => a.name);
      // Both directions of one real grant round-trip (grant/revoke drop the caches).
      const gate = async (granted: boolean) => {
        granted ? await items().grant(group!, 'Delete') : await items().revoke(group!, 'Delete');
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
      // Open (the platform's default row action), History and Copy link are
      // always offered; Share needs BOTH a row-mode table and Share on the row.
      for (const n of ['Open', 'History', 'Copy link'])
        expect(allowed.actions.includes(n), true, `${n} action missing: ${JSON.stringify(allowed.actions)}`);
      expect(allowed.actions.includes('Share...'), allowed.perms.share === true,
        `Share... must follow the row's Share permission: ${JSON.stringify(allowed.perms)}`);
      // Unsaved rows have no securing entity, hence no permissions — for admins
      // too — and therefore no actions at all (no address, history or link).
      const p = await handler.newRow().permissions();
      expect(p.edit || p.delete || p.share, false,
        `an unsaved row must hold no permissions: ${JSON.stringify(p)}`);
      expect((await handler.getRibbonActions(handler.newRow() as any)).length, 0,
        'an unsaved row must offer no actions');
      expect(handler.deepLink(handler.newRow() as any), null,
        'an unsaved row must have no deep link');
    } finally {
      // Delete needs the grant; restore the baseline only afterwards.
      if (group != null)
        try {
          await items().grant(group, 'Delete');
        } catch (e) {
          // Never silent: a failed restore is exactly how a baseline grant gets lost.
          console.error(`restoring the Delete grant on apitests.item failed: ${e}`);
        }
      await items().delete(ins.id);
      if (group != null && !hadDelete)
        await items().revoke(group, 'Delete');
      grok.dapi.domains.invalidateUiCaches();
    }
  });

  // ─────────── DG.DomainView + the dialog openers (ui-js-api WO-5) ───────────
  // Every opener enters the SAME Dart flow the built-in surfaces use, so these
  // tests drive the real platform dialogs (find them by title, click their
  // buttons) instead of asserting that a method exists.

  /** The open dialog whose title contains [fragment]. */
  const openDialog = async (fragment: string): Promise<_DG.Dialog> => {
    const find = () => DG.Dialog.getOpenDialogs().find((d) => `${d.title}`.includes(fragment));
    await awaitCheck(() => find() != null, `dialog "${fragment}" did not open`, 10000);
    return find()!;
  };
  const click = (dlg: _DG.Dialog, text: string) => dlg.getButton(text).click();
  /** Closes anything still on screen — a half-open dialog fails every later test. */
  const closeDialogs = () => DG.Dialog.getOpenDialogs().forEach((d) => d.close());

  test('DomainView: lists rows, honors permanentFilter, reports its query', async () => {
    const mine = `SKU-DV-${stamp()}`;
    const other = `SKU-DV-${stamp()}`;
    // One insert inside the try so no partial pair can leak; sku and name carry
    // the same stamp, so whichever the table's name column is, the rendered row
    // text contains it.
    let inserted: {id: string}[] = [];
    let view: _DG.DomainView | null = null;
    try {
      inserted = await items().insert([{sku: mine, name: mine}, {sku: other, name: other}]);
      // Embedded: a brief DOM list instead of the canvas grid, so the rendered
      // rows are readable from the test (and the render mode is not remembered).
      view = DG.DomainView.create({schema: 'apitests', table: 'item',
        permanentFilter: `sku = "${mine}"`, embedded: true});
      grok.shell.addView(view);
      expect(view instanceof DG.DomainView, true, 'create must resolve a DG.DomainView');
      expect(view.permanentFilter, `sku = "${mine}"`, 'permanentFilter did not round-trip');
      await awaitCheck(() => view!.root.textContent!.includes(mine),
        'the permanently filtered row never appeared in the view', 15000);
      expect(view.root.textContent!.includes(other), false,
        'a row outside permanentFilter is listed');
      const q = view.query;
      expect(q.schema, 'apitests');
      expect(q.table, 'item');
      expect(q.limit! > 0, true, `query must carry the view's row cap: ${JSON.stringify(q)}`);
      // The reported query is the WHOLE current state: the invisible
      // permanentFilter is merged into the filter elements like everything else.
      expect((q.filters ?? []).some((f) => f.includes(mine)), true,
        `permanentFilter missing from the reported query: ${JSON.stringify(q)}`);
      // ...and so is what the user typed into the search box.
      view.searchValue = 'quantity > 0';
      await awaitCheck(() => (view!.query.filters ?? []).some((f) => f.includes('quantity > 0')),
        `searchValue never reached the reported query: ${JSON.stringify(view.query)}`, 15000);
      view.searchValue = '';
      // The inherited CardView refresh() re-runs the query behind the view. The
      // row is changed on the SERVER first, so only a re-query can show the new
      // value — asserting the row is still listed would pass without one.
      const renamed = `Renamed ${stamp()}`;
      await items().update(inserted[0].id, {name: renamed});
      view.refresh();
      await awaitCheck(() => view!.root.textContent!.includes(renamed),
        'refresh() did not re-run the query', 15000);
      view.showFilters();
      await awaitCheck(() => view!.root.querySelector('.grok-domains-filter-panel') != null,
        'showFilters did not dock the filter panel', 15000);
      view.showFilters();      // idempotent
      // The registry is the view's meta source — assigning one is refused by name.
      let assigned = '';
      try {
        (view as any).meta = null;
      } catch (e: any) {
        assigned = e.message ?? `${e}`;
      }
      expect(assigned.includes('domain registry'), true,
        `DomainView.meta must fail with a named error, got: '${assigned}'`);
    } finally {
      if (view != null)
        view.close();
      for (const r of inserted)
        await items().delete(r.id);
    }
  });

  test('openers: create dialog opens, cancel resolves false', async () => {
    const info = await grok.dapi.domains.registry.tableInfo('apitests.item');
    const pending = DG.DomainObjectHandler.createRow('apitests.item');
    try {
      const dlg = await openDialog(`New ${info.singularName}`);
      click(dlg, 'CANCEL');
      expect(await pending, false, 'a cancelled create must resolve false');
    } finally {
      closeDialogs();
    }
  });

  test('openers: edit dialog saves and resolves true', async () => {
    const handler = new DG.DomainObjectHandler('apitests.item');
    const [ins] = await items().insert({sku: `SKU-ED-${stamp()}`, name: 'Editor probe'});
    try {
      const row = await handler.getById(ins.id) as _DG.DomainRow;
      // The instance member delegates to the static: exercising it covers both.
      const pending = handler.editRow(row as any);
      const dlg = await openDialog('Edit apitests.item');
      // Captions are the platform's column labels, not the wire names.
      const input = dlg.inputs.find((i) => `${i.caption}`.toLowerCase() === 'name');
      expect(input != null, true,
        `the editor has no input for the name column: ${JSON.stringify(dlg.inputs.map((i) => i.caption))}`);
      const edited = `Edited ${stamp()}`;
      input!.value = edited;
      click(dlg, 'OK');
      expect(await pending, true, 'saving the edit dialog must resolve true');
      // A real versioned update went through the platform's own save path.
      const fresh = await items().get(ins.id);
      expect(fresh!.name, edited, 'the edited value did not reach the server');
      expect(fresh!.version > ins.version!, true,
        `version not bumped (${ins.version} -> ${fresh!.version})`);
    } finally {
      closeDialogs();
      await items().delete(ins.id);
    }
  });

  test('openers: picker hosts the target Domain View, cancel resolves null', async () => {
    const pending = DG.DomainObjectHandler.pickRow('apitests.item');
    try {
      const dlg = await openDialog('Select apitests.item');
      click(dlg, 'CANCEL');
      expect(await pending, null, 'a cancelled pick must resolve null');
    } finally {
      closeDialogs();
    }
  });

  test('openers: conflict dialog resolves the decision', async () => {
    for (const [button, expected] of [['RELOAD', 'reload'], ['OVERWRITE', 'overwrite'],
      ['CANCEL', null]] as [string, string | null][]) {
      const pending = DG.DomainObjectHandler.showConflictDialog('SKU-CONFLICT');
      try {
        click(await openDialog('Conflicting changes'), button);
        expect(await pending, expected, `${button} must resolve ${expected}`);
      } finally {
        closeDialogs();
      }
    }
  });

  test('openers: delete confirmation, audit pane, grants pane, navigation', async () => {
    const handler = new DG.DomainObjectHandler('apitests.item');
    const [ins] = await items().insert({sku: `SKU-OP-${stamp()}`, name: 'Openers probe'});
    const previous = grok.shell.v;
    try {
      const row = await handler.getById(ins.id) as _DG.DomainRow;
      const pending = handler.deleteRow(row as any);
      click(await openDialog('Are you sure?'), 'CANCEL');
      expect(await pending, false, 'a cancelled delete must resolve false');

      // The audit trail of a just-inserted row: one 'insert' event.
      const audit = handler.auditPane(row as any);
      expect(audit instanceof HTMLElement, true, 'auditPane must return an element');
      await awaitCheck(() => audit.textContent!.includes('insert'),
        `the audit pane did not render the insert event: "${audit.textContent}"`, 10000);

      // Grants of the table's REGISTRY entity (rows share through shareRow).
      const schema = await grok.dapi.domains.schemas.include('tables').filter('name = "apitests"').first();
      const table = schema.tables.find((t) => t.name === 'item')!;
      const grants = DG.DomainObjectHandler.grantsPane(table.id, 'apitests.item', {readOnly: true});
      expect(grants instanceof HTMLElement, true, 'grantsPane must return an element');
      await awaitCheck(() => grants.textContent!.length > 0, 'the grants pane stayed empty', 10000);

      // Navigation goes through the platform's own routing (one spelling).
      expect(DG.DomainObjectHandler.deepLink(row), handler.deepLink(row as any),
        'the static and instance deep links disagree');
      handler.openRow(row as any);
      await awaitCheck(() => `${grok.shell.v?.path ?? ''}`.startsWith('/domains/apitests/item/'),
        `openRow did not open the Entity View (at "${grok.shell.v?.path}")`, 15000);
    } finally {
      closeDialogs();
      if (grok.shell.v !== previous) {
        grok.shell.v.close();
        if (previous != null)
          grok.shell.v = previous;
      }
      await items().delete(ins.id);
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

    // rowFrom takes ONE row of a query() result. Anything else is read key by
    // key into a row of nonsense values that only surfaces much later, as a card
    // with no identity — so it is rejected by name, right where it was passed.
    const handler = new DG.DomainObjectHandler('apitests.item');
    for (const garbage of ['garbage', 42, [{id: 1}], new Date()] as any[]) {
      let rejected: any = null;
      try {
        handler.rowFrom(garbage);
      } catch (e) {
        rejected = e;
      }
      expect(rejected instanceof Error, true, `rowFrom accepted ${JSON.stringify(garbage)}`);
      expect(`${rejected?.message}`.includes('apitests.item'), true,
        `the rejection does not name the table: ${rejected?.message}`);
    }
    // The legal shapes still pass: a values map, and null (a new row).
    expect(handler.rowFrom({sku: 'x'}).values['sku'], 'x', 'a values map was rejected');
    expect(handler.newRow().id == null, true, 'newRow() stopped working');

    const unknown = new DG.DomainObjectHandler('apitests.nosuch');
    for (const [member, action] of [['getProperties', () => unknown.getProperties()],
      ['capabilities', () => unknown.capabilities()]] as [string, () => Promise<any>][]) {
      const e = await thrown(action);
      expect(e instanceof DG.DomainValidationError, true,
        `${member} on an unknown table must reject with DomainValidationError, got ${e?.constructor?.name}: ${e?.message}`);
    }
  });
}, {owner: 'askalkin@datagrok.ai'});
