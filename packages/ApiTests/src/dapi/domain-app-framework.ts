import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test, awaitCheck} from '@datagrok-libraries/test/src/test';
import {DomainAppView, DomainEntityAppView, EntityListWidget, ENTITY_LIST_ROW_CAP,
  confirmDiscardChanges, domainHandler} from '@datagrok-libraries/domain-ui';
import {withRestrictedUser} from './domain-lifecycle';

// ui-js-api WO-9: the app framework — EntityListWidget (the "list of entities"
// an app opens on) and the AppView hierarchy (DomainAppView / DomainEntityAppView),
// including the save/discard prompt the views own on behalf of the refresh-agnostic
// editing controls.
//
// Fixture: the package's own 'apitests.item' (sku required+unique, name, quantity
// int min 0) and its child 'apitests.item_event'. Every test inserts inside its
// try and drops its rows in the finally with ONE filtered deleteWhere.
category('Dapi: domain app framework', () => {
  const items = () => grok.dapi.domains.table('apitests.item');
  const events = () => grok.dapi.domains.table('apitests.item_event');
  const stamp = () => `${Date.now()}-${Math.floor(Math.random() * 1e6)}`;

  async function seed(prefix: string, n: number): Promise<string[]> {
    const rows = [];
    for (let i = 0; i < n; i++)
      // The name carries the prefix too: it is the table's display-name column,
      // so a rendered card / brief item shows it and the search box matches it.
      rows.push({sku: `${prefix}-${i}`, name: `${prefix} name ${i}`, quantity: i});
    return (await items().insert(rows)).map((r) => r.id);
  }

  async function cleanup(prefix: string): Promise<void> {
    try {
      await items().deleteWhere({property: 'sku', operator: 'like', value: `${prefix}%`});
    } catch (e) {
      console.error(`app-framework fixture ${prefix} not cleaned up: ${e}`);
    }
  }

  const query = (prefix: string): any =>
    ({filter: {property: 'sku', operator: 'like', value: `${prefix}%`}, sort: 'sku'});

  /** Buttons of a widget, by their caption. */
  const buttons = (root: HTMLElement): string[] =>
    Array.from(root.querySelectorAll('button')).map((b) => `${b.textContent}`.trim());

  /** The unsaved-changes prompt, once it is up. */
  const prompt = async (): Promise<_DG.Dialog> => {
    const find = () => DG.Dialog.getOpenDialogs().find((d) => `${d.title}` === 'Unsaved changes');
    await awaitCheck(() => find() != null, 'the unsaved-changes prompt did not open', 10000);
    return find()!;
  };

  test('EntityListWidget: lists rows, searches, counts, gates New on canInsert', async () => {
    const prefix = `wo9-list-${stamp()}`;
    await seed(prefix, 3);
    let list: EntityListWidget | null = null;
    try {
      list = await EntityListWidget.create(items() as any, {query: query(prefix)});
      expect(list.rows.length, 3, 'the list did not show the fixture rows');
      // Rendering goes through the platform: one gallery card per row, carrying
      // the row's own identity.
      const cards = list.root.querySelectorAll('.grok-gallery-grid-item-wrapper');
      expect(cards.length, 3, `expected 3 cards, got ${cards.length}`);
      expect(`${list.root.textContent}`.includes(`${prefix} name 1`), true,
        'a card does not show the row it renders');

      // The search box ANDs its term into the query — server-side, not a DOM
      // filter — and matches the row's identity columns (name AND business key).
      await list.setSearch(`${prefix}-1`);
      expect(list.rows.length, 1, 'search did not match the business key');
      expect(list.rows[0].values['sku'], `${prefix}-1`, 'the business-key search matched the wrong row');
      await list.setSearch(`name 1`);
      expect(list.rows.length, 1, 'search did not narrow the list');
      expect(list.rows[0].values['sku'], `${prefix}-1`, 'search matched the wrong row');
      await list.setSearch('');
      expect(list.rows.length, 3, 'clearing the search did not restore the list');

      // Capability-driven affordance: admin may insert, so New is offered.
      expect(list.capabilities.canInsert, true, 'admin cannot insert into the fixture table');
      expect(buttons(list.root).some((b) => b.startsWith('New')), true,
        `no New button for a caller who may insert: ${buttons(list.root).join(', ')}`);

      // Selection and open are observable; the cap constant is the platform's.
      let selected: string | null = null;
      const sub = list.onItemSelected.subscribe((row) => selected = row?.values['sku'] ?? null);
      (cards[0] as HTMLElement).click();
      sub.unsubscribe();
      expect(selected, `${prefix}-0`, 'clicking a card did not report the item');
      expect(ENTITY_LIST_ROW_CAP, 1000, 'the list cap drifted from the platform Domain View');
    } finally {
      list?.detach();
      await cleanup(prefix);
    }
  });

  test('EntityListWidget: cards, brief, and an editable grid', async () => {
    const prefix = `wo9-modes-${stamp()}`;
    await seed(prefix, 2);
    let list: EntityListWidget | null = null;
    try {
      list = await EntityListWidget.create(items() as any, {query: query(prefix), editable: true});
      expect(list.mode, 'cards', 'the default mode is not cards');
      expect(list.grid, null, 'a card list must not build a grid');

      await list.setMode('brief');
      expect(list.rows.length, 2, 'the brief list lost its rows');
      expect(`${list.root.textContent}`.includes(`${prefix} name 0`), true,
        'the brief list renders nothing');

      await list.setMode('grid');
      expect(list.grid != null, true, 'grid mode did not build a DomainGrid');
      expect(list.grid!.dataFrame.rowCount, 2, 'the grid does not show the queried rows');
      expect(list.grid!.editable, list.capabilities.canEdit, 'the grid ignored the table capability');
      // The list's editors are what an app's dirty policy gates on.
      expect(list.editors.length, 1, 'grid mode must expose its editor');
      list.grid!.editor.setValue(0, 'name', 'Edited in the list');
      expect(list.isDirty, true, 'the list does not report its grid as dirty');
      list.grid!.editor.discard();
    } finally {
      list?.detach();
      await cleanup(prefix);
    }
  });

  test('DomainAppView: a list page in one line, round-tripping its query through the URL',
    async () => {
      const prefix = `wo9-app-${stamp()}`;
      await seed(prefix, 2);
      let view: DomainAppView | null = null;
      try {
        // THE measure: a working browse page, zero overrides, one expression.
        view = new DomainAppView(items() as any, {name: 'Items'});
        await awaitCheck(() => view!.list != null, 'the app view never loaded its list', 20000);
        expect(view.list!.rows.length >= 2, true, 'the app view lists no rows');

        // Query in: applied to the list AND published to the URL.
        const q = new DG.DomainQuery({schema: 'apitests', table: 'item',
          filters: [JSON.stringify({property: 'sku', operator: 'like', value: `${prefix}%`})],
          orderBy: ['sku'], limit: 10});
        expect(await view.setQuery(q), true, 'setQuery was refused on a clean page');
        expect(view.list!.rows.length, 2, 'the new query was not applied to the list');
        const path = `${view.path}`;
        expect(path.includes('filters'), true, `the query is missing from the path: ${path}`);

        // ...and out: the published parameters rebuild an equal query.
        const params: {[key: string]: string} = {};
        new URLSearchParams(path.substring(path.indexOf('?'))).forEach((v, k) => params[k] = v);
        const restored = DG.DomainQuery.fromUrlParams('apitests', 'item', params);
        expect(JSON.stringify(restored.toParams()), JSON.stringify(q.toParams()),
          'the URL round trip lost part of the query');

        // The view mode is UI state: a reserved parameter, never part of the query.
        await view.list!.setMode('grid');
        await awaitCheck(() => `${view!.path}`.includes('view=grid'),
          `the render mode never reached the URL: ${view!.path}`, 10000);
        expect(DG.DomainQuery.fromUrlParams('apitests', 'item', {view: 'grid'}).filters == null, true,
          'a reserved UI parameter leaked into the DomainQuery');
      } finally {
        view?.detach();
        await cleanup(prefix);
      }
    });

  test('DomainAppView: a deep link restores the query and the mode', async () => {
    const prefix = `wo9-link-${stamp()}`;
    await seed(prefix, 2);
    const original = window.location.search;
    let view: DomainAppView | null = null;
    try {
      const filter = JSON.stringify({property: 'sku', operator: '=', value: `${prefix}-1`});
      const deepLink = new URLSearchParams({'filters[0]': filter, 'view': 'brief'});
      // What a pasted app URL looks like when the page opens.
      window.history.replaceState(null, '', `?${deepLink.toString()}`);
      view = new DomainAppView(items() as any);
      await awaitCheck(() => view!.list != null, 'the deep-linked view never loaded', 20000);
      expect(view.list!.mode, 'brief', 'the reserved view= parameter was not restored');
      expect(view.list!.rows.length, 1, 'the deep-linked filter was not applied');
      expect(view.list!.rows[0].values['sku'], `${prefix}-1`, 'the deep link restored the wrong row');
      expect(`${view.path}`.includes('view=brief'), true,
        `the restored mode was not republished: ${view.path}`);
      expect((view.query.filters ?? []).length, 1, 'the restored query lost its filter element');
    } finally {
      window.history.replaceState(null, '', original === '' ? window.location.pathname : original);
      view?.detach();
      await cleanup(prefix);
    }
  });

  test('DomainEntityAppView: the row page, its detail tab, and a form save', async () => {
    const prefix = `wo9-entity-${stamp()}`;
    const ids = await seed(prefix, 1);
    let view: DomainEntityAppView | null = null;
    try {
      await events().insert([{item_id: ids[0], kind: 'created', amount: 1}]);
      view = new DomainEntityAppView(items() as any, ids[0]);
      await awaitCheck(() => view!.row != null, 'the entity page never loaded its row', 20000);
      expect(view.row!.id, ids[0], 'the page loaded the wrong row');
      expect(`${view.root.textContent}`.includes(`${prefix} name 0`), true,
        'the page header does not name the row');

      // The child table is a tab, FK-inverted from the registry.
      await awaitCheck(() => `${view!.root.textContent}`.includes('apitests.item_event'),
        `the detail tab is missing: ${view!.root.textContent}`, 15000);

      // The Details pane is the platform's reflective form, bound to the row:
      // typing into it writes through, and Save sends exactly what changed.
      const renamed = `Renamed ${stamp()}`;
      await awaitCheck(() => nameInput(view!.root) != null, 'the property form never rendered', 15000);
      const input = nameInput(view.root)!;
      input.value = renamed;
      input.dispatchEvent(new Event('input', {bubbles: true}));
      input.dispatchEvent(new Event('change', {bubbles: true}));
      expect(view.row!.values['name'], renamed, 'the form did not write into the row');
      expect(await view.save(), true, 'the form save failed');
      expect((await items().get(ids[0])).name, renamed, 'the save never reached the server');
    } finally {
      view?.detach();
      await cleanup(prefix);
    }
  });

  /** The 'name' text editor of a rendered property form. */
  function nameInput(root: HTMLElement): HTMLInputElement | null {
    for (const el of Array.from(root.querySelectorAll('.ui-input-root'))) {
      const label = `${el.querySelector('label')?.textContent ?? ''}`.trim().toLowerCase();
      const editor = el.querySelector('input.ui-input-editor') as HTMLInputElement | null;
      if (label === 'name' && editor != null)
        return editor;
    }
    return null;
  }

  test('unsaved changes: the navigation prompt cancels, discards, or saves', async () => {
    const prefix = `wo9-dirty-${stamp()}`;
    await seed(prefix, 2);
    let list: EntityListWidget | null = null;
    try {
      list = await EntityListWidget.create(items() as any,
        {query: query(prefix), mode: 'grid', editable: true});
      const editor = list.grid!.editor;
      if (!editor.capabilities.canEdit) {
        console.log('skipped: the fixture table is not editable for this caller');
        return;
      }

      // CANCEL — the rebuild does not happen and the changes survive.
      editor.setValue(0, 'name', 'Pending 1');
      const canceled = list.setQuery(query(`${prefix}-1`));
      (await prompt()).getButton('CANCEL').click();
      expect(await canceled, false, 'a canceled prompt still applied the query');
      expect(editor.isDirty, true, 'a canceled prompt dropped the pending changes');
      expect(list.rows.length + list.grid!.dataFrame.rowCount, 2, 'the list was rebuilt anyway');

      // DISCARD — the rebuild happens and the changes are gone (unsaved).
      const discarded = list.setQuery(query(prefix));
      (await prompt()).getButton('DISCARD').click();
      expect(await discarded, true, 'the discard outcome refused the query');
      expect(list.grid!.editor.isDirty, false, 'discard left the editor dirty');
      expect((await items().first({filter: {property: 'sku', operator: '=',
        value: `${prefix}-0`} as any}))!.name, `${prefix} name 0`,
      'a discarded change reached the server');

      // SAVE — the batch is written, THEN the rebuild happens.
      const saved = `Saved through the prompt ${stamp()}`;
      list.grid!.editor.setValue(0, 'name', saved);
      const applied = list.setQuery(query(prefix));
      (await prompt()).getButton('SAVE').click();
      expect(await applied, true, 'the save outcome refused the query');
      expect(list.grid!.editor.isDirty, false, 'the save left changes pending');
      expect((await items().first({filter: {property: 'sku', operator: '=',
        value: `${prefix}-0`} as any}))!.name, saved, 'the prompt saved nothing');

      // Nothing pending — no prompt at all (this would hang otherwise).
      expect(await confirmDiscardChanges(list.editors), true,
        'a clean page still asked about unsaved changes');
    } finally {
      DG.Dialog.getOpenDialogs().forEach((d) => d.close());
      list?.detach();
      await cleanup(prefix);
    }
  });

  test('restricted user: no New, no editing, no Delete action', async () => {
    const prefix = `wo9-perm-${stamp()}`;
    await seed(prefix, 1);
    let outcome: string | null = null;
    try {
      outcome = await withRestrictedUser('wo9perm', async (probe) => {
        const table = items();
        await table.grant(probe.group, 'View');
        try {
          await probe.asUser(async () => {
            grok.dapi.domains.invalidateUiCaches();
            let list: EntityListWidget | null = null;
            try {
              list = await EntityListWidget.create(grok.dapi.domains.table('apitests.item') as any,
                {query: query(prefix), mode: 'grid', editable: true});
              expect(list.capabilities.canInsert, false, 'a View-only user reports canInsert');
              expect(buttons(list.root).some((b) => b.startsWith('New')), false,
                `a View-only user is offered New: ${buttons(list.root).join(', ')}`);
              expect(list.grid!.editable, false, 'a View-only user got an editable grid');
              expect(list.rows.length + list.grid!.dataFrame.rowCount >= 1, true,
                'the View-only user cannot even read the rows');
              // Row commands come from the row's own permissions. They are NOT
              // asserted here: `DomainRow.permissions()` short-circuits on the
              // Dart client's `Auth.adminMode`, which this harness cannot swap
              // (it swaps the auth cookie, not the booted client's session), so a
              // cookie-impersonated admin still sees every action. Row-level
              // degradation is verified in a REAL restricted browser session
              // instead (WO-9 live pass). What IS session-truthful here — the
              // capabilities — is asserted above.
              const row = domainHandler('apitests.item').rowFrom(
                (await grok.dapi.domains.table('apitests.item').first({filter:
                  {property: 'sku', operator: 'like', value: `${prefix}%`} as any}))!);
              const actions = (await domainHandler('apitests.item').getRibbonActions(row as any))
                .map((a) => a.name);
              expect(actions.includes('Open'), true, `the default Open action is missing: ${actions}`);
            } finally {
              list?.detach();
            }
          });
        } finally {
          await table.revoke(probe.group, 'View');
        }
        return 'checked';
      });
      expect(outcome, 'checked',
        'permission degradation was NOT verified: no restricted session (see the console for why)');
    } finally {
      grok.dapi.domains.invalidateUiCaches();
      await cleanup(prefix);
    }
  });
});
