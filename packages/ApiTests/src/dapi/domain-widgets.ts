import type * as _grok from 'datagrok-api/grok';
import type * as _ui from 'datagrok-api/ui';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok;
declare let ui: typeof _ui;
declare let DG: typeof _DG;

import {category, expect, test, awaitCheck} from '@datagrok-libraries/test/src/test';
import {DomainAppView, DomainForm, DomainFrameEditor, DomainGrid,
  domainHandler, domains} from '@datagrok-libraries/domain-ui';

// ui-js-api WO-8R: the `domains` facade — the prefetched table handle, `DomainForm`
// on a one-row frame editor, the view/dialog hosts, and the MACHINE SURFACE
// (`props` / `getWidgetStatus()` / `getFunctions()`), which is what these tests
// drive: no DOM scraping, no direct editor calls.
//
// The tests are the UI-FACADE.md §1 examples, run against the package's own
// 'apitests.item' (sku required + unique, name, quantity int min 0) and its child
// 'apitests.item_event' instead of the doc's 'grit.issue'. Every test inserts inside
// its try and drops its rows in the finally with ONE filtered deleteWhere.
category('Dapi: domain widgets (smoke)', () => {
  const items = () => domains.table('apitests.item');
  const events = () => domains.table('apitests.item_event');
  const stamp = () => `${Date.now()}-${Math.floor(Math.random() * 1e6)}`;

  async function cleanup(prefix: string): Promise<void> {
    try {
      await grok.dapi.domains.table('apitests.item')
        .deleteWhere({property: 'sku', operator: 'like', value: `${prefix}%`});
    } catch (e) {
      console.error(`domain-widgets fixture ${prefix} not cleaned up: ${e}`);
    }
  }

  // The dialog that opened after the action started — the user's next gesture.
  async function openedDialog(before: _DG.Dialog[]): Promise<_DG.Dialog> {
    const fresh = () => DG.Dialog.getOpenDialogs().filter((d) => !before.includes(d));
    await awaitCheck(() => fresh().length > 0, 'the dialog never opened', 15000);
    return fresh()[0];
  }

  // Captions of everything on the view's ribbon.
  const ribbon = (view: any): string[] => (view.getRibbonPanels() as HTMLElement[][])
    .reduce((all: HTMLElement[], panel) => all.concat(panel), [])
    .map((e) => `${e.textContent}`.trim());

  test('UI-FACADE 1.2: an insert dialog saves, a dismissed one discards silently', async () => {
    const prefix = `wo8r-dlg-${stamp()}`;
    try {
      // §1.2, with the fixture standing in for `project`. The example awaits the
      // dialog; a test has to start it, answer it, and only then await.
      const issues = await items();
      const opened = DG.Dialog.getOpenDialogs();
      const saved = issues.formDialog({values: {sku: `${prefix}-0`, name: `${prefix} name`, quantity: 1}});
      (await openedDialog(opened)).getButton('OK').click();
      expect(await saved, true, 'the insert dialog did not report a save');
      const rows = await issues.query({filter: DG.cond('sku', '=', `${prefix}-0`)});
      expect(rows.length, 1, 'the dialog saved nothing');
      expect(rows[0]['quantity'], 1, 'the prefilled value did not reach the server');

      // Editing an existing row: `await issues.formDialog({row})` — dismissed here,
      // which is a SILENT discard resolving false (the dialog ruling).
      const row = domainHandler('apitests.item').rowFrom(rows[0]);
      const before = DG.Dialog.getOpenDialogs();
      const editing = issues.formDialog({row: row});
      const dialog = await openedDialog(before);
      dialog.close();
      expect(await editing, false, 'a dismissed dialog did not resolve false');
      expect((await issues.get(rows[0]['id']))['name'], `${prefix} name`,
        'a dismissed dialog changed the row');
    } finally {
      DG.Dialog.getOpenDialogs().forEach((d) => d.close());
      await cleanup(prefix);
    }
  });

  test('UI-FACADE 1.3: a routed form page carries the form actions and its dirty state',
    async () => {
      const prefix = `wo8r-page-${stamp()}`;
      let view: any = null;
      try {
        // §1.3 verbatim, modulo the fixture and the `as any` the typing of
        // `grok.shell.preview` forces (it is `View | null`, and every AppView is a
        // `ViewBase` — recorded as an API defect of the example).
        const issues = await items();
        const form = issues.form({values: {name: `${prefix} seed`}});
        view = domains.formView(form);
        grok.shell.preview = view as any;
        await form.ready;

        expect(view.widgets.length, 1, 'the page does not host the form');
        const captions = ribbon(view);
        expect(captions.includes('Save'), true, `Save is not on the ribbon: ${captions.join(', ')}`);
        expect(captions.includes('Reset'), true, `Reset is not on the ribbon: ${captions.join(', ')}`);

        // The page answers for the form's editor — one dirty state, one gate. An
        // insert form is PRISTINE UNTIL TOUCHED: its unsaved row (prefilled or
        // not) is not a pending change until something writes into it, so an
        // untouched New page closes without a prompt.
        expect(view.editors.length, 1, 'the form editor did not roll up to the page');
        expect(view.isDirty, false, 'an untouched New form armed the unsaved-changes gate');
        form.props['sku'] = `${prefix}-0`;
        expect(view.isDirty, true, 'the first write did not make the form dirty');
        expect(`${view.statusBarPanels[0].textContent}`.includes('unsaved change'), true,
          `the pending count is not on the status bar: ${view.statusBarPanels[0].textContent}`);

        // The page's status nests the form's inputs under the child's name.
        const status = view.getWidgetStatus() as _DG.IWidgetStatus;
        expect((status.inputs ?? []).some((i) => i.name === 'DomainForm[0].sku'), true,
          `the page did not nest the form inputs: ${(status.inputs ?? []).map((i) => i.name)}`);
        form.discard();
      } finally {
        view?.close();
        await cleanup(prefix);
      }
    });

  test('UI-FACADE 1.6c/1.6d: one replaced input, one extra validator', async () => {
    const prefix = `wo8r-custom-${stamp()}`;
    let form: DomainForm | null = null;
    try {
      const issues = await items();
      form = issues.form();
      // §1.6c — the doc replaces 'description' on grit.issue; the fixture's free-text
      // column is 'name'.
      form.replaceInput('name', (p) => ui.input.textArea(p.caption));
      expect(form.getInput('name')!.root.querySelector('textarea') != null, true,
        'replaceInput did not put the custom editor on the form');

      // §1.6d verbatim.
      form.addValidator('name', (s) => s.trim().length < 5 ? 'At least 5 characters' : null);
      await form.ready;
      form.props['name'] = 'x';
      const short = form.getWidgetStatus().inputs!.find((i) => i.name === 'name')!;
      expect(`${short.error}`.includes('At least 5 characters'), true,
        `the validator did not mark the field: ${short.error}`);
      expect(short.valid, false, 'the marked input still reports itself valid');

      form.props['name'] = `${prefix} long enough`;
      expect(form.getWidgetStatus().inputs!.find((i) => i.name === 'name')!.error == null, true,
        'the validator did not clear');
      // The replaced input received the value through the editor, like a keystroke.
      expect(form.getInput('name')!.value, `${prefix} long enough`,
        'the custom input was not synced from the row');
    } finally {
      form?.detach();
      await cleanup(prefix);
    }
  });

  test('UI-FACADE 1.11: the form is discoverable, describable and invocable', async () => {
    const prefix = `wo8r-machine-${stamp()}`;
    let form: DomainForm | null = null;
    try {
      const issues = await items();
      form = issues.form();
      await form.ready;

      // Discovery is the platform's own registry of widgets.
      const all = DG.Widget.getAll();
      expect(all.includes(form as any), true,
        'DG.Widget.getAll() does not report the form');
      expect(all.find((x) => x.type === 'DomainForm') != null, true,
        'no widget reports itself as a DomainForm');
      // §1.11's `find` takes the FIRST registered DomainForm. A JS widget stays in
      // the platform's widget map until its ROOT is killed with the DOM, so a run
      // that built forms in earlier tests has to address its own.
      const w = all.filter((x) => x.type === 'DomainForm').find((x) => x === (form as any))!;

      const status = w.getWidgetStatus();
      expect(`${status.description}`.includes('apitests.item'), true,
        `the description does not name the table: ${status.description}`);
      expect(Object.keys(status.parts).includes('sku'), true,
        `the input roots are not named parts: ${Object.keys(status.parts)}`);
      const sku = status.inputs!.find((i) => i.name === 'sku')!;
      expect(sku.required, true, 'a required column does not report itself required');
      expect(sku.type, 'string', 'the input type is not the registry type');

      // A property write takes the single-writer path a keystroke takes.
      w.props['quantity'] = 7;
      expect(form.getValue('quantity'), 7, 'a props write did not reach the row');
      expect(form.editor!.isChanged(0, 'quantity'), true,
        'a props write did not go through the frame editor');
      expect(w.getWidgetStatus().inputs!.find((i) => i.name === 'quantity')!.value, 7,
        'the status does not report the value that was just written');

      // The actions are REAL platform Funcs.
      const save = w.getFunctions().find((f) => f.name === 'Save')!;
      expect(save != null, true, 'the form offers no Save function');
      expect(save.inputs.map((p) => p.name).includes('widget'), true,
        'the Save function does not take the widget it acts on');
    } finally {
      form?.detach();
      await cleanup(prefix);
    }
  });

  // ─────────────────────── WO-9R ─────────────────────────

  test('an insert form is pristine until touched (the gate arms on the first write)',
    async () => {
      const prefix = `wo9r-pristine-${stamp()}`;
      let form: DomainForm | null = null;
      try {
        // Prefilled BY THE DEVELOPER, not by the user: nothing has been written
        // into the form, so there is nothing to prompt about when it closes.
        const table = await items();
        form = table.form({values: {name: `${prefix} seeded`}});
        await form.ready;
        expect(form.editor!.stateOf(0), 'new', 'the insert form has no unsaved row');
        expect(form.editor!.changeCount, 0, 'a prefilled but untouched form counts as a change');
        expect(form.isDirty, false, 'an untouched New form is dirty');
        expect(form.getValue('name'), `${prefix} seeded`, 'the seeded value is not in the row');

        // The first write — user or machine, same path — arms it.
        form.props['sku'] = `${prefix}-0`;
        expect(form.isDirty, true, 'the first write did not arm the gate');
        expect(form.editor!.changeCount, 1, 'a touched new row counts more than once');

        // ...and a reset puts it back to pristine, exactly as it opened.
        form.reset();
        expect(form.isDirty, false, 'a reset form stayed dirty');
        expect(form.getValue('name'), `${prefix} seeded`, 'the reset dropped the seeded value');
        expect(form.getValue('sku'), null, 'the reset kept the written value');

        // An 'Add row' in a GRID is a user gesture: pending immediately, as before.
        const grid = table.grid({query: {filter: DG.cond('sku', 'like', `${prefix}%`)}});
        try {
          await grid.ready;
          grid.addRow();
          expect(grid.editor.changeCount, 1, 'a row the user added is not pending');
        } finally {
          grid.detach();
        }
      } finally {
        form?.detach();
        await cleanup(prefix);
      }
    });

  test('UI-FACADE 1.5: an editable grid over a filtered subset, through its status',
    async () => {
      const prefix = `wo9r-grid-${stamp()}`;
      let grid: DomainGrid | null = null;
      try {
        await grok.dapi.domains.table('apitests.item')
          .insert([{sku: `${prefix}-0`, name: `${prefix} grid item`, quantity: 2}]);
        // §1.5 verbatim, modulo the fixture: the factory is SYNCHRONOUS on the
        // handle, and `ready` is the one place the loading boundary shows.
        const table = await items();
        grid = table.grid({query: {filter: DG.cond('sku', 'like', `${prefix}%`)},
          defaults: {sku: `${prefix}-default`}});
        grok.shell.newView('My items', [grid.root]);
        await grid.ready;
        expect(grid.dataFrame.rowCount, 1, 'the grid did not run the query it was given');

        // A grid's machine surface is status + functions + the dataFrame — no
        // per-cell inputs.
        let status = grid.getWidgetStatus();
        expect(status.inputs === undefined, true, 'a grid reported per-cell inputs');
        expect(`${status.description}`.includes('0 unsaved changes'), true,
          `the grid does not report its pending count: ${status.description}`);
        expect(status.error, null, `a freshly loaded grid reports an error: ${status.error}`);
        expect(Object.keys(status.parts).includes('grid'), true,
          `the grid is not a named part: ${Object.keys(status.parts)}`);

        // An invalid edit flips both counts and the error — read live.
        grid.editor.setValue(0, 'quantity', -3);
        status = grid.getWidgetStatus();
        expect(`${status.description}`.includes('1 unsaved change,'), true,
          `the pending count did not flip: ${status.description}`);
        expect(`${status.error}`.includes('quantity'), true,
          `the invalid cell is not reported: ${status.error}`);

        // ...and the Save function writes the batch (readback), through the same
        // Func the AI would call.
        grid.editor.setValue(0, 'quantity', 5);
        expect(grid.getWidgetStatus().error, null,
          `the grid stayed invalid: ${grid.getWidgetStatus().error}`);
        const save = grid.getFunctions().find((f) => f.name === 'Save')!;
        expect(save != null, true, `no Save function on an editable grid: ` +
          `${grid.getFunctions().map((f) => f.name)}`);
        await save.apply({widget: grid});
        expect((await table.query({filter: DG.cond('sku', '=', `${prefix}-0`)}))[0]['quantity'], 5,
          'the Save function did not write the grid batch');
        expect(grid.getWidgetStatus().error, null, 'the saved grid still reports an error');
      } finally {
        grid?.detach();
        await cleanup(prefix);
      }
    });

  test('UI-FACADE 1.8/1.9: a composed page — one dirty state, one gate, one status',
    async () => {
      const prefix = `wo9r-page-${stamp()}`;
      let view: any = null;
      try {
        // §1.8: a form and a grid on ONE page. §1.9's chart is left out — a
        // DG.Viewer is already a widget, and the point being pinned here is the
        // aggregation, not the chart.
        const table = await items();
        const form = table.form({values: {name: `${prefix} composed`}});
        const grid = table.grid({query: {filter: DG.cond('sku', 'like', `${prefix}%`)}});
        view = domains.view([form, grid], {name: 'Composed'});
        grok.shell.addView(view);
        await Promise.all([form.ready, grid.ready]);

        expect(view.widgets.length, 2, 'the page does not host both widgets');
        expect(view.editors.length, 2, 'the page did not roll both editors up');
        expect(view.isDirty, false, 'an untouched composed page is dirty');
        form.props['sku'] = `${prefix}-0`;
        expect(view.isDirty, true, 'the page does not see the form as dirty');

        // ONE status, with the children nested under their own names, and the
        // functions of both merged.
        const status = view.getWidgetStatus() as _DG.IWidgetStatus;
        expect((status.inputs ?? []).some((i) => i.name === 'DomainForm[0].sku'), true,
          `the form inputs are not nested: ${(status.inputs ?? []).map((i) => i.name)}`);
        expect(Object.keys(status.parts).some((p) => p.startsWith('DomainGrid[1].')), true,
          `the grid parts are not nested: ${Object.keys(status.parts)}`);
        expect(`${status.description}`.includes('hosts 2 widgets'), true,
          `the page does not describe what it hosts: ${status.description}`);
        const names = (view.getFunctions() as _DG.Func[]).map((f) => f.name);
        expect(names.includes('Save') && names.includes('AddRow'), true,
          `the page did not merge its children's functions: ${names}`);
        view.discard();
      } finally {
        view?.close();
        await cleanup(prefix);
      }
    });

  test('routing: one router per app — back onto ?entity= reuses the row page', async () => {
    const prefix = `wo9r-route-${stamp()}`;
    const base = `/apps/ApiTests/DomainWidgets-${stamp()}`;
    const original = window.location.search;
    let view: DomainAppView | null = null;
    try {
      const [row] = await grok.dapi.domains.table('apitests.item')
        .insert([{sku: `${prefix}-0`, name: `${prefix} routed`, quantity: 1}]);
      const table = await items();
      view = table.listView({name: 'Routed', query: {filter: DG.cond('sku', 'like', `${prefix}%`)}});
      grok.shell.addView(view);
      await awaitCheck(() => view!.list != null, 'the list page never loaded', 20000);
      // A page opened outside an app call has no URL prefix of its own; basePath
      // is what the platform composes one from.
      view.basePath = base;

      // Only the ROOT page claims the app prefix — a row page never does, which
      // is what makes back / forward independent of the order of the views.
      const page = view.open(row.id)!;
      await awaitCheck(() => page.row != null, 'the row page never loaded', 20000);
      expect(view.acceptsPath(base), true, 'the app root does not accept its own prefix');
      expect(page.acceptsPath(base), false, 'a row page claims the app prefix');

      // Three back/forward round trips onto ?entity= — the SAME page every time,
      // and no new view is opened.
      window.history.replaceState(null, '', `?entity=${row.id}`);
      const count = () => Array.from(grok.shell.views).length;
      const views = count();
      for (let i = 0; i < 3; i++) {
        view.handlePath(base);
        await DG.delay(50);
        expect(view.open(row.id) === page, true, `round trip ${i + 1} spawned another row page`);
      }
      expect(count(), views, 'back / forward onto ?entity= stacked duplicate row pages');
    } finally {
      window.history.replaceState(null, '', original === '' ? window.location.pathname : original);
      view?.detach();
      view?.close();
      await cleanup(prefix);
    }
  });

  test('WO-9F2: stale edit snapshots, cancelled actions, and the cap-banner box', async () => {
    const prefix = `wo9r-fixes-${stamp()}`;
    let editor: DomainFrameEditor | null = null;
    let host: HTMLElement | null = null;
    try {
      const client = grok.dapi.domains.table('apitests.item');
      await client.insert([0, 1, 2].map((i) =>
        ({sku: `${prefix}-${i}`, name: `${prefix} name ${i}`, quantity: i})));
      editor = await DomainFrameEditor.create(client as any,
        {query: {filter: {property: 'sku', operator: 'like', value: `${prefix}%`} as any,
          sort: 'sku'}});
      expect(editor.dataFrame.rowCount, 3, 'the fixture rows did not load');

      // An EXTERNAL row removal shifts every index below it — including the keys
      // of the pre-edit snapshots. A stale one made `beginEdit` early-return and
      // the next `commitEdit` record ANOTHER row's original value.
      editor.beginEdit(1);
      editor.dataFrame.rows.removeAt(0, 1);
      editor.beginEdit(1);
      editor.dataFrame.set('name', 1, 'Edited after the shift');
      editor.commitEdit(1, 'name');
      expect(editor.changesOf(1)['name'], `${prefix} name 2`,
        'the stale snapshot recorded another row\'s original value');

      // A CANCELLED action must not reload the row page (which would pop the
      // unsaved-changes prompt right behind the dialog the user dismissed).
      const rows = await client.query({filter: {property: 'sku', operator: 'like',
        value: `${prefix}%`} as any, limit: 1}) as any[];
      const table = await items();
      const page = table.entityView(rows[0]['id'], {
        actions: (_row, defaults) => defaults.concat([
          {name: 'Cancelled', icon: 'pencil', run: () => false},
          {name: 'Applied', icon: 'pencil', changesRow: true, run: () => true}]),
      });
      try {
        await awaitCheck(() => page.form != null, 'the row page never loaded', 20000);
        const form = page.form;
        const actions = await domainHandler('apitests.item')
          .getRibbonActions(page.row as any);
        expect(actions.some((a) => a.name === 'Open'), true,
          'the fixture has no Open action to drop');
        await page.runAction({name: 'Cancelled', icon: 'pencil', run: () => false});
        expect(page.form === form, true, 'a cancelled action reloaded the page');
        await page.runAction({name: 'Applied', icon: 'pencil', changesRow: true,
          run: () => true});
        await awaitCheck(() => page.form !== form, 'changesRow did not reload the page', 20000);
      } finally {
        page.detach();
      }

      // The cap-banner wrapper of a child pane: without its own rule the box
      // under it takes ui.css's default 400x300, whatever the pane is wide.
      host = ui.div([], 'ui-box');
      host.style.width = '900px';
      host.style.height = '400px';
      document.body.appendChild(host);
      const inner = ui.div([]);
      host.appendChild(ui.divV([ui.divText('Showing the first N of M rows'), ui.box(inner)],
        'domain-ui-detail-pane'));
      await DG.delay(50);
      expect((inner.parentElement as HTMLElement).getBoundingClientRect().width > 400, true,
        'the cap-banner branch collapses its grid to the default 400px box: ' +
        `${(inner.parentElement as HTMLElement).getBoundingClientRect().width}`);
    } finally {
      host?.remove();
      editor?.detach();
      await cleanup(prefix);
    }
  });

  test('the machine surface drives a form end to end (props, status, functions)', async () => {
    const prefix = `wo8r-e2e-${stamp()}`;
    let form: DomainForm | null = null;
    let child: DomainForm | null = null;
    try {
      const table = await items();
      form = table.form();
      await form.ready;

      // Invalid — the registry's own constraint, reported through the status.
      form.props['quantity'] = -1;
      const quantity = form.getWidgetStatus().inputs!.find((i) => i.name === 'quantity')!;
      expect(`${quantity.error}`.includes('less than'), true,
        `an out-of-range value was not reported: ${quantity.error}`);
      expect(form.getWidgetStatus().error != null, true, 'the form reports itself valid');

      // Valid — and saved through the function, not through the editor.
      form.props['sku'] = `${prefix}-0`;
      form.props['name'] = `${prefix} machine item`;
      form.props['quantity'] = 3;
      expect(form.getWidgetStatus().error, null,
        `the form is still invalid: ${form.getWidgetStatus().error}`);
      await form.getFunctions().find((f) => f.name === 'Save')!.apply({widget: form});
      const saved = await table.query({filter: DG.cond('sku', '=', `${prefix}-0`)});
      expect(saved.length, 1, 'the Save function saved nothing');
      expect(saved[0]['quantity'], 3, 'the saved row lost a value');

      // A reference written as DISPLAY TEXT resolves through the same suggestion
      // source the picker's type-ahead uses.
      const eventTable = await events();
      child = eventTable.form();
      await child.ready;
      child.props['item_id'] = `${prefix} machine item`;
      child.props['kind'] = 'created';
      await child.settle();
      expect(child.getWidgetStatus().error, null,
        `the reference did not resolve: ${child.getWidgetStatus().error}`);
      await child.getFunctions().find((f) => f.name === 'Save')!.apply({widget: child});
      const written = await eventTable.query({filter: DG.cond('item_id', '=', saved[0]['id'])});
      expect(written.length, 1, 'the reference resolved from its display name did not save');

      // An action the form does not offer is simply absent — no throw, no guessing.
      expect(form.getFunctions().find((f) => f.name === 'Escalate') == null, true,
        'the form claims an action it does not have');

      // Ambiguity is a validation error, never a throw.
      child.props['item_id'] = 'no such item at all';
      await child.settle();
      expect(`${child.getWidgetStatus().inputs!.find((i) => i.name === 'item_id')!.error}`
        .includes('No apitests.item matches'), true,
      'an unresolvable reference was not reported as a validation error');
    } finally {
      child?.detach();
      form?.detach();
      await cleanup(prefix);
    }
  });
});
