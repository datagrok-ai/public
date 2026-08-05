import type * as _grok from 'datagrok-api/grok';
import type * as _ui from 'datagrok-api/ui';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok;
declare let ui: typeof _ui;
declare let DG: typeof _DG;

import {category, expect, test, awaitCheck} from '@datagrok-libraries/test/src/test';
import {DomainForm, domainHandler, domains} from '@datagrok-libraries/domain-ui';

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
        // INSERT form is pending from birth: its unsaved row is a real pending
        // change, exactly as an 'Add row' in the grid is.
        expect(view.editors.length, 1, 'the form editor did not roll up to the page');
        expect(view.isDirty, true, 'the page does not see the pending new row');
        form.props['sku'] = `${prefix}-0`;
        expect(view.isDirty, true, 'the page does not see the form as dirty');
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
