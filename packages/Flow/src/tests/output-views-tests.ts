/** Output-view tabs (Spotfire-style internal pages): the status-bar strip
 *  renders Canvas + one tab per table output (Table Output, dataframe-typed
 *  Value Output, table-carrying SetVar terminal); an empty tab shows the
 *  "run the flow" message; a value materializes a detached DG.TableView
 *  lazily on activation (`_onAdded` against a live, laid-out pane); a re-run
 *  refreshes the SAME TableView in place; ribbon/toolbox swap to the active
 *  tab's (MultiView pattern) and restore on Canvas; layouts persist in the
 *  `.ffjson` keyed by paramName and survive node-id remapping across a
 *  save → load round-trip. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/utils/src/test';

import {FuncFlowView} from '../funcflow-view';
import {FlowEditor} from '../rete/flow-editor';
import {ensureFuncNodeType} from '../rete/node-factory';
import {parseFlowBody} from '../serialization/flow-script-format';
import {FuncFlowDocument} from '../serialization/flow-schema';
import {until, addNode} from './test-utils';

interface ViewHarness {
  view: FuncFlowView;
  host: HTMLElement;
  flow: FlowEditor;
}

/** A FuncFlowView mounted offscreen but laid out with a real size — the tab
 *  panes must be non-zero-sized when a TableView `_onAdded`s into them. */
async function makeView(): Promise<ViewHarness> {
  const view = new FuncFlowView();
  const host = ui.div([view.root], {style: {
    width: '1100px', height: '700px', position: 'absolute', left: '-10000px',
  }});
  document.body.appendChild(host);
  const ok = await until(() => (view as unknown as {flow?: FlowEditor}).flow != null, 10000);
  if (!ok) throw new Error('editor did not initialize');
  return {view, host, flow: (view as unknown as {flow: FlowEditor}).flow};
}

function destroyView(h: ViewHarness): void {
  try {
    h.view.outputViews?.destroy();
  } catch {/* already gone */}
  try {
    (h.view as unknown as {flow?: FlowEditor}).flow?.destroy();
  } catch {/* already gone */}
  h.host.remove();
}

const chips = (h: ViewHarness): HTMLElement[] =>
  Array.from(h.view.root.querySelectorAll('[data-testid="ff-view-tabs"] .ff-view-tab'));
const outputChips = (h: ViewHarness): HTMLElement[] =>
  chips(h).filter((c) => c.dataset.nodeId != null);
const chipByParam = (h: ViewHarness, param: string): HTMLElement | undefined =>
  outputChips(h).find((c) => c.dataset.param === param);

function numericDf(name: string, rows: number): DG.DataFrame {
  const xs = Array.from({length: rows}, (_, i) => i);
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList('int', 'x', xs),
    DG.Column.fromList('double', 'y', xs.map((v) => v * 1.5)),
  ]);
  df.name = name;
  return df;
}

category('Flow: output views', () => {
  test('the strip renders Canvas plus one tab per table output; scalar outputs stay out', async () => {
    const h = await makeView();
    try {
      expect(chips(h).length, 1, 'an empty flow shows just the Canvas tab');
      expect(chips(h)[0].dataset.param, 'canvas');

      const tableOut = await addNode(h.flow, 'Outputs/Table Output', 100, 0);
      tableOut.properties['paramName'] = 'result';
      h.flow.notifyNodeParamsChanged(tableOut.id);
      expect(outputChips(h).length, 1, 'a Table Output gets a tab');

      const valueOut = await addNode(h.flow, 'Outputs/Value Output', 100, 100);
      valueOut.properties['paramName'] = 'vals';
      h.flow.notifyNodeParamsChanged(valueOut.id);
      expect(outputChips(h).length, 1, 'a double-typed Value Output stays out');

      valueOut.properties['outputType'] = 'dataframe';
      h.flow.notifyNodeParamsChanged(valueOut.id);
      expect(outputChips(h).length, 2, 'a dataframe-typed Value Output gets a tab');
      expect(!!chipByParam(h, 'vals'), true, 'the tab is keyed/labeled by paramName');
      expect(chipByParam(h, 'vals')!.dataset.state, 'empty', 'no value yet → empty state');
    } finally {
      destroyView(h);
    }
  }, {timeout: 60000});

  test('node add / remove / rename updates the tabs; removal falls back to Canvas', async () => {
    const h = await makeView();
    try {
      const out = await addNode(h.flow, 'Outputs/Table Output', 100, 0);
      out.properties['paramName'] = 'first';
      h.flow.notifyNodeParamsChanged(out.id);
      expect(!!chipByParam(h, 'first'), true);

      out.properties['paramName'] = 'renamed';
      h.flow.notifyNodeParamsChanged(out.id);
      expect(!!chipByParam(h, 'renamed'), true, 'rename tracks the paramName');
      expect(chipByParam(h, 'renamed')!.textContent!.includes('renamed'), true, 'label follows');

      h.view.outputViews.activate(out.id);
      expect(h.view.outputViews.activeKey, out.id, 'the tab activated');
      await h.flow.removeNode(out.id);
      expect(outputChips(h).length, 0, 'removing the node removes its tab');
      expect(h.view.outputViews.activeKey, 'canvas', 'an active removed tab falls back to Canvas');
    } finally {
      destroyView(h);
    }
  }, {timeout: 60000});

  test('an empty tab shows the run message; Canvas re-activation restores the canvas', async () => {
    const h = await makeView();
    try {
      const out = await addNode(h.flow, 'Outputs/Table Output', 100, 0);
      const canvasPane = h.view.root.querySelector('[data-testid="ff-root"]') as HTMLElement;
      const tab = h.view.outputViews.getTab(out.id)!;

      h.view.outputViews.activate(out.id);
      expect(tab.pane.style.display, 'flex', 'the tab pane is shown');
      expect(canvasPane.style.display, 'none', 'the canvas is hidden');
      const empty = tab.pane.querySelector('[data-testid="ff-output-view-empty"]') as HTMLElement;
      expect(!!empty, true, 'the empty state is rendered');
      expect((empty.textContent ?? '').includes('Run the flow'), true, 'it says to run the flow');
      expect(!!empty.querySelector('button'), true, 'with a Run button');

      h.view.outputViews.activate('canvas');
      expect(tab.pane.style.display, 'none', 'the tab pane hides');
      expect(canvasPane.style.display === 'none', false, 'the canvas is back');
    } finally {
      destroyView(h);
    }
  }, {timeout: 60000});

  test('a SetVar terminal carrying a table gets a tab, named by its variable', async () => {
    const setVarFunc = DG.Func.find({name: 'SetVar'})[0];
    if (!setVarFunc) {
      expect(true, true); // no live backend — nothing to check
      return;
    }
    const h = await makeView();
    try {
      const input = await addNode(h.flow, 'Inputs/Table Input', 0, 0);
      const setVar = await addNode(h.flow, ensureFuncNodeType(setVarFunc), 240, 0);
      setVar.inputValues['variableName'] = 'MyResult';
      await h.flow.addConnectionByKeys(input.id, 'table', setVar.id, 'value');
      h.flow.notifyNodeParamsChanged(setVar.id);
      expect(!!chipByParam(h, 'MyResult'), true, 'SetVar acts exactly as an output (Q2)');

      // A scalar-fed SetVar stays out.
      const constNode = await addNode(h.flow, 'Constants/Int', 0, 200);
      const setVar2 = await addNode(h.flow, ensureFuncNodeType(setVarFunc), 240, 200);
      setVar2.inputValues['variableName'] = 'MyNumber';
      await h.flow.addConnectionByKeys(constNode.id, 'value', setVar2.id, 'value');
      h.flow.notifyNodeParamsChanged(setVar2.id);
      expect(!!chipByParam(h, 'MyNumber'), false, 'a scalar SetVar gets no tab');
    } finally {
      destroyView(h);
    }
  }, {timeout: 60000});

  test('a value materializes the TableView lazily; a new value refreshes the SAME view', async () => {
    const h = await makeView();
    try {
      const out = await addNode(h.flow, 'Outputs/Table Output', 100, 0);
      out.properties['paramName'] = 'result';
      h.flow.notifyNodeParamsChanged(out.id);
      const tab = h.view.outputViews.getTab(out.id)!;

      h.view.outputViews.setValue(out.id, numericDf('people', 3));
      expect(tab.tv == null, true, 'no TableView until the tab is activated (lazy)');
      const chip = chipByParam(h, 'result')!;
      expect(chip.dataset.state, 'ready');
      expect(chip.textContent!.includes('people'), true, 'the tab is named as the table');

      h.view.outputViews.activate(out.id);
      expect(await until(() => tab.tv != null, 15000), true, 'activation creates the TableView');
      expect(await until(() => tab.pane.querySelector('.d4-grid') != null, 15000), true,
        'the grid rendered inside the pane');
      const tv = tab.tv!;
      const root = tv.root;

      h.view.outputViews.setValue(out.id, numericDf('people2', 5));
      expect(tab.tv === tv, true, 'refresh rebinds, never recreates');
      expect(tv.root === root, true, 'same live DOM');
      expect(tv.dataFrame.rowCount, 5, 'the new value landed');
      expect(chipByParam(h, 'result')!.textContent!.includes('people2'), true, 'label tracks the table name');

      h.view.outputViews.markStale(out.id);
      expect(chipByParam(h, 'result')!.dataset.state, 'stale', 'invalidation shows the amber dot');
      expect(tab.pane.querySelector('.d4-grid') != null, true, 'the last table stays visible (Q3)');
    } finally {
      destroyView(h);
    }
  }, {timeout: 90000});

  test('switching tabs swaps the ribbon and toolbox and restores them on Canvas', async () => {
    const h = await makeView();
    try {
      const out = await addNode(h.flow, 'Outputs/Table Output', 100, 0);
      h.view.outputViews.setValue(out.id, numericDf('swap', 4));

      // The Dart toolbox getter may wrap the set element — assert by
      // containment of the function browser's root, not reference identity.
      const browserRoot = (h.view as unknown as {functionBrowser: {root: HTMLElement}}).functionBrowser.root;
      const toolboxHasBrowser = (): boolean => {
        const tb = h.view.toolbox;
        return tb === browserRoot || tb.contains(browserRoot);
      };
      const saveButton = (h.view as unknown as {saveButton: HTMLElement}).saveButton;
      const inRibbon = (el: Element | null): boolean =>
        el != null && h.view.getRibbonPanels().flat().some((p) => p === el || p.contains(el));
      const ribbonHasTid = (tid: string): boolean =>
        h.view.getRibbonPanels().flat().some((p) =>
          p.getAttribute('data-testid') === tid || p.querySelector(`[data-testid="${tid}"]`) != null);
      expect(toolboxHasBrowser(), true, 'Flow\'s toolbox starts as the function browser');
      expect(inRibbon(saveButton), true, 'Flow\'s Save button starts in the ribbon');
      expect(ribbonHasTid('ff-ribbon-run'), true, 'Flow\'s Run icon starts in the ribbon');

      h.view.outputViews.activate(out.id);
      const tab = h.view.outputViews.getTab(out.id)!;
      expect(await until(() => tab.tv != null, 15000), true, 'TableView created');
      const tvItems = tab.tv!.getRibbonPanels().flat().length;
      expect(tvItems > 0, true, `a detached TableView has ribbon items (got ${tvItems})`);
      expect(await until(() => !toolboxHasBrowser(), 5000), true,
        'the toolbox swapped to the TableView\'s');
      expect(ribbonHasTid('ff-ribbon-run'), false, 'Flow\'s icons left the ribbon');
      expect(inRibbon(saveButton), true, 'Flow\'s Save pill stays on the table tab (saves the layout)');
      const swapped = h.view.getRibbonPanels().flat().length;
      expect(swapped > 1, true, `the TableView's items landed on the host ribbon (got ${swapped})`);
      // The TableView's own (core-hidden) Save button must not be copied.
      const hiddenCopied = h.view.getRibbonPanels().flat().some((p) => {
        const inner = p.firstElementChild as HTMLElement | null;
        return (p as HTMLElement).style.display === 'none' || inner?.style?.display === 'none';
      });
      expect(hiddenCopied, false, 'no vestigial hidden items copied from the TableView');

      h.view.outputViews.activate('canvas');
      expect(toolboxHasBrowser(), true, 'Canvas restores Flow\'s toolbox');
      expect(inRibbon(saveButton), true, 'Flow\'s Save is back home');
      expect(ribbonHasTid('ff-ribbon-run'), true, 'Flow\'s ribbon panels are restored (same references)');

      // Regression: setRibbonPanels MOVES elements — a second activation must
      // re-use the captured inner elements, not re-read empty tv wrappers.
      h.view.outputViews.activate(out.id);
      expect(await until(() => !ribbonHasTid('ff-ribbon-run'), 5000), true,
        'second switch swaps again');
      const secondCount = h.view.getRibbonPanels().flat().length;
      expect(secondCount, swapped, `second switch shows the same items (got ${secondCount}, was ${swapped})`);
      const emptyItems = h.view.getRibbonPanels().flat()
        .filter((p) => p.children.length === 0 && (p.textContent ?? '') === '').length;
      expect(emptyItems, 0, 'no empty husks on the second switch');

      h.view.outputViews.activate('canvas');
      expect(ribbonHasTid('ff-ribbon-run'), true, 'Canvas restores again after the second round-trip');
    } finally {
      destroyView(h);
    }
  }, {timeout: 90000});

  test('the core Save-project dialog opens seeded with a flow table and cancels to null', async () => {
    // End-to-end probe of the publish path: the Dart-side interop
    // (grok_Project_OpenSaveDialog → ProjectMeta.publishTables) must be
    // registered, the dialog must show the passed table (workspace NOT
    // scanned), and closing without saving must resolve null.
    const showSaveDialog = (DG.Project as unknown as {
      showSaveDialog?: (o: object) => Promise<unknown>;
    }).showSaveDialog;
    expect(typeof showSaveDialog, 'function', 'DG.Project.showSaveDialog exists (js-api)');
    expect(typeof (window as unknown as Record<string, unknown>).grok_Project_OpenSaveDialog,
      'function', 'the Dart interop is registered');

    const df = numericDf('ff pub probe', 4);
    let dartError: unknown = null;
    const promise = (showSaveDialog!.call(DG.Project, {tables: [df], name: 'FlowPublishProbe'}) as
      Promise<unknown>).catch((e) => {
      dartError = e;
      return 'ERROR';
    });
    // The project name sits in an input VALUE (not textContent) — identify the
    // dialog by the offered table row instead.
    const dialogEl = (): HTMLElement | null => {
      for (const d of Array.from(document.querySelectorAll('.d4-dialog')) as HTMLElement[]) {
        if ((d.textContent ?? '').includes('ff pub probe')) return d;
      }
      return null;
    };
    const opened = await until(() => dialogEl() != null, 15000);
    expect(opened, true, `the Save-project dialog opened with the passed table` +
      (dartError != null ? ` (dart error: ${dartError})` : ''));
    const dlg = dialogEl()!;
    const nameBox = Array.from(dlg.querySelectorAll('input[type="text"], input#name'))
      .find((i) => (i as HTMLInputElement).value === 'FlowPublishProbe');
    expect(!!nameBox, true, 'the project name is prefilled');

    // Cancel: press Escape (Modal closes on ESC), fall back to the close icon.
    document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', keyCode: 27, bubbles: true} as
      KeyboardEventInit));
    if (dialogEl() != null)
      (dlg.querySelector('.d4-dialog-header i[class*="close"], .d4-dialog-header .fa-times') as
        HTMLElement | null)?.click();
    expect(await until(() => dialogEl() == null, 10000), true, 'the dialog closed');
    const result = await Promise.race([promise, new Promise((r) => setTimeout(() => r('TIMEOUT'), 10000))]);
    expect(result == null, true, `a cancelled dialog resolves null (got ${result})`);
  }, {timeout: 60000});

  test('publishing rebinds the same project and ships layouts without a live view', async () => {
    const showSaveDialog = (DG.Project as unknown as {
      showSaveDialog: (o: object) => Promise<unknown>;
    }).showSaveDialog;
    const dialogEl = (marker: string): HTMLElement | null => {
      for (const d of Array.from(document.querySelectorAll('.d4-dialog')) as HTMLElement[]) {
        if ((d.textContent ?? '').includes(marker)) return d;
      }
      return null;
    };
    const openAndOk = async (opts: object, marker: string): Promise<DG.Project | null> => {
      const promise = showSaveDialog.call(DG.Project, opts) as Promise<unknown>;
      const ok = await until(() => dialogEl(marker) != null, 15000);
      expect(ok, true, `the dialog opened (${marker})`);
      const btn = Array.from(dialogEl(marker)!.querySelectorAll('button'))
        .find((b) => (b.textContent ?? '').trim().toUpperCase() === 'OK') as HTMLElement | undefined;
      expect(!!btn, true, 'the OK button found');
      btn!.click();
      return await Promise.race([
        promise,
        new Promise((r) => setTimeout(() => r('TIMEOUT'), 60000)),
      ]) as DG.Project | null;
    };

    // A real layout state string, captured from a scratch view.
    const tmp = grok.shell.addTableView(numericDf('ff bind layout src', 3));
    const layoutState = tmp.saveLayout().viewState;
    tmp.close();

    let project: DG.Project | null = null;
    try {
      const df = numericDf('ff pub bind', 3);
      project = await openAndOk(
        {tables: [df], layouts: [layoutState], name: 'FlowBindProbe'}, 'ff pub bind');
      expect(project != null && project !== ('TIMEOUT' as unknown), true, 'first publish saved a project');
      const children = project!.children;
      expect(children.some((c) => c instanceof DG.ViewInfo), true,
        'the layout shipped as a ViewInfo without a live view');
      expect(children.some((c) => c instanceof DG.TableInfo), true, 'the table is a child');

      // Re-publish with the binding: the SAME project is updated, not a new one.
      const df2 = numericDf('ff pub bind', 5);
      const again = await openAndOk(
        {tables: [df2], name: 'FlowBindProbe', project: project!.id}, 'ff pub bind');
      expect(again != null && again !== ('TIMEOUT' as unknown), true, 'second publish saved');
      expect(again!.id, project!.id, 're-publishing updates the bound project (no new project per save)');
      project = again;
    } finally {
      if (project != null && (project as unknown) !== 'TIMEOUT')
        await grok.dapi.projects.delete(project).catch(() => {});
    }
  }, {timeout: 180000});

  test('tab layouts persist by paramName and survive a save → load round-trip', async () => {
    const h = await makeView();
    let h2: ViewHarness | null = null;
    try {
      const out = await addNode(h.flow, 'Outputs/Table Output', 100, 0);
      out.properties['paramName'] = 'result';
      h.flow.notifyNodeParamsChanged(out.id);
      h.view.outputViews.setValue(out.id, numericDf('layout-df', 6));
      h.view.outputViews.activate(out.id);
      const tab = h.view.outputViews.getTab(out.id)!;
      expect(await until(() => tab.tv != null, 15000), true, 'TableView created');

      tab.tv!.addViewer('Scatter plot');
      expect(await until(() => tab.tv!.viewers.length >= 2, 15000), true, 'scatter plot added');

      const layouts = h.view.outputViews.captureLayouts();
      expect(typeof layouts['result'] === 'string' && layouts['result'].length > 0, true,
        'the layout captured under the paramName');

      // The saved entity body carries the layouts; the graph snapshot used for
      // dirty tracking (`serializeFlow`) does not.
      const body = (h.view as unknown as {entityBodyText(): string}).entityBodyText();
      const doc: FuncFlowDocument = parseFlowBody(body).doc;
      expect(!!doc.outputViews?.['result']?.layout, true, 'the .ffjson carries outputViews');

      // Fresh view + load: node ids remap, the layout still finds its tab.
      h2 = await makeView();
      await h2.view.loadFromDoc(doc);
      const out2 = h2.flow.getNodes().find((n) => n.dgTypeName === 'Outputs/Table Output')!;
      expect(out2.id !== out.id, true, 'node id was remapped on load');
      h2.view.outputViews.setValue(out2.id, numericDf('layout-df', 6));
      h2.view.outputViews.activate(out2.id);
      const tab2 = h2.view.outputViews.getTab(out2.id)!;
      expect(await until(() => tab2.tv != null, 15000), true, 'TableView created in the fresh view');
      expect(await until(() => tab2.tv!.viewers.length >= 2, 15000), true,
        'the saved layout re-created the scatter plot');
    } finally {
      if (h2) destroyView(h2);
      destroyView(h);
    }
  }, {timeout: 120000});
});
