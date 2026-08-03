/** Widget/Viewer outputs: a node whose output is a DG.Widget or DG.Viewer keeps
 *  the live object (captured by reference during the in-tab run) and renders its
 *  `.root` in the bottom output panel; the context-panel meta shows
 *  "widget"/"viewer" instead of "[object Object]". Also covers the panel's
 *  in-view splitter behavior: hidden → expanded on first renderable output,
 *  minimize is remembered (content updates never pop it back up), `clear()`
 *  hides but keeps the preference, disabled panels (embedded hosts) never show.
 *  Pure DOM — no live backend needed. */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {category, test, expect} from '@datagrok-libraries/utils/src/test';
import {NodeExecState, NodeExecStatus, ValueSummary} from '../execution/execution-state';
import {connectionCountLabel} from '../execution/execution-controller';
import {
  buildExecutionMeta, buildValuePreviews, hasRenderablePreview, buildPreview,
} from '../execution/value-inspector';
import {OutputPreviewPanel, OutputPanelState} from '../execution/output-preview';
import {FuncFlowView} from '../funcflow-view';
import {until} from './test-utils';

function renderableState(): NodeExecState {
  const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('x', ['a', 'b'])]);
  return {
    status: NodeExecStatus.completed,
    outputs: {result: {type: 'dataframe', rows: 2, cols: 1, clone: df}},
  };
}

function widgetState(type: 'widget' | 'viewer'): {state: NodeExecState; root: HTMLElement} {
  const root = document.createElement('div');
  root.textContent = 'live-content';
  const state: NodeExecState = {
    status: NodeExecStatus.completed,
    startTime: 0, endTime: 9,
    outputs: {result: {type, value: {root}}},
  };
  return {state, root};
}

category('Flow: execution preview', () => {
  test('edge count labels read "rows × cols" for tables; __pt entries stay out of panels', async () => {
    // A table with both dims → "N × K"; dims-only summaries (pass-throughs)
    // count too; columns keep the bare count; nothing countable → null.
    expect(connectionCountLabel({type: 'dataframe', rows: 1204, cols: 8}), '1,204 × 8');
    expect(connectionCountLabel({type: 'dataframe', rows: 5} as ValueSummary), '5 rows', 'no cols → rows fallback');
    expect(connectionCountLabel({type: 'column', length: 42} as ValueSummary), '42');
    expect(connectionCountLabel(null), null);
    expect(connectionCountLabel({type: 'primitive', value: 3} as ValueSummary), null);

    // Panel listings skip the `__pt` bookkeeping entries (and tolerate null).
    const state: NodeExecState = {
      status: NodeExecStatus.completed,
      outputs: {
        'result': {type: 'primitive', value: 7},
        'table__pt': {type: 'dataframe', rows: 3, cols: 2},
        'other__pt': null as unknown as ValueSummary,
      },
    };
    expect(hasRenderablePreview(state), false, 'dims-only entries are not renderable');
    const meta = buildExecutionMeta(state);
    expect((meta.textContent ?? '').includes('__pt'), false, 'no __pt rows in the context-panel meta');
  });

  test('widget output is renderable and mounts its live root', async () => {
    const {state, root} = widgetState('widget');
    expect(hasRenderablePreview(state), true, 'widget counts as renderable');
    const previews = buildValuePreviews(state);
    expect(previews.contains(root), true, 'the live widget root is mounted in the preview');
  });

  test('viewer output is renderable and mounts its live root', async () => {
    const {state, root} = widgetState('viewer');
    expect(hasRenderablePreview(state), true, 'viewer counts as renderable');
    const preview = buildPreview('result', state.outputs!.result);
    expect(!!preview && preview.contains(root), true, 'the live viewer root is mounted');
  });

  test('context-panel meta names the kind, not "[object Object]"', async () => {
    const {state} = widgetState('widget');
    const meta = buildExecutionMeta(state);
    const text = meta.textContent ?? '';
    expect(text.includes('widget'), true, 'meta says "widget"');
    expect(text.includes('[object Object]'), false, 'no [object Object] leak');
  });

  test('a column output renders as a one-column DataFrame grid, not a text sample', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('x', ['a', 'b', 'c'])]);
    const state: NodeExecState = {
      status: NodeExecStatus.completed,
      outputs: {result: {type: 'column', name: 'x', length: 3, sample: ['a', 'b', 'c'], clone: df}},
    };
    expect(hasRenderablePreview(state), true, 'a column with a captured DataFrame is renderable');
    const preview = buildPreview('result', state.outputs!.result);
    expect(!!preview, true, 'preview built');
    // Grid mode → the fallback text-sample <table> must NOT be used.
    expect(preview!.querySelector('table'), null, 'renders a grid, not the text sample table');
    expect(preview!.childElementCount > 0, true, 'grid root mounted');
  });

  test('the threaded "(modified)" table is suppressed when a column is output', async () => {
    const col = DG.DataFrame.fromColumns([DG.Column.fromStrings('x', ['a', 'b'])]);
    const tbl = DG.DataFrame.fromColumns([DG.Column.fromStrings('y', ['p', 'q'])]);
    const state: NodeExecState = {
      status: NodeExecStatus.completed,
      outputs: {
        result: {type: 'column', name: 'x', length: 2, clone: col},
        'table (modified)': {type: 'dataframe', rows: 2, cols: 1, clone: tbl},
      },
    };
    const previews = buildValuePreviews(state);
    const blocks = previews.querySelectorAll('[data-testid^="ff-preview-block"]');
    expect(blocks.length, 1, 'only the column preview shows; the (modified) table is suppressed');
  });

  test('a column with no captured DataFrame falls back to the text sample', async () => {
    const state: NodeExecState = {
      status: NodeExecStatus.completed,
      outputs: {result: {type: 'column', name: 'x', length: 2, sample: ['a', 'b']}},
    };
    const preview = buildPreview('result', state.outputs!.result);
    expect(!!preview && !!preview.querySelector('table'), true, 'sample table rendered when no clone');
  });

  test('the panel opens on the first renderable output; minimize is remembered', async () => {
    const panel = new OutputPreviewPanel();
    expect(panel.panelState, 'hidden', 'starts hidden');

    // No renderable values → stays hidden.
    panel.showForNode({id: 'n0', label: 'empty'}, {status: NodeExecStatus.completed, outputs: {}});
    expect(panel.panelState, 'hidden', 'a value-less node never opens the panel');

    panel.showForNode({id: 'n1', label: 'a'}, renderableState());
    expect(panel.panelState, 'expanded', 'first renderable output expands the panel');

    panel.minimize();
    expect(panel.panelState, 'minimized');
    panel.showForNode({id: 'n2', label: 'b'}, renderableState());
    expect(panel.panelState, 'minimized', 'content updates never pop a minimized panel back up');
    expect(panel.root.querySelector('[data-testid="ff-output-panel-node"]')!.textContent, 'b',
      'the content still updates in place while minimized');

    panel.expand();
    expect(panel.panelState, 'expanded', 'restoring is explicit');
  });

  test('clear() hides the panel but keeps the minimized preference', async () => {
    const panel = new OutputPreviewPanel();
    panel.showForNode({id: 'n1', label: 'a'}, renderableState());
    panel.minimize();

    panel.clear();
    expect(panel.panelState, 'hidden', 'clear hides the panel');
    expect(panel.root.querySelector('[data-testid="ff-output-panel-content"]')!.childElementCount, 0,
      'stale content is dropped');

    panel.showForNode({id: 'n2', label: 'b'}, renderableState());
    expect(panel.panelState, 'minimized', 'after a clear, new content respects the remembered preference');
  });

  test('the caret click toggles minimized/expanded and reports state changes', async () => {
    const panel = new OutputPreviewPanel();
    const states: OutputPanelState[] = [];
    panel.onStateChanged = (s): void => {states.push(s);};

    panel.showForNode({id: 'n1', label: 'a'}, renderableState());
    const caret = panel.root.querySelector('[data-testid="ff-output-panel-caret"]') as HTMLElement;
    caret.click();
    caret.click();
    expect(states.join(','), 'expanded,minimized,expanded', 'every transition is reported (divider sync)');

    // Only the caret toggles — the rest of the header sits right under the
    // splitter divider, and a near-miss resize click must not collapse the panel.
    const header = panel.root.querySelector('[data-testid="ff-output-panel-header"]') as HTMLElement;
    header.dispatchEvent(new MouseEvent('click', {bubbles: false}));
    expect(panel.panelState, 'expanded', 'clicking the header body does not toggle');
  });

  test('a disabled panel never shows — embedded hosts (dialogs)', async () => {
    const panel = new OutputPreviewPanel({enabled: false});
    panel.showForNode({id: 'n1', label: 'a'}, renderableState());
    expect(panel.panelState, 'hidden', 'disabled → renderable output is ignored');

    panel.setEnabled(true);
    panel.showForNode({id: 'n1', label: 'a'}, renderableState());
    expect(panel.panelState, 'expanded', 'enabling turns it back into a real panel');
  });

  test('re-clicking the same node with the same state does not rebuild the preview', async () => {
    // Every rebuild re-mounts the grids and makes the panel jump — same node +
    // same captured state (state objects are replaced, never mutated, so
    // reference identity IS value identity) must keep the existing DOM.
    const panel = new OutputPreviewPanel();
    const state = renderableState();
    const contentEl = (): Element | null =>
      panel.root.querySelector('[data-testid="ff-output-panel-content"]')?.firstElementChild ?? null;

    panel.showForNode({id: 'n1', label: 'a'}, state);
    const first = contentEl();
    expect(!!first, true, 'preview rendered');

    // Same node, same state object → untouched DOM.
    panel.showForNode({id: 'n1', label: 'a'}, state);
    expect(contentEl() === first, true, 'same node + same state → no rebuild');

    // A fresh state object (e.g. after a re-run) → rebuilt.
    panel.showForNode({id: 'n1', label: 'a'}, renderableState());
    const rebuilt = contentEl();
    expect(!!rebuilt && rebuilt !== first, true, 'a new captured state rebuilds the preview');

    // A different node → rebuilt as well.
    panel.showForNode({id: 'n2', label: 'b'}, state);
    expect(contentEl() !== rebuilt, true, 'another node rebuilds the preview');
  });

  test('pin freezes the preview on the shown node; unpin releases it', async () => {
    const panel = new OutputPreviewPanel();
    const label = (): string =>
      panel.root.querySelector('[data-testid="ff-output-panel-node"]')!.textContent ?? '';
    const pin = panel.root.querySelector('[data-testid="ff-output-panel-pin"]') as HTMLElement;

    expect(pin.style.display, 'none', 'no content → nothing to pin, icon hidden');

    panel.showForNode({id: 'n1', label: 'a'}, renderableState());
    expect(pin.style.display === 'none', false, 'the icon appears with content');
    pin.click();
    expect(panel.pinnedNodeId, 'n1', 'pin captures the shown node');
    expect(pin.dataset.pinned, 'true', 'pinned state on the icon');

    panel.showForNode({id: 'n2', label: 'b'}, renderableState());
    expect(label(), 'a', 'other nodes cannot replace a pinned preview');
    expect(panel.currentNodeId, 'n1', 'still showing the pinned node');

    // Fresh state for the pinned node itself still re-renders (the controller
    // re-shows it on node-complete — the "watch the chart while editing" case).
    panel.showForNode({id: 'n1', label: 'a2'}, renderableState());
    expect(label(), 'a2', 'the pinned node itself updates in place');

    pin.click();
    expect(panel.pinnedNodeId, null, 'second click unpins');
    panel.showForNode({id: 'n2', label: 'b'}, renderableState());
    expect(label(), 'b', 'unpinned → node clicks switch the preview again');
  });

  test('unpinning via the icon reports onUnpinned; programmatic unpin stays silent', async () => {
    // The controller listens to onUnpinned to re-show the currently selected
    // node — re-clicking an already-selected node fires no selection event, so
    // without this "unpin to see what I clicked while pinned" showed stale content.
    const panel = new OutputPreviewPanel();
    let fires = 0;
    panel.onUnpinned = (): void => {fires++;};
    panel.showForNode({id: 'n1', label: 'a'}, renderableState());
    panel.togglePin();
    expect(fires, 0, 'pinning does not fire');
    panel.togglePin();
    expect(fires, 1, 'a user unpin fires');
    panel.togglePin();
    panel.unpin();
    expect(fires, 1, 'programmatic unpin (node deleted / graph replaced) stays silent');
  });

  test('markUpdating overlays a spinner on the kept content; the next render releases it', async () => {
    // The pinned-preview recompute path: instead of hiding + re-docking (which
    // reads as jumping), the stale content stays with a "Recalculating…"
    // indicator until the fresh render lands.
    const panel = new OutputPreviewPanel();
    const content = panel.root.querySelector('[data-testid="ff-output-panel-content"]') as HTMLElement;
    const shadow = (): Element | null => content.querySelector('.d4-update-shadow');

    panel.markUpdating();
    expect(shadow(), null, 'no content → no indicator (nothing to overlay)');

    panel.showForNode({id: 'n1', label: 'a'}, renderableState());
    panel.togglePin();
    panel.markUpdating();
    expect(panel.panelState, 'expanded', 'the panel stays visible — no hide');
    expect(!!shadow(), true, 'indicator overlays the kept content');
    expect(content.childElementCount >= 2, true, 'the stale preview is still mounted under it');

    panel.showForNode({id: 'n1', label: 'a'}, renderableState());
    expect(shadow(), null, 'a fresh render releases the indicator');

    panel.markUpdating();
    panel.clearUpdating();
    expect(shadow(), null, 'clearUpdating releases it without a render (error / skipped run)');

    panel.markUpdating();
    panel.clear();
    expect(shadow(), null, 'clear() drops the indicator with the content');
  });

  test('the pin survives clear(); unpin() drops it', async () => {
    const panel = new OutputPreviewPanel();
    panel.showForNode({id: 'n1', label: 'a'}, renderableState());
    panel.togglePin();

    panel.clear();
    expect(panel.pinnedNodeId, 'n1', 'clear (invalidation / new run) keeps the pin');
    panel.showForNode({id: 'n2', label: 'b'}, renderableState());
    expect(panel.panelState, 'hidden', 'other nodes stay gated even after a clear');
    panel.showForNode({id: 'n1', label: 'a'}, renderableState());
    expect(panel.panelState, 'expanded', 'the pinned node\'s fresh result re-opens the panel');

    panel.unpin();
    panel.showForNode({id: 'n2', label: 'b'}, renderableState());
    expect(panel.root.querySelector('[data-testid="ff-output-panel-node"]')!.textContent, 'b',
      'unpin releases the gate');
  });

  test('the panel is a pane of the view splitter; embedded views create it disabled', async () => {
    const view = new FuncFlowView();
    try {
      expect(view.root.contains(view.outputPreview.root), true, 'panel root lives inside the view root');
      expect(!!view.root.querySelector('.ui-split-v-divider'), true, 'the vertical splitter divider exists');
      expect(view.outputPreview.panelState, 'hidden', 'no outputs yet → hidden');
      expect(view.outputPreview.isEnabled, true, 'the real editor view has the panel on');
    } finally {
      await new Promise((r) => setTimeout(r, 120)); // let the deferred initEditor run before teardown
      ((view as any).flow)?.destroy?.();
      view.root.remove();
    }

    const embedded = new FuncFlowView([], {outputPanel: false});
    try {
      expect(embedded.outputPreview.isEnabled, false, 'embedded hosts get a disabled panel');
      embedded.outputPreview.showForNode({id: 'n1', label: 'a'}, renderableState());
      expect(embedded.outputPreview.panelState, 'hidden', 'it never shows in a dialog');
      embedded.enableOutputPanel();
      expect(embedded.outputPreview.isEnabled, true, 'Open In Editor re-enables it');
    } finally {
      await new Promise((r) => setTimeout(r, 120));
      ((embedded as any).flow)?.destroy?.();
      embedded.root.remove();
    }
  }, {timeout: 30000});

  test('opening the output preview minimizes the overview minimap', async () => {
    // Regression: the minimap and the bottom output panel crowd the same corner,
    // so the minimap auto-minimizes to its header the first time the preview
    // opens. This wiring was dropped in the dock → splitter rework.
    const view = new FuncFlowView();
    const host = ui.div([view.root], {style: {
      width: '900px', height: '600px', position: 'absolute', left: '-10000px',
    }});
    document.body.appendChild(host);
    try {
      // Wait for the deferred editor build to mount the minimap.
      await until(() => view.root.querySelector('.ff-minimap') != null);
      const mm = view.root.querySelector('.ff-minimap') as HTMLElement;
      expect(!!mm, true, 'minimap mounted');
      expect(mm.dataset.collapsed, 'false', 'minimap starts expanded');
      expect(view.outputPreview.panelState, 'hidden', 'preview starts hidden');

      // Opening the preview (hidden → expanded) minimizes the minimap.
      view.outputPreview.showForNode({id: 'n1', label: 'a'}, renderableState());
      expect(view.outputPreview.panelState, 'expanded', 'preview opened');
      expect(mm.dataset.collapsed, 'true', 'minimap minimized when the preview opened');

      // One-shot: manually reopening the minimap while the preview stays open sticks.
      view.setMinimapCollapsed(false);
      view.outputPreview.showForNode({id: 'n2', label: 'b'}, renderableState());
      expect(mm.dataset.collapsed, 'false', 'no re-collapse while the preview was already open');
    } finally {
      ((view as any).flow)?.destroy?.();
      host.remove();
    }
  }, {timeout: 30000});

  test('the canvas container clips — no scrollbars around the transformed canvas', async () => {
    // Regression: as a `div.ui-div` directly inside the splitter's `.ui-box`
    // pane, the canvas matched the core rule `div.ui-box > div.ui-div
    // { overflow: auto !important }`, which overrode the package's
    // `overflow: hidden` and grew giant scrollbars around the Rete canvas
    // (node views live at large canvas coordinates).
    const view = new FuncFlowView();
    const host = ui.div([view.root], {style: {
      width: '900px', height: '600px', position: 'absolute', left: '-10000px',
    }});
    document.body.appendChild(host);
    try {
      await new Promise((r) => setTimeout(r, 150)); // let the deferred initEditor build the canvas
      const canvas = view.root.querySelector('.funcflow-canvas-container') as HTMLElement;
      expect(!!canvas, true, 'canvas container exists');
      const style = getComputedStyle(canvas);
      const diag = `self="${canvas.className}" style="${canvas.getAttribute('style') ?? ''}" ` +
        `parent="${canvas.parentElement?.className}" grandparent="${canvas.parentElement?.parentElement?.className}"`;
      expect(style.overflowX, 'hidden', `canvas clips horizontally — ${diag}`);
      expect(style.overflowY, 'hidden', `canvas clips vertically — ${diag}`);
    } finally {
      ((view as any).flow)?.destroy?.();
      host.remove();
    }
  }, {timeout: 30000});

  test('a widget without a real root is not renderable', async () => {
    const state: NodeExecState = {
      status: NodeExecStatus.completed,
      outputs: {result: {type: 'widget', value: {}}},
    };
    expect(hasRenderablePreview(state), false, 'no .root → nothing to render');
    expect(buildPreview('result', state.outputs!.result), null, 'buildPreview returns null');
  });

  const wsBtn = (el: HTMLElement | null): Element | null =>
    el?.querySelector('[data-testid="ff-add-to-workspace"]') ?? null;

  test('a dataframe preview carries an "Add to workspace" button', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('x', ['a', 'b'])]);
    const preview = buildPreview('result', {type: 'dataframe', rows: 2, cols: 1, clone: df});
    expect(!!wsBtn(preview), true, 'the button is overlaid on the grid');
  });

  test('a column preview carries an "Add to workspace" button', async () => {
    const col = DG.DataFrame.fromColumns([DG.Column.fromStrings('x', ['a', 'b', 'c'])]);
    const preview = buildPreview('result', {type: 'column', name: 'x', length: 3, clone: col});
    expect(!!wsBtn(preview), true, 'the one-column grid gets the button too');
  });

  test('an in-place column output previews the whole table (scrolled), not a lone column', async () => {
    // The output column belongs to the node's single input table (Add New Column
    // idiom): the summary carries the full table + the column to scroll to, so
    // the column is shown in the context of its table.
    const table = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('name', ['a', 'b']),
      DG.Column.fromFloat32Array('extra', new Float32Array([1, 2])),
    ]);
    const oneCol = DG.DataFrame.fromColumns([DG.Column.fromFloat32Array('extra', new Float32Array([1, 2]))]);
    const summary = {
      type: 'column', name: 'extra', length: 2, clone: oneCol,
      tableClone: table, scrollToColumn: 'extra',
    } as const;
    const preview = buildPreview('result', summary);
    expect(!!preview, true, 'preview built');
    // Grid mode → the fallback text-sample <table> must NOT appear.
    expect(preview!.querySelector('table'), null, 'renders a grid, not a text sample');
    expect(preview!.childElementCount > 0, true, 'grid mounted');
    expect(!!wsBtn(preview), true, 'the in-place table preview is addable to the workspace');
  });

  test('a viewer preview has the workspace button, plus the gear when editable', async () => {
    const {state} = widgetState('viewer');
    const summary = state.outputs!.result;

    const noGear = buildPreview('result', summary);
    expect(!!wsBtn(noGear), true, 'the workspace button is always present on a viewer');
    expect(noGear!.querySelector('[data-testid="ff-viewer-edit"]'), null, 'no gear without an edit handler');

    const withGear = buildPreview('result', summary, () => {});
    expect(!!wsBtn(withGear), true, 'the workspace button coexists with the gear');
    expect(!!withGear!.querySelector('[data-testid="ff-viewer-edit"]'), true, 'the edit gear is present');
  });

  test('a single output fills the panel; multiple outputs split side by side', async () => {
    const df = (name: string): DG.DataFrame =>
      DG.DataFrame.fromColumns([DG.Column.fromStrings(name, ['a', 'b'])]);

    // One output → the block is mounted directly (no splitter).
    const one = buildValuePreviews({
      status: NodeExecStatus.completed,
      outputs: {result: {type: 'dataframe', rows: 2, cols: 1, clone: df('x')}},
    });
    expect(one.querySelectorAll('.funcflow-preview-block').length, 1, 'one preview block');
    expect(one.querySelector('.ui-split-h'), null, 'a lone output is not wrapped in a splitter');

    // Two renderable outputs (e.g. a multi-output func) → side-by-side splitH.
    const two = buildValuePreviews({
      status: NodeExecStatus.completed,
      outputs: {
        result1: {type: 'dataframe', rows: 2, cols: 1, clone: df('x')},
        result2: {type: 'dataframe', rows: 2, cols: 1, clone: df('y')},
      },
    });
    expect(two.querySelectorAll('.funcflow-preview-block').length, 2, 'both previews built');
    expect(!!two.querySelector('.ui-split-h'), true, 'multiple outputs share a horizontal splitter');
  });
});
