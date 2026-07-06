/** Widget/Viewer outputs: a node whose output is a DG.Widget or DG.Viewer keeps
 *  the live object (captured by reference during the in-tab run) and renders its
 *  `.root` in the docked preview; the context-panel meta shows "widget"/"viewer"
 *  instead of "[object Object]". Pure DOM — no live backend needed. */
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/utils/src/test';
import {NodeExecState, NodeExecStatus} from '../execution/execution-state';
import {
  buildExecutionMeta, buildValuePreviews, hasRenderablePreview, buildPreview,
} from '../execution/value-inspector';
import {OutputPreviewPanel} from '../execution/output-preview';

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

  test('onDocked fires when the bottom panel is first created, not on updates', async () => {
    // The view wires this to minimize the minimap (they share the same
    // corner). It must fire once per dock creation — an in-place content
    // update reuses the open dock and stays silent.
    const panel = new OutputPreviewPanel();
    let docked = 0;
    panel.onDocked = (): void => {docked++;};
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('x', ['a', 'b'])]);
    const state: NodeExecState = {
      status: NodeExecStatus.completed,
      outputs: {result: {type: 'dataframe', rows: 2, cols: 1, clone: df}},
    };
    try {
      // No renderable values → no dock, no callback.
      panel.showForNode({id: 'n0', label: 'empty'}, {status: NodeExecStatus.completed, outputs: {}});
      expect(docked, 0, 'no dock for a value-less node');

      panel.showForNode({id: 'n1', label: 'a'}, state);
      expect(docked, 1, 'dock created → onDocked fired');
      panel.showForNode({id: 'n2', label: 'b'}, state);
      expect(docked, 1, 'in-place update → not fired again');

      // Re-docking after a close fires again (the dock is newly created).
      panel.close();
      panel.showForNode({id: 'n3', label: 'c'}, state);
      expect(docked, 2, 'a fresh dock after close fires again');
    } finally {
      panel.close();
    }
  });

  test('re-clicking the same node with the same state does not rebuild the preview', async () => {
    // Every rebuild re-mounts the grids and makes the panel jump — same node +
    // same captured state (state objects are replaced, never mutated, so
    // reference identity IS value identity) must keep the existing DOM.
    const panel = new OutputPreviewPanel();
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('x', ['a', 'b'])]);
    const mkState = (): NodeExecState => ({
      status: NodeExecStatus.completed,
      outputs: {result: {type: 'dataframe', rows: 2, cols: 1, clone: df}},
    });
    const state = mkState();
    const contentEl = (): Element | null =>
      document.querySelector('[data-testid="ff-output-panel"]')?.firstElementChild ?? null;
    try {
      panel.showForNode({id: 'n1', label: 'a'}, state);
      const first = contentEl();
      expect(!!first, true, 'preview rendered');

      // Same node, same state object → untouched DOM.
      panel.showForNode({id: 'n1', label: 'a'}, state);
      expect(contentEl() === first, true, 'same node + same state → no rebuild');

      // A fresh state object (e.g. after a re-run) → rebuilt.
      panel.showForNode({id: 'n1', label: 'a'}, mkState());
      const rebuilt = contentEl();
      expect(!!rebuilt && rebuilt !== first, true, 'a new captured state rebuilds the preview');

      // A different node → rebuilt as well.
      panel.showForNode({id: 'n2', label: 'b'}, state);
      expect(contentEl() !== rebuilt, true, 'another node rebuilds the preview');
    } finally {
      panel.close();
    }
  });

  test('a widget without a real root is not renderable', async () => {
    const state: NodeExecState = {
      status: NodeExecStatus.completed,
      outputs: {result: {type: 'widget', value: {}}},
    };
    expect(hasRenderablePreview(state), false, 'no .root → nothing to render');
    expect(buildPreview('result', state.outputs!.result), null, 'buildPreview returns null');
  });
});
