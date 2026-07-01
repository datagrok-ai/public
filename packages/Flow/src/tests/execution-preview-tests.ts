/** Widget/Viewer outputs: a node whose output is a DG.Widget or DG.Viewer keeps
 *  the live object (captured by reference during the in-tab run) and renders its
 *  `.root` in the docked preview; the context-panel meta shows "widget"/"viewer"
 *  instead of "[object Object]". Pure DOM — no live backend needed. */
import {category, test, expect} from '@datagrok-libraries/utils/src/test';
import {NodeExecState, NodeExecStatus} from '../execution/execution-state';
import {
  buildExecutionMeta, buildValuePreviews, hasRenderablePreview, buildPreview,
} from '../execution/value-inspector';

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

  test('a widget without a real root is not renderable', async () => {
    const state: NodeExecState = {
      status: NodeExecStatus.completed,
      outputs: {result: {type: 'widget', value: {}}},
    };
    expect(hasRenderablePreview(state), false, 'no .root → nothing to render');
    expect(buildPreview('result', state.outputs!.result), null, 'buildPreview returns null');
  });
});
