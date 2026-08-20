/** In-node viewer/widget preview: the title-bar toggle, the resizable container
 *  that mounts the captured live root, serialization of the choice/size, and the
 *  bottom-panel ownership handoff. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {ClassicPreset} from 'rete';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions, createNode, getRegisteredFuncs} from '../rete/node-factory';
import {FlowEditor} from '../rete/flow-editor';
import {
  FlowNode, supportsInlinePreview, inlinePreviewEnabled, inlinePreviewSize,
  INLINE_PREVIEW_PROP, INLINE_PREVIEW_SIZE_PROP, INLINE_PREVIEW_DEFAULT_SIZE, INLINE_HOSTED_DATA_KEY,
} from '../rete/scheme';
import {getSocket} from '../rete/sockets';
import {estimateNodeHeight, estimateNodeWidth} from '../rete/graph-layout';
import {serializeFlow, deserializeFlow} from '../serialization/flow-serializer';
import {buildPreview} from '../execution/value-inspector';
import {OutputPreviewPanel} from '../execution/output-preview';
import {ExecutionController} from '../execution/execution-controller';
import {NodeExecStatus} from '../execution/execution-state';
import {makeEditor, destroyEditor, addNode, until, TestEditor} from './test-utils';

const SETTINGS = {name: 'T', description: '', tags: []};

function outputNodeOf(dgType: string): FlowNode {
  const n = new FlowNode(dgType);
  n.addOutput('result', new ClassicPreset.Output(getSocket(dgType), 'result'));
  return n;
}

const SVG_SAMPLE =
  '<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10"><rect width="10" height="10"/></svg>';

function funcTypeName(name: string): string | null {
  const lc = name.toLowerCase();
  return getRegisteredFuncs().find((f) => f.func.name.toLowerCase() === lc)?.nodeTypeName ?? null;
}

/** Editor whose in-node previews resolve through an ExecutionController, wired
 *  the same way FuncFlowView wires the real one. */
function makeWiredEditor(): TestEditor & {ctrl: ExecutionController} {
  const container = ui.div([], {style: {width: '1000px', height: '700px', position: 'absolute', left: '-10000px'}});
  document.body.appendChild(container);
  const box: {ctrl?: ExecutionController} = {};
  const flow = new FlowEditor(container, {
    getInlinePreviewContent: (nodeId) => box.ctrl?.inlinePreviewRoot(nodeId) ?? null,
    onInlinePreviewToggled: (nodeId) => box.ctrl?.syncInlinePreviewOwnership(nodeId),
  });
  const ctrl = new ExecutionController(flow);
  box.ctrl = ctrl;
  return {flow, container, ctrl};
}

function previewEl(container: HTMLElement, nodeId: string): HTMLElement | null {
  return container.querySelector(`[data-node-id="${nodeId}"] [data-testid="ff-node-preview"]`);
}

category('Flow: inline preview', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('only viewer-, widget- and graphics-producing nodes support the in-node preview', async () => {
    const viewer = createNode('Viewers/Scatter Plot')!;
    expect(supportsInlinePreview(viewer), true, 'a viewer node supports it');

    expect(supportsInlinePreview(outputNodeOf('widget')), true, 'a widget output supports it');
    expect(supportsInlinePreview(outputNodeOf('graphics')), true,
      'a graphics output supports it (Gasteiger-style scripts)');

    const tableInput = createNode('Inputs/Table Input')!;
    expect(supportsInlinePreview(tableInput), false, 'a table input does not');

    const selTable = createNode('Utilities/Select Table')!;
    expect(supportsInlinePreview(selTable), false, 'a table-producing utility does not');
  });

  test('default off; enabling sets the property; size falls back to 300×300', async () => {
    const viewer = createNode('Viewers/Scatter Plot')!;
    expect(inlinePreviewEnabled(viewer), false, 'off by default');
    viewer.properties[INLINE_PREVIEW_PROP] = true;
    expect(inlinePreviewEnabled(viewer), true, 'the property turns it on');
    expect(inlinePreviewSize(viewer).width, INLINE_PREVIEW_DEFAULT_SIZE, 'default width');
    expect(inlinePreviewSize(viewer).height, INLINE_PREVIEW_DEFAULT_SIZE, 'default height');
    viewer.properties[INLINE_PREVIEW_SIZE_PROP] = {width: -5, height: 'x'};
    expect(inlinePreviewSize(viewer).width, INLINE_PREVIEW_DEFAULT_SIZE, 'malformed width falls back');
    expect(inlinePreviewSize(viewer).height, INLINE_PREVIEW_DEFAULT_SIZE, 'malformed height falls back');

    // A non-capable node never reads as enabled, even with a stray property.
    const tableInput = createNode('Inputs/Table Input')!;
    tableInput.properties[INLINE_PREVIEW_PROP] = true;
    expect(inlinePreviewEnabled(tableInput), false, 'gated on capability');
  });

  test('the choice and the size survive a save/load round-trip', async () => {
    const e = makeEditor();
    const e2 = makeEditor();
    try {
      const viewer = await addNode(e.flow, 'Viewers/Scatter Plot');
      await e.flow.setInlinePreview(viewer.id, true);
      viewer.properties[INLINE_PREVIEW_SIZE_PROP] = {width: 420, height: 260};

      const doc = serializeFlow(e.flow, {scriptName: 'T', scriptDescription: '', tags: []});
      await deserializeFlow(doc, e2.flow);
      const loaded = e2.flow.getNodes().find((n) => n.dgTypeName === 'Viewers/Scatter Plot');
      expect(loaded != null, true, 'the viewer node loaded');
      expect(inlinePreviewEnabled(loaded!), true, 'the in-node preview choice is remembered');
      expect(inlinePreviewSize(loaded!).width, 420, 'width is remembered');
      expect(inlinePreviewSize(loaded!).height, 260, 'height is remembered');
    } finally {
      destroyEditor(e2);
      destroyEditor(e);
    }
  });

  test('the title-bar toggle renders only on capable nodes and flips the preview', async () => {
    const e = makeEditor();
    try {
      const viewer = await addNode(e.flow, 'Viewers/Scatter Plot');
      const tableInput = await addNode(e.flow, 'Inputs/Table Input', 0, 200);
      const toggleOf = (id: string): HTMLElement | null =>
        e.container.querySelector(`[data-node-id="${id}"] [data-testid="ff-node-preview-toggle"]`);

      expect(await until(() => toggleOf(viewer.id) != null), true, 'the viewer node has the toggle');
      expect(toggleOf(tableInput.id), null, 'a non-capable node has no toggle');
      expect(toggleOf(viewer.id)!.dataset.on, 'false', 'reads off before the first click');

      toggleOf(viewer.id)!.click();
      expect(await until(() => inlinePreviewEnabled(viewer)), true, 'the click enables the preview');
      expect(await until(() => previewEl(e.container, viewer.id) != null), true, 'the container appears');
      const el = previewEl(e.container, viewer.id)!;
      expect(el.style.width, '300px', 'default width 300');
      expect(el.style.height, '300px', 'default height 300');
      expect(await until(() => toggleOf(viewer.id)?.dataset.on === 'true'), true, 'the toggle reads on');
      const nodeEl = e.container.querySelector(`[data-node-id="${viewer.id}"]`) as HTMLElement;
      expect(nodeEl.dataset.inlinePreview, 'true', 'the node carries the preview attribute');
      expect(!!el.querySelector('[data-testid="ff-node-preview-placeholder"]'), true,
        'no captured value yet → placeholder');

      toggleOf(viewer.id)!.click();
      expect(await until(() => !inlinePreviewEnabled(viewer)), true, 'the second click disables');
      expect(await until(() => previewEl(e.container, viewer.id) == null), true, 'the container is gone');
      expect(viewer.properties[INLINE_PREVIEW_PROP] === undefined, true, 'off deletes the property (tidy saves)');
    } finally {
      destroyEditor(e);
    }
  });

  test('the node context menu offers Show/Hide preview on node', async () => {
    const e = makeEditor();
    try {
      const viewer = await addNode(e.flow, 'Viewers/Scatter Plot');
      await until(() => !!e.container.querySelector('.ff-node'));
      const nodeEl = e.container.querySelector(`[data-node-id="${viewer.id}"]`) as HTMLElement;
      nodeEl.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, clientX: 10, clientY: 10}));
      const item = (label: string): HTMLElement | null =>
        Array.from(document.querySelectorAll<HTMLElement>('.d4-menu-item-label'))
          .find((el) => el.textContent?.trim() === label) ?? null;
      expect(await until(() => item('Show preview on node') != null), true, 'the menu offers the preview');
      item('Show preview on node')!.click();
      expect(await until(() => inlinePreviewEnabled(viewer)), true, 'the menu item enables it');
    } finally {
      for (const el of Array.from(document.querySelectorAll('.d4-menu-popup, .d4-menu-dropdown')))
        el.remove();
      destroyEditor(e);
    }
  });

  test('the container mounts the captured live root; toggling off releases it', async () => {
    const container = ui.div([], {style: {width: '1000px', height: '700px', position: 'absolute', left: '-10000px'}});
    document.body.appendChild(container);
    const fakeRoot = ui.divText('live-content');
    let content: HTMLElement | null = fakeRoot;
    const flow = new FlowEditor(container, {getInlinePreviewContent: () => content});
    try {
      const viewer = createNode('Viewers/Scatter Plot')!;
      await flow.addNodeAt(viewer, 0, 0);
      await flow.setInlinePreview(viewer.id, true);

      expect(await until(() => previewEl(container, viewer.id)?.contains(fakeRoot) === true), true,
        'the live root is mounted inside the node');
      expect(fakeRoot.dataset[INLINE_HOSTED_DATA_KEY], 'true', 'the root is marked as node-hosted');

      // A collapsed node folds the preview away and releases the root.
      await flow.toggleCollapsed(viewer.id);
      expect(await until(() => fakeRoot.dataset[INLINE_HOSTED_DATA_KEY] === undefined), true,
        'collapse releases the hosted marker');
      await flow.toggleCollapsed(viewer.id);
      expect(await until(() => fakeRoot.dataset[INLINE_HOSTED_DATA_KEY] === 'true'), true,
        'expand claims it back');

      await flow.setInlinePreview(viewer.id, false);
      expect(await until(() => fakeRoot.dataset[INLINE_HOSTED_DATA_KEY] === undefined), true,
        'toggle-off releases the hosted marker');
      expect(previewEl(container, viewer.id), null, 'the container is gone');

      // No content → placeholder, and stale content swaps out cleanly.
      content = null;
      await flow.setInlinePreview(viewer.id, true);
      expect(await until(() =>
        !!previewEl(container, viewer.id)?.querySelector('[data-testid="ff-node-preview-placeholder"]')), true,
      'no captured value → placeholder');
    } finally {
      flow.destroy();
      container.remove();
    }
  });

  test('a drag-resize persists the size into the node properties', async () => {
    const e = makeEditor();
    try {
      const viewer = await addNode(e.flow, 'Viewers/Scatter Plot');
      await e.flow.setInlinePreview(viewer.id, true);
      expect(await until(() => previewEl(e.container, viewer.id) != null), true, 'container mounted');
      const el = previewEl(e.container, viewer.id)!;
      // CSS `resize: both` writes inline width/height — simulate the result.
      el.style.width = '400px';
      el.style.height = '350px';
      expect(await until(() =>
        inlinePreviewSize(viewer).width === 400 && inlinePreviewSize(viewer).height === 350), true,
      'the ResizeObserver persisted the new size');
    } finally {
      destroyEditor(e);
    }
  });

  test('layout estimates account for the preview size', async () => {
    const viewer = createNode('Viewers/Scatter Plot')!;
    const hOff = estimateNodeHeight(viewer);
    const wOff = estimateNodeWidth(viewer);
    viewer.properties[INLINE_PREVIEW_PROP] = true;
    viewer.properties[INLINE_PREVIEW_SIZE_PROP] = {width: 500, height: 400};
    expect(estimateNodeHeight(viewer) >= hOff + 400, true, 'height grows by the preview height');
    expect(estimateNodeWidth(viewer) >= 500, true, 'width grows past the 280 cap');
    viewer.collapsed = true;
    expect(estimateNodeHeight(viewer), 30, 'collapsed stays a title bar');
    expect(estimateNodeWidth(viewer) <= 280, true, 'collapsed width ignores the preview');
    expect(wOff <= 280, true, 'preview-off width keeps the CSS cap');
  });

  test('the bottom panel shows a note instead of stealing a node-hosted root', async () => {
    const root = ui.divText('live-content');
    const summary = {type: 'viewer' as const, value: {root}};

    root.dataset[INLINE_HOSTED_DATA_KEY] = 'true';
    const note = buildPreview('result', summary);
    expect(note != null, true, 'a hosted root still yields a block');
    expect(note!.dataset.testid, 'ff-preview-inline-note', 'the block is the note');
    expect(note!.contains(root), false, 'the live root is NOT stolen from the node');

    delete root.dataset[INLINE_HOSTED_DATA_KEY];
    const mounted = buildPreview('result', summary);
    expect(!!mounted && mounted.contains(root), true, 'an unhosted root mounts in the panel as before');
  });

  test('OutputPreviewPanel.refresh rebuilds the shown preview in place', async () => {
    const panel = new OutputPreviewPanel();
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('x', ['a', 'b'])]);
    const state = {
      status: NodeExecStatus.completed,
      outputs: {result: {type: 'dataframe' as const, rows: 2, cols: 1, clone: df}},
    };
    const contentEl = (): Element | null =>
      panel.root.querySelector('[data-testid="ff-output-panel-content"]')?.firstElementChild ?? null;

    panel.refresh();
    expect(panel.panelState, 'hidden', 'refresh on an empty panel is a no-op');

    panel.showForNode({id: 'n1', label: 'a'}, state);
    const first = contentEl();
    panel.refresh();
    const rebuilt = contentEl();
    expect(!!rebuilt && rebuilt !== first, true, 'refresh rebuilds even with the same state');
    expect(panel.root.querySelector('[data-testid="ff-output-panel-node"]')!.textContent, 'a',
      'the node label is kept');
  });

  test('controller: inlinePreviewRoot finds the captured root; ownership follows the toggle', async () => {
    const e = makeEditor();
    try {
      const viewer = await addNode(e.flow, 'Viewers/Scatter Plot');
      const ctrl = new ExecutionController(e.flow);
      const root = ui.divText('live-content');

      expect(ctrl.inlinePreviewRoot(viewer.id), null, 'no run → no root');
      ctrl.state.setNodeStatus(viewer.id, NodeExecStatus.completed, {
        outputs: {viewer: {type: 'viewer', value: {root}}},
      });
      expect(ctrl.inlinePreviewRoot(viewer.id) === root, true, 'the captured root is found');

      ctrl.syncInlinePreviewOwnership(viewer.id);
      expect(root.dataset[INLINE_HOSTED_DATA_KEY] === undefined, true, 'preview off → not claimed');

      viewer.properties[INLINE_PREVIEW_PROP] = true;
      ctrl.syncInlinePreviewOwnership(viewer.id);
      expect(root.dataset[INLINE_HOSTED_DATA_KEY], 'true', 'preview on → claimed for the node');

      viewer.collapsed = true;
      ctrl.syncInlinePreviewOwnership(viewer.id);
      expect(root.dataset[INLINE_HOSTED_DATA_KEY] === undefined, true, 'collapsed → released to the panel');

      // Stale keeps the last result visible on the node, like the bottom panel.
      viewer.collapsed = false;
      ctrl.state.markStale([viewer.id]);
      expect(ctrl.inlinePreviewRoot(viewer.id) === root, true, 'a stale result still previews');
    } finally {
      destroyEditor(e);
    }
  });

  test('controller: a graphics output builds a cached element (SVG and PNG)', async () => {
    const e = makeEditor();
    try {
      const node = outputNodeOf('graphics');
      await e.flow.addNodeAt(node, 0, 0);
      const ctrl = new ExecutionController(e.flow);

      ctrl.state.setNodeStatus(node.id, NodeExecStatus.completed, {
        outputs: {charges: {type: 'graphics', value: SVG_SAMPLE}},
      });
      const el = ctrl.inlinePreviewRoot(node.id);
      expect(el != null, true, 'an element is built from the SVG');
      expect(!!el!.querySelector('svg'), true, 'the SVG markup is rendered inline');
      expect(ctrl.inlinePreviewRoot(node.id) === el, true,
        'the same captured summary returns the same element (React re-attach, no rebuild)');

      ctrl.state.setNodeStatus(node.id, NodeExecStatus.completed, {
        outputs: {charges: {type: 'graphics', value: 'aGVsbG8='}},
      });
      const png = ctrl.inlinePreviewRoot(node.id);
      expect(png !== el, true, 'a new run builds a fresh element');
      expect(png!.style.backgroundImage.includes('data:image/png;base64'), true,
        'a non-SVG payload renders as a base64 PNG');

      // Graphics is copied data — ownership stamping never touches it, so the
      // bottom panel keeps its own copy alongside the node's.
      node.properties[INLINE_PREVIEW_PROP] = true;
      ctrl.syncInlinePreviewOwnership(node.id);
      expect(png!.dataset[INLINE_HOSTED_DATA_KEY] === undefined, true,
        'no hosted marker on a graphics element');
      const panelBlock = buildPreview('charges',
        ctrl.state.getNodeState(node.id)!.outputs!['charges']);
      expect(panelBlock != null, true, 'the bottom panel still renders graphics normally');
      expect(panelBlock!.dataset.testid !== 'ff-preview-inline-note', true,
        'no yield-note for graphics — both previews can coexist');
    } finally {
      destroyEditor(e);
    }
  });

  test('a graphics node mounts its rendered image inside the node', async () => {
    const e = makeWiredEditor();
    try {
      const node = outputNodeOf('graphics');
      await e.flow.addNodeAt(node, 0, 0);
      e.ctrl.state.setNodeStatus(node.id, NodeExecStatus.completed, {
        outputs: {charges: {type: 'graphics', value: SVG_SAMPLE}},
      });
      await e.flow.setInlinePreview(node.id, true);
      expect(await until(() =>
        !!previewEl(e.container, node.id)?.querySelector('svg')), true,
      'the graphics element is mounted inside the node container');
      expect(previewEl(e.container, node.id)!.contains(e.ctrl.inlinePreviewRoot(node.id)!), true,
        'it is the controller-cached element');
    } finally {
      destroyEditor(e);
    }
  });

  test('a live Gasteiger Partial Charges run renders its graphics in the node (live)', async () => {
    const typeName = funcTypeName('ChemistryGasteigerPartialCharges');
    if (!typeName) return; // Chem not installed on this stand
    // Chem ships TestPythonRunning "to be used in tests to ensure JKG is up and
    // running" — when the compute VM is down, a script call hangs into a 504
    // minutes later, so probe first and skip with the cause named.
    const probe = await Promise.race([
      grok.functions.call('Chem:TestPythonRunning', {x: 2, y: 3}).then((v) => v, () => null),
      new Promise<'timeout'>((r) => setTimeout(() => r('timeout'), 30000)),
    ]);
    if (probe !== 5) {
      console.warn('Flow: inline preview: skipping the Gasteiger live run — ' +
        `Python scripting is not answering on this stand (probe: ${String(probe)})`);
      return;
    }
    const e = makeWiredEditor();
    try {
      const node = await addNode(e.flow, typeName);
      // The mol/contours inputs carry declared defaults, so the node is ready as dropped.
      expect(supportsInlinePreview(node), true, 'the graphics-outputting func is preview-capable');
      await e.flow.setInlinePreview(node.id, true);

      expect(e.ctrl.runAutorun(new Set(), SETTINGS), 'started', 'the run starts');
      await until(() => {
        const s = e.ctrl.state.getNodeState(node.id)?.status;
        return s === NodeExecStatus.completed || s === NodeExecStatus.errored;
      }, 90000);
      const st = e.ctrl.state.getNodeState(node.id);
      expect(st?.status, NodeExecStatus.completed,
        `the Python script completed (status=${st?.status ?? 'none'}, error=${st?.error ?? ''})`);

      const root = (): HTMLElement | null => e.ctrl.inlinePreviewRoot(node.id);
      expect(await until(() => root() != null, 5000), true, 'a graphics element was built');
      expect(await until(() => previewEl(e.container, node.id)?.contains(root()!) === true, 5000), true,
        'the charges image is mounted inside the node');
      const el = root()!;
      expect(!!el.querySelector('svg') || el.style.backgroundImage.includes('data:image'), true,
        'the element carries a real rendered image');
    } finally {
      destroyEditor(e);
    }
  }, {timeout: 120000});

  test('a live run mounts the viewer inside the node and the panel yields (live)', async () => {
    const e = makeWiredEditor();
    const df = DG.DataFrame.fromCsv('x,y\n1,2\n3,4\n5,6\n');
    df.name = 'ffInlinePreviewLive';
    const shellTable = grok.shell.addTable(df);
    try {
      const sel = await addNode(e.flow, 'Utilities/Select Table');
      sel.properties['tableName'] = 'ffInlinePreviewLive';
      const viewer = await addNode(e.flow, 'Viewers/Scatter Plot', 300, 0);
      await e.flow.addConnectionByKeys(sel.id, 'table', viewer.id, 'table');
      await e.flow.setInlinePreview(viewer.id, true);

      expect(e.ctrl.runAutorun(new Set(), SETTINGS), 'started', 'the run starts');
      expect(await until(() =>
        e.ctrl.state.getNodeState(viewer.id)?.status === NodeExecStatus.completed, 15000), true,
      'the viewer node completed');

      const root = (): HTMLElement | null => e.ctrl.inlinePreviewRoot(viewer.id);
      expect(await until(() => root() != null, 5000), true, 'a live viewer root was captured');
      expect(await until(() => previewEl(e.container, viewer.id)?.contains(root()!) === true, 5000), true,
        'the live viewer is mounted inside the node');
      expect(root()!.dataset[INLINE_HOSTED_DATA_KEY], 'true', 'the root is claimed by the node');

      const block = buildPreview('viewer',
        e.ctrl.state.getNodeState(viewer.id)!.outputs!['viewer']);
      expect(block?.dataset.testid, 'ff-preview-inline-note',
        'the bottom panel yields with a note while the node hosts the viewer');

      await e.flow.setInlinePreview(viewer.id, false);
      expect(await until(() => root()!.dataset[INLINE_HOSTED_DATA_KEY] === undefined), true,
        'turning the preview off hands the root back to the panel');
    } finally {
      try {
        grok.shell.closeTable(shellTable);
      } catch {/* best effort */}
      destroyEditor(e);
    }
  }, {timeout: 60000});
});
