/** In-node viewer/widget preview: the title-bar toggle, the resizable container
 *  that mounts the captured live root, serialization of the choice/size, and the
 *  bottom-panel ownership handoff. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {ClassicPreset} from 'rete';
import {category, test, expect, expectFloat, before} from '@datagrok-libraries/utils/src/test';

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
import {FuncFlowView} from '../funcflow-view';

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
function makeWiredEditor(opts: {onScreen?: boolean} = {}): TestEditor & {ctrl: ExecutionController} {
  const container = ui.div([], {style: opts.onScreen ?
    {width: '1000px', height: '700px', position: 'fixed', left: '0', top: '0',
      zIndex: '5000', background: '#fff'} :
    {width: '1000px', height: '700px', position: 'absolute', left: '-10000px'}});
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

/** The screen-space portal actually hosting a node's preview content. */
function portalEl(container: HTMLElement, nodeId: string): HTMLElement | null {
  return container.querySelector(`.ff-node-preview-portal[data-node-id="${nodeId}"]`);
}

/** Last measured portal/host rects — appended to failure messages. */
let lastPortalGeom = '';

/** True once the portal's on-screen box matches the in-card container's. */
function portalAligned(container: HTMLElement, nodeId: string): boolean {
  const portal = portalEl(container, nodeId);
  const host = previewEl(container, nodeId);
  if (!portal || !host) {
    lastPortalGeom = `portal=${!!portal} host=${!!host}`;
    return false;
  }
  const p = portal.getBoundingClientRect();
  const h = host.getBoundingClientRect();
  lastPortalGeom = `portal=[${[p.left, p.top, p.width, p.height].map(Math.round)}] ` +
    `host=[${[h.left, h.top, h.width, h.height].map(Math.round)}] ` +
    `style=[${portal.style.left},${portal.style.top},${portal.style.width},${portal.style.display}]`;
  return Math.abs(p.left - h.left) < 2 && Math.abs(p.top - h.top) < 2 &&
    Math.abs(p.width - h.width) < 2 && Math.abs(p.height - h.height) < 2;
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

      // The content mounts once into the screen-space portal and never moves
      // again — no re-mount or resize on hover, nothing visibly "happens".
      expect(await until(() => portalEl(container, viewer.id)?.contains(fakeRoot) === true), true,
        'the live root is mounted in the portal');
      expect(fakeRoot.closest('.ff-node'), null,
        'the content has NO transformed node ancestor (fixed-position popups stay true)');
      expect(await until(() => portalAligned(container, viewer.id)), true,
        `the portal covers the in-card container — ${lastPortalGeom}`);
      expect(fakeRoot.dataset[INLINE_HOSTED_DATA_KEY], 'true', 'the root is marked as node-hosted');
      expect(!!portalEl(container, viewer.id)!.querySelector('[data-testid="ff-node-preview-grip"]'), true,
        'the portal carries the resize grip');
      const mountedPortal = portalEl(container, viewer.id);
      previewEl(container, viewer.id)!.dispatchEvent(new PointerEvent('pointerenter', {buttons: 0}));
      mountedPortal!.dispatchEvent(new PointerEvent('pointerenter', {buttons: 0}));
      await new Promise((r) => setTimeout(r, 200));
      expect(portalEl(container, viewer.id) === mountedPortal, true,
        'hover changes nothing — same portal, no re-mount');
      expect(mountedPortal!.contains(fakeRoot), true, 'the content never moved');

      // A collapsed node folds the preview away and releases the root.
      await flow.toggleCollapsed(viewer.id);
      expect(await until(() => fakeRoot.dataset[INLINE_HOSTED_DATA_KEY] === undefined), true,
        'collapse releases the hosted marker');
      expect(portalEl(container, viewer.id), null, 'collapse removes the portal');
      await flow.toggleCollapsed(viewer.id);
      expect(await until(() => fakeRoot.dataset[INLINE_HOSTED_DATA_KEY] === 'true'), true,
        'expand claims it back');

      await flow.setInlinePreview(viewer.id, false);
      expect(await until(() => fakeRoot.dataset[INLINE_HOSTED_DATA_KEY] === undefined), true,
        'toggle-off releases the hosted marker');
      expect(portalEl(container, viewer.id), null, 'toggle-off removes the portal');
      expect(previewEl(container, viewer.id), null, 'the container is gone');

      // No content → placeholder, and stale content swaps out cleanly.
      content = null;
      await flow.setInlinePreview(viewer.id, true);
      expect(await until(() =>
        !!previewEl(container, viewer.id)?.querySelector('[data-testid="ff-node-preview-placeholder"]')), true,
      'no captured value → placeholder');
      expect(portalEl(container, viewer.id), null, 'no content → no portal');
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
        !!portalEl(e.container, node.id)?.querySelector('svg')), true,
      'the graphics element is mounted in the node preview portal');
      expect(portalEl(e.container, node.id)!.contains(e.ctrl.inlinePreviewRoot(node.id)!), true,
        'it is the controller-cached element');
    } finally {
      destroyEditor(e);
    }
  });

  test('overlapping cards interleave with previews: portal order + clip holes', async () => {
    // The reported artifact: previews painted above every card. Portals now
    // follow their nodes' paint order, and a card ABOVE a portal's node cuts a
    // click-through hole into it — visual and interactive interleaving.
    const e = makeWiredEditor({onScreen: true});
    try {
      const a = outputNodeOf('graphics');
      const b = outputNodeOf('graphics');
      await e.flow.addNodeAt(a, 20, 20);
      // B is added later → painted above A; place its card over A's preview box.
      await e.flow.addNodeAt(b, 180, 120);
      for (const n of [a, b]) {
        e.ctrl.state.setNodeStatus(n.id, NodeExecStatus.completed, {
          outputs: {charges: {type: 'graphics', value: SVG_SAMPLE}},
        });
        await e.flow.setInlinePreview(n.id, true);
      }
      expect(await until(() => portalAligned(e.container, a.id) && portalAligned(e.container, b.id)),
        true, 'both portals mounted and tracking');

      const pa = portalEl(e.container, a.id)!;
      const pb = portalEl(e.container, b.id)!;
      // eslint-disable-next-line no-bitwise
      expect((pa.compareDocumentPosition(pb) & Node.DOCUMENT_POSITION_FOLLOWING) !== 0, true,
        'portal order follows node paint order (B above A)');

      // A's portal is cut where B's card covers it; B's portal is uncut.
      const diag = (): string => {
        const aCardEl = e.container.querySelector(`[data-node-id="${a.id}"]`) as HTMLElement;
        const wrap = aCardEl?.parentElement;
        const sibs: string[] = [];
        for (let sib = wrap?.nextElementSibling; sib; sib = sib.nextElementSibling) {
          sibs.push(`[${sib.className.toString().slice(0, 20)}|hasCard=${
            !!sib.querySelector(':scope > .ff-node')}]`);
        }
        return `wrapParent=${wrap?.parentElement?.className?.toString()?.slice(0, 30)} ` +
          `sibsAfterA=${sibs.join(',')} clipA='${pa.style.clipPath.slice(0, 60)}'`;
      };
      expect(await until(() => pa.style.clipPath.includes('evenodd')), true,
        `the covered portal carries an evenodd clip with a hole — ${diag()}`);
      expect(pb.style.clipPath, '', 'the top portal is not clipped');

      // The hole is real: a point inside the overlap hits B (card or portal),
      // never A's portal — clicks reach the covering card.
      const bCard = e.container.querySelector(`[data-node-id="${b.id}"]`) as HTMLElement;
      const br = bCard.getBoundingClientRect();
      const ar = pa.getBoundingClientRect();
      const x = Math.max(br.left, ar.left) + 8;
      const y = Math.max(br.top, ar.top) + 8;
      expect(x < Math.min(br.right, ar.right) && y < Math.min(br.bottom, ar.bottom), true,
        'the test geometry actually overlaps');
      const hit = document.elementFromPoint(x, y);
      expect(!!hit && !pa.contains(hit), true,
        `the overlap point hits the covering card, not the clipped portal (hit=${
          (hit as HTMLElement)?.className?.toString()?.slice(0, 40)})`);
      expect(!!hit && (bCard.contains(hit) || pb.contains(hit)), true,
        'what it hits belongs to the covering node');
    } finally {
      destroyEditor(e);
    }
  });

  test('repro + fix: fixed-position popups break in the canvas, stay true in the portal', async () => {
    const e = makeWiredEditor();
    try {
      const node = outputNodeOf('graphics');
      await e.flow.addNodeAt(node, 120, 60);
      e.ctrl.state.setNodeStatus(node.id, NodeExecStatus.completed, {
        outputs: {charges: {type: 'graphics', value: SVG_SAMPLE}},
      });
      await e.flow.setInlinePreview(node.id, true);
      expect(await until(() => portalEl(e.container, node.id) != null), true, 'portal mounted');

      // A viewer's axis/color selector popup is `position: fixed` at viewport
      // coordinates — this probe stands in for it.
      const probeAt = (parent: Element): DOMRect => {
        const probe = document.createElement('div');
        probe.style.cssText = 'position:fixed;left:50px;top:60px;width:10px;height:10px;';
        parent.appendChild(probe);
        const r = probe.getBoundingClientRect();
        probe.remove();
        return r;
      };
      for (const k of [1, 2]) {
        e.flow.setZoom(k);
        expect(await until(() => portalAligned(e.container, node.id), 3000), true,
          `the portal tracks the container at k=${k} — ${lastPortalGeom}`);
        // REPRO of the bug: inside the transformed canvas subtree, the transform
        // becomes the fixed-position containing block — the "viewport" popup
        // lands relative to the node instead, far from where DG placed it.
        const nodeBody = e.container.querySelector(
          `[data-node-id="${node.id}"] .ff-node-body`) as HTMLElement;
        const broken = probeAt(nodeBody);
        expect(Math.abs(broken.left - 50) > 20 || Math.abs(broken.top - 60) > 20, true,
          `k=${k}: a fixed popup inside the node lands far from its viewport coords ` +
          `(got ${Math.round(broken.left)},${Math.round(broken.top)} for 50,60) — the reported bug`);
        // THE FIX: the portal has no transformed ancestor, so fixed coordinates
        // are true viewport coordinates — popups land where DG puts them.
        const fixedOk = probeAt(portalEl(e.container, node.id)!);
        expectFloat(fixedOk.left, 50, 1.5, `k=${k}: a fixed popup in the portal lands at its coords (x)`);
        expectFloat(fixedOk.top, 60, 1.5, `k=${k}: ...and (y)`);
        expectFloat(fixedOk.width, 10, 1.5, `k=${k}: ...unscaled`);
        // And the content itself renders untransformed — canvas hit-tests true.
        const content = e.ctrl.inlinePreviewRoot(node.id)!;
        expectFloat(content.getBoundingClientRect().width, content.clientWidth, 2,
          `k=${k}: content client px == layout px`);
      }
      e.flow.setZoom(1);
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
      grok.functions.call('Chem:TestPythonRunning', {x: 2, y: 3}).then((v) => v as number, () => null),
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
      expect(await until(() => portalEl(e.container, node.id)?.contains(root()!) === true, 5000), true,
        'the charges image is mounted in the node preview portal');
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
      expect(await until(() => portalEl(e.container, viewer.id)?.contains(root()!) === true, 5000), true,
        'the live viewer is mounted in the node preview portal');
      expect(root()!.dataset[INLINE_HOSTED_DATA_KEY], 'true', 'the root is claimed by the node');

      const block = buildPreview('viewer',
        e.ctrl.state.getNodeState(viewer.id)!.outputs!['viewer']);
      expect(block?.dataset.testid, 'ff-preview-inline-note',
        'the bottom panel yields with a note while the node hosts the viewer');

      // The scatter plot hit-tests on its own canvas in client px — in the
      // portal the REAL viewer renders untransformed at every zoom.
      for (const k of [0.5, 2]) {
        e.flow.setZoom(k);
        expect(await until(() => portalAligned(e.container, viewer.id), 3000), true,
          `the portal re-tracks the container at k=${k} — ${lastPortalGeom}`);
        const r = root()!;
        expectFloat(r.getBoundingClientRect().width, r.clientWidth, 2,
          `the live scatter plot's client geometry is unscaled at k=${k}`);
        expectFloat(r.getBoundingClientRect().height, r.clientHeight, 2,
          `...vertically too at k=${k}`);
      }
      e.flow.setZoom(1);
      await until(() => portalAligned(e.container, viewer.id), 3000);

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

  test('right-click inside the preview reaches the viewer, not the node menu', async () => {
    const container = ui.div([], {style: {width: '1000px', height: '700px', position: 'absolute', left: '-10000px'}});
    document.body.appendChild(container);
    const fakeRoot = ui.divText('live-content');
    const flow = new FlowEditor(container, {getInlinePreviewContent: () => fakeRoot});
    const flowMenuItem = (label: string): HTMLElement | null =>
      Array.from(document.querySelectorAll<HTMLElement>('.d4-menu-item-label'))
        .find((el) => el.textContent?.trim() === label) ?? null;
    try {
      const viewer = createNode('Viewers/Scatter Plot')!;
      await flow.addNodeAt(viewer, 0, 0);
      await flow.setInlinePreview(viewer.id, true);
      expect(await until(() => portalEl(container, viewer.id)?.contains(fakeRoot) === true), true,
        'content mounted in the portal');

      // Right-click on the content: the portal swallows it before the canvas —
      // no node menu, and the event is left for the viewer's own menu.
      const ev = new MouseEvent('contextmenu', {bubbles: true, cancelable: true, clientX: 5, clientY: 5});
      fakeRoot.dispatchEvent(ev);
      await new Promise((r) => setTimeout(r, 300));
      expect(flowMenuItem('Duplicate'), null, 'no node context menu over the preview content');
      expect(ev.defaultPrevented, false, 'the event is not defaultPrevented — the viewer may handle it');

      // Right-click on the node title still opens the node menu.
      const title = container.querySelector(
        `[data-node-id="${viewer.id}"] [data-testid="ff-node-title"]`) as HTMLElement;
      title.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, clientX: 5, clientY: 5}));
      expect(await until(() => flowMenuItem('Duplicate') != null), true, 'the node menu still works on the card');
    } finally {
      for (const el of Array.from(document.querySelectorAll('.d4-menu-popup, .d4-menu-dropdown')))
        el.remove();
      flow.destroy();
      container.remove();
    }
  });

  test('axis selector popup opens AT the viewer; right-click gets the viewer menu (live)', async () => {
    // Faithful environment: the Flow view attached to the shell, on screen.
    const view = new FuncFlowView();
    grok.shell.addView(view);
    const df = DG.DataFrame.fromCsv('height,weight,age\n160,60,30\n170,70,40\n180,80,50\n175,72,35\n');
    df.name = 'ffInlinePopupLive';
    const shellTable = grok.shell.addTable(df);
    grok.shell.v = view;
    try {
      await until(() => (view as any).flow != null, 10000);
      const flow = (view as any).flow as FlowEditor;
      const ctrl = (view as any).executionController as ExecutionController;
      const sel = await addNode(flow, 'Utilities/Select Table');
      sel.properties['tableName'] = 'ffInlinePopupLive';
      const viewer = await addNode(flow, 'Viewers/Scatter Plot', 300, 0);
      await flow.addConnectionByKeys(sel.id, 'table', viewer.id, 'table');
      await flow.setInlinePreview(viewer.id, true);
      expect(ctrl.runAutorun(new Set(), SETTINGS), 'started', 'the run starts');
      await until(() => ctrl.state.getNodeState(viewer.id)?.status === NodeExecStatus.completed, 15000);
      const root = (): HTMLElement | null => ctrl.inlinePreviewRoot(viewer.id);
      expect(await until(() => root() != null &&
        !!root()!.querySelector('.d4-column-selector-column'), 8000), true,
      'the live scatter plot with its DOM axis selectors is mounted');
      expect(await until(() => portalEl(view.root, viewer.id)?.contains(root()!) === true, 5000), true,
        'the viewer is mounted in the preview portal');

      // DG may reuse a previously-built popup element on later clicks — remember
      // the ones we've seen and accept a known popup that became visible again.
      const knownPopups = new Set<HTMLElement>();
      const clickSelector = async (): Promise<HTMLElement[]> => {
        const selEl = root()!.querySelector('.d4-column-selector-column') as HTMLElement;
        const r = selEl.getBoundingClientRect();
        const px = r.left + r.width / 2;
        const py = r.top + r.height / 2;
        const added: HTMLElement[] = [];
        const mo = new MutationObserver((muts) => {
          for (const m of muts) {
            for (const n of Array.from(m.addedNodes))
              if (n instanceof HTMLElement) added.push(n);
          }
        });
        mo.observe(document.documentElement, {childList: true, subtree: true});
        for (const t of ['pointermove', 'mousemove', 'pointerdown', 'mousedown', 'pointerup', 'mouseup', 'click']) {
          const Ev = t.startsWith('pointer') ? PointerEvent : MouseEvent;
          selEl.dispatchEvent(new Ev(t, {bubbles: true, cancelable: true, view: window,
            clientX: px, clientY: py, button: 0} as MouseEventInit));
        }
        await new Promise((r2) => setTimeout(r2, 600));
        mo.disconnect();
        const isShown = (el: HTMLElement): boolean => el.isConnected &&
          getComputedStyle(el).position === 'fixed' && el.getBoundingClientRect().width > 10;
        const found = added.filter(isShown);
        for (const el of found) knownPopups.add(el);
        return found.length > 0 ? found : Array.from(knownPopups).filter(isShown);
      };
      const assertNear = (popup: HTMLElement, label: string): void => {
        const p = popup.getBoundingClientRect();
        const v = root()!.getBoundingClientRect();
        const near = p.left > v.left - 200 && p.left < v.right + 200 &&
          p.top > v.top - 200 && p.top < v.bottom + 200;
        expect(near, true, `${label}: the popup [${Math.round(p.left)},${Math.round(p.top)}] opens near ` +
          `the viewer [${Math.round(v.left)},${Math.round(v.top)},${Math.round(v.width)}x${Math.round(v.height)}]`);
        const center = document.elementFromPoint(p.left + p.width / 2, p.top + p.height / 2);
        expect(!!center && (popup === center || popup.contains(center)), true,
          `${label}: the popup is actually visible (nothing covers it)`);
      };

      // A click can TOGGLE a popup Escape failed to dismiss — retry once.
      const openSelectorPopup = async (): Promise<HTMLElement[]> => {
        const first = await clickSelector();
        return first.length > 0 ? first : clickSelector();
      };

      // Zoom 1: the popup used to land ~200px off (transformed containing block).
      let popups = await openSelectorPopup();
      expect(popups.length > 0, true, 'clicking the axis selector opens its popup');
      assertNear(popups[0], 'k=1');
      document.body.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
      await new Promise((r) => setTimeout(r, 300));

      // Zoomed: same guarantee — the portal never scales or displaces the viewer.
      flow.setZoom(1.5);
      await until(() => portalAligned(view.root, viewer.id), 3000);
      popups = await openSelectorPopup();
      expect(popups.length > 0, true, 'the selector still opens its popup at zoom 1.5');
      assertNear(popups[0], 'k=1.5');
      document.body.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
      await new Promise((r) => setTimeout(r, 300));
      flow.setZoom(1);

      // Right-click on the plot: the VIEWER's menu, not the node's.
      const canvas = root()!.querySelector('canvas') as HTMLElement;
      const cr = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, view: window,
        clientX: cr.left + cr.width / 2, clientY: cr.top + cr.height / 2}));
      const menuLabels = (): string[] => Array.from(
        document.querySelectorAll<HTMLElement>('.d4-menu-item-label'))
        .map((el) => el.textContent?.trim() ?? '');
      expect(await until(() => menuLabels().length > 0, 3000), true,
        'right-clicking the plot opens a context menu (the viewer\'s own)');
      expect(menuLabels().includes('Run up to here & preview'), false,
        'it is NOT the flow node menu');
      document.body.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
    } finally {
      for (const el of Array.from(document.querySelectorAll('.d4-menu-popup, .d4-menu-dropdown')))
        el.remove();
      try {
        grok.shell.closeTable(shellTable);
      } catch {/* best effort */}
      view.close();
    }
  }, {timeout: 90000});
});
