/** Datagrok selection semantics on the canvas (mirrors d4 `selectRows` /
 *  `areaSelector`): plain click selects exclusively, Shift+click adds,
 *  Ctrl+click toggles, Ctrl+Shift+click removes; Shift+drag draws a marquee
 *  that ADDS the covered nodes (Ctrl at mouse-up removes them instead);
 *  Ctrl+A / Ctrl+Shift+A select / deselect everything. Exercised through
 *  real PointerEvent / KeyboardEvent dispatch on the live editor DOM. */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';
import * as ui from 'datagrok-api/ui';

import {registerBuiltinNodes} from '../rete/node-factory';
import {FlowEditor} from '../rete/flow-editor';
import {FlowNode} from '../rete/scheme';
import {makeEditor, destroyEditor, addNode, until, TestEditor} from './test-utils';

function isSelected(e: TestEditor, node: FlowNode): boolean {
  return (e.flow.editor.getNode(node.id) as {selected?: boolean})?.selected === true;
}

function nodeEl(e: TestEditor, node: FlowNode): HTMLElement {
  const el = e.container.querySelector<HTMLElement>(`.ff-node[data-node-id="${node.id}"]`);
  if (!el) throw new Error('node element not rendered');
  return el;
}

function pointer(target: EventTarget, type: string, init: PointerEventInit): void {
  target.dispatchEvent(new PointerEvent(type, {bubbles: true, cancelable: true, button: 0, ...init}));
}

const tick = (ms = 30): Promise<void> => new Promise((res) => setTimeout(res, ms));

/** The rete area holder — where real empty-canvas pointerdowns land. Events
 *  dispatched on the container ITSELF would also fire the AreaPlugin's own
 *  same-element listeners (stopPropagation can't suppress those), starting a
 *  phantom pan whose click-release clears the selection. */
function canvasBg(e: TestEditor): HTMLElement {
  return e.container.firstElementChild as HTMLElement;
}

/** Click a node's title area with the given modifiers. The tick between down
 *  and up lets rete's async `nodepicked` add settle, as a real click would. */
async function clickNode(e: TestEditor, node: FlowNode, init: PointerEventInit = {}): Promise<void> {
  const el = nodeEl(e, node);
  const r = el.getBoundingClientRect();
  const at = {clientX: r.left + r.width / 2, clientY: r.top + 8};
  pointer(el, 'pointerdown', {...at, ...init});
  await tick();
  pointer(el, 'pointerup', {...at, ...init});
  await tick();
}

/** Shift+drag a marquee between two client points on the empty canvas. All
 *  three events target the canvas background (not `window`): real events
 *  reach window listeners through the capture phase, so the editor's
 *  window-capture handlers must run before rete's window-bubble ones. */
async function dragMarquee(e: TestEditor, from: {x: number; y: number}, to: {x: number; y: number},
  upInit: PointerEventInit = {}): Promise<void> {
  pointer(canvasBg(e), 'pointerdown', {clientX: from.x, clientY: from.y, shiftKey: true});
  await tick();
  pointer(canvasBg(e), 'pointermove', {clientX: to.x, clientY: to.y, shiftKey: true});
  await tick();
  pointer(canvasBg(e), 'pointerup', {clientX: to.x, clientY: to.y, shiftKey: true, ...upInit});
  await tick();
}

/** Client-space rectangle that covers the given nodes with a margin. */
function marqueeOver(e: TestEditor, nodes: FlowNode[]): {from: {x: number; y: number}; to: {x: number; y: number}} {
  const rects = nodes.map((n) => nodeEl(e, n).getBoundingClientRect());
  return {
    from: {x: Math.min(...rects.map((r) => r.left)) - 10, y: Math.min(...rects.map((r) => r.top)) - 10},
    to: {x: Math.max(...rects.map((r) => r.right)) + 10, y: Math.max(...rects.map((r) => r.bottom)) + 10},
  };
}

async function threeNodes(e: TestEditor): Promise<[FlowNode, FlowNode, FlowNode]> {
  const a = await addNode(e.flow, 'Inputs/String Input', 0, 0);
  const b = await addNode(e.flow, 'Inputs/String Input', 250, 0);
  const c = await addNode(e.flow, 'Inputs/String Input', 0, 450);
  await until(() => !!e.container.querySelector(`.ff-node[data-node-id="${c.id}"]`));
  return [a, b, c];
}

category('Flow: selection', () => {
  before(async () => {
    registerBuiltinNodes();
  });

  test('re-clicking or grabbing the already-selected node fires onNodeSelected once', async () => {
    const container = ui.div([], {style: {width: '1000px', height: '700px', position: 'absolute', left: '-10000px'}});
    document.body.appendChild(container);
    let fires = 0;
    const flow = new FlowEditor(container, {onNodeSelected: () => fires++});
    const e: TestEditor = {flow, container};
    try {
      const a = await addNode(flow, 'Inputs/String Input', 0, 0);
      await until(() => !!container.querySelector(`.ff-node[data-node-id="${a.id}"]`));

      await clickNode(e, a);
      expect(await until(() => isSelected(e, a)), true, 'first click selects');
      expect(fires, 1, 'first click fires once');

      // Re-clicks and grabs of the current node change nothing — no re-fire,
      // so the host never rebuilds its panels for a no-op.
      await clickNode(e, a);
      await clickNode(e, a);
      expect(fires, 1, 're-clicks do not re-fire');

      // Deselect all, then click again → a real change fires again.
      await flow.unselectAllNodes();
      expect(await until(() => !isSelected(e, a)), true, 'deselected');
      await clickNode(e, a);
      expect(await until(() => isSelected(e, a)), true, 're-selected after deselect-all');
      expect(fires, 2, 'a real selection change fires again');
    } finally {
      flow.destroy();
      container.remove();
    }
  });

  test('shift+drag marquee adds the covered nodes to the selection', async () => {
    const e = makeEditor();
    try {
      const [a, b, c] = await threeNodes(e);
      await e.flow.selectNode(c.id); // pre-selected, outside the marquee
      const {from, to} = marqueeOver(e, [a, b]);
      await dragMarquee(e, from, to);
      expect(await until(() => isSelected(e, a) && isSelected(e, b)), true, 'a and b selected');
      expect(isSelected(e, c), true, 'existing selection kept (additive, never replaces)');
    } finally {
      destroyEditor(e);
    }
  });

  test('ctrl at mouse-up turns the marquee into remove-from-selection', async () => {
    const e = makeEditor();
    try {
      const [a, b, c] = await threeNodes(e);
      for (const n of [a, b, c]) await e.flow.selectNode(n.id, true);
      const {from, to} = marqueeOver(e, [a, b]);
      await dragMarquee(e, from, to, {ctrlKey: true});
      expect(await until(() => !isSelected(e, a) && !isSelected(e, b)), true, 'a and b removed');
      expect(isSelected(e, c), true, 'node outside the marquee stays selected');
    } finally {
      destroyEditor(e);
    }
  });

  test('ctrl+drag does not start a marquee anymore', async () => {
    const e = makeEditor();
    try {
      await threeNodes(e);
      pointer(canvasBg(e), 'pointerdown', {clientX: -9990, clientY: 10, ctrlKey: true});
      pointer(canvasBg(e), 'pointermove', {clientX: -9600, clientY: 300, ctrlKey: true});
      expect(e.container.querySelector('.ff-rect-select') == null, true, 'no marquee element');
      pointer(canvasBg(e), 'pointerup', {clientX: -9600, clientY: 300, ctrlKey: true});
    } finally {
      destroyEditor(e);
    }
  });

  test('click modifiers mirror selectRows: exclusive / add / toggle / remove', async () => {
    const e = makeEditor();
    try {
      const [a, b] = await threeNodes(e);
      await clickNode(e, a);
      expect(await until(() => isSelected(e, a)), true, 'plain click selects');

      await clickNode(e, b, {shiftKey: true});
      expect(await until(() => isSelected(e, b)), true, 'shift+click adds');
      expect(isSelected(e, a), true, 'shift+click keeps the selection');

      await clickNode(e, b, {ctrlKey: true});
      expect(await until(() => !isSelected(e, b)), true, 'ctrl+click toggles off');
      expect(isSelected(e, a), true, 'the rest of the selection survives');

      await clickNode(e, b, {ctrlKey: true});
      expect(await until(() => isSelected(e, b)), true, 'ctrl+click toggles back on');

      await clickNode(e, a, {ctrlKey: true, shiftKey: true});
      expect(await until(() => !isSelected(e, a)), true, 'ctrl+shift+click removes');
      expect(isSelected(e, b), true, 'other nodes untouched');
    } finally {
      destroyEditor(e);
    }
  });

  test('plain click on a multi-selected node collapses the selection on release', async () => {
    const e = makeEditor();
    try {
      const [a, b, c] = await threeNodes(e);
      for (const n of [a, b, c]) await e.flow.selectNode(n.id, true);
      await clickNode(e, a);
      expect(await until(() => !isSelected(e, b) && !isSelected(e, c)), true, 'others deselected');
      expect(isSelected(e, a), true, 'clicked node stays selected');
    } finally {
      destroyEditor(e);
    }
  });

  test('marquee, Ctrl+A, and modifier clicks report onSelectionChanged', async () => {
    // The marquee's release is swallowed (stopImmediatePropagation) before any
    // container/window listener sees it — hosts that track the selection (the
    // suggestion pane) rely on this callback instead.
    const container = ui.div([], {style: {width: '1000px', height: '700px', position: 'absolute', left: '-10000px'}});
    document.body.appendChild(container);
    let changes = 0;
    const e: TestEditor = {flow: new FlowEditor(container, {onSelectionChanged: () => changes++}), container};
    try {
      const [a, b] = await threeNodes(e);
      const {from, to} = marqueeOver(e, [a, b]);
      const before = changes;
      await dragMarquee(e, from, to);
      expect(await until(() => changes > before), true, 'marquee reports a selection change');
      expect(isSelected(e, a) && isSelected(e, b), true, 'marquee selected the nodes');

      const beforeKeys = changes;
      window.dispatchEvent(new KeyboardEvent('keydown', {key: 'a', ctrlKey: true, bubbles: true, cancelable: true}));
      expect(await until(() => changes > beforeKeys), true, 'Ctrl+A reports a selection change');

      const beforeToggle = changes;
      await clickNode(e, a, {ctrlKey: true});
      expect(await until(() => changes > beforeToggle), true, 'ctrl+click toggle-off reports');
    } finally {
      destroyEditor(e);
    }
  });

  test('Ctrl+A selects every node, Ctrl+Shift+A deselects all', async () => {
    const e = makeEditor();
    try {
      const [a, b, c] = await threeNodes(e);
      window.dispatchEvent(new KeyboardEvent('keydown', {key: 'a', ctrlKey: true, bubbles: true, cancelable: true}));
      expect(await until(() => [a, b, c].every((n) => isSelected(e, n))), true, 'all selected');

      window.dispatchEvent(new KeyboardEvent('keydown',
        {key: 'A', ctrlKey: true, shiftKey: true, bubbles: true, cancelable: true}));
      expect(await until(() => [a, b, c].every((n) => !isSelected(e, n))), true, 'all deselected');
    } finally {
      destroyEditor(e);
    }
  });
});
