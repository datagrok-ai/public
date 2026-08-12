/** Datagrok selection semantics on the canvas (mirrors d4 `selectRows` / `areaSelector`),
 *  exercised through real PointerEvent / KeyboardEvent dispatch on the live editor DOM. */
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

/** Dispatch on the area holder, not the container — the container's own AreaPlugin listeners would start a phantom pan whose click-release clears the selection. */
function canvasBg(e: TestEditor): HTMLElement {
  return e.container.firstElementChild as HTMLElement;
}

/** The tick between down and up lets rete's async `nodepicked` add settle, as a real click would. */
async function clickNode(e: TestEditor, node: FlowNode, init: PointerEventInit = {}): Promise<void> {
  const el = nodeEl(e, node);
  const r = el.getBoundingClientRect();
  const at = {clientX: r.left + r.width / 2, clientY: r.top + 8};
  pointer(el, 'pointerdown', {...at, ...init});
  await tick();
  pointer(el, 'pointerup', {...at, ...init});
  await tick();
}

/** Events target the canvas background (not `window`) so the editor's window-capture handlers run before rete's window-bubble ones, as with real events. */
async function dragMarquee(e: TestEditor, from: {x: number; y: number}, to: {x: number; y: number},
  upInit: PointerEventInit = {}): Promise<void> {
  pointer(canvasBg(e), 'pointerdown', {clientX: from.x, clientY: from.y, shiftKey: true});
  await tick();
  pointer(canvasBg(e), 'pointermove', {clientX: to.x, clientY: to.y, shiftKey: true});
  await tick();
  pointer(canvasBg(e), 'pointerup', {clientX: to.x, clientY: to.y, shiftKey: true, ...upInit});
  await tick();
}

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

      await clickNode(e, a);
      await clickNode(e, a);
      expect(fires, 1, 're-clicks do not re-fire');

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

  test('isNodeContextCurrent veto: a stale host context re-fires on re-click', async () => {
    const container = ui.div([], {style: {width: '1000px', height: '700px', position: 'absolute', left: '-10000px'}});
    document.body.appendChild(container);
    let fires = 0;
    let contextCurrent = true;
    const flow = new FlowEditor(container, {
      onNodeSelected: () => fires++,
      isNodeContextCurrent: () => contextCurrent,
    });
    const e: TestEditor = {flow, container};
    try {
      const a = await addNode(flow, 'Inputs/String Input', 0, 0);
      await until(() => !!container.querySelector(`.ff-node[data-node-id="${a.id}"]`));

      await clickNode(e, a);
      expect(await until(() => isSelected(e, a)), true, 'first click selects');
      expect(fires, 1, 'first click fires once');

      await clickNode(e, a);
      expect(fires, 1, 'context current → re-click stays deduped');

      contextCurrent = false;
      await clickNode(e, a);
      expect(fires, 2, 'context stale → re-click re-fires');

      contextCurrent = true;
      await clickNode(e, a);
      expect(fires, 2, 'context current again → deduped again');
    } finally {
      flow.destroy();
      container.remove();
    }
  });

  test('shift+drag marquee adds the covered nodes to the selection', async () => {
    const e = makeEditor();
    try {
      const [a, b, c] = await threeNodes(e);
      await e.flow.selectNode(c.id);
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
    // The marquee's release is swallowed by stopImmediatePropagation — hosts must rely on this callback.
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

  test('arrow keys nudge the selection; Ctrl+arrows pan the canvas instead', async () => {
    const e = makeEditor();
    try {
      const [a, b] = await threeNodes(e);
      const key = (k: string, ctrl = false): boolean => window.dispatchEvent(
        new KeyboardEvent('keydown', {key: k, ctrlKey: ctrl, bubbles: true, cancelable: true}));

      const ax = a.pos.x;
      key('ArrowRight');
      await new Promise((r) => setTimeout(r, 50));
      expect(a.pos.x, ax, 'an unselected node does not move');

      await e.flow.selectNode(a.id);
      await e.flow.selectNode(b.id, true);
      const [bx, by] = [b.pos.x, b.pos.y];
      key('ArrowRight');
      expect(await until(() => a.pos.x === ax + 10 && b.pos.x === bx + 10), true,
        'ArrowRight nudged every selected node by 10');
      key('ArrowUp');
      expect(await until(() => b.pos.y === by - 10), true, 'ArrowUp nudged up');

      const t0x = e.flow.area.area.transform.x;
      const [nx, ny] = [a.pos.x, a.pos.y];
      key('ArrowLeft', true);
      expect(await until(() => e.flow.area.area.transform.x !== t0x), true, 'Ctrl+arrow panned the canvas');
      await new Promise((r) => setTimeout(r, 50));
      expect(a.pos.x === nx && a.pos.y === ny, true, 'Ctrl+arrow left the selection in place');
    } finally {
      destroyEditor(e);
    }
  });

  test('a nudge is never snapped — a clicked (picked) off-grid node moves by exactly 10', async () => {
    const e = makeEditor();
    try {
      const [a] = await threeNodes(e);
      await e.flow.translate(a.id, 537, 231); // off-grid, BEFORE picking
      await clickNode(e, a);
      expect(await until(() => isSelected(e, a)), true, 'the click selected the node');
      expect(a.pos.x, 537, 'the click itself did not move it');
      window.dispatchEvent(new KeyboardEvent('keydown',
        {key: 'ArrowRight', bubbles: true, cancelable: true}));
      expect(await until(() => a.pos.x === 547), true,
        `the nudge landed on exactly 547 (got ${a.pos.x})`);
    } finally {
      destroyEditor(e);
    }
  });
});
