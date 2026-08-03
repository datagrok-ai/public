/** Node groups: create/ungroup, frame drag carries members, minimize hides
 *  members + internal wires and re-anchors boundary wires to the card's edge
 *  dots, maximize restores everything, membership survives serialization.
 *  Groups are editor-level (the graph stays flat) — see rete/node-group.ts. */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes} from '../rete/node-factory';
import {serializeFlow, deserializeFlow} from '../serialization/flow-serializer';
import {GROUP_DOT_TOP} from '../rete/node-group';
import {makeEditor, destroyEditor, addNode, until, TestEditor} from './test-utils';
import {FlowNode} from '../rete/scheme';

function key(init: KeyboardEventInit): void {
  window.dispatchEvent(new KeyboardEvent('keydown', {bubbles: true, cancelable: true, ...init}));
}

/** Body drag on a group element (down → move → up); handler reads client deltas. */
function drag(el: HTMLElement, dx: number, dy: number): void {
  const at = (x: number, y: number): PointerEventInit =>
    ({bubbles: true, cancelable: true, button: 0, clientX: 400 + x, clientY: 400 + y});
  el.dispatchEvent(new PointerEvent('pointerdown', at(0, 0)));
  el.dispatchEvent(new PointerEvent('pointermove', at(dx, dy)));
  el.dispatchEvent(new PointerEvent('pointerup', at(dx, dy)));
}

/** Whether a node's canvas view is hidden inside a minimized group. */
function isHidden(e: TestEditor, nodeId: string): boolean {
  const el = e.container.querySelector(`.ff-node[data-node-id="${nodeId}"]`);
  return !!el && el.closest('.ff-group-hidden') !== null;
}

/** Start/end coords of a rendered connection's path. */
function pathEndpoints(e: TestEditor, connId: string): {start: {x: number; y: number}; end: {x: number; y: number}} | null {
  const wrap = e.container.querySelector(`[data-connection-id="${connId}"]`);
  const d = wrap?.querySelector('path')?.getAttribute('d') ?? '';
  const nums = d.match(/-?\d+(\.\d+)?/g)?.map(Number) ?? [];
  if (nums.length < 4) return null;
  return {
    start: {x: nums[0], y: nums[1]},
    end: {x: nums[nums.length - 2], y: nums[nums.length - 1]},
  };
}

/** A(outside) → C(member); D(member) → B(member) internal; B → V(output row).
 *  So the group [D, B, C] has one boundary input (C.table) and one boundary
 *  output (B.column). */
async function fixture(e: TestEditor): Promise<{a: FlowNode; d: FlowNode; b: FlowNode; c: FlowNode; v: FlowNode}> {
  const a = await addNode(e.flow, 'Inputs/Table Input', 0, 200);
  const d = await addNode(e.flow, 'Inputs/Table Input', 300, 0);
  const b = await addNode(e.flow, 'Utilities/Select Column', 550, 0);
  const c = await addNode(e.flow, 'Utilities/Select Column', 550, 150);
  const v = await addNode(e.flow, 'Outputs/Value Output', 900, 0);
  await e.flow.addConnectionByKeys(a.id, 'table', c.id, 'table');
  await e.flow.addConnectionByKeys(d.id, 'table', b.id, 'table');
  await e.flow.addConnectionByKeys(b.id, 'column', v.id, 'value');
  await until(() => !!e.container.querySelector(`.ff-node[data-node-id="${c.id}"]`));
  return {a, d, b, c, v};
}

category('Flow: groups', () => {
  before(async () => {
    registerBuiltinNodes();
  });

  test('createGroup wraps the members in a frame; strip outputs are excluded', async () => {
    const e = makeEditor();
    try {
      const {d, b, c, v} = await fixture(e);
      const g = e.flow.createGroup([d.id, b.id, c.id, v.id])!;
      expect(g.memberIds.size, 3, 'the output-strip node was filtered out');
      expect(g.memberIds.has(v.id), false);
      expect(e.flow.groupOf(b.id)?.id, g.id);
      expect(g.element.dataset.minimized, 'false');
      expect(await until(() => Math.abs(g.pos.x - (300 - 14)) < 1), true,
        'frame refits to the member bbox (left edge = min x - pad)');
      expect(e.flow.createGroup([b.id]), null, 'a grouped node cannot join another group');
    } finally {
      destroyEditor(e);
    }
  });

  test('dragging the frame carries the members', async () => {
    const e = makeEditor();
    try {
      const {d, b, c} = await fixture(e);
      const g = e.flow.createGroup([d.id, b.id, c.id])!;
      await until(() => Math.abs(g.pos.x - 286) < 1);
      drag(g.element, 100, 60);
      expect(await until(() => b.pos.x === 650 && b.pos.y === 60), true, 'member b followed');
      expect(await until(() => d.pos.x === 400), true, 'member d followed');
      expect(await until(() => Math.abs(g.pos.x - 386) < 1), true, 'frame refitted around the moved members');
    } finally {
      destroyEditor(e);
    }
  });

  test('minimize hides members + internal wires; boundary wires re-anchor to card dots', async () => {
    const e = makeEditor();
    try {
      const {a, d, b, c} = await fixture(e);
      const g = e.flow.createGroup([d.id, b.id, c.id])!;
      await until(() => Math.abs(g.pos.x - 286) < 1);
      await e.flow.minimizeGroup(g.id);

      expect(g.element.dataset.minimized, 'true');
      expect(isHidden(e, d.id), true, 'member views hidden');
      expect(isHidden(e, b.id), true);
      expect(isHidden(e, a.id), false, 'outside node stays visible');

      const conns = e.flow.getConnections();
      const internal = conns.find((x) => x.source === d.id && x.target === b.id)!;
      const boundaryIn = conns.find((x) => x.source === a.id && x.target === c.id)!;
      expect(await until(() => {
        const el = e.container.querySelector(`[data-connection-id="${internal.id}"]`);
        return !!el && el.classList.contains('ff-group-hidden');
      }), true, 'internal wire hidden');

      expect(g.element.querySelectorAll('.ff-group-dots-in .ff-group-dot').length, 1, 'one input dot');
      expect(g.element.querySelectorAll('.ff-group-dots-out .ff-group-dot').length, 1, 'one output dot');

      // The A → C wire's tail must land on the card's left-edge dot.
      expect(await until(() => {
        const ep = pathEndpoints(e, boundaryIn.id);
        return !!ep && Math.abs(ep.end.x - g.pos.x) < 6 &&
          Math.abs(ep.end.y - (g.pos.y + GROUP_DOT_TOP)) < 6;
      }), true, 'boundary wire re-anchored to the card edge');
    } finally {
      destroyEditor(e);
    }
  });

  test('dragging the card moves the hidden members and their wire anchors', async () => {
    const e = makeEditor();
    try {
      const {a, d, b, c} = await fixture(e);
      const g = e.flow.createGroup([d.id, b.id, c.id])!;
      await until(() => Math.abs(g.pos.x - 286) < 1);
      await e.flow.minimizeGroup(g.id);
      const bBefore = {...b.pos};
      const cardBefore = {...g.pos};

      drag(g.element, 70, 30);
      expect(g.pos.x, cardBefore.x + 70, 'card moved');
      expect(await until(() => b.pos.x === bBefore.x + 70 && b.pos.y === bBefore.y + 30), true,
        'hidden members traveled with the card');

      const boundaryIn = e.flow.getConnections().find((x) => x.source === a.id && x.target === c.id)!;
      expect(await until(() => {
        const ep = pathEndpoints(e, boundaryIn.id);
        return !!ep && Math.abs(ep.end.x - g.pos.x) < 6;
      }), true, 'wire anchor followed the card');
    } finally {
      destroyEditor(e);
    }
  });

  test('maximize restores members, wires, and the frame around them', async () => {
    const e = makeEditor();
    try {
      const {d, b, c} = await fixture(e);
      const g = e.flow.createGroup([d.id, b.id, c.id])!;
      await until(() => Math.abs(g.pos.x - 286) < 1);
      await e.flow.minimizeGroup(g.id);
      drag(g.element, 200, 0); // travel while minimized
      await until(() => d.pos.x === 500);
      await e.flow.maximizeGroup(g.id);

      expect(g.element.dataset.minimized, 'false');
      expect(isHidden(e, d.id), false, 'members visible again');
      const internal = e.flow.getConnections().find((x) => x.source === d.id && x.target === b.id)!;
      const el = e.container.querySelector(`[data-connection-id="${internal.id}"]`);
      expect(el?.classList.contains('ff-group-hidden'), false, 'internal wire visible again');
      expect(await until(() => Math.abs(g.pos.x - 486) < 1), true,
        'frame reappears around the traveled members');
    } finally {
      destroyEditor(e);
    }
  });

  test('ungroup dissolves the frame; members and wires stay', async () => {
    const e = makeEditor();
    try {
      const {d, b, c} = await fixture(e);
      const nodeCount = e.flow.getNodeCount();
      const connCount = e.flow.getConnectionCount();
      const g = e.flow.createGroup([d.id, b.id, c.id])!;
      await e.flow.minimizeGroup(g.id);
      e.flow.ungroup(g.id);

      expect(e.flow.getGroups().length, 0);
      expect(e.flow.groupOf(b.id) === undefined, true, 'membership gone');
      expect(e.flow.getNodeCount(), nodeCount, 'no node lost');
      expect(e.flow.getConnectionCount(), connCount, 'no wire lost');
      expect(g.element.isConnected, false, 'frame element removed');
      expect(isHidden(e, d.id), false, 'members un-hidden on ungroup of a minimized group');
    } finally {
      destroyEditor(e);
    }
  });

  test('deleting members shrinks the group; the last one dissolves it', async () => {
    const e = makeEditor();
    try {
      const {d, b} = await fixture(e);
      const g = e.flow.createGroup([d.id, b.id])!;
      await e.flow.removeNode(d.id);
      expect(g.memberIds.size, 1);
      expect(e.flow.getGroups().length, 1);
      await e.flow.removeNode(b.id);
      expect(e.flow.getGroups().length, 0, 'empty group dissolved');
    } finally {
      destroyEditor(e);
    }
  });

  test('Ctrl+G groups the selection, Ctrl+Shift+G ungroups it', async () => {
    const e = makeEditor();
    try {
      const {d, b} = await fixture(e);
      await e.flow.selectNode(d.id, true);
      await e.flow.selectNode(b.id, true);
      key({key: 'g', ctrlKey: true});
      expect(await until(() => e.flow.getGroups().length === 1), true, 'Ctrl+G grouped');

      await e.flow.selectNode(d.id);
      key({key: 'g', ctrlKey: true, shiftKey: true});
      expect(await until(() => e.flow.getGroups().length === 0), true, 'Ctrl+Shift+G ungrouped');
    } finally {
      destroyEditor(e);
    }
  });

  test('Ctrl+A skips members hidden inside a minimized group', async () => {
    const e = makeEditor();
    try {
      const {a, d, b, c, v} = await fixture(e);
      const g = e.flow.createGroup([d.id, b.id])!;
      await e.flow.minimizeGroup(g.id);
      key({key: 'a', ctrlKey: true});
      expect(await until(() => {
        const sel = new Set(e.flow.getSelectedNodeIds());
        return sel.has(a.id) && sel.has(c.id) && sel.has(v.id) && !sel.has(d.id) && !sel.has(b.id);
      }), true, 'visible nodes selected, hidden members skipped');
    } finally {
      destroyEditor(e);
    }
  });

  test('the card aggregates member run statuses', async () => {
    const e = makeEditor();
    try {
      const {d, b} = await fixture(e);
      const g = e.flow.createGroup([d.id, b.id])!;
      await e.flow.minimizeGroup(g.id);

      (d as FlowNode & {dgStatus?: string}).dgStatus = 'running';
      await e.flow.updateNode(d.id);
      expect(g.element.dataset.status, 'running');

      (d as FlowNode & {dgStatus?: string}).dgStatus = 'completed';
      (b as FlowNode & {dgStatus?: string}).dgStatus = 'completed';
      await e.flow.updateNode(d.id);
      expect(g.element.dataset.status, 'completed', 'all done → completed');

      (b as FlowNode & {dgStatus?: string}).dgStatus = 'errored';
      await e.flow.updateNode(b.id);
      expect(g.element.dataset.status, 'errored', 'any error wins over done');
    } finally {
      destroyEditor(e);
    }
  });

  test('groups survive a serialize/deserialize round-trip', async () => {
    const e = makeEditor();
    const e2 = makeEditor();
    try {
      const {d, b} = await fixture(e);
      const g = e.flow.createGroup([d.id, b.id], {title: 'Prep', description: 'load & pick'})!;
      await e.flow.minimizeGroup(g.id);
      const doc = serializeFlow(e.flow, {scriptName: 'x', scriptDescription: '', tags: []});
      expect(doc.groups!.length, 1);

      await deserializeFlow(doc, e2.flow);
      const g2 = e2.flow.getGroups()[0];
      expect(g2.title, 'Prep');
      expect(g2.description, 'load & pick');
      expect(g2.minimized, true);
      expect(g2.memberIds.size, 2, 'member ids remapped through the loader idMap');
      for (const id of g2.memberIds)
        expect(!!e2.flow.getNodeById(id), true, 'remapped member exists');
      const hiddenMember = Array.from(g2.memberIds)[0];
      expect(await until(() => isHidden(e2, hiddenMember)), true,
        'loaded minimized group hides its members');
    } finally {
      destroyEditor(e);
      destroyEditor(e2);
    }
  });
});
