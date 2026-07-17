/** Annotation carry: dragging an annotation moves whatever sits inside it —
 *  nodes whose center is within the frame, annotations fully contained in it,
 *  and the waypoints of connections between carried nodes. The cargo is
 *  captured at drag START (stateless): a node dragged out of the frame simply
 *  isn't inside at the next grab. Strip-pinned output rows are never carried. */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes} from '../rete/node-factory';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';

/** Simulate a body drag on an annotation element (down → move → up). The drag
 *  handler only reads client deltas, so absolute coords are arbitrary. */
function drag(el: HTMLElement, dx: number, dy: number): void {
  const at = (x: number, y: number): PointerEventInit =>
    ({bubbles: true, cancelable: true, button: 0, clientX: 500 + x, clientY: 500 + y});
  el.dispatchEvent(new PointerEvent('pointerdown', at(0, 0)));
  el.dispatchEvent(new PointerEvent('pointermove', at(dx, dy)));
  el.dispatchEvent(new PointerEvent('pointerup', at(dx, dy)));
}

category('Flow: annotations', () => {
  before(async () => {
    registerBuiltinNodes();
  });

  test('dragging an annotation carries the nodes inside it', async () => {
    const e = makeEditor();
    try {
      const inside = await addNode(e.flow, 'Inputs/String Input', 20, 20);
      const outside = await addNode(e.flow, 'Inputs/String Input', 600, 20);
      await until(() => !!e.container.querySelector(`.ff-node[data-node-id="${outside.id}"]`));
      const ann = e.flow.addAnnotation({pos: {x: -20, y: -20}, size: {w: 320, h: 220}});

      drag(ann.element, 100, 50);
      expect(ann.pos.x, 80, 'annotation moved');
      expect(await until(() => inside.pos.x === 120 && inside.pos.y === 70), true,
        'the node inside moved by the same delta');
      expect(outside.pos.x, 600, 'the node outside did not move');
    } finally {
      destroyEditor(e);
    }
  });

  test('a node dragged out of the annotation detaches from the next drag', async () => {
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, 'Inputs/String Input', 20, 20);
      await until(() => !!e.container.querySelector(`.ff-node[data-node-id="${node.id}"]`));
      const ann = e.flow.addAnnotation({pos: {x: -20, y: -20}, size: {w: 320, h: 220}});

      drag(ann.element, 40, 0);
      expect(await until(() => node.pos.x === 60), true, 'carried while inside');

      await e.flow.translate(node.id, 2000, 2000); // user drags the node out
      drag(ann.element, 40, 0);
      await new Promise((r) => setTimeout(r, 100));
      expect(node.pos.x, 2000, 'a node outside the frame is no longer carried');
    } finally {
      destroyEditor(e);
    }
  });

  test('nested annotations move along; overlapping ones do not', async () => {
    const e = makeEditor();
    try {
      const outer = e.flow.addAnnotation({pos: {x: 0, y: 0}, size: {w: 400, h: 300}});
      const nested = e.flow.addAnnotation({pos: {x: 40, y: 40}, size: {w: 100, h: 80}});
      const straddling = e.flow.addAnnotation({pos: {x: 350, y: 40}, size: {w: 200, h: 80}});

      drag(outer.element, 60, 30);
      expect(nested.pos.x, 100, 'fully-contained annotation carried');
      expect(nested.pos.y, 70);
      expect(straddling.pos.x, 350, 'partially-overlapping annotation stays');
    } finally {
      destroyEditor(e);
    }
  });

  test('waypoints of connections between carried nodes move too', async () => {
    const e = makeEditor();
    try {
      const a = await addNode(e.flow, 'Inputs/Table Input', 20, 20);
      const b = await addNode(e.flow, 'Utilities/Select Column', 20, 160);
      await until(() => !!e.container.querySelector(`.ff-node[data-node-id="${b.id}"]`));
      await e.flow.addConnectionByKeys(a.id, 'table', b.id, 'table');
      const conn = e.flow.getConnections()[0];
      e.flow.addWaypoint(conn, {x: 150, y: 100});
      const ann = e.flow.addAnnotation({pos: {x: -20, y: -20}, size: {w: 400, h: 320}});

      drag(ann.element, 25, 35);
      expect(await until(() => conn.waypoints?.[0].x === 175), true, 'waypoint x carried');
      expect(conn.waypoints![0].y, 135, 'waypoint y carried');
    } finally {
      destroyEditor(e);
    }
  });

  test('strip-pinned output rows are never carried', async () => {
    const e = makeEditor();
    try {
      const out = await addNode(e.flow, 'Outputs/Table Output', 20, 20);
      const ann = e.flow.addAnnotation({pos: {x: -50, y: -50}, size: {w: 500, h: 400}});
      drag(ann.element, 80, 0);
      await new Promise((r) => setTimeout(r, 100));
      expect(out.pos.x, 20, 'output row position untouched');
    } finally {
      destroyEditor(e);
    }
  });
});
