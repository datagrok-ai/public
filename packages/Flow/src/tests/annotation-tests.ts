/** Annotation carry: dragging an annotation moves the nodes, contained annotations,
 *  and waypoints inside it; the cargo is captured at drag start. Plus styling:
 *  the color palette and the serialized title font size. */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';
import * as ui from 'datagrok-api/ui';

import {registerBuiltinNodes} from '../rete/node-factory';
import {FlowEditor} from '../rete/flow-editor';
import {ANNOTATION_COLORS, ANNOTATION_TITLE_SIZES, ANN_DEFAULT_FONT_SIZE} from '../rete/annotation';
import {serializeFlow, deserializeFlow} from '../serialization/flow-serializer';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';

const SETTINGS = {scriptName: 'AnnStyle', scriptDescription: '', tags: []};

function hexToRgb(hex: string): string {
  const n = parseInt(hex.slice(1), 16);
  return `rgb(${(n >> 16) & 255}, ${(n >> 8) & 255}, ${n & 255})`;
}

/** The drag handler only reads client deltas, so absolute coords are arbitrary. */
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

      await e.flow.translate(node.id, 2000, 2000);
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

  test('Delete removes the last-clicked annotation; clicking elsewhere disarms', async () => {
    const e = makeEditor();
    try {
      const ann = e.flow.addAnnotation({pos: {x: 10, y: 10}, size: {w: 200, h: 120}, text: 'doomed'});
      const down: PointerEventInit = {bubbles: true, cancelable: true, button: 0, clientX: 0, clientY: 0};
      ann.element.dispatchEvent(new PointerEvent('pointerdown', down));
      ann.element.dispatchEvent(new PointerEvent('pointerup', down));
      expect(ann.element.classList.contains('ff-annotation-active'), true, 'clicked → armed and marked');

      e.container.querySelector('.ff-canvas')!.dispatchEvent(new PointerEvent('pointerdown', down));
      window.dispatchEvent(new KeyboardEvent('keydown', {key: 'Delete', bubbles: true, cancelable: true}));
      await new Promise((r) => setTimeout(r, 50));
      expect(e.flow.getAnnotations().length, 1, 'a disarmed annotation survives Delete');
      expect(ann.element.classList.contains('ff-annotation-active'), false, 'the mark is gone');

      ann.element.dispatchEvent(new PointerEvent('pointerdown', down));
      ann.element.dispatchEvent(new PointerEvent('pointerup', down));
      window.dispatchEvent(new KeyboardEvent('keydown', {key: 'Delete', bubbles: true, cancelable: true}));
      expect(await until(() => e.flow.getAnnotations().length === 0), true,
        'Delete removed the clicked annotation');
    } finally {
      destroyEditor(e);
    }
  });

  test('with nodes selected, Delete removes the nodes, not the armed annotation', async () => {
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, 'Inputs/String Input', 600, 20);
      const ann = e.flow.addAnnotation({pos: {x: 10, y: 10}, size: {w: 100, h: 80}});
      const down: PointerEventInit = {bubbles: true, cancelable: true, button: 0, clientX: 0, clientY: 0};
      ann.element.dispatchEvent(new PointerEvent('pointerdown', down));
      ann.element.dispatchEvent(new PointerEvent('pointerup', down));
      await e.flow.selectNode(node.id);
      window.dispatchEvent(new KeyboardEvent('keydown', {key: 'Delete', bubbles: true, cancelable: true}));
      expect(await until(() => e.flow.getNodes().length === 0), true, 'the selected node went');
      expect(e.flow.getAnnotations().length, 1, 'the annotation stayed — selection wins');
    } finally {
      destroyEditor(e);
    }
  });

  test('the color palette is broad, unique, and border-matched', async () => {
    expect(ANNOTATION_COLORS.length >= 12, true, `only ${ANNOTATION_COLORS.length} colors`);
    expect(new Set(ANNOTATION_COLORS.map((c) => c.bg.toLowerCase())).size, ANNOTATION_COLORS.length,
      'backgrounds are unique');
    expect(new Set(ANNOTATION_COLORS.map((c) => c.name)).size, ANNOTATION_COLORS.length, 'names are unique');
    for (const c of ANNOTATION_COLORS) {
      expect(/^#[0-9A-Fa-f]{6}$/.test(c.bg) && /^#[0-9A-Fa-f]{6}$/.test(c.border), true,
        `${c.name} has hex bg+border`);
    }
    const e = makeEditor();
    try {
      // A post-expansion color must resolve its own border, not the blue fallback.
      const orange = ANNOTATION_COLORS.find((c) => c.name === 'Orange')!;
      const ann = e.flow.addAnnotation({color: orange.bg});
      expect(ann.element.style.background, hexToRgb(orange.bg), 'background applied');
      expect(ann.element.style.borderColor, hexToRgb(orange.border), 'matching border resolved');
    } finally {
      destroyEditor(e);
    }
  });

  test('title font size: default, custom, malformed fallback, tidy toDoc', async () => {
    const e = makeEditor();
    try {
      const plain = e.flow.addAnnotation({});
      expect(plain.fontSize, ANN_DEFAULT_FONT_SIZE, 'default size');
      expect(plain.titleEl.style.fontSize, `${ANN_DEFAULT_FONT_SIZE}px`, 'default applied to the title');
      expect('fontSize' in plain.toDoc(), false, 'default size stays out of the doc');

      const big = e.flow.addAnnotation({fontSize: 20});
      expect(big.titleEl.style.fontSize, '20px', 'custom size applied');
      expect(big.toDoc().fontSize, 20, 'custom size serialized');

      const bad = e.flow.addAnnotation({fontSize: Number.NaN});
      expect(bad.fontSize, ANN_DEFAULT_FONT_SIZE, 'malformed size falls back');

      expect(ANNOTATION_TITLE_SIZES.some((s) => s.size === ANN_DEFAULT_FONT_SIZE), true,
        'the menu offers the default size');
    } finally {
      destroyEditor(e);
    }
  });

  test('color and title size round-trip through the .flow doc', async () => {
    const e = makeEditor();
    const e2 = makeEditor();
    try {
      const green = ANNOTATION_COLORS.find((c) => c.name === 'Green')!;
      e.flow.addAnnotation({text: 'styled', color: green.bg, fontSize: 26});
      e.flow.addAnnotation({text: 'plain'});
      const doc = serializeFlow(e.flow, SETTINGS);
      const styled = doc.annotations!.find((a) => a.text === 'styled')!;
      expect(styled.fontSize, 26, 'doc carries the size');
      expect(styled.color, green.bg, 'doc carries the color');
      expect('fontSize' in doc.annotations!.find((a) => a.text === 'plain')!, false,
        'the untouched annotation saves without the key');

      await deserializeFlow(doc, e2.flow);
      const loaded = e2.flow.getAnnotations().find((a) => a.text === 'styled')!;
      expect(loaded.fontSize, 26, 'size survives the load');
      expect(loaded.titleEl.style.fontSize, '26px', 'and is applied');
      expect(loaded.element.style.borderColor, hexToRgb(green.border), 'color border survives too');
    } finally {
      destroyEditor(e2);
      destroyEditor(e);
    }
  });

  test('the editor setters apply the style and mark the graph changed', async () => {
    const container = ui.div([], {style: {width: '1000px', height: '700px', position: 'absolute', left: '-10000px'}});
    document.body.appendChild(container);
    let changes = 0;
    const flow = new FlowEditor(container, {onGraphChanged: () => changes++});
    try {
      const ann = flow.addAnnotation({});
      const red = ANNOTATION_COLORS.find((c) => c.name === 'Red')!;
      const base = changes;
      flow.setAnnotationColor(ann.id, red.bg);
      expect(ann.element.style.background, hexToRgb(red.bg), 'recolored');
      expect(changes, base + 1, 'recolor is an unsaved change');
      flow.setAnnotationColor(ann.id, red.bg);
      expect(changes, base + 1, 'same color again is not');

      flow.setAnnotationFontSize(ann.id, 20);
      expect(ann.titleEl.style.fontSize, '20px', 'resized');
      expect(changes, base + 2, 'resize is an unsaved change');
      flow.setAnnotationFontSize(ann.id, 20);
      expect(changes, base + 2, 'same size again is not');
    } finally {
      try {
        flow.destroy();
      } finally {
        container.remove();
      }
    }
  });

  test('right-click offers Color and Title size groups', async () => {
    const e = makeEditor();
    try {
      const ann = e.flow.addAnnotation({pos: {x: 10, y: 10}, size: {w: 200, h: 120}});
      ann.element.dispatchEvent(new MouseEvent('contextmenu',
        {bubbles: true, cancelable: true, clientX: 20, clientY: 20}));
      const labels = (): string[] => Array.from(document.querySelectorAll<HTMLElement>('.d4-menu-item-label'))
        .map((el) => el.textContent?.trim() ?? '');
      expect(await until(() => labels().includes('Color') && labels().includes('Title size') &&
        labels().includes('Delete')), true, `menu groups missing; saw: ${labels().join(', ')}`);
    } finally {
      for (const el of Array.from(document.querySelectorAll('.d4-menu-popup, .d4-menu-dropdown')))
        el.remove();
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
