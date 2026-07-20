/** Tests for the Outputs strip — the thin column OUTSIDE the canvas viewport
 *  hosting every output node as a screen-space chip (zoom-independent), with
 *  wires ending analytically at the canvas' right edge. */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions} from '../rete/node-factory';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';

/** The chip element for a node id. */
function chipEl(container: HTMLElement, id: string): HTMLElement | null {
  return container.querySelector(`.ff-output-row[data-node-id="${id}"]`);
}

/** Whether the chip exists and lies fully inside the strip column. */
function inStrip(container: HTMLElement, id: string): boolean {
  const el = chipEl(container, id);
  const strip = container.querySelector('.ff-output-strip');
  if (!el || !strip) return false;
  const c = el.getBoundingClientRect();
  const s = strip.getBoundingClientRect();
  return c.left >= s.left - 1 && c.right <= s.right + 1 && c.top >= s.top - 1 && c.bottom <= s.bottom + 1;
}

category('Flow: output strip', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('strip column mounts OUTSIDE the canvas viewport with a vertical label', async () => {
    const e = makeEditor();
    try {
      const strip = e.container.querySelector('.ff-output-strip') as HTMLElement;
      expect(strip != null, true, 'strip mounted');
      expect(strip.querySelector('.ff-output-strip-header')?.textContent, 'Outputs');
      const canvasR = e.flow.canvasEl.getBoundingClientRect();
      const stripR = strip.getBoundingClientRect();
      expect(stripR.left >= canvasR.right - 1, true,
        `strip sits right of the canvas (canvas right ${canvasR.right}, strip left ${stripR.left})`);
      expect(stripR.width <= 50, true, `thin column (${stripR.width}px)`);
      await until(() => strip.dataset.empty === 'true');
      expect(strip.dataset.empty, 'true', 'empty until an output exists');
    } finally {
      destroyEditor(e);
    }
  });

  test('an output node renders as a chip INSIDE the strip; its canvas view is hidden', async () => {
    const e = makeEditor();
    try {
      const out = await addNode(e.flow, 'Outputs/Table Output', 100, 100);
      out.properties['paramName'] = 'finalResult';
      await e.flow.updateNode(out.id);
      expect(await until(() => inStrip(e.container, out.id)), true, 'chip inside the strip column');
      const el = chipEl(e.container, out.id)!;
      const r = el.getBoundingClientRect();
      expect(Math.round(r.width), 40, 'fixed chip width');
      expect(Math.round(r.height), 24, 'fixed chip height');
      expect(el.dataset.bound, 'false', 'unbound chip marked');
      expect(el.dataset.nodeTypeName, 'Outputs/Table Output', 'guide-addressable type name');
      expect(el.querySelector('[data-testid="ff-output-row-letter"]')?.textContent?.length, 1, 'one-letter type');
      expect(el.querySelector('[data-testid="ff-socket-input-table"]') != null, true, 'socket dot rendered');
      expect((el.title ?? '').includes('finalResult'), true, 'tooltip names the output');
      // The canvas node view is an empty hidden placeholder — no card anywhere.
      const canvasCard = e.flow.canvasEl.querySelector(`.ff-node[data-node-id="${out.id}"]`);
      expect(canvasCard == null, true, 'no canvas card for an output node');
    } finally {
      destroyEditor(e);
    }
  });

  test('chips are screen-space: pan and zoom change neither size nor position', async () => {
    const e = makeEditor();
    try {
      const a = await addNode(e.flow, 'Outputs/Table Output', 0, 0);
      const b = await addNode(e.flow, 'Outputs/Value Output', 0, 0);
      expect(await until(() => inStrip(e.container, a.id) && inStrip(e.container, b.id)), true, 'both in strip');
      const topOf = (id: string): number => chipEl(e.container, id)!.getBoundingClientRect().top;
      expect(topOf(a.id) < topOf(b.id), true, 'chips stack in insertion order');
      // The group is vertically centered: its center of mass sits at the strip
      // middle (each added chip re-balances the rest).
      const strip = e.container.querySelector('.ff-output-strip')!.getBoundingClientRect();
      const rA = chipEl(e.container, a.id)!.getBoundingClientRect();
      const rB = chipEl(e.container, b.id)!.getBoundingClientRect();
      const groupCenter = (rA.top + rB.bottom) / 2;
      expect(Math.abs(groupCenter - (strip.top + strip.height / 2)) < 2, true,
        `chip group centered in the strip (group ${groupCenter}, strip mid ${strip.top + strip.height / 2})`);

      const before = chipEl(e.container, a.id)!.getBoundingClientRect();
      await e.flow.area.area.translate(-380, 240); // pan
      await e.flow.area.area.zoom(2.2); // zoom
      await new Promise((r) => setTimeout(r, 50));
      const after = chipEl(e.container, a.id)!.getBoundingClientRect();
      expect(Math.abs(after.left - before.left) < 1 && Math.abs(after.top - before.top) < 1, true,
        'position unchanged by pan/zoom');
      expect(Math.abs(after.width - before.width) < 1 && Math.abs(after.height - before.height) < 1, true,
        'size unchanged by pan/zoom');
    } finally {
      destroyEditor(e);
    }
  });

  test('wire endpoints track the canvas edge at the chip row across pan and zoom', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input', 0, 0);
      const out = await addNode(e.flow, 'Outputs/Table Output', 0, 0);
      await e.flow.addConnectionByKeys(input.id, 'table', out.id, 'table');

      const flow = e.flow as unknown as {
        listenChipSocket(id: string, cb: (p: {x: number; y: number}) => void): () => void;
      };
      let last: {x: number; y: number} | null = null;
      const unsub = flow.listenChipSocket(out.id, (p) => {last = p;});
      const screenX = (): number => {
        const t = e.flow.area.area.transform;
        return t.x + last!.x * t.k;
      };
      const screenY = (): number => {
        const t = e.flow.area.area.transform;
        return t.y + last!.y * t.k;
      };
      expect(last != null, true, 'endpoint emitted on subscribe');
      expect(Math.abs(screenX() - e.flow.canvasEl.clientWidth) < 1, true, 'endpoint at the canvas right edge');
      // Level with the chip: the sole chip is centered → endpoint at mid-height.
      expect(Math.abs(screenY() - e.flow.canvasEl.clientHeight / 2) < 1, true,
        'endpoint level with the centered chip');

      await e.flow.area.area.translate(-250, 130);
      await e.flow.area.area.zoom(1.8);
      const settled = await until(() => Math.abs(screenX() - e.flow.canvasEl.clientWidth) < 1);
      expect(settled, true, 'endpoint stays at the edge after pan/zoom');
      unsub();
    } finally {
      destroyEditor(e);
    }
  });

  test('clicking a chip selects its node (full pointer gesture)', async () => {
    const e = makeEditor();
    try {
      const out = await addNode(e.flow, 'Outputs/Value Output', 0, 0);
      expect(await until(() => chipEl(e.container, out.id) != null), true);
      // The REAL gesture — pointerdown, pointerup, click on the same element.
      // A strip that rebuilds its chips on pointerup kills the browser's click
      // dispatch (the pressed element vanishes mid-gesture), so the chip must
      // survive the whole sequence.
      const chip = chipEl(e.container, out.id)!;
      chip.dispatchEvent(new PointerEvent('pointerdown', {bubbles: true, button: 0}));
      chip.dispatchEvent(new PointerEvent('pointerup', {bubbles: true, button: 0}));
      await new Promise((r) => setTimeout(r, 0)); // the pointerup microtasks run
      expect(chipEl(e.container, out.id) === chip, true, 'chip element survives the gesture');
      chip.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      const selected = await until(() => e.flow.getSelectedNodeIds().includes(out.id) &&
        chipEl(e.container, out.id)?.dataset.selected === 'true');
      expect(selected, true, 'chip click selects the node and marks the chip');
    } finally {
      destroyEditor(e);
    }
  });

  test('dropping a dataframe output on the strip creates a wired Table Output chip', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input', 0, 0);
      const bind = (e.flow as unknown as {
        bindOutputToStrip(src: {nodeId: string; outputKey: string; dgType: string}): Promise<void>;
      }).bindOutputToStrip.bind(e.flow);
      await bind({nodeId: input.id, outputKey: 'table', dgType: 'dataframe'});

      const out = e.flow.getNodes().find((n) => n.dgTypeName === 'Outputs/Table Output');
      expect(out != null, true, 'Table Output created');
      expect(String(out!.properties['paramName']), 'table', 'named after the source slot');
      expect(e.flow.isInputConnected(out!.id, 'table'), true, 'wired to the source');
      expect(await until(() => chipEl(e.container, out!.id)?.dataset.bound === 'true'), true, 'chip shows bound');
      const boundTitle = await until(() => (chipEl(e.container, out!.id)?.title ?? '').includes('Table Input'));
      expect(boundTitle, true, `tooltip names the source (got: ${chipEl(e.container, out!.id)?.title})`);

      // A second drop from the same slot gets a unique name.
      await bind({nodeId: input.id, outputKey: 'table', dgType: 'dataframe'});
      const names = e.flow.getNodes().filter((n) => n.dgNodeType === 'output')
        .map((n) => String(n.properties['paramName'])).sort().join(',');
      expect(names, 'table,table2', 'param names stay unique');
    } finally {
      destroyEditor(e);
    }
  });

  test('dropping a scalar output on the strip creates a Value Output typed after the source', async () => {
    const e = makeEditor();
    try {
      const c = await addNode(e.flow, 'Constants/String', 0, 0);
      await (e.flow as unknown as {
        bindOutputToStrip(src: {nodeId: string; outputKey: string; dgType: string}): Promise<void>;
      }).bindOutputToStrip({nodeId: c.id, outputKey: 'value', dgType: 'string'});

      const out = e.flow.getNodes().find((n) => n.dgTypeName === 'Outputs/Value Output');
      expect(out != null, true, 'Value Output created');
      expect(e.flow.isInputConnected(out!.id, 'value'), true, 'wired');
      expect(String(out!.properties['outputType']), 'string', 'declared type auto-set from the source');
      const letterOk = await until(() =>
        chipEl(e.container, out!.id)?.querySelector('[data-testid="ff-output-row-letter"]')?.textContent === 's');
      expect(letterOk, true, 'chip letter reflects the detected type');
    } finally {
      destroyEditor(e);
    }
  });

  test('a data-output drag lights the strip as a drop target', async () => {
    const e = makeEditor();
    try {
      const input = await addNode(e.flow, 'Inputs/Table Input', 0, 0);
      const flow = e.flow as unknown as {
        beginConnectHints(n: string, k: string, s: string): void;
        endConnectHints(): void;
      };
      const strip = e.container.querySelector('.ff-output-strip') as HTMLElement;
      flow.beginConnectHints(input.id, 'table', 'output');
      expect(strip.classList.contains('ff-strip-droptarget'), true, 'lit during the drag');
      flow.endConnectHints();
      expect(strip.classList.contains('ff-strip-droptarget'), false, 'cleared on drop');
    } finally {
      destroyEditor(e);
    }
  });

  test('minimap and zoom-to-fit ignore output nodes', async () => {
    const e = makeEditor();
    try {
      await addNode(e.flow, 'Constants/String', 0, 0);
      await addNode(e.flow, 'Outputs/Value Output', 0, 0);
      const drawn = await until(() => e.container.querySelectorAll('.ff-minimap-node').length === 1);
      expect(drawn, true, 'one minimap rect — the chip is not part of the overview');

      await e.flow.zoomToFit();
      const k = e.flow.area.area.transform.k;
      expect(Number.isFinite(k) && k > 0, true, 'sane transform after fit');
    } finally {
      destroyEditor(e);
    }
  });

  test('destroy removes the strip and canvas wrapper', async () => {
    const e = makeEditor();
    const container = e.container;
    try {
      expect(container.querySelector('.ff-output-strip') != null, true);
    } finally {
      destroyEditor(e);
    }
    expect(container.querySelector('.ff-output-strip') == null, true, 'strip removed on destroy');
    expect(container.querySelector('.ff-canvas-wrap') == null, true, 'wrapper removed on destroy');
  });
});
