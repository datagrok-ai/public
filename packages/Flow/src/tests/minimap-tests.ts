import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions} from '../rete/node-factory';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';

category('Flow: minimap', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('minimap draws a rect per node plus the viewport rectangle', async () => {
    const e = makeEditor();
    try {
      expect(e.container.querySelector('.ff-minimap') != null, true, 'minimap mounted');
      await addNode(e.flow, 'Constants/String', 0, 0);
      await addNode(e.flow, 'Utilities/ToString', 300, 0);
      // Redraw is rAF-coalesced — poll until both node rects are drawn.
      const drawn = await until(() => e.container.querySelectorAll('.ff-minimap-node').length === 2);
      expect(drawn, true, 'one minimap rect per node');
      expect(e.container.querySelector('.ff-minimap-viewport') != null, true, 'viewport rectangle drawn');
    } finally {
      destroyEditor(e);
    }
  });

  test('setMinimapCollapsed collapses and restores the minimap', async () => {
    const e = makeEditor();
    try {
      // Needs a node — the overview is hidden on an empty canvas (below).
      await addNode(e.flow, 'Constants/String', 0, 0);
      const mm = e.container.querySelector('.ff-minimap') as HTMLElement;
      expect(mm != null, true);
      expect(mm.dataset.collapsed, 'false', 'starts expanded');
      e.flow.setMinimapCollapsed(true);
      expect(mm.dataset.collapsed, 'true');
      e.flow.setMinimapCollapsed(false);
      expect(mm.dataset.collapsed, 'false');
      // Collapse only minimizes to the header bar — it does not hide outright.
      await until(() => mm.style.display !== 'none');
      expect(mm.style.display !== 'none', true, 'present while the canvas has nodes');
    } finally {
      destroyEditor(e);
    }
  });

  test('overview is hidden on an empty canvas, shown once a node is added', async () => {
    const e = makeEditor();
    try {
      const mm = e.container.querySelector('.ff-minimap') as HTMLElement;
      // Empty canvas: nothing to overview → hidden.
      const hidden = await until(() => mm.style.display === 'none');
      expect(hidden, true, 'hidden while empty');
      // First node lands → it reappears.
      await addNode(e.flow, 'Constants/String', 0, 0);
      const shown = await until(() => mm.style.display !== 'none');
      expect(shown, true, 'shown once a node exists');
    } finally {
      destroyEditor(e);
    }
  });

  test('clicking a spot on the minimap centers the viewport on that graph point', async () => {
    // The inverse of the draw transform: offsetX/offsetY already fold in
    // `-min * scale`, so the nav mapping must NOT add min again — that shifted
    // every pan target down/right by the graph origin (the "click the preview,
    // land far below" bug, worst for graphs placed away from (0,0)).
    const e = makeEditor();
    try {
      // Deliberately far from the origin so a double-counted min is loud.
      const a = await addNode(e.flow, 'Constants/String', 600, 500);
      await addNode(e.flow, 'Utilities/ToString', 950, 520);
      await until(() => e.container.querySelectorAll('.ff-minimap-node').length === 2);

      const svg = e.container.querySelector('.ff-minimap-svg') as SVGSVGElement & {dataset: DOMStringMap};
      const fit = {
        scale: parseFloat(svg.dataset.scale!),
        offsetX: parseFloat(svg.dataset.offsetX!),
        offsetY: parseFloat(svg.dataset.offsetY!),
      };
      // Click the minimap pixel where node A's top-left is drawn.
      const target = {x: a.pos.x, y: a.pos.y};
      const rect = svg.getBoundingClientRect();
      const body = e.container.querySelector('.ff-minimap-body') as HTMLElement;
      body.dispatchEvent(new PointerEvent('pointerdown', {
        bubbles: true, cancelable: true, button: 0,
        clientX: rect.left + target.x * fit.scale + fit.offsetX,
        clientY: rect.top + target.y * fit.scale + fit.offsetY,
      }));
      body.dispatchEvent(new PointerEvent('pointerup', {bubbles: true, cancelable: true, button: 0}));

      // The pan centers the clicked graph point in the container.
      const cont = e.container.getBoundingClientRect();
      const area = (e.flow as unknown as {area: {area: {transform: {x: number; y: number; k: number}}}}).area;
      const ok = await until(() => {
        const t = area.area.transform;
        return Math.abs(t.x - (cont.width / 2 - target.x * t.k)) < 5 &&
               Math.abs(t.y - (cont.height / 2 - target.y * t.k)) < 5;
      });
      const t = area.area.transform;
      expect(ok, true, `pan centered the clicked point (t=${t.x},${t.y} k=${t.k}, ` +
        `expected ${cont.width / 2 - target.x * t.k},${cont.height / 2 - target.y * t.k})`);
    } finally {
      destroyEditor(e);
    }
  });

  test('clicking the minimap header toggles collapse', async () => {
    const e = makeEditor();
    try {
      await addNode(e.flow, 'Constants/String', 0, 0);
      const mm = e.container.querySelector('.ff-minimap') as HTMLElement;
      const header = mm.querySelector('.ff-minimap-header') as HTMLElement;
      header.click();
      expect(mm.dataset.collapsed, 'true', 'collapsed after clicking the header');
      header.click();
      expect(mm.dataset.collapsed, 'false', 'restored on second header click');
    } finally {
      destroyEditor(e);
    }
  });
});
