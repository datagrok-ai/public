/* ---
sub_features_covered: [legend.corner.hysteresis-band, legend.corner.priority,
  legend.occupancy.row-cap, legend.placement.zoom-stable, legend.items.icon-no-shrink]
--- */
// Corner placement stability and its supporting invariants. Pins: repeated shrinks never
// flap a legend corner↔docked (the enter/stay checks form a two-sided band around
// cornerStickyPaddingPx — entering needs clearance, staying tolerates creep); the corner
// ladder is pure priority order, so a free rightTop always wins; zooming/panning a
// filter-by-zoom scatter never moves the legend even as categories vanish (item counts are
// not placement inputs); a narrowed docked legend ellipsizes labels but never crushes the
// marker icons (flex: 0 0 auto); above 20K rows corner auto-placement and occupancy
// marking are skipped entirely while an explicit corner still works; and a floating pie
// resized slowly at a fractional devicePixelRatio never flaps docked↔corner (the pie's
// marks and corner offsets describe the same viewer-frame future-circle geometry, so the
// decision cannot destroy its own justification).

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';

async function openDemog(page: any) {
  await page.evaluate(async (path: string) => {
    const grok = (window as any).grok;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
  }, demogPath);
  await page.waitForFunction(() =>
    (window as any).grok?.shell?.tv?.dataFrame?.rowCount > 0, null, {timeout: 30000});
}

async function recordPlacements(page: any) {
  await page.evaluate(() => {
    (window as any).__legendMoves = [];
    if ((window as any).__legendSub) return;
    (window as any).__legendSub = (window as any).grok.events
      .onEvent('d4-legend-placement-changed')
      .subscribe((a: any) => (window as any).__legendMoves.push(
        `${a.args.width}x${a.args.height} ${a.args.fromMode}/${a.args.fromSlot}->${a.args.toMode}/${a.args.toSlot}`));
  });
}

const moves = (page: any): Promise<string[]> =>
  page.evaluate(() => (window as any).__legendMoves);

async function legendState(page: any, viewerType: string) {
  return await page.evaluate((vt: string) => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === vt);
    const root = x?.root.querySelector('[name="legend"]') as HTMLElement;
    if (!root) return {noLegend: true};
    return {mode: root.getAttribute('data-legend-mode'), slot: root.getAttribute('data-legend-slot'),
      items: Number(root.getAttribute('data-legend-items'))};
  }, viewerType);
}

function expectNoOscillation(seen: string[], label: string) {
  expect(seen.length, `${label}: too many placement moves:\n${seen.join('\n')}`)
    .toBeLessThanOrEqual(1);
  const slots = seen.map((s) => s.split('->')[1]);
  const repeats = slots.filter((s, i) => i >= 2 && slots[i - 2] === s);
  expect(repeats, `${label}: returned to a slot it had just left:\n${seen.join('\n')}`).toEqual([]);
}

test('Pie and bar corner stability on repeated shrinks (demog)', async ({page}) => {
  test.setTimeout(900_000);
  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await openDemog(page);

  await softStep('A fresh pie legend settles directly, with no docked flash before the corner', async () => {
    await recordPlacements(page);
    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Pie chart', {categoryColumnName: 'RACE'});
    });
    await v.waitForLegendIdle(page, 'Pie chart');
    const seen = await moves(page);
    const jumped = seen.some((m, i) => m.includes('->docked') &&
      seen.slice(i + 1).some((n) => n.includes('->corner')));
    expect(jumped, `legend showed docked, then jumped to a corner:\n${seen.join('\n')}`).toBe(false);
  });

  await softStep('Pie: shrinking height repeatedly never flaps corner<->docked', async () => {
    await page.evaluate(() => { (window as any).__legendMoves = []; });
    for (const h of [700, 600, 500, 420, 340, 280])
      await v.resizeViewer(page, 'Pie chart', 1100, h);
    expectNoOscillation(await moves(page), 'pie height shrink');
  });

  await softStep('Pie: a square viewer prefers the free rightTop corner', async () => {
    await v.resizeViewer(page, 'Pie chart', 620, 620);
    const s = await legendState(page, 'Pie chart');
    expect(`${s.mode}/${s.slot}`).toBe('corner/rightTop');
    await page.evaluate(() => {
      (window as any).grok.shell.tv.viewers.find((q: any) => q.type === 'Pie chart')?.close();
    });
  });

  await softStep('Bar: full-height bars — no oscillation, corner never lands on the long top bars', async () => {
    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Bar chart',
        {splitColumnName: 'DEMOG', stackColumnName: 'RACE'});
    });
    await v.waitForLegendIdle(page, 'Bar chart');
    await recordPlacements(page);
    await page.evaluate(() => { (window as any).__legendMoves = []; });
    for (const [w, h] of [[1100, 620], [900, 620], [700, 620], [500, 620],
      [1100, 620], [1100, 480], [1100, 340], [1100, 280]] as any) {
      await v.resizeViewer(page, 'Bar chart', w, h);
      const s = await legendState(page, 'Bar chart');
      expect(s.slot, `legend on the longest bars at ${w}x${h}`).not.toBe('rightTop');
    }
    expectNoOscillation(await moves(page), 'bar shrink');
  });

  expect(errors, `page errors:\n${errors.join('\n')}`).toEqual([]);
  await v.cleanupShell(page);
  v.finishSpec('Corner stability failures');
});

test('Zoom stability, icon clipping, large-dataset gate (demog)', async ({page}) => {
  test.setTimeout(900_000);
  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await openDemog(page);

  await page.evaluate(() => {
    (window as any).grok.shell.tv.addViewer('Scatter plot',
      {xColumnName: 'HEIGHT', yColumnName: 'WEIGHT', colorColumnName: 'RACE',
       zoomAndFilter: 'filter by zoom'});
  });
  await v.waitForLegendIdle(page, 'Scatter plot');
  await v.resizeViewer(page, 'Scatter plot', 950, 620);

  await softStep('Zoom and pan under filter-by-zoom never move the legend, even as categories vanish', async () => {
    await recordPlacements(page);
    const res = await page.evaluate(async () => {
      const x = (window as any).grok.shell.tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      const itemCounts: number[] = [];
      const step = async () => {
        await new Promise((r) => setTimeout(r, 800));
        itemCounts.push(Number(root.getAttribute('data-legend-items')));
      };
      for (let i = 0; i < 5; i++) {
        const vp = x.viewport;
        x.viewport = {x: vp.x + vp.width * 0.15, y: vp.y + vp.height * 0.15,
          width: vp.width * 0.7, height: vp.height * 0.7};
        await step();
      }
      for (let i = 0; i < 4; i++) {
        const vp = x.viewport;
        x.viewport = {x: vp.x - vp.width * 0.6, y: vp.y, width: vp.width, height: vp.height};
        await step();
      }
      return {itemCounts};
    });
    // without a real category change the rest proves nothing
    expect(new Set(res.itemCounts).size, `item counts never changed: ${res.itemCounts}`)
      .toBeGreaterThan(1);
    expect(await moves(page), 'zoom/pan re-resolved the placement').toEqual([]);
  });

  await softStep('Narrowing a docked dual legend truncates labels, never the marker icons', async () => {
    await page.evaluate(() => {
      const x = (window as any).grok.shell.tv.viewers.find((q: any) => q.type === 'Scatter plot');
      x.props.markersColumnName = 'SEX';
      x.props.legendPosition = 'Right';
    });
    await v.waitForLegendIdle(page, 'Scatter plot');
    const res = await page.evaluate(async () => {
      const x = (window as any).grok.shell.tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      const splitter = x.root.querySelector('.splitbar-vertical') as HTMLElement;
      const r = splitter.getBoundingClientRect();
      const x0 = r.x + r.width / 2, y0 = r.y + r.height / 2;
      const ev = (type: string, cx: number) => new MouseEvent(type,
        {bubbles: true, cancelable: true, clientX: cx, clientY: y0});
      const startWidth = root.offsetWidth;
      splitter.dispatchEvent(ev('mousedown', x0));
      document.dispatchEvent(ev('mousemove', x0 + (startWidth - 42)));
      await new Promise((r2) => setTimeout(r2, 200));
      document.dispatchEvent(ev('mouseup', x0 + (startWidth - 42)));
      await new Promise((r2) => setTimeout(r2, 800));
      const icons = Array.from(root.querySelectorAll('.d4-legend-item'))
        .map((el: any) => el.querySelector('i, .grok-icon') as HTMLElement)
        .filter((ic) => !!ic)
        .map((ic) => ({w: ic.offsetWidth, shrink: getComputedStyle(ic).flexShrink}));
      return {endWidth: root.offsetWidth, icons};
    });
    expect(res.endWidth, 'the drag never narrowed the legend').toBeLessThan(60);
    expect(res.icons.length, 'no marker icons found').toBeGreaterThan(0);
    for (const ic of res.icons) {
      expect(ic.shrink).toBe('0');
      expect(ic.w, `icon crushed to ${ic.w}px`).toBeGreaterThanOrEqual(20);
    }
  });

  await softStep('Above 20K rows the auto ladder never corners; an explicit corner still works', async () => {
    await page.evaluate(async () => {
      const grok = (window as any).grok;
      const df = grok.data.demo.demog(25000);
      grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      (window as any).grok.shell.tv.addViewer('Pie chart', {categoryColumnName: 'race'});
    });
    await v.waitForLegendIdle(page, 'Pie chart');
    await v.resizeViewer(page, 'Pie chart', 1100, 500);
    const auto = await legendState(page, 'Pie chart');
    expect(auto.mode, `25K pie went to ${auto.mode}/${auto.slot}`).toBe('docked');
    await v.setViewerProps(page, 'Pie chart', [{set: {legendPosition: 'RightTop'}, wait: 1500}]);
    const explicit = await legendState(page, 'Pie chart');
    expect(`${explicit.mode}/${explicit.slot}`).toBe('corner/rightTop');
  });

  expect(errors, `page errors:\n${errors.join('\n')}`).toEqual([]);
  await v.cleanupShell(page);
  v.finishSpec('Zoom stability / icon clipping / row-cap failures');
});

test('Floating pie: slow height growth at fractional DPR never flaps (SPGI)', async ({page}) => {
  test.setTimeout(900_000);
  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await v.openTable(page);

  await softStep('Dragging the floating dialog taller moves the legend at most once', async () => {
    await page.evaluate(() => {
      Object.defineProperty(window, 'devicePixelRatio', {get: () => 0.8999999, configurable: true});
      (window as any).grok.shell.tv.addViewer('Pie chart', {categoryColumnName: 'Primary Series Name'});
    });
    await v.waitForLegendIdle(page, 'Pie chart');
    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Pie chart');
      tv.dockManager.findNode(x.root).container.float();
    });
    await v.waitForLegendIdle(page, 'Pie chart');

    // drag the SE corner handle so the floating dialog becomes exactly 395x513
    await page.evaluate(async () => {
      const dlg = document.querySelector('.dialog-floating') as HTMLElement;
      const se = dlg.querySelector('.resize-handle-se') as HTMLElement;
      const r0 = dlg.getBoundingClientRect();
      const hr = se.getBoundingClientRect();
      const ev = (type: string, x: number, y: number) => new MouseEvent(type,
        {bubbles: true, cancelable: true, view: window, clientX: x, clientY: y, button: 0});
      const sx = hr.x + hr.width / 2;
      const sy = hr.y + hr.height / 2;
      const dx = 395 - r0.width;
      const dy = 513 - r0.height;
      se.dispatchEvent(ev('mousedown', sx, sy));
      for (let i = 1; i <= 10; i++) {
        document.dispatchEvent(ev('mousemove', sx + dx * i / 10, sy + dy * i / 10));
        await new Promise((r) => setTimeout(r, 30));
      }
      document.dispatchEvent(ev('mouseup', sx + dx, sy + dy));
    });
    await v.waitForLegendIdle(page, 'Pie chart');
    await recordPlacements(page);
    await page.evaluate(() => { (window as any).__legendMoves = []; });

    // slowly drag the south handle 190px down, 2px per move
    await page.evaluate(async () => {
      const dlg = document.querySelector('.dialog-floating') as HTMLElement;
      const s = dlg.querySelector('.resize-handle-s') as HTMLElement;
      const hr = s.getBoundingClientRect();
      const ev = (type: string, cx: number, cy: number) => new MouseEvent(type,
        {bubbles: true, cancelable: true, view: window, clientX: cx, clientY: cy, button: 0});
      const sx = hr.x + hr.width / 2;
      const sy = hr.y + hr.height / 2;
      s.dispatchEvent(ev('mousedown', sx, sy));
      for (let d = 2; d <= 190; d += 2) {
        document.dispatchEvent(ev('mousemove', sx, sy + d));
        await new Promise((r) => setTimeout(r, 130));
      }
      document.dispatchEvent(ev('mouseup', sx, sy + 190));
    });
    await v.waitForLegendIdle(page, 'Pie chart');
    expectNoOscillation(await moves(page), 'floating pie height growth');
  });

  expect(errors, `page errors:\n${errors.join('\n')}`).toEqual([]);
  await v.cleanupShell(page);
  v.finishSpec('Floating pie stability failures');
});
