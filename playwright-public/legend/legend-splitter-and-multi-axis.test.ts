/* ---
sub_features_covered: [legend.occupancy.box-plot-points, legend.splitter.live-drag,
  legend.multi-axis.selectors, legend.multi-axis.split-labels]
--- */
// Splitter live-drag feedback, box plot point occupancy, and the multi-axis line chart
// legend. Pins: the box plot marks every drawn point into the occupancy grid, not just its
// box/whisker rects, so the auto corner never lands on outside-value points; dragging the
// legend splitter resizes the legend before mouseup even on viewers with a debounced
// render (trellis); multi-axis legend items keep their live y-selector widgets across
// frames (per-frame refresh must not kill widgets that live inside legend items); and
// multi-axis split items are prefixed with their y column (the visible label is
// LegendItemInfo.title, not the raw category).

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

async function legendState(page: any, viewerType: string) {
  return await page.evaluate((vt: string) => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === vt);
    if (!x) return {noViewer: true};
    const root = x.root.querySelector('[name="legend"]') as HTMLElement;
    if (!root) return {noLegend: true};
    const vr = x.root.getBoundingClientRect();
    const r = root.getBoundingClientRect();
    return {mode: root.getAttribute('data-legend-mode'), slot: root.getAttribute('data-legend-slot'),
      items: Number(root.getAttribute('data-legend-items')),
      box: {x: Math.round(r.x - vr.x), y: Math.round(r.y - vr.y),
        w: Math.round(r.width), h: Math.round(r.height)},
      viewer: {w: Math.round(vr.width), h: Math.round(vr.height)}};
  }, viewerType);
}

async function closeViewer(page: any, viewerType: string) {
  await page.evaluate((vt: string) => {
    const tv = (window as any).grok.shell.tv;
    tv.viewers.find((q: any) => q.type === vt)?.close();
  }, viewerType);
  await page.waitForTimeout(500);
}

test('Box plot point occupancy, trellis splitter live drag (SPGI)', async ({page}) => {
  test.setTimeout(900_000);
  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await v.openTable(page);

  await softStep('Box plot with points everywhere keeps the legend docked, not on the points', async () => {
    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Box plot');
    });
    await page.waitForTimeout(2000);
    // set after adding — the box plot resets a color passed at construction to its category
    await v.setViewerProps(page, 'Box plot', [{set: {markerColorColumnName: 'Series'}, wait: 2000}]);
    await v.resizeViewer(page, 'Box plot', 950, 780);
    await page.waitForTimeout(1500);
    const s = await legendState(page, 'Box plot');
    expect(s.mode, `legend went to a corner over the points: ${JSON.stringify(s)}`).toBe('docked');
    await closeViewer(page, 'Box plot');
  });

  await softStep('Dragging the trellis legend splitter resizes the legend before mouseup', async () => {
    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Trellis plot');
    });
    await page.waitForTimeout(4000);
    await v.resizeViewer(page, 'Trellis plot', 1000, 700);
    await page.waitForTimeout(2000);
    const before = await legendState(page, 'Trellis plot');
    expect(before.mode).toBe('docked');

    const res = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Trellis plot');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      const splitter = x.root.querySelector('.splitbar-vertical') as HTMLElement;
      if (!splitter) return {noSplitter: true};
      const r = splitter.getBoundingClientRect();
      const x0 = r.x + r.width / 2, y0 = r.y + r.height / 2;
      const ev = (type: string, cx: number) => new MouseEvent(type,
        {bubbles: true, cancelable: true, clientX: cx, clientY: y0});
      const startWidth = root.offsetWidth;
      splitter.dispatchEvent(ev('mousedown', x0));
      const widths: number[] = [];
      for (let dx = 15; dx <= 90; dx += 15) {
        document.dispatchEvent(ev('mousemove', x0 - dx));
        await new Promise((r2) => setTimeout(r2, 30));
        widths.push(root.offsetWidth);
      }
      const duringDrag = Math.max(...widths);
      document.dispatchEvent(ev('mouseup', x0 - 90));
      return {startWidth, widths, duringDrag};
    });
    expect(res.noSplitter, 'trellis legend splitter not found').toBeFalsy();
    expect(res.duringDrag - res.startWidth,
      `no live feedback during drag: ${JSON.stringify(res)}`).toBeGreaterThanOrEqual(40);
    await closeViewer(page, 'Trellis plot');
  });

  expect(errors, `page errors:\n${errors.join('\n')}`).toEqual([]);
  await v.cleanupShell(page);
  v.finishSpec('Box plot occupancy / trellis splitter failures');
});

test('Multi-axis legend selectors and split labels (SPGI)', async ({page}) => {
  test.setTimeout(900_000);
  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await v.openTable(page);

  const yCols = ['Average Mass', 'TPSA', 'NIBR logP', 'Chemical Space X'];
  await page.evaluate((cols: string[]) => {
    (window as any).grok.shell.tv.addViewer('Line chart',
      {xColumnName: 'Competition assay Date', yColumnNames: cols, multiAxis: true});
  }, yCols);
  await page.waitForTimeout(3000);
  await v.resizeViewer(page, 'Line chart', 950, 620);
  await page.waitForTimeout(2000);

  const selectorItems = async () => await page.evaluate(() => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === 'Line chart');
    const items = Array.from(x.root.querySelectorAll(
      '[name="legend"] .d4-legend-axis-selector')) as HTMLElement[];
    return items.map((el) => {
      const selector = el.querySelector('.d4-line-chart-y-selector') as HTMLElement;
      return {children: el.children.length,
        selectorInDom: !!selector,
        selectorVisible: !!selector && selector.offsetWidth > 0,
        text: (el.textContent ?? '').trim().slice(0, 40)};
    });
  });

  await softStep('Every multi-axis legend item shows its y-axis selector, not just the close button', async () => {
    const items = await selectorItems();
    expect(items.length, 'no selector items in the legend').toBe(yCols.length);
    for (const it of items)
      expect(it.selectorVisible, `selector missing from item: ${JSON.stringify(it)}`).toBe(true);
  });

  await softStep('Selectors survive a resize (no per-frame kill, no flicker)', async () => {
    await v.resizeViewer(page, 'Line chart', 780, 540);
    await page.waitForTimeout(1500);
    const items = await selectorItems();
    expect(items.length).toBe(yCols.length);
    for (const it of items)
      expect(it.selectorVisible, `selector gone after resize: ${JSON.stringify(it)}`).toBe(true);
  });

  await softStep('Multi-axis + split items are prefixed with their y column', async () => {
    await v.setViewerProps(page, 'Line chart', [{set: {splitColumnNames: ['Series']}, wait: 2500}]);
    const labels: string[] = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Line chart');
      return Array.from(x.root.querySelectorAll('[name="legend"] .d4-legend-item'))
        .map((el: any) => (el.textContent ?? '').trim());
    });
    expect(labels.length, 'no split legend items').toBeGreaterThan(0);
    const prefixed = labels.filter((l) => yCols.some((c) => l.startsWith(`${c} / `)));
    expect(prefixed.length, `labels lost their y-column prefix:\n${labels.join('\n')}`)
      .toBe(labels.length);
    expect(new Set(labels).size, `duplicate labels — every block shows the same values:\n${labels.join('\n')}`)
      .toBe(labels.length);
  });

  expect(errors, `page errors:\n${errors.join('\n')}`).toEqual([]);
  await v.cleanupShell(page);
  v.finishSpec('Multi-axis legend failures');
});
