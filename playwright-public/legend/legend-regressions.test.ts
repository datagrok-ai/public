/* ---
sub_features_covered: [legend.placement.stability, legend.sizing.estimate,
  legend.sections.split, legend.item.click-filter, legend.item.color-picker]
--- */
// Placement and sizing regressions of the frame-bound legend:
// placement recomputed while panning/zooming or resizing the legend, a dual legend sized
// to its 40% ceiling regardless of its labels, the second section squeezed to nothing,
// a crash clicking a scatter-plot marker item, and a line-chart split legend that neither
// filters nor recolours.
//
// `data-legend-passes` counts placement resolutions and `data-legend-model` model rebuilds:
// both must stay put while only the viewport moves.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

interface Counters {passes: number; model: number; mode: string; slot: string}

async function counters(page: any, viewerType: string): Promise<Counters> {
  return await page.evaluate((vt: string) => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === vt);
    const root = x.root.querySelector('[name="legend"]') as HTMLElement;
    return {passes: Number(root.getAttribute('data-legend-passes') ?? '-1'),
      model: Number(root.getAttribute('data-legend-model') ?? '-1'),
      mode: root.getAttribute('data-legend-mode') ?? '',
      slot: root.getAttribute('data-legend-slot') ?? ''};
  }, viewerType);
}

test('Legend placement — only the viewer size moves it', async ({page}) => {
  test.setTimeout(900_000);

  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {column: 'Stereo Category', viewers: [
    {type: 'Scatter plot', column: 'Stereo Category'}]});
  await page.evaluate(() => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
    x.props.markersColumnName = 'Series';
  });
  await v.resizeViewer(page, 'Scatter plot', 900, 600);
  await v.waitForLegendIdle(page, 'Scatter plot');

  const zoomAndPan = async (): Promise<boolean> => await page.evaluate(async () => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
    const vp0 = JSON.stringify(x.viewport);
    for (let i = 0; i < 4; i++) {
      const vp = x.viewport;
      x.viewport = {x: vp.x + vp.width * 0.1, y: vp.y + vp.height * 0.1,
        width: vp.width * 0.8, height: vp.height * 0.8};
      await new Promise((r) => setTimeout(r, 500));
    }
    for (let i = 0; i < 3; i++) {
      const vp = x.viewport;
      x.viewport = {x: vp.x + vp.width * 0.15, y: vp.y, width: vp.width, height: vp.height};
      await new Promise((r) => setTimeout(r, 500));
    }
    return vp0 !== JSON.stringify(x.viewport);
  });

  await softStep('Zooming and panning rebuild nothing when the data stands still', async () => {
    await v.setViewerProps(page, 'Scatter plot',
      [{set: {zoomAndFilter: 'no action'}, wait: 1500}]);
    await v.waitForLegendIdle(page, 'Scatter plot');
    const before = await counters(page, 'Scatter plot');
    const moved = await zoomAndPan();
    const after = await counters(page, 'Scatter plot');

    expect(moved, 'the viewport never changed — the rest proves nothing').toBe(true);
    expect(after.passes, 'zoom/pan re-resolved the placement').toBe(before.passes);
    expect(after.model, 'zoom/pan rebuilt the legend model').toBe(before.model);
  });

  await softStep('Filter-by-zoom refilters the data but still does not move the legend', async () => {
    // the default mode: zooming really does change which rows the legend sees, so its model
    // must rebuild — its *placement* must not follow
    await v.setViewerProps(page, 'Scatter plot',
      [{set: {zoomAndFilter: 'filter by zoom'}, wait: 1500}]);
    await v.waitForLegendIdle(page, 'Scatter plot');
    const rect = async () => await page.evaluate(() => {
      const x = (window as any).grok.shell.tv.viewers
        .find((q: any) => q.type === 'Scatter plot');
      const r = (x.root.querySelector('[name="legend"]') as HTMLElement)
        .getBoundingClientRect();
      const vr = x.root.getBoundingClientRect();
      return {rows: x.filter.trueCount, left: Math.round(r.left - vr.left),
        top: Math.round(r.top - vr.top), width: Math.round(r.width),
        height: Math.round(r.height)};
    });

    const before = await counters(page, 'Scatter plot');
    const r0 = await rect();
    const moved = await zoomAndPan();
    const after = await counters(page, 'Scatter plot');
    const r1 = await rect();

    expect(moved, 'the viewport never changed').toBe(true);
    expect(r1.rows, 'filter-by-zoom filtered nothing — the rest proves nothing')
      .toBeLessThan(r0.rows);
    expect([after.mode, after.slot], 'the legend moved while zooming')
      .toEqual([before.mode, before.slot]);
    expect([r1.left, r1.top, r1.width, r1.height], 'the legend box moved while zooming')
      .toEqual([r0.left, r0.top, r0.width, r0.height]);
  });

  await softStep('Dragging the legend splitter resizes it without re-placing it', async () => {
    await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'Right'}, wait: 1500}]);
    await v.waitForLegendIdle(page, 'Scatter plot');
    const before = await counters(page, 'Scatter plot');
    const width0 = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      return (x.root.querySelector('[name="legend"]') as HTMLElement).getBoundingClientRect().width;
    });

    await v.dragLegendSplitter(page, 'Scatter plot', -60, 0);
    const after = await counters(page, 'Scatter plot');
    const width1 = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      return (x.root.querySelector('[name="legend"]') as HTMLElement).getBoundingClientRect().width;
    });

    expect(Math.round(width1 - width0), 'the drag did not resize the legend')
      .toBeGreaterThan(20);
    expect(after.slot, 'a manual resize moved the legend').toBe(before.slot);
    expect(after.passes, 'a manual resize re-resolved the placement').toBe(before.passes);
  });

  expect(errors, 'the page threw while panning, zooming or resizing').toEqual([]);
  await v.cleanupShell(page);
  v.finishSpec('Legend placement stability failures');
});

test('Legend sizing — a many-category legend is measured, not maxed out', async ({page}) => {
  test.setTimeout(900_000);

  await loginToDatagrok(page);
  await v.openTable(page);
  await page.evaluate(async () => {
    const df = (window as any).grok.shell.tv.dataFrame;
    df.columns.addNewString('ShortCats').init((i: number) => 'c' + (i % 300));
    df.columns.addNewString('WideCats').init((i: number) =>
      'a-very-long-category-label-number-' + (i % 300));
    await new Promise((r) => setTimeout(r, 800));
  });

  await v.addLegendViewers(page, {column: 'ShortCats', viewers: ['Scatter plot']});
  await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'Right'}, wait: 1500}]);
  await v.resizeViewer(page, 'Scatter plot', 1000, 600);
  await v.waitForLegendIdle(page, 'Scatter plot');

  const measure = async (): Promise<{legend: number; viewer: number; label: number}> =>
    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      const labels = Array.from(root.querySelectorAll('.d4-legend-value')) as HTMLElement[];
      return {legend: root.getBoundingClientRect().width,
        viewer: x.root.getBoundingClientRect().width,
        label: Math.max(...labels.map((l) => l.scrollWidth))};
    });

  const short = await measure();
  await softStep('Short labels get a narrow legend, not the 40% ceiling', async () => {
    expect(short.legend, 'a 300-category legend of "c123" labels took its full ceiling')
      .toBeLessThan(short.viewer * 0.25);
    expect(short.legend, 'the legend is narrower than its own widest label')
      .toBeGreaterThan(short.label);
  });

  await v.setViewerProps(page, 'Scatter plot',
    [{set: {colorColumnName: 'WideCats'}, wait: 2000}]);
  await v.waitForLegendIdle(page, 'Scatter plot');
  const wide = await measure();

  await softStep('Long labels widen the legend, still under the 40% ceiling', async () => {
    expect(wide.legend, 'label length does not affect the width at all')
      .toBeGreaterThan(short.legend + 30);
    expect(wide.legend, 'the legend broke through its 40% ceiling')
      .toBeLessThanOrEqual(wide.viewer * 0.4 + 2);
  });

  await v.cleanupShell(page);
  v.finishSpec('Legend sizing failures');
});

test('Legend sections — the shorter section keeps the rows it needs', async ({page}) => {
  test.setTimeout(900_000);

  await loginToDatagrok(page);
  await v.openTable(page);
  await page.evaluate(async () => {
    const df = (window as any).grok.shell.tv.dataFrame;
    df.columns.addNewString('ManyCats').init((i: number) => 'cat-' + (i % 120));
    df.columns.addNewString('FewShapes').init((i: number) => 'shape-' + (i % 4));
    await new Promise((r) => setTimeout(r, 800));
  });

  await v.addLegendViewers(page, {column: 'ManyCats', viewers: ['Scatter plot']});
  await page.evaluate(() => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
    x.props.markersColumnName = 'FewShapes';
  });
  await v.setViewerProps(page, 'Scatter plot', [{set: {legendPosition: 'Right'}, wait: 2000}]);
  await v.resizeViewer(page, 'Scatter plot', 1000, 600);
  await v.waitForLegendIdle(page, 'Scatter plot');

  const layout = await page.evaluate(() => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
    const root = x.root.querySelector('[name="legend"]') as HTMLElement;
    const read = (id: string) => {
      const s = root.querySelector(`[data-section="${id}"]`) as HTMLElement;
      const list = s.querySelector('.d4-legend-list') as HTMLElement;
      const rows = Array.from(list.querySelectorAll('[data-item-key]')) as HTMLElement[];
      const lr = list.getBoundingClientRect();
      return {items: Number(s.getAttribute('data-section-items') ?? '0'),
        listHeight: list.clientHeight, rendered: rows.length,
        // the list positions rows absolutely, so "does it fit" is the last row's bottom
        lowestRow: rows.length === 0 ? 0
          : Math.max(...rows.map((r) => r.getBoundingClientRect().bottom - lr.top))};
    };
    return {main: read('main'), extra: read('extra'), height: root.clientHeight};
  });

  expect(layout.extra.items, 'the markers section did not build its items').toBe(4);
  expect(layout.extra.rendered, 'the markers section renders fewer rows than it has')
    .toBe(layout.extra.items);
  expect(layout.extra.lowestRow, 'the markers section has to scroll to show four rows')
    .toBeLessThanOrEqual(layout.extra.listHeight);
  expect(layout.extra.listHeight, 'the markers section took more than half the legend')
    .toBeLessThanOrEqual(layout.height / 2 + 2);
  expect(layout.main.listHeight, 'the categories section got no room')
    .toBeGreaterThan(layout.height / 3);

  await v.cleanupShell(page);
});

test('Legend markers — clicking any marker item filters instead of throwing', async ({page}) => {
  test.setTimeout(900_000);

  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {column: 'Stereo Category', viewers: ['Scatter plot']});
  await page.evaluate(() => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
    x.props.markersColumnName = 'Series';
  });
  await v.resizeViewer(page, 'Scatter plot', 1000, 600);
  await v.waitForLegendIdle(page, 'Scatter plot');

  const items = page.locator('[name="viewer-Scatter-plot"] [data-section="extra"] [data-item-key]');
  const count = await items.count();
  expect(count, 'no marker items to click').toBeGreaterThanOrEqual(2);

  const filtered: number[] = [];
  for (let i = 0; i < Math.min(count, 3); i++) {
    await items.nth(i).click();
    await page.waitForTimeout(900);
    filtered.push(await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Scatter plot');
      return x.filter.trueCount;
    }));
    // clicking the same item again clears it, so the next click starts from the same state
    await items.nth(i).click();
    await page.waitForTimeout(900);
  }

  expect(errors, 'clicking a marker legend item threw').toEqual([]);
  const total = await page.evaluate(() => (window as any).grok.shell.tv.dataFrame.rowCount);
  for (let i = 0; i < filtered.length; i++)
    expect(filtered[i], `marker item ${i} did not filter the viewer`).toBeLessThan(total);

  await v.cleanupShell(page);
});

test('Line chart — a split legend filters, recolours and keeps its colour', async ({page}) => {
  test.setTimeout(900_000);

  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {column: 'Stereo Category', viewers: ['Line chart']});
  await v.resizeViewer(page, 'Line chart', 1000, 600);
  await v.waitForLegendIdle(page, 'Line chart');

  const items = page.locator('[name="viewer-Line-chart"] [data-section="extra"] [data-item-key]');
  expect(await items.count(), 'the split legend built no items').toBeGreaterThanOrEqual(2);

  await softStep('Clicking a split item filters the chart', async () => {
    const before = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      return tv.viewers.find((q: any) => q.type === 'Line chart').filter.trueCount;
    });

    await items.nth(1).click();
    await page.waitForTimeout(1200);
    const after = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      return tv.viewers.find((q: any) => q.type === 'Line chart').filter.trueCount;
    });

    expect(after, 'clicking a split legend item did nothing').toBeLessThan(before);
    await items.nth(1).click();
    await page.waitForTimeout(1200);
  });

  await softStep('Right-clicking a split item opens the colour picker', async () => {
    await items.nth(1).click({button: 'right'});
    await page.locator('.d4-dialog .d4-color-bar').first().waitFor({timeout: 5000});
    expect(await page.locator('.d4-dialog').count(), 'right-click stacked colour pickers')
      .toBe(1);
  });

  await softStep('The recoloured item keeps its colour across renders', async () => {
    const key = await items.nth(1).getAttribute('data-item-key');
    const chosen = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog')!;
      const bars = Array.from(dlg.querySelectorAll('.d4-color-bar')) as HTMLElement[];
      const item = (window as any).grok.shell.tv.viewers
        .find((q: any) => q.type === 'Line chart').root
        .querySelectorAll('[data-section="extra"] [data-item-key]')[1] as HTMLElement;
      const current = getComputedStyle(item).color;
      const sw = bars.find((s) => s.style.backgroundColor && s.style.backgroundColor !== current);
      if (!sw) throw new Error('no swatch to pick in the colour picker');
      const o = {bubbles: true, cancelable: true, view: window, button: 0};
      sw.dispatchEvent(new MouseEvent('mousedown', o));
      sw.dispatchEvent(new MouseEvent('mouseup', o));
      sw.dispatchEvent(new MouseEvent('click', o));
      return sw.style.backgroundColor;
    });
    await page.locator('.d4-dialog [name="button-OK"]').click({timeout: 5000});
    await page.waitForTimeout(900);

    const after = await page.evaluate(async (k) => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Line chart');
      const read = () => {
        const el = x.root.querySelector(
          `[data-section="extra"] [data-item-key="${k}"]`) as HTMLElement;
        return el ? getComputedStyle(el).color : null;
      };
      const immediate = read();
      // the whole point: a repaint must not throw the colour away
      x.invalidateCanvas?.();
      await new Promise((r) => setTimeout(r, 1500));
      return {immediate, later: read()};
    }, key);

    expect(after.immediate, 'the legend item did not take the chosen colour').toBe(chosen);
    expect(after.later, 'the legend item lost its colour on the next render').toBe(chosen);
  });

  expect(errors, 'the line chart legend threw').toEqual([]);
  await v.cleanupShell(page);
  v.finishSpec('Line chart legend failures');
});
