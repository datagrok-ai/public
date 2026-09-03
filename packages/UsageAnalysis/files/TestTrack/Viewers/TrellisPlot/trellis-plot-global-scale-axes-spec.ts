/* ---
realizes: [trellisplot.cp.global-scale-inner-axes, trellisplot.int.global-scale-range-slider-sync]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const isBenignError = (text: string) =>
  /Failed to load resource/.test(text) || /404 \(\)/.test(text) || /favicon/.test(text) ||
  /Unable to find element in cloned iframe/.test(text);

// a cell whose canvas has not been painted yet hashes to 0; a baseline read there compares a blank
// canvas with the painted one and books the paint as a zoom
async function paintedHash(page: Page, idx: number): Promise<number | null> {
  return v.pollValue(async () => (await cellHashes(page, [idx]))[0], (h) => h !== null && h !== 0, 2000, 50);
}

async function cellHashes(page: Page, idxs: number[]): Promise<(number | null)[]> {
  return page.evaluate((idxs) => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    function hash(cellIdx: number): number | null {
      const cell = root.querySelectorAll('.d4-trellis-plot-cell')[cellIdx];
      const cv = cell?.querySelector('canvas') as HTMLCanvasElement | null;
      if (!cv) return null;
      try {
        const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
        let h = 0;
        for (let i = 0; i < img.length; i += 4)
          h = (h * 31 + ((img[i] << 16) | (img[i + 1] << 8) | img[i + 2])) % 2147483647;
        return h;
      } catch { return null; }
    }
    return idxs.map(hash);
  }, idxs);
}

async function trellisMenuLabels(page: Page): Promise<string[]> {
  const labels = await page.evaluate(async () => {
    const w = window as any;
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    const grid = (root.querySelector('.d4-trellis-plot-charts-grid') as HTMLElement) ?? root;
    const gr = grid.getBoundingClientRect();
    const items = () => Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
      .map((e) => (e as HTMLElement).innerText.trim());
    grid.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true,
      clientX: gr.left + 4, clientY: gr.top + 4}));
    // the menu is built in stages (inner viewer items first, trellis items after), so settle on
    // the label list holding still rather than on its first item
    await w.__poll(() => items().length, (n: number) => n > 0, 800, 40);
    await w.__settledFor(() => items().join('|'), 250, 800, 25);
    const out = items();
    document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
    document.body.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
    await w.__poll(() => document.querySelectorAll('.d4-menu-popup').length, (n: number) => n === 0, 400, 40);
    return out;
  });
  return labels;
}

async function rangeSliderCount(page: Page): Promise<number> {
  return page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    return root.querySelectorAll('[type="range-slider"]').length;
  });
}

async function cellCenter(page: Page, idx: number): Promise<{x: number; y: number} | null> {
  return page.evaluate((idx) => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    const cell = root.querySelectorAll('.d4-trellis-plot-cell')[idx] as HTMLElement | undefined;
    if (!cell) return null;
    const r = cell.getBoundingClientRect();
    return {x: r.left + r.width / 2, y: r.top + r.height / 2};
  }, idx);
}

async function setInnerType(page: Page, viewerType: string): Promise<void> {
  await page.evaluate((vt) => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    tp.props.viewerType = vt;
  }, viewerType);
  await v.waitForViewerRendered(page, 'Trellis plot', 900);
}

async function assertInnerTypeSwitched(page: Page, viewerType: string,
  beforeSwitch: number | null): Promise<void> {
  const applied = await page.evaluate(() => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    return tp.props.viewerType as string;
  });
  const [afterSwitch] = await cellHashes(page, [0]);
  console.log(`[inner type -> ${viewerType}] readBack=${applied} hash ${beforeSwitch} -> ${afterSwitch}`);
  expect(applied, `inner viewer type did not switch to ${viewerType}`).toBe(viewerType);
  expect(beforeSwitch).not.toBeNull();
  expect(afterSwitch).not.toBeNull();
  expect(afterSwitch,
    `cells were not repainted by the switch to ${viewerType} — the arm below would probe the previous inner type`)
    .not.toBe(beforeSwitch);
}

// The trellis keeps ONE active inner slider per axis (features/inner_viewer_axes.dart) plus hidden
// per-column copies whose handles are never shown; the active one is revealed by cell.onMouseEnter
// and stays revealed on the axis strip, and takes mousedown on its handle elements only.
async function dragInnerRangeSlider(page: Page, axis: 'x' | 'y', rootIndex = 0): Promise<boolean> {
  const cell = await page.evaluate((rootIdx) => {
    const root = document.querySelectorAll('[name="viewer-Trellis-plot"]')[rootIdx];
    const c = Array.from(root?.querySelectorAll('.d4-trellis-plot-cell') ?? []).find((x) => x.querySelector('canvas'));
    if (!c) return null;
    const b = c.getBoundingClientRect();
    return {x: b.x + b.width / 2, y: b.y + b.height / 2};
  }, rootIndex);
  if (!cell) return false;
  // the reveal is cell.onMouseEnter, so the pointer has to come from outside the cell
  await page.mouse.move(2, 2);
  await page.mouse.move(cell.x, cell.y);
  const geo = await v.pollValue(() => page.evaluate(({rootIdx, ax}) => {
    const root = document.querySelectorAll('[name="viewer-Trellis-plot"]')[rootIdx];
    const svgs = Array.from(root?.querySelectorAll(`.d4-range-selector > svg[type="range-slider"][name="${ax}-slider"]`) ?? []) as SVGElement[];
    const shownHandles = (svg: SVGElement) => ['min-handle', 'max-handle']
      .map((n) => svg.querySelector(`[name="${n}"]`) as SVGElement | null)
      .filter((h) => !!h && getComputedStyle(h).display !== 'none' && h.getBoundingClientRect().width > 0)
      .map((h) => { const b = h!.getBoundingClientRect(); return {x: b.x + b.width / 2, y: b.y + b.height / 2}; });
    for (const svg of svgs) {
      const wrap = svg.closest('.d4-range-selector') as HTMLElement | null;
      if (!wrap || getComputedStyle(wrap).visibility === 'hidden' || getComputedStyle(svg).visibility === 'hidden') continue;
      const hs = shownHandles(svg);
      if (hs.length !== 2) continue;
      const b = svg.getBoundingClientRect();
      const end = hs.sort((p, q) => ax === 'x' ? q.x - p.x : q.y - p.y)[0];
      return {end, svg: {x: b.x, y: b.y, w: b.width, h: b.height}};
    }
    return null;
  }, {rootIdx: rootIndex, ax: axis}), (g) => g !== null, 1500, 30);
  if (!geo) return false;
  await page.mouse.move(geo.end.x, geo.end.y, {steps: 4});
  await page.mouse.down();
  if (axis === 'x') await page.mouse.move(geo.svg.x + geo.svg.w * 0.45, geo.end.y, {steps: 12});
  else await page.mouse.move(geo.end.x, geo.svg.y + geo.svg.h * 0.45, {steps: 12});
  await page.mouse.up();
  return true;
}

const innerTypeTabs = ['Scatter plot', 'Bar chart', 'Box plot', 'Histogram', 'Line chart', 'Pie chart'];

async function openInnerViewerTab(page: Page, tabName?: string): Promise<void> {

  await v.openViewerGear(page, 'Trellis plot');

  await page.evaluate(({name, names}) => {
    const headers = Array.from(document.querySelectorAll('.d4-tab-header')) as HTMLElement[];
    const tab = (name ? headers.find((h) => h.innerText.trim() === name) : undefined) ??
      headers.find((h) => names.includes(h.innerText.trim()));
    tab?.click();
  }, {name: tabName ?? '', names: innerTypeTabs});
  const row = page.locator('.property-grid tr[name="prop-allow-zoom"]');
  await row.first().waitFor({state: 'attached', timeout: 800}).catch(() => {});
  if (await row.count() > 0 && !(await row.isVisible())) {
    await page.locator('.property-grid tr[name="prop-category-misc"]').first().click();
    await row.first().waitFor({state: 'visible', timeout: 500}).catch(() => {});
  }
}

async function allowZoomState(page: Page): Promise<boolean | null> {
  return page.evaluate(() => {
    const row = document.querySelector('.property-grid tr[name="prop-allow-zoom"]');
    const cb = row?.querySelector('input[type="checkbox"]') as HTMLInputElement | null;
    return cb ? cb.checked : null;
  });
}

async function setAllowZoom(page: Page, desired: boolean, tabName = 'Scatter plot'): Promise<boolean | null> {
  await openInnerViewerTab(page, tabName);
  const row = page.locator('.property-grid tr[name="prop-allow-zoom"]');
  if (await row.count() === 0) return null;
  await row.waitFor({state: 'visible', timeout: 10000});
  const box = row.locator('input[type="checkbox"]').first();
  if (await box.isChecked() !== desired) {
    await box.click();

    await v.waitForViewerRendered(page, 'Trellis plot', 900);
  }
  return box.isChecked();
}

// the post-wheel window is the assertion for the no-zoom arms (nothing may repaint), so it
// stays as a capped poll for the hash leaving `before`; the zoom arm exits on the first change
async function wheelOver(page: Page, pt: {x: number; y: number}, before: number | null, steps = 5): Promise<number | null> {
  await page.mouse.move(pt.x, pt.y);
  await page.waitForTimeout(200);
  for (let i = 0; i < steps; i++) {
    await page.mouse.wheel(0, 120);
    await page.waitForTimeout(60);
  }
  return v.pollValue(async () => (await cellHashes(page, [0]))[0], (h) => h !== before, 600, 50);
}

test('Trellis plot: global scale, inner axes, range slider reset', async ({page}) => {
  test.setTimeout(240_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => {
    if (m.type() === 'error' && !isBenignError(m.text()) && !isLocalBootNoise(m.text()))
      consoleErrors.push(m.text());
  });

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  const setup = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    return {rowCount: df.rowCount, sex: df.col('SEX').categories.length, race: df.col('RACE').categories.length};
  });
  expect(setup).toEqual({rowCount: 5850, sex: 2, race: 4});

  await v.addViewerByIcon(page, 'trellis-plot', 'Trellis-plot', 15000);
  await page.evaluate(() => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    tp.props.xColumnNames = ['SEX'];
    tp.props.yColumnNames = ['RACE'];
    tp.props.viewerType = 'Scatter plot';
    tp.props.showRangeSliders = true;
    tp.props.showXAxes = 'Always';
    tp.props.showYAxes = 'Always';
  });
  await v.waitForViewerRendered(page, 'Trellis plot', 900);
  const cellLocator = page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell');
  await expect(cellLocator).toHaveCount(8);

  await softStep('Scenario 1 Step 5', async () => {

    const probes = [0, 2];
    const errBefore = consoleErrors.length;
    const pageErrBefore = pageErrors.length;
    const before = await cellHashes(page, probes);
    await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.globalScale = true;
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    const after = await cellHashes(page, probes);

    expect(before.every((h) => h !== null)).toBe(true);
    expect(after.every((h) => h !== null)).toBe(true);
    console.log(`[Scenario 1 Step 5] before=${JSON.stringify(before)} after=${JSON.stringify(after)}`);
    expect(after[0]).not.toBe(before[0]);
    expect(after[1]).not.toBe(before[1]);
    expect(after[0]).not.toBe(after[1]);
    expect(consoleErrors.slice(errBefore)).toEqual([]);
    expect(pageErrors.slice(pageErrBefore)).toEqual([]);
  });

  await softStep('Scenario 1 Step 7', async () => {
    const labels = await trellisMenuLabels(page);
    expect(labels).toContain('Reset Inner Range Sliders');
  });

  await softStep('Scenario 1 Step 9', async () => {

    await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.showXAxes = 'Never';
      tp.props.showYAxes = 'Never';
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    const labels = await trellisMenuLabels(page);

    expect(labels).toContain('Properties...');
    expect(labels).not.toContain('Reset Inner Range Sliders');
  });

  await softStep('Scenario 1 Step 10', async () => {
    await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.showXAxes = 'Always';
      tp.props.showYAxes = 'Always';
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    const labels = await trellisMenuLabels(page);
    expect(labels).toContain('Reset Inner Range Sliders');
  });

  await softStep('Scenario 1 Step 11', async () => {
    const countAxesShown = await rangeSliderCount(page);
    await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.showXAxes = 'Never';
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    const countXHidden = await rangeSliderCount(page);
    console.log(`[Scenario 1 Step 11] count_axes_shown=${countAxesShown} count_x_hidden=${countXHidden}`);
    expect(countAxesShown).toBeGreaterThan(countXHidden);

    await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.showXAxes = 'Always';
      tp.props.showYAxes = 'Always';
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    await expect(cellLocator).toHaveCount(8);
  });

  await softStep('Scenario 2 Step 4', async () => {

    await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.globalScale = true;
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    const probes = [0, 2];

    const errBefore = consoleErrors.length;
    const pageErrBefore = pageErrors.length;
    const baseline = await cellHashes(page, probes);
    expect(baseline.every((h) => h !== null)).toBe(true);
    expect(baseline[0]).not.toBe(baseline[1]);

    const dragged = await dragInnerRangeSlider(page, 'x');
    expect(dragged).toBe(true);
    const after = await v.pollValue(() => cellHashes(page, probes),
      (h) => h[0] !== baseline[0] && h[1] !== baseline[1], 2000, 50);
    console.log(`[Scenario 2 Step 4] baseline=${JSON.stringify(baseline)} after=${JSON.stringify(after)}`);
    expect(after.every((h) => h !== null)).toBe(true);

    expect(after[0]).not.toBe(baseline[0]);
    expect(after[1]).not.toBe(baseline[1]);
    expect(after[0]).not.toBe(after[1]);
    expect(consoleErrors.slice(errBefore)).toEqual([]);
    expect(pageErrors.slice(pageErrBefore)).toEqual([]);

    await page.evaluate(({narrowed, base}) => {
      (window as any).__narrowed = narrowed;
      (window as any).__baseline = base;
    }, {narrowed: after, base: baseline});
  });

  await softStep('Scenario 2 Step 5', async () => {

    const probes = [0, 2];
    const narrowed = await page.evaluate(() => (window as any).__narrowed as (number | null)[]);
    const baseline = await page.evaluate(() => (window as any).__baseline as (number | null)[]);
    const errBefore = consoleErrors.length;
    const pageErrBefore = pageErrors.length;
    await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const grid = (root.querySelector('.d4-trellis-plot-charts-grid') as HTMLElement) ?? root;
      const gr = grid.getBoundingClientRect();
      grid.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true,
        clientX: gr.left + 4, clientY: gr.top + 4}));
      const find = () => Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
        .find((e) => (e as HTMLElement).innerText.trim() === 'Reset Inner Range Sliders');
      const target = await (window as any).__poll(find, (e: Element | undefined) => !!e, 800, 40);
      (target?.closest('.d4-menu-item') as HTMLElement | null)?.click();
    });
    const after = await v.pollValue(() => cellHashes(page, probes),
      (h) => h[0] === baseline[0] && h[1] === baseline[1], 2000, 50);
    console.log(`[Scenario 2 Step 5] baseline=${JSON.stringify(baseline)} narrowed=${JSON.stringify(narrowed)} ` +
      `after=${JSON.stringify(after)}`);
    expect(after.every((h) => h !== null)).toBe(true);
    expect(after[0]).not.toBe(narrowed[0]);
    expect(after[1]).not.toBe(narrowed[1]);
    expect(after[0], 'cell 0 did not return exactly to its pre-drag baseline').toBe(baseline[0]);
    expect(after[1], 'cell 2 did not return exactly to its pre-drag baseline').toBe(baseline[1]);
    expect(consoleErrors.slice(errBefore)).toEqual([]);
    expect(pageErrors.slice(pageErrBefore)).toEqual([]);
  });

  await page.evaluate(() => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    tp.props.globalScale = false;
    tp.props.showRangeSliders = false;
    tp.props.xColumnNames = ['SEX'];
    tp.props.yColumnNames = ['RACE'];
    tp.props.viewerType = 'Scatter plot';
  });
  await v.waitForViewerRendered(page, 'Trellis plot', 900);
  await expect(cellLocator).toHaveCount(8);

  try {

    await softStep('Scenario 3 Step 2', async () => {
      await setInnerType(page, 'Scatter plot');

      await openInnerViewerTab(page, 'Scatter plot');
      const defaultAllowZoom = await allowZoomState(page);
      console.log(`[Scenario 3 Step 2] untouched default allowZoom=${defaultAllowZoom}`);
      expect(defaultAllowZoom).toBe(false);
      const pt = await cellCenter(page, 0);
      expect(pt).not.toBeNull();
      const errBefore = consoleErrors.length;
      const pageErrBefore = pageErrors.length;
      const before = await paintedHash(page, 0);
      expect(before).not.toBeNull();
      const after = await wheelOver(page, pt!, before);
      console.log(`[Scenario 3 Step 2] scatter before=${before} after=${after}`);
      expect(after).not.toBeNull();
      expect(after).toBe(before);
      expect(consoleErrors.slice(errBefore)).toEqual([]);
      expect(pageErrors.slice(pageErrBefore)).toEqual([]);
    });

    await softStep('Scenario 3 Step 3', async () => {

      const [beforeSwitch] = await cellHashes(page, [0]);
      await setInnerType(page, 'Bar chart');

      await expect(cellLocator).toHaveCount(8);
      await assertInnerTypeSwitched(page, 'Bar chart', beforeSwitch);
      const pt = await cellCenter(page, 0);
      expect(pt).not.toBeNull();
      const errBefore = consoleErrors.length;
      const pageErrBefore = pageErrors.length;
      const before = await paintedHash(page, 0);
      expect(before).not.toBeNull();
      const after = await wheelOver(page, pt!, before);
      console.log(`[Scenario 3 Step 3] bar before=${before} after=${after}`);
      expect(after).not.toBeNull();
      expect(after).toBe(before);
      expect(consoleErrors.slice(errBefore)).toEqual([]);
      expect(pageErrors.slice(pageErrBefore)).toEqual([]);
    });

    await softStep('Scenario 3 Step 4', async () => {

      const [beforeSwitch] = await cellHashes(page, [0]);
      await setInnerType(page, 'Box plot');
      await expect(cellLocator).toHaveCount(8);
      await assertInnerTypeSwitched(page, 'Box plot', beforeSwitch);
      const pt = await cellCenter(page, 0);
      expect(pt).not.toBeNull();
      const errBefore = consoleErrors.length;
      const pageErrBefore = pageErrors.length;
      const before = await paintedHash(page, 0);
      expect(before).not.toBeNull();
      const after = await wheelOver(page, pt!, before);
      console.log(`[Scenario 3 Step 4] box before=${before} after=${after}`);
      expect(after).not.toBeNull();
      expect(after).toBe(before);
      expect(consoleErrors.slice(errBefore)).toEqual([]);
      expect(pageErrors.slice(pageErrBefore)).toEqual([]);
    });

    await softStep('Scenario 3 Step 5', async () => {

      await setInnerType(page, 'Scatter plot');
      const on = await setAllowZoom(page, true);
      expect(on).toBe(true);
      const pt = await cellCenter(page, 0);
      expect(pt).not.toBeNull();
      const errBefore = consoleErrors.length;
      const pageErrBefore = pageErrors.length;
      const before = await paintedHash(page, 0);
      expect(before).not.toBeNull();
      const after = await wheelOver(page, pt!, before);
      console.log(`[Scenario 3 Step 5] allowZoom=true before=${before} after=${after}`);
      expect(after).not.toBeNull();

      expect(after).not.toBe(before);
      expect(consoleErrors.slice(errBefore)).toEqual([]);
      expect(pageErrors.slice(pageErrBefore)).toEqual([]);
    });

    await softStep('Scenario 3 Step 6', async () => {

      const off = await setAllowZoom(page, false);
      expect(off).toBe(false);
      const pt = await cellCenter(page, 0);
      expect(pt).not.toBeNull();
      const errBefore = consoleErrors.length;
      const pageErrBefore = pageErrors.length;
      const before = await paintedHash(page, 0);
      expect(before).not.toBeNull();
      const after = await wheelOver(page, pt!, before);
      console.log(`[Scenario 3 Step 6] allowZoom=false before=${before} after=${after}`);
      expect(after).not.toBeNull();
      expect(after).toBe(before);
      expect(consoleErrors.slice(errBefore)).toEqual([]);
      expect(pageErrors.slice(pageErrBefore)).toEqual([]);
    });
  } finally {

    await setAllowZoom(page, false).catch(() => {});
  }

  await v.closeAllAndWait(page);
  v.finishSpec();
});
