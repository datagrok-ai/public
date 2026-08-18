/* ---
realizes: [trellisplot.cp.global-scale-inner-axes, trellisplot.int.global-scale-range-slider-sync]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const isBenignError = (text: string) =>
  /Failed to load resource/.test(text) || /404 \(\)/.test(text) || /favicon/.test(text) ||
  /Unable to find element in cloned iframe/.test(text);

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
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    const grid = (root.querySelector('.d4-trellis-plot-charts-grid') as HTMLElement) ?? root;
    const gr = grid.getBoundingClientRect();
    grid.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true,
      clientX: gr.left + 4, clientY: gr.top + 4}));
    await new Promise((r) => setTimeout(r, 800)); 

    const out = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
      .map((e) => (e as HTMLElement).innerText.trim());

    document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
    document.body.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
    return out;
  });
  await page.waitForTimeout(400); 
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

const innerTypeTabs = ['Scatter plot', 'Bar chart', 'Box plot', 'Histogram', 'Line chart', 'Pie chart'];

async function openInnerViewerTab(page: Page, tabName?: string): Promise<void> {

  await v.openViewerGear(page, 'Trellis plot');

  await page.evaluate(({name, names}) => {
    const headers = Array.from(document.querySelectorAll('.d4-tab-header')) as HTMLElement[];
    const tab = (name ? headers.find((h) => h.innerText.trim() === name) : undefined) ??
      headers.find((h) => names.includes(h.innerText.trim()));
    tab?.click();
  }, {name: tabName ?? '', names: innerTypeTabs});
  await page.waitForTimeout(800); 

  const row = page.locator('.property-grid tr[name="prop-allow-zoom"]');
  if (await row.count() > 0 && !(await row.isVisible())) {
    await page.locator('.property-grid tr[name="prop-category-misc"]').first().click();
    await page.waitForTimeout(500); 
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

async function wheelOver(page: Page, pt: {x: number; y: number}, steps = 5): Promise<void> {
  await page.mouse.move(pt.x, pt.y);
  await page.waitForTimeout(200); 
  for (let i = 0; i < steps; i++) {
    await page.mouse.wheel(0, 120);
    await page.waitForTimeout(60); 
  }

  await page.waitForTimeout(600);
}

test('Trellis plot: global scale, inner axes, range slider reset', async ({page}) => {
  test.setTimeout(900_000);
  page.setDefaultTimeout(120_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error' && !isBenignError(m.text())) consoleErrors.push(m.text()); });

  await loginToDatagrok(page);

  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch {}
    try { grok.shell.windows.simpleMode = true; } catch {}
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    await new Promise((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

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

    const track = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const sliders = Array.from(root.querySelectorAll(
        '.d4-range-selector > svg[type="range-slider"][name="x-slider"]')) as SVGElement[];
      if (sliders.length === 0) return null;
      const s = sliders[0];
      const r = s.getBoundingClientRect();
      s.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: r.left + r.width / 2, clientY: r.top + 5}));
      return {left: r.left, top: r.top, width: r.width, height: r.height};
    });
    expect(track).not.toBeNull();
    const t = track!;
    const y = t.top + t.height / 2;
    await page.mouse.move(t.left + t.width - 6, y);
    await page.mouse.down();
    await page.mouse.move(t.left + t.width * 0.75, y, {steps: 6});
    await page.mouse.move(t.left + t.width * 0.5, y, {steps: 8});
    await page.mouse.up();
    await v.waitForViewerRendered(page, 'Trellis plot', 900);

    const after = await cellHashes(page, probes);
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
      await new Promise((r) => setTimeout(r, 800)); 

      const target = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
        .find((e) => (e as HTMLElement).innerText.trim() === 'Reset Inner Range Sliders');
      (target?.closest('.d4-menu-item') as HTMLElement | null)?.click();
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    const after = await cellHashes(page, probes);
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
      const [before] = await cellHashes(page, [0]);
      expect(before).not.toBeNull();
      await wheelOver(page, pt!);
      const [after] = await cellHashes(page, [0]);
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
      const [before] = await cellHashes(page, [0]);
      expect(before).not.toBeNull();
      await wheelOver(page, pt!);
      const [after] = await cellHashes(page, [0]);
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
      const [before] = await cellHashes(page, [0]);
      expect(before).not.toBeNull();
      await wheelOver(page, pt!);
      const [after] = await cellHashes(page, [0]);
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
      const [before] = await cellHashes(page, [0]);
      expect(before).not.toBeNull();
      await wheelOver(page, pt!);
      const [after] = await cellHashes(page, [0]);
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
      const [before] = await cellHashes(page, [0]);
      expect(before).not.toBeNull();
      await wheelOver(page, pt!);
      const [after] = await cellHashes(page, [0]);
      console.log(`[Scenario 3 Step 6] allowZoom=false before=${before} after=${after}`);
      expect(after).not.toBeNull();
      expect(after).toBe(before);
      expect(consoleErrors.slice(errBefore)).toEqual([]);
      expect(pageErrors.slice(pageErrBefore)).toEqual([]);
    });
  } finally {

    await setAllowZoom(page, false).catch(() => {});
  }

  v.finishSpec();
});
