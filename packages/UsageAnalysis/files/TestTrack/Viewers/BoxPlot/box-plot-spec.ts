/* ---
realizes: [boxplot.cp.property-surface-smoke, boxplot.int.inside-outside-values, boxplot.int.auto-layout-hides-chrome, boxplot.int.showmarkers-gates-marker-props]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

// A stratified 1000-row subset of demog (same SEX / RACE / DIS_POP proportions): every paint of
// the plot draws one marker per row, and nothing here asserts on the row count.
const datasetPath = 'System:DemoFiles/demog-1000.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';

async function bpProp(page: Page, prop: string): Promise<any> {
  return page.evaluate((p) => {
    const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    return bp?.props?.[p];
  }, prop);
}

async function setBpProp(page: Page, prop: string, value: any, settleMs = 800): Promise<void> {
  await v.setViewerProps(page, 'Box plot', [{set: {[prop]: value}}], settleMs);
}

// One settle for the whole group: the sets coalesce into a single repaint.
async function setBpProps(page: Page, props: Record<string, any>, settleMs = 800): Promise<void> {
  await v.setViewerProps(page, 'Box plot', [{set: props}], settleMs);
}

async function dismissMenu(page: Page): Promise<void> {
  await page.evaluate(() => (window as any).__closeContextMenu());
}

async function canvasRect(page: Page): Promise<{x: number; y: number; w: number; h: number}> {
  return page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Box-plot"]')!;
    const c = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
    return {x: c.x, y: c.y, w: c.width, h: c.height};
  });
}

async function canvasInk(page: Page): Promise<number> {
  return page.evaluate(() => {
    const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    const cv = bp.root.querySelector('canvas[name="canvas"]') as HTMLCanvasElement;
    const data = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
    let n = 0;
    for (let i = 0; i < data.length; i += 4) {
      const r = data[i], g = data[i + 1], b = data[i + 2], a = data[i + 3];
      if (a !== 0 && !(r >= 250 && g >= 250 && b >= 250)) n++;
    }
    return n;
  });
}

// Clicks a leaf of the main context menu and waits for the repaint the click causes.
async function clickMainMenuLeaf(page: Page, leafName: string, renderCapMs = 900): Promise<boolean> {
  const r = await canvasRect(page);
  const clicked = await page.evaluate(async ({name, cx, cy, cap}) => {
    const w = window as any;
    const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    const cv = bp.root.querySelector('canvas[name="canvas"]') as HTMLCanvasElement;
    await w.__openContextMenu(cv, cx, cy);
    const leaf: HTMLElement | null = await w.__poll(
      () => document.querySelector(`[name="${name}"]`), (el: HTMLElement | null) => el !== null, 2000, 50);
    if (leaf) await w.__settled('viewer:Box plot.onViewerRendered', () => leaf.click(), cap);
    return !!leaf;
  }, {name: leafName, cx: r.x + r.w * 0.5, cy: r.y + r.h * 0.5, cap: renderCapMs});
  await dismissMenu(page);
  return clicked;
}

async function menuItemsAt(page: Page, fx: number, fy: number): Promise<
  {name: string | null; d4name: string | null; opacity: string}[]> {
  const r = await canvasRect(page);
  const items = await page.evaluate(async ({cx, cy}) => {
    const w = window as any;
    const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    const cv = bp.root.querySelector('canvas[name="canvas"]') as HTMLCanvasElement;
    await w.__openContextMenu(cv, cx, cy);
    return Array.from(document.querySelectorAll('.d4-menu-item')).map((i) => ({
      name: i.getAttribute('name'),
      d4name: i.getAttribute('d4-name'),
      opacity: getComputedStyle(i).opacity,
    }));
  }, {cx: r.x + r.w * fx, cy: r.y + r.h * fy});
  await dismissMenu(page);
  return items;
}

// The reveal icons appear on hover over the top strip of the canvas; the icon's laid-out width is
// the signal (no tooltip event fires for the strip). Each poll hovers a fresh point so a
// mousemove is actually delivered.
async function clickRevealIcon(page: Page, iconName: string): Promise<void> {
  const origin = await page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Box-plot"]')!;
    const c = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
    return {x: c.x, y: c.y};
  });
  let hover = 0;
  const pt = await v.pollValue(async () => {
    await page.mouse.move(origin.x + 35 + (hover++ % 4) * 5, origin.y + 15);
    return page.evaluate((name) => {
      const el = document.querySelector(`[name="${name}"]`) as HTMLElement | null;
      const r = el?.getBoundingClientRect();
      return r && r.width > 0 ? {x: r.x + r.width / 2, y: r.y + r.height / 2} : null;
    }, iconName);
  }, (p) => p !== null, 3000, 100);
  expect(pt).not.toBeNull();
  await page.mouse.click(pt!.x, pt!.y);
  await v.waitForViewerRendered(page, 'Box plot', 1200);
}

async function selectorState(page: Page, suffix: string): Promise<{display: string; w: number; h: number}> {
  return page.evaluate((sfx) => {
    const root = document.querySelector('[name="viewer-Box-plot"]')!;
    const el = root.querySelector(`[name="div-column-combobox-${sfx}"]`) as HTMLElement | null;
    if (!el) return {display: 'absent', w: 0, h: 0};
    const b = el.getBoundingClientRect();
    return {display: getComputedStyle(el).display, w: Math.round(b.width), h: Math.round(b.height)};
  }, suffix);
}

async function viewportRect(page: Page): Promise<{top: number; bottom: number; height: number}> {
  return page.evaluate(() => {
    const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    const vp = bp.viewport;
    return {top: vp.top, bottom: vp.bottom, height: vp.height};
  });
}

async function verticalSliderHandles(
  page: Page,
): Promise<{top: {x: number; y: number}; bottom: {x: number; y: number}}> {
  const readSlider = () => page.evaluate(() => {
    const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    const slider = Array.from(bp.root.querySelectorAll('svg[type="range-slider"]'))
      .find((s: any) => {
        const r = s.getBoundingClientRect();
        return r.height > r.width;
      }) as SVGElement | undefined;
    if (!slider) return null;
    const sr = slider.getBoundingClientRect();
    const circles = Array.from(slider.querySelectorAll('circle'));
    if (circles.length < 2) return null;
    const centers = circles.slice(0, 2)
      .map((c) => {
        const r = c.getBoundingClientRect();
        return {x: r.x + r.width / 2, y: r.y + r.height / 2};
      })
      .sort((a, b) => a.y - b.y);
    return {rect: {x: sr.x, y: sr.y, w: sr.width, h: sr.height},
      top: centers[0], bottom: centers[1]};
  });
  for (let i = 0; i < 20; i++) {
    const s = await readSlider();
    if (s && s.bottom.y - s.top.y > 50) return {top: s.top, bottom: s.bottom};
    if (s) {
      const fy = 0.15 + (i % 5) * 0.15;
      await page.mouse.move(s.rect.x + s.rect.w / 2, s.rect.y + s.rect.h * fy, {steps: 4});
    }
    await v.pollStable(readSlider,
      (a, b) => JSON.stringify(a) === JSON.stringify(b), 300, 100);
  }
  throw new Error('vertical range slider handles did not lay out to a usable span');
}

test('Box plot property surface smoke', async ({page}) => {
  test.setTimeout(600_000);

  const consoleErrors: string[] = [];
  const pageErrors: string[] = [];
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });
  page.on('pageerror', (e) => { pageErrors.push(String(e)); });

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await v.installEventWaits(page, {immediateRendering: true});
  await page.evaluate(() => {
    grok.shell.tv.dataFrame.name = 'demog';
    const bp = grok.shell.tv.addViewer('Box plot');
    bp.props.valueColumnName = 'AGE';
    bp.props.category1ColumnName = 'SEX';
  });
  await page.locator('[name="viewer-Box-plot"]').waitFor({timeout: 10000});
  await v.waitForViewerRendered(page, 'Box plot', 1500);
  await v.waitForViewerQuiet(page, 'Box plot');

  await softStep('[anchor: Context menus as property paths] Misc menu Show Inside/Outside Values flip the prop AND the drawn points; Markers menu Size grays with a size column', async () => {
    expect(await bpProp(page, 'showInsideValues')).toBe(true);
    const inkBoth = await canvasInk(page);
    const insideClicked = await clickMainMenuLeaf(page, 'div-Misc---Show-Inside-Values');
    expect(insideClicked).toBe(true);
    expect(await bpProp(page, 'showInsideValues')).toBe(false);
    const inkNoInside = await v.pollValue(() => canvasInk(page), (n) => n < inkBoth * 0.9, 3000, 150);
    console.log('Misc Show Inside Values ink both/off:', inkBoth, inkNoInside);
    expect(inkNoInside).toBeLessThan(inkBoth * 0.9);
    const outsideClicked = await clickMainMenuLeaf(page, 'div-Misc---Show-Outside-Values');
    expect(outsideClicked).toBe(true);
    expect(await bpProp(page, 'showOutsideValues')).toBe(false);
    const inkNeither = await v.pollValue(() => canvasInk(page), (n) => n < inkNoInside, 3000, 150);
    console.log('Misc Show Outside Values ink off:', inkNeither);
    expect(inkNeither).toBeLessThan(inkNoInside);
    await clickMainMenuLeaf(page, 'div-Misc---Show-Inside-Values', 700);
    await clickMainMenuLeaf(page, 'div-Misc---Show-Outside-Values');
    expect(await bpProp(page, 'showInsideValues')).toBe(true);
    expect(await bpProp(page, 'showOutsideValues')).toBe(true);
    const inkRestored = await v.pollValue(() => canvasInk(page), (n) => n > inkNoInside, 3000, 150);
    console.log('Misc round-trip restored ink:', inkRestored);
    expect(inkRestored).toBeGreaterThan(inkNoInside);
    await setBpProp(page, 'markerSizeColumnName', 'WEIGHT', 700);
    const markerItems = await menuItemsAt(page, 0.5, 0.5);
    const sizeItem = markerItems.find((i) => i.name === 'div-Markers---Size');
    console.log('Markers menu Size item:', JSON.stringify(sizeItem));
    expect(sizeItem).toBeTruthy();
    expect(parseFloat(sizeItem!.opacity)).toBeLessThan(1);
    await setBpProp(page, 'markerSizeColumnName', '', 500);
  });

  await softStep('[anchor: Statistics and group-comparison menu regions] Stats-region menu grays Group Comparison items while off; enabling GC ungrays them and adds an "Add ... Table" item', async () => {
    await setBpProps(page, {showStatistics: true, showGroupComparison: false, showPValue: true}, 600);

    const width = await page.evaluate(() =>
      document.querySelector('[name="viewer-Box-plot"]')!.getBoundingClientRect().width);
    expect(width).toBeGreaterThan(300);
    const statsItems = await menuItemsAt(page, 0.5, 0.93);
    const grayedGc = statsItems.find((i) => i.name === 'div-Group-Comparison---Show-Assumption-Checks');
    console.log('Stats-region Group Comparison Show Assumption Checks (GC off):', JSON.stringify(grayedGc));
    expect(grayedGc).toBeTruthy();
    expect(parseFloat(grayedGc!.opacity)).toBeLessThan(1);

    const rs = await canvasRect(page);
    const groupPt = await page.evaluate(async ({cx, cy}) => {
      const w = window as any;
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      const cv = bp.root.querySelector('canvas[name="canvas"]') as HTMLCanvasElement;
      await w.__openContextMenu(cv, cx, cy);
      const b: DOMRect | undefined = await w.__poll(
        () => document.querySelector('[name="div-Group-Comparison"]')?.getBoundingClientRect(),
        (r: DOMRect | undefined) => !!r && r.width > 0, 2000, 50);
      return b && b.width > 0 ? {x: b.x + b.width / 2, y: b.y + b.height / 2} : null;
    }, {cx: rs.x + rs.w * 0.5, cy: rs.y + rs.h * 0.93});
    expect(groupPt).not.toBeNull();
    await page.mouse.move(groupPt!.x, groupPt!.y);
    const leafPt = await page.evaluate(async () => {
      const b: DOMRect | undefined = await (window as any).__poll(
        () => document.querySelector('[name="div-Group-Comparison---Show-Assumption-Checks"]')
          ?.getBoundingClientRect(),
        (r: DOMRect | undefined) => !!r && r.width > 0, 4500, 100);
      return b && b.width > 0 ? {x: b.x + b.width / 2, y: b.y + b.height / 2} : null;
    });
    expect(leafPt).not.toBeNull();
    const hintShown = await v.armEvent(page, 'grok.events.onTooltipShown', 3000);
    await page.mouse.move(leafPt!.x, leafPt!.y);
    await hintShown();
    const hint = await v.pollValue(
      () => page.evaluate(() => (document.querySelector('.d4-tooltip')?.textContent ?? '').trim()),
      (t) => /Show Group Comparison/i.test(t), 1500, 100);
    console.log('Gated-item hover hint:', JSON.stringify(hint));

    expect(hint).toMatch(/Show Group Comparison/i);
    await dismissMenu(page);

    let pMenu: {statsFormatCount: number} | null = null;
    for (const [dx, dy] of [[60, 20], [70, 10], [64, 14], [56, 18]]) {
      const rp = await canvasRect(page);
      const shape = await page.evaluate(async ({cx, cy}) => {
        const w = window as any;
        const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
        const cv = bp.root.querySelector('canvas[name="canvas"]') as HTMLCanvasElement;
        await w.__openContextMenu(cv, cx, cy);
        const b: DOMRect | undefined = await w.__poll(
          () => document.querySelector('[name="div-Show-P-Value"]')?.getBoundingClientRect(),
          (r: DOMRect | undefined) => !!r && r.width > 0, 500, 50);
        return b && b.width > 0
          ? {statsFormatCount: document.querySelectorAll(
            '[name="div-Statistics-Format"], [name="div-Statistics---Statistics-Format"]').length}
          : null;
      }, {cx: rp.x + dx, cy: rp.y + dy});
      await dismissMenu(page);
      if (shape) { pMenu = shape; break; }
    }
    console.log('P-value region exclusive menu (Statistics Format count):', JSON.stringify(pMenu));
    expect(pMenu).not.toBeNull();
    expect(pMenu!.statsFormatCount).toBe(0);
    await clickRevealIcon(page, 'show-group-stats');
    expect(await bpProp(page, 'showGroupComparison')).toBe(true);

    let addItem: string | null = null;
    let ungrayed: {name: string | null; d4name: string | null; opacity: string} | undefined;
    for (const [dx, dy] of [[42, 16], [50, 16], [60, 14], [80, 14], [42, 30]]) {
      const rr = await canvasRect(page);
      const stripItems = await menuItemsAt(page, dx / rr.w, dy / rr.h);
      const add = stripItems.find((i) => /^Add .*Table$/i.test(i.d4name ?? ''));
      ungrayed = stripItems.find((i) => i.name === 'div-Show-Assumption-Checks')
        ?? stripItems.find((i) => i.d4name === 'Show Assumption Checks');
      if (add) { addItem = add.d4name; break; }
    }
    console.log('Comparison-strip Add-Table item / ungrayed Assumption Checks:', addItem, JSON.stringify(ungrayed));
    expect(addItem).toBeTruthy();
    expect(ungrayed).toBeTruthy();
    expect(parseFloat(ungrayed!.opacity)).toBe(1);
    await setBpProp(page, 'showGroupComparison', false, 600);
  });

  await softStep('[anchor: Resize and auto layout] Auto Layout hides the column selectors at a small size and restores them; a narrow resize with a coloring raises no error (GROK-18677)', async () => {
    await setBpProps(page, {valueColumnName: 'AGE', category1ColumnName: 'SEX', autoLayout: true}, 500);
    const visibleSelectors = () => page.evaluate(() =>
      Array.from(document.querySelectorAll('[name="viewer-Box-plot"] [name^="div-column-combobox-"]'))
        .filter((s) => { const b = s.getBoundingClientRect(); return b.width > 0 && b.height > 0; }).length);
    const largeCount = await v.pollValue(visibleSelectors, (n) => n > 1, 3000, 150);
    console.log('Auto-layout visible selectors (large):', largeCount);
    expect(largeCount).toBeGreaterThan(1);
    await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Box-plot"]') as HTMLElement;
      (window as any).__bpOrigSize = {w: root.style.width, h: root.style.height};
      root.style.width = '170px'; root.style.height = '150px';
      window.dispatchEvent(new Event('resize'));
    });
    await v.waitForViewerRendered(page, 'Box plot', 1500);
    const smallCount = await v.pollValue(visibleSelectors, (n) => n < largeCount, 3000, 150);
    console.log('Auto-layout visible selectors (small):', smallCount);
    expect(smallCount).toBeLessThan(largeCount);
    await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Box-plot"]') as HTMLElement;
      const s = (window as any).__bpOrigSize ?? {w: '', h: ''};
      root.style.width = s.w; root.style.height = s.h;
      window.dispatchEvent(new Event('resize'));
    });
    await v.waitForViewerRendered(page, 'Box plot', 1200);
    const restoredCount = await v.pollValue(visibleSelectors, (n) => n > smallCount, 3000, 150);
    console.log('Auto-layout visible selectors (restored):', restoredCount);
    expect(restoredCount).toBeGreaterThan(smallCount);
    await setBpProp(page, 'markerColorColumnName', 'SEX', 800);
    const errBefore = consoleErrors.length;
    const pageErrBefore = pageErrors.length;
    await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Box-plot"]') as HTMLElement;
      root.style.width = '120px';
      window.dispatchEvent(new Event('resize'));
    });
    await v.waitForViewerRendered(page, 'Box plot', 700);
    await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Box-plot"]') as HTMLElement;
      const s = (window as any).__bpOrigSize ?? {w: '', h: ''};
      root.style.width = s.w;
      window.dispatchEvent(new Event('resize'));
    });
    await v.waitForViewerRendered(page, 'Box plot', 1000);
    const errDelta = consoleErrors.slice(errBefore);
    const pageErrDelta = pageErrors.slice(pageErrBefore);
    console.log('GROK-18677 narrow-resize error deltas:', JSON.stringify(errDelta), JSON.stringify(pageErrDelta));
    expect(errDelta).toEqual([]);
    expect(pageErrDelta).toEqual([]);
    await setBpProp(page, 'markerColorColumnName', '', 500);
  });

  // a settled set (setViewerProps) leaves the canvas final: one property change is one repaint
  await softStep('[anchor: Marker gate and size scaling] Disabling Show Markers removes the points and grays the Marker group; Size Scaling linear/log repaints the markers', async () => {
    const inkWithMarkers = await canvasInk(page);
    await v.snapshotCanvasColors(page, 'Box plot');
    await setBpProp(page, 'showMarkers', false, 300);
    const markerGateDelta = await v.waitForCanvasChange(page, 'Box plot', {minDelta: 2000, timeoutMs: 15000});
    console.log('Show Markers off canvas delta:', markerGateDelta);
    expect(markerGateDelta).toBeGreaterThanOrEqual(0);
    expect(markerGateDelta).toBeGreaterThan(2000);
    const inkNoMarkers = await canvasInk(page);
    expect(inkNoMarkers).toBeLessThan(inkWithMarkers);
    await page.evaluate(() => { grok.shell.o = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot'); });
    const markerRowOpacity = await page.evaluate(async () => {
      const row: HTMLElement | null = await (window as any).__poll(
        () => document.querySelector('.property-grid tr[name="prop-marker-type"]'),
        (el: HTMLElement | null) => el !== null, 10000, 250);
      return row ? getComputedStyle(row).opacity : null;
    });
    console.log('Marker group prop-marker-type opacity (showMarkers off):', markerRowOpacity);
    expect(markerRowOpacity).not.toBeNull();
    expect(parseFloat(markerRowOpacity as string)).toBeLessThan(1);
    await setBpProp(page, 'showMarkers', true, 800);
    const inkReturned = await v.pollValue(() => canvasInk(page), (n) => n > inkNoMarkers, 3000, 150);
    expect(inkReturned).toBeGreaterThan(inkNoMarkers);
    await setBpProp(page, 'markerSizeColumnName', 'WEIGHT', 800);
    await v.snapshotCanvasColors(page, 'Box plot');
    await setBpProp(page, 'markerSizeScaling', 'logarithmic', 300);
    const scalingDelta = await v.waitForCanvasChange(page, 'Box plot', {minDelta: 1, timeoutMs: 15000});
    console.log('Size Scaling linear→log canvas delta:', scalingDelta);
    expect(scalingDelta).toBeGreaterThanOrEqual(0);
    expect(scalingDelta).toBeGreaterThan(0);
    await v.snapshotCanvasColors(page, 'Box plot');
    await setBpProp(page, 'markerSizeScaling', 'linear', 300);
    const scalingBackDelta = await v.waitForCanvasChange(page, 'Box plot', {minDelta: 1, timeoutMs: 15000});
    console.log('Size Scaling log→linear canvas delta:', scalingBackDelta);
    expect(scalingBackDelta).toBeGreaterThanOrEqual(0);
    expect(scalingBackDelta).toBeGreaterThan(0);
    await setBpProp(page, 'markerSizeColumnName', '', 500);
  });

  await softStep('[anchor: Whisker and control-band style] Whisker line width / width ratio each repaint; Control Band Color sets without error', async () => {
    await v.snapshotCanvasColors(page, 'Box plot');
    await setBpProp(page, 'whiskerLineWidth', 4, 300);
    const lwDelta = await v.waitForCanvasChange(page, 'Box plot', {minDelta: 1, timeoutMs: 15000});
    console.log('Whisker Line Width canvas delta:', lwDelta);
    expect(lwDelta).toBeGreaterThanOrEqual(0);
    expect(lwDelta).toBeGreaterThan(0);
    await v.snapshotCanvasColors(page, 'Box plot');
    await setBpProp(page, 'whiskerWidthRatio', 0.3, 300);
    const wrDelta = await v.waitForCanvasChange(page, 'Box plot', {minDelta: 1, timeoutMs: 15000});
    console.log('Whisker Width Ratio canvas delta:', wrDelta);
    expect(wrDelta).toBeGreaterThanOrEqual(0);
    expect(wrDelta).toBeGreaterThan(0);

    const errBefore = consoleErrors.length;
    const pageErrBefore = pageErrors.length;
    await setBpProp(page, 'controlBandColor', 0xFF00AA00, 600);
    const errDelta = consoleErrors.slice(errBefore);
    const pageErrDelta = pageErrors.slice(pageErrBefore);
    console.log('Control Band Color set error deltas:', JSON.stringify(errDelta), JSON.stringify(pageErrDelta));
    expect(errDelta).toEqual([]);
    expect(pageErrDelta).toEqual([]);
    await setBpProps(page, {whiskerLineWidth: 2, whiskerWidthRatio: 0.5}, 400);
  });

  await softStep('[anchor: Controls visibility] Size selector absent by default; each visibility toggle adds/removes its chrome (DOM for the selectors, canvas repaint for the canvas-drawn chrome); the round trip restores the default baseline', async () => {
    console.log('showSizeSelector default:', await bpProp(page, 'showSizeSelector'));
    const sizeBaseline = await selectorState(page, 'marker--size');
    console.log('Size selector baseline state:', JSON.stringify(sizeBaseline));
    expect(sizeBaseline.display).toBe('none');
    await setBpProp(page, 'showSizeSelector', true, 700);
    const sizeOn = await v.pollValue(() => selectorState(page, 'marker--size'),
      (s) => s.display !== 'none' && s.w > 0, 3000, 150);
    console.log('Size selector after enable:', JSON.stringify(sizeOn));
    expect(sizeOn.display).not.toBe('none');
    expect(sizeOn.w).toBeGreaterThan(0);
    await setBpProp(page, 'showValueSelector', false, 700);
    const valueOff = await v.pollValue(() => selectorState(page, 'value'),
      (s) => s.display === 'none', 3000, 150);
    console.log('Value selector after disable:', JSON.stringify(valueOff));
    expect(valueOff.display).toBe('none');
    await setBpProp(page, 'showColorSelector', false, 700);
    const colorOff = await v.pollValue(() => selectorState(page, 'marker--color'),
      (s) => s.display === 'none', 3000, 150);
    console.log('Color selector after disable:', JSON.stringify(colorOff));
    expect(colorOff.display).toBe('none');

    const canvasToggle = async (prop: string, value: boolean) => {
      await v.snapshotCanvasColors(page, 'Box plot');
      await setBpProp(page, prop, value, 300);
      const delta = await v.waitForCanvasChange(page, 'Box plot', {minDelta: 2000, timeoutMs: 15000});
      console.log(`${prop}=${value} canvas delta:`, delta);
      expect(delta).toBeGreaterThanOrEqual(0);
      expect(delta).toBeGreaterThan(2000);
    };
    await canvasToggle('showCategorySelector', false);
    await canvasToggle('showValueAxis', false);
    await canvasToggle('showCategoryAxis', false);
    await canvasToggle('showCategoryAxis', true);
    await canvasToggle('showValueAxis', true);
    await canvasToggle('showCategorySelector', true);
    await setBpProps(page, {showValueSelector: true, showColorSelector: true, showSizeSelector: false}, 700);
    const sizeRestored = await v.pollValue(() => selectorState(page, 'marker--size'),
      (s) => s.display === 'none', 3000, 150);
    const valueRestored = await v.pollValue(() => selectorState(page, 'value'),
      (s) => s.display !== 'none', 3000, 150);
    const colorRestored = await v.pollValue(() => selectorState(page, 'marker--color'),
      (s) => s.display !== 'none', 3000, 150);
    console.log('Controls restored size/value/color:', JSON.stringify(sizeRestored),
      JSON.stringify(valueRestored), JSON.stringify(colorRestored));
    expect(sizeRestored.display).toBe('none');
    expect(valueRestored.display).not.toBe('none');
    expect(colorRestored.display).not.toBe('none');
  });

  await softStep('[anchor: Title and description] Title text appears in the panel titlebar; description appears while Always, moves to Bottom, and disappears while Never', async () => {
    await setBpProps(page, {showTitle: true, title: 'Age by Race'}, 800);
    const readTitle = () => page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Box-plot"]')!;
      const panel = root.closest('.panel-base');
      return (panel?.querySelector('.panel-titlebar-text')?.textContent ?? '').trim();
    });
    const titleText = await v.pollValue(readTitle, (t) => t === 'Age by Race', 3000, 150);
    console.log('Panel titlebar text:', JSON.stringify(titleText));
    expect(titleText).toBe('Age by Race');
    await setBpProps(page, {description: 'Box plot of patient ages', descriptionVisibilityMode: 'Always'}, 700);
    const readDesc = () => page.evaluate(() => {
      const el = document.querySelector('[name="viewer-Box-plot"] .d4-viewer-description');
      return el ? (el.textContent ?? '').trim() : null;
    });
    const descAlways = await v.pollValue(readDesc, (d) => d === 'Box plot of patient ages', 3000, 150);
    console.log('Description (Always):', JSON.stringify(descAlways));
    expect(descAlways).toBe('Box plot of patient ages');
    await setBpProp(page, 'descriptionPosition', 'Bottom', 700);
    const descBottom = await v.pollValue(readDesc, (d) => d !== null, 3000, 150) !== null;
    expect(descBottom).toBe(true);
    await setBpProp(page, 'descriptionVisibilityMode', 'Never', 700);
    const descNever = await v.pollValue(readDesc, (d) => d === null, 3000, 150) !== null;
    console.log('Description host present when Never:', descNever);
    expect(descNever).toBe(false);
    await setBpProps(page, {showTitle: false, title: '', description: '', descriptionVisibilityMode: 'Auto',
      descriptionPosition: 'Top'}, 400);
  });

  await softStep('[anchor: Axis font] Changing Axis Font repaints the labels with no error — no Infinity.floor (GROK-19297); restoring completes without error', async () => {
    const original = await bpProp(page, 'axisFont');
    const errBefore = consoleErrors.length;
    const pageErrBefore = pageErrors.length;
    await setBpProp(page, 'axisFont', 'normal normal 16px "Roboto"', 800);
    expect(await bpProp(page, 'axisFont')).toBe('normal normal 16px "Roboto"');
    await setBpProp(page, 'axisFont', original, 700);
    const errDelta = consoleErrors.slice(errBefore);
    const pageErrDelta = pageErrors.slice(pageErrBefore);
    console.log('GROK-19297 axis-font error deltas:', JSON.stringify(errDelta), JSON.stringify(pageErrDelta));
    expect(errDelta.filter((e) => /Infinity\.floor/i.test(e))).toEqual([]);
    expect(errDelta).toEqual([]);
    expect(pageErrDelta).toEqual([]);
  });

  await softStep('[anchor: Date category mapping] Category 1 Map Month then Quarter restructures the datetime category axis; returning to a categorical column restores plain categories', async () => {
    await setBpProp(page, 'category1ColumnName', 'STARTED', 1200);
    expect(await bpProp(page, 'category1ColumnName')).toBe('STARTED');
    await v.snapshotCanvasColors(page, 'Box plot');
    await setBpProp(page, 'category1Map', 'month', 300);
    const monthDelta = await v.waitForCanvasChange(page, 'Box plot', {minDelta: 1, timeoutMs: 15000});
    console.log('Category1Map month canvas delta:', monthDelta);
    expect(monthDelta).toBeGreaterThanOrEqual(0);
    expect(monthDelta).toBeGreaterThan(0);
    expect(await bpProp(page, 'category1Map')).toBe('month');
    await v.snapshotCanvasColors(page, 'Box plot');
    await setBpProp(page, 'category1Map', 'quarter', 300);
    const quarterDelta = await v.waitForCanvasChange(page, 'Box plot', {minDelta: 1, timeoutMs: 15000});
    console.log('Category1Map quarter canvas delta:', quarterDelta);
    expect(quarterDelta).toBeGreaterThanOrEqual(0);
    expect(quarterDelta).toBeGreaterThan(0);
    expect(await bpProp(page, 'category1Map')).toBe('quarter');
    const errBefore = consoleErrors.length;
    await setBpProp(page, 'category1ColumnName', 'RACE', 1000);
    expect(await bpProp(page, 'category1ColumnName')).toBe('RACE');
    console.log('Date-mapping restore error delta:', JSON.stringify(consoleErrors.slice(errBefore)));
    expect(consoleErrors.slice(errBefore)).toEqual([]);
    await setBpProp(page, 'category1ColumnName', 'SEX', 800);
  });

  await softStep('[anchor: Custom tooltip] Row Tooltip AGE, SEX, WEIGHT shows exactly those three columns on marker hover; resetting to inherit restores the default tooltip', async () => {
    await setBpProps(page, {valueColumnName: 'AGE', category1ColumnName: 'RACE', markerSize: 10,
      rowTooltip: 'AGE\nSEX\nWEIGHT', showTooltip: 'show custom tooltip'}, 600);
    expect(await bpProp(page, 'rowTooltip')).toBe('AGE\nSEX\nWEIGHT');
    const readTipCols = () => page.evaluate(() => {
      const tip = document.querySelector('.d4-tooltip table.d4-row-tooltip-table');
      if (!tip) return [];
      return Array.from(tip.querySelectorAll('tr'))
        .map((tr) => ((tr as HTMLTableRowElement).cells[0]?.textContent ?? '').trim())
        .filter((t) => t.length > 0);
    });
    const r = await canvasRect(page);

    // The tooltip host keeps its table in the DOM and only hides, so "gone" is the platform's
    // onTooltipClosed on leaving the viewer, confirmed by the host's display.
    const tooltipHidden = () => page.evaluate(() => {
      const t = document.querySelector('.d4-tooltip') as HTMLElement | null;
      return !t || getComputedStyle(t).display === 'none';
    });
    const moveAwayAndHideTooltip = async () => {
      const closed = await v.armEvent(page, 'grok.events.onTooltipClosed', 1000);
      await page.mouse.move(r.x + r.w * 0.5, r.y - 40);
      await closed();
      expect(await v.pollValue(tooltipHidden, (hidden) => hidden, 500, 50)).toBe(true);
    };
    let tipCols: string[] = [];
    for (const [fx, fy] of [[0.62, 0.55], [0.3, 0.5], [0.5, 0.55], [0.4, 0.45], [0.6, 0.6], [0.5, 0.4]]) {
      await page.mouse.move(r.x + r.w * fx, r.y + r.h * fy);
      tipCols = await v.pollValue(readTipCols, (c) => c.length > 0, 500, 100);
      if (tipCols.length > 0) break;
    }
    console.log('Custom tooltip columns:', JSON.stringify(tipCols));
    const uniqueCols = Array.from(new Set(tipCols.map((c) => c.toUpperCase())));
    expect(uniqueCols.sort()).toEqual(['AGE', 'SEX', 'WEIGHT']);
    await moveAwayAndHideTooltip();
    await setBpProps(page, {showTooltip: 'inherit from table', rowTooltip: ''}, 500);
    expect(await bpProp(page, 'rowTooltip')).toBe('');
    let defaultCols: string[] = [];
    for (const [fx, fy] of [[0.62, 0.55], [0.3, 0.5], [0.5, 0.55], [0.4, 0.45], [0.6, 0.6]]) {
      await page.mouse.move(r.x + r.w * fx, r.y + r.h * fy);
      defaultCols = await v.pollValue(readTipCols, (c) => c.length > 0, 500, 100);
      if (defaultCols.length > 0) break;
    }
    const defaultUnique = Array.from(new Set(defaultCols.map((c) => c.toUpperCase())));
    console.log('Default (inherited) tooltip columns:', JSON.stringify(defaultUnique));
    expect(defaultUnique.length).toBeGreaterThan(0);
    expect(defaultUnique.sort()).not.toEqual(['AGE', 'SEX', 'WEIGHT']);
    await moveAwayAndHideTooltip();
  });

  await softStep('[anchor: Table switching resets Category 2] Switching Table demog→spgi-100 with a two-level category resets Category 2 to a consistent state — no stale demog column (GROK-18361)', async () => {
    await page.evaluate(async (path) => {
      const df = await grok.dapi.files.readCsv(path);
      df.name = 'spgi-100';
      grok.shell.addTableView(df);
    }, spgiPath);
    await v.pollValue(() => page.evaluate(() =>
      Array.from(grok.shell.tableViews).some((vw: any) => vw.dataFrame?.name === 'spgi-100')),
    (present) => present, 1500, 100);
    await page.evaluate(() => {
      const home = Array.from(grok.shell.tableViews).find((vw: any) => vw.dataFrame?.name === 'demog');
      grok.shell.v = home;
    });
    await page.locator('[name="viewer-Box-plot"]').waitFor({timeout: 10000});
    await v.pollValue(() => page.evaluate(() =>
      grok.shell.tv?.dataFrame?.name === 'demog'
      && !!grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot')),
    (ready) => ready, 600, 100);
    await setBpProps(page, {valueColumnName: 'AGE', category1ColumnName: 'SEX', category2ColumnName: 'RACE'}, 800);
    expect(await bpProp(page, 'category2ColumnName')).toBe('RACE');
    const errBefore = consoleErrors.length;
    const pageErrBefore = pageErrors.length;
    await page.evaluate(async () => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      await (window as any).__settled('viewer:Box plot.onViewerRendered', () => { bp.props.table = 'spgi-100'; }, 2000);
    });

    const afterCat2 = await bpProp(page, 'category2ColumnName');
    const spgiHasCat2 = await page.evaluate((c) => {
      const t = grok.shell.tables.find((tb: any) => tb.name === 'spgi-100');
      return c == null || c === '' || t.columns.names().includes(c);
    }, afterCat2);
    console.log('GROK-18361 Category2 after table switch:', JSON.stringify(afterCat2));
    expect(afterCat2).not.toBe('RACE');
    expect(spgiHasCat2).toBe(true);
    const errDelta = consoleErrors.slice(errBefore);
    const pageErrDelta = pageErrors.slice(pageErrBefore);
    console.log('GROK-18361 table-switch error deltas:', JSON.stringify(errDelta), JSON.stringify(pageErrDelta));
    expect(errDelta).toEqual([]);
    expect(pageErrDelta).toEqual([]);
    await setBpProps(page, {valueColumnName: 'Average Mass', category1ColumnName: 'Series'}, 800);
    const afterDf = await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      return bp.dataFrame?.name;
    });
    expect(afterDf).toBe('spgi-100');
    await setBpProp(page, 'table', 'demog', 1500);
    await setBpProps(page, {valueColumnName: 'AGE', category1ColumnName: 'SEX', category2ColumnName: ''}, 600);
  });

  await softStep('[anchor: Legend minimum under coloring] A legend-bearing coloring keeps the render valid: canvas keeps ink, warnings delta zero, markerColorColumnName stays applied', async () => {
    const errBefore = consoleErrors.length;
    const pageErrBefore = pageErrors.length;
    await setBpProp(page, 'markerColorColumnName', 'RACE', 1000);
    expect(await bpProp(page, 'markerColorColumnName')).toBe('RACE');
    const px = await v.pollValue(() => v.countCanvasPixels(page, 'Box plot'),
      (c) => c.total > 0, 3000, 150);
    console.log('Legend-minimum canvas pixels:', px.total);
    expect(px.total).toBeGreaterThan(0);
    const errDelta = consoleErrors.slice(errBefore);
    const pageErrDelta = pageErrors.slice(pageErrBefore);
    console.log('Legend-minimum error deltas:', JSON.stringify(errDelta), JSON.stringify(pageErrDelta));
    expect(errDelta).toEqual([]);
    expect(pageErrDelta).toEqual([]);
    await setBpProp(page, 'markerColorColumnName', '', 600);
  });

  await softStep('[anchor: Double-click resets the view] Range-slider zoom narrows the viewport; double-clicking empty plot space fires d4-boxplot-reset-view AND restores the full range', async () => {
    const vpFull = await viewportRect(page);
    const handles = await verticalSliderHandles(page);
    await page.mouse.move(handles.top.x, handles.top.y);
    await page.mouse.down();
    // the slider repaints the plot on every mousemove, so each step costs a paint
    await page.mouse.move(handles.top.x,
      handles.top.y + (handles.bottom.y - handles.top.y) * 0.4, {steps: 4});
    await page.mouse.up();
    const vpZoomed = await v.pollValue(() => viewportRect(page),
      (vp) => vp.height < vpFull.height * 0.95, 900, 100);
    console.log('Viewport full/zoomed:', JSON.stringify(vpFull), JSON.stringify(vpZoomed));
    expect(vpZoomed.height).toBeLessThan(vpFull.height * 0.95);
    await page.evaluate(() => {
      const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
      (window as any).__bpResetFired = false;
      (window as any).__bpResetSub = bp.onEvent('d4-boxplot-reset-view')
        .subscribe(() => { (window as any).__bpResetFired = true; });
    });

    const r = await canvasRect(page);
    await page.mouse.dblclick(r.x + r.w * 0.5, r.y + r.h * 0.05);
    await v.pollValue(() => page.evaluate(() => !!(window as any).__bpResetFired),
      (fired) => fired, 700, 100);
    const resetFired = await page.evaluate(() => {
      const fired = (window as any).__bpResetFired;
      try { (window as any).__bpResetSub?.unsubscribe(); } catch (_) {  }
      return fired;
    });
    const vpReset = await viewportRect(page);
    console.log('d4-boxplot-reset-view fired on dblclick:', resetFired,
      'viewport after reset:', JSON.stringify(vpReset));
    expect(resetFired).toBe(true);
    expect(vpReset.height).toBeCloseTo(vpFull.height, 1);
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});
