/* ---
realizes: [trellisplot.cp.split-and-pick-inner, trellisplot.cp.click-to-filter, trellisplot.cp.global-scale-inner-axes, trellisplot.cp.scroll-categories, trellisplot.cp.tiled-single-column, trellisplot.int.selectors-labels-visibility-coupling, trellisplot.int.undo-redo-viewer-lifecycle]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const curvesPath = 'System:DemoFiles/curves.csv';

function axisViewportCount(n: number, oneColumnOnly: boolean): number {
  return Math.min(oneColumnOnly ? 5 : (n < 5 * 1.5 ? n : 5), n);
}

async function cellIndexFor(page: Page, xValue: string, yValue: string): Promise<number> {
  return page.evaluate(({xValue, yValue}) => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    const df = grok.shell.tv.dataFrame;
    const xCol = df.col(tp.props.xColumnNames[0]);
    const yCol = df.col(tp.props.yColumnNames[0]);
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
    const byTop: Record<string, number> = {};
    for (const c of Array.from(root?.querySelectorAll('.d4-trellis-plot-cell') ?? [])) {
      const t = String(Math.round(c.getBoundingClientRect().top));
      byTop[t] = (byTop[t] || 0) + 1;
    }
    const perRow = Object.keys(byTop).length ? Math.max(...Object.values(byTop)) : 0;
    const stride = perRow > 0 ? Math.min(perRow, tp.xCategoriesCount) : tp.xCategoriesCount;
    return yCol.categories.indexOf(yValue) * stride + xCol.categories.indexOf(xValue);
  }, {xValue, yValue});
}

async function comboRowCount(page: Page, xCol: string, xValue: string, yCol: string, yValue: string): Promise<number> {
  return page.evaluate(({xCol, xValue, yCol, yValue}) => {
    const df = grok.shell.tv.dataFrame;
    const xc = df.col(xCol), yc = df.col(yCol);
    let n = 0;
    for (let i = 0; i < df.rowCount; i++)
      if (xc.get(i) === xValue && yc.get(i) === yValue) n++;
    return n;
  }, {xCol, xValue, yCol, yValue});
}

// Both delegate to the installed in-page __menuLeaf (the same navigator drivePanelMenuLeaf /
// driveContextMenuLeaf use) — a submenu tree is DOM-present but zero-size until a real pointer
// hover expands each level, which __menuLeaf performs. The menu is opened by the caller's real
// right-click; these only navigate the already-open .d4-menu-popup.
async function clickMenuItemInGroup(page: Page, group: string, item: string): Promise<void> {
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 10000});
  await page.evaluate(({g, l}) => (window as any).__menuLeaf(g, l), {g: group, l: item});
}

async function clickTopLevelMenuItem(page: Page, item: string): Promise<void> {
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 10000});
  await page.evaluate((l) => (window as any).__menuLeaf(null, l), item);
}

async function ensureStyleExpanded(page: Page): Promise<void> {
  // The Style section is collapsed by default and its rows are zero-height until expanded
  // (refdoc: pitfall 16); ensurePropertyCategory expands the owning category and settles on the
  // probe row becoming visible, replacing the fixed 600ms wait after a header click.
  await v.ensurePropertyCategory(page, 'Trellis-plot', 'style', 'show-all-categories').catch(() => {});
}

async function showAllCategoriesState(page: Page): Promise<boolean | null> {
  await ensureStyleExpanded(page);
  const val = await v.propertyGridValue(page, 'show-all-categories', 'style');
  return val === '' ? null : val === 'true';
}

async function openBoxPlotTab(page: Page): Promise<void> {
  await v.openViewerGear(page, 'Trellis plot');
  await page.locator('.property-grid').first().waitFor({state: 'visible', timeout: 10000}).catch(() => {});
  await page.evaluate(() => {
    const tab = Array.from(document.querySelectorAll('.d4-tab-header'))
      .find((h) => (h as HTMLElement).innerText.trim() === 'Box plot') as HTMLElement | undefined;
    tab?.click();
  });
  // the Box plot tab renders its property grid asynchronously; settle on the show-all-categories
  // row appearing rather than a fixed 800ms
  await page.locator('.property-grid tr[name="prop-show-all-categories"]').first()
    .waitFor({state: 'attached', timeout: 10000}).catch(() => {});
}

async function setShowAllCategories(page: Page, desired: boolean): Promise<boolean | null> {
  await ensureStyleExpanded(page);
  const box = page.locator('.property-grid tr[name="prop-show-all-categories"] input[type="checkbox"]');
  if (await box.count() === 0) return null;
  // CLICK the checkbox (never setOptions — a props write lands in the saved look and makes the
  // GROK-20432 guard true by construction); setPropertyGridCheckbox clicks and settles on the
  // read-back matching desired, replacing the 900ms sleep.
  await v.setPropertyGridCheckbox(page, 'show-all-categories', desired, 'style').catch(() => {});
  return box.isChecked();
}

async function legendCategories(page: Page): Promise<string[]> {
  await page.waitForFunction(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
    return !!root && root.querySelectorAll('[name="legend"] .d4-legend-item .d4-legend-cross').length > 0;
  }, {timeout: 20000});
  return page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    return Array.from(root.querySelectorAll('[name="legend"] .d4-legend-item'))
      .map((it) => (it.querySelector('.d4-legend-value')?.textContent ?? '').trim())
      .filter((v) => v.length > 0);
  });
}

interface LegendItemState {clicked: boolean; present: boolean; active: boolean; opacity: number}

async function legendItemState(page: Page, category: string): Promise<LegendItemState> {
  return page.evaluate((cat) => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
    const item = root && Array.from(root.querySelectorAll('[name="legend"] .d4-legend-item'))
      .find((it) => (it.querySelector('.d4-legend-value')?.textContent ?? '').trim() === cat) as HTMLElement | undefined;
    if (!item) return {clicked: false, present: false, active: false, opacity: -1};
    return {
      clicked: false,
      present: true,
      active: item.classList.contains('d4-legend-item-current'),
      opacity: parseFloat(getComputedStyle(item).opacity),
    };
  }, category);
}

async function uncheckLegendCategory(page: Page, category: string): Promise<LegendItemState> {
  const clicked = await page.evaluate((cat) => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    const item = Array.from(root.querySelectorAll('[name="legend"] .d4-legend-item'))
      .find((it) => (it.querySelector('.d4-legend-value')?.textContent ?? '').trim() === cat);
    const cross = item?.querySelector('.d4-legend-cross') as HTMLElement | null;
    if (!cross) return false;
    cross.dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true, view: window, button: 0}));
    return true;
  }, category);
  if (!clicked) return {clicked: false, present: false, active: false, opacity: -1};
  await page.waitForFunction((cat) => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
    const item = root && Array.from(root.querySelectorAll('[name="legend"] .d4-legend-item'))
      .find((it) => (it.querySelector('.d4-legend-value')?.textContent ?? '').trim() === cat) as HTMLElement | undefined;
    return !!item && parseFloat(getComputedStyle(item).opacity) < 1;
  }, category, {timeout: 15000});
  return {...(await legendItemState(page, category)), clicked: true};
}

async function selectOnlyLegendCategory(page: Page, category: string): Promise<LegendItemState> {
  const clicked = await page.evaluate((cat) => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    const item = Array.from(root.querySelectorAll('[name="legend"] .d4-legend-item'))
      .find((it) => (it.querySelector('.d4-legend-value')?.textContent ?? '').trim() === cat);
    const label = item?.querySelector('.d4-legend-value') as HTMLElement | null;
    if (!label) return false;
    label.dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true, view: window, button: 0}));
    return true;
  }, category);
  if (!clicked) return {clicked: false, present: false, active: false, opacity: -1};
  await page.waitForFunction((cat) => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
    if (!root) return false;
    const items = Array.from(root.querySelectorAll('[name="legend"] .d4-legend-item')) as HTMLElement[];
    const name = (it: HTMLElement) => (it.querySelector('.d4-legend-value')?.textContent ?? '').trim();
    const shown = items.filter((it) => parseFloat(getComputedStyle(it).opacity) >= 1);
    return shown.length === 1 && name(shown[0]) === cat &&
      shown[0].classList.contains('d4-legend-item-current');
  }, category, {timeout: 15000});
  return {...(await legendItemState(page, category)), clicked: true};
}

async function legendSnapshot(page: Page): Promise<{name: string; opacity: number; current: boolean}[]> {
  return page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
    if (!root) return [];

    return Array.from(root.querySelectorAll('[name="legend"] .d4-legend-item')).map((it) => ({
      name: (it.querySelector('.d4-legend-value')?.textContent ?? '').trim(),
      opacity: parseFloat(getComputedStyle(it).opacity),
      current: it.classList.contains('d4-legend-item-current'),
    })).filter((s) => s.name.length > 0);
  });
}

async function resetLegendViaLastCross(page: Page, category: string): Promise<boolean> {
  const clicked = await page.evaluate((cat) => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    const item = Array.from(root.querySelectorAll('[name="legend"] .d4-legend-item'))
      .find((it) => (it.querySelector('.d4-legend-value')?.textContent ?? '').trim() === cat);
    const cross = item?.querySelector('.d4-legend-cross') as HTMLElement | null;
    if (!cross) return false;
    cross.dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true, view: window, button: 0}));
    return true;
  }, category);
  if (!clicked) return false;
  await page.waitForFunction(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
    if (!root) return false;
    const items = Array.from(root.querySelectorAll('[name="legend"] .d4-legend-item')) as HTMLElement[];
    return items.length > 0 && items.every((it) =>
      parseFloat(getComputedStyle(it).opacity) >= 1 && !it.classList.contains('d4-legend-item-current'));
  }, undefined, {timeout: 15000});
  return true;
}

async function trellisCount(page: Page): Promise<number> {
  return page.evaluate(() =>
    grok.shell.tv.viewers.filter((x: any) => x.type === 'Trellis plot').length);
}

async function waitForTrellisCount(page: Page, target: number, capMs: number): Promise<void> {
  await page.waitForFunction(
    (t) => grok.shell.tv.viewers.filter((x: any) => x.type === 'Trellis plot').length === t,
    target, {timeout: capMs}).catch(() => {});
}

async function closeTrellis(page: Page): Promise<void> {
  await page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]');
    const panel = root?.closest('.panel-base');
    const btn = panel?.querySelector('.panel-titlebar [name="Close"]') as HTMLElement | null;
    btn?.click();
  });
}

async function gridGeometry(page: Page): Promise<{count: number; rows: number; cols: number; maxPerRow: number}> {
  return page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
    if (!root) return {count: -1, rows: -1, cols: -1, maxPerRow: -1};
    const cells = Array.from(root.querySelectorAll('.d4-trellis-plot-cell'));
    const byTop: Record<string, number> = {};
    const byLeft: Record<string, number> = {};
    for (const c of cells) {
      const b = c.getBoundingClientRect();
      const t = String(Math.round(b.top));
      const l = String(Math.round(b.left));
      byTop[t] = (byTop[t] || 0) + 1;
      byLeft[l] = (byLeft[l] || 0) + 1;
    }
    const rowCounts = Object.values(byTop);
    return {
      count: cells.length,
      rows: rowCounts.length,
      cols: Object.keys(byLeft).length,
      maxPerRow: rowCounts.length ? Math.max(...rowCounts) : 0,
    };
  });
}

async function scrollSliderExtent(page: Page, axis: 'x' | 'y'): Promise<{present: number; ratio: number}> {
  return page.evaluate((ax) => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
    if (!root) return {present: -1, ratio: -1};
    const svgs = Array.from(root.querySelectorAll(
      `.d4-layout-center > svg[type="range-slider"][name="${ax}-slider"]`));
    if (svgs.length === 0) return {present: 0, ratio: -1};
    const bar = svgs[0].querySelector('[name="pan-handle"]');
    if (!bar) return {present: svgs.length, ratio: -1};
    const tb = svgs[0].getBoundingClientRect();
    const bb = bar.getBoundingClientRect();
    const track = ax === 'x' ? tb.width : tb.height;
    const span = ax === 'x' ? bb.width : bb.height;
    return {present: svgs.length, ratio: track > 0 ? span / track : -1};
  }, axis);
}

async function categoryLabels(page: Page): Promise<{x: string[]; y: string[]; xAngles: number[]; yAngles: number[]}> {
  return page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
    const read = (cls: string) => {
      const nodes = root ? Array.from(root.querySelectorAll('.' + cls)) : [];
      return nodes.filter((n) => n.tagName.toLowerCase() === 'text');
    };
    const angleOf = (el: Element) => {
      const m = /rotate\(\s*(-?[\d.]+)/.exec(el.getAttribute('transform') ?? '');
      return m ? parseFloat(m[1]) : NaN;
    };
    const xs = read('d4-trellis-plot-cat-item-horz');
    const ys = read('d4-trellis-plot-cat-item-vert');
    return {
      x: xs.map((n) => (n.textContent ?? '').trim()),
      y: ys.map((n) => (n.textContent ?? '').trim()),
      xAngles: xs.map(angleOf),
      yAngles: ys.map(angleOf),
    };
  });
}

async function currentCellIndex(page: Page): Promise<number> {
  return page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
    if (!root) return -1;
    return Array.from(root.querySelectorAll('.d4-trellis-plot-cell'))
      .findIndex((c) => c.classList.contains('d4-trellis-cell-current'));
  });
}

async function focusChartsGrid(page: Page): Promise<void> {
  await page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
    const grid = (root?.querySelector('.d4-trellis-plot-charts-grid') ??
      root?.querySelector('div[tabindex="-1"]')) as HTMLElement | null;
    grid?.focus();
  });
}

async function comboFromLastCellChange(page: Page): Promise<Record<string, string>> {
  return page.evaluate(() => {
    const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
    const df = grok.shell.tv.dataFrame;
    const cols: string[] = [...tp.props.xColumnNames, ...tp.props.yColumnNames];
    const cats: string[][] = cols.map((c) => Array.from(df.col(c).categories));
    const out: Record<string, string> = {};
    const seen = new Set<any>();
    const walk = (node: any, depth: number) => {
      if (node == null || depth > 6) return;
      if (typeof node === 'string') {
        const owners = cols.filter((_c, i) => cats[i].indexOf(node) >= 0);
        if (owners.length === 1) out[owners[0]] = node;
        return;
      }
      if (typeof node !== 'object' || seen.has(node)) return;
      seen.add(node);
      for (const k of Object.keys(node)) walk(node[k], depth + 1);
    };
    const list = ((window as any).__escCc ?? []) as any[];
    walk(list[list.length - 1], 0);
    return out;
  });
}

async function armDfEvent(page: Page, event: 'onRowsFiltered' | 'onSelectionChanged' | 'onFilterChanged'): Promise<void> {
  await page.evaluate((ev) => {
    const df = grok.shell.tv.dataFrame;
    const w = window as any;
    w.__dfEvt = {fired: false, sub: null};
    try { w.__dfEvt.sub = df[ev].subscribe(() => { w.__dfEvt.fired = true; }); } catch (_) {}
  }, event);
}

async function awaitArmedDfEvent(page: Page, capMs: number): Promise<void> {
  await page.waitForFunction(() => (window as any).__dfEvt?.fired === true, null, {timeout: capMs}).catch(() => {});
  await page.evaluate(() => { try { (window as any).__dfEvt?.sub?.unsubscribe(); } catch (_) {} });
}

test('Trellis plot tests', async ({page}) => {

  test.setTimeout(900_000);
  page.setDefaultTimeout(120_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath});
  await v.installEventWaits(page);

  const setup = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    return {rowCount: df.rowCount, sex: df.col('SEX').categories.length, race: df.col('RACE').categories.length};
  });
  const fullRowCount = setup.rowCount;
  expect(setup).toEqual({rowCount: 5850, sex: 2, race: 4});

  const canonicalCellCount = axisViewportCount(setup.sex, false) * axisViewportCount(setup.race, false);

  await v.addViewerByIcon(page, 'trellis-plot', 'Trellis-plot', 15000);

  const restoreCanonical = async () => {
    await page.evaluate(async () => {
      try {
        for (const vw of Array.from(grok.shell.tv.viewers) as any[])
          if (vw.type !== 'Grid' && vw.type !== 'Trellis plot') vw.close();
        const trellises = Array.from(grok.shell.tv.viewers).filter((v: any) => v.type === 'Trellis plot') as any[];
        for (let i = 1; i < trellises.length; i++) trellises[i].close();
        let tp = trellises[0];
        if (!tp) tp = grok.shell.tv.addViewer('Trellis plot');
        tp.props.globalScale = false;
        tp.props.showXAxes = 'Auto';
        tp.props.showYAxes = 'Auto';
        tp.props.onClick = 'None';
        tp.props.viewerType = 'Scatter plot';
        tp.props.xColumnNames = ['SEX'];
        tp.props.yColumnNames = ['RACE'];
        await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1500);
      } catch (_) {  }
    });
  };

  await softStep('Inner viewer types', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;

      const settle = (cap: number) => (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, cap);
      const w = window as any;
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      await settle(1500);
      const r: any[] = [];
      const cellsHaveCanvas = () => {
        const cells = root.querySelectorAll('.d4-trellis-plot-cell');
        let withCanvas = 0;
        for (const c of Array.from(cells)) if (c.querySelector('canvas')) withCanvas++;
        return withCanvas;
      };

      const uiSwitch = async (icon: string) => {
        const vs = root.querySelector('[name="viewer selector"]') as HTMLElement;
        // open the combo and settle on its drop-down attaching rather than a fixed 600ms
        vs.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
        await w.__poll(() => document.querySelector('.d4-combo-drop-down'), (e: Element | null) => !!e, 1000, 40);
        const item = document.querySelector(`.d4-combo-drop-down [name="${icon}"]`);
        // switching the inner type is a heavy rebuild — race its repaint burst going quiet
        await w.__settled('viewer:Trellis plot.onViewerRendered',
          () => (item?.closest('.d4-list-item') as HTMLElement | null)?.click(), 1600);
        await settle(1600);
      };
      const types: Array<[string, string, any]> = [
        ['Scatter plot', 'icon-scatter-plot', {xColumnName: 'WEIGHT', yColumnName: 'HEIGHT', colorColumnName: 'RACE'}],
        ['Bar chart', 'icon-bar-chart', {splitColumnName: 'RACE', valueColumnName: 'AGE', valueAggrType: 'avg'}],
        ['Histogram', 'icon-histogram', {valueColumnName: 'AGE', splitColumnName: 'RACE'}],
        ['Line chart', 'icon-line-chart', {xColumnName: 'STARTED', yColumnName: 'AGE', splitColumnName: 'RACE'}],
        ['Box plot', 'icon-box-plot', {categoryColumnName: 'SEX', valueColumnName: 'AGE'}],
        ['Pie chart', 'icon-pie-chart', {categoryColumnName: 'RACE'}],
        ['Density plot', 'icon-density-plot', {xColumnName: 'WEIGHT', yColumnName: 'HEIGHT'}],
        ['Summary', 'icon-summary', {visualization: 'bars'}],
        ['Sparklines', 'icon-sparklines', {sparklineType: 'Bar Chart'}],
        ['PC Plot', 'icon-pc-plot', {colorColumnName: 'SEX'}],

        ['Heatmap', 'icon-heat-map', {columnNames: ['AGE', 'HEIGHT', 'WEIGHT']}],
      ];
      for (const [, icon, look] of types) {
        await uiSwitch(icon);
        try { tp.setOptions({innerViewerLook: look}); } catch {}
        await settle(900);
        r.push({type: tp.props.viewerType, cellsWithCanvas: cellsHaveCanvas()});
      }
      return r;
    });

    const expectedTypes = ['Scatter plot', 'Bar chart', 'Histogram', 'Line chart', 'Box plot',
      'Pie chart', 'Density plot', 'Summary', 'Sparklines', 'PC Plot', 'Heatmap'];
    expect(result.map((x: any) => x.type)).toEqual(expectedTypes);
    for (const x of result) expect(x.cellsWithCanvas).toBeGreaterThan(0);

    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      await (window as any).__settled('viewer:Trellis plot.onViewerRendered',
        () => tp.setOptions({innerViewerLook: {columnNames: ['AGE', 'HEIGHT', 'WEIGHT']}}), 1400);
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1400);
    });
    const pointsBefore = await v.trellisCellHashes(page, [0, 1]);
    const appliedAfter = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      await (window as any).__settled('viewer:Trellis plot.onViewerRendered',
        () => tp.setOptions({innerViewerLook: {columnNames: ['AGE']}}), 1400);
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1400);
      return tp.getOptions().look.innerViewerLook.columnNames;
    });
    const pointsAfter = await v.trellisCellHashes(page, [0, 1]);
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1400);
    });
    expect(appliedAfter).toEqual(['AGE']);

    expect(pointsBefore[0]).not.toBeNull();
    expect(pointsBefore[1]).not.toBeNull();
    expect(pointsAfter[0]).not.toBeNull();
    expect(pointsAfter[1]).not.toBeNull();
    expect(pointsAfter[0]).not.toBe(pointsBefore[0]);
    expect(pointsAfter[1]).not.toBe(pointsBefore[1]);
  });

  await softStep('Global scale', async () => {
   try {

    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      tp.setOptions({innerViewerLook: {xColumnName: 'WEIGHT', yColumnName: 'HEIGHT'}});
      tp.props.globalScale = false;
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 3000);
    });
    const idxA = await cellIndexFor(page, 'F', 'Caucasian');
    const idxB = await cellIndexFor(page, 'M', 'Caucasian');
    const cellHashes = async () => {
      const h = await v.trellisCellHashes(page, [idxA, idxB]);
      return {a: h[0], b: h[1]};
    };

    const idleBefore = await cellHashes();
    // the idle driven-guard: with nothing touched the cells must NOT repaint on their own, so a
    // subsequent delta cannot be booked by ambient repaint. There is no event to await for the
    // absence of a repaint — __quiet returns 0 when nothing fires, which is exactly the guard.
    const idleRenders = await page.evaluate(() =>
      (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 400, 1500));
    const idleAfter = await cellHashes();
    expect(idleRenders).toBe(0);
    expect(idleBefore.a).not.toBeNull();
    expect(idleBefore.b).not.toBeNull();
    expect(idleAfter.a).toBe(idleBefore.a);
    expect(idleAfter.b).toBe(idleBefore.b);

    const flip = async (value: boolean) => {
      await page.evaluate(async (val) => {
        const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
        await (window as any).__settled('viewer:Trellis plot.onViewerRendered',
          () => { tp.props.globalScale = val; }, 1600);
        await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1600);
      }, value);
      return cellHashes();
    };

    const on = await flip(true);
    expect(on.a).not.toBeNull();
    expect(on.b).not.toBeNull();
    expect(on.a).not.toBe(idleAfter.a);
    expect(on.b).not.toBe(idleAfter.b);

    const off = await flip(false);
    expect(off.a).not.toBeNull();
    expect(off.b).not.toBeNull();
    expect(off.a).not.toBe(on.a);
    expect(off.b).not.toBe(on.b);

    const onAgain = await flip(true);
    expect(onAgain.a).not.toBeNull();
    expect(onAgain.b).not.toBeNull();
    expect(onAgain.a).not.toBe(off.a);
    expect(onAgain.b).not.toBe(off.b);
   } finally {

    await restoreCanonical();
   }
  });

  await softStep('Axes visibility', async () => {
   try {

    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;

      const settle = (cap: number) => (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, cap);
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      tp.props.viewerType = 'Scatter plot';
      tp.props.globalScale = true;
      await settle(1500);

      const innerSliderCount = () => root.querySelectorAll('.d4-range-selector > svg[type="range-slider"]').length;

      const axisSliders = (ax: string) =>
        Array.from(root.querySelectorAll(`.d4-range-selector > svg[type="range-slider"][name="${ax}-slider"]`))
          .filter((el) => {
            const b = el.getBoundingClientRect();
            return b.width > 0 && b.height > 0;
          }).length;
      tp.props.showRangeSliders = true;
      tp.props.showYAxes = 'Always';
      await settle(1200);
      const r: any[] = [];
      const xByMode: Record<string, number> = {};
      for (const val of ['Always', 'Never', 'Auto']) {
        tp.props.showXAxes = val;
        await settle(1200);
        r.push(tp.props.showXAxes);
        xByMode[val] = axisSliders('x');
      }
      tp.props.showXAxes = 'Always';
      await settle(1200);
      const yByMode: Record<string, number> = {};
      for (const val of ['Always', 'Never', 'Auto']) {
        tp.props.showYAxes = val;
        await settle(1200);
        r.push(tp.props.showYAxes);
        yByMode[val] = axisSliders('y');
      }
      tp.props.showXAxes = 'Always';
      tp.props.showYAxes = 'Always';
      tp.props.showRangeSliders = false;
      await settle(1200);
      const slidersOff = innerSliderCount();
      tp.props.showRangeSliders = true;
      await settle(1200);
      const slidersOn = innerSliderCount();
      return {modes: r, xByMode, yByMode, slidersOff, slidersOn};
    });
    expect(result.modes).toEqual(['Always', 'Never', 'Auto', 'Always', 'Never', 'Auto']);

    expect(result.xByMode.Always).toBeGreaterThan(0);
    expect(result.xByMode.Never).toBe(0);
    expect(result.yByMode.Always).toBeGreaterThan(0);
    expect(result.yByMode.Never).toBe(0);

    expect(result.xByMode.Auto).toBe(result.xByMode.Always);
    expect(result.yByMode.Auto).toBe(result.yByMode.Always);

    expect(result.slidersOff).toBe(0);
    expect(result.slidersOn).toBeGreaterThan(0);
   } finally {

    await page.evaluate(async () => {
      try {
        const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
        tp.props.showRangeSliders = true;
        await new Promise((res) => setTimeout(res, 400));
      } catch (_) {  }
    });
    await restoreCanonical();
   }
  });

  await softStep('Range sliders with global scale', async () => {
   try {

    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      tp.props.viewerType = 'Scatter plot';
      tp.setOptions({innerViewerLook: {xColumnName: 'WEIGHT', yColumnName: 'HEIGHT'}});
      tp.props.globalScale = true;
      tp.props.showRangeSliders = true;
      tp.props.showXAxes = 'Always';
      tp.props.showYAxes = 'Always';
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 2200);
    });

    const sliderBox = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const innerX = root.querySelector('.d4-range-selector > svg[type="range-slider"][name="x-slider"]') as SVGElement | null;
      if (!innerX) return null;
      const wrap = innerX.closest('.d4-range-selector') as HTMLElement;
      const wb = wrap.getBoundingClientRect();
      wrap.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: wb.left + wb.width / 2, clientY: wb.top + wb.height / 2}));
      await new Promise((r) => setTimeout(r, 400));
      const b = innerX.getBoundingClientRect();
      return {x: b.x, y: b.y, w: b.width, h: b.height};
    });
    expect(sliderBox).not.toBeNull();
    expect(sliderBox!.w).toBeGreaterThan(0);

    const idxA = await cellIndexFor(page, 'F', 'Caucasian');
    const idxB = await cellIndexFor(page, 'M', 'Caucasian');
    const beforeH = await v.trellisCellHashes(page, [idxA, idxB]);
    const before = {a: beforeH[0], b: beforeH[1]};
    await page.mouse.move(sliderBox!.x + sliderBox!.w - 4, sliderBox!.y + sliderBox!.h / 2);
    await page.mouse.down();
    await page.mouse.move(sliderBox!.x + sliderBox!.w * 0.45, sliderBox!.y + sliderBox!.h / 2, {steps: 12});
    await page.mouse.up();
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    const afterH = await v.trellisCellHashes(page, [idxA, idxB]);
    const after = {a: afterH[0], b: afterH[1]};

    expect(before.a).not.toBeNull();
    expect(before.b).not.toBeNull();
    expect(after.a).not.toBeNull();
    expect(after.b).not.toBeNull();
    expect(after.a).not.toBe(before.a);
    expect(after.b).not.toBe(before.b);

    await page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell').first().click({button: 'right', position: {x: 6, y: 6}});
    await clickTopLevelMenuItem(page, 'Reset Inner Range Sliders');
    await v.waitForViewerRendered(page, 'Trellis plot', 1200);
    const resetH = await v.trellisCellHashes(page, [idxA, idxB]);
    const reset = {restoredA: resetH[0] === before.a, restoredB: resetH[1] === before.b};
    expect(reset.restoredA).toBe(true);
    expect(reset.restoredB).toBe(true);

    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      await (window as any).__settled('viewer:Trellis plot.onViewerRendered',
        () => { tp.props.showYAxes = 'Always'; }, 1200);
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1200);
    });
    const ySliderBox = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const innerY = root.querySelector('.d4-range-selector > svg[type="range-slider"][name="y-slider"]') as SVGElement | null;
      if (!innerY) return null;
      const wrap = innerY.closest('.d4-range-selector') as HTMLElement;
      const wb = wrap.getBoundingClientRect();
      wrap.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: wb.left + wb.width / 2, clientY: wb.top + wb.height / 2}));
      await new Promise((r) => setTimeout(r, 400));
      const b = innerY.getBoundingClientRect();
      return {x: b.x, y: b.y, w: b.width, h: b.height};
    });
    expect(ySliderBox).not.toBeNull();
    expect(ySliderBox!.h).toBeGreaterThan(0);
    const yBeforeH = await v.trellisCellHashes(page, [idxA, idxB]);
    const yBefore = {a: yBeforeH[0], b: yBeforeH[1]};
    await page.mouse.move(ySliderBox!.x + ySliderBox!.w / 2, ySliderBox!.y + ySliderBox!.h - 4);
    await page.mouse.down();
    await page.mouse.move(ySliderBox!.x + ySliderBox!.w / 2, ySliderBox!.y + ySliderBox!.h * 0.45, {steps: 12});
    await page.mouse.up();
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    const yAfterH = await v.trellisCellHashes(page, [idxA, idxB]);
    const yAfter = {a: yAfterH[0], b: yAfterH[1]};
    expect(yBefore.a).not.toBeNull();
    expect(yBefore.b).not.toBeNull();
    expect(yAfter.a).not.toBeNull();
    expect(yAfter.b).not.toBeNull();
    expect(yAfter.a).not.toBe(yBefore.a);
    expect(yAfter.b).not.toBe(yBefore.b);
   } finally {

    await restoreCanonical();
   }
  });

  await softStep('Gridlines', async () => {

    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const settle = (cap: number) => (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 250, cap);
      const hasGrid = () => !!root.querySelector('.d4-trellis-plot-charts-grid');
      const r: {mode: string; prop: string; grid: boolean}[] = [];
      for (const val of ['always', 'never', 'auto']) {
        tp.props.showGridlines = val;
        await settle(600);
        r.push({mode: val, prop: tp.props.showGridlines, grid: hasGrid()});
      }
      tp.props.viewerType = 'Bar chart';
      await settle(1200);
      const autoNonScatter = hasGrid();
      tp.props.viewerType = 'Scatter plot';
      await settle(1200);
      return {r, autoNonScatter, autoScatter: hasGrid()};
    });
    expect(result.r.map((x) => x.prop)).toEqual(['always', 'never', 'auto']);
    expect(result.r.map((x) => x.grid)).toEqual([true, false, true]);
    expect(result.autoNonScatter).toBe(false);
    expect(result.autoScatter).toBe(true);
  });

  await softStep('Tiles mode', async () => {
   try {

    const setup = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const settle = (cap: number) => (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, cap);
      tp.props.globalScale = false;
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['RACE'];
      tp.props.yColumnNames = [];
      await settle(1500);
      tp.props.useTiledView = true;
      await settle(1800);
      return {raceCats: grok.shell.tv.dataFrame.col('RACE').categories.length,
        tilesPerRow: tp.props.tilesPerRow};
    });

    expect(setup.raceCats).toBeGreaterThanOrEqual(3);
    expect(setup.raceCats).toBeLessThanOrEqual(6);

    const tiledGeometry = (n: number, tilesPerRow: number) => {
      const xRaw = Math.min(tilesPerRow, n);
      const w = axisViewportCount(xRaw, true);
      const h = axisViewportCount(Math.ceil(n / xRaw), true);
      return {count: w * h, rows: h, maxPerRow: w};
    };
    const tilesOn = await gridGeometry(page);
    expect({count: tilesOn.count, rows: tilesOn.rows, maxPerRow: tilesOn.maxPerRow})
      .toEqual(tiledGeometry(setup.raceCats, setup.tilesPerRow));

    const setTiles = async (n: number) => {
      await page.evaluate(async (val) => {
        const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
        tp.props.tilesPerRow = val;
        await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1600);
      }, n);
      return gridGeometry(page);
    };

    const atTwo = await setTiles(2);
    expect({count: atTwo.count, rows: atTwo.rows, maxPerRow: atTwo.maxPerRow})
      .toEqual(tiledGeometry(setup.raceCats, 2));
    expect(atTwo.maxPerRow).toBeLessThanOrEqual(2);
    expect(atTwo.rows).toBeGreaterThan(1);

    const atSix = await setTiles(6);
    expect({count: atSix.count, rows: atSix.rows, maxPerRow: atSix.maxPerRow})
      .toEqual(tiledGeometry(setup.raceCats, 6));
    expect(atSix.rows).toBe(1);

    const backToTwo = await setTiles(2);
    expect(backToTwo.rows).toBeGreaterThan(1);
    const untiled = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.useTiledView = false;
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1800);
      return tp.props.useTiledView;
    });
    const offGeometry = await gridGeometry(page);
    expect(untiled).toBe(false);

    expect(offGeometry.count).toBe(axisViewportCount(setup.raceCats, true));
    expect({rows: offGeometry.rows, maxPerRow: offGeometry.maxPerRow})
      .not.toEqual({rows: backToTwo.rows, maxPerRow: backToTwo.maxPerRow});
   } finally {
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      try {
        tp.props.useTiledView = true;
        tp.props.tilesPerRow = 4;
        await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 600);
      } catch (_) {  }
    });
    await restoreCanonical();
   }
  });

  await softStep('Category management', async () => {

    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const cells = () => root.querySelectorAll('.d4-trellis-plot-cell').length;
      const r: any[] = [];

      const settle = (cap: number) => (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, cap);
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      await settle(1500);
      r.push({x: [...tp.props.xColumnNames], cells: cells()});

      tp.props.xColumnNames = ['SEX', 'DIS_POP'];
      await settle(1500);
      r.push({x: [...tp.props.xColumnNames], cells: cells()});

      tp.props.xColumnNames = ['SEX'];
      await settle(1500);
      r.push({x: [...tp.props.xColumnNames], cells: cells()});

      return r;
    });

    const cardinalities = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {sex: df.col('SEX').categories.length, race: df.col('RACE').categories.length,
        disPop: df.col('DIS_POP').categories.length};
    });
    const baseCells = axisViewportCount(cardinalities.sex, false) *
      axisViewportCount(cardinalities.race, false);
    const grownCells = axisViewportCount(cardinalities.sex * cardinalities.disPop, false) *
      axisViewportCount(cardinalities.race, false);
    expect(grownCells).toBeGreaterThan(baseCells);
    expect(result[0].x).toEqual(['SEX']);
    expect(result[0].cells).toBe(baseCells);
    expect(result[1].x).toEqual(['SEX', 'DIS_POP']);
    expect(result[1].cells).toBe(grownCells);
    expect(result[1].cells).toBeGreaterThan(result[0].cells);
    expect(result[2].x).toEqual(['SEX']);
    expect(result[2].cells).toBe(result[0].cells);

    const labelsOn = await categoryLabels(page);
    expect(labelsOn.x.length).toBeGreaterThan(0);
    expect(labelsOn.y.length).toBeGreaterThan(0);

    const setLabels = async (x: boolean, y: boolean) => {
      await page.evaluate(async (v) => {
        const tp = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Trellis plot') as any;
        tp.props.showXLabels = v.x;
        tp.props.showYLabels = v.y;
        await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1600);
      }, {x, y});
      return categoryLabels(page);
    };

    const labelsOff = await setLabels(false, false);
    expect(labelsOff.x.length).toBe(0);
    expect(labelsOff.y.length).toBe(0);

    const labelsBack = await setLabels(true, true);
    expect(labelsBack.x.length).toBe(labelsOn.x.length);
    expect(labelsBack.y.length).toBe(labelsOn.y.length);
  });

  await softStep('Pack categories', async () => {

    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const df = grok.shell.tv.dataFrame;
      const cells = () => root.querySelectorAll('.d4-trellis-plot-cell').length;
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      const w = window as any;
      const settle = (cap: number) => w.__quiet('viewer:Trellis plot.onViewerRendered', 300, cap);
      tp.props.packCategories = true;
      await settle(1600);
      const baseCells = cells();

      grok.shell.tv.getFiltersGroup();
      await w.__quiet('viewer:Trellis plot.onViewerRendered', 300, 1500);
      const fg = grok.shell.tv.getFiltersGroup();
      const cats = df.col('RACE').categories;
      // filter out one category — settle on the rows-filtered event, then the trellis repaint
      await w.__settled('df.onRowsFiltered',
        () => fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: cats.filter((c: string) => c !== 'Asian')}), 1500);
      await settle(1500);
      const packedOnCells = cells();

      tp.props.packCategories = false;
      await settle(1200);
      const packedOffCells = cells();

      tp.props.packCategories = true;
      await settle(1200);
      const packedOnAgainCells = cells();

      await w.__settled('df.onRowsFiltered',
        () => fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: cats}), 800);
      return {baseCells, packedOnCells, packedOffCells, packedOnAgainCells};
    });
    expect(result.packedOnCells).toBeLessThan(result.baseCells);
    expect(result.packedOffCells).toBe(result.baseCells);
    expect(result.packedOnAgainCells).toBe(result.packedOnCells);
  });

  await softStep('On Click functionality', async () => {
    const cellLocator = page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell');
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      df.filter.setAll(true); df.selection.setAll(false); df.rows.requestFilter();
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1500);
    });
    await expect(cellLocator).toHaveCount(canonicalCellCount);

    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      await (window as any).__settled('viewer:Trellis plot.onPropertyValueChanged',
        () => { tp.props.onClick = 'Select'; }, 600);
    });
    const selExpected = await comboRowCount(page, 'SEX', 'F', 'RACE', 'Caucasian');
    let idx = await cellIndexFor(page, 'F', 'Caucasian');
    await armDfEvent(page, 'onSelectionChanged');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await awaitArmedDfEvent(page, 900);
    const selAfterClick = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(selAfterClick).toBe(selExpected);

    const selAfterTypeChange = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const settle = (cap: number) => (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, cap);
      tp.props.viewerType = 'Bar chart';
      tp.setOptions({innerViewerLook: {splitColumnName: 'RACE', valueColumnName: 'AGE'}});
      await settle(1200);
      const sel = grok.shell.tv.dataFrame.selection.trueCount;
      tp.props.viewerType = 'Scatter plot';
      await settle(1000);
      return sel;
    });
    expect(selAfterTypeChange).toBe(selExpected);

    const mBlackExpected = await comboRowCount(page, 'SEX', 'M', 'RACE', 'Black');
    idx = await cellIndexFor(page, 'M', 'Black');
    await armDfEvent(page, 'onSelectionChanged');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await awaitArmedDfEvent(page, 900);
    const selAfterOther = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(selAfterOther).toBe(mBlackExpected);
    expect(selAfterOther).not.toBe(selExpected);

    const selAfterAxis = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const settle = (cap: number) => (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, cap);
      tp.props.yColumnNames = ['SEVERITY'];
      await settle(1200);
      const sel = grok.shell.tv.dataFrame.selection.trueCount;
      tp.props.yColumnNames = ['RACE'];
      await settle(1000);
      return sel;
    });
    expect(selAfterAxis).toBe(mBlackExpected);

    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true); df.selection.setAll(false); df.rows.requestFilter();
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      await (window as any).__settled('viewer:Trellis plot.onPropertyValueChanged',
        () => { tp.props.onClick = 'Filter'; }, 700);
    });
    const filterExpected = await comboRowCount(page, 'SEX', 'F', 'RACE', 'Caucasian');
    idx = await cellIndexFor(page, 'F', 'Caucasian');
    await armDfEvent(page, 'onRowsFiltered');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await awaitArmedDfEvent(page, 900);
    const filterAfterClick = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const filters: string[] = [];
      for (const f of df.rows.filters) filters.push(f);
      return {count: df.filter.trueCount, filters};
    });
    expect(filterAfterClick.count).toBe(filterExpected);
    expect(filterAfterClick.count).toBeLessThan(fullRowCount);
    expect(filterAfterClick.filters).toContain('SEX: F');
    expect(filterAfterClick.filters).toContain('RACE: Caucasian');

    const filterAfterType = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const settle = (cap: number) => (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, cap);
      tp.props.viewerType = 'Pie chart';
      tp.setOptions({innerViewerLook: {categoryColumnName: 'RACE'}});
      await settle(1200);
      const c = grok.shell.tv.dataFrame.filter.trueCount;
      tp.props.viewerType = 'Scatter plot';
      await settle(1000);
      return c;
    });
    expect(filterAfterType).toBe(filterExpected);

    const mBlackFilter = await comboRowCount(page, 'SEX', 'M', 'RACE', 'Black');
    idx = await cellIndexFor(page, 'M', 'Black');
    await armDfEvent(page, 'onRowsFiltered');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await awaitArmedDfEvent(page, 900);
    const filterAfterOther = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const filters: string[] = [];
      for (const f of df.rows.filters) filters.push(f);
      return {count: df.filter.trueCount, filters};
    });
    expect(filterAfterOther.count).toBe(mBlackFilter);
    expect(filterAfterOther.filters).toContain('SEX: M');
    expect(filterAfterOther.filters).toContain('RACE: Black');

    const filterAfterAxis = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const settle = (cap: number) => (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, cap);
      tp.props.xColumnNames = ['CONTROL'];
      await settle(1500);
      const c = grok.shell.tv.dataFrame.filter.trueCount;
      tp.props.xColumnNames = ['SEX'];
      await settle(1500);
      return c;
    });
    expect(filterAfterAxis).toBe(fullRowCount);

    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true); df.selection.setAll(false); df.rows.requestFilter();
      await new Promise((r) => setTimeout(r, 500));
    });
    idx = await cellIndexFor(page, 'F', 'Caucasian');
    await armDfEvent(page, 'onRowsFiltered');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await awaitArmedDfEvent(page, 700);
    const beforeEsc = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const filters: string[] = [];
      for (const f of df.rows.filters) filters.push(f);
      return {count: df.filter.trueCount, filters};
    });
    expect(beforeEsc.count).toBeLessThan(fullRowCount);
    expect(beforeEsc.filters.length).toBeGreaterThan(0);
    await focusChartsGrid(page);
    await armDfEvent(page, 'onRowsFiltered');
    await page.keyboard.press('Escape');
    await awaitArmedDfEvent(page, 900);
    const afterEsc = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const filters: string[] = [];
      for (const f of df.rows.filters) filters.push(f);
      return {count: df.filter.trueCount, filters};
    });
    expect(afterEsc.count).toBe(fullRowCount);
    expect(afterEsc.filters).toEqual([]);

    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.onClick = 'Select';
      (window as any).__escCc = [];
      (window as any).__escCcSub = tp.onEvent('d4-trellis-plot-current-cell-changed').subscribe((a: any) => {
        const mc = (a && a.args && a.args.matchCondition) ? a.args.matchCondition : (a && a.matchCondition ? a.matchCondition : a);
        (window as any).__escCc.push(mc);
      });
      await (window as any).__settled('df.onRowsFiltered',
        () => { df.filter.setAll(true); df.selection.setAll(false); df.rows.requestFilter(); }, 700);
    });
    idx = await cellIndexFor(page, 'F', 'Caucasian');
    await armDfEvent(page, 'onSelectionChanged');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await awaitArmedDfEvent(page, 900);

    expect(await currentCellIndex(page)).toBe(idx);
    const selComboRows = await comboRowCount(page, 'SEX', 'F', 'RACE', 'Caucasian');
    expect(selComboRows).toBeGreaterThan(0);
    const selBeforeEsc = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(selBeforeEsc).toBe(selComboRows);
    await focusChartsGrid(page);
    await armDfEvent(page, 'onSelectionChanged');
    await page.keyboard.press('Escape');
    await awaitArmedDfEvent(page, 900);
    const selAfterEsc = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(selAfterEsc).toBe(0);

    expect(await comboFromLastCellChange(page)).toEqual({SEX: 'F', RACE: 'Caucasian'});

    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      (window as any).__escCcSub?.unsubscribe?.();
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.onClick = 'Filter';
      await (window as any).__settled('df.onRowsFiltered',
        () => { df.filter.setAll(true); df.selection.setAll(false); df.rows.requestFilter(); }, 700);
    });

    const panelOnly = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true); df.rows.requestFilter();
      const ageCol = df.col('AGE');
      (window as any).__panelSub = df.onRowsFiltering.subscribe(() => {
        for (let i = 0; i < df.rowCount; i++) if (ageCol.get(i) < 40) df.filter.set(i, false, false);
      });
      await (window as any).__settled('df.onRowsFiltered', () => df.rows.requestFilter(), 800);
      return df.filter.trueCount;
    });
    expect(panelOnly).toBeLessThan(fullRowCount);
    const mAsianAge = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const s = df.col('SEX'), r = df.col('RACE'), a = df.col('AGE');
      let n = 0;
      for (let i = 0; i < df.rowCount; i++) if (s.get(i) === 'M' && r.get(i) === 'Asian' && a.get(i) >= 40) n++;
      return n;
    });
    idx = await cellIndexFor(page, 'M', 'Asian');
    await armDfEvent(page, 'onRowsFiltered');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await awaitArmedDfEvent(page, 900);
    const bothActive = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(bothActive).toBe(mAsianAge);
    expect(bothActive).toBeLessThan(panelOnly);

    const noneResult = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      (window as any).__panelSub?.unsubscribe?.();
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      await (window as any).__settled('df.onRowsFiltered',
        () => { df.filter.setAll(true); df.selection.setAll(false); df.rows.requestFilter(); }, 700);
      tp.props.onClick = 'None';
      return {filterBefore: df.filter.trueCount, selBefore: df.selection.trueCount};
    });
    idx = await cellIndexFor(page, 'F', 'Caucasian');
    // On Click = None: the click must change nothing. Arm both df channels, click, and confirm
    // neither fired within the window — a stronger negative than a blind sleep-then-read.
    const noneEvents = await v.armEvent(page, 'df.onSelectionChanged', 900);
    const noneFilterEvt = await v.armEvent(page, 'df.onRowsFiltered', 900);
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    const selEvt = await noneEvents();
    const filtEvt = await noneFilterEvt();
    const noneAfter = await page.evaluate(() => ({
      filterAfter: grok.shell.tv.dataFrame.filter.trueCount,
      selAfter: grok.shell.tv.dataFrame.selection.trueCount,
    }));
    expect(selEvt).toBeNull();
    expect(filtEvt).toBeNull();
    expect(noneAfter.filterAfter).toBe(noneResult.filterBefore);
    expect(noneAfter.selAfter).toBe(noneResult.selBefore);
  });

  await softStep('Selectors', async () => {
   try {

    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;

      const settle = (cap: number) => (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, cap);
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      tp.props.showControlPanel = true;
      await settle(1500);
      const vsVisible = () => {
        const el = root.querySelector('[name="viewer selector"]') as HTMLElement | null;
        if (!el) return false;
        const b = el.getBoundingClientRect();
        return b.width > 0 && b.height > 0;
      };
      tp.props.showXSelectors = false;
      tp.props.showYSelectors = false;
      tp.props.showControlPanel = false;
      await settle(1200);
      const off = {props: {x: tp.props.showXSelectors, y: tp.props.showYSelectors, cp: tp.props.showControlPanel}, vs: vsVisible()};
      tp.props.showXSelectors = true;
      tp.props.showYSelectors = true;
      tp.props.showControlPanel = true;
      await settle(1200);
      const on = {props: {x: tp.props.showXSelectors, y: tp.props.showYSelectors, cp: tp.props.showControlPanel}, vs: vsVisible()};
      return {off, on};
    });
    expect(result.off.props).toEqual({x: false, y: false, cp: false});
    expect(result.on.props).toEqual({x: true, y: true, cp: true});
    expect(result.off.vs).toBe(false);
    expect(result.on.vs).toBe(true);

    const columnSelectorCount = () => page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
      if (!root) return -1;
      const reallyVisible = (el: Element) => {
        const b = el.getBoundingClientRect();
        if (b.width <= 0 || b.height <= 0) return false;
        for (let n: Element | null = el; n && n !== document.documentElement; n = n.parentElement) {
          const s = getComputedStyle(n);
          if (s.visibility === 'hidden' || s.visibility === 'collapse' || s.display === 'none') return false;
        }
        return true;
      };

      return Array.from(root.querySelectorAll('[name="div-column-combobox-"]'))
        .filter(reallyVisible).length;
    });
    const controlPanelVisible = () => page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
      const el = root?.querySelector('[name="viewer selector"]') as HTMLElement | null;
      if (!el) return false;
      const b = el.getBoundingClientRect();
      return b.width > 0 && b.height > 0;
    });
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1200);
      tp.props.autoLayout = true;
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1200);
    });

    const bothSelectors = await columnSelectorCount();
    expect(bothSelectors).toBe(2);

    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.showXSelectors = false;
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1500);
    });
    const xTurnedOff = await columnSelectorCount();
    expect(xTurnedOff).toBe(1);

    await page.setViewportSize({width: 500, height: 400});

    await page.waitForFunction(() => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
      const el = root?.querySelector('[name="viewer selector"]') as HTMLElement | null;
      if (!el) return false;
      const b = el.getBoundingClientRect();
      return !(b.width > 0 && b.height > 0);
    }, null, {timeout: 2000}).catch(() => {});
    const shrunkPanel = await controlPanelVisible();
    const shrunkSelectors = await columnSelectorCount();
    await page.setViewportSize({width: 1920, height: 1080});
    await page.waitForFunction(() => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
      const el = root?.querySelector('[name="viewer selector"]') as HTMLElement | null;
      if (!el) return false;
      const b = el.getBoundingClientRect();
      return b.width > 0 && b.height > 0;
    }, null, {timeout: 2000}).catch(() => {});
    const restoredPanel = await controlPanelVisible();
    const afterCycle = await columnSelectorCount();

    expect(shrunkPanel).toBe(false);
    expect(shrunkSelectors).toBe(0);
    expect(restoredPanel).toBe(true);

    expect(afterCycle).toBe(xTurnedOff);

    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.showXSelectors = true;
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1200);
    });
    expect(await columnSelectorCount()).toBe(bothSelectors);
   } finally {

    await page.setViewportSize({width: 1920, height: 1080});
    await v.waitForViewerRendered(page, 'Trellis plot', 1000);
    await page.evaluate(async () => {
      try {
        const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
        tp.props.showXSelectors = true;
        tp.props.showYSelectors = true;
        tp.props.showControlPanel = true;
        await new Promise((r) => setTimeout(r, 500));
      } catch (_) {  }
    });
   }
  });

  await softStep('Allow viewer full screen', async () => {
   try {

    await page.mouse.move(5, 5);
    await page.waitForTimeout(400); // technical: park pointer, no event to await before the baseline
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;

      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      tp.props.allowViewerFullScreen = true;
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1500);

      const rootEl = () => document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
      const cells = () => Array.from(rootEl()?.querySelectorAll('.d4-trellis-plot-cell') ?? []);
      const iconState = () => {
        const list = Array.from(rootEl()?.querySelectorAll('.d4-viewer-icon[name="icon-expand-arrows"]') ?? []) as HTMLElement[];
        if (list.length !== 1) return {count: list.length, parent: -1, top: -1, left: -1};
        const b = list[0].getBoundingClientRect();
        return {count: 1, parent: cells().indexOf(list[0].parentElement as Element),
          top: Math.round(b.top), left: Math.round(b.left)};
      };
      const w = window as any;
      const iconCount = () =>
        (rootEl()?.querySelectorAll('.d4-viewer-icon[name="icon-expand-arrows"]') ?? []).length;
      const pointer = async (idx: number, into: boolean) => {
        const cell = cells()[idx];
        if (!cell) return;
        const r = cell.getBoundingClientRect();
        const o = {bubbles: true, cancelable: true, clientX: r.left + r.width / 2, clientY: r.top + r.height / 2};
        for (const t of into ? ['mouseenter', 'mouseover', 'mousemove'] : ['mouseout', 'mouseleave'])
          cell.dispatchEvent(new MouseEvent(t, o));
        // settle on the hover icon appearing (into) or clearing (out), capped — not a blind 800ms
        await w.__poll(iconCount, (n: number) => (into ? n >= 1 : n === 0), 800, 40);
      };

      const all = cells();
      const firstIdx = all.findIndex((c) => !!c.querySelector('canvas'));
      const firstTop = firstIdx >= 0 ? Math.round(all[firstIdx].getBoundingClientRect().top) : -1;

      const secondIdx = all.findIndex((c, i) => i !== firstIdx && !!c.querySelector('canvas') &&
        Math.round(c.getBoundingClientRect().top) !== firstTop);
      const out: any = {cellsFound: false, firstIdx, secondIdx,
        beforeHover: iconState().count, onFirst: {count: -1, parent: -1, top: -1, left: -1},
        onSecond: {count: -1, parent: -1, top: -1, left: -1}, afterLeave: -1,
        modalOpened: false, modalTitle: null as string | null, iconWhenOff: true};
      if (firstIdx < 0 || secondIdx < 0) return out;
      out.cellsFound = true;

      await pointer(firstIdx, true);
      out.onFirst = iconState();

      await pointer(firstIdx, false);
      await pointer(secondIdx, true);
      out.onSecond = iconState();

      await pointer(secondIdx, false);
      out.afterLeave = iconState().count;

      await pointer(firstIdx, true);
      const icon = rootEl()?.querySelector('.d4-viewer-icon[name="icon-expand-arrows"]') as HTMLElement | null;
      if (icon) {
        const ib = icon.getBoundingClientRect();
        const io = {bubbles: true, cancelable: true, button: 0, clientX: ib.left + ib.width / 2, clientY: ib.top + ib.height / 2};
        icon.dispatchEvent(new MouseEvent('mousedown', io));
        icon.dispatchEvent(new MouseEvent('mouseup', io));
        icon.dispatchEvent(new MouseEvent('click', io));
        // settle on the full-screen dialog attaching, capped
        const dlg = await w.__poll(() => document.querySelector('.d4-dialog'),
          (e: Element | null) => !!e, 1200, 50) as Element | null;
        out.modalOpened = !!dlg;
        if (dlg) {
          out.modalTitle = (dlg.querySelector('.d4-dialog-title') as HTMLElement | null)?.textContent?.trim() ?? null;
          const cancel = (dlg.querySelector('[name="button-CANCEL"]') as HTMLElement | null)
            ?? (dlg.querySelector('.d4-dialog-header [name="icon-times"]') as HTMLElement | null);
          if (cancel) cancel.click(); else document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
          await w.__poll(() => document.querySelector('.d4-dialog'), (e: Element | null) => !e, 500, 40);
        }
      }
      await pointer(firstIdx, false);

      await w.__settled('viewer:Trellis plot.onPropertyValueChanged',
        () => { tp.props.allowViewerFullScreen = false; }, 800);
      await pointer(firstIdx, true);
      out.iconWhenOff = iconState().count > 0;
      await pointer(firstIdx, false);
      tp.props.allowViewerFullScreen = true;
      return out;
    });
    expect(result.cellsFound).toBe(true);

    expect(result.beforeHover).toBe(0);
    expect(result.onFirst.count).toBe(1);
    expect(result.onFirst.parent).toBe(result.firstIdx);
    expect(result.onSecond.count).toBe(1);
    expect(result.onSecond.parent).toBe(result.secondIdx);
    expect(result.onSecond.top === result.onFirst.top && result.onSecond.left === result.onFirst.left).toBe(false);
    expect(result.afterLeave).toBe(0);
    expect(result.modalOpened).toBe(true);

    expect(result.modalTitle).toContain('SEX');
    expect(result.iconWhenOff).toBe(false);
   } finally {

    await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
      for (const c of Array.from(root?.querySelectorAll('.d4-trellis-plot-cell') ?? [])) {
        const r = c.getBoundingClientRect();
        const o = {bubbles: true, cancelable: true, clientX: r.left + r.width / 2, clientY: r.top + r.height / 2};
        c.dispatchEvent(new MouseEvent('mouseout', o));
        c.dispatchEvent(new MouseEvent('mouseleave', o));
      }
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      if (tp) tp.props.allowViewerFullScreen = true;
      await new Promise((r) => setTimeout(r, 400));
    });
    await page.mouse.move(5, 5);

    await page.waitForTimeout(400);
   }
  });

  await softStep('Scrolling', async () => {
   try {

    const cardinality = (cols: string[]) => page.evaluate((cs) => {
      const df = grok.shell.tv.dataFrame;
      return cs.reduce((n: number, c: string) => n * df.col(c).categories.length, 1);
    }, cols);
    const configure = (x: string[], y: string[]) => page.evaluate(async (cols) => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.packCategories = false;
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = cols.x;
      tp.props.yColumnNames = cols.y;
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 400, 2500);
      return {x: [...tp.props.xColumnNames], y: [...tp.props.yColumnNames]};
    }, {x, y});

    const fitting = await configure(['SEX'], ['RACE']);
    expect(fitting.x).toEqual(['SEX']);
    const fitsRawX = await cardinality(['SEX']);
    const fitsRawY = await cardinality(['RACE']);
    const fitsGeom = await gridGeometry(page);
    expect(fitsGeom.cols).toBe(axisViewportCount(fitsRawX, false));
    expect(fitsGeom.cols).toBe(fitsRawX);
    expect(fitsGeom.rows).toBe(fitsRawY);
    const fitsX = await scrollSliderExtent(page, 'x');
    const fitsY = await scrollSliderExtent(page, 'y');
    expect(fitsX.present).toBeGreaterThan(0);
    expect(fitsY.present).toBeGreaterThan(0);
    expect(fitsX.ratio).toBeGreaterThan(0.5);
    expect(fitsY.ratio).toBeGreaterThan(0.5);

    const overX = await configure(['SEX', 'DIS_POP', 'RACE'], ['SEVERITY']);
    expect(overX.x).toEqual(['SEX', 'DIS_POP', 'RACE']);
    const overRawX = await cardinality(['SEX', 'DIS_POP', 'RACE']);
    const overGeomX = await gridGeometry(page);

    expect(overGeomX.cols).toBe(axisViewportCount(overRawX, false));
    expect(overGeomX.cols).toBeLessThan(overRawX);
    const scrolledX = await scrollSliderExtent(page, 'x');
    expect(scrolledX.present).toBeGreaterThan(0);

    expect(scrolledX.ratio).toBeGreaterThan(0);
    expect(scrolledX.ratio).toBeLessThan(fitsX.ratio);

    const overY = await configure(['SEX'], ['RACE', 'DIS_POP', 'SEVERITY']);
    expect(overY.y).toEqual(['RACE', 'DIS_POP', 'SEVERITY']);
    const overRawY = await cardinality(['RACE', 'DIS_POP', 'SEVERITY']);
    const overGeomY = await gridGeometry(page);
    expect(overGeomY.rows).toBe(axisViewportCount(overRawY, false));
    expect(overGeomY.rows).toBeLessThan(overRawY);
    const scrolledY = await scrollSliderExtent(page, 'y');
    expect(scrolledY.present).toBeGreaterThan(0);
    expect(scrolledY.ratio).toBeGreaterThan(0);
    expect(scrolledY.ratio).toBeLessThan(fitsY.ratio);
   } finally {

    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      try {
        tp.props.packCategories = true;
        await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 600);
      } catch (_) {  }
    });
    await restoreCanonical();
   }
  });

  await softStep('Legend', async () => {
   try {

    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const legendVisible = () => {
        const el = tp.root.querySelector('[name="legend"]') as HTMLElement | null;
        if (!el) return false;
        const cs = getComputedStyle(el);
        const rect = el.getBoundingClientRect();
        return cs.display !== 'none' && cs.visibility !== 'hidden' &&
          el.offsetParent !== null && rect.width > 0 && rect.height > 0;
      };
      const settle = (cap: number) => (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, cap);
      tp.props.viewerType = 'Scatter plot';
      tp.setOptions({innerViewerLook: {colorColumnName: 'SEX'}});
      await settle(1000);

      tp.props.legendVisibility = 'Always';
      await settle(600);
      const always = {vis: tp.props.legendVisibility, rendered: legendVisible()};

      const legendSlot = () => {
        const el = tp.root.querySelector('[name="legend"]');
        const side = el && el.closest('.d4-layout-left, .d4-layout-right, .d4-layout-top, .d4-layout-bottom');
        if (!side) return null;
        return ['left', 'right', 'top', 'bottom'].find((s) => side.classList.contains(`d4-layout-${s}`)) || null;
      };
      const positions: string[] = [];
      const slots: (string | null)[] = [];
      for (const pos of ['Left', 'Right', 'Top', 'Bottom']) {
        tp.props.legendPosition = pos;
        await settle(900);
        positions.push(tp.props.legendPosition);
        slots.push(legendSlot());
      }

      tp.props.legendVisibility = 'Never';
      await settle(600);
      const never = {vis: tp.props.legendVisibility, rendered: legendVisible()};
      return {always, positions, slots, never};
    });
    expect(result.always.vis).toBe('Always');
    expect(result.always.rendered).toBe(true);
    expect(result.positions).toEqual(['Left', 'Right', 'Top', 'Bottom']);
    expect(result.slots).toEqual(['left', 'right', 'top', 'bottom']);
    expect(result.never.vis).toBe('Never');
    expect(result.never.rendered).toBe(false);

    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      tp.props.viewerType = 'Box plot';
      tp.props.legendVisibility = 'Always';
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 2500);
    });
    await expect(page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell')).toHaveCount(canonicalCellCount);

    await openBoxPlotTab(page);
    await expect(
      page.locator('[name="viewer-Trellis-plot"] [name="legend"] .d4-legend-item .d4-legend-cross').first(),
    ).toBeAttached({timeout: 15000});

    const setAll = await setShowAllCategories(page, true);
    expect(setAll).toBe(true);

    const cats = await legendCategories(page);
    console.log(`[Legend] legend categories = ${JSON.stringify(cats)}`);

    expect(cats.length).toBeGreaterThanOrEqual(3);
    const pristine = await legendSnapshot(page);
    expect(pristine.every((it) => it.opacity >= 1)).toBe(true);
    expect(pristine.some((it) => it.current)).toBe(false);
    const firstBefore = await legendItemState(page, cats[0]);
    expect(firstBefore.present).toBe(true);
    expect(firstBefore.opacity).toBe(1);
    expect(firstBefore.active).toBe(false);
    const firstUnchecked = await uncheckLegendCategory(page, cats[0]);
    expect(firstUnchecked.clicked).toBe(true);
    expect(firstUnchecked.opacity).toBeLessThan(1);

    const afterFirstSnap = await legendSnapshot(page);
    expect(afterFirstSnap.filter((it) => it.name !== cats[0]).every((it) => it.current && it.opacity >= 1)).toBe(true);

    const afterFirst = await showAllCategoriesState(page);
    console.log(`[Legend] showAllCategories after single uncheck = ${afterFirst}`);
    expect(afterFirst).toBe(true);

    const errBeforeSeq = consoleErrors.length;
    const pageErrBeforeSeq = pageErrors.length;
    const secondBefore = await legendItemState(page, cats[1]);

    expect(secondBefore.opacity).toBe(1);
    expect(secondBefore.active).toBe(true);
    const secondUnchecked = await uncheckLegendCategory(page, cats[1]);
    expect(secondUnchecked.clicked).toBe(true);
    expect(secondUnchecked.opacity).toBeLessThan(1);

    const firstStillOff = await legendItemState(page, cats[0]);
    expect(firstStillOff.present).toBe(true);
    expect(firstStillOff.opacity).toBeLessThan(1);
    const afterSecond = await showAllCategoriesState(page);
    console.log(`[Legend] showAllCategories after two sequential unchecks = ${afterSecond}`);
    expect(afterSecond).toBe(true);
    expect(consoleErrors.length).toBe(errBeforeSeq);
    expect(pageErrors.length).toBe(pageErrBeforeSeq);

    const onlyFirst = await selectOnlyLegendCategory(page, cats[0]);
    expect(onlyFirst.clicked).toBe(true);
    expect(onlyFirst.opacity).toBeGreaterThanOrEqual(1);
    expect(onlyFirst.active).toBe(true);
    const exclusiveSnap = await legendSnapshot(page);
    expect(exclusiveSnap.filter((it) => it.opacity >= 1).map((it) => it.name)).toEqual([cats[0]]);
    expect(exclusiveSnap.find((it) => it.name === cats[2])!.opacity).toBeLessThan(1);

    expect(await resetLegendViaLastCross(page, cats[0])).toBe(true);
    const resetSnap = await legendSnapshot(page);
    expect(resetSnap.length).toBe(cats.length);
    expect(resetSnap.every((it) => it.opacity >= 1)).toBe(true);
    expect(resetSnap.some((it) => it.current)).toBe(false);

    const afterRecheck = await showAllCategoriesState(page);
    console.log(`[Legend] showAllCategories after exclusive select + reset = ${afterRecheck}`);
    expect(afterRecheck).toBe(true);
    expect(consoleErrors.length).toBe(errBeforeSeq);
    expect(pageErrors.length).toBe(pageErrBeforeSeq);
   } finally {

    await page.evaluate(async () => {
      try {
        const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
        tp.props.legendVisibility = 'Never';
        await new Promise((r) => setTimeout(r, 400));
      } catch (_) {  }
    });
    await restoreCanonical();
   }
  });

  await softStep('Context menu', async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.viewerType = 'Scatter plot';
      await w.__quiet('viewer:Trellis plot.onViewerRendered', 300, 800);
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const cell = root.querySelectorAll('.d4-trellis-plot-cell')[0];
      const r = cell.getBoundingClientRect();
      cell.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2, clientX: r.left + r.width / 2, clientY: r.top + r.height / 2}));
      // settle on the context menu popup attaching, capped
      await w.__poll(() => document.querySelector('.d4-menu-popup .d4-menu-item-label'),
        (e: Element | null) => !!e, 900, 40);

      const labels = Array.from(document.querySelectorAll('.d4-menu-item-label')).map((el) => (el as HTMLElement).textContent?.trim());
      document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
      return {
        hasInnerGroup: labels.includes('Scatter plot'),
        hasProperties: labels.includes('Properties...'),

        hasLasso: labels.includes('Lasso Tool'),
        hasRegression: labels.includes('Show Regression Line'),
      };
    });

    expect(result.hasInnerGroup).toBe(true);
    expect(result.hasProperties).toBe(true);
    expect(result.hasLasso).toBe(true);
    expect(result.hasRegression).toBe(true);
  });

  await softStep('Inner viewer properties', async () => {
   try {

    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      tp.setOptions({innerViewerLook: {xColumnName: 'WEIGHT', yColumnName: 'HEIGHT'}});
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1200);
    });

    await v.openViewerGear(page, 'Trellis plot');
    await page.locator('.property-grid').first().waitFor({timeout: 10000});

    const tabs = await page.evaluate(() => ({
      hasTrellisTab: !!document.querySelector('.d4-tab-header[name="Trellis"]'),
      hasInnerTab: !!document.querySelector('.d4-tab-header[name="Scatter plot"]'),
    }));
    expect(tabs.hasTrellisTab).toBe(true);
    expect(tabs.hasInnerTab).toBe(true);

    const idxs = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const cells = root.querySelectorAll('.d4-trellis-plot-cell');
      const out: number[] = [];
      for (let i = 0; i < cells.length && out.length < 2; i++) if (cells[i].querySelector('canvas')) out.push(i);
      return out;
    });
    const before = await v.trellisCellHashes(page, idxs);
    const after = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      await (window as any).__settled('viewer:Trellis plot.onViewerRendered',
        () => tp.setOptions({innerViewerLook: {xColumnName: 'AGE', yColumnName: 'WEIGHT'}}), 1500);
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1500);
      return {type: tp.props.viewerType};
    });
    const afterHashes = await v.trellisCellHashes(page, idxs);
    expect(after.type).toBe('Scatter plot');
    expect(idxs.length).toBeGreaterThan(0);
    for (let i = 0; i < idxs.length; i++) {
      expect(before[i]).not.toBeNull();
      expect(afterHashes[i]).not.toBeNull();
      expect(afterHashes[i]).not.toBe(before[i]);
    }
   } finally {
    await restoreCanonical();
   }
  });

  await softStep('Use in Trellis', async () => {

    try {

      await page.evaluate(async () => {
        const w = window as any;
        const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
        if (tp) tp.close();
        await w.__poll(() => document.querySelector('[name="viewer-Trellis-plot"]'), (e: Element | null) => !e, 600, 40);
        const sp = grok.shell.tv.addViewer('Scatter plot') as any;
        sp.setOptions({xColumnName: 'AGE', yColumnName: 'HEIGHT', colorColumnName: 'SEX'});
        // settle on the scatter plot attaching, capped
        await w.__poll(() => document.querySelector('[name="viewer-Scatter-plot"]'), (e: Element | null) => !!e, 1200, 60);
      });
      await page.locator('[name="viewer-Scatter-plot"]').first().click({button: 'right'});
      await clickMenuItemInGroup(page, 'General', 'Use in Trellis');
      const scatterResult = await page.evaluate(async () => {
        const w = window as any;
        let newTp: any = null;
        for (let i = 0; i < 40 && !newTp; i++) {
          newTp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot');
          if (!newTp) await new Promise((res) => setTimeout(res, 250));
        }
        // settle on the new trellis finishing its first paint, capped
        await w.__quiet('viewer:Trellis plot.onViewerRendered', 300, 800);

        let ivl: any = null;
        try {
          const opts = newTp ? newTp.getOptions(true) : null;
          const raw = opts?.look?.innerViewerLook;
          ivl = typeof raw === 'string' ? JSON.parse(raw) : raw;
        } catch {}
        for (const vw of Array.from(grok.shell.tv.viewers) as any[]) if (vw.type === 'Scatter plot') vw.close();
        return {trellisCreated: !!newTp, innerType: newTp?.props.viewerType,
          preservedX: ivl?.xColumnName, preservedY: ivl?.yColumnName, preservedColor: ivl?.colorColumnName};
      });
      expect(scatterResult.trellisCreated).toBe(true);
      expect(scatterResult.innerType).toBe('Scatter plot');
      expect(scatterResult.preservedX).toBe('AGE');
      expect(scatterResult.preservedY).toBe('HEIGHT');
      expect(scatterResult.preservedColor).toBe('SEX');

      const useInTrellis = async (viewerType: string, look: any) => {
        await page.evaluate(async ({viewerType, look}) => {
          const w = window as any;
          for (const vw of Array.from(grok.shell.tv.viewers) as any[]) if (vw.type === 'Trellis plot') vw.close();
          await w.__poll(() => document.querySelector('[name="viewer-Trellis-plot"]'), (e: Element | null) => !e, 500, 40);
          const src = grok.shell.tv.addViewer(viewerType) as any;
          src.setOptions(look);
          const sel = '[name="viewer-' + viewerType.replace(/\s+/g, '-') + '"]';
          await w.__poll(() => document.querySelector(sel), (e: Element | null) => !!e, 1200, 60);
        }, {viewerType, look});
        await page.locator('[name="viewer-' + viewerType.replace(/\s+/g, '-') + '"]').first().click({button: 'right'});
        await clickMenuItemInGroup(page, 'General', 'Use in Trellis');
        return page.evaluate(async (viewerType) => {
          const w = window as any;
          let newTp: any = null;
          for (let i = 0; i < 40 && !newTp; i++) {
            newTp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot');
            if (!newTp) await new Promise((res) => setTimeout(res, 250));
          }
          await w.__quiet('viewer:Trellis plot.onViewerRendered', 300, 800);
          const innerType = newTp?.props.viewerType;
          for (const vw of Array.from(grok.shell.tv.viewers) as any[]) if (vw.type === viewerType) vw.close();
          return {trellisCreated: !!newTp, innerType};
        }, viewerType);
      };

      const barResult: any = await useInTrellis('Bar chart', {splitColumnName: 'RACE', valueColumnName: 'AGE'});
      expect(barResult.trellisCreated).toBe(true);
      expect(barResult.innerType).toBe('Bar chart');

      const histoResult: any = await useInTrellis('Histogram', {valueColumnName: 'AGE', splitColumnName: 'RACE'});
      expect(histoResult.trellisCreated).toBe(true);
      expect(histoResult.innerType).toBe('Histogram');

      const lineResult: any = await useInTrellis('Line chart', {xColumnName: 'STARTED', yColumnName: 'AGE'});
      expect(lineResult.trellisCreated).toBe(true);
      expect(lineResult.innerType).toBe('Line chart');

      const boxResult: any = await useInTrellis('Box plot', {categoryColumnName: 'SEX', valueColumnName: 'AGE'});
      expect(boxResult.trellisCreated).toBe(true);
      expect(boxResult.innerType).toBe('Box plot');
    } finally {

      await page.keyboard.press('Escape').catch(() => {});
      await restoreCanonical();
    }
  });

  await softStep('Auto layout', async () => {
   try {

    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.autoLayout = true;
      tp.props.showXLabels = true;
      tp.props.showYLabels = true;
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 600);
    });
    const vsVisible = () => page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const el = root?.querySelector('[name="viewer selector"]') as HTMLElement | null;
      if (!el) return false;
      const b = el.getBoundingClientRect();
      return b.width > 0 && b.height > 0;
    });
    const viewerWidth = () => page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
      return root ? Math.round(root.getBoundingClientRect().width) : -1;
    });

    const cards = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {sex: df.col('SEX').categories.length, race: df.col('RACE').categories.length};
    });
    const expectedX = axisViewportCount(cards.sex, false);
    const expectedY = axisViewportCount(cards.race, false);

    expect(await vsVisible()).toBe(true);
    const wideLabels = await categoryLabels(page);
    expect(wideLabels.x.length).toBe(expectedX);
    expect(wideLabels.y.length).toBe(expectedY);

    await page.setViewportSize({width: 500, height: 400});

    await page.waitForFunction(() => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
      if (!root) return false;
      return Array.from(root.querySelectorAll('.d4-trellis-plot-cat-item-horz'))
        .filter((n) => n.tagName.toLowerCase() === 'text').length === 0;
    }, null, {timeout: 1500}).catch(() => {});
    const smallVisible = await vsVisible();
    const smallLabels = await categoryLabels(page);
    expect(smallVisible).toBe(false);
    expect(smallLabels.x.length).toBe(0);
    expect(smallLabels.y.length).toBe(0);

    await page.setViewportSize({width: 1920, height: 1080});

    await page.waitForFunction((n) => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
      if (!root) return false;
      return Array.from(root.querySelectorAll('.d4-trellis-plot-cat-item-horz'))
        .filter((el) => el.tagName.toLowerCase() === 'text').length === n;
    }, expectedX, {timeout: 1500}).catch(() => {});
    const largeVisible = await vsVisible();
    const restoredLabels = await categoryLabels(page);
    expect(largeVisible).toBe(true);
    expect(restoredLabels.x.length).toBe(expectedX);
    expect(restoredLabels.y.length).toBe(expectedY);

    // The band search keeps HEAD's fixed settles: waitForViewerRendered returns before the
    // resize relayout finishes, so an event-waited probe reads an intermediate viewer width
    // and the ~20px band is never observed. Measured 2026-08-31 (2/8 passes event-waited).
    const wideWidth = await viewerWidth();
    await page.setViewportSize({width: 1420, height: 1080});

    await page.waitForTimeout(1200);
    const midWidth = await viewerWidth();
    const slope = Math.max((wideWidth - midWidth) / 500, 0.05);

    console.log(`[Auto layout] band calibration: wideWidth=${wideWidth} midWidth=${midWidth} slope=${slope.toFixed(3)}`);
    let band: {window: number; viewer: number; x: number; y: number} | null = null;
    for (let target = 235; target >= 150 && !band; target -= 5) {
      const win = Math.round(Math.min(1900, Math.max(420, 1420 - (midWidth - target) / slope)));
      await page.setViewportSize({width: win, height: 1080});

      await page.waitForTimeout(1000);
      const seen = await categoryLabels(page);
      console.log(`[Auto layout] band probe: targetViewer=${target} window=${win} viewer=${await viewerWidth()} ` +
        `x=${seen.x.length} y=${seen.y.length}`);
      if (seen.x.length > 0 || seen.y.length === 0) continue;

      await page.waitForTimeout(900);
      const confirmed = await categoryLabels(page);
      console.log(`[Auto layout] band candidate re-read: window=${win} viewer=${await viewerWidth()} ` +
        `x=${confirmed.x.length} y=${confirmed.y.length}`);
      if (confirmed.x.length === 0 && confirmed.y.length > 0)
        band = {window: win, viewer: await viewerWidth(), x: confirmed.x.length, y: confirmed.y.length};
    }
    console.log(`[Auto layout] band search finished: ${JSON.stringify(band)}`);
    expect(band).not.toBeNull();
    expect(band!.x).toBe(0);
    expect(band!.y).toBe(expectedY);

    await page.setViewportSize({width: 1920, height: 1080});
    await v.waitForViewerRendered(page, 'Trellis plot', 1200);

    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.autoLayout = false;
      await new Promise((r) => setTimeout(r, 500));
    });
    await page.setViewportSize({width: 500, height: 400});

    await page.waitForFunction((n) => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
      if (!root) return false;
      return Array.from(root.querySelectorAll('.d4-trellis-plot-cat-item-horz'))
        .filter((el) => el.tagName.toLowerCase() === 'text').length === n;
    }, expectedX, {timeout: 1500}).catch(() => {});
    const offSmallVisible = await vsVisible();
    const offSmallLabels = await categoryLabels(page);
    expect(offSmallVisible).toBe(true);
    expect(offSmallLabels.x.length).toBe(expectedX);
    expect(offSmallLabels.y.length).toBe(expectedY);
   } finally {

    await page.setViewportSize({width: 1920, height: 1080});
    await v.waitForViewerRendered(page, 'Trellis plot', 1000);
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      if (tp) {
        tp.props.autoLayout = true;
        tp.props.showXLabels = true;
        tp.props.showYLabels = true;
      }
      await new Promise((r) => setTimeout(r, 500));
    });
   }
  });

  await softStep('Title and description', async () => {
   try {

    const titleShown = (text: string) => page.evaluate((t) => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
      if (!root) return false;
      const el = root.querySelector('.d4-viewer-title') as HTMLElement | null;
      if (el) {
        const editor = el.querySelector('textarea, input') as HTMLInputElement | null;
        const shown = ((editor ? editor.value : el.textContent) ?? '').trim();
        return shown === t && el.getBoundingClientRect().height > 0;
      }
      return Array.from(root.querySelectorAll('*')).some((n) =>
        n.children.length === 0 && (n.textContent ?? '').trim() === t);
    }, text);
    const descriptionSlot = () => page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
      const el = root?.querySelector('.d4-viewer-description');
      const side = el && el.closest('.d4-layout-left, .d4-layout-right, .d4-layout-top, .d4-layout-bottom');
      if (!side) return null;
      return ['left', 'right', 'top', 'bottom'].find((s) => side.classList.contains(`d4-layout-${s}`)) ?? null;
    });

    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.showTitle = true;
      tp.props.title = 'My Trellis';
      tp.props.description = 'Test description';
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1500);
    });
    expect(await titleShown('My Trellis')).toBe(true);

    const slots: (string | null)[] = [];
    for (const pos of ['Bottom', 'Top', 'Left', 'Right']) {
      await page.evaluate(async (p) => {
        const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
        tp.props.descriptionPosition = p;
        await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1200);
      }, pos);
      slots.push(await descriptionSlot());
    }
    expect(slots).toEqual(['bottom', 'top', 'left', 'right']);

    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.showTitle = false;
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1200);
    });
    expect(await titleShown('My Trellis')).toBe(false);
   } finally {
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      try {
        tp.props.description = '';
        tp.props.title = '';
        tp.props.showTitle = false;
        await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 600);
      } catch (_) {  }
    });
   }
  });

  await softStep('Label orientation', async () => {
   try {

    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['RACE'];
      tp.props.yColumnNames = ['SEX'];
      tp.props.showXLabels = true;
      tp.props.showYLabels = true;
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1800);
    });
    const setOrientation = async (axis: 'x' | 'y', value: string) => {
      await page.evaluate(async (o) => {
        const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
        if (o.axis === 'x') tp.props.xLabelsOrientation = o.value;
        else tp.props.yLabelsOrientation = o.value;
        await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1500);
      }, {axis, value});
      return categoryLabels(page);
    };

    const xHorz = await setOrientation('x', 'Horz');
    expect(xHorz.xAngles.length).toBeGreaterThan(0);
    expect(xHorz.xAngles.every((a) => a === 0)).toBe(true);
    const xVert = await setOrientation('x', 'Vert');
    expect(xVert.xAngles.length).toBeGreaterThan(0);
    expect(xVert.xAngles.every((a) => a === -90)).toBe(true);

    const yHorz = await setOrientation('y', 'Horz');
    expect(yHorz.yAngles.length).toBeGreaterThan(0);
    expect(yHorz.yAngles.every((a) => a === 0)).toBe(true);
    const yVert = await setOrientation('y', 'Vert');
    expect(yVert.yAngles.length).toBeGreaterThan(0);
    expect(yVert.yAngles.every((a) => a === -90)).toBe(true);

    const xAuto = await setOrientation('x', 'Auto');
    expect(xAuto.xAngles.length).toBeGreaterThan(0);
    expect(xAuto.xAngles.every((a) => a === 0 || a === -90)).toBe(true);
    const yAuto = await setOrientation('y', 'Auto');
    expect(yAuto.yAngles.length).toBeGreaterThan(0);
    expect(yAuto.yAngles.every((a) => a === 0 || a === -90)).toBe(true);
   } finally {
    await restoreCanonical();
   }
  });

  await softStep('Pick Up / Apply', async () => {
   try {

    await page.evaluate(async () => {
      const tp1 = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      grok.shell.tv.addViewer('Trellis plot');
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1200);

      tp1.props.yColumnNames = ['RACE'];
      tp1.props.viewerType = 'Bar chart';
      tp1.setOptions({innerViewerLook: {splitColumnName: 'RACE', valueColumnName: 'AGE'}});
      tp1.props.legendVisibility = 'Always';
      tp1.props.legendPosition = 'Top';
      tp1.props.showTitle = true;
      tp1.props.title = 'First Trellis';
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1000);
    });

    const trellisCell = (n: number) => page.locator('[name="viewer-Trellis-plot"]').nth(n).locator('.d4-trellis-plot-cell').first();
    await trellisCell(0).click({button: 'right', position: {x: 6, y: 6}});
    await clickMenuItemInGroup(page, 'Pick Up / Apply', 'Pick Up');
    // the picked-up look is held in memory; settle on the menu closing before re-opening it
    await page.locator('.d4-menu-popup').last().waitFor({state: 'detached', timeout: 800}).catch(() => {});
    await trellisCell(1).click({button: 'right', position: {x: 6, y: 6}});
    await clickMenuItemInGroup(page, 'Pick Up / Apply', 'Apply');
    await v.waitForViewerRendered(page, 'Trellis plot', 900);

    const result: any = await page.evaluate(async () => {
      const tps = Array.from(grok.shell.tv.viewers).filter((v: any) => v.type === 'Trellis plot') as any[];
      const applied = {type: tps[1]?.props.viewerType, title: tps[1]?.props.title, legendPos: tps[1]?.props.legendPosition,
        x: [...(tps[1]?.props.xColumnNames ?? [])], y: [...(tps[1]?.props.yColumnNames ?? [])]};
      const source = {x: [...tps[0].props.xColumnNames], y: [...tps[0].props.yColumnNames]};

      const tp2TypeBefore = tps[1]?.props.viewerType;
      tps[0].props.xColumnNames = ['CONTROL'];
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 1000);
      const tp2XAfterFirstChange = [...tps[1].props.xColumnNames];
      const tp2Independent = tps[1]?.props.viewerType === tp2TypeBefore;

      return {applied, source, tp2XAfterFirstChange, tp2Independent};
    });
    expect(result.applied.type).toBe('Bar chart');
    expect(result.applied.title).toBe('First Trellis');
    expect(result.applied.legendPos).toBe('Top');
    expect(result.applied.x).toEqual(result.source.x);
    expect(result.applied.y).toEqual(result.source.y);
    expect(result.tp2Independent).toBe(true);
    expect(result.tp2XAfterFirstChange).not.toContain('CONTROL');

    const step8Setup = await page.evaluate(async () => {
      const tps = Array.from(grok.shell.tv.viewers).filter((v: any) => v.type === 'Trellis plot') as any[];
      const tp1 = tps[0], tp2 = tps[1];
      tp1.props.viewerType = 'Scatter plot';
      tp1.props.xColumnNames = ['SEX'];
      tp1.props.yColumnNames = ['RACE'];

      tp2.props.viewerType = 'Scatter plot';
      tp2.props.xColumnNames = ['SEX'];
      tp2.props.yColumnNames = ['RACE'];
      tp2.setOptions({innerViewerLook: {xColumnName: 'WEIGHT', yColumnName: 'HEIGHT'}});
      tp2.props.globalScale = true;
      tp2.props.showRangeSliders = true;
      tp2.props.showXAxes = 'Always';
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 2200);
      const root2 = document.querySelectorAll('[name="viewer-Trellis-plot"]')[1] as HTMLElement;
      const innerX = root2.querySelector('.d4-range-selector > svg[type="range-slider"][name="x-slider"]') as SVGElement | null;
      let box: {x: number; y: number; w: number; h: number} | null = null;
      if (innerX) {
        const wrap = innerX.closest('.d4-range-selector') as HTMLElement;
        const wb = wrap.getBoundingClientRect();
        wrap.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: wb.left + wb.width / 2, clientY: wb.top + wb.height / 2}));
        await new Promise((r) => setTimeout(r, 200));
        const b = innerX.getBoundingClientRect();
        box = {x: b.x, y: b.y, w: b.width, h: b.height};
      }
      return {box, look1: JSON.stringify(tp1.getOptions(true)?.look ?? null)};
    });
    const setupHash1 = await v.trellisAllCellHashes(page, {rootIndex: 0});
    const setupHash2 = await v.trellisAllCellHashes(page, {rootIndex: 1});
    expect(step8Setup.box).not.toBeNull();
    expect(step8Setup.box!.w).toBeGreaterThan(0);
    await page.mouse.move(step8Setup.box!.x + step8Setup.box!.w - 4, step8Setup.box!.y + step8Setup.box!.h / 2);
    await page.mouse.down();
    await page.mouse.move(step8Setup.box!.x + step8Setup.box!.w * 0.45, step8Setup.box!.y + step8Setup.box!.h / 2, {steps: 12});
    await page.mouse.up();

    await v.waitForViewerRendered(page, 'Trellis plot', 1500);
    const afterHash1 = await v.trellisAllCellHashes(page, {rootIndex: 0});
    const afterHash2 = await v.trellisAllCellHashes(page, {rootIndex: 1});
    const step8After = {
      hash1: afterHash1,
      hash2: afterHash2,
      look1: await page.evaluate(() => {
        const tps = Array.from(grok.shell.tv.viewers).filter((v: any) => v.type === 'Trellis plot') as any[];
        return JSON.stringify(tps[0].getOptions(true)?.look ?? null);
      }),
    };
    const setupHashes = {hash1: setupHash1, hash2: setupHash2, look1: step8Setup.look1};
    const secondChanged = setupHashes.hash2.some((h, i) => h !== null && step8After.hash2[i] !== null && h !== step8After.hash2[i]);
    expect(secondChanged).toBe(true);

    const comparable1 = setupHashes.hash1.filter((h, i) => h !== null && step8After.hash1[i] !== null).length;
    const firstMoved = setupHashes.hash1
      .map((h, i) => (h !== null && step8After.hash1[i] !== null && h !== step8After.hash1[i] ? i : -1))
      .filter((i) => i >= 0);
    expect(comparable1).toBeGreaterThan(0);
    expect(firstMoved).toEqual([]);

    expect(step8After.look1).toBe(setupHashes.look1);
   } finally {

    await restoreCanonical();
   }
  });

  await softStep('Layout and Project save/restore', async () => {
   try {

    const result = await page.evaluate(async () => {
      const r: any = {};
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      r.viewersAtSave = Array.from(grok.shell.tv.viewers).map((v: any) => v.type).sort();

      const beforeAddCount = grok.shell.tv.viewers.length;
      grok.shell.tv.addViewer('Histogram');
      grok.shell.tv.addViewer('Bar chart');
      // settle on both viewers actually attaching, capped — not a blind 1200ms
      await (window as any).__poll(() => grok.shell.tv.viewers.length,
        (n: number) => n >= beforeAddCount + 2, 1200, 60);
      r.viewersBefore = Array.from(grok.shell.tv.viewers).map((v: any) => v.type).sort();

      const saved = await grok.dapi.layouts.find(layoutId);

      const applied = new Promise<void>((res) => {
        let sub: any = null;
        try { sub = grok.events.onViewLayoutApplied.subscribe(() => { sub.unsubscribe(); res(); }); }
        catch (_) {  }
        setTimeout(() => { try { sub?.unsubscribe(); } catch (_) {} res(); }, 3000);
      });
      grok.shell.tv.loadLayout(saved);
      await applied;
      r.viewersAfter = Array.from(grok.shell.tv.viewers).map((v: any) => v.type).sort();

      await grok.dapi.layouts.delete(saved);
      return r;
    });
    expect(result.viewersBefore.length).toBeGreaterThan(result.viewersAfter.length);
    expect(result.viewersAfter).toEqual(result.viewersAtSave);
    expect(result.viewersAfter).toContain('Trellis plot');

    const errBeforeSave = consoleErrors.length;
    const pageErrBeforeSave = pageErrors.length;
    const saveBtn = page.locator('[name="button-Save"]').first();
    await saveBtn.waitFor({state: 'visible', timeout: 15000});
    await saveBtn.click();

    await page.locator('.d4-dialog').first().waitFor({state: 'visible', timeout: 2000}).catch(() => {});
    const dialogOpen = await page.evaluate(() => !!document.querySelector('.d4-dialog'));
    const errAfterSave = consoleErrors.length;
    const pageErrAfterSave = pageErrors.length;

    await page.locator('[name="button-CANCEL"]').first().click({timeout: 5000}).catch(() => {});
    await page.locator('.d4-dialog').first().waitFor({state: 'detached', timeout: 500}).catch(() => {});
    expect(dialogOpen).toBe(true);
    expect(errAfterSave).toBe(errBeforeSave);
    expect(pageErrAfterSave).toBe(pageErrBeforeSave);
   } finally {

    await restoreCanonical();
   }
  });

  await softStep('Viewer filter formula', async () => {

    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true); df.rows.requestFilter();
      await new Promise((r) => setTimeout(r, 500));
      const dfBefore = df.filter.trueCount;
      const r: any[] = [];
      tp.props.filter = '${AGE} > 40';
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 600);
      r.push({filter: tp.props.filter, dfCount: df.filter.trueCount});
      tp.props.filter = '';
      await (window as any).__quiet('viewer:Trellis plot.onViewerRendered', 300, 600);
      r.push({filter: tp.props.filter, dfCount: df.filter.trueCount});
      return {r, dfBefore};
    });
    expect(result.r[0].filter).toBe('${AGE} > 40');
    expect(result.r[1].filter).toBe('');
    expect(result.r[0].dfCount).toBe(result.dfBefore);
  });

  await softStep('Multi Curve inner viewer', async () => {
    const result = await page.evaluate(async (cPath) => {

      const hasSexRace = (view: any) => { try { const df = view.dataFrame; return !!(df && df.col('SEX') && df.col('RACE')); } catch { return false; } };
      const hasTrellis = (view: any) => { try { return Array.from(view.viewers ?? []).some((x: any) => x.type === 'Trellis plot'); } catch { return false; } };
      const anchorDemog = (): any => {
        const views = Array.from(grok.shell.views ?? []) as any[];
        const withTrellis = views.find((v) => hasSexRace(v) && hasTrellis(v));
        const target = withTrellis ?? views.find(hasSexRace);
        if (target) { grok.shell.v = target; return target; }
        return grok.shell.tv;
      };
      const w = window as any;
      const settle = (cap: number) => w.__quiet('viewer:Trellis plot.onViewerRendered', 300, cap);
      const demogView = anchorDemog();
      await settle(600);
      const demogName = grok.shell.tv.dataFrame.name;
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      if (!tp) return {error: 'Trellis plot not found on demog view'};
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      await settle(1200);

      const dfCurves = await grok.dapi.files.readCsv(cPath);
      grok.shell.addTableView(dfCurves);
      // settle on the new grid attaching, capped
      await w.__poll(() => document.querySelector('.d4-grid[name="viewer-Grid"]'),
        (e: Element | null) => !!e, 1800, 80);
      grok.shell.v = demogView;
      await settle(800);

      let switchError: string | null = null;
      // rebind the trellis to the curves table — settle on its repaint, capped
      await w.__settled('viewer:Trellis plot.onViewerRendered', () => {
        try { tp.props.table = dfCurves.name; } catch (e) { switchError = String(e); }
      }, 1500);
      await settle(1500);

      const boundToCurves = (() => {
        try { const d = tp.dataFrame; return {name: d?.name ?? null, rows: d?.rowCount ?? -1}; }
        catch { return {name: null, rows: -1}; }
      })();

      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const vs = root.querySelector('[name="viewer selector"]') as HTMLElement;
      vs.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
      await w.__poll(() => document.querySelector('.d4-combo-drop-down'), (e: Element | null) => !!e, 600, 40);
      const mc = document.querySelector('.d4-combo-drop-down [name="icon-multicurveviewer"]');
      const mcClicked = !!mc;
      await w.__settled('viewer:Trellis plot.onViewerRendered',
        () => (mc?.closest('.d4-list-item') as HTMLElement | null)?.click(), 1800);
      await settle(1800);
      const vt = tp.props.viewerType;

      let gearOpenedPropGrid = false;
      const panel = root.closest('.panel-base') as HTMLElement | null;
      const gear = panel?.querySelector('.panel-titlebar [name="icon-font-icon-settings"]') as HTMLElement | null;
      if (gear) {
        gear.click();
        await w.__poll(() => document.querySelector('.property-grid'), (e: Element | null) => !!e, 1200, 60);
        gearOpenedPropGrid = !!document.querySelector('.property-grid');
      }

      let restoreError: string | null = null;
      await w.__settled('viewer:Trellis plot.onViewerRendered', () => {
        try { tp.props.table = demogName; } catch (e) { restoreError = String(e); }
      }, 1800);
      await settle(1800);
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      await settle(1800);
      const boundBack = (() => { try { return tp.dataFrame?.name ?? null; } catch { return null; } })();
      const restoredCells = document.querySelectorAll('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell').length;
      return {viewerType: vt, mcClicked, switchError, boundToCurves, restoreError, boundBack,
        curvesName: dfCurves.name, curvesRows: dfCurves.rowCount, demogName,
        restoredCells, gearOpenedPropGrid};
    }, curvesPath);
    expect(result.mcClicked).toBe(true);
    expect(['MultiCurveViewer', 'Multi curve viewer', 'Curves'].includes(result.viewerType)).toBe(true);
    expect(result.gearOpenedPropGrid).toBe(true);

    console.log(`[Multi Curve] table switch: error=${result.switchError} ` +
      `bound=${JSON.stringify(result.boundToCurves)} expected={"name":"${result.curvesName}","rows":${result.curvesRows}}`);
    console.log(`[Multi Curve] restore: error=${result.restoreError} boundBack=${result.boundBack} expected=${result.demogName}`);
    expect(result.switchError).toBeNull();
    expect(result.boundToCurves.name).toBe(result.curvesName);
    expect(result.boundToCurves.rows).toBe(result.curvesRows);
    expect(result.restoreError).toBeNull();
    expect(result.boundBack).toBe(result.demogName);
    expect(result.restoredCells).toBe(canonicalCellCount);
  });

  await softStep('To Script', async () => {

    await page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell').first().click({button: 'right', position: {x: 6, y: 6}});
    await clickMenuItemInGroup(page, 'To Script', 'To JavaScript');
    const result = await page.evaluate(async () => {
      // settle on the To-Script balloon appearing, capped — not a blind 1500ms
      const balloon = await (window as any).__poll(() => document.querySelector('.d4-balloon'),
        (e: Element | null) => !!e, 1500, 60) as Element | null;
      const generated = !!balloon;
      if (balloon) {
        const close = balloon.querySelector('.close') || balloon.querySelector('[name="icon-times"]');
        if (close) (close as HTMLElement).click();
      }
      return {scriptGenerated: generated};
    });
    expect(result.scriptGenerated).toBe(true);
  });

  await softStep('Keyboard navigation', async () => {
    const cellLocator = page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell');

    await page.evaluate(async () => {
      const hasSexRace = (view: any) => { try { const df = view.dataFrame; return !!(df && df.col('SEX') && df.col('RACE')); } catch { return false; } };
      const hasTrellis = (view: any) => { try { return Array.from(view.viewers ?? []).some((x: any) => x.type === 'Trellis plot'); } catch { return false; } };
      const w = window as any;
      const views = Array.from(grok.shell.views ?? []) as any[];
      const target = views.find((v) => hasSexRace(v) && hasTrellis(v)) ?? views.find(hasSexRace);
      if (target) { try { grok.shell.v = target; } catch {} }
      await new Promise((r) => setTimeout(r, 400));

      for (const vw of Array.from(grok.shell.tv.viewers) as any[]) if (vw.type === 'Trellis plot') vw.close();
      await w.__poll(() => document.querySelector('[name="viewer-Trellis-plot"]'), (e: Element | null) => !e, 600, 40);
      const tp = grok.shell.tv.addViewer('Trellis plot') as any;
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      tp.props.onClick = 'None';
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true); df.selection.setAll(false); df.rows.requestFilter();
      await w.__quiet('viewer:Trellis plot.onViewerRendered', 300, 1800);
    });
    await expect(cellLocator).toHaveCount(canonicalCellCount);

    const navCats = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {x: [...df.col('SEX').categories], y: [...df.col('RACE').categories]};
    });
    expect(navCats.x.length).toBeGreaterThanOrEqual(2);
    expect(navCats.y.length).toBeGreaterThanOrEqual(2);
    const idx = await cellIndexFor(page, navCats.x[0], navCats.y[0]);
    const idxRight = await cellIndexFor(page, navCats.x[1], navCats.y[0]);
    const idxDown = await cellIndexFor(page, navCats.x[0], navCats.y[1]);
    expect(idxRight).not.toBe(idx);
    expect(idxDown).not.toBe(idx);
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});

    await page.waitForFunction((i) => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
      if (!root) return false;
      return Array.from(root.querySelectorAll('.d4-trellis-plot-cell'))
        .findIndex((c) => c.classList.contains('d4-trellis-cell-current')) === i;
    }, idx, {timeout: 700}).catch(() => {});
    expect(await currentCellIndex(page)).toBe(idx);
    await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      (window as any).__cc = [];
      (window as any).__ccSub = tp.onEvent('d4-trellis-plot-current-cell-changed').subscribe((a: any) => {
        const mc = (a && a.args && a.args.matchCondition) ? a.args.matchCondition : (a && a.matchCondition ? a.matchCondition : a);
        (window as any).__cc.push(mc);
      });
    });
    const focusGrid = () => focusChartsGrid(page);

    const arrow = async (key: string) => {
      await focusGrid();
      const before = await page.evaluate(() => ((window as any).__cc?.length ?? 0));
      await page.keyboard.press(key);

      await page.waitForFunction((n) => ((window as any).__cc?.length ?? 0) > n, before, {timeout: 500}).catch(() => {});
      return currentCellIndex(page);
    };
    const afterRight = await arrow('ArrowRight');
    const afterLeft = await arrow('ArrowLeft');
    const afterDown = await arrow('ArrowDown');
    const afterUp = await arrow('ArrowUp');
    const events = await page.evaluate(() => {
      (window as any).__ccSub?.unsubscribe?.();
      return (window as any).__cc as any[];
    });
    expect(afterRight).toBe(idxRight);
    expect(afterLeft).toBe(idx);
    expect(afterDown).toBe(idxDown);
    expect(afterUp).toBe(idx);

    expect(events.length).toBeGreaterThanOrEqual(4);
    expect(JSON.stringify(events[1])).toBe(JSON.stringify(events[3]));
    expect(JSON.stringify(events[0])).not.toBe(JSON.stringify(events[1]));
    expect(JSON.stringify(events[2])).not.toBe(JSON.stringify(events[1]));
    expect(JSON.stringify(events[0])).not.toBe(JSON.stringify(events[2]));

    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.onClick = 'Filter';
      df.filter.setAll(true); df.rows.requestFilter();
      await new Promise((r) => setTimeout(r, 500));
    });
    await armDfEvent(page, 'onRowsFiltered');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await awaitArmedDfEvent(page, 700);
    const filteredBeforeEsc = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filteredBeforeEsc).toBeLessThan(fullRowCount);

    await focusGrid();
    await armDfEvent(page, 'onRowsFiltered');
    await page.keyboard.press('Escape');
    await awaitArmedDfEvent(page, 900);
    const filteredAfterEsc = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filteredAfterEsc).toBe(fullRowCount);
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.onClick = 'None';
      await new Promise((r) => setTimeout(r, 400));
    });
  });

  await softStep('Undo/redo', async () => {
   try {

    await restoreCanonical();
    await expect(page.locator('[name="viewer-Trellis-plot"]')).toBeVisible();
    expect(await page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell').count()).toBeGreaterThan(0);
    expect(await trellisCount(page)).toBe(1);

    await closeTrellis(page);
    await waitForTrellisCount(page, 0, 1500);
    expect(await trellisCount(page)).toBe(0);
    await expect(page.locator('[name="viewer-Trellis-plot"]')).toHaveCount(0);

    let errBefore = consoleErrors.length;
    let pageErrBefore = pageErrors.length;
    await page.keyboard.press('Control+z');
    await waitForTrellisCount(page, 1, 1800);
    expect(await trellisCount(page)).toBe(1);
    await expect(page.locator('[name="viewer-Trellis-plot"]')).toBeVisible();
    expect(consoleErrors.length).toBe(errBefore);
    expect(pageErrors.length).toBe(pageErrBefore);

    errBefore = consoleErrors.length;
    pageErrBefore = pageErrors.length;
    await page.keyboard.press('Control+Shift+z');
    await waitForTrellisCount(page, 0, 1800);
    expect(await trellisCount(page)).toBe(0);
    expect(consoleErrors.length).toBe(errBefore);
    expect(pageErrors.length).toBe(pageErrBefore);
    expect(await page.locator('.d4-balloon.error').count()).toBe(0);

    await page.keyboard.press('Control+z');
    await waitForTrellisCount(page, 1, 1800);
    expect(await trellisCount(page)).toBe(1);

    errBefore = consoleErrors.length;
    pageErrBefore = pageErrors.length;
    await page.keyboard.press('Control+Shift+z');
    await waitForTrellisCount(page, 0, 1800);
    expect(await trellisCount(page)).toBe(0);
    expect(consoleErrors.length).toBe(errBefore);
    expect(pageErrors.length).toBe(pageErrBefore);
    expect(await page.locator('.d4-balloon.error').count()).toBe(0);
   } finally {

    await restoreCanonical();
   }
  });

  v.finishSpec();
});
