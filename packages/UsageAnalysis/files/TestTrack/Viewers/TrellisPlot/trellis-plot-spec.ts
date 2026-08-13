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

// PER-AXIS CATEGORY VIEWPORT CLAMP — never assume one cell per category: the rendered
// window, not n, is what an axis puts on screen (trellis_plot_core.dart:855-857, read back
// by getViewport() :887-889); oneColumnOnly (:43) = one column serves both [DOM 2026-08-11].
function axisViewportCount(n: number, oneColumnOnly: boolean): number {
  return Math.min(oneColumnOnly ? 5 : (n < 5 * 1.5 ? n : 5), n);
}

// Row-major grid: a cell's DOM index is yIndex * stride + xIndex. The stride is the RENDERED
// row width — the clamped X viewport — which coincides with xCategoriesCount only while the
// X axis is unclamped, and that is what the min() covers.
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

// Structural row count of an X x Y combination (dataset-level, ignores filter).
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

// Exact-text menu-label locator restricted to labels actually on screen: Playwright's
// :visible is the same predicate a real click enforces, so a hidden item does not match.
function visibleMenuLabel(page: Page, text: string) {
  return page.locator('.d4-menu-item-label:visible').filter({hasText: new RegExp('^\\s*' + text.replace(/[.*+?^${}()|[\]\\/]/g, '\\$&') + '\\s*$')});
}

// Open a context-menu GROUP and click one of its children the way a person does. MENU
// REACHABILITY IS NOT DOM PRESENCE [DOM 2026-08-12]: a collapsed child is in the DOM but not
// on screen. A group opens on HOVER, and the pointer must travel along the header's own row.
async function clickMenuItemInGroup(page: Page, group: string, item: string): Promise<void> {
  const header = visibleMenuLabel(page, group).last();
  await expect(header, `context-menu group "${group}" is not on screen`).toBeVisible({timeout: 10000});
  const box = await header.boundingBox();
  expect(box, `context-menu group "${group}" has no box`).not.toBeNull();
  const y = box!.y + box!.height / 2;
  await page.mouse.move(box!.x + 2, y);
  await page.mouse.move(box!.x + box!.width - 2, y, {steps: 10});
  const child = visibleMenuLabel(page, item).last();
  await expect(child, `"${item}" never became visible after hovering "${group}"`).toBeVisible({timeout: 10000});
  await child.click();
}

// Click a menu item that must be reachable WITHOUT opening any group; same reachability
// rule. When it is unreachable the failure names the collapsed container that hides it.
async function clickTopLevelMenuItem(page: Page, item: string): Promise<void> {
  const target = visibleMenuLabel(page, item).last();
  await target.waitFor({state: 'visible', timeout: 10000}).catch(() => {});
  if (await visibleMenuLabel(page, item).count() === 0) {
    const where = await page.evaluate((item) => {
      const el = Array.from(document.querySelectorAll('.d4-menu-item-label')).find((n) => n.textContent?.trim() === item) as HTMLElement | undefined;
      if (!el) return 'it is not in the menu at all';
      for (let p = el.parentElement; p; p = p.parentElement)
        if (getComputedStyle(p).display === 'none') return `it is collapsed inside "${p.className}" — that group has to be hovered open first`;
      return 'it is in the DOM but has no box';
    }, item);
    throw new Error(`context-menu item "${item}" is not reachable: ${where}`);
  }
  await target.click();
}

// Class-2 selectors for the Legend section's GROK-20432 guard [DOM 2026-08-12]: the Box plot
// tab's prop-show-all-categories is the ONLY reachable read channel for the inner box plot's
// look state — the trellis exposes no getter for its inner-viewer instance.

// Expand the Box plot tab's collapsed "Style" section. While it is collapsed the
// prop-show-all-categories row is zero-size (offsetParent null) and a click on it hangs.
// Idempotent, so every caller can front-load it.
async function ensureStyleExpanded(page: Page): Promise<void> {
  await page.evaluate(() => {
    const target = document.querySelector('.property-grid tr[name="prop-show-all-categories"]') as HTMLElement | null;
    const laidOut = !!target && target.getBoundingClientRect().height > 0 && target.offsetParent !== null;
    if (laidOut) return;
    const header = document.querySelector('.property-grid tr[name="prop-category-style"] td') as HTMLElement | null;
    header?.dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true, view: window}));
  });
  await page.waitForTimeout(600);
}

// Read the inner Box plot's "Show All Categories" off the Box plot tab of the trellis
// Context Panel; null when that row is absent.
async function showAllCategoriesState(page: Page): Promise<boolean | null> {
  await ensureStyleExpanded(page);
  return page.evaluate(() => {
    const tr = document.querySelector('.property-grid tr[name="prop-show-all-categories"]') as HTMLElement | null;
    const box = tr?.querySelector('input[type="checkbox"]') as HTMLInputElement | null;
    return box ? box.checked : null;
  });
}

async function openBoxPlotTab(page: Page): Promise<void> {
  await page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    const panel = root.closest('.panel-base') as HTMLElement | null;
    const gear = panel?.querySelector('.panel-titlebar [name="icon-font-icon-settings"]') as HTMLElement | null;
    gear?.click();
  });
  await page.waitForTimeout(1000);
  await page.evaluate(() => {
    const tab = Array.from(document.querySelectorAll('.d4-tab-header'))
      .find((h) => (h as HTMLElement).innerText.trim() === 'Box plot') as HTMLElement | undefined;
    tab?.click();
  });
  await page.waitForTimeout(800);
}

// Set "Show All Categories" through the real property-grid editor rather than the props
// channel, so the guard exercises the same path a user takes.
async function setShowAllCategories(page: Page, desired: boolean): Promise<boolean | null> {
  await ensureStyleExpanded(page);
  const box = page.locator('.property-grid tr[name="prop-show-all-categories"] input[type="checkbox"]');
  if (await box.count() === 0) return null;
  if (await box.isChecked() !== desired) {
    await box.click();
    await page.waitForTimeout(900);
  }
  return box.isChecked();
}

// Value labels of the inner Box plot's color legend, once its items have rendered. The
// set is configuration-dependent, so it is read from the live DOM and NEVER hardcoded.
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

// Observable state of ONE legend item [DOM 2026-08-12]. `opacity` is the check-state — 1
// shown, 0.5 switched off — and the only field graded on. `active` reports
// .d4-legend-item-current: subset membership, NOT "checked", absent from a pristine legend.
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

// Switch a legend category OFF via its "x" cross, settling on the OBSERVED dimming — fatal
// on purpose, since the mere presence of a cross element would otherwise stand in for the
// uncheck. Never the LAST still-shown category: that click RESETS — use the reset helper.
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
  await page.waitForTimeout(600);
  return {...(await legendItemState(page, category)), clicked: true};
}

// Click a legend category's VALUE LABEL: an EXCLUSIVE SELECT — "show only this one" — and not
// the inverse of the cross [DOM 2026-08-12]. The settle waits for BOTH halves: waiting only on
// the clicked item is what made the retired "re-check" claim look true.
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
  await page.waitForTimeout(600);
  return {...(await legendItemState(page, category)), clicked: true};
}

// Whole-legend snapshot: every entry's opacity and subset-marker class, in DOM order. The
// exclusive-select and reset steps are graded on the WHOLE legend — an assertion that only
// looks at the clicked item cannot tell an exclusive select from a re-check.
async function legendSnapshot(page: Page): Promise<{name: string; opacity: number; current: boolean}[]> {
  return page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
    if (!root) return [];
    // Same filter legendCategories() applies, so a snapshot and a category list are
    // comparable by length.
    return Array.from(root.querySelectorAll('[name="legend"] .d4-legend-item')).map((it) => ({
      name: (it.querySelector('.d4-legend-value')?.textContent ?? '').trim(),
      opacity: parseFloat(getComputedStyle(it).opacity),
      current: it.classList.contains('d4-legend-item-current'),
    })).filter((s) => s.name.length > 0);
  });
}

// Cross-click the LAST still-shown category: this does not hide the final series, it RESETS
// the legend to the pristine "everything shown" state [DOM 2026-08-12]. The section's only
// restore channel — the value label narrows the subset further and never widens it back.
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
  await page.waitForTimeout(600);
  return true;
}

async function trellisCount(page: Page): Promise<number> {
  return page.evaluate(() =>
    grok.shell.tv.viewers.filter((x: any) => x.type === 'Trellis plot').length);
}

// Close the Trellis plot through its dock title-bar close button — the title bar lives on
// the .panel-base ancestor, not inside the viewer root.
async function closeTrellis(page: Page): Promise<void> {
  await page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]');
    const panel = root?.closest('.panel-base');
    const btn = panel?.querySelector('.panel-titlebar [name="Close"]') as HTMLElement | null;
    btn?.click();
  });
}

// Cell geometry of the rendered grid. There is no DOM row/column container, so a row is a
// group of cells sharing a rounded top coordinate and a column a rounded left one.
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

// Category scroll slider: how much of its track the pan handle spans. The sliders sit in the
// DOM even when every category fits, so PRESENCE carries no signal — the EXTENT does. SENTINEL:
// ratio -1 means NOT MEASURED, and every caller must REJECT it before comparing extents.
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

// Category labels are SVG <text> nodes classed .d4-trellis-plot-cat-item-horz (X strip) and
// -vert (Y strip); the same classes land on the invisible hit-test rects, hence the tag
// filter. The rotate() transform — 0 horizontal, -90 vertical — serves the orientation reads.
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

// DOM index of the cell carrying the current-cell marker, row-major and therefore directly
// comparable with cellIndexFor(); -1 when no cell is current. The index is what makes a
// keyboard move DIRECTIONAL: "some cell is current" is already true before any arrow key.
async function currentCellIndex(page: Page): Promise<number> {
  return page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
    if (!root) return -1;
    return Array.from(root.querySelectorAll('.d4-trellis-plot-cell'))
      .findIndex((c) => c.classList.contains('d4-trellis-cell-current'));
  });
}

// Route keyboard input to the trellis: its key handler runs only while the focused element is
// the charts div or sits inside it, and a cell click focuses the inner canvas instead. The
// .d4-trellis-plot-charts-grid class is there only with gridlines, hence the tabindex fallback.
async function focusChartsGrid(page: Page): Promise<void> {
  await page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
    const grid = (root?.querySelector('.d4-trellis-plot-charts-grid') ??
      root?.querySelector('div[tabindex="-1"]')) as HTMLElement | null;
    grid?.focus();
  });
}

// The category combination the trellis reports for the newly current cell, from the last
// d4-trellis-plot-current-cell-changed payload parked on window.__escCc. matchCondition arrives
// as an opaque Dart handle, so it is walked for strings owned by exactly one split column.
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

test('Trellis plot tests', async ({page}) => {
  // Many viewer attaches + layout round-trips + a browser resize need a large budget.
  test.setTimeout(900_000);
  page.setDefaultTimeout(120_000);

  // Console/page-error floor, attached before any section runs (grok.shell.warnings is not
  // exposed by the client). Sections assert a DELTA across their own actuation: the monolith
  // accumulates ambient noise that would make a global zero-error floor false-red.
  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });

  await loginToDatagrok(page);
  await page.waitForTimeout(5000);

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

  // The live DemoFiles demog on dev has 5850 rows, SEX(2), RACE(4). Read from the
  // live dataframe — never hardcoded downstream.
  const setup = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    return {rowCount: df.rowCount, sex: df.col('SEX').categories.length, race: df.col('RACE').categories.length};
  });
  const fullRowCount = setup.rowCount;
  expect(setup).toEqual({rowCount: 5850, sex: 2, race: 4});
  // Cell count of the canonical SEX x RACE grid, computed through the clamp rather than
  // written down as 8: a literal states the product instead of the rule and would keep
  // passing on a build whose category viewport clamp changed.
  const canonicalCellCount = axisViewportCount(setup.sex, false) * axisViewportCount(setup.race, false);

  await v.addViewerByIcon(page, 'trellis-plot', 'Trellis-plot', 15000);

  // The canonical state every section starts from: one Trellis plot on the demog view, Scatter
  // SEX x RACE, neutral axes/onClick, no stray viewers. Sections with risky UI actuation call
  // it from a finally so a thrown assertion cannot cascade into the next section.
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
        await new Promise((r) => setTimeout(r, 1500));
      } catch (_) { /* best-effort restore */ }
    });
  };

  // #### Inner viewer types
  await softStep('Inner viewer types', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      // Prescribed small grid (SEX x RACE) so every populated cell hosts a canvas.
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      await new Promise((res) => setTimeout(res, 1500));
      const r: any[] = [];
      const cellsHaveCanvas = () => {
        const cells = root.querySelectorAll('.d4-trellis-plot-cell');
        let withCanvas = 0;
        for (const c of Array.from(cells)) if (c.querySelector('canvas')) withCanvas++;
        return withCanvas;
      };
      // Switch the inner type through the control-panel selector combo: the props channel does
      // NOT reliably rebuild the heavier types (Summary / Sparklines / PC Plot) and echoes an
      // invalid type while rendering nothing, whereas the selector drives setViewerType.
      const uiSwitch = async (icon: string) => {
        const vs = root.querySelector('[name="viewer selector"]') as HTMLElement;
        vs.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
        await new Promise((res) => setTimeout(res, 600));
        const item = document.querySelector(`.d4-combo-drop-down [name="${icon}"]`);
        (item?.closest('.d4-list-item') as HTMLElement | null)?.click();
        await new Promise((res) => setTimeout(res, 1600));
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
        // Points viewer — the 11th trellisable Dart type. Its descriptor name is 'Heatmap' and
        // it reuses the heat-map icon, so it lists as 'Heatmap' / icon-heat-map in the selector
        // (the separate 'Heat map' grid viewer is not trellisable).
        ['Heatmap', 'icon-heat-map', {columnNames: ['AGE', 'HEIGHT', 'WEIGHT']}],
      ];
      for (const [, icon, look] of types) {
        await uiSwitch(icon);
        try { tp.setOptions({innerViewerLook: look}); } catch {}
        await new Promise((res) => setTimeout(res, 900));
        r.push({type: tp.props.viewerType, cellsWithCanvas: cellsHaveCanvas()});
      }
      return r;
    });
    // Not just a prop echo: every switch must also rebuild the per-cell inner viewers,
    // so each populated cell has to host a canvas.
    const expectedTypes = ['Scatter plot', 'Bar chart', 'Histogram', 'Line chart', 'Box plot',
      'Pie chart', 'Density plot', 'Summary', 'Sparklines', 'PC Plot', 'Heatmap'];
    expect(result.map((x: any) => x.type)).toEqual(expectedTypes);
    for (const x of result) expect(x.cellsWithCanvas).toBeGreaterThan(0);

    // The canvas-presence floor above is weak, so prove the Points viewer's per-cell
    // render RESPONDS to its 'Columns' setting: shrink the column set and require a real
    // ink delta. Ends on Scatter SEX x RACE so downstream sections start canonical.
    const points = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const hash = (i: number) => {
        const cv = root.querySelectorAll('.d4-trellis-plot-cell')[i]?.querySelector('canvas') as HTMLCanvasElement | null;
        if (!cv) return null;
        try {
          const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
          let h = 0;
          for (let k = 0; k < img.length; k += 4) h = (h * 31 + ((img[k] << 16) | (img[k + 1] << 8) | img[k + 2])) % 2147483647;
          return h;
        } catch { return null; }
      };
      tp.setOptions({innerViewerLook: {columnNames: ['AGE', 'HEIGHT', 'WEIGHT']}});
      await new Promise((r) => setTimeout(r, 1400));
      const before = {a: hash(0), b: hash(1)};
      tp.setOptions({innerViewerLook: {columnNames: ['AGE']}});
      await new Promise((r) => setTimeout(r, 1400));
      const after = {a: hash(0), b: hash(1)};
      const applied = tp.getOptions().look.innerViewerLook.columnNames;
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      await new Promise((r) => setTimeout(r, 1400));
      return {before, after, applied};
    });
    expect(points.applied).toEqual(['AGE']);
    // BOTH endpoints of BOTH cells: a vanished canvas hashes to null and would satisfy the
    // inequality without any repaint, so an unguarded `before` end admits exactly the same
    // hole as an unguarded `after` end.
    expect(points.before.a).not.toBeNull();
    expect(points.before.b).not.toBeNull();
    expect(points.after.a).not.toBeNull();
    expect(points.after.b).not.toBeNull();
    expect(points.after.a).not.toBe(points.before.a);
    expect(points.after.b).not.toBe(points.before.b);
  });

  // #### Global scale
  await softStep('Global scale', async () => {
   try {
    // Graded on the shared re-bounding of the inner axes, NOT on reading globalScale back: a
    // set-then-read through the same channel proves nothing. With per-cell WEIGHT/HEIGHT ranges
    // it re-bounds every cell at once, so two cells in different columns must both repaint.
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      tp.setOptions({innerViewerLook: {xColumnName: 'WEIGHT', yColumnName: 'HEIGHT'}});
      tp.props.globalScale = false;
      await new Promise((r) => setTimeout(r, 3000));
    });
    const idxA = await cellIndexFor(page, 'F', 'Caucasian');
    const idxB = await cellIndexFor(page, 'M', 'Caucasian');
    const cellHashes = () => page.evaluate(({idxA, idxB}) => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const hash = (i: number) => {
        const cv = root.querySelectorAll('.d4-trellis-plot-cell')[i]?.querySelector('canvas') as HTMLCanvasElement | null;
        if (!cv) return null;
        try {
          const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
          let h = 0;
          for (let k = 0; k < img.length; k += 4) h = (h * 31 + ((img[k] << 16) | (img[k + 1] << 8) | img[k + 2])) % 2147483647;
          return h;
        } catch { return null; }
      };
      return {a: hash(idxA), b: hash(idxB)};
    }, {idxA, idxB});

    // Driven guard: across an idle settle window of the same length, with nothing actuated,
    // both cells must hash identically — otherwise a viewer that repainted on its own would
    // satisfy every "changed" assertion below by itself.
    const idleBefore = await cellHashes();
    await page.waitForTimeout(1500);
    const idleAfter = await cellHashes();
    expect(idleBefore.a).not.toBeNull();
    expect(idleBefore.b).not.toBeNull();
    expect(idleAfter.a).toBe(idleBefore.a);
    expect(idleAfter.b).toBe(idleBefore.b);

    const flip = async (value: boolean) => {
      await page.evaluate(async (val) => {
        const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
        tp.props.globalScale = val;
        await new Promise((r) => setTimeout(r, 1600));
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
    // Global scale and the WEIGHT/HEIGHT inner look must not leak into the sections below
    // even if a hash read threw.
    await restoreCanonical();
   }
  });

  // #### Axes visibility
  await softStep('Axes visibility', async () => {
   try {
    // PRECONDITION — GLOBAL SCALE ON, and a Scatter inner [DOM 2026-08-12]: Show X/Y Axes is a
    // NO-OP otherwise (scatterplot_meta.dart:14,16 x trellis_plot_core.dart:1021-1029 ->
    // inner_viewer_axes.dart:227-241) — measured 0 sliders in every mode, vs 8/3/5 with it on.
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      // Prescribed Scatter inner + SEX x RACE so inner axes (and their sliders) exist.
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      tp.props.viewerType = 'Scatter plot';
      tp.props.globalScale = true;
      await new Promise((res) => setTimeout(res, 1500));
      // Inner-axis sliders are wrapped in .d4-range-selector — the class-1
      // discriminator vs the two always-present category scroll sliders (bare in
      // .d4-layout-center); count only the inner-axis kind so the toggle is a real delta.
      const innerSliderCount = () => root.querySelectorAll('.d4-range-selector > svg[type="range-slider"]').length;
      // Per-axis counts are the DOM witness for Show X/Y Axes: the axis canvases carry no class
      // of their own, but an inner-axis slider exists only while its axis is displayed. Laid-out
      // only — the wrappers are visibility:hidden until hovered, which still leaves a box.
      const axisSliders = (ax: string) =>
        Array.from(root.querySelectorAll(`.d4-range-selector > svg[type="range-slider"][name="${ax}-slider"]`))
          .filter((el) => {
            const b = el.getBoundingClientRect();
            return b.width > 0 && b.height > 0;
          }).length;
      tp.props.showRangeSliders = true;
      tp.props.showYAxes = 'Always';
      await new Promise((res) => setTimeout(res, 1200));
      const r: any[] = [];
      const xByMode: Record<string, number> = {};
      for (const val of ['Always', 'Never', 'Auto']) {
        tp.props.showXAxes = val;
        await new Promise((res) => setTimeout(res, 1200));
        r.push(tp.props.showXAxes);
        xByMode[val] = axisSliders('x');
      }
      tp.props.showXAxes = 'Always';
      await new Promise((res) => setTimeout(res, 1200));
      const yByMode: Record<string, number> = {};
      for (const val of ['Always', 'Never', 'Auto']) {
        tp.props.showYAxes = val;
        await new Promise((res) => setTimeout(res, 1200));
        r.push(tp.props.showYAxes);
        yByMode[val] = axisSliders('y');
      }
      tp.props.showXAxes = 'Always';
      tp.props.showYAxes = 'Always';
      tp.props.showRangeSliders = false;
      await new Promise((res) => setTimeout(res, 1200));
      const slidersOff = innerSliderCount();
      tp.props.showRangeSliders = true;
      await new Promise((res) => setTimeout(res, 1200));
      const slidersOn = innerSliderCount();
      return {modes: r, xByMode, yByMode, slidersOff, slidersOn};
    });
    expect(result.modes).toEqual(['Always', 'Never', 'Auto', 'Always', 'Never', 'Auto']);
    // Always vs Never is graded on the axis's own slider population, not on the mode string read
    // back. Never is EXACTLY ZERO, not merely fewer — measured 0 x-sliders at X = Never and 0 y
    // at Y = Never [DOM 2026-08-12] — which catches a partial hide.
    expect(result.xByMode.Always).toBeGreaterThan(0);
    expect(result.xByMode.Never).toBe(0);
    expect(result.yByMode.Always).toBeGreaterThan(0);
    expect(result.yByMode.Never).toBe(0);
    // AUTO IS NOT A THIRD STATE [DOM 2026-08-12]: trellis_plot_core.dart:1021-1029 gates the axis
    // on "showXAxes != Never", so Auto and Always take the same branch and measure identical —
    // graded as distinguishable from Never, indistinguishable from Always.
    expect(result.xByMode.Auto).toBe(result.xByMode.Always);
    expect(result.yByMode.Auto).toBe(result.yByMode.Always);
    // Show Range Sliders creates/destroys the inner-axis slider DOM. OFF IS EXACTLY ZERO, not
    // merely fewer — measured 0 wrappers / 0 x / 0 y [DOM 2026-08-12] — which catches a partial
    // teardown; the ON reading is the anti-vacuity control that sliders were ever there.
    expect(result.slidersOff).toBe(0);
    expect(result.slidersOn).toBeGreaterThan(0);
   } finally {
    // restoreCanonical() puts globalScale and both axis modes back but does not touch
    // showRangeSliders, so that one is reset by hand.
    await page.evaluate(async () => {
      try {
        const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
        tp.props.showRangeSliders = true;
        await new Promise((res) => setTimeout(res, 400));
      } catch (_) { /* best-effort restore */ }
    });
    await restoreCanonical();
   }
  });

  // #### Range sliders with global scale
  await softStep('Range sliders with global scale', async () => {
   try {
    // Prime the viewer so an inner-axis slider (the .d4-range-selector kind) exists.
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
      await new Promise((r) => setTimeout(r, 2200));
    });

    // Step 4: the inner-axis slider is revealed by hovering its axis area — measure the
    // box only after the hover.
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

    // Step 5: a TRUSTED CDP drag from the max end of the track narrows the shared axis;
    // under global scale that re-bounds ALL cells at once, so two cells must repaint.
    const idxA = await cellIndexFor(page, 'F', 'Caucasian');
    const idxB = await cellIndexFor(page, 'M', 'Caucasian');
    const before = await page.evaluate(({idxA, idxB}) => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const hash = (i: number) => {
        const cv = root.querySelectorAll('.d4-trellis-plot-cell')[i]?.querySelector('canvas') as HTMLCanvasElement | null;
        if (!cv) return null;
        try {
          const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
          let h = 0;
          for (let i = 0; i < img.length; i += 4) h = (h * 31 + ((img[i] << 16) | (img[i + 1] << 8) | img[i + 2])) % 2147483647;
          return h;
        } catch { return null; }
      };
      return {a: hash(idxA), b: hash(idxB)};
    }, {idxA, idxB});
    await page.mouse.move(sliderBox!.x + sliderBox!.w - 4, sliderBox!.y + sliderBox!.h / 2);
    await page.mouse.down();
    await page.mouse.move(sliderBox!.x + sliderBox!.w * 0.45, sliderBox!.y + sliderBox!.h / 2, {steps: 12});
    await page.mouse.up();
    await page.waitForTimeout(1500);
    const after = await page.evaluate(({idxA, idxB}) => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const hash = (i: number) => {
        const cv = root.querySelectorAll('.d4-trellis-plot-cell')[i]?.querySelector('canvas') as HTMLCanvasElement | null;
        if (!cv) return null;
        try {
          const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
          let h = 0;
          for (let i = 0; i < img.length; i += 4) h = (h * 31 + ((img[i] << 16) | (img[i + 1] << 8) | img[i + 2])) % 2147483647;
          return h;
        } catch { return null; }
      };
      return {a: hash(idxA), b: hash(idxB)};
    }, {idxA, idxB});
    // `before.b` is also the baseline step 6 restores back to, so it must be a real hash —
    // otherwise that equality would be null against null.
    expect(before.a).not.toBeNull();
    expect(before.b).not.toBeNull();
    expect(after.a).not.toBeNull();
    expect(after.b).not.toBeNull();
    expect(after.a).not.toBe(before.a);
    expect(after.b).not.toBe(before.b);

    // Step 6: 'Reset Inner Range Sliders' must restore every cell to the pre-drag hash. The item
    // is conditional on the inner viewer having axis sliders — absent with a Pie chart inner
    // [DOM 2026-08-11] — so the Scatter inner primed above is a precondition, not decoration.
    await page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell').first().click({button: 'right', position: {x: 6, y: 6}});
    await clickTopLevelMenuItem(page, 'Reset Inner Range Sliders');
    const reset = await page.evaluate(async ({idxA, idxB, beforeA, beforeB}) => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      await new Promise((res) => setTimeout(res, 1200));
      const hash = (i: number) => {
        const cv = root.querySelectorAll('.d4-trellis-plot-cell')[i]?.querySelector('canvas') as HTMLCanvasElement | null;
        if (!cv) return null;
        try {
          const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
          let h = 0;
          for (let i = 0; i < img.length; i += 4) h = (h * 31 + ((img[i] << 16) | (img[i + 1] << 8) | img[i + 2])) % 2147483647;
          return h;
        } catch { return null; }
      };
      return {restoredA: hash(idxA) === beforeA, restoredB: hash(idxB) === beforeB};
    }, {idxA, idxB, beforeA: before.a, beforeB: before.b});
    expect(reset.restoredA).toBe(true);
    expect(reset.restoredB).toBe(true);

    // Step 7: the Y-axis mirror of the X flow above (name="y-slider").
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.showYAxes = 'Always';
      await new Promise((r) => setTimeout(r, 1200));
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
    const yBefore = await page.evaluate(({idxA, idxB}) => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const hash = (i: number) => {
        const cv = root.querySelectorAll('.d4-trellis-plot-cell')[i]?.querySelector('canvas') as HTMLCanvasElement | null;
        if (!cv) return null;
        try {
          const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
          let h = 0;
          for (let i = 0; i < img.length; i += 4) h = (h * 31 + ((img[i] << 16) | (img[i + 1] << 8) | img[i + 2])) % 2147483647;
          return h;
        } catch { return null; }
      };
      return {a: hash(idxA), b: hash(idxB)};
    }, {idxA, idxB});
    await page.mouse.move(ySliderBox!.x + ySliderBox!.w / 2, ySliderBox!.y + ySliderBox!.h - 4);
    await page.mouse.down();
    await page.mouse.move(ySliderBox!.x + ySliderBox!.w / 2, ySliderBox!.y + ySliderBox!.h * 0.45, {steps: 12});
    await page.mouse.up();
    await page.waitForTimeout(1500);
    const yAfter = await page.evaluate(({idxA, idxB}) => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const hash = (i: number) => {
        const cv = root.querySelectorAll('.d4-trellis-plot-cell')[i]?.querySelector('canvas') as HTMLCanvasElement | null;
        if (!cv) return null;
        try {
          const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
          let h = 0;
          for (let i = 0; i < img.length; i += 4) h = (h * 31 + ((img[i] << 16) | (img[i + 1] << 8) | img[i + 2])) % 2147483647;
          return h;
        } catch { return null; }
      };
      return {a: hash(idxA), b: hash(idxB)};
    }, {idxA, idxB});
    expect(yBefore.a).not.toBeNull();
    expect(yBefore.b).not.toBeNull();
    expect(yAfter.a).not.toBeNull();
    expect(yAfter.b).not.toBeNull();
    expect(yAfter.a).not.toBe(yBefore.a);
    expect(yAfter.b).not.toBe(yBefore.b);
   } finally {
    // A failure must not leave global scale on for the following sections.
    await restoreCanonical();
   }
  });

  // #### Gridlines
  await softStep('Gridlines', async () => {
    // The structural signal is the d4-trellis-plot-charts-grid class on the charts
    // host, not the property read-back. In 'auto' the class follows the inner type:
    // gridlines only for the Scatter plot.
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const hasGrid = () => !!root.querySelector('.d4-trellis-plot-charts-grid');
      const r: {mode: string; prop: string; grid: boolean}[] = [];
      for (const val of ['always', 'never', 'auto']) {
        tp.props.showGridlines = val;
        await new Promise((res) => setTimeout(res, 600));
        r.push({mode: val, prop: tp.props.showGridlines, grid: hasGrid()});
      }
      tp.props.viewerType = 'Bar chart';
      await new Promise((res) => setTimeout(res, 1200));
      const autoNonScatter = hasGrid();
      tp.props.viewerType = 'Scatter plot';
      await new Promise((res) => setTimeout(res, 1200));
      return {r, autoNonScatter, autoScatter: hasGrid()};
    });
    expect(result.r.map((x) => x.prop)).toEqual(['always', 'never', 'auto']);
    expect(result.r.map((x) => x.grid)).toEqual([true, false, true]);
    expect(result.autoNonScatter).toBe(false);
    expect(result.autoScatter).toBe(true);
  });

  // #### Tiles mode
  await softStep('Tiles mode', async () => {
   try {
    // Single X column so the tiled view is meaningful. Tiles Per Row is graded on the
    // GEOMETRY it re-flows — how many cells share a row — not on reading tilesPerRow back
    // out of the property it was written to.
    const setup = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.globalScale = false;
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['RACE'];
      tp.props.yColumnNames = [];
      await new Promise((r) => setTimeout(r, 1500));
      tp.props.useTiledView = true;
      await new Promise((r) => setTimeout(r, 1800));
      return {raceCats: grok.shell.tv.dataFrame.col('RACE').categories.length,
        tilesPerRow: tp.props.tilesPerRow};
    });
    // RACE must carry between three and six categories: below three the 2-per-row and 6-per-row
    // layouts coincide and nothing is measured; above six the 6-per-row step's premise — more
    // tiles per row than there are tiles — stops being true.
    expect(setup.raceCats).toBeGreaterThanOrEqual(3);
    expect(setup.raceCats).toBeLessThanOrEqual(6);
    // Tiled single-column geometry: row width min(tilesPerRow, N) (trellis_plot_core.dart:81-88),
    // row count ceil(N / rowWidth) (:89-95), each under the per-axis clamp. The grid is a PADDED
    // rectangle, so "one cell per category" holds only when the row width divides N [DOM 2026-08-11].
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
        await new Promise((r) => setTimeout(r, 1600));
      }, n);
      return gridGeometry(page);
    };

    // Tiles Per Row = 2 -> no row wider than two tiles, and more than one row.
    const atTwo = await setTiles(2);
    expect({count: atTwo.count, rows: atTwo.rows, maxPerRow: atTwo.maxPerRow})
      .toEqual(tiledGeometry(setup.raceCats, 2));
    expect(atTwo.maxPerRow).toBeLessThanOrEqual(2);
    expect(atTwo.rows).toBeGreaterThan(1);

    // Tiles Per Row = 6 exceeds the tile count, so every tile collapses into one row. The width
    // saturates at min(6, N) and is then capped at 5 by the clamp, so "six per row" differs from
    // "five per row" only while N stays at or below five — which the computed expectation carries.
    const atSix = await setTiles(6);
    expect({count: atSix.count, rows: atSix.rows, maxPerRow: atSix.maxPerRow})
      .toEqual(tiledGeometry(setup.raceCats, 6));
    expect(atSix.rows).toBe(1);

    // Turning Tiled View off has to change the geometry, so it is measured from the 2-per-row
    // layout rather than the 6-per-row one, which is already a single row. The direction of the
    // untiled strip is not pinned here — that is the focused tiled scenario's contract.
    const backToTwo = await setTiles(2);
    expect(backToTwo.rows).toBeGreaterThan(1);
    const untiled = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.useTiledView = false;
      await new Promise((r) => setTimeout(r, 1800));
      return tp.props.useTiledView;
    });
    const offGeometry = await gridGeometry(page);
    expect(untiled).toBe(false);
    // Untiled, a single split column lays its categories out along ONE axis, so the strip
    // holds the clamped viewport of that axis — five cells at most, whichever direction it
    // runs in. A 6-category column would render 5, not 6.
    expect(offGeometry.count).toBe(axisViewportCount(setup.raceCats, true));
    expect({rows: offGeometry.rows, maxPerRow: offGeometry.maxPerRow})
      .not.toEqual({rows: backToTwo.rows, maxPerRow: backToTwo.maxPerRow});
   } finally {
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      try {
        tp.props.useTiledView = true;
        tp.props.tilesPerRow = 4;
        await new Promise((r) => setTimeout(r, 600));
      } catch (_) { /* best-effort restore */ }
    });
    await restoreCanonical();
   }
  });

  // #### Category management
  await softStep('Category management', async () => {
    // Adding/removing an X category column changes the DOM cell count (product of the
    // split cardinalities under the viewport clamp), not just the xColumnNames prop.
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const cells = () => root.querySelectorAll('.d4-trellis-plot-cell').length;
      const r: any[] = [];

      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      await new Promise((res) => setTimeout(res, 1500));
      r.push({x: [...tp.props.xColumnNames], cells: cells()});

      tp.props.xColumnNames = ['SEX', 'DIS_POP'];
      await new Promise((res) => setTimeout(res, 1500));
      r.push({x: [...tp.props.xColumnNames], cells: cells()});

      tp.props.xColumnNames = ['SEX'];
      await new Promise((res) => setTimeout(res, 1500));
      r.push({x: [...tp.props.xColumnNames], cells: cells()});

      return r;
    });
    // The grid is the product of the two axes' CLAMPED viewports, not of their raw cardinalities:
    // adding DIS_POP takes X to 12 raw slots, of which the clamp renders 5. Both states are pinned
    // exactly from the live cardinalities — "more cells than before" would pass on any growth.
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

    // Steps 4-6: the label toggles are graded on the rendered category labels, not on the
    // properties they were written to — turning a toggle off has to empty its strip and
    // turning it back on has to repopulate it.
    const labelsOn = await categoryLabels(page);
    expect(labelsOn.x.length).toBeGreaterThan(0);
    expect(labelsOn.y.length).toBeGreaterThan(0);

    const setLabels = async (x: boolean, y: boolean) => {
      await page.evaluate(async (v) => {
        const tp = Array.from(grok.shell.tv.viewers).find((vw: any) => vw.type === 'Trellis plot') as any;
        tp.props.showXLabels = v.x;
        tp.props.showYLabels = v.y;
        await new Promise((r) => setTimeout(r, 1600));
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

  // #### Pack categories
  await softStep('Pack categories', async () => {
    // Pack Categories on: a filter that empties a RACE category drops its whole row from the
    // grid, a real cell-count delta; off: the empty row reappears at the base count.
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const df = grok.shell.tv.dataFrame;
      const cells = () => root.querySelectorAll('.d4-trellis-plot-cell').length;
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      tp.props.packCategories = true;
      await new Promise((r) => setTimeout(r, 1600));
      const baseCells = cells();

      grok.shell.tv.getFiltersGroup();
      await new Promise((r) => setTimeout(r, 1500));
      const fg = grok.shell.tv.getFiltersGroup();
      const cats = df.col('RACE').categories;
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: cats.filter((c: string) => c !== 'Asian')});
      await new Promise((r) => setTimeout(r, 1500));
      const packedOnCells = cells();

      tp.props.packCategories = false;
      await new Promise((r) => setTimeout(r, 1200));
      const packedOffCells = cells();

      tp.props.packCategories = true;
      await new Promise((r) => setTimeout(r, 1200));
      const packedOnAgainCells = cells();

      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: cats});
      await new Promise((r) => setTimeout(r, 800));
      return {baseCells, packedOnCells, packedOffCells, packedOnAgainCells};
    });
    expect(result.packedOnCells).toBeLessThan(result.baseCells);
    expect(result.packedOffCells).toBe(result.baseCells);
    expect(result.packedOnAgainCells).toBe(result.packedOnCells);
  });

  // #### On Click functionality
  await softStep('On Click functionality', async () => {
    const cellLocator = page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell');
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      df.filter.setAll(true); df.selection.setAll(false); df.rows.requestFilter();
      await new Promise((r) => setTimeout(r, 1500));
    });
    await expect(cellLocator).toHaveCount(canonicalCellCount);

    // Step 2: On Click = Select, corner-click F|Caucasian -> selection == that cell's rows.
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.onClick = 'Select';
      await new Promise((r) => setTimeout(r, 600));
    });
    const selExpected = await comboRowCount(page, 'SEX', 'F', 'RACE', 'Caucasian');
    let idx = await cellIndexFor(page, 'F', 'Caucasian');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(900);
    const selAfterClick = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(selAfterClick).toBe(selExpected);

    // Step 3: changing the inner viewer type does NOT alter the current selection.
    const selAfterTypeChange = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.viewerType = 'Bar chart';
      tp.setOptions({innerViewerLook: {splitColumnName: 'RACE', valueColumnName: 'AGE'}});
      await new Promise((r) => setTimeout(r, 1200));
      const sel = grok.shell.tv.dataFrame.selection.trueCount;
      tp.props.viewerType = 'Scatter plot';
      await new Promise((r) => setTimeout(r, 1000));
      return sel;
    });
    expect(selAfterTypeChange).toBe(selExpected);

    // Step 4: clicking another non-empty cell changes the selection to that cell's rows.
    const mBlackExpected = await comboRowCount(page, 'SEX', 'M', 'RACE', 'Black');
    idx = await cellIndexFor(page, 'M', 'Black');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(900);
    const selAfterOther = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(selAfterOther).toBe(mBlackExpected);
    expect(selAfterOther).not.toBe(selExpected);

    // Step 5: changing an axis column does NOT change the current selection.
    const selAfterAxis = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.yColumnNames = ['SEVERITY'];
      await new Promise((r) => setTimeout(r, 1200));
      const sel = grok.shell.tv.dataFrame.selection.trueCount;
      tp.props.yColumnNames = ['RACE'];
      await new Promise((r) => setTimeout(r, 1000));
      return sel;
    });
    expect(selAfterAxis).toBe(mBlackExpected);

    // Step 7: On Click = Filter, corner-click F|Caucasian -> df.filter drops to that cell's rows.
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true); df.selection.setAll(false); df.rows.requestFilter();
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.onClick = 'Filter';
      await new Promise((r) => setTimeout(r, 700));
    });
    const filterExpected = await comboRowCount(page, 'SEX', 'F', 'RACE', 'Caucasian');
    idx = await cellIndexFor(page, 'F', 'Caucasian');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(900);
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

    // Step 8: changing the inner viewer does NOT alter the active trellis filter.
    const filterAfterType = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.viewerType = 'Pie chart';
      tp.setOptions({innerViewerLook: {categoryColumnName: 'RACE'}});
      await new Promise((r) => setTimeout(r, 1200));
      const c = grok.shell.tv.dataFrame.filter.trueCount;
      tp.props.viewerType = 'Scatter plot';
      await new Promise((r) => setTimeout(r, 1000));
      return c;
    });
    expect(filterAfterType).toBe(filterExpected);

    // Step 9: clicking another non-empty cell updates the filter to that cell's categories.
    const mBlackFilter = await comboRowCount(page, 'SEX', 'M', 'RACE', 'Black');
    idx = await cellIndexFor(page, 'M', 'Black');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(900);
    const filterAfterOther = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const filters: string[] = [];
      for (const f of df.rows.filters) filters.push(f);
      return {count: df.filter.trueCount, filters};
    });
    expect(filterAfterOther.count).toBe(mBlackFilter);
    expect(filterAfterOther.filters).toContain('SEX: M');
    expect(filterAfterOther.filters).toContain('RACE: Black');

    // Step 10: changing an axis column resets the trellis filter to the unfiltered value.
    const filterAfterAxis = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['CONTROL'];
      await new Promise((r) => setTimeout(r, 1500));
      const c = grok.shell.tv.dataFrame.filter.trueCount;
      tp.props.xColumnNames = ['SEX'];
      await new Promise((r) => setTimeout(r, 1500));
      return c;
    });
    expect(filterAfterAxis).toBe(fullRowCount);

    // Step 12, filter half: ESC resets filtering. The PRE-ESC reading makes the reset
    // attributable — a click that silently did nothing leaves the count at full and sails through
    // the post-ESC assertion. The trellis's own df.rows.filters entries are graded alongside.
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true); df.selection.setAll(false); df.rows.requestFilter();
      await new Promise((r) => setTimeout(r, 500));
    });
    idx = await cellIndexFor(page, 'F', 'Caucasian');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(700);
    const beforeEsc = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const filters: string[] = [];
      for (const f of df.rows.filters) filters.push(f);
      return {count: df.filter.trueCount, filters};
    });
    expect(beforeEsc.count).toBeLessThan(fullRowCount);
    expect(beforeEsc.filters.length).toBeGreaterThan(0);
    await focusChartsGrid(page);
    await page.keyboard.press('Escape');
    await page.waitForTimeout(900);
    const afterEsc = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const filters: string[] = [];
      for (const f of df.rows.filters) filters.push(f);
      return {count: df.filter.trueCount, filters};
    });
    expect(afterEsc.count).toBe(fullRowCount);
    expect(afterEsc.filters).toEqual([]);

    // Step 12, selection half — checked under On Click = Select, the only mode where a cell
    // click selects. Under Filter the click selects nothing, so a post-ESC "nothing is
    // selected" cannot tell a cleared selection from one that was never made.
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.onClick = 'Select';
      df.filter.setAll(true); df.selection.setAll(false); df.rows.requestFilter();
      (window as any).__escCc = [];
      (window as any).__escCcSub = tp.onEvent('d4-trellis-plot-current-cell-changed').subscribe((a: any) => {
        const mc = (a && a.args && a.args.matchCondition) ? a.args.matchCondition : (a && a.matchCondition ? a.matchCondition : a);
        (window as any).__escCc.push(mc);
      });
      await new Promise((r) => setTimeout(r, 700));
    });
    idx = await cellIndexFor(page, 'F', 'Caucasian');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(900);
    // Witness the click BEFORE the ESC, or a silent no-op passes the post-ESC assertion
    // untouched: the current-cell marker on the clicked cell and exactly that
    // combination's rows — a bare "something is selected" would pass on the wrong cell.
    expect(await currentCellIndex(page)).toBe(idx);
    const selComboRows = await comboRowCount(page, 'SEX', 'F', 'RACE', 'Caucasian');
    expect(selComboRows).toBeGreaterThan(0);
    const selBeforeEsc = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(selBeforeEsc).toBe(selComboRows);
    await focusChartsGrid(page);
    await page.keyboard.press('Escape');
    await page.waitForTimeout(900);
    const selAfterEsc = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(selAfterEsc).toBe(0);
    // Second product-side confirmation of WHICH combination was selected. Asserted AFTER
    // the ESC contract on purpose: the payload crosses the Dart/JS boundary and an interop
    // shape change must not be able to mask the contract this step exists for.
    expect(await comboFromLastCellChange(page)).toEqual({SEX: 'F', RACE: 'Caucasian'});

    // Back to Filter (and nothing selected) for Step 14, whose AND-composition needs the
    // cell click to filter again.
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      (window as any).__escCcSub?.unsubscribe?.();
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.onClick = 'Filter';
      df.filter.setAll(true); df.selection.setAll(false); df.rows.requestFilter();
      await new Promise((r) => setTimeout(r, 700));
    });

    // Step 14: a Filter-Panel filter ANDs with a trellis cell filter (strictly below
    // the panel-only value). Use a persistent onRowsFiltering subscriber as neutral
    // setup (the real AND-composition is produced by the corner-click below).
    const panelOnly = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true); df.rows.requestFilter();
      const ageCol = df.col('AGE');
      (window as any).__panelSub = df.onRowsFiltering.subscribe(() => {
        for (let i = 0; i < df.rowCount; i++) if (ageCol.get(i) < 40) df.filter.set(i, false, false);
      });
      df.rows.requestFilter();
      await new Promise((r) => setTimeout(r, 800));
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
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(900);
    const bothActive = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(bothActive).toBe(mAsianAge);
    expect(bothActive).toBeLessThan(panelOnly);

    // Step 15 ("Right-click > Reset view") is intentionally NOT implemented: the live trellis
    // context menu has no such item (refdoc Context Menu), and substituting onClick='None' +
    // setAll(true) would prove nothing about a menu path that does not exist.

    // Step 16: On Click = None -> clicking any cell does not change filtering or selection.
    const noneResult = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      (window as any).__panelSub?.unsubscribe?.();
      df.filter.setAll(true); df.selection.setAll(false); df.rows.requestFilter();
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.onClick = 'None';
      await new Promise((r) => setTimeout(r, 700));
      return {filterBefore: df.filter.trueCount, selBefore: df.selection.trueCount};
    });
    idx = await cellIndexFor(page, 'F', 'Caucasian');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(900);
    const noneAfter = await page.evaluate(() => ({
      filterAfter: grok.shell.tv.dataFrame.filter.trueCount,
      selAfter: grok.shell.tv.dataFrame.selection.trueCount,
    }));
    expect(noneAfter.filterAfter).toBe(noneResult.filterBefore);
    expect(noneAfter.selAfter).toBe(noneResult.selBefore);
  });

  // #### Selectors
  await softStep('Selectors', async () => {
   try {
    // showXSelectors/showYSelectors/showControlPanel govern DOM visibility of the
    // selector hosts; assert the viewer-selector's DOM visibility, not just the prop.
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      // Clean SEX x RACE Scatter grid so the control panel + selectors are present.
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      tp.props.showControlPanel = true;
      await new Promise((r) => setTimeout(r, 1500));
      const vsVisible = () => {
        const el = root.querySelector('[name="viewer selector"]') as HTMLElement | null;
        if (!el) return false;
        const b = el.getBoundingClientRect();
        return b.width > 0 && b.height > 0;
      };
      tp.props.showXSelectors = false;
      tp.props.showYSelectors = false;
      tp.props.showControlPanel = false;
      await new Promise((r) => setTimeout(r, 1200));
      const off = {props: {x: tp.props.showXSelectors, y: tp.props.showYSelectors, cp: tp.props.showControlPanel}, vs: vsVisible()};
      tp.props.showXSelectors = true;
      tp.props.showYSelectors = true;
      tp.props.showControlPanel = true;
      await new Promise((r) => setTimeout(r, 1200));
      const on = {props: {x: tp.props.showXSelectors, y: tp.props.showYSelectors, cp: tp.props.showControlPanel}, vs: vsVisible()};
      return {off, on};
    });
    expect(result.off.props).toEqual({x: false, y: false, cp: false});
    expect(result.on.props).toEqual({x: true, y: true, cp: true});
    expect(result.off.vs).toBe(false);
    expect(result.on.vs).toBe(true);

    // Step 5: the CONJUNCTION int.selectors-labels-visibility-coupling claims — an explicit
    // selector-off outlives the reversible Auto Layout hide. VISIBILITY, NOT BOX [DOM 2026-08-12]:
    // hiding sets visibility:hidden on the host, whose box stays — a box count reads 2 in every state.
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
      // The exact attribute value excludes [name="div-column-combobox-category"], the
      // legend's category picker, which mirrors the split column's label (refdoc TRAP).
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
      tp.props.autoLayout = true;
      await new Promise((r) => setTimeout(r, 1200));
    });
    // Exactly two pickers in the SEX x RACE configuration — one X, one Y — pinned rather
    // than "> 0" so the step below reads as a ladder and not as a direction.
    const bothSelectors = await columnSelectorCount();
    expect(bothSelectors).toBe(2);

    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.showXSelectors = false;
      await new Promise((r) => setTimeout(r, 1500));
    });
    const xTurnedOff = await columnSelectorCount();
    expect(xTurnedOff).toBe(1);

    await page.setViewportSize({width: 500, height: 400});
    await page.waitForTimeout(2000);
    const shrunkPanel = await controlPanelVisible();
    const shrunkSelectors = await columnSelectorCount();
    await page.setViewportSize({width: 1920, height: 1080});
    await page.waitForTimeout(2000);
    const restoredPanel = await controlPanelVisible();
    const afterCycle = await columnSelectorCount();
    // The auto-hide really fired and really reverted — on the control panel and on the
    // remaining Y picker alike ([DOM 2026-08-12]: at 500x400 both hosts go
    // visibility:hidden).
    expect(shrunkPanel).toBe(false);
    expect(shrunkSelectors).toBe(0);
    expect(restoredPanel).toBe(true);
    // ...and the explicit off survived it: the restore brings back only what the user did
    // not switch off by hand.
    expect(afterCycle).toBe(xTurnedOff);

    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.showXSelectors = true;
      await new Promise((r) => setTimeout(r, 1200));
    });
    expect(await columnSelectorCount()).toBe(bothSelectors);
   } finally {
    // A throw between the two resizes would hand every later section a 500x400 viewport,
    // with the explicit selector-off leaking along with it.
    await page.setViewportSize({width: 1920, height: 1080});
    await page.waitForTimeout(1000);
    await page.evaluate(async () => {
      try {
        const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
        tp.props.showXSelectors = true;
        tp.props.showYSelectors = true;
        tp.props.showControlPanel = true;
        await new Promise((r) => setTimeout(r, 500));
      } catch (_) { /* best-effort restore */ }
    });
   }
  });

  // #### Allow viewer full screen
  await softStep('Allow viewer full screen', async () => {
   try {
    // THE ICON IS A STATE OF THE HOVER, NOT A CONSTANT [DOM 2026-08-12]: none exists before a
    // hover, one is created and re-parented between cells, and mouseleave removes it
    // (trellis_plot_core.dart:1649-1710) — so the real pointer is parked outside first.
    await page.mouse.move(5, 5);
    await page.waitForTimeout(800);
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      // Wide SEX x RACE cells (the full-screen icon is not rendered on narrow cells).
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      tp.props.allowViewerFullScreen = true;
      await new Promise((r) => setTimeout(r, 1500));

      // Every read re-queries the viewer root and its cells: the full-screen modal
      // re-renders the grid, so a node list captured up front would go stale and every
      // later dispatch would land on detached nodes.
      const rootEl = () => document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
      const cells = () => Array.from(rootEl()?.querySelectorAll('.d4-trellis-plot-cell') ?? []);
      const iconState = () => {
        const list = Array.from(rootEl()?.querySelectorAll('.d4-viewer-icon[name="icon-expand-arrows"]') ?? []) as HTMLElement[];
        if (list.length !== 1) return {count: list.length, parent: -1, top: -1, left: -1};
        const b = list[0].getBoundingClientRect();
        return {count: 1, parent: cells().indexOf(list[0].parentElement as Element),
          top: Math.round(b.top), left: Math.round(b.left)};
      };
      const pointer = async (idx: number, into: boolean) => {
        const cell = cells()[idx];
        if (!cell) return;
        const r = cell.getBoundingClientRect();
        const o = {bubbles: true, cancelable: true, clientX: r.left + r.width / 2, clientY: r.top + r.height / 2};
        for (const t of into ? ['mouseenter', 'mouseover', 'mousemove'] : ['mouseout', 'mouseleave'])
          cell.dispatchEvent(new MouseEvent(t, o));
        await new Promise((res) => setTimeout(res, 800));
      };

      const all = cells();
      const firstIdx = all.findIndex((c) => !!c.querySelector('canvas'));
      const firstTop = firstIdx >= 0 ? Math.round(all[firstIdx].getBoundingClientRect().top) : -1;
      // A cell in another ROW, so the move is visible in both coordinates.
      const secondIdx = all.findIndex((c, i) => i !== firstIdx && !!c.querySelector('canvas') &&
        Math.round(c.getBoundingClientRect().top) !== firstTop);
      const out: any = {cellsFound: false, firstIdx, secondIdx,
        beforeHover: iconState().count, onFirst: {count: -1, parent: -1, top: -1, left: -1},
        onSecond: {count: -1, parent: -1, top: -1, left: -1}, afterLeave: -1,
        modalOpened: false, modalTitle: null as string | null, iconWhenOff: true};
      if (firstIdx < 0 || secondIdx < 0) return out;
      out.cellsFound = true;

      // Step 1: hover the first cell — one icon, inside that very cell.
      await pointer(firstIdx, true);
      out.onFirst = iconState();

      // Step 2: leave it and hover a cell in another row — still one icon, now in that
      // cell and at different coordinates.
      await pointer(firstIdx, false);
      await pointer(secondIdx, true);
      out.onSecond = iconState();

      // Step 3: the pointer leaves and the icon goes with it.
      await pointer(secondIdx, false);
      out.afterLeave = iconState().count;

      // Step 4: hover again and click the icon -> the inner viewer opens full screen.
      await pointer(firstIdx, true);
      const icon = rootEl()?.querySelector('.d4-viewer-icon[name="icon-expand-arrows"]') as HTMLElement | null;
      if (icon) {
        const ib = icon.getBoundingClientRect();
        const io = {bubbles: true, cancelable: true, button: 0, clientX: ib.left + ib.width / 2, clientY: ib.top + ib.height / 2};
        icon.dispatchEvent(new MouseEvent('mousedown', io));
        icon.dispatchEvent(new MouseEvent('mouseup', io));
        icon.dispatchEvent(new MouseEvent('click', io));
        await new Promise((res) => setTimeout(res, 1200));
        const dlg = document.querySelector('.d4-dialog');
        out.modalOpened = !!dlg;
        if (dlg) {
          out.modalTitle = (dlg.querySelector('.d4-dialog-title') as HTMLElement | null)?.textContent?.trim() ?? null;
          const cancel = (dlg.querySelector('[name="button-CANCEL"]') as HTMLElement | null)
            ?? (dlg.querySelector('.d4-dialog-header [name="icon-times"]') as HTMLElement | null);
          if (cancel) cancel.click(); else document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
          await new Promise((res) => setTimeout(res, 500));
        }
      }
      await pointer(firstIdx, false);

      // Step 6: with the setting off, no icon appears on hover.
      tp.props.allowViewerFullScreen = false;
      await new Promise((res) => setTimeout(res, 800));
      await pointer(firstIdx, true);
      out.iconWhenOff = iconState().count > 0;
      await pointer(firstIdx, false);
      tp.props.allowViewerFullScreen = true;
      return out;
    });
    expect(result.cellsFound).toBe(true);
    // Baseline: nothing on screen before the first dispatch, so the reads below grade the
    // hover and not the mere existence of an icon.
    expect(result.beforeHover).toBe(0);
    expect(result.onFirst.count).toBe(1);
    expect(result.onFirst.parent).toBe(result.firstIdx);
    expect(result.onSecond.count).toBe(1);
    expect(result.onSecond.parent).toBe(result.secondIdx);
    expect(result.onSecond.top === result.onFirst.top && result.onSecond.left === result.onFirst.left).toBe(false);
    expect(result.afterLeave).toBe(0);
    expect(result.modalOpened).toBe(true);
    // The modal is titled with the cell's category combination.
    expect(result.modalTitle).toContain('SEX');
    expect(result.iconWhenOff).toBe(false);
   } finally {
    // An icon left parented to a cell survives until the next repaint and pollutes whatever
    // the following sections read out of the grid.
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

  // #### Scrolling
  await softStep('Scrolling', async () => {
   try {
    // The category scroll sliders are in the DOM even when everything fits (refdoc
    // trellis_plot.md), so "a slider exists" grades nothing — the handle's EXTENT does. Packing is
    // off so the overflow is structural, and the rendered counts are asserted alongside the ratio.
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
      await new Promise((r) => setTimeout(r, 2500));
      return {x: [...tp.props.xColumnNames], y: [...tp.props.yColumnNames]};
    }, {x, y});

    // Everything fits: SEX (2) x RACE (4) stays well under the clamp, so both axes put
    // every slot on screen and both handles span (almost) their whole track.
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

    // Horizontal overflow: three X split columns put far more slots on the axis than the clamp
    // renders, so the horizontal handle must contract. No column is reused across the two axes —
    // a column split against itself would make most combinations structurally impossible.
    const overX = await configure(['SEX', 'DIS_POP', 'RACE'], ['SEVERITY']);
    expect(overX.x).toEqual(['SEX', 'DIS_POP', 'RACE']);
    const overRawX = await cardinality(['SEX', 'DIS_POP', 'RACE']);
    const overGeomX = await gridGeometry(page);
    // The window really is smaller than the axis: what is off screen is what the
    // shortened handle offers to scroll to.
    expect(overGeomX.cols).toBe(axisViewportCount(overRawX, false));
    expect(overGeomX.cols).toBeLessThan(overRawX);
    const scrolledX = await scrollSliderExtent(page, 'x');
    expect(scrolledX.present).toBeGreaterThan(0);
    // Cut the NOT-MEASURED sentinel off BEFORE the comparison: a handle that vanished reads
    // -1 and "satisfies" the inequality without shrinking. The fitting ends are already
    // fenced by their 0.5 floor.
    expect(scrolledX.ratio).toBeGreaterThan(0);
    expect(scrolledX.ratio).toBeLessThan(fitsX.ratio);

    // Vertical overflow: the same construction on the Y axis.
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
    // packCategories is not part of restoreCanonical, so it is put back by hand.
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      try {
        tp.props.packCategories = true;
        await new Promise((r) => setTimeout(r, 600));
      } catch (_) { /* best-effort restore */ }
    });
    await restoreCanonical();
   }
  });

  // #### Legend
  await softStep('Legend', async () => {
   try {
    // Legend props (visibility/position) come from LegendMixinLook; assert the legend
    // DOM element presence as the non-prop witness for Always, absence for Never.
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const legendRendered = () => !!tp.root.querySelector('[name="legend"]');
      tp.props.viewerType = 'Scatter plot';
      tp.setOptions({innerViewerLook: {colorColumnName: 'SEX'}});
      await new Promise((r) => setTimeout(r, 1000));

      tp.props.legendVisibility = 'Always';
      await new Promise((r) => setTimeout(r, 600));
      const always = {vis: tp.props.legendVisibility, rendered: legendRendered()};

      // Legend Position is graded on WHERE the legend lands, not on the property echo: a position
      // re-parents the legend root into the matching side panel of the viewer's flex layout
      // (.d4-layout-left / -right / -top / -bottom), a witness with no pixel constants in it.
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
        await new Promise((r) => setTimeout(r, 900));
        positions.push(tp.props.legendPosition);
        slots.push(legendSlot());
      }

      tp.props.legendVisibility = 'Never';
      await new Promise((r) => setTimeout(r, 600));
      const never = {vis: tp.props.legendVisibility, rendered: legendRendered()};
      return {always, positions, slots, never};
    });
    expect(result.always.vis).toBe('Always');
    expect(result.always.rendered).toBe(true);
    expect(result.positions).toEqual(['Left', 'Right', 'Top', 'Bottom']);
    expect(result.slots).toEqual(['left', 'right', 'top', 'bottom']);
    expect(result.never.vis).toBe('Never');
    expect(result.never.rendered).toBe(false);

    // ===== GROK-20432 guard: a legend-driven refresh must not clobber the inner Box plot's look
    // state. The read channel is the Box plot tab of the Context Panel (the trellis exposes no
    // getter for its inner-viewer instance) — a cross-channel product read, not a prop echo. =====

    // Step 5: Box plot inner viewer with its color legend back on screen.
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      tp.props.viewerType = 'Box plot';
      tp.props.legendVisibility = 'Always';
      await new Promise((r) => setTimeout(r, 2500));
    });
    await expect(page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell')).toHaveCount(canonicalCellCount);

    // Step 6: open the Box plot tab of the Context Panel.
    await openBoxPlotTab(page);
    await expect(
      page.locator('[name="viewer-Trellis-plot"] [name="legend"] .d4-legend-item .d4-legend-cross').first(),
    ).toBeAttached({timeout: 15000});

    // Step 7: Show All Categories to its non-default value (default is unchecked).
    const setAll = await setShowAllCategories(page, true);
    expect(setAll).toBe(true);

    // Step 8: switch ONE real legend category off, graded on the observed dimming, never on a
    // cross element being present. The PRE-click witness is opacity alone: .d4-legend-item-current
    // is absent while the legend is pristine [DOM 2026-08-12], so asserting it there could only fail.
    const cats = await legendCategories(page);
    console.log(`[Legend] legend categories = ${JSON.stringify(cats)}`);
    // THREE, not two: step 11 grades the value-label click as an EXCLUSIVE select, which needs a
    // third entry that was shown and is pushed off by a click aimed elsewhere — with two entries
    // that is indistinguishable from "the other one was switched off in step 10".
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
    // The switch-off produced a non-empty subset: the survivors now carry the marker that
    // nobody carried a moment ago.
    const afterFirstSnap = await legendSnapshot(page);
    expect(afterFirstSnap.filter((it) => it.name !== cats[0]).every((it) => it.current && it.opacity >= 1)).toBe(true);

    // Step 9: Show All Categories survived the legend-driven refresh (GROK-20432).
    const afterFirst = await showAllCategoriesState(page);
    console.log(`[Legend] showAllCategories after single uncheck = ${afterFirst}`);
    expect(afterFirst).toBe(true);

    // Step 10: a SECOND uncheck in sequence, with no panel re-open in between, leaves the value
    // alone and raises no new console error. The second category must be active going in and the
    // first still off coming out — otherwise this would be one uncheck plus a no-op.
    const errBeforeSeq = consoleErrors.length;
    const pageErrBeforeSeq = pageErrors.length;
    const secondBefore = await legendItemState(page, cats[1]);
    // Both halves are asserted, so a second click that landed on an already-dimmed entry
    // could not pass as a second uncheck.
    expect(secondBefore.opacity).toBe(1);
    expect(secondBefore.active).toBe(true);
    const secondUnchecked = await uncheckLegendCategory(page, cats[1]);
    expect(secondUnchecked.clicked).toBe(true);
    expect(secondUnchecked.opacity).toBeLessThan(1);
    // Presence first, dimming second: legendItemState() reports opacity -1 for an entry that
    // is not in the legend at all, and -1 satisfies "still dimmed" without the entry ever
    // being there.
    const firstStillOff = await legendItemState(page, cats[0]);
    expect(firstStillOff.present).toBe(true);
    expect(firstStillOff.opacity).toBeLessThan(1);
    const afterSecond = await showAllCategoriesState(page);
    console.log(`[Legend] showAllCategories after two sequential unchecks = ${afterSecond}`);
    expect(afterSecond).toBe(true);
    expect(consoleErrors.length).toBe(errBeforeSeq);
    expect(pageErrors.length).toBe(pageErrBeforeSeq);

    // Step 11: the OTHER legend affordance, and the way back to "everything shown". A value-label
    // click is an exclusive select, not a re-check [DOM 2026-08-12]. 11a: select ONLY cats[0] —
    // the point is cats[2], which WAS shown and is pushed off by a click that never touched it.
    const onlyFirst = await selectOnlyLegendCategory(page, cats[0]);
    expect(onlyFirst.clicked).toBe(true);
    expect(onlyFirst.opacity).toBeGreaterThanOrEqual(1);
    expect(onlyFirst.active).toBe(true);
    const exclusiveSnap = await legendSnapshot(page);
    expect(exclusiveSnap.filter((it) => it.opacity >= 1).map((it) => it.name)).toEqual([cats[0]]);
    expect(exclusiveSnap.find((it) => it.name === cats[2])!.opacity).toBeLessThan(1);
    // 11b: back to the pristine legend through the observed reset — crossing out the last
    // still-shown category restores "everything shown" rather than hiding the final series.
    // This is the only channel that widens the subset back to all.
    expect(await resetLegendViaLastCross(page, cats[0])).toBe(true);
    const resetSnap = await legendSnapshot(page);
    expect(resetSnap.length).toBe(cats.length);
    expect(resetSnap.every((it) => it.opacity >= 1)).toBe(true);
    expect(resetSnap.some((it) => it.current)).toBe(false);
    // Neither the exclusive select nor the reset disturbs the inner Box plot's look state
    // (GROK-20432) or the error floor.
    const afterRecheck = await showAllCategoriesState(page);
    console.log(`[Legend] showAllCategories after exclusive select + reset = ${afterRecheck}`);
    expect(afterRecheck).toBe(true);
    expect(consoleErrors.length).toBe(errBeforeSeq);
    expect(pageErrors.length).toBe(pageErrBeforeSeq);
   } finally {
    // Step 12: the guard left a Box plot inner viewer with a live legend — put the legend
    // back to the section's terminal Never state and restore the canonical trellis.
    await page.evaluate(async () => {
      try {
        const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
        tp.props.legendVisibility = 'Never';
        await new Promise((r) => setTimeout(r, 400));
      } catch (_) { /* best-effort restore */ }
    });
    await restoreCanonical();
   }
  });

  // #### Context menu
  await softStep('Context menu', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.viewerType = 'Scatter plot';
      await new Promise((r) => setTimeout(r, 800));
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const cell = root.querySelectorAll('.d4-trellis-plot-cell')[0];
      const r = cell.getBoundingClientRect();
      cell.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2, clientX: r.left + r.width / 2, clientY: r.top + r.height / 2}));
      await new Promise((res) => setTimeout(res, 900));
      // This step grades menu COMPOSITION — which commands the trellis and its inner viewer put
      // in the menu — and deliberately not reachability, so the read is a textContent sweep.
      // Presence here is NOT permission to actuate: commands go through the click helpers.
      const labels = Array.from(document.querySelectorAll('.d4-menu-item-label')).map((el) => (el as HTMLElement).textContent?.trim());
      document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
      return {
        hasInnerGroup: labels.includes('Scatter plot'),
        hasProperties: labels.includes('Properties...'),
        // Scatter-plot-specific inner items — present only because the inner viewer's
        // own context menu is delegated into this menu.
        hasLasso: labels.includes('Lasso Tool'),
        hasRegression: labels.includes('Show Regression Line'),
      };
    });
    // Step 1: a group named after the inner viewer type plus the standard items.
    expect(result.hasInnerGroup).toBe(true);
    expect(result.hasProperties).toBe(true);
    expect(result.hasLasso).toBe(true);
    expect(result.hasRegression).toBe(true);
  });

  // #### Inner viewer properties
  await softStep('Inner viewer properties', async () => {
   try {
    // Step 1: Scatter inner (with a known X/Y look) so the inner-viewer property tab exists.
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      tp.setOptions({innerViewerLook: {xColumnName: 'WEIGHT', yColumnName: 'HEIGHT'}});
      await new Promise((r) => setTimeout(r, 1200));
    });

    // Step 2: open the properties Context Panel via the title-bar gear.
    await v.openViewerGear(page, 'Trellis plot');
    await page.locator('.property-grid').first().waitFor({timeout: 10000});

    // Step 3: the Context Panel exposes two tabs — 'Trellis' and one named after the
    // inner viewer type. The tab-header name selectors are class-1 (live-recon 2026-08-06).
    const tabs = await page.evaluate(() => ({
      hasTrellisTab: !!document.querySelector('.d4-tab-header[name="Trellis"]'),
      hasInnerTab: !!document.querySelector('.d4-tab-header[name="Scatter plot"]'),
    }));
    expect(tabs.hasTrellisTab).toBe(true);
    expect(tabs.hasInnerTab).toBe(true);

    // Steps 4-5: change an inner viewer property (X/Y column) — the change applies to
    // ALL cells (every sampled populated cell's canvas repaints off its prior frame),
    // not merely "a canvas is present".
    const idxs = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const cells = root.querySelectorAll('.d4-trellis-plot-cell');
      const out: number[] = [];
      for (let i = 0; i < cells.length && out.length < 2; i++) if (cells[i].querySelector('canvas')) out.push(i);
      return out;
    });
    const cellHashes = (idxsArg: number[]) => page.evaluate((idxs) => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const hash = (i: number) => {
        const cv = root.querySelectorAll('.d4-trellis-plot-cell')[i]?.querySelector('canvas') as HTMLCanvasElement | null;
        if (!cv) return null;
        try {
          const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
          let h = 0;
          for (let i = 0; i < img.length; i += 4) h = (h * 31 + ((img[i] << 16) | (img[i + 1] << 8) | img[i + 2])) % 2147483647;
          return h;
        } catch { return null; }
      };
      return idxs.map(hash);
    }, idxsArg);
    const before = await cellHashes(idxs);
    const after = await page.evaluate(async (idxsArg) => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.setOptions({innerViewerLook: {xColumnName: 'AGE', yColumnName: 'WEIGHT'}});
      await new Promise((r) => setTimeout(r, 1500));
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const hash = (i: number) => {
        const cv = root.querySelectorAll('.d4-trellis-plot-cell')[i]?.querySelector('canvas') as HTMLCanvasElement | null;
        if (!cv) return null;
        try {
          const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
          let h = 0;
          for (let i = 0; i < img.length; i += 4) h = (h * 31 + ((img[i] << 16) | (img[i + 1] << 8) | img[i + 2])) % 2147483647;
          return h;
        } catch { return null; }
      };
      return {type: tp.props.viewerType, hashes: idxsArg.map(hash)};
    }, idxs);
    expect(after.type).toBe('Scatter plot');
    expect(idxs.length).toBeGreaterThan(0);
    for (let i = 0; i < idxs.length; i++) {
      expect(before[i]).not.toBeNull();
      expect(after.hashes[i]).not.toBeNull();
      expect(after.hashes[i]).not.toBe(before[i]);
    }
   } finally {
    await restoreCanonical();
   }
  });

  // #### Use in Trellis
  await softStep('Use in Trellis', async () => {
    // 'Use in Trellis' is NOT a flat top-level item: it sits in the source viewer's 'General'
    // group (viewer_context_commands.dart:61-63), and a real click without opening the group is
    // refused [DOM 2026-08-12]. The finally restores a clean trellis for the sections below.
    try {
      // Step 3: right-click a Scatter plot > General | Use in Trellis -> a trellis with
      // the scatter as inner viewer and its X/Y/color settings preserved.
      await page.evaluate(async () => {
        const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
        if (tp) tp.close();
        await new Promise((res) => setTimeout(res, 600));
        const sp = grok.shell.tv.addViewer('Scatter plot') as any;
        sp.setOptions({xColumnName: 'AGE', yColumnName: 'HEIGHT', colorColumnName: 'SEX'});
        await new Promise((res) => setTimeout(res, 1200));
      });
      await page.locator('[name="viewer-Scatter-plot"]').first().click({button: 'right'});
      await clickMenuItemInGroup(page, 'General', 'Use in Trellis');
      const scatterResult = await page.evaluate(async () => {
        let newTp: any = null;
        for (let i = 0; i < 40 && !newTp; i++) {
          newTp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot');
          if (!newTp) await new Promise((res) => setTimeout(res, 250));
        }
        await new Promise((res) => setTimeout(res, 800));
        // props.innerViewerLook is an obfuscated Dart object; the preserved inner settings
        // are readable as plain JSON only via getOptions(true).
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

      // Steps 6-15: the same flow from the Bar chart, Histogram, Line chart and Box plot
      // menus.
      const useInTrellis = async (viewerType: string, look: any) => {
        await page.evaluate(async ({viewerType, look}) => {
          for (const vw of Array.from(grok.shell.tv.viewers) as any[]) if (vw.type === 'Trellis plot') vw.close();
          await new Promise((res) => setTimeout(res, 500));
          const src = grok.shell.tv.addViewer(viewerType) as any;
          src.setOptions(look);
          await new Promise((res) => setTimeout(res, 1200));
        }, {viewerType, look});
        await page.locator('[name="viewer-' + viewerType.replace(/\s+/g, '-') + '"]').first().click({button: 'right'});
        await clickMenuItemInGroup(page, 'General', 'Use in Trellis');
        return page.evaluate(async (viewerType) => {
          let newTp: any = null;
          for (let i = 0; i < 40 && !newTp; i++) {
            newTp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot');
            if (!newTp) await new Promise((res) => setTimeout(res, 250));
          }
          await new Promise((res) => setTimeout(res, 800));
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
      // A leg that throws can leave the popup open; ESC first, so restoreCanonical's viewer
      // closes are not swallowed by a modal menu.
      await page.keyboard.press('Escape').catch(() => {});
      await restoreCanonical();
    }
  });

  // #### Auto layout
  await softStep('Auto layout', async () => {
   try {
    // The resize is a real page.setViewportSize, not a CSS poke: Auto Layout reacts to the
    // actual viewport, so nothing short of a genuine resize exercises it.
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.autoLayout = true;
      tp.props.showXLabels = true;
      tp.props.showYLabels = true;
      await new Promise((r) => setTimeout(r, 600));
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
    // THE LABEL HALF IS GRADED BY NODE COUNT, NEVER BY BOX [DOM 2026-08-12]. Auto Layout hides a
    // strip by not creating it — `_showXLabels` (trellis_plot_core.dart:1019) / `_showYLabels`
    // (:982) gate render calls that run AFTER clear() (:1866, :1872) — and SVG text keeps a box.
    const cards = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {sex: df.col('SEX').categories.length, race: df.col('RACE').categories.length};
    });
    const expectedX = axisViewportCount(cards.sex, false);
    const expectedY = axisViewportCount(cards.race, false);

    // Step 1: baseline at the full window — the panel is on screen and both strips carry one
    // label per rendered category, counted off the live cardinalities through the clamp.
    expect(await vsVisible()).toBe(true);
    const wideLabels = await categoryLabels(page);
    expect(wideLabels.x.length).toBe(expectedX);
    expect(wideLabels.y.length).toBe(expectedY);

    // Step 2: shrink the viewport -> the panel goes and BOTH strips empty.
    await page.setViewportSize({width: 500, height: 400});
    await page.waitForTimeout(1500);
    const smallVisible = await vsVisible();
    const smallLabels = await categoryLabels(page);
    expect(smallVisible).toBe(false);
    expect(smallLabels.x.length).toBe(0);
    expect(smallLabels.y.length).toBe(0);

    // Step 3: restore the large viewport -> panel and both strips come back in full.
    await page.setViewportSize({width: 1920, height: 1080});
    await page.waitForTimeout(1500);
    const largeVisible = await vsVisible();
    const restoredLabels = await categoryLabels(page);
    expect(largeVisible).toBe(true);
    expect(restoredLabels.x.length).toBe(expectedX);
    expect(restoredLabels.y.length).toBe(expectedY);

    // Step 4: the two axes are governed SEPARATELY — X goes at viewer width 200, Y only at ~140 +
    // the Y control's width — so there is a band where the horizontal labels are gone and the
    // vertical ones are still drawn. It is ~25 viewer px wide, hence calibrated, never hardcoded.
    const wideWidth = await viewerWidth();
    await page.setViewportSize({width: 1420, height: 1080});
    await page.waitForTimeout(1200);
    const midWidth = await viewerWidth();
    const slope = Math.max((wideWidth - midWidth) / 500, 0.05);
    // NON-ASSERTIVE DIAGNOSTICS: the per-rung readings are the only thing that tells "the sweep
    // never entered the band" (widen the range) from "this layout has no band at all" (take the
    // rung out). Those are opposite repairs and nothing else in the run distinguishes them.
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
      // Re-read after a further settle: a relayout that lagged the first read must not be
      // able to book a band the settled DOM does not show.
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
    await page.waitForTimeout(1200);

    // Step 5: THE ANTI-VACUITY CONTROL for step 2. With Auto Layout off the very same 500x400
    // window keeps the panel and both strips fully populated — every label survives viewer widths
    // of 113, 53 and even 12px [DOM 2026-08-12] — so the zeroes above are Auto Layout hiding them.
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.autoLayout = false;
      await new Promise((r) => setTimeout(r, 500));
    });
    await page.setViewportSize({width: 500, height: 400});
    await page.waitForTimeout(1500);
    const offSmallVisible = await vsVisible();
    const offSmallLabels = await categoryLabels(page);
    expect(offSmallVisible).toBe(true);
    expect(offSmallLabels.x.length).toBe(expectedX);
    expect(offSmallLabels.y.length).toBe(expectedY);
   } finally {
    // A thrown assertion must not leave every later section running in a 500x400 window.
    await page.setViewportSize({width: 1920, height: 1080});
    await page.waitForTimeout(1000);
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

  // #### Title and description
  await softStep('Title and description', async () => {
   try {
    // Graded on what reaches the screen: the title text in the viewer's panel, and Description
    // Position on the flex side panel .d4-viewer-description lands in. Scoped to the VIEWER ROOT —
    // the dock panel's title bar shows the viewer name whatever Show Title is set to.
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
      await new Promise((r) => setTimeout(r, 1500));
    });
    expect(await titleShown('My Trellis')).toBe(true);

    const slots: (string | null)[] = [];
    for (const pos of ['Bottom', 'Top', 'Left', 'Right']) {
      await page.evaluate(async (p) => {
        const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
        tp.props.descriptionPosition = p;
        await new Promise((r) => setTimeout(r, 1200));
      }, pos);
      slots.push(await descriptionSlot());
    }
    expect(slots).toEqual(['bottom', 'top', 'left', 'right']);

    // Round-trip: switching Show Title off takes the text back off the screen, so the
    // presence reading above is a real state and not a match on something permanent.
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.showTitle = false;
      await new Promise((r) => setTimeout(r, 1200));
    });
    expect(await titleShown('My Trellis')).toBe(false);
   } finally {
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      try {
        tp.props.description = '';
        tp.props.title = '';
        tp.props.showTitle = false;
        await new Promise((r) => setTimeout(r, 600));
      } catch (_) { /* best-effort restore */ }
    });
   }
  });

  // #### Label orientation
  await softStep('Label orientation', async () => {
   try {
    // Graded on the rendered labels, not on the property: the rotate() transform on each
    // label <text> is 0 horizontal / -90 vertical, so Horz and Vert are two DOM states.
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['RACE'];
      tp.props.yColumnNames = ['SEX'];
      tp.props.showXLabels = true;
      tp.props.showYLabels = true;
      await new Promise((r) => setTimeout(r, 1800));
    });
    const setOrientation = async (axis: 'x' | 'y', value: string) => {
      await page.evaluate(async (o) => {
        const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
        if (o.axis === 'x') tp.props.xLabelsOrientation = o.value;
        else tp.props.yLabelsOrientation = o.value;
        await new Promise((r) => setTimeout(r, 1500));
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

    // Auto is an HONEST FLOOR: the chosen angle depends on the label widths against the
    // available cell box, so neither value is a correct expectation — only that labels are
    // still rendered and the angle is one of the two the renderer can emit.
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

  // #### Pick Up / Apply
  await softStep('Pick Up / Apply', async () => {
   try {
    // Step 6: pick up the first trellis's settings, apply to the second -> the second
    // matches (inner type, axes, legend position, title). Steps 7/8: post-Apply the two
    // plots are independent (axis / range-slider changes on one do not affect the other).
    await page.evaluate(async () => {
      const tp1 = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      grok.shell.tv.addViewer('Trellis plot');
      await new Promise((r) => setTimeout(r, 1200));

      tp1.props.yColumnNames = ['RACE'];
      tp1.props.viewerType = 'Bar chart';
      tp1.setOptions({innerViewerLook: {splitColumnName: 'RACE', valueColumnName: 'AGE'}});
      tp1.props.legendVisibility = 'Always';
      tp1.props.legendPosition = 'Top';
      tp1.props.showTitle = true;
      tp1.props.title = 'First Trellis';
      await new Promise((r) => setTimeout(r, 1000));
    });

    // 'Pick Up' and 'Apply' are children of the 'Pick Up / Apply' group, so both legs hover
    // the group header with the real mouse, wait for the child to become visible and click
    // it for real — the same actuation 'Use in Trellis' uses.
    const trellisCell = (n: number) => page.locator('[name="viewer-Trellis-plot"]').nth(n).locator('.d4-trellis-plot-cell').first();
    await trellisCell(0).click({button: 'right', position: {x: 6, y: 6}});
    await clickMenuItemInGroup(page, 'Pick Up / Apply', 'Pick Up');
    await page.waitForTimeout(800);
    await trellisCell(1).click({button: 'right', position: {x: 6, y: 6}});
    await clickMenuItemInGroup(page, 'Pick Up / Apply', 'Apply');
    await page.waitForTimeout(1200);

    const result: any = await page.evaluate(async () => {
      const tps = Array.from(grok.shell.tv.viewers).filter((v: any) => v.type === 'Trellis plot') as any[];
      const applied = {type: tps[1]?.props.viewerType, title: tps[1]?.props.title, legendPos: tps[1]?.props.legendPosition,
        x: [...(tps[1]?.props.xColumnNames ?? [])], y: [...(tps[1]?.props.yColumnNames ?? [])]};
      const source = {x: [...tps[0].props.xColumnNames], y: [...tps[0].props.yColumnNames]};

      // Step 7: change the X axis on the first -> the second is unaffected.
      const tp2TypeBefore = tps[1]?.props.viewerType;
      tps[0].props.xColumnNames = ['CONTROL'];
      await new Promise((r) => setTimeout(r, 1000));
      const tp2XAfterFirstChange = [...tps[1].props.xColumnNames];
      const tp2Independent = tps[1]?.props.viewerType === tp2TypeBefore;

      // Leave BOTH plots open — Step 8 below drives a TRUSTED range-slider drag on
      // the second plot, which needs real page.mouse input outside this evaluate.
      return {applied, source, tp2XAfterFirstChange, tp2Independent};
    });
    expect(result.applied.type).toBe('Bar chart');
    expect(result.applied.title).toBe('First Trellis');
    expect(result.applied.legendPos).toBe('Top');
    expect(result.applied.x).toEqual(result.source.x);
    expect(result.applied.y).toEqual(result.source.y);
    expect(result.tp2Independent).toBe(true);
    expect(result.tp2XAfterFirstChange).not.toContain('CONTROL');

    // Step 8: a range-slider drag on the SECOND plot must not touch the FIRST. Both halves are
    // graded on the cell canvases — the second's own repaint makes the drag attributable. The
    // first plot's look JSON is a companion only: the axis range window is not in the look.
    const step8Setup = await page.evaluate(async () => {
      const tps = Array.from(grok.shell.tv.viewers).filter((v: any) => v.type === 'Trellis plot') as any[];
      const tp1 = tps[0], tp2 = tps[1];
      tp1.props.viewerType = 'Scatter plot';
      tp1.props.xColumnNames = ['SEX'];
      tp1.props.yColumnNames = ['RACE'];
      // Second plot: an inner-axis slider (in a .d4-range-selector wrapper) on ITS root.
      tp2.props.viewerType = 'Scatter plot';
      tp2.props.xColumnNames = ['SEX'];
      tp2.props.yColumnNames = ['RACE'];
      tp2.setOptions({innerViewerLook: {xColumnName: 'WEIGHT', yColumnName: 'HEIGHT'}});
      tp2.props.globalScale = true;
      tp2.props.showRangeSliders = true;
      tp2.props.showXAxes = 'Always';
      await new Promise((r) => setTimeout(r, 2200));
      const root1 = document.querySelectorAll('[name="viewer-Trellis-plot"]')[0] as HTMLElement;
      const root2 = document.querySelectorAll('[name="viewer-Trellis-plot"]')[1] as HTMLElement;
      const hashRoot = (root: HTMLElement) => {
        const hs: (number | null)[] = [];
        for (const cell of Array.from(root.querySelectorAll('.d4-trellis-plot-cell'))) {
          const cv = cell.querySelector('canvas') as HTMLCanvasElement | null;
          if (!cv) { hs.push(null); continue; }
          try {
            const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
            let h = 0;
            for (let i = 0; i < img.length; i += 4) h = (h * 31 + ((img[i] << 16) | (img[i + 1] << 8) | img[i + 2])) % 2147483647;
            hs.push(h);
          } catch { hs.push(null); }
        }
        return hs;
      };
      const innerX = root2.querySelector('.d4-range-selector > svg[type="range-slider"][name="x-slider"]') as SVGElement | null;
      let box: {x: number; y: number; w: number; h: number} | null = null;
      if (innerX) {
        const wrap = innerX.closest('.d4-range-selector') as HTMLElement;
        const wb = wrap.getBoundingClientRect();
        wrap.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: wb.left + wb.width / 2, clientY: wb.top + wb.height / 2}));
        await new Promise((r) => setTimeout(r, 400));
        const b = innerX.getBoundingClientRect();
        box = {x: b.x, y: b.y, w: b.width, h: b.height};
      }
      return {box, hash1: hashRoot(root1), hash2: hashRoot(root2), look1: JSON.stringify(tp1.getOptions(true)?.look ?? null)};
    });
    expect(step8Setup.box).not.toBeNull();
    expect(step8Setup.box!.w).toBeGreaterThan(0);
    await page.mouse.move(step8Setup.box!.x + step8Setup.box!.w - 4, step8Setup.box!.y + step8Setup.box!.h / 2);
    await page.mouse.down();
    await page.mouse.move(step8Setup.box!.x + step8Setup.box!.w * 0.45, step8Setup.box!.y + step8Setup.box!.h / 2, {steps: 12});
    await page.mouse.up();
    await page.waitForTimeout(1500);
    const step8After = await page.evaluate(() => {
      const tps = Array.from(grok.shell.tv.viewers).filter((v: any) => v.type === 'Trellis plot') as any[];
      const roots = document.querySelectorAll('[name="viewer-Trellis-plot"]');
      const hashRoot = (root: HTMLElement) => {
        const hs: (number | null)[] = [];
        for (const cell of Array.from(root.querySelectorAll('.d4-trellis-plot-cell'))) {
          const cv = cell.querySelector('canvas') as HTMLCanvasElement | null;
          if (!cv) { hs.push(null); continue; }
          try {
            const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
            let h = 0;
            for (let i = 0; i < img.length; i += 4) h = (h * 31 + ((img[i] << 16) | (img[i + 1] << 8) | img[i + 2])) % 2147483647;
            hs.push(h);
          } catch { hs.push(null); }
        }
        return hs;
      };
      return {
        hash1: hashRoot(roots[0] as HTMLElement),
        hash2: hashRoot(roots[1] as HTMLElement),
        look1: JSON.stringify(tps[0].getOptions(true)?.look ?? null),
      };
    });
    const secondChanged = step8Setup.hash2.some((h, i) => h !== null && step8After.hash2[i] !== null && h !== step8After.hash2[i]);
    expect(secondChanged).toBe(true);
    // Mirror of `secondChanged` on the first plot: comparable = the cells whose hash came
    // back on BOTH sides, so a vanished canvas can neither book a difference nor stand in
    // for a match. At least one such cell must exist, and none of them may have moved.
    const comparable1 = step8Setup.hash1.filter((h, i) => h !== null && step8After.hash1[i] !== null).length;
    const firstMoved = step8Setup.hash1
      .map((h, i) => (h !== null && step8After.hash1[i] !== null && h !== step8After.hash1[i] ? i : -1))
      .filter((i) => i >= 0);
    expect(comparable1).toBeGreaterThan(0);
    expect(firstMoved).toEqual([]);
    // Companion, not the witness: the first plot's configuration must not move either.
    expect(step8After.look1).toBe(step8Setup.look1);
   } finally {
    // A leftover second plot would break the Layout section below.
    await restoreCanonical();
   }
  });

  // #### Layout and Project save/restore
  await softStep('Layout and Project save/restore', async () => {
   try {
    // Step 3: applying the saved layout must restore EXACTLY the saved viewer set — the
    // extras added in between have to be gone, not merely "fewer than before".
    const result = await page.evaluate(async () => {
      const r: any = {};
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      r.viewersAtSave = Array.from(grok.shell.tv.viewers).map((v: any) => v.type).sort();
      await new Promise((res) => setTimeout(res, 1000));

      grok.shell.tv.addViewer('Histogram');
      grok.shell.tv.addViewer('Bar chart');
      await new Promise((res) => setTimeout(res, 1200));
      r.viewersBefore = Array.from(grok.shell.tv.viewers).map((v: any) => v.type).sort();

      const saved = await grok.dapi.layouts.find(layoutId);
      grok.shell.tv.loadLayout(saved);
      await new Promise((res) => setTimeout(res, 3000));
      r.viewersAfter = Array.from(grok.shell.tv.viewers).map((v: any) => v.type).sort();

      await grok.dapi.layouts.delete(saved);
      return r;
    });
    expect(result.viewersBefore.length).toBeGreaterThan(result.viewersAfter.length);
    expect(result.viewersAfter).toEqual(result.viewersAtSave);
    expect(result.viewersAfter).toContain('Trellis plot');

    // Step 6 (persist entry-point floor): only the entry point lives here — the real ribbon Save
    // button opens the project dialog and raises nothing on the error channels. The
    // reopen-and-read half is closed by the focused trellis-plot-split-and-pick-inner.md.
    const errBeforeSave = consoleErrors.length;
    const pageErrBeforeSave = pageErrors.length;
    const saveBtn = page.locator('[name="button-Save"]').first();
    await saveBtn.waitFor({state: 'visible', timeout: 15000});
    await saveBtn.click();
    await page.waitForTimeout(2000);
    const dialogOpen = await page.evaluate(() => !!document.querySelector('.d4-dialog'));
    const errAfterSave = consoleErrors.length;
    const pageErrAfterSave = pageErrors.length;
    // Dismiss the dialog if one opened (CANCEL before OK -> no project persisted).
    await page.locator('[name="button-CANCEL"]').first().click({timeout: 5000}).catch(() => {});
    await page.waitForTimeout(500);
    expect(dialogOpen).toBe(true);
    expect(errAfterSave).toBe(errBeforeSave);
    expect(pageErrAfterSave).toBe(pageErrBeforeSave);
   } finally {
    // The Histogram and Bar chart added above must not leak into the following sections.
    await restoreCanonical();
   }
  });

  // #### Viewer filter formula
  await softStep('Viewer filter formula', async () => {
    // look.filter is viewer-local — df.filter.trueCount stays invariant when it changes.
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true); df.rows.requestFilter();
      await new Promise((r) => setTimeout(r, 500));
      const dfBefore = df.filter.trueCount;
      const r: any[] = [];
      tp.props.filter = '${AGE} > 40';
      await new Promise((res) => setTimeout(res, 600));
      r.push({filter: tp.props.filter, dfCount: df.filter.trueCount});
      tp.props.filter = '';
      await new Promise((res) => setTimeout(res, 600));
      r.push({filter: tp.props.filter, dfCount: df.filter.trueCount});
      return {r, dfBefore};
    });
    expect(result.r[0].filter).toBe('${AGE} > 40');
    expect(result.r[1].filter).toBe('');
    expect(result.r[0].dfCount).toBe(result.dfBefore);
  });

  // #### Multi Curve inner viewer (and table switching)
  await softStep('Multi Curve inner viewer', async () => {
    const result = await page.evaluate(async (cPath) => {
      // Re-anchor the demog view by CONTENT (the table view whose dataFrame has SEX &
      // RACE): an earlier section can close the original demog view, so a stored view
      // object goes stale — a content lookup does not.
      const hasSexRace = (view: any) => { try { const df = view.dataFrame; return !!(df && df.col('SEX') && df.col('RACE')); } catch { return false; } };
      const hasTrellis = (view: any) => { try { return Array.from(view.viewers ?? []).some((x: any) => x.type === 'Trellis plot'); } catch { return false; } };
      const anchorDemog = (): any => {
        const views = Array.from(grok.shell.views ?? []) as any[];
        const withTrellis = views.find((v) => hasSexRace(v) && hasTrellis(v));
        const target = withTrellis ?? views.find(hasSexRace);
        if (target) { grok.shell.v = target; return target; }
        return grok.shell.tv;
      };
      const demogView = anchorDemog();
      await new Promise((r) => setTimeout(r, 600));
      const demogName = grok.shell.tv.dataFrame.name;
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      if (!tp) return {error: 'Trellis plot not found on demog view'};
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      await new Promise((r) => setTimeout(r, 1200));

      const dfCurves = await grok.dapi.files.readCsv(cPath);
      grok.shell.addTableView(dfCurves);
      await new Promise((r) => setTimeout(r, 1800));
      grok.shell.v = demogView;
      await new Promise((r) => setTimeout(r, 800));
      // TABLE SWITCHING IS HALF THIS SECTION'S TITLE and used to be witnessed by nothing: the
      // assignment sat in an empty `catch` and the curves row count was returned without ever
      // being read. The error is captured instead of dropped, and asserted null below.
      let switchError: string | null = null;
      try { tp.props.table = dfCurves.name; } catch (e) { switchError = String(e); }
      await new Promise((r) => setTimeout(r, 1500));
      // Product-side witness, read off the viewer rather than off the handle the test opened
      // the file with. Name AND row count together: a rebind to a same-named empty frame
      // cannot pass for a real one, and neither can a name echo on a stale binding.
      const boundToCurves = (() => {
        try { const d = tp.dataFrame; return {name: d?.name ?? null, rows: d?.rowCount ?? -1}; }
        catch { return {name: null, rows: -1}; }
      })();

      // Switch to the multi-curve viewer through the control-panel selector: its runtime
      // type is 'MultiCurveViewer' and the props channel with a display-name string does
      // not drive the switch.
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const vs = root.querySelector('[name="viewer selector"]') as HTMLElement;
      vs.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
      await new Promise((r) => setTimeout(r, 600));
      const mc = document.querySelector('.d4-combo-drop-down [name="icon-multicurveviewer"]');
      const mcClicked = !!mc;
      (mc?.closest('.d4-list-item') as HTMLElement | null)?.click();
      await new Promise((r) => setTimeout(r, 1800));
      const vt = tp.props.viewerType;

      // Step 8: the title-bar gear must still open the property grid while the multi-curve
      // inner viewer is active. A synthetic .click() suffices here only because the page
      // runs with body.selenium.
      let gearOpenedPropGrid = false;
      const panel = root.closest('.panel-base') as HTMLElement | null;
      const gear = panel?.querySelector('.panel-titlebar [name="icon-font-icon-settings"]') as HTMLElement | null;
      if (gear) {
        gear.click();
        await new Promise((r) => setTimeout(r, 1200));
        gearOpenedPropGrid = !!document.querySelector('.property-grid');
      }

      // Restore: rebind to the demog table by its real name, back to a Scatter SEX x RACE
      // grid. The way back is graded too — it is the round-trip half of the same switch.
      let restoreError: string | null = null;
      try { tp.props.table = demogName; } catch (e) { restoreError = String(e); }
      await new Promise((r) => setTimeout(r, 1800));
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      await new Promise((r) => setTimeout(r, 1800));
      const boundBack = (() => { try { return tp.dataFrame?.name ?? null; } catch { return null; } })();
      const restoredCells = document.querySelectorAll('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell').length;
      return {viewerType: vt, mcClicked, switchError, boundToCurves, restoreError, boundBack,
        curvesName: dfCurves.name, curvesRows: dfCurves.rowCount, demogName,
        restoredCells, gearOpenedPropGrid};
    }, curvesPath);
    expect(result.mcClicked).toBe(true);
    expect(['MultiCurveViewer', 'Multi curve viewer', 'Curves'].includes(result.viewerType)).toBe(true);
    expect(result.gearOpenedPropGrid).toBe(true);
    // Printed before it is asserted: if the switch is refused on some build, the log says
    // whether it threw, what the viewer ended up bound to and what it should have been —
    // which is what tells a product refusal apart from a stale binding.
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

  // #### To Script
  await softStep('To Script', async () => {
    // 'To JavaScript' is a child of the 'To Script' group — reached by hovering the group
    // header with the real mouse and clicking the child for real, like every other menu
    // actuation in this file.
    await page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell').first().click({button: 'right', position: {x: 6, y: 6}});
    await clickMenuItemInGroup(page, 'To Script', 'To JavaScript');
    const result = await page.evaluate(async () => {
      await new Promise((r) => setTimeout(r, 1500));
      const balloon = document.querySelector('.d4-balloon');
      const generated = !!balloon;
      if (balloon) {
        const close = balloon.querySelector('.close') || balloon.querySelector('[name="icon-times"]');
        if (close) (close as HTMLElement).click();
      }
      return {scriptGenerated: generated};
    });
    expect(result.scriptGenerated).toBe(true);
  });

  // #### Keyboard navigation
  await softStep('Keyboard navigation', async () => {
    const cellLocator = page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell');
    // Re-anchor the demog view by content first — an earlier section can leave the active
    // view or the trellis table elsewhere.
    await page.evaluate(async () => {
      const hasSexRace = (view: any) => { try { const df = view.dataFrame; return !!(df && df.col('SEX') && df.col('RACE')); } catch { return false; } };
      const hasTrellis = (view: any) => { try { return Array.from(view.viewers ?? []).some((x: any) => x.type === 'Trellis plot'); } catch { return false; } };
      const views = Array.from(grok.shell.views ?? []) as any[];
      const target = views.find((v) => hasSexRace(v) && hasTrellis(v)) ?? views.find(hasSexRace);
      if (target) { try { grok.shell.v = target; } catch {} }
      await new Promise((r) => setTimeout(r, 600));
      // Replace any existing trellis with a FRESH one: a viewer whose table was rebound
      // earlier (Multi Curve) keeps a detached keyboard/event wiring, so arrow
      // navigation stops firing current-cell-changed.
      for (const vw of Array.from(grok.shell.tv.viewers) as any[]) if (vw.type === 'Trellis plot') vw.close();
      await new Promise((r) => setTimeout(r, 600));
      const tp = grok.shell.tv.addViewer('Trellis plot') as any;
      tp.props.viewerType = 'Scatter plot';
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      tp.props.onClick = 'None';
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true); df.selection.setAll(false); df.rows.requestFilter();
      await new Promise((r) => setTimeout(r, 1800));
    });
    await expect(cellLocator).toHaveCount(canonicalCellCount);

    // Step 1: a real corner-click sets the current cell; the charts grid is focused
    // explicitly before each arrow. The start cell and its neighbours are computed from the
    // live category order, so each arrow has a NAMED destination to be checked against.
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
    await page.waitForTimeout(700);
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
    // Steps 2-5: Right, Left, Down, Up each move the current cell to a SPECIFIC neighbour.
    // "Some cell carries .d4-trellis-cell-current" is already true after the click above and
    // could not fail whatever the arrows did; the INDEX is what carries the direction.
    const arrow = async (key: string) => {
      await focusGrid();
      await page.keyboard.press(key);
      await page.waitForTimeout(500);
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
    // The event stream has to agree with the DOM: four moves, the two back-moves landing on the
    // same combination. The matchCondition payload crosses the Dart/JS boundary, so it is compared
    // structurally — the exact combination is already pinned by the index assertions above.
    expect(events.length).toBeGreaterThanOrEqual(4);
    expect(JSON.stringify(events[1])).toBe(JSON.stringify(events[3]));
    expect(JSON.stringify(events[0])).not.toBe(JSON.stringify(events[1]));
    expect(JSON.stringify(events[2])).not.toBe(JSON.stringify(events[1]));
    expect(JSON.stringify(events[0])).not.toBe(JSON.stringify(events[2]));

    // Step 6: ESC resets the trellis filter.
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.onClick = 'Filter';
      df.filter.setAll(true); df.rows.requestFilter();
      await new Promise((r) => setTimeout(r, 500));
    });
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(700);
    const filteredBeforeEsc = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filteredBeforeEsc).toBeLessThan(fullRowCount);
    // ESC is handled by the focused charts-grid — focus it before the keypress.
    await focusGrid();
    await page.keyboard.press('Escape');
    await page.waitForTimeout(900);
    const filteredAfterEsc = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filteredAfterEsc).toBe(fullRowCount);
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.onClick = 'None';
      await new Promise((r) => setTimeout(r, 400));
    });
  });

  // #### Undo/redo
  // Stands LAST in the file on purpose: the section closes the trellis and re-opens it
  // through the platform undo stack, the largest cascade risk here.
  await softStep('Undo/redo', async () => {
   try {
    // Steps 1-2: exactly one Trellis plot on the demog view, cells rendered.
    await restoreCanonical();
    await expect(page.locator('[name="viewer-Trellis-plot"]')).toBeVisible();
    expect(await page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell').count()).toBeGreaterThan(0);
    expect(await trellisCount(page)).toBe(1);

    // Step 3: close through the dock title-bar X (the lightest actuation path).
    await closeTrellis(page);
    await page.waitForTimeout(1500);
    expect(await trellisCount(page)).toBe(0);
    await expect(page.locator('[name="viewer-Trellis-plot"]')).toHaveCount(0);

    // Step 4: Ctrl+Z restores the viewer; the reopen itself raises nothing.
    let errBefore = consoleErrors.length;
    let pageErrBefore = pageErrors.length;
    await page.keyboard.press('Control+z');
    await page.waitForTimeout(1800);
    expect(await trellisCount(page)).toBe(1);
    await expect(page.locator('[name="viewer-Trellis-plot"]')).toBeVisible();
    expect(consoleErrors.length).toBe(errBefore);
    expect(pageErrors.length).toBe(pageErrBefore);

    // Step 5: GROK-16560 — the redo after an undo-reopen must not throw. The viewer
    // round-trips back to closed with zero new console/page errors and no error balloon.
    errBefore = consoleErrors.length;
    pageErrBefore = pageErrors.length;
    await page.keyboard.press('Control+Shift+z');
    await page.waitForTimeout(1800);
    expect(await trellisCount(page)).toBe(0);
    expect(consoleErrors.length).toBe(errBefore);
    expect(pageErrors.length).toBe(pageErrBefore);
    expect(await page.locator('.d4-balloon.error').count()).toBe(0);

    // Step 6: a second undo brings the viewer back again.
    await page.keyboard.press('Control+z');
    await page.waitForTimeout(1800);
    expect(await trellisCount(page)).toBe(1);

    // Step 7: the second redo is clean too — the undo manager left no dangling state
    // that makes only the first redo safe.
    errBefore = consoleErrors.length;
    pageErrBefore = pageErrors.length;
    await page.keyboard.press('Control+Shift+z');
    await page.waitForTimeout(1800);
    expect(await trellisCount(page)).toBe(0);
    expect(consoleErrors.length).toBe(errBefore);
    expect(pageErrors.length).toBe(pageErrBefore);
    expect(await page.locator('.d4-balloon.error').count()).toBe(0);
   } finally {
    // Step 8: the section ends with the trellis closed by design — restore the canonical
    // trellis.
    await restoreCanonical();
   }
  });

  v.finishSpec();
});
