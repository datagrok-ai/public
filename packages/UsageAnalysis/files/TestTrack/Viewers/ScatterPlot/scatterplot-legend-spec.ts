/* ---
realizes: [scatterplot.cp.legend, viewers.scatter-plot]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaUI, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';
const RACE_CATEGORIES = ['Asian', 'Black', 'Caucasian', 'Other'];
const SEX_CATEGORIES = ['F', 'M'];
const CONDITIONAL_RANGES = ['634783-634820', '634820-634885'];
const EMPTY_PROBE_COLUMN = 'ZZ_EMPTY_PROBE';
const NO_VALUE_LABEL = '(no value)';

const SELECTION_INK_MAX_FRACTION = 0.5;
const PANEL_RESPECT_MARGIN_PX = 400;
const EMPTY_SATURATED_MAX_FRACTION = 0.1;

const isAmbientError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text) || /Failed to connect to Claude runtime/.test(text) ||
  /powerPreference option is currently ignored/.test(text) ||
  /willReadFrequently/.test(text);

let inProjectSaveWindow = false;
const isBenignError = (text: string) => {
  if (isAmbientError(text)) return true;
  if (inProjectSaveWindow)
    return /Unable to find element in cloned iframe/.test(text) ||
      /Stack trace [A-Za-z]+/.test(text) ||
      /NullError: method not found: '\w+' on null/.test(text);
  return false;
};

const canvasRect = (page: Page) => page.evaluate(() => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'))!;
  const r = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
  return {x: r.x, y: r.y, width: r.width, height: r.height};
});

const waitPlotCanvas = (page: Page) => page.waitForFunction(() => {
  return [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .filter((e) => !e.closest('.d4-dialog'))
    .some((e) => {
      const c = e.querySelector('canvas[name="canvas"]');
      return !!c && c.getBoundingClientRect().width > 0 && c.getBoundingClientRect().height > 0;
    });
}, null, {timeout: 20000});

async function waitBackdrop(page: Page, timeout = 5000): Promise<boolean> {
  return await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
    null, {timeout}).then(() => true).catch(() => false);
}

async function waitFor(page: Page, probe: () => Promise<boolean>, capMs: number): Promise<boolean> {
  return v.pollValue(probe, (ok) => ok, capMs, 100);
}

const menuPopupReady = (page: Page) => page.evaluate(() => {
  const popups = [...document.querySelectorAll('.d4-menu-popup')];
  const last = popups[popups.length - 1];
  if (!last) return false;
  return [...last.querySelectorAll('[name]')].some((e) => {
    const b = e.getBoundingClientRect();
    return b.width > 0 && b.height > 0;
  });
});

const menuItemSized = (page: Page, name: string) => page.evaluate((n: string) => {
  const els = [...document.querySelectorAll(`.d4-menu-popup [name="${n}"]`)];
  const el = els[els.length - 1];
  if (!el) return false;
  const b = el.getBoundingClientRect();
  return b.width > 0 && b.height > 0;
}, name);

async function commitColumn(page: Page, column: string): Promise<void> {
  const text = column.toLowerCase();
  await page.keyboard.press(text[0]);
  await page.waitForTimeout(150); 
  if (text.length > 1) await page.keyboard.type(text.slice(1));
  await page.waitForTimeout(200); 
  await page.keyboard.press('Enter');
  await waitFor(page, async () => (await viewerColumns(page)).markers === column, 4000);
  await settledLegend(page);
}

const pickOnViewer = (page: Page, role: string, column: string) =>
  v.pickColumnViaSelectorTrusted(page, {role, columnName: column});

const MARKER_SELECTOR = '[name="prop-markers"] [name="div-column-combobox-markers"]';

const markerSelector = (page: Page) => page.locator(MARKER_SELECTOR);

async function pickMarkerColumn(page: Page, column: string): Promise<void> {
  await v.openViewerGear(page, 'Scatter plot');
  if (!await revealPropEditor(page, MARKER_SELECTOR, 'marker'))
    throw new Error('the Markers column row never became reachable');
  const sel = markerSelector(page);
  await sel.scrollIntoViewIfNeeded();
  await sel.click();
  if (!await waitBackdrop(page)) throw new Error('marker column popup did not open');
  await commitColumn(page, column);
}

async function clearMarkerColumn(page: Page): Promise<void> {
  await v.openViewerGear(page, 'Scatter plot');
  if (!await revealPropEditor(page, MARKER_SELECTOR, 'marker'))
    throw new Error('the Markers column row never became reachable');
  const sel = markerSelector(page);
  await sel.scrollIntoViewIfNeeded();
  await sel.click();
  if (!await waitBackdrop(page)) throw new Error('marker column popup did not open');
  const pt = await page.evaluate(() => {
    const b = document.querySelector('.d4-column-selector-backdrop')!.getBoundingClientRect();
    return {x: b.x + b.width / 2, y: b.y + 30};
  });
  await page.mouse.click(pt.x, pt.y);
  await waitFor(page, async () => (await viewerColumns(page)).markers === '', 4000);
  await settledLegend(page);
}

interface Legend {
  present: boolean;
  all: number;
  coloring: number;
  extra: number;
  labels: string[];
  colorLabels: string[];

  glyphLabels: string[];
  current: string[];
  dimmed: string[];
}

const readLegend = (page: Page): Promise<Legend> => page.evaluate(() => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'));
  const empty = {present: false, all: 0, coloring: 0, extra: 0, labels: [] as string[],
    colorLabels: [] as string[], glyphLabels: [] as string[], current: [] as string[], dimmed: [] as string[]};
  if (!root) return empty;
  const items = [...root.querySelectorAll('[name="legend"] .d4-legend-item')] as HTMLElement[];
  const txt = (i: Element) => (i.querySelector('.d4-legend-value')?.textContent ?? '').trim();
  const colorItems = items.filter((i) => i.classList.contains('d4-legend-item-coloring'));
  return {
    present: !!root.querySelector('[name="legend"]'),
    all: items.length,
    coloring: colorItems.length,
    extra: items.filter((i) => i.classList.contains('d4-legend-item-extra')).length,
    labels: items.map(txt),
    colorLabels: colorItems.map(txt),
    glyphLabels: items.filter((i) => i.querySelector('i[name="legend-icon-color-picker"]')).map(txt),
    current: items.filter((i) => i.classList.contains('d4-legend-item-current')).map(txt),

    dimmed: items.filter((i) => parseFloat(getComputedStyle(i).opacity) < 0.9).map(txt),
  };
});

async function settledLegend(page: Page): Promise<Legend> {
  let prev = await readLegend(page);
  for (let i = 0; i < 30; i++) {
    await page.waitForTimeout(100);
    const cur = await readLegend(page);
    if (JSON.stringify(cur) === JSON.stringify(prev)) return cur;
    prev = cur;
  }
  return prev;
}

interface Ink {
  nonWhite: number;
  saturated: number;

  pale: number;
}

const readInk = (page: Page): Promise<Ink> => page.evaluate(() => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'))!;
  const canvas = root.querySelector('canvas[name="canvas"]') as HTMLCanvasElement;
  const ctx = canvas.getContext('2d', {willReadFrequently: true})!;
  const d = ctx.getImageData(0, 0, canvas.width, canvas.height).data;
  let nonWhite = 0; let saturated = 0; let pale = 0;
  for (let i = 0; i < d.length; i += 4) {
    if (d[i + 3] === 0) continue;
    const r = d[i]; const g = d[i + 1]; const b = d[i + 2];
    if (r >= 250 && g >= 250 && b >= 250) continue;
    nonWhite++;
    const spread = Math.max(r, g, b) - Math.min(r, g, b);
    if (spread >= 12) saturated++;
    else if (Math.max(r, g, b) >= 200) pale++;
  }
  return {nonWhite, saturated, pale};
});

async function parkPointer(page: Page): Promise<void> {
  await page.mouse.move(4, 4);
  await v.waitForViewerRendered(page, 'Scatter plot', 250);
}

const sameInk = (a: Ink, b: Ink) =>
  a.nonWhite === b.nonWhite && a.saturated === b.saturated && a.pale === b.pale;

async function settledInk(page: Page): Promise<Ink> {
  await parkPointer(page);
  let prev = await readInk(page);
  for (let i = 0; i < 40; i++) {
    await page.waitForTimeout(250);
    const cur = await readInk(page);
    if (sameInk(cur, prev)) return cur;
    prev = cur;
  }
  return prev;
}

async function settledInkAfterChange(page: Page, from: Ink): Promise<Ink> {
  await parkPointer(page);
  let prev = await readInk(page);
  for (let i = 0; i < 60; i++) {
    await page.waitForTimeout(250);
    const cur = await readInk(page);
    if (!sameInk(cur, from) && sameInk(cur, prev)) return cur;
    prev = cur;
  }
  return prev;
}

const legendSelection = async (page: Page) => (await readLegend(page)).current.join('|');

async function clickLegendEntry(page: Page, label: string): Promise<void> {
  const before = await legendSelection(page);
  const pt = await legendPoint(page, label);
  expect(pt, `legend entry ${label}`).not.toBeNull();
  await page.mouse.click(pt!.x, pt!.y);
  await waitFor(page, async () => await legendSelection(page) !== before, 4000);
}

async function clearLegendSelection(page: Page): Promise<void> {
  if (!await legendFiltering(page)) return;
  const current = (await readLegend(page)).current;
  if (!current.length) return;
  const pt = await legendPoint(page, current[0]);
  if (!pt) return;
  await page.mouse.click(pt.x, pt.y);
  await waitFor(page, async () => !await legendFiltering(page), 4000);
}

const legendFiltering = (page: Page) => page.evaluate(() => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'));
  const host = root?.querySelector('[name="legend"]');
  return !!host?.classList.contains('d4-legend-filtering');
});

const legendPoint = (page: Page, label: string) => page.evaluate((l: string) => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'));
  const item = [...(root?.querySelectorAll('[name="legend"] .d4-legend-item') ?? [])]
    .find((i) => (i.querySelector('.d4-legend-value')?.textContent ?? '').trim() === l);
  if (!item) return null;
  const el = item.querySelector('.d4-legend-value') ?? item;
  const b = el.getBoundingClientRect();
  return {x: b.x + b.width / 2, y: b.y + b.height / 2};
}, label);

const filterCount = (page: Page) =>
  page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount as number);

async function settledFilterCountAfterChange(page: Page, from: number): Promise<number> {
  let prev = -1;
  let stable = 0;
  for (let i = 0; i < 70; i++) {
    const c = await filterCount(page);
    if (c === prev) {
      stable++;
      if (c !== from && stable >= 4) return c;
    } else {
      stable = 0;
      prev = c;
    }
    await page.waitForTimeout(250);
  }
  return prev;
}

async function settledFilterCountUnchanged(page: Page, floorMs = 2000): Promise<number> {

  await page.waitForTimeout(floorMs);
  let last = -1;
  let stable = 0;
  for (let i = 0; i < 25; i++) {
    const c = await filterCount(page);
    if (c === last) {
      stable++;
      if (stable >= 3) return c;
    } else {
      stable = 0;
      last = c;
    }
    await page.waitForTimeout(400);
  }
  return last;
}

const viewerColumns = (page: Page) => page.evaluate(() => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return {color: sp.props.colorColumnName, markers: sp.props.markersColumnName};
});

const legendHostCount = (page: Page) => page.evaluate(() => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'));
  return root ? root.querySelectorAll('[name="legend"]').length : 0;
});

const legendVisibility = (page: Page) => page.evaluate(() => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return sp.props.legendVisibility as string;
});

const VISIBILITY_MENU = ['div-Properties...', 'div-Properties...---Legend',
  'div-Properties...---Legend---Legend-Visibility'];

async function settingsPanelReady(page: Page): Promise<boolean> {
  const panelUp = () => page.evaluate(() => !!document.querySelector('[name="prop-category-legend"]'));
  for (let i = 0; i < 4; i++) {
    if (await panelUp()) return true;
    await v.openViewerGear(page, 'Scatter plot');
    await waitFor(page, panelUp, 1000);
  }
  return false;
}

async function revealPropEditor(page: Page, editorSelector: string, category: string): Promise<boolean> {
  const sized = () => page.evaluate((sel: string) => {
    const el = document.querySelector(sel) as HTMLElement | null;
    if (!el || !el.offsetParent) return false;
    const b = el.getBoundingClientRect();
    return b.width > 0 && b.height > 0;
  }, editorSelector);
  for (let i = 0; i < 8; i++) {
    if (await sized()) return true;
    const header = page.locator(`[name="prop-category-${category}"]`);
    if (await header.count() > 0 && await header.isVisible()) await header.click();
    await waitFor(page, sized, 800);
  }
  return false;
}

async function openPlotContextMenu(page: Page): Promise<void> {
  const r = await canvasRect(page);
  await page.mouse.click(r.x + r.width / 2, r.y + r.height / 2, {button: 'right'});
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 8000});
  await waitFor(page, () => menuPopupReady(page), 500);
}

async function setLegendVisibility(page: Page, value: string): Promise<void> {
  let driven = false;
  const editor = '[name="prop-view-legend-visibility"]';
  if (await settingsPanelReady(page) && await revealPropEditor(page, editor, 'legend')) {
    const cell = page.locator(editor);
    await cell.scrollIntoViewIfNeeded();
    await cell.click();
    const select = page.locator('[name="prop-legend-visibility"] select.property-grid-item-editor-spinner');
    await select.waitFor({timeout: 3000}).catch(() => {});
    if (await select.count() > 0) {
      const options = (await select.locator('option').allTextContents()).map((o) => o.trim());
      expect(options).toContain(value);
      await select.selectOption(value);
      await waitFor(page, async () => await legendVisibility(page) === value, 4000);
      driven = true;
    }
  }
  if (!driven) {
    await openPlotContextMenu(page);
    const leafName = `${VISIBILITY_MENU[2]}---${value}`;
    for (let i = 0; i < VISIBILITY_MENU.length; i++) {
      await page.locator(`.d4-menu-popup [name="${VISIBILITY_MENU[i]}"]`).last().hover();

      const next = i + 1 < VISIBILITY_MENU.length ? VISIBILITY_MENU[i + 1] : leafName;
      await waitFor(page, () => menuItemSized(page, next), 700);
    }
    const leaf = page.locator(`.d4-menu-popup [name="${leafName}"]`);
    expect(await leaf.count()).toBeGreaterThan(0);
    await leaf.last().click();
    await waitFor(page, async () => await legendVisibility(page) === value, 4000);
    await page.keyboard.press('Escape');
    await waitFor(page, async () => !await menuPopupReady(page), 500);
  }
  expect(await legendVisibility(page)).toBe(value);
}

async function openColumnHeaderMenu(page: Page, column: string): Promise<void> {
  const pt = await page.evaluate((name: string) => {
    const grid = grok.shell.tv.grid;
    const r = grid.root.getBoundingClientRect();
    const gc = grid.columns.byName(name);
    return {x: r.x + gc.left + gc.width / 2, y: r.y + (grid.colHeaderH || 20) / 2};
  }, column);
  await page.mouse.click(pt.x, pt.y, {button: 'right'});
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 8000});
  await waitFor(page, () => menuPopupReady(page), 400);
}

async function clickHeaderMenuLeaf(page: Page, group: string, leaf: string): Promise<void> {
  await page.locator(`.d4-menu-popup [name="${group}"]`).last().hover();
  await waitFor(page, () => menuItemSized(page, leaf), 600);
  await page.locator(`.d4-menu-popup [name="${leaf}"]`).last().click();

  await waitFor(page, async () => !await menuPopupReady(page), 1200);
}

const conditionalRanges = (page: Page, column: string) => page.evaluate((c: string) => {
  const col = grok.shell.tv.dataFrame.col(c);
  const raw = col.tags['.color-coding-conditional'];
  return {type: col.tags['.color-coding-type'] ?? null, ranges: raw ? Object.keys(JSON.parse(raw)) : []};
}, column);

test('Scatter Plot — Legend Lifecycle, Filter Interplay, Persistence', async ({page}: {page: Page}) => {
  test.setTimeout(1_500_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => { if (!isBenignError(String(e))) pageErrors.push(String(e)); });
  const consoleErrors: string[] = [];
  page.on('console', (m) => {
    if (m.type() === 'error' && !isBenignError(m.text())) consoleErrors.push(m.text());
  });
  const errCount = () => pageErrors.length + consoleErrors.length;

  await loginToDatagrok(page);
  await v.openTable(page, {path: demogPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'scatter-plot', 'Scatter-plot');
  await waitPlotCanvas(page);
  await pickOnViewer(page, 'x', 'WEIGHT');
  await pickOnViewer(page, 'y', 'HEIGHT');
  const fullRowCount = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount as number);
  expect(fullRowCount).toBeGreaterThan(0);

  await softStep('Color legend and the joint Color plus Marker legend', async () => {
    const errBefore = errCount();

    await pickOnViewer(page, 'color', 'RACE');
    const colorOnly = await readLegend(page);
    expect(colorOnly.present).toBe(true);
    expect(colorOnly.colorLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
    expect(colorOnly.all).toBe(RACE_CATEGORIES.length);

    await pickMarkerColumn(page, 'SEX');
    const joint = await readLegend(page);
    expect(joint.coloring).toBe(RACE_CATEGORIES.length);
    expect(joint.extra).toBe(SEX_CATEGORIES.length);
    expect(joint.glyphLabels.sort()).toEqual([...SEX_CATEGORIES].sort());
    expect(joint.all).toBe(RACE_CATEGORIES.length + SEX_CATEGORIES.length);

    await pickOnViewer(page, 'color', 'AGE');
    const numericColor = await readLegend(page);
    expect(numericColor.glyphLabels.sort()).toEqual([...SEX_CATEGORIES].sort());
    expect(numericColor.extra).toBe(SEX_CATEGORIES.length);

    await pickOnViewer(page, 'color', 'RACE');
    const back = await readLegend(page);
    expect(back.coloring).toBe(RACE_CATEGORIES.length);
    expect(back.extra).toBe(SEX_CATEGORIES.length);
    expect(back.all).toBe(RACE_CATEGORIES.length + SEX_CATEGORIES.length);
    expect(back.colorLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
    expect(new Set(back.labels).size).toBe(back.labels.length);

    await pickMarkerColumn(page, 'RACE');
    const same = await readLegend(page);
    expect(same.colorLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
    expect(same.glyphLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
    expect(same.all).toBe(RACE_CATEGORIES.length);

    await pickMarkerColumn(page, 'SEX');
    expect((await viewerColumns(page)).markers).toBe('SEX');
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Clearing the Marker column removes its glyph entries', async () => {
    const errBefore = errCount();
    await pickOnViewer(page, 'color', 'SEX');
    await pickMarkerColumn(page, 'SEX');
    const before = await readLegend(page);
    expect(before.colorLabels.sort()).toEqual([...SEX_CATEGORIES].sort());
    expect(before.glyphLabels.sort()).toEqual([...SEX_CATEGORIES].sort());

    await clearMarkerColumn(page);
    const after = await readLegend(page);

    expect((await viewerColumns(page)).markers).toBe('');
    expect(after.glyphLabels).toEqual([]);
    expect(after.colorLabels.sort()).toEqual([...SEX_CATEGORIES].sort());
    expect(after.present).toBe(true);

    await pickOnViewer(page, 'color', 'RACE');
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Filtered-out categories are absent from the marker legend', async () => {
    const errBefore = errCount();
    await v.openFilterPanel(page);

    await v.applyCategoricalFilter(page, 'RACE', ['Asian', 'Caucasian'], 600);
    const narrowed = await settledFilterCountAfterChange(page, fullRowCount);
    expect(narrowed).toBeLessThan(fullRowCount);
    await pickMarkerColumn(page, 'RACE');
    const afterFilter = await readLegend(page);
    expect(afterFilter.glyphLabels.sort()).toEqual(['Asian', 'Caucasian']);

    await v.resetFilters(page);
    await settledFilterCountAfterChange(page, narrowed);
    await clearMarkerColumn(page);
    await pickMarkerColumn(page, 'RACE');
    const unfiltered = await readLegend(page);
    expect(unfiltered.glyphLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
    await v.applyCategoricalFilter(page, 'RACE', ['Asian', 'Caucasian'], 600);
    await settledFilterCountAfterChange(page, fullRowCount);
    await settledLegend(page);
    const beforeFilter = await readLegend(page);
    expect(beforeFilter.glyphLabels.sort()).toEqual(['Asian', 'Caucasian']);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Filtered-out categories are absent from the marker legend — filter cleared', async () => {
    const before = await filterCount(page);
    await v.resetFilters(page);
    const restored = await settledFilterCountAfterChange(page, before);
    expect(restored).toBe(fullRowCount);
    const legend = await settledLegend(page);
    expect(legend.glyphLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
  });

  await softStep('Clicking a legend category hides the other categories on the plot', async () => {
    const errBefore = errCount();
    await pickOnViewer(page, 'color', 'RACE');
    await pickMarkerColumn(page, 'SEX');
    try {
      const baseline = await settledInk(page);
      expect(baseline.nonWhite).toBeGreaterThan(0);
      const baseCount = await settledFilterCountUnchanged(page, 500);
      expect(baseCount).toBe(fullRowCount);

      await clickLegendEntry(page, 'Asian');
      const legendOnly = await settledInkAfterChange(page, baseline);
      expect(legendOnly.nonWhite)
        .toBeLessThan(Math.round(baseline.nonWhite * SELECTION_INK_MAX_FRACTION));
      expect(await legendFiltering(page)).toBe(true);
      const marked = await readLegend(page);
      expect(marked.current).toContain('Asian');
      expect(marked.dimmed.length).toBeGreaterThan(0);
      expect(await settledFilterCountUnchanged(page)).toBe(baseCount);

      await v.applyCategoricalFilter(page, 'SEX', ['F'], 600);
      const narrowed = await settledFilterCountAfterChange(page, fullRowCount);
      expect(narrowed).toBeLessThan(fullRowCount);
      const legendThenPanel = await settledInkAfterChange(page, legendOnly);
      expect(legendThenPanel.nonWhite)
        .toBeLessThan(legendOnly.nonWhite - PANEL_RESPECT_MARGIN_PX);
      expect((await readLegend(page)).labels).not.toContain('M');

      await clearLegendSelection(page);
      await v.resetFilters(page);
      await settledFilterCountAfterChange(page, narrowed);
      await v.applyCategoricalFilter(page, 'SEX', ['F'], 600);
      await settledFilterCountAfterChange(page, fullRowCount);
      const panelOnly = await settledInk(page);
      await clickLegendEntry(page, 'Asian');
      const panelThenLegend = await settledInkAfterChange(page, panelOnly);
      expect(panelThenLegend).toEqual(legendThenPanel);
      expect(errCount()).toBe(errBefore);
    } finally {
      await clearLegendSelection(page);
      await v.resetFilters(page);
      await settledFilterCountUnchanged(page, 500);
    }
  });

  await softStep('Clicking a legend category hides the other categories on the plot — second click restores the plot exactly',
    async () => {
      const errBefore = errCount();
      const baseline = await settledInk(page);
      try {
        await clickLegendEntry(page, 'Asian');
        expect(await legendFiltering(page)).toBe(true);
        const selected = await settledInkAfterChange(page, baseline);
        await clickLegendEntry(page, 'Asian');

        expect(await settledInkAfterChange(page, selected)).toEqual(baseline);
        const legend = await readLegend(page);
        expect(legend.current).toEqual([]);
        expect(legend.dimmed).toEqual([]);
        expect(await legendFiltering(page)).toBe(false);
        expect(errCount()).toBe(errBefore);
      } finally {
        await clearLegendSelection(page);
      }
    });

  await softStep('Clicking a legend category hides the other categories on the plot — the entry for empty values',
    async () => {
      const errBefore = errCount();
      const baseCount = await settledFilterCountUnchanged(page, 500);
      try {

        await page.evaluate((name: string) => {
          const df = grok.shell.tv.dataFrame;
          const col = df.columns.addNewString(name);
          for (let i = 0; i < df.rowCount; i++)
            col.set(i, i % 3 === 0 ? null : (i % 3 === 1 ? 'alpha' : 'beta'));
        }, EMPTY_PROBE_COLUMN);
        await v.waitForViewerRendered(page, 'Scatter plot', 1200);
        await pickOnViewer(page, 'color', EMPTY_PROBE_COLUMN);
        expect((await settledLegend(page)).colorLabels).toContain(NO_VALUE_LABEL);

        const baseline = await settledInk(page);
        expect(baseline.saturated).toBeGreaterThan(0);

        await clickLegendEntry(page, NO_VALUE_LABEL);
        const selected = await settledInkAfterChange(page, baseline);
        expect(selected.saturated)
          .toBeLessThan(Math.round(baseline.saturated * EMPTY_SATURATED_MAX_FRACTION));
        expect(selected.pale).toBeGreaterThan(baseline.pale);
        expect(await legendFiltering(page)).toBe(true);
        expect((await readLegend(page)).current).toContain(NO_VALUE_LABEL);
        expect(await settledFilterCountUnchanged(page)).toBe(baseCount);

        await v.applyCategoricalFilter(page, 'SEX', ['F'], 600);
        expect(await settledFilterCountAfterChange(page, baseCount)).toBeLessThan(baseCount);
        const withPanel = await settledInkAfterChange(page, selected);
        expect(withPanel.nonWhite).toBeLessThan(selected.nonWhite - PANEL_RESPECT_MARGIN_PX);

        await clearLegendSelection(page);
        const filtered = await filterCount(page);
        await v.resetFilters(page);
        await settledFilterCountAfterChange(page, filtered);
        expect(await settledInkAfterChange(page, withPanel)).toEqual(baseline);
        expect(errCount()).toBe(errBefore);
      } finally {
        await clearLegendSelection(page);
        await v.resetFilters(page);
        await page.evaluate((name: string) => {
          const df = grok.shell.tv.dataFrame;
          if (df.col(name)) df.columns.remove(name);
        }, EMPTY_PROBE_COLUMN);
        await pickOnViewer(page, 'color', 'RACE');
        await pickMarkerColumn(page, 'SEX');
      }
      const peak = await viewerColumns(page);
      expect(peak.color).toBe('RACE');
      expect(peak.markers).toBe('SEX');
    });

  await softStep('A column with conditional color coding still produces a legend', async () => {
    const errBefore = errCount();
    try {
      await page.evaluate(async (path: string) => {
        const df = await grok.dapi.files.readCsv(path);
        grok.shell.addTableView(df);
        await new Promise((resolve) => {
          const s = df.onSemanticTypeDetected.subscribe(() => { s.unsubscribe(); resolve(null); });
          setTimeout(resolve, 5000);
        });
      }, spgiPath);
      await page.waitForFunction(() =>
        !!grok.shell.tv?.dataFrame?.col('CAST Idea ID'), null, {timeout: 30000});

      await openColumnHeaderMenu(page, 'CAST Idea ID');
      await clickHeaderMenuLeaf(page, 'div-Color-Coding', 'div-Color-Coding---Conditional');
      await openColumnHeaderMenu(page, 'CAST Idea ID');
      await clickHeaderMenuLeaf(page, 'div-Color-Coding', 'div-Color-Coding---Edit...');

      const dialog = page.locator('[name^="dialog-Color-coding"]');
      await dialog.waitFor({timeout: 10000});
      const rangeInputs = dialog.locator('input.ui-input-editor');
      for (let guard = 0; guard < 10 && await rangeInputs.count() > CONDITIONAL_RANGES.length; guard++) {
        const before = await rangeInputs.count();
        await dialog.locator('[name="button-Remove-row"]').last().click();
        await waitFor(page, async () => await rangeInputs.count() < before, 600);
      }
      expect(await rangeInputs.count()).toBe(CONDITIONAL_RANGES.length);
      for (let i = 0; i < CONDITIONAL_RANGES.length; i++) {
        await rangeInputs.nth(i).click();
        await page.keyboard.press('Control+A');
        await page.keyboard.type(CONDITIONAL_RANGES[i]);
        await page.keyboard.press('Enter');
        await waitFor(page, async () => await rangeInputs.nth(i).inputValue() === CONDITIONAL_RANGES[i], 800);
      }
      await dialog.locator('[name="button-CLOSE"]').click();
      await waitFor(page, async () => await dialog.count() === 0, 1200);

      const coding = await conditionalRanges(page, 'CAST Idea ID');
      expect(coding.type).toBe('Conditional');
      expect(coding.ranges).toEqual(CONDITIONAL_RANGES);

      await v.addViewerByIcon(page, 'scatter-plot', 'Scatter-plot');
      await waitPlotCanvas(page);
      const axes = await page.evaluate(() => {
        const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
        return {x: sp.props.xColumnName, y: sp.props.yColumnName};
      });
      expect(axes.x).toBe('CAST Idea ID');
      expect(axes.y).toBe('Idea ID');

      await pickOnViewer(page, 'color', 'CAST Idea ID');
      const legend = await readLegend(page);
      expect(legend.present).toBe(true);
      expect(legend.colorLabels).toEqual(CONDITIONAL_RANGES);
      expect(errCount()).toBe(errBefore);
    } finally {

      await page.evaluate(() => grok.shell.v.close());
      await waitFor(page, () => page.evaluate(() =>
        [...grok.shell.tableViews].some((x: any) => x.dataFrame.columns.contains('RACE'))), 1500);
      await page.evaluate(() => {
        for (const view of grok.shell.tableViews) {
          if (view.dataFrame.columns.contains('RACE')) {
            grok.shell.v = view;
            break;
          }
        }
      });

      await waitFor(page, () => page.evaluate(() => {
        const tv = grok.shell.tv;
        return !!tv?.dataFrame?.columns.contains('RACE') &&
          tv.viewers.some((x: any) => x.type === 'Scatter plot');
      }), 1500);
    }
    const peak = await viewerColumns(page);
    expect(peak.color).toBe('RACE');
    expect(peak.markers).toBe('SEX');
  });

  await softStep('Legend visibility hides and restores the legend', async () => {
    const errBefore = errCount();
    const initial = await legendVisibility(page);
    const before = await readLegend(page);
    expect(before.present).toBe(true);
    expect(await legendHostCount(page)).toBe(1);

    try {
      await setLegendVisibility(page, 'Never');
      expect(await legendHostCount(page)).toBe(0);
      expect((await readLegend(page)).present).toBe(false);

      await setLegendVisibility(page, 'Always');
      expect(await legendHostCount(page)).toBe(1);
      const shown = await readLegend(page);
      expect(shown.present).toBe(true);
      expect(shown.all).toBe(before.all);
      expect(shown.colorLabels.sort()).toEqual([...before.colorLabels].sort());
      expect(shown.glyphLabels.sort()).toEqual([...before.glyphLabels].sort());

      await setLegendVisibility(page, initial);
      const restored = await readLegend(page);
      expect(restored.present).toBe(true);
      expect(restored.all).toBe(before.all);
      expect(restored.colorLabels.sort()).toEqual([...before.colorLabels].sort());
      expect(errCount()).toBe(errBefore);
    } finally {

      await page.evaluate((val: string) => {
        const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
        if (sp && sp.props.legendVisibility !== val) sp.props.legendVisibility = val;
      }, initial);
      await v.waitForViewerRendered(page, 'Scatter plot', 1200);
    }
  });

  await softStep('Layout and project persistence at peak configuration — layout round-trip', async () => {
    const peak = await viewerColumns(page);
    expect(peak.color).toBe('RACE');
    expect(peak.markers).toBe('SEX');
    const peakLegend = await readLegend(page);
    expect(peakLegend.coloring).toBe(RACE_CATEGORIES.length);
    expect(peakLegend.extra).toBe(SEX_CATEGORIES.length);

    const layoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      return layout.id as string;
    });

    await waitFor(page, () => page.evaluate(async (id: string) =>
      !!await grok.dapi.layouts.find(id), layoutId), 1500);
    try {
      await page.evaluate(() => grok.shell.tv.addViewer('Histogram'));
      await waitFor(page, () => page.evaluate(() =>
        grok.shell.tv.viewers.some((x: any) => x.type === 'Histogram')), 1200);
      const [clearedColor] = await v.setViewerProps(page, 'Scatter plot',
        [{set: {colorColumnName: ''}, wait: 1000, read: 'colorColumnName'}]);
      await page.evaluate(async (id: string) => {
        grok.shell.tv.loadLayout(await grok.dapi.layouts.find(id));
      }, layoutId);

      await waitFor(page, () => page.evaluate(() => {
        const tv = grok.shell.tv;
        if (!tv) return false;
        const sp = tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
        return !!sp && sp.props.colorColumnName === 'RACE' &&
          !tv.viewers.some((x: any) => x.type === 'Histogram');
      }), 4000);
      const result = await page.evaluate((cleared: string) => {
        const tv = grok.shell.tv;
        const restored = tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
        return {
          clearedColor: cleared,
          hasScatter: tv.viewers.some((x: any) => x.type === 'Scatter plot'),
          hasHistogram: tv.viewers.some((x: any) => x.type === 'Histogram'),
          color: restored?.props.colorColumnName, markers: restored?.props.markersColumnName,
        };
      }, clearedColor);

      expect(result.clearedColor).toBe('');

      expect(result.hasScatter).toBe(true);
      expect(result.hasHistogram).toBe(false);
      expect(result.color).toBe('RACE');
      expect(result.markers).toBe('SEX');

      const legend = await settledLegend(page);
      expect(legend.present).toBe(true);
      expect(legend.colorLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
      expect(legend.glyphLabels.sort()).toEqual([...SEX_CATEGORIES].sort());
    } finally {
      await page.evaluate(async (id: string) => {
        try {
          const saved = await grok.dapi.layouts.find(id);
          if (saved) await grok.dapi.layouts.delete(saved);
        } catch (_) {  }
      }, layoutId);
    }
  });

  await softStep('Layout and project persistence at peak configuration — project save / Close All / reopen',
    async () => {
      const projName = 'zz-scatterplot-legend-probe-' + Date.now();
      let projectId: string | null = null;
      inProjectSaveWindow = true;
      try {
        const saved = await saveProjectViaUI(page, projName);
        projectId = saved.projectId;
        expect(projectId).toBeTruthy();

        await v.closeAllAndWait(page);
        const result = await page.evaluate(async (id: string) => {
          const full = await grok.dapi.projects.find(id);
          await full.open();
          let sp: any = null;
          for (let t = 0; t < 20 && !sp; t++) {
            await new Promise((r) => setTimeout(r, 1000));
            for (const view of grok.shell.tableViews)
              for (const vw of view.viewers)
                if (vw.type === 'Scatter plot') sp = vw;
          }
          const items = sp ? [...sp.root.querySelectorAll('[name="legend"] .d4-legend-item')] : [];
          const txt = (i: Element) => (i.querySelector('.d4-legend-value')?.textContent ?? '').trim();
          return {
            found: !!sp,
            color: sp?.props.colorColumnName, markers: sp?.props.markersColumnName,
            legendPresent: !!sp?.root.querySelector('[name="legend"]'),
            colorLabels: items.filter((i) => i.classList.contains('d4-legend-item-coloring')).map(txt),
            glyphLabels: items.filter((i) => i.querySelector('i[name="legend-icon-color-picker"]')).map(txt),
          };
        }, projectId);

        expect(result.found).toBe(true);
        expect(result.color).toBe('RACE');
        expect(result.markers).toBe('SEX');
        expect(result.legendPresent).toBe(true);
        expect(result.colorLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
        expect(result.glyphLabels.sort()).toEqual([...SEX_CATEGORIES].sort());
      } finally {
        inProjectSaveWindow = false;
        if (projectId) await deleteProjectWithCleanup(page, {projectId});
      }
    });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});
