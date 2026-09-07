/* ---
realizes: [scatterplot.cp.legend, viewers.scatter-plot]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import * as sp from './scatterplot-shared';

declare const grok: any;
declare const DG: any;

// The legend ladder on the local lane; the layout and project round-trips of the same
// scenario live in scatterplot-legend-server-spec.ts.
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

const MARKERS_ROW = 'prop-markers';
const MARKERS_COMBO = 'div-column-combobox-markers';

async function pickMarkerColumn(page: Page, column: string): Promise<void> {
  await sp.pickPanelColumn(page, MARKERS_ROW, MARKERS_COMBO, 'marker', column, 'markersColumnName');
  await settledLegend(page);
}

async function clearMarkerColumn(page: Page): Promise<void> {
  await sp.clearPanelColumn(page, MARKERS_ROW, MARKERS_COMBO, 'marker', 'markersColumnName');
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
  const root = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot')?.root as HTMLElement | undefined;
  const empty = {present: false, all: 0, coloring: 0, extra: 0, labels: [] as string[],
    colorLabels: [] as string[], glyphLabels: [] as string[], current: [] as string[], dimmed: [] as string[]};
  if (!root) return empty;
  // legendVisibility: Never sets display:none on the host rather than removing it (LegendHost.apply)
  const host = root.querySelector('[name="legend"]') as HTMLElement | null;
  const hostVisible = !!host && getComputedStyle(host).display !== 'none';
  const items = hostVisible ? [...root.querySelectorAll('[name="legend"] .d4-legend-item')] as HTMLElement[] : [];
  const txt = (i: Element) => (i.querySelector('.d4-legend-value')?.textContent ?? '').trim();
  const colorItems = items.filter((i) => !i.classList.contains('d4-legend-item-extra'));
  return {
    present: hostVisible,
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

const settledLegend = (page: Page): Promise<Legend> =>
  v.pollStable(() => readLegend(page), (a, b) => JSON.stringify(a) === JSON.stringify(b), 3000, 100);

interface Ink { nonWhite: number; saturated: number; pale: number; }

const readInk = (page: Page): Promise<Ink> => page.evaluate(() => {
  const root = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot').root as HTMLElement;
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

const sameInk = (a: Ink, b: Ink) =>
  a.nonWhite === b.nonWhite && a.saturated === b.saturated && a.pale === b.pale;

async function settledInk(page: Page): Promise<Ink> {
  await sp.parkPointer(page);
  return v.pollStable(() => readInk(page), sameInk, 5000, 100);
}

async function settledInkAfterChange(page: Page, from: Ink): Promise<Ink> {
  await sp.parkPointer(page);
  await v.pollValue(() => readInk(page), (cur) => !sameInk(cur, from), 6000, 100);
  return v.pollStable(() => readInk(page), sameInk, 4000, 100);
}

const legendSelection = async (page: Page) => (await readLegend(page)).current.join('|');

const legendFiltering = (page: Page) => page.evaluate(() => {
  const root = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot')?.root as HTMLElement | undefined;
  return !!root?.querySelector('[name="legend"]')?.classList.contains('d4-legend-filtering');
});

const legendPoint = (page: Page, label: string) => page.evaluate((l: string) => {
  const root = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot')?.root as HTMLElement | undefined;
  const item = [...(root?.querySelectorAll('[name="legend"] .d4-legend-item') ?? [])]
    .find((i) => (i.querySelector('.d4-legend-value')?.textContent ?? '').trim() === l);
  if (!item) return null;
  const b = (item.querySelector('.d4-legend-value') ?? item).getBoundingClientRect();
  return {x: b.x + b.width / 2, y: b.y + b.height / 2};
}, label);

async function clickLegendEntry(page: Page, label: string): Promise<void> {
  const before = await legendSelection(page);
  const pt = await legendPoint(page, label);
  expect(pt, `legend entry ${label}`).not.toBeNull();
  await page.mouse.click(pt!.x, pt!.y);
  await v.pollValue(() => legendSelection(page), (cur) => cur !== before, 4000, 50);
}

async function clearLegendSelection(page: Page): Promise<void> {
  if (!await legendFiltering(page)) return;
  const current = (await readLegend(page)).current;
  if (!current.length) return;
  const pt = await legendPoint(page, current[0]);
  if (!pt) return;
  await page.mouse.click(pt.x, pt.y);
  await v.pollValue(() => legendFiltering(page), (on) => !on, 4000, 50);
}

const viewerColumns = (page: Page) => page.evaluate(() => {
  const s = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return {color: s.props.colorColumnName, markers: s.props.markersColumnName};
});

const legendHostCount = (page: Page) => page.evaluate(() => {
  const root = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot')?.root as HTMLElement | undefined;
  return [...(root?.querySelectorAll('[name="legend"]') ?? [])]
    .filter((el) => getComputedStyle(el as HTMLElement).display !== 'none').length;
});

const legendVisibility = (page: Page) => sp.readProp(page, 'legendVisibility');

const VISIBILITY_MENU = ['div-Properties...', 'div-Properties...---Legend',
  'div-Properties...---Legend---Legend-Visibility'];

async function setLegendVisibility(page: Page, value: string): Promise<void> {
  let driven = false;
  const editor = '[name="prop-view-legend-visibility"]';
  await sp.openSettings(page);
  await sp.revealPropEditor(page, editor, 'legend');
  const cell = page.locator(editor);
  await cell.scrollIntoViewIfNeeded();
  await cell.click();
  const select = page.locator('[name="prop-legend-visibility"] select.property-grid-item-editor-spinner');
  await select.waitFor({timeout: 3000}).catch(() => {});
  if (await select.count() > 0) {
    const options = (await select.locator('option').allTextContents()).map((o) => o.trim());
    expect(options).toContain(value);
    await select.selectOption(value);
    await sp.propIs(page, 'legendVisibility', value, 4000);
    driven = true;
  }
  if (!driven) {
    const leafName = `${VISIBILITY_MENU[2]}---${value}`;
    await sp.clickContextMenuLeaf(page, [...VISIBILITY_MENU, leafName]);
    await sp.propIs(page, 'legendVisibility', value, 4000);
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
  await v.pollValue(() => sp.menuLeafNames(page), (names) => names.includes('div-Color-Coding'), 3000, 50);
}

async function clickHeaderMenuLeaf(page: Page, group: string, leaf: string): Promise<void> {
  await page.locator(`.d4-menu-popup [name="${group}"]`).last().hover();
  await page.locator(`.d4-menu-popup [name="${leaf}"]`).last().waitFor({timeout: 5000});
  await page.locator(`.d4-menu-popup [name="${leaf}"]`).last().click();
  await v.pollValue(() => page.locator('.d4-menu-popup:visible').count(), (n) => n === 0, 2000, 50);
}

const conditionalRanges = (page: Page, column: string) => page.evaluate((c: string) => {
  const col = grok.shell.tv.dataFrame.col(c);
  const raw = col.tags['.color-coding-conditional'];
  return {type: col.tags['.color-coding-type'] ?? null, ranges: raw ? Object.keys(JSON.parse(raw)) : []};
}, column);

test('Scatter Plot — Legend Lifecycle, Filter Interplay', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  const errors = sp.trackErrors(page);
  const errCount = errors.count;

  await openDatagrok(page);
  await v.openTable(page, {path: demogPath, semTypeTimeoutMs: 3000});
  await sp.addScatterPlot(page);
  await sp.pickOnViewer(page, 'x', 'WEIGHT');
  await sp.pickOnViewer(page, 'y', 'HEIGHT');
  const fullRowCount = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount as number);
  expect(fullRowCount).toBeGreaterThan(0);

  await softStep('Color legend and the joint Color plus Marker legend', async () => {
    const errBefore = errCount();

    await sp.pickOnViewer(page, 'color', 'RACE');
    const colorOnly = await settledLegend(page);
    expect(colorOnly.present).toBe(true);
    expect(colorOnly.colorLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
    expect(colorOnly.all).toBe(RACE_CATEGORIES.length);

    await pickMarkerColumn(page, 'SEX');
    const joint = await settledLegend(page);
    expect(joint.coloring).toBe(RACE_CATEGORIES.length);
    expect(joint.extra).toBe(SEX_CATEGORIES.length);
    expect(joint.glyphLabels.sort()).toEqual([...SEX_CATEGORIES].sort());
    expect(joint.all).toBe(RACE_CATEGORIES.length + SEX_CATEGORIES.length);

    await sp.pickOnViewer(page, 'color', 'AGE');
    const numericColor = await settledLegend(page);
    expect(numericColor.glyphLabels.sort()).toEqual([...SEX_CATEGORIES].sort());
    expect(numericColor.extra).toBe(SEX_CATEGORIES.length);

    await sp.pickOnViewer(page, 'color', 'RACE');
    const back = await settledLegend(page);
    expect(back.coloring).toBe(RACE_CATEGORIES.length);
    expect(back.extra).toBe(SEX_CATEGORIES.length);
    expect(back.all).toBe(RACE_CATEGORIES.length + SEX_CATEGORIES.length);
    expect(back.colorLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
    expect(new Set(back.labels).size).toBe(back.labels.length);

    await pickMarkerColumn(page, 'RACE');
    const same = await settledLegend(page);
    expect(same.colorLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
    expect(same.glyphLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
    expect(same.all).toBe(RACE_CATEGORIES.length);

    await pickMarkerColumn(page, 'SEX');
    expect((await viewerColumns(page)).markers).toBe('SEX');
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Clearing the Marker column removes its glyph entries', async () => {
    const errBefore = errCount();
    await sp.pickOnViewer(page, 'color', 'SEX');
    await pickMarkerColumn(page, 'SEX');
    const before = await settledLegend(page);
    expect(before.colorLabels.sort()).toEqual([...SEX_CATEGORIES].sort());
    expect(before.glyphLabels.sort()).toEqual([...SEX_CATEGORIES].sort());

    await clearMarkerColumn(page);
    const after = await v.pollValue(() => readLegend(page), (l) => l.glyphLabels.length === 0, 3000, 100);

    expect((await viewerColumns(page)).markers).toBe('');
    expect(after.glyphLabels).toEqual([]);
    expect(after.colorLabels.sort()).toEqual([...SEX_CATEGORIES].sort());
    expect(after.present).toBe(true);

    await sp.pickOnViewer(page, 'color', 'RACE');
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Filtered-out categories are absent from the marker legend', async () => {
    const errBefore = errCount();
    await sp.openFilterPanel(page);

    await v.applyCategoricalFilter(page, 'RACE', ['Asian', 'Caucasian'], 600);
    const narrowed = await sp.filterMoved(page, fullRowCount);
    expect(narrowed).toBeLessThan(fullRowCount);
    await pickMarkerColumn(page, 'RACE');
    const afterFilter = await readLegend(page);
    expect(afterFilter.glyphLabels.sort()).toEqual(['Asian', 'Caucasian']);

    await v.resetFilters(page);
    await sp.filterMoved(page, narrowed);
    await clearMarkerColumn(page);
    await pickMarkerColumn(page, 'RACE');
    const unfiltered = await readLegend(page);
    expect(unfiltered.glyphLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
    await v.applyCategoricalFilter(page, 'RACE', ['Asian', 'Caucasian'], 600);
    await sp.filterMoved(page, fullRowCount);
    const beforeFilter = await v.pollValue(() => readLegend(page), (l) => l.glyphLabels.length === 2, 3000, 100);
    expect(beforeFilter.glyphLabels.sort()).toEqual(['Asian', 'Caucasian']);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Filtered-out categories are absent from the marker legend — filter cleared', async () => {
    const before = await sp.filterCount(page);
    await v.resetFilters(page);
    const restored = await sp.filterMoved(page, before);
    expect(restored).toBe(fullRowCount);
    const legend = await v.pollValue(() => readLegend(page),
      (l) => l.glyphLabels.length === RACE_CATEGORIES.length, 3000, 100);
    expect(legend.glyphLabels.sort()).toEqual([...RACE_CATEGORIES].sort());
  });

  await softStep('Clicking a legend category hides the other categories on the plot', async () => {
    const errBefore = errCount();
    await sp.pickOnViewer(page, 'color', 'RACE');
    await pickMarkerColumn(page, 'SEX');
    try {
      const baseline = await settledInk(page);
      expect(baseline.nonWhite).toBeGreaterThan(0);
      const baseCount = await sp.filterHeld(page, 500);
      expect(baseCount).toBe(fullRowCount);

      await clickLegendEntry(page, 'Asian');
      const legendOnly = await settledInkAfterChange(page, baseline);
      expect(legendOnly.nonWhite)
        .toBeLessThan(Math.round(baseline.nonWhite * SELECTION_INK_MAX_FRACTION));
      expect(await legendFiltering(page)).toBe(true);
      const marked = await readLegend(page);
      expect(marked.current).toContain('Asian');
      expect(marked.dimmed.length).toBeGreaterThan(0);
      expect(await sp.filterHeld(page)).toBe(baseCount);

      await v.applyCategoricalFilter(page, 'SEX', ['F'], 600);
      const narrowed = await sp.filterMoved(page, fullRowCount);
      expect(narrowed).toBeLessThan(fullRowCount);
      const legendThenPanel = await settledInkAfterChange(page, legendOnly);
      expect(legendThenPanel.nonWhite)
        .toBeLessThan(legendOnly.nonWhite - PANEL_RESPECT_MARGIN_PX);
      expect((await readLegend(page)).labels).not.toContain('M');

      await clearLegendSelection(page);
      await v.resetFilters(page);
      await sp.filterMoved(page, narrowed);
      await v.applyCategoricalFilter(page, 'SEX', ['F'], 600);
      await sp.filterMoved(page, fullRowCount);
      const panelOnly = await settledInk(page);
      await clickLegendEntry(page, 'Asian');
      const panelThenLegend = await settledInkAfterChange(page, panelOnly);
      expect(panelThenLegend).toEqual(legendThenPanel);
      expect(errCount()).toBe(errBefore);
    } finally {
      await clearLegendSelection(page);
      await v.resetFilters(page);
      await sp.filterHeld(page, 500);
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
      const baseCount = await sp.filterHeld(page, 500);
      try {
        await page.evaluate((name: string) => {
          const df = grok.shell.tv.dataFrame;
          const values: (string | null)[] = [];
          for (let i = 0; i < df.rowCount; i++) values.push(i % 3 === 0 ? null : (i % 3 === 1 ? 'alpha' : 'beta'));
          df.columns.add(DG.Column.fromStrings(name, values));
        }, EMPTY_PROBE_COLUMN);
        await v.waitForViewerRendered(page, sp.SP_TYPE, 1200);
        await sp.pickOnViewer(page, 'color', EMPTY_PROBE_COLUMN);
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
        expect(await sp.filterHeld(page)).toBe(baseCount);

        await v.applyCategoricalFilter(page, 'SEX', ['F'], 600);
        expect(await sp.filterMoved(page, baseCount)).toBeLessThan(baseCount);
        const withPanel = await settledInkAfterChange(page, selected);
        expect(withPanel.nonWhite).toBeLessThan(selected.nonWhite - PANEL_RESPECT_MARGIN_PX);

        await clearLegendSelection(page);
        const filtered = await sp.filterCount(page);
        await v.resetFilters(page);
        await sp.filterMoved(page, filtered);
        expect(await settledInkAfterChange(page, withPanel)).toEqual(baseline);
        expect(errCount()).toBe(errBefore);
      } finally {
        await clearLegendSelection(page);
        await v.resetFilters(page);
        await page.evaluate((name: string) => {
          const df = grok.shell.tv.dataFrame;
          if (df.col(name)) df.columns.remove(name);
        }, EMPTY_PROBE_COLUMN);
        await sp.pickOnViewer(page, 'color', 'RACE');
        await pickMarkerColumn(page, 'SEX');
      }
      const peak = await viewerColumns(page);
      expect(peak.color).toBe('RACE');
      expect(peak.markers).toBe('SEX');
    });

  await softStep('A column with conditional color coding still produces a legend', async () => {
    const errBefore = errCount();
    try {
      await sp.addTableView(page, spgiPath, 'CAST Idea ID');
      await page.locator('.d4-grid[name="viewer-Grid"]').last().waitFor({timeout: 15_000});
      await v.pollValue(() => page.evaluate(() => {
        const grid = grok.shell.tv.grid;
        return grid.root.getBoundingClientRect().height > 0 && grid.columns.byName('CAST Idea ID')?.width > 0;
      }), (ok) => ok, 5000, 50);

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
        await v.pollValue(() => rangeInputs.count(), (n) => n < before, 1000, 50);
      }
      expect(await rangeInputs.count()).toBe(CONDITIONAL_RANGES.length);
      for (let i = 0; i < CONDITIONAL_RANGES.length; i++) {
        await rangeInputs.nth(i).click();
        await page.keyboard.press('Control+A');
        await page.keyboard.type(CONDITIONAL_RANGES[i]);
        await page.keyboard.press('Enter');
        await v.pollValue(() => rangeInputs.nth(i).inputValue(), (val) => val === CONDITIONAL_RANGES[i], 1000, 50);
      }
      await dialog.locator('[name="button-CLOSE"]').click();
      await dialog.waitFor({state: 'detached', timeout: 5000});

      const coding = await conditionalRanges(page, 'CAST Idea ID');
      expect(coding.type).toBe('Conditional');
      expect(coding.ranges).toEqual(CONDITIONAL_RANGES);

      await sp.addScatterPlot(page);
      const axes = await page.evaluate(() => {
        const s = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
        return {x: s.props.xColumnName, y: s.props.yColumnName};
      });
      expect(axes.x).toBe('CAST Idea ID');
      expect(axes.y).toBe('Idea ID');

      await sp.pickOnViewer(page, 'color', 'CAST Idea ID');
      const legend = await v.pollValue(() => readLegend(page), (l) => l.present && l.colorLabels.length > 0, 3000, 100);
      expect(legend.present).toBe(true);
      expect(legend.colorLabels).toEqual(CONDITIONAL_RANGES);
      expect(errCount()).toBe(errBefore);
    } finally {
      await page.evaluate(() => {
        grok.shell.v.close();
        for (const view of grok.shell.tableViews)
          if (view.dataFrame.columns.contains('RACE')) grok.shell.v = view;
      });
      await v.pollValue(() => page.evaluate(() => {
        const tv = grok.shell.tv;
        return !!tv?.dataFrame?.columns.contains('RACE') && tv.viewers.some((x: any) => x.type === 'Scatter plot');
      }), (ok) => ok, 3000, 50);
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
      const hidden = await v.pollValue(() => readLegend(page), (l) => l.present === false, 3000, 100);
      expect(hidden.present).toBe(false);

      await setLegendVisibility(page, 'Always');
      const shown = await v.pollValue(() => readLegend(page),
        (l) => l.present && l.all === before.all, 3000, 100);
      expect(shown.present).toBe(true);
      expect(shown.all).toBe(before.all);
      expect(shown.colorLabels.sort()).toEqual([...before.colorLabels].sort());
      expect(shown.glyphLabels.sort()).toEqual([...before.glyphLabels].sort());

      await setLegendVisibility(page, initial);
      const restored = await v.pollValue(() => readLegend(page),
        (l) => l.present && l.all === before.all, 3000, 100);
      expect(restored.present).toBe(true);
      expect(restored.all).toBe(before.all);
      expect(restored.colorLabels.sort()).toEqual([...before.colorLabels].sort());
      expect(errCount()).toBe(errBefore);
    } finally {
      await v.setViewerProps(page, sp.SP_TYPE, [{set: {legendVisibility: initial}, wait: 1200}]);
    }
  });

  await v.cleanupShell(page);
  v.finishSpec();
});
