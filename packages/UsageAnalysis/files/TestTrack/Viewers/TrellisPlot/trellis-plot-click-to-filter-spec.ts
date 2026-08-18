/* ---
realizes: [trellisplot.cp.click-to-filter, trellisplot.int.click-cell-filter-select-events, trellisplot.int.onclick-filter-panel-collab]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const PANEL_COLUMN = 'DIS_POP';

const ROW_SOURCE_OPTIONS = ['All', 'CurrentRow', 'Filtered', 'FilteredSelected',
  'MouseOverGroup', 'MouseOverRow', 'Selected', 'SelectedOrCurrent'];

const isBenignError = (text: string) =>
  /Failed to load resource/.test(text) || /404 \(\)/.test(text) || /favicon/.test(text) ||
  /Unable to find element in cloned iframe/.test(text);

async function cellIndexFor(page: Page, xValue: string, yValue: string): Promise<number> {
  return page.evaluate(({xValue, yValue}) => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    const df = grok.shell.tv.dataFrame;
    const xCol = df.col(tp.props.xColumnNames[0]);
    const yCol = df.col(tp.props.yColumnNames[0]);
    const xi = xCol.categories.indexOf(xValue);
    const yi = yCol.categories.indexOf(yValue);
    return yi * tp.xCategoriesCount + xi;
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

async function panelOnlyRowCount(page: Page, column: string, selected: string[]): Promise<number> {
  return page.evaluate(({column, selected}) => {
    const df = grok.shell.tv.dataFrame;
    const c = df.col(column);
    let n = 0;
    for (let i = 0; i < df.rowCount; i++)
      if (selected.indexOf(c.get(i)) >= 0) n++;
    return n;
  }, {column, selected});
}

async function waitForFilterCount(page: Page, expected: number, capMs = 3000): Promise<void> {
  await page.waitForFunction((exp) => grok.shell.tv.dataFrame.filter.trueCount === exp,
    expected, {timeout: capMs}).catch(() => {});
}
async function waitForFilterBelow(page: Page, ceiling: number, capMs = 3000): Promise<void> {
  await page.waitForFunction((c) => grok.shell.tv.dataFrame.filter.trueCount < c,
    ceiling, {timeout: capMs}).catch(() => {});
}
async function waitForSelectionCount(page: Page, expected: number, capMs = 3000): Promise<void> {
  await page.waitForFunction((exp) => grok.shell.tv.dataFrame.selection.trueCount === exp,
    expected, {timeout: capMs}).catch(() => {});
}

async function openTrellisPropertyGrid(page: Page, probeProp: string): Promise<void> {
  if (await page.locator(`.property-grid tr[name="prop-${probeProp}"]`).count() === 0)
    await v.openViewerGear(page, 'Trellis-plot');
  const tab = page.locator('.d4-tab-header[name="Trellis"]');
  if (await tab.count() > 0) {
    await tab.first().click();
    await page.waitForTimeout(300); 
  }
  await page.locator(`.property-grid tr[name="prop-${probeProp}"]`).first()
    .waitFor({state: 'attached', timeout: 15000});
}

async function propCategoryOf(page: Page, prop: string): Promise<string | undefined> {
  const cat = await page.evaluate((p) => {

    for (const grid of Array.from(document.querySelectorAll('.property-grid'))) {
      const rows = Array.from(grid.querySelectorAll('tr'));
      const i = rows.findIndex((r) => r.getAttribute('name') === `prop-${p}`);
      if (i < 0) continue;
      for (let j = i - 1; j >= 0; j--) {
        if (!rows[j].className.includes('property-grid-category')) continue;
        const name = rows[j].getAttribute('name') ?? '';

        return name.startsWith('prop-category-') ? name.substring('prop-category-'.length) : '';
      }
      return '';
    }
    return '';
  }, prop);
  return cat === '' ? undefined : cat;
}

async function setTrellisChoiceViaPanel(page: Page, prop: string, value: string): Promise<string> {
  await openTrellisPropertyGrid(page, prop);
  const category = await propCategoryOf(page, prop);
  if (category) await v.ensurePropertyCategory(page, 'Trellis-plot', category, prop);
  await v.selectPropertyGridChoice(page, prop, value, category);
  return v.propertyGridValue(page, prop, category);
}

async function propertyGridChoices(page: Page, prop: string): Promise<string[]> {
  await openTrellisPropertyGrid(page, prop);
  const category = await propCategoryOf(page, prop);

  if (category) await v.ensurePropertyCategory(page, 'Trellis-plot', category, prop);
  const row = page.locator(`.property-grid tr[name="prop-${prop}"]`).first();
  await row.locator('td').last().click();

  await row.locator('select').waitFor({state: 'attached', timeout: 3000});
  const options = await row.locator('select option').allTextContents();
  await page.keyboard.press('Escape');
  await page.waitForTimeout(300); 
  return options.map((s) => s.trim()).filter((s) => s.length > 0);
}

test('Trellis plot: click-to-filter, click-to-select, events, keyboard navigation', async ({page}) => {
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
    return {
      rowCount: df.rowCount,
      sex: df.col('SEX').categories.length,
      race: df.col('RACE').categories.length,
    };
  });
  expect(setup).toEqual({rowCount: 5850, sex: 2, race: 4});
  const fullRowCount = setup.rowCount;

  let panelSelected: string[] = [];

  await v.addViewerByIcon(page, 'trellis-plot', 'Trellis-plot', 15000);
  await page.evaluate(() => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    tp.props.xColumnNames = ['SEX'];
    tp.props.yColumnNames = ['RACE'];
    tp.props.viewerType = 'Scatter plot';
  });
  await v.waitForViewerRendered(page, 'Trellis plot', 900);
  const cellLocator = page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell');
  await expect(cellLocator).toHaveCount(8);

  await softStep('Scenario 1 Step 7', async () => {

    const idxA = await cellIndexFor(page, 'F', 'Caucasian');
    const idxB = await cellIndexFor(page, 'M', 'Caucasian');
    const result = await page.evaluate(async ({idxA, idxB}) => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;

      const rendered = (capMs: number) => new Promise<void>((resolve) => {
        let sub: any = null;
        try { sub = tp.onViewerRendered.subscribe(() => { sub.unsubscribe(); resolve(); }); }
        catch {  }
        setTimeout(() => { try { sub?.unsubscribe(); } catch {} resolve(); }, capMs);
      });
      function cellHash(cellIdx: number): number | null {
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
      const noneSettled = rendered(2500);
      tp.props.onClick = 'None';
      await noneSettled;
      const beforeA = cellHash(idxA), beforeB = cellHash(idxB);

      const filterSettled = rendered(2500);
      tp.props.onClick = 'Filter';
      await filterSettled;
      const afterA = cellHash(idxA), afterB = cellHash(idxB);
      return {
        read: beforeA !== null && beforeB !== null && afterA !== null && afterB !== null,

        distinct: beforeA !== beforeB,
        unchangedA: beforeA === afterA,
        unchangedB: beforeB === afterB,
      };
    }, {idxA, idxB});
    expect(result.read).toBe(true);
    expect(result.distinct).toBe(true);
    expect(result.unchangedA).toBe(true);
    expect(result.unchangedB).toBe(true);
  });

  await softStep('Scenario 1 Step 5', async () => {
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true); df.rows.requestFilter();
    });
    const expected = await comboRowCount(page, 'SEX', 'F', 'RACE', 'Caucasian');
    const idx = await cellIndexFor(page, 'F', 'Caucasian');

    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});

    await waitForFilterCount(page, expected, 3000);
    const after = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const filters: string[] = [];
      for (const f of df.rows.filters) filters.push(f);
      return {filterCount: df.filter.trueCount, filters};
    });
    expect(after.filterCount).toBe(expected);
    expect(after.filterCount).toBeLessThan(fullRowCount);
    expect(after.filters).toContain('SEX: F');
    expect(after.filters).toContain('RACE: Caucasian');
  });

  await softStep('Scenario 1 Step 6', async () => {
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true); df.rows.requestFilter(); 
    });

    const idx = await cellIndexFor(page, 'F', 'Caucasian');
    await page.evaluate((idx) => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      (window as any).__ev = {inner: null as any, cell: null as any};
      (window as any).__subs = [
        tp.onEvent('d4-trellis-plot-inner-viewer-clicked').subscribe((a: any) => {
          const arg = (a && a.args !== undefined) ? a.args : a;
          (window as any).__ev.inner = arg;
        }),
        tp.onEvent('d4-trellis-plot-current-cell-changed').subscribe((a: any) => {
          const mc = (a && a.args && a.args.matchCondition) ? a.args.matchCondition
            : (a && a.matchCondition ? a.matchCondition : a);
          (window as any).__ev.cell = mc;
        }),
      ];
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const cell = root.querySelectorAll('.d4-trellis-plot-cell')[idx];
      const cv = cell.querySelector('canvas') as HTMLCanvasElement;
      const r = cv.getBoundingClientRect();
      cv.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, cancelable: true, button: 0,
        clientX: r.left + r.width / 2, clientY: r.top + r.height / 2}));
    }, idx);

    await page.waitForFunction(() => (window as any).__ev?.inner != null, null, {timeout: 3000}).catch(() => {});
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});

    await page.waitForFunction(() => (window as any).__ev?.cell != null, null, {timeout: 3000}).catch(() => {});
    const ev = await page.evaluate(() => {
      ((window as any).__subs || []).forEach((s: any) => s?.unsubscribe?.());
      return (window as any).__ev;
    });

    expect(Array.isArray(ev.inner)).toBe(true);
    expect(ev.inner.length).toBeGreaterThan(0);
    expect(ev.cell).toMatchObject({SEX: 'F', RACE: 'Caucasian'});
  });

  await softStep('Scenario 1 Step 9', async () => {

    const seedIdx = await cellIndexFor(page, 'F', 'Caucasian');
    await cellLocator.nth(seedIdx).click({position: {x: 6, y: 6}});
    await waitForFilterBelow(page, fullRowCount, 3000);
    await page.keyboard.press('Escape');

    await waitForFilterCount(page, fullRowCount, 3000);
    const cleared = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const filters: string[] = [];
      for (const f of df.rows.filters) filters.push(f);
      return {trueCount: df.filter.trueCount,
        trellisEntries: filters.filter((s) => s.startsWith('SEX:') || s.startsWith('RACE:')).length};
    });

    expect(cleared).toEqual({trueCount: fullRowCount, trellisEntries: 0});

    await v.openFilterPanel(page);
    const cats: string[] = await page.evaluate((col) => grok.shell.tv.dataFrame.col(col).categories, PANEL_COLUMN);
    expect(cats.length).toBeGreaterThan(1);

    panelSelected = cats.filter((c) => c !== cats[0]);
    const expectedPanelOnly = await panelOnlyRowCount(page, PANEL_COLUMN, panelSelected);
    const {filteredCount: panelOnly} = await v.applyCategoricalFilter(page, PANEL_COLUMN, panelSelected);
    expect(panelOnly).toBe(expectedPanelOnly);
    expect(panelOnly).toBeLessThan(fullRowCount);
    expect(panelOnly).toBeGreaterThan(0);

    await expect(cellLocator).toHaveCount(8);

    const geom = await page.evaluate(({col, sel}) => {
      const df = grok.shell.tv.dataFrame;
      const s = df.col('SEX'), r = df.col('RACE'), p = df.col(col);
      let mAsian = 0, mAsianPanel = 0;
      for (let i = 0; i < df.rowCount; i++)
        if (s.get(i) === 'M' && r.get(i) === 'Asian') { mAsian++; if (sel.indexOf(p.get(i)) >= 0) mAsianPanel++; }
      return {mAsian, mAsianPanel};
    }, {col: PANEL_COLUMN, sel: panelSelected});
    expect(geom.mAsianPanel).toBeGreaterThan(0);
    const idx = await cellIndexFor(page, 'M', 'Asian');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await waitForFilterCount(page, geom.mAsianPanel, 3000);
    const after = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);

    expect(after).toBe(geom.mAsianPanel);
    expect(after).toBeLessThan(panelOnly);
    expect(after).toBeLessThan(geom.mAsian);
  });

  await softStep('Scenario 1 Step 10', async () => {
    const panelOnly = await panelOnlyRowCount(page, PANEL_COLUMN, panelSelected);

    const idx = await cellIndexFor(page, 'M', 'Asian');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await waitForFilterBelow(page, panelOnly, 3000);

    const beforeEsc = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(beforeEsc).toBeLessThan(panelOnly);
    await page.keyboard.press('Escape');
    await waitForFilterCount(page, panelOnly, 3000);
    const after = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const filters: string[] = [];
      for (const f of df.rows.filters) filters.push(f);
      return {filterCount: df.filter.trueCount, hasSex: filters.some((s) => s.startsWith('SEX:')),
        hasRace: filters.some((s) => s.startsWith('RACE:'))};
    });
    expect(after.filterCount).toBe(panelOnly);
    expect(after.hasSex).toBe(false);
    expect(after.hasRace).toBe(false);
  });

  await softStep('Scenario 1 Step 11', async () => {
    const panelOnly = await panelOnlyRowCount(page, PANEL_COLUMN, panelSelected);
    const idx = await cellIndexFor(page, 'M', 'Asian');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await waitForFilterBelow(page, panelOnly, 3000);

    const beforeChange = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(beforeChange).toBeLessThan(panelOnly);
    await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['CONTROL'];
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    await waitForFilterCount(page, panelOnly, 3000);
    const after = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(after).toBe(panelOnly);

    await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['SEX'];
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
  });

  await softStep('Scenario 1 Step 12', async () => {

    const before = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(before).toBeLessThan(fullRowCount);

    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) fg.remove(f);
    });

    await page.waitForFunction((col) => {
      const filters: string[] = [];
      for (const f of grok.shell.tv.dataFrame.rows.filters) filters.push(f);
      return filters.filter((s) => s.startsWith(col + ':')).length === 0;
    }, PANEL_COLUMN, {timeout: 3000}).catch(() => {});
    const after = await page.evaluate(() => {
      const tv = grok.shell.tv;
      tv.dataFrame.rows.requestFilter();
      const filters: string[] = [];
      for (const f of tv.dataFrame.rows.filters) filters.push(f);
      return {filterCount: tv.dataFrame.filter.trueCount, filters};
    });
    expect(after.filterCount).toBe(fullRowCount);
    expect(after.filters.filter((s) => s.startsWith(PANEL_COLUMN + ':'))).toEqual([]);
  });

  await softStep('Scenario 2 Step 3', async () => {

    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true); df.selection.setAll(false); df.rows.requestFilter();
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.onClick = 'None';
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);

    const parked = await setTrellisChoiceViaPanel(page, 'row-source', 'Filtered');
    expect(parked).toBe('Filtered');
    const beforeCommit = await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      return tp.props.rowSource;
    });
    expect(beforeCommit).toBe('Filtered');

    const shown = await setTrellisChoiceViaPanel(page, 'on-click', 'Filter');

    expect(shown).toBe('Filter');
    const committed = await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      return {onClick: tp.props.onClick, rowSource: tp.props.rowSource};
    });

    expect(committed.onClick).toBe('Filter');
    expect(committed.rowSource).toBe('All');

    const expected = await comboRowCount(page, 'SEX', 'F', 'RACE', 'Black');
    const idx = await cellIndexFor(page, 'F', 'Black');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await waitForFilterCount(page, expected, 3000);
    const after = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const filters: string[] = [];
      for (const f of df.rows.filters) filters.push(f);
      return {filterCount: df.filter.trueCount, filters};
    });
    expect(after.filterCount).toBe(expected);
    expect(after.filters).toContain('SEX: F');
    expect(after.filters).toContain('RACE: Black');

    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.onClick = 'None';
      df.filter.setAll(true); df.rows.requestFilter();
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
  });

  await softStep('Scenario 2 Step 5', async () => {
    const errBefore = consoleErrors.length;
    const pageErrBefore = pageErrors.length;

    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      df.filter.setAll(true); df.rows.requestFilter();
      tp.props.packCategories = false;
      df.selection.init((i: number) => i % 2 === 0);
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    await setTrellisChoiceViaPanel(page, 'on-click', 'Filter');

    await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['CONTROL'];
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['SEX'];
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);

    await page.evaluate(() => grok.shell.tv.dataFrame.rows.requestFilter());
    const start = await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      const df = grok.shell.tv.dataFrame;
      const filters: string[] = [];
      for (const f of df.rows.filters) filters.push(f);
      return {onClick: tp.props.onClick, trueCount: df.filter.trueCount,
        trellisEntries: filters.filter((s) => s.startsWith('SEX:') || s.startsWith('RACE:')).length};
    });

    expect(start.onClick).toBe('Filter');
    expect(start.trueCount).toBe(fullRowCount);
    expect(start.trellisEntries).toBe(0);

    const choices = await propertyGridChoices(page, 'row-source');

    expect([...choices].sort()).toEqual([...ROW_SOURCE_OPTIONS].sort());

    const rungs: {source: string; shown: string; rowSource: string; onClick: string;
      cells: number; trueCount: number; trellisEntries: number; mapped: boolean; wide: boolean;
      fedCells: number; blankFedCells: number[]; distinctHashes: number; signature: string}[] = [];
    for (const source of choices) {
      const shown = await setTrellisChoiceViaPanel(page, 'row-source', source);

      await v.waitForViewerRendered(page, 'Trellis plot', 900);
      const state = await page.evaluate(async () => {
        const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
        const df = grok.shell.tv.dataFrame;
        const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
        const filters: string[] = [];
        for (const f of df.rows.filters) filters.push(f);

        function cellProbe(cell: Element): {hash: number | null; painted: boolean} {
          const cv = cell.querySelector('canvas') as HTMLCanvasElement | null;
          if (!cv) return {hash: null, painted: false};
          try {
            const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
            let h = 0, first = -1, painted = false;
            for (let i = 0; i < img.length; i += 4) {
              const rgb = (img[i] << 16) | (img[i + 1] << 8) | img[i + 2];
              h = (h * 31 + rgb) % 2147483647;
              if (first === -1) first = rgb;
              else if (rgb !== first) painted = true;
            }
            return {hash: h, painted};
          } catch { return {hash: null, painted: false}; }
        }

        const rowCount = df.rowCount;
        const maskOf = (bs: any) => {
          const m = new Uint8Array(rowCount);
          const idxs = Array.from(bs.getSelectedIndexes()) as number[];
          for (const i of idxs) m[i] = 1;
          return m;
        };
        const filt = maskOf(df.filter), sel = maskOf(df.selection);
        const cur = df.currentRowIdx, mouse = df.mouseOverRowIdx;
        const sources: Record<string, (i: number) => boolean> = {
          All: () => true,
          Filtered: (i) => filt[i] === 1,
          Selected: (i) => sel[i] === 1,
          FilteredSelected: (i) => filt[i] === 1 && sel[i] === 1,
          SelectedOrCurrent: (i) => sel[i] === 1 || i === cur,
          CurrentRow: (i) => i === cur,
          MouseOverRow: (i) => i === mouse,

          MouseOverGroup: () => false,
        };
        const pred = sources[tp.props.rowSource];

        const cells = root.querySelectorAll('.d4-trellis-plot-cell');
        const probes = Array.from(cells).map(cellProbe);
        const xCol = df.col(tp.props.xColumnNames[0]), yCol = df.col(tp.props.yColumnNames[0]);
        const xAt: Record<string, number> = {}, yAt: Record<string, number> = {};
        (xCol.categories as string[]).forEach((c: string, i: number) => xAt[c] = i);
        (yCol.categories as string[]).forEach((c: string, i: number) => yAt[c] = i);
        const fedRows: number[] = new Array(cells.length).fill(0);
        if (pred)
          for (let i = 0; i < rowCount; i++) {
            if (!pred(i)) continue;
            const cx = xAt[xCol.get(i)], cy = yAt[yCol.get(i)];
            if (cx === undefined || cy === undefined) continue;
            const idx = cy * tp.xCategoriesCount + cx;
            if (idx < fedRows.length) fedRows[idx]++;
          }
        const fed: number[] = [];
        for (let i = 0; i < fedRows.length; i++) if (fedRows[i] > 0) fed.push(i);

        return {
          rowSource: tp.props.rowSource, onClick: tp.props.onClick,
          cells: cells.length,
          trueCount: df.filter.trueCount,
          trellisEntries: filters.filter((s) => s.startsWith('SEX:') || s.startsWith('RACE:')).length,
          mapped: !!pred,

          wide: !!pred && fed.length >= 2,
          fedCells: fed.length,
          blankFedCells: fed.filter((i) => !probes[i].painted),
          distinctHashes: Array.from(new Set(fed.map((i) => probes[i].hash))).length,
          signature: JSON.stringify(probes.map((p) => p.hash)),
        };
      });
      rungs.push({source, shown, ...state});
    }
    expect(consoleErrors.slice(errBefore)).toEqual([]);
    expect(pageErrors.slice(pageErrBefore)).toEqual([]);
    for (const r of rungs) {

      expect(r.shown).toBe(r.source);

      expect(r.onClick === 'Filter' && r.rowSource === 'Filtered').toBe(false);

      expect(r.cells).toBe(8);

      expect(r.trellisEntries).toBe(0);
      expect(r.trueCount).toBe(fullRowCount);
    }

    const filteredIdx = rungs.findIndex((r) => r.rowSource === 'Filtered');
    expect(filteredIdx).toBeGreaterThanOrEqual(0);

    const beforeFiltered = filteredIdx > 0 ? rungs[filteredIdx - 1].onClick : start.onClick;
    expect(beforeFiltered).toBe('Filter');
    expect(rungs[filteredIdx].onClick).toBe('None');

    expect(rungs.filter((r) => !r.mapped).map((r) => r.rowSource)).toEqual([]);

    expect(rungs.some((r) => r.rowSource === 'All' && r.wide)).toBe(true);
    const wideRungs = rungs.filter((r) => r.wide);

    expect(wideRungs.filter((r) => r.distinctHashes < 2)
      .map((r) => `${r.rowSource}: ${r.distinctHashes} hash(es) over ${r.fedCells} fed cells`)).toEqual([]);

    expect(wideRungs.filter((r) => r.blankFedCells.length > 0)
      .map((r) => `${r.rowSource}: ${r.blankFedCells.join(',')}`)).toEqual([]);

    expect(Array.from(new Set(rungs.map((r) => r.signature))).length).toBeGreaterThanOrEqual(2);

    await setTrellisChoiceViaPanel(page, 'row-source', rungs.find((r) => r.rowSource === 'All')?.source ?? 'All');
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.onClick = 'None';
      tp.props.packCategories = true;
      df.selection.setAll(false);
      df.filter.setAll(true); df.rows.requestFilter();
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    const restored = await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      return {rowSource: tp.props.rowSource, onClick: tp.props.onClick};
    });
    expect(restored).toEqual({rowSource: 'All', onClick: 'None'});
    await expect(cellLocator).toHaveCount(8);
  });

  await softStep('Scenario 3 Step 4', async () => {
    const viewportBefore = await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      const before = {xCat: tp.xCategoriesCount, yCat: tp.yCategoriesCount,
        cells: document.querySelectorAll('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell').length};
      tp.props.onClick = 'Select';
      return before;
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    const expected = await comboRowCount(page, 'SEX', 'F', 'RACE', 'Caucasian');
    const idx = await cellIndexFor(page, 'F', 'Caucasian');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});

    await waitForSelectionCount(page, expected, 3000);
    const after = await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      return {sel: grok.shell.tv.dataFrame.selection.trueCount,
        xCat: tp.xCategoriesCount, yCat: tp.yCategoriesCount,
        cells: document.querySelectorAll('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell').length};
    });
    expect(after.sel).toBe(expected);

    expect({xCat: after.xCat, yCat: after.yCat, cells: after.cells})
      .toEqual({xCat: viewportBefore.xCat, yCat: viewportBefore.yCat, cells: viewportBefore.cells});
  });

  await softStep('Scenario 3 Step 6', async () => {

    const fCaucasian = await comboRowCount(page, 'SEX', 'F', 'RACE', 'Caucasian');
    const before = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;

      df.selection.setAll(false);
      const s = df.col('SEX'), r = df.col('RACE');
      for (let i = 0; i < df.rowCount; i++)
        if (s.get(i) === 'M' && r.get(i) === 'Black') df.selection.set(i, true, false);
      return df.selection.trueCount; 
    });
    expect(before).toBeGreaterThan(0);

    const {filteredCount} = await v.applyCategoricalFilter(page, 'SEX', ['M']);
    await v.openFilterPanel(page);
    expect(filteredCount).toBeLessThan(fullRowCount);
    expect(filteredCount).toBeGreaterThan(0);

    await expect(cellLocator).toHaveCount(8);
    const idx = await cellIndexFor(page, 'F', 'Caucasian');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}, modifiers: ['Control']});
    await waitForSelectionCount(page, before + fCaucasian, 3000);
    const after = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(after).toBe(before + fCaucasian);

    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) fg.remove(f);
    });
    await waitForFilterCount(page, fullRowCount, 3000);
    const cleared = await page.evaluate(() => {
      grok.shell.tv.dataFrame.rows.requestFilter();
      return grok.shell.tv.dataFrame.filter.trueCount;
    });
    expect(cleared).toBe(fullRowCount);
  });

  await softStep('Scenario 3 Step 7', async () => {

    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      df.selection.setAll(false);
      tp.props.packCategories = false;
      tp.props.xColumnNames = ['RACE'];
      tp.props.yColumnNames = ['SEVERITY'];
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    const before = await page.evaluate(() => {

      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const cells = root.querySelectorAll('.d4-trellis-plot-cell');
      let popIdx = -1;
      for (let i = 0; i < cells.length; i++) if (cells[i].querySelector('canvas')) { popIdx = i; break; }
      return {popIdx};
    });
    await cellLocator.nth(before.popIdx).click({position: {x: 6, y: 6}});

    await page.waitForFunction(() => grok.shell.tv.dataFrame.selection.trueCount > 0,
      null, {timeout: 3000}).catch(() => {});
    const emptyIdx = await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      const df = grok.shell.tv.dataFrame;
      const xCol = df.col(tp.props.xColumnNames[0]);
      const yCol = df.col(tp.props.yColumnNames[0]);
      const xi = xCol.categories.indexOf('Asian');
      const yi = yCol.categories.indexOf('Critical');
      return yi * tp.xCategoriesCount + xi;
    });
    const selBefore = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(selBefore).toBeGreaterThan(0);

    expect(emptyIdx).not.toBe(before.popIdx);
    await cellLocator.nth(emptyIdx).click({position: {x: 6, y: 6}, modifiers: ['Control']});

    await page.waitForFunction((idx) => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const cells = Array.from(root.querySelectorAll('.d4-trellis-plot-cell'));
      const marked = cells.filter((c) => c.classList.contains('d4-trellis-cell-current'));
      return marked.length === 1 && cells.indexOf(marked[0]) === idx;
    }, emptyIdx, {timeout: 3000}).catch(() => {});
    const after = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const cells = Array.from(root.querySelectorAll('.d4-trellis-plot-cell'));
      const marked = cells.filter((c) => c.classList.contains('d4-trellis-cell-current'));
      return {sel: grok.shell.tv.dataFrame.selection.trueCount,
        currentCount: marked.length, currentIdx: cells.indexOf(marked[0])};
    });
    expect(after.sel).toBe(selBefore);

    expect({count: after.currentCount, idx: after.currentIdx}).toEqual({count: 1, idx: emptyIdx});
  });

  await softStep('Scenario 3 Step 8', async () => {

    const target = await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      const df = grok.shell.tv.dataFrame;
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const cells = root.querySelectorAll('.d4-trellis-plot-cell');
      const xCol = df.col(tp.props.xColumnNames[0]);
      const yCol = df.col(tp.props.yColumnNames[0]);
      const xAt: Record<string, number> = {}, yAt: Record<string, number> = {};
      (xCol.categories as string[]).forEach((c: string, i: number) => xAt[c] = i);
      (yCol.categories as string[]).forEach((c: string, i: number) => yAt[c] = i);
      const selMask = new Uint8Array(df.rowCount);
      for (const i of Array.from(df.selection.getSelectedIndexes()) as number[]) selMask[i] = 1;
      const total: number[] = new Array(cells.length).fill(0);
      const selected: number[] = new Array(cells.length).fill(0);
      for (let i = 0; i < df.rowCount; i++) {
        const cx = xAt[xCol.get(i)], cy = yAt[yCol.get(i)];
        if (cx === undefined || cy === undefined) continue;
        const k = cy * tp.xCategoriesCount + cx;
        if (k >= total.length) continue;
        total[k]++;
        if (selMask[i]) selected[k]++;
      }
      const before = df.selection.trueCount;
      const seedIdx = selected.findIndex((n) => n > 0);
      let idx = -1;
      for (let i = 0; i < cells.length; i++) {
        if (i === seedIdx || selected[i] > 0) continue;
        if (!cells[i].querySelector('canvas') || total[i] === 0 || total[i] === before) continue;
        idx = i; break;
      }
      return {idx, seedIdx, before, expected: idx >= 0 ? total[idx] : -1};
    });

    expect(target.before).toBeGreaterThan(0);
    expect(target.idx).toBeGreaterThanOrEqual(0);
    expect(target.idx).not.toBe(target.seedIdx);
    expect(target.expected).not.toBe(target.before);
    await cellLocator.nth(target.idx).click({position: {x: 6, y: 6}});

    await waitForSelectionCount(page, target.expected, 3000);
    const after = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(after).toBe(target.expected);

    expect(after).not.toBe(target.before);
    expect(after).not.toBe(target.before + target.expected);

    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      df.selection.setAll(false);
      tp.props.packCategories = true;
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      tp.props.onClick = 'None';
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    await expect(cellLocator).toHaveCount(8);
  });

  await softStep('Scenario 4 Step 3', async () => {

    const idx = await cellIndexFor(page, 'F', 'Asian');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});

    await page.waitForFunction((i) => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const cells = Array.from(root.querySelectorAll('.d4-trellis-plot-cell'));
      const marked = cells.filter((c) => c.classList.contains('d4-trellis-cell-current'));
      return marked.length === 1 && cells.indexOf(marked[0]) === i;
    }, idx, {timeout: 3000}).catch(() => {});

    const currentBefore = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const cells = Array.from(root.querySelectorAll('.d4-trellis-plot-cell'));
      const marked = cells.filter((c) => c.classList.contains('d4-trellis-cell-current'));
      return {count: marked.length, idx: cells.indexOf(marked[0])};
    });
    expect(currentBefore).toEqual({count: 1, idx});
    await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      (window as any).__cc = [];
      (window as any).__ccSub = tp.onEvent('d4-trellis-plot-current-cell-changed').subscribe((a: any) => {
        const mc = (a && a.args && a.args.matchCondition) ? a.args.matchCondition
          : (a && a.matchCondition ? a.matchCondition : a);
        (window as any).__cc.push(mc);
      });
    });

    const awaitCc = (n: number) => page.waitForFunction((c) => (window as any).__cc.length >= c,
      n, {timeout: 3000}).catch(() => {});
    await page.keyboard.press('ArrowRight');
    await awaitCc(1);
    await page.keyboard.press('ArrowDown');
    await awaitCc(2);
    await page.keyboard.press('ArrowLeft');
    await awaitCc(3);
    const events = await page.evaluate(() => {
      (window as any).__ccSub?.unsubscribe?.();
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const cells = Array.from(root.querySelectorAll('.d4-trellis-plot-cell'));
      const marked = cells.filter((c) => c.classList.contains('d4-trellis-cell-current'));
      return {cc: (window as any).__cc, currentCount: marked.length, currentIdx: cells.indexOf(marked[0])};
    });
    expect(events.cc.length).toBeGreaterThanOrEqual(3);
    for (const mc of events.cc) {
      expect(mc).toHaveProperty('SEX');
      expect(mc).toHaveProperty('RACE');
    }

    const combos = Array.from(new Set((events.cc as any[])
      .map((mc) => JSON.stringify([mc.SEX, mc.RACE]))));
    expect(combos.length).toBeGreaterThanOrEqual(3);

    expect(events.currentCount).toBe(1);
    expect(events.currentIdx).not.toBe(currentBefore.idx);
  });

  await softStep('Scenario 4 Step 4', async () => {

    const idx = await cellIndexFor(page, 'F', 'Asian');
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true); df.rows.requestFilter();
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.onClick = 'None';
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});

    await page.waitForFunction((i) => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const cells = Array.from(root.querySelectorAll('.d4-trellis-plot-cell'));
      const marked = cells.filter((c) => c.classList.contains('d4-trellis-cell-current'));
      return marked.length === 1 && cells.indexOf(marked[0]) === i;
    }, idx, {timeout: 3000}).catch(() => {});
    await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.onClick = 'Filter';
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);

    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await waitForFilterBelow(page, fullRowCount, 3000);
    const beforeArrow = await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      (window as any).__ccMove = [];
      (window as any).__ccMoveSub = tp.onEvent('d4-trellis-plot-current-cell-changed').subscribe((a: any) => {
        const mc = (a && a.args && a.args.matchCondition) ? a.args.matchCondition
          : (a && a.matchCondition ? a.matchCondition : a);
        (window as any).__ccMove.push(mc);
      });
      return grok.shell.tv.dataFrame.filter.trueCount;
    });

    expect(beforeArrow).toBeLessThan(fullRowCount);
    await page.keyboard.press('ArrowDown');

    await page.waitForFunction(() => (window as any).__ccMove.length > 0, null, {timeout: 3000}).catch(() => {});
    const mc = await page.evaluate(() => {
      (window as any).__ccMoveSub?.unsubscribe?.();
      const cc = (window as any).__ccMove;
      return cc.length > 0 ? cc[cc.length - 1] : null;
    });

    expect(typeof mc?.SEX).toBe('string');
    expect(typeof mc?.RACE).toBe('string');

    expect({SEX: mc.SEX, RACE: mc.RACE}).not.toEqual({SEX: 'F', RACE: 'Asian'});
    const expectedAfterArrow = await comboRowCount(page, 'SEX', mc.SEX, 'RACE', mc.RACE);
    expect(expectedAfterArrow).toBeGreaterThan(0);

    await waitForFilterCount(page, expectedAfterArrow, 3000);
    const movedFilterCount = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(movedFilterCount).toBe(expectedAfterArrow);
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.onClick = 'None';
      df.filter.setAll(true); df.rows.requestFilter();
    });
    await v.waitForViewerRendered(page, 'Trellis plot', 900);
  });

  v.finishSpec();
});
