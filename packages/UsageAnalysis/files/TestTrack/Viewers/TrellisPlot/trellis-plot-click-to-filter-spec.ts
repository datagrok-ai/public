/* ---
realizes: [trellisplot.cp.click-to-filter, trellisplot.int.click-cell-filter-select-events, trellisplot.int.onclick-filter-panel-collab]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

// Must stay INDEPENDENT of the trellis split (SEX x RACE): a panel filter on a split
// column is swallowed by the cell constraint, and the AND below could then hold with one
// of the two contributors silently dropped.
const PANEL_COLUMN = 'DIS_POP';

// The eight Row Source options the property-grid `<select>` offers [DOM 2026-08-12]: the RAW
// RowSet enum names (ddt utils.dart:466-475), NOT the humanised context-menu captions the
// refdoc pins. Compared as a SET (both sides sorted): the option ORDER is not guaranteed.
const ROW_SOURCE_OPTIONS = ['All', 'CurrentRow', 'Filtered', 'FilteredSelected',
  'MouseOverGroup', 'MouseOverRow', 'Selected', 'SelectedOrCurrent'];

// Ambient noise the packed dev build emits regardless of what is actuated. Nothing
// product-specific is filtered out, so a NullError or a Dart stack trace still counts.
const isBenignError = (text: string) =>
  /Failed to load resource/.test(text) || /404 \(\)/.test(text) || /favicon/.test(text) ||
  /Unable to find element in cloned iframe/.test(text);

// Row-major grid (X categories = columns, Y = rows). The category order is read LIVE: a
// hardcoded index silently clicks the wrong combination once categories reorder.
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

// Row count of a column-value combination over the WHOLE dataset — deliberately
// blind to df.filter, so it stays a stable baseline the filter asserts compare against.
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

// Rows surviving the Filter Panel constraint ALONE — like comboRowCount, computed
// outside df.filter so the "panel-only" asserts never measure the thing under test.
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

// Bring the trellis's OWN property grid on screen: gear, then the `Trellis` tab — the tab
// named after the inner viewer type does not carry the trellis rows [refdoc trellis_plot.md
// HTML Structure]. Waits for `attached`, not `visible`: a collapsed row has zero size.
async function openTrellisPropertyGrid(page: Page, probeProp: string): Promise<void> {
  if (await page.locator(`.property-grid tr[name="prop-${probeProp}"]`).count() === 0)
    await v.openViewerGear(page, 'Trellis-plot');
  const tab = page.locator('.d4-tab-header[name="Trellis"]');
  if (await tab.count() > 0) {
    await tab.first().click();
    await page.waitForTimeout(600);
  }
  await page.locator(`.property-grid tr[name="prop-${probeProp}"]`).first()
    .waitFor({state: 'attached', timeout: 15000});
}

// The owning category is READ off the live grid, by CLASS (`property-grid-category`) and never
// by a `prop-category-` NAME prefix: the ORDINARY row `prop-category-label-font` sits above
// `prop-on-click` [DOM 2026-08-12], and `ensurePropertyCategory` hangs on that non-category.
async function propCategoryOf(page: Page, prop: string): Promise<string | undefined> {
  const cat = await page.evaluate((p) => {
    // Three `.property-grid` tables are on screen at once — two from the inner
    // ScatterPlotLook, one from TrellisPlotLook [DOM 2026-08-12] — so the walk is scoped
    // to the grid that OWNS the row and can never run off into a neighbouring one.
    for (const grid of Array.from(document.querySelectorAll('.property-grid'))) {
      const rows = Array.from(grid.querySelectorAll('tr'));
      const i = rows.findIndex((r) => r.getAttribute('name') === `prop-${p}`);
      if (i < 0) continue;
      for (let j = i - 1; j >= 0; j--) {
        if (!rows[j].className.includes('property-grid-category')) continue;
        const name = rows[j].getAttribute('name') ?? '';
        // A header not in the `prop-category-<id>` form cannot be addressed by
        // `ensurePropertyCategory`, so report "no category" rather than hand it a
        // selector that will never match.
        return name.startsWith('prop-category-') ? name.substring('prop-category-'.length) : '';
      }
      return '';
    }
    return '';
  }, prop);
  return cat === '' ? undefined : cat;
}

// Drive a trellis choice property through the REAL Context Panel editor and return the value
// the panel then DISPLAYS — the second actuation path of GROK-17711. A scenario using it must
// never fall back to `tp.props.<name> = ...`, or both paths collapse into one channel.
async function setTrellisChoiceViaPanel(page: Page, prop: string, value: string): Promise<string> {
  await openTrellisPropertyGrid(page, prop);
  const category = await propCategoryOf(page, prop);
  if (category) await v.ensurePropertyCategory(page, 'Trellis-plot', category, prop);
  await v.selectPropertyGridChoice(page, prop, value, category);
  return v.propertyGridValue(page, prop, category);
}

// The values the live choice editor OFFERS, read as `<select>` option TEXT, so a drifted
// spelling surfaces as the product statement it is rather than as a missing-option throw.
// Mind the surface: PROPERTY GRID (raw enum names), not the context menu (humanised captions).
async function propertyGridChoices(page: Page, prop: string): Promise<string[]> {
  await openTrellisPropertyGrid(page, prop);
  const category = await propCategoryOf(page, prop);
  // Expand first, like the setter: a row in a collapsed category has zero size and the
  // click below would hang instead of opening the editor.
  if (category) await v.ensurePropertyCategory(page, 'Trellis-plot', category, prop);
  const row = page.locator(`.property-grid tr[name="prop-${prop}"]`).first();
  await row.locator('td').last().click();
  await page.waitForTimeout(400);
  const options = await row.locator('select option').allTextContents();
  await page.keyboard.press('Escape');
  await page.waitForTimeout(300);
  return options.map((s) => s.trim()).filter((s) => s.length > 0);
}

test('Trellis plot: click-to-filter, click-to-select, events, keyboard navigation', async ({page}) => {
  test.setTimeout(900_000);
  page.setDefaultTimeout(120_000);

  // Attached before login so it covers every actuation below. These two are the ONLY error
  // channels available (no `grok.shell.warnings` symbol). The GROK-13205 ladder asserts a
  // DELTA over its own window, never a global zero — a whole-run floor reads noise as red.
  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error' && !isBenignError(m.text())) consoleErrors.push(m.text()); });

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

  // The live DemoFiles demog on dev has 5850 rows (NOT the 1588 the older scenarios
  // quote). Pinned once here so a dataset swap fails loudly at setup; every later
  // expectation derives from the live dataframe rather than restating a literal.
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
  // Categories the Filter Panel keeps in Step 9; Steps 10-12 recompute the panel-only
  // count from it.
  let panelSelected: string[] = [];

  await v.addViewerByIcon(page, 'trellis-plot', 'Trellis-plot', 15000);
  await page.evaluate(async () => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    tp.props.xColumnNames = ['SEX'];
    tp.props.yColumnNames = ['RACE'];
    tp.props.viewerType = 'Scatter plot';
    await new Promise((r) => setTimeout(r, 1800));
  });
  const cellLocator = page.locator('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell');
  await expect(cellLocator).toHaveCount(8);

  // ================= Scenario 1: On Click = Filter + Filter Panel collab + ESC =================

  // #### Scenario 1 Steps 2-4/7: mode switch alone does not alter the inner canvases (GROK-17708)
  await softStep('Scenario 1 Step 7', async () => {
    // GROK-17708 guard: flipping On Click must not re-render cell content, so the diff
    // is taken BEFORE any click.
    const idxA = await cellIndexFor(page, 'F', 'Caucasian');
    const idxB = await cellIndexFor(page, 'M', 'Caucasian');
    const result = await page.evaluate(async ({idxA, idxB}) => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
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
      tp.props.onClick = 'None';
      await new Promise((r) => setTimeout(r, 1500));
      const beforeA = cellHash(idxA), beforeB = cellHash(idxB);
      // The settle wait after the flip is what makes the comparison meaningful: a
      // spurious re-render needs time to land, else "unchanged" only proves we were early.
      tp.props.onClick = 'Filter';
      await new Promise((r) => setTimeout(r, 1500));
      const afterA = cellHash(idxA), afterB = cellHash(idxB);
      return {
        read: beforeA !== null && beforeB !== null && afterA !== null && afterB !== null,
        // Anti-vacuity: two cells hashing DIFFERENTLY proves the probes carry real
        // content — a blank canvas hashes equal and "unchanged" would pass vacuously.
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

  // #### Scenario 1 Step 5: corner-click F|Caucasian drops df.filter to the combo row count
  await softStep('Scenario 1 Step 5', async () => {
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true); df.rows.requestFilter();
      await new Promise((r) => setTimeout(r, 500));
    });
    const expected = await comboRowCount(page, 'SEX', 'F', 'RACE', 'Caucasian');
    const idx = await cellIndexFor(page, 'F', 'Caucasian');
    // Every cell click in this spec is a CORNER click (~6px inset), and must stay one:
    // the cell CENTER is intercepted by the inner scatter canvas and never reaches the
    // trellis handler, so a centered click would silently run no On Click action at all.
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(900);
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

  // #### Scenario 1 Step 6: the two trellis events fire with their payloads
  await softStep('Scenario 1 Step 6', async () => {
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true); df.rows.requestFilter();
      await new Promise((r) => setTimeout(r, 500));
    });
    // Recon fact [DOM 2026-08-06]: the CORNER click fires only current-cell-changed;
    // d4-trellis-plot-inner-viewer-clicked fires on MOUSE-DOWN over the inner CANVAS (center).
    // Neither gesture raises both events, which is why one cell gets two gestures here.
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
    await page.waitForTimeout(400);
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(900);
    const ev = await page.evaluate(() => {
      ((window as any).__subs || []).forEach((s: any) => s?.unsubscribe?.());
      return (window as any).__ev;
    });
    // An event that fires with an EMPTY payload is not a pass.
    expect(Array.isArray(ev.inner)).toBe(true);
    expect(ev.inner.length).toBeGreaterThan(0);
    expect(ev.cell).toMatchObject({SEX: 'F', RACE: 'Caucasian'});
  });

  // #### Scenario 1 Steps 8-9: a persistent panel filter ANDs with the trellis cell click
  await softStep('Scenario 1 Step 9', async () => {
    // The trellis contribution left by Steps 5-6 has to go the PRODUCT way before the
    // panel-only value is measured: `setAll(true)` + `requestFilter()` does NOT remove it
    // [DOM 2026-08-12] — the trellis re-applies its cell and panel-only reads 2771, not 5573.
    const seedIdx = await cellIndexFor(page, 'F', 'Caucasian');
    await cellLocator.nth(seedIdx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(500);
    await page.keyboard.press('Escape');
    await page.waitForTimeout(900);
    const cleared = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const filters: string[] = [];
      for (const f of df.rows.filters) filters.push(f);
      return {trueCount: df.filter.trueCount,
        trellisEntries: filters.filter((s) => s.startsWith('SEX:') || s.startsWith('RACE:')).length};
    });
    // Load-bearing witness: if the cell survived, every number below is silently an AND.
    expect(cleared).toEqual({trueCount: fullRowCount, trellisEntries: 0});
    // The REAL Filter Panel: a private df.onRowsFiltering subscription would prove only that
    // the bitset arithmetic ANDs, leaving the panel channel — the actual claim of
    // int.onclick-filter-panel-collab — untested.
    await v.openFilterPanel(page);
    const cats: string[] = await page.evaluate((col) => grok.shell.tv.dataFrame.col(col).categories, PANEL_COLUMN);
    expect(cats.length).toBeGreaterThan(1);
    // Drop exactly one category: the panel-only count stays far above any single cell's
    // row count, so the AND below lands strictly under BOTH contributors.
    panelSelected = cats.filter((c) => c !== cats[0]);
    const expectedPanelOnly = await panelOnlyRowCount(page, PANEL_COLUMN, panelSelected);
    const {filteredCount: panelOnly} = await v.applyCategoricalFilter(page, PANEL_COLUMN, panelSelected);
    expect(panelOnly).toBe(expectedPanelOnly);
    expect(panelOnly).toBeLessThan(fullRowCount);
    expect(panelOnly).toBeGreaterThan(0);
    // The cell index math below is only valid while the grid still spans the full 2x4
    // product — true under rowSource = All, not if packing repacked it.
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
    await page.waitForTimeout(900);
    const after = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    // The equality alone could hold with one of the two filters silently dropped, so
    // both inequalities pin the result strictly BELOW each filter taken on its own.
    expect(after).toBe(geom.mAsianPanel);
    expect(after).toBeLessThan(panelOnly);
    expect(after).toBeLessThan(geom.mAsian);
  });

  // #### Scenario 1 Step 10: ESC removes only the trellis contribution; panel survives
  await softStep('Scenario 1 Step 10', async () => {
    const panelOnly = await panelOnlyRowCount(page, PANEL_COLUMN, panelSelected);
    // ESC only reaches the charts grid when keyboard focus is already there, and only a
    // REAL click puts it there (a synthetic focus() does not). Dropping the click before
    // the key press turns this step into a no-op.
    const idx = await cellIndexFor(page, 'M', 'Asian');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(500);
    // Witness that there IS a contribution to remove: without it a click that silently
    // did nothing leaves the count at panel-only and ESC "passes".
    const beforeEsc = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(beforeEsc).toBeLessThan(panelOnly);
    await page.keyboard.press('Escape');
    await page.waitForTimeout(900);
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

  // #### Scenario 1 Step 11: changing the X split column drops the trellis contribution
  await softStep('Scenario 1 Step 11', async () => {
    const panelOnly = await panelOnlyRowCount(page, PANEL_COLUMN, panelSelected);
    const idx = await cellIndexFor(page, 'M', 'Asian');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(600);
    // Same witness as the ESC step: a DROP can only be shown against a contribution that
    // was demonstrably there a moment earlier.
    const beforeChange = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(beforeChange).toBeLessThan(panelOnly);
    const after = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['CONTROL'];
      await new Promise((r) => setTimeout(r, 1800));
      return grok.shell.tv.dataFrame.filter.trueCount;
    });
    expect(after).toBe(panelOnly);
    // Restore X = SEX: the later scenarios assume the SEX x RACE split.
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['SEX'];
      await new Promise((r) => setTimeout(r, 1500));
    });
  });

  // #### Scenario 1 Step 12: resetting all filters restores the full row count
  await softStep('Scenario 1 Step 12', async () => {
    // A reset asserted from an already-full frame proves nothing, hence this witness.
    const before = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(before).toBeLessThan(fullRowCount);
    // Reset the PRODUCT way — NOT df.filter.setAll(true), NOT v.resetFilters (which ends in
    // setAll): asserting that an all-true bitset is all-true is a tautology. This leaves the
    // filters group EMPTIED, which constrains every later panel use — see Scenario 3 Step 6.
    const after = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const fg = tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) fg.remove(f);
      await new Promise((r) => setTimeout(r, 500));
      tv.dataFrame.rows.requestFilter();
      await new Promise((r) => setTimeout(r, 900));
      const filters: string[] = [];
      for (const f of tv.dataFrame.rows.filters) filters.push(f);
      return {filterCount: tv.dataFrame.filter.trueCount, filters};
    });
    expect(after.filterCount).toBe(fullRowCount);
    expect(after.filters.filter((s) => s.startsWith(PANEL_COLUMN + ':'))).toEqual([]);
  });

  // ================= Scenario 2: On Click = Filter set via Context Panel (GROK-17711) =================

  // #### Scenario 2 Step 3/4: the Context Panel actuation commits the filter subscription
  await softStep('Scenario 2 Step 3', async () => {
    // GROK-17711 (actuation-path parity): On Click = Filter set through the Context Panel must
    // engage the filter subscription just like the property channel does — the bug was a SILENT
    // no-op that left the property looking correct while df.filter stayed at the full count.
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true); df.selection.setAll(false); df.rows.requestFilter();
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.onClick = 'None';
      await new Promise((r) => setTimeout(r, 500));
    });
    // Park Row Source on `Filtered` FIRST, or the auto-correct graded below has nothing to
    // correct: Scenario 1 already knocked it off. Through the PANEL, the only actuation channel
    // this scenario may use; the literal is safe here — `Filtered` is alike on both surfaces.
    const parked = await setTrellisChoiceViaPanel(page, 'row-source', 'Filtered');
    expect(parked).toBe('Filtered');
    const beforeCommit = await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      return tp.props.rowSource;
    });
    expect(beforeCommit).toBe('Filtered');
    // THE OTHER ACTUATION PATH: the real Context Panel property-grid row, never
    // `tp.props.onClick`. Assigning the prop here would make this scenario a duplicate of
    // Scenario 1 on a different cell and leave GROK-17711 with no coverage at all.
    const shown = await setTrellisChoiceViaPanel(page, 'on-click', 'Filter');
    // Actuation guard, not the signal: the panel echoing 'Filter' proves the editor
    // committed. The corner-click below is what proves the subscription engaged.
    expect(shown).toBe('Filter');
    const committed = await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      return {onClick: tp.props.onClick, rowSource: tp.props.rowSource};
    });
    // rowSource auto-correcting away from Filtered (trellis_plot_core.dart:582) is the SIDE
    // EFFECT of the real commit path. Graded as the exact landing on `All` [refdoc
    // trellis_plot.md Interaction table, DOM 2026-08-05], read off props.rowSource (raw names).
    expect(committed.onClick).toBe('Filter');
    expect(committed.rowSource).toBe('All');

    const expected = await comboRowCount(page, 'SEX', 'F', 'RACE', 'Black');
    const idx = await cellIndexFor(page, 'F', 'Black');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(900);
    const after = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const filters: string[] = [];
      for (const f of df.rows.filters) filters.push(f);
      return {filterCount: df.filter.trueCount, filters};
    });
    expect(after.filterCount).toBe(expected);
    expect(after.filters).toContain('SEX: F');
    expect(after.filters).toContain('RACE: Black');
    // Resets the FILTER, not the CURRENT CELL — the cell clicked above outlives both ESC
    // and On Click = None [DOM 2026-08-12]; the ladder below clears the cell itself
    // rather than trusting this reset.
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.onClick = 'None';
      df.filter.setAll(true); df.rows.requestFilter();
      await new Promise((r) => setTimeout(r, 600));
    });
  });

  // #### Scenario 2 Step 5: the Row Source ladder — every row source keeps the trellis
  // filter contract and raises no error (GROK-13205)
  await softStep('Scenario 2 Step 5', async () => {
    const errBefore = consoleErrors.length;
    const pageErrBefore = pageErrors.length;
    // Packing follows the ROW SOURCE's combined filter (trellis_plot_core.dart:163/182), so a
    // narrow source legitimately collapses a PACKED grid; pack off pins the full 2x4 product
    // for every rung. The broad selection keeps the selection-based sources from being empty.
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      df.filter.setAll(true); df.rows.requestFilter();
      tp.props.packCategories = false;
      df.selection.init((i: number) => i % 2 === 0);
      await new Promise((r) => setTimeout(r, 1500));
    });
    await setTrellisChoiceViaPanel(page, 'on-click', 'Filter');
    // Clear the trellis's CURRENT CELL: the F|Black click seated one, it survives both ESC and
    // On Click = None, and returning to Filter re-applies its constraint with no click
    // [DOM 2026-08-12]. A SPLIT COLUMN change is the one product path observed to drop it.
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['CONTROL'];
      await new Promise((r) => setTimeout(r, 1500));
      tp.props.xColumnNames = ['SEX'];
      await new Promise((r) => setTimeout(r, 1800));
      // Force every contributor to re-announce itself: a cell that secretly survived the
      // split change would re-apply here and be caught by the entry assert below.
      grok.shell.tv.dataFrame.rows.requestFilter();
      await new Promise((r) => setTimeout(r, 800));
    });
    const start = await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      const df = grok.shell.tv.dataFrame;
      const filters: string[] = [];
      for (const f of df.rows.filters) filters.push(f);
      return {onClick: tp.props.onClick, trueCount: df.filter.trueCount,
        trellisEntries: filters.filter((s) => s.startsWith('SEX:') || s.startsWith('RACE:')).length};
    });
    // Entry state the whole ladder is graded against: On Click IS Filter (so the Filtered
    // rung below has an auto-correct to make) and nothing is filtered yet — asserted on
    // both channels, the same two reads every rung repeats.
    expect(start.onClick).toBe('Filter');
    expect(start.trueCount).toBe(fullRowCount);
    expect(start.trellisEntries).toBe(0);

    const choices = await propertyGridChoices(page, 'row-source');
    // Graded as an exact SET, never a floor of three: a build that dropped or renamed one of
    // the eight RowSet values (ddt utils.dart:466-475, offered under their RAW names — see
    // ROW_SOURCE_OPTIONS) would otherwise walk a shorter ladder and still pass.
    expect([...choices].sort()).toEqual([...ROW_SOURCE_OPTIONS].sort());

    const rungs: {source: string; shown: string; rowSource: string; onClick: string;
      cells: number; trueCount: number; trellisEntries: number; mapped: boolean; wide: boolean;
      fedCells: number; blankFedCells: number[]; distinctHashes: number; signature: string}[] = [];
    for (const source of choices) {
      const shown = await setTrellisChoiceViaPanel(page, 'row-source', source);
      const state = await page.evaluate(async () => {
        // A RENDER settle, not just a property settle: the cell canvases are read below,
        // and a re-render that has not landed yet would be graded as "blank".
        await new Promise((r) => setTimeout(r, 1500));
        const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
        const df = grok.shell.tv.dataFrame;
        const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
        const filters: string[] = [];
        for (const f of df.rows.filters) filters.push(f);

        // `painted` is this section's emptiness definition at pixel level: a cell is empty
        // when it has no canvas child (the bare DIV of Scenario 3 Step 7) or its canvas
        // carries a single flat colour. Same getImageData probe as the GROK-17708 guard.
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

        // Which rows this row source feeds, restated over the LIVE frame: the js-api exposes
        // no `getRowsBitSet`, so each RowSet value (ddt utils.dart:466-475) is expressed
        // through the frame state its NAME denotes. Width is not declared — it falls out below.
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
          // The ladder never puts a pointer over the grid, so the mouse-over group is
          // empty by construction and can carry no render signal here.
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
          // WIDE = the row source feeds MORE THAN ONE category combination. Only then is
          // "every cell draws the same picture" / "a fed cell is blank" a defect: a source
          // that feeds one row or none legitimately renders one cell, or none at all.
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
      // The editor committed this rung; without it the row grades the previous source.
      expect(r.shown).toBe(r.source);
      // The platform never lets (On Click = Filter, Row Source = Filtered) stand:
      // picking either one rewrites the other (trellis_plot_core.dart:589, 610-613).
      expect(r.onClick === 'Filter' && r.rowSource === 'Filtered').toBe(false);
      // Pack off: every category combination keeps its cell whatever the row source is.
      expect(r.cells).toBe(8);
      // No cell is clicked anywhere in this ladder, so no row source may leave a trellis
      // contribution behind — a row source that quietly filters on its own is the
      // GROK-13205 symptom class.
      expect(r.trellisEntries).toBe(0);
      expect(r.trueCount).toBe(fullRowCount);
    }
    // The auto-correct asserted positively, not as an absence: the rung landing on
    // Filtered must be the one that pushed On Click off Filter.
    const filteredIdx = rungs.findIndex((r) => r.rowSource === 'Filtered');
    expect(filteredIdx).toBeGreaterThanOrEqual(0);
    // ...witnessed LIVE on the rung IMMEDIATELY BEFORE it: the entry assert is several rungs
    // back, so `None` could equally be inherited. Only `Filtered` rewrites On Click
    // (trellis_plot_core.dart:589, 609-613); if Filtered is FIRST, the entry state is witness.
    const beforeFiltered = filteredIdx > 0 ? rungs[filteredIdx - 1].onClick : start.onClick;
    expect(beforeFiltered).toBe('Filter');
    expect(rungs[filteredIdx].onClick).toBe('None');

    // ---- GROK-13205's own symptom: IDENTICAL or EMPTY cells ---- The floor above cannot see
    // it: a ladder can raise no error, keep all 8 cells and leave df.filter untouched while
    // every cell draws one repeated picture. An inexpressible row source fails loudly here.
    expect(rungs.filter((r) => !r.mapped).map((r) => r.rowSource)).toEqual([]);
    // The classifier itself has to work: with the frame unfiltered and half the rows
    // selected, All cannot be anything but wide. Without this, a classifier that called
    // everything narrow would grade nothing.
    expect(rungs.some((r) => r.rowSource === 'All' && r.wide)).toBe(true);
    const wideRungs = rungs.filter((r) => r.wide);
    // (a) "identical data": the cells fed by this row source do not all carry one hash.
    expect(wideRungs.filter((r) => r.distinctHashes < 2)
      .map((r) => `${r.rowSource}: ${r.distinctHashes} hash(es) over ${r.fedCells} fed cells`)).toEqual([]);
    // (b) "empty data": no cell that still receives rows is blank. Cells the row source
    // feeds NO rows are not graded — a blank cell is the correct rendering there.
    expect(wideRungs.filter((r) => r.blankFedCells.length > 0)
      .map((r) => `${r.rowSource}: ${r.blankFedCells.join(',')}`)).toEqual([]);
    // Driven-guard for (a) and (b) together: the per-cell hash signature must differ
    // between at least two rungs. If the row source never reached the render, every rung
    // would report the SAME frame and "the cells differ" would grade an untouched grid.
    expect(Array.from(new Set(rungs.map((r) => r.signature))).length).toBeGreaterThanOrEqual(2);

    // rowSource MUST end on All: Scenario 3 filters SEX = M and still needs the F column
    // present, which only holds while the trellis reads all rows — the shipped default
    // `Filtered` repacks the grid. Restored through the option text the ladder OBSERVED.
    await setTrellisChoiceViaPanel(page, 'row-source', rungs.find((r) => r.rowSource === 'All')?.source ?? 'All');
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.onClick = 'None';
      tp.props.packCategories = true;
      df.selection.setAll(false);
      df.filter.setAll(true); df.rows.requestFilter();
      await new Promise((r) => setTimeout(r, 1500));
    });
    const restored = await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      return {rowSource: tp.props.rowSource, onClick: tp.props.onClick};
    });
    expect(restored).toEqual({rowSource: 'All', onClick: 'None'});
    await expect(cellLocator).toHaveCount(8);
  });

  // ================= Scenario 3: On Click = Select — additive Ctrl+click + GROK-19790 =================

  // #### Scenario 3 Steps 3-4: plain select equals the cell structural row count (GROK-17714)
  await softStep('Scenario 3 Step 4', async () => {
    const viewportBefore = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      const before = {xCat: tp.xCategoriesCount, yCat: tp.yCategoriesCount,
        cells: document.querySelectorAll('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell').length};
      tp.props.onClick = 'Select';
      await new Promise((r) => setTimeout(r, 700));
      return before;
    });
    const expected = await comboRowCount(page, 'SEX', 'F', 'RACE', 'Caucasian');
    const idx = await cellIndexFor(page, 'F', 'Caucasian');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(900);
    const after = await page.evaluate(() => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      return {sel: grok.shell.tv.dataFrame.selection.trueCount,
        xCat: tp.xCategoriesCount, yCat: tp.yCategoriesCount,
        cells: document.querySelectorAll('[name="viewer-Trellis-plot"] .d4-trellis-plot-cell').length};
    });
    expect(after.sel).toBe(expected);
    // GROK-17714: the +/- viewport configuration must survive the On Click = Select switch
    // and the first click — the bug rebuilt the grid out from under the user.
    expect({xCat: after.xCat, yCat: after.yCat, cells: after.cells})
      .toEqual({xCat: viewportBefore.xCat, yCat: viewportBefore.yCat, cells: viewportBefore.cells});
  });

  // #### Scenario 3 Step 6: Ctrl+click a filter-emptied cell ADDS its structural rows
  await softStep('Scenario 3 Step 6', async () => {
    // A filter-emptied cell: SEX = M hides F|Caucasian's rows while the structural combination
    // still exists — not the structurally empty cell of the next step. Selection is
    // dataset-level, so the Ctrl+click must still ADD the cell's full structural row count.
    const fCaucasian = await comboRowCount(page, 'SEX', 'F', 'RACE', 'Caucasian');
    const before = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      // Seed a DIFFERENT cell than the one clicked below, so additive growth is
      // distinguishable from a replace.
      df.selection.setAll(false);
      const s = df.col('SEX'), r = df.col('RACE');
      for (let i = 0; i < df.rowCount; i++)
        if (s.get(i) === 'M' && r.get(i) === 'Black') df.selection.set(i, true, false);
      await new Promise((r) => setTimeout(r, 300));
      return df.selection.trueCount;
    });
    expect(before).toBeGreaterThan(0);
    // ORDER MATTERS: Scenario 1 Step 12 emptied the filters group and `getFiltersGroup()` does
    // not refill it, so `openFilterPanel` FIRST would wait out its whole 15s `.d4-filter`
    // timeout [DOM 2026-08-12]. Adding the filter rebuilds one, and that wait is the proof.
    const {filteredCount} = await v.applyCategoricalFilter(page, 'SEX', ['M']);
    await v.openFilterPanel(page);
    expect(filteredCount).toBeLessThan(fullRowCount);
    expect(filteredCount).toBeGreaterThan(0);
    // The F column has to SURVIVE the SEX = M filter for this step to be about a
    // filter-emptied cell at all — it does because the trellis reads all rows (rowSource =
    // All); packing over filtered data would drop it and the index would address another cell.
    await expect(cellLocator).toHaveCount(8);
    const idx = await cellIndexFor(page, 'F', 'Caucasian');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}, modifiers: ['Control']});
    await page.waitForTimeout(900);
    const after = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(after).toBe(before + fCaucasian);
    // Clear through the panel itself: the next step re-splits the grid from scratch and
    // must not inherit a SEX constraint.
    const cleared = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const fg = tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) fg.remove(f);
      await new Promise((r) => setTimeout(r, 500));
      tv.dataFrame.rows.requestFilter();
      await new Promise((r) => setTimeout(r, 700));
      return tv.dataFrame.filter.trueCount;
    });
    expect(cleared).toBe(fullRowCount);
  });

  // #### Scenario 3 Step 7: Ctrl+click a STRUCTURALLY empty cell must NOT clear (GROK-19790)
  await softStep('Scenario 3 Step 7', async () => {
    // Pack Categories must be OFF for a structurally empty cell to exist at all;
    // Asian|Critical is such a cell (a DIV with no canvas). GROK-19790 guard: Ctrl+click on
    // it must add nothing rather than wipe the selection built on a populated cell.
    const before = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      df.selection.setAll(false);
      tp.props.packCategories = false;
      tp.props.xColumnNames = ['RACE'];
      tp.props.yColumnNames = ['SEVERITY'];
      await new Promise((r) => setTimeout(r, 2000));
      // A canvas child is the populated-cell marker: empty cells render a bare DIV.
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const cells = root.querySelectorAll('.d4-trellis-plot-cell');
      let popIdx = -1;
      for (let i = 0; i < cells.length; i++) if (cells[i].querySelector('canvas')) { popIdx = i; break; }
      return {popIdx};
    });
    await cellLocator.nth(before.popIdx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(700);
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
    // A DIFFERENT cell from the one that seeded the selection — otherwise the highlight
    // already sits on it and the arrival witness below is true before the click is sent.
    expect(emptyIdx).not.toBe(before.popIdx);
    await cellLocator.nth(emptyIdx).click({position: {x: 6, y: 6}, modifiers: ['Control']});
    await page.waitForTimeout(900);
    const after = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const cells = Array.from(root.querySelectorAll('.d4-trellis-plot-cell'));
      const marked = cells.filter((c) => c.classList.contains('d4-trellis-cell-current'));
      return {sel: grok.shell.tv.dataFrame.selection.trueCount,
        currentCount: marked.length, currentIdx: cells.indexOf(marked[0])};
    });
    expect(after.sel).toBe(selBefore);
    // DRIVEN-GUARD: "the selection did not change" is ALSO what a click that never reached the
    // trellis produces, so the relocated highlight witnesses arrival. The CLASS, not
    // `matchCondition` — that payload is built from the very field GROK-19790 degenerates.
    expect({count: after.currentCount, idx: after.currentIdx}).toEqual({count: 1, idx: emptyIdx});
  });

  // #### Scenario 3 Step 8: a plain click on a different populated cell REPLACES the selection
  await softStep('Scenario 3 Step 8', async () => {
    // The target must NOT be the cell that SEEDED the standing selection, and must hold a
    // DIFFERENT number of rows: on a re-click "replaces" reads identically to "adds" and to
    // "the click did nothing". Which is why Step 7's first-populated-cell scan cannot pick it.
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
    // Entry conditions that make the grading possible at all: a grid with no qualifying
    // candidate fails LOUDLY here instead of degenerating into a re-click of the seeded
    // cell, the shape that made this step unfalsifiable before.
    expect(target.before).toBeGreaterThan(0);
    expect(target.idx).toBeGreaterThanOrEqual(0);
    expect(target.idx).not.toBe(target.seedIdx);
    expect(target.expected).not.toBe(target.before);
    await cellLocator.nth(target.idx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(900);
    const after = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(after).toBe(target.expected);
    // Neither what it was before (a click that did nothing) nor the sum of the two
    // disjoint cells (a click that ADDED instead of replacing).
    expect(after).not.toBe(target.before);
    expect(after).not.toBe(target.before + target.expected);
    // Restore the SEX x RACE split + Pack on: Scenario 4 navigates the 8-cell grid.
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      df.selection.setAll(false);
      tp.props.packCategories = true;
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      tp.props.onClick = 'None';
      await new Promise((r) => setTimeout(r, 1800));
    });
    await expect(cellLocator).toHaveCount(8);
  });

  // ================= Scenario 4: Arrow-key navigation — current-cell-changed events =================

  // #### Scenario 4 Step 3: each arrow press fires current-cell-changed with a new matchCondition
  await softStep('Scenario 4 Step 3', async () => {
    // The click moves keyboard focus onto the charts grid (see the ESC step); without it
    // the arrow presses go nowhere and no event is ever recorded.
    const idx = await cellIndexFor(page, 'F', 'Asian');
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(600);
    // WHERE the highlight sits before the arrows, not merely THAT one exists: the click above
    // already seated a current cell, so a presence read stays true even for the regression it
    // names. Selector pinned in refdoc trellis_plot.md.
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
    await page.keyboard.press('ArrowRight');
    await page.waitForTimeout(400);
    await page.keyboard.press('ArrowDown');
    await page.waitForTimeout(400);
    await page.keyboard.press('ArrowLeft');
    await page.waitForTimeout(400);
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
    // A set of THREE, not a first-pair difference: the presses land on (x+1, y), (x+1, y'),
    // (x, y') — pairwise different by geometry — so cc[0] vs cc[1] would admit a build where
    // ONE arrow stopped moving the cell while still announcing itself. Key order cannot split.
    const combos = Array.from(new Set((events.cc as any[])
      .map((mc) => JSON.stringify([mc.SEX, mc.RACE]))));
    expect(combos.length).toBeGreaterThanOrEqual(3);
    // The highlight MOVED, graded as a relocation of a SINGLE class: two highlighted cells is
    // the "the previous cell kept it" regression, and Right -> Down -> Left ends one row below
    // F|Asian by geometry, so an index differing from the seeded one is the expectation.
    expect(events.currentCount).toBe(1);
    expect(events.currentIdx).not.toBe(currentBefore.idx);
  });

  // #### Scenario 4 Step 4: with Filter mode, an arrow press re-applies the filter to the new cell
  await softStep('Scenario 4 Step 4', async () => {
    // Under On Click = Filter an arrow press both moves the current cell and re-applies
    // the trellis filter to it [DOM 2026-08-06] — so the count after the press must
    // EQUAL the newly focused cell's row count, which is what this step grades.
    const idx = await cellIndexFor(page, 'F', 'Asian');
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true); df.rows.requestFilter();
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.onClick = 'None';
      await new Promise((r) => setTimeout(r, 400));
    });
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(500);
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.onClick = 'Filter';
      await new Promise((r) => setTimeout(r, 500));
    });
    // The first click ran under On Click = None, so it seated no filter — a second
    // click after the mode switch is what gives the arrow press something to move.
    await cellLocator.nth(idx).click({position: {x: 6, y: 6}});
    await page.waitForTimeout(600);
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
    // The click must have seated a filter, else "the filter follows the arrow" has nothing
    // to follow and every reading below is degenerate.
    expect(beforeArrow).toBeLessThan(fullRowCount);
    await page.keyboard.press('ArrowDown');
    await page.waitForTimeout(800);
    const moved = await page.evaluate(() => {
      (window as any).__ccMoveSub?.unsubscribe?.();
      const cc = (window as any).__ccMove;
      return {mc: cc.length > 0 ? cc[cc.length - 1] : null, filterCount: grok.shell.tv.dataFrame.filter.trueCount};
    });
    // The filter must equal THAT cell's row count. A "changed and non-zero" hedge would
    // also be satisfied by the regression itself: an arrow that RESETS the filter instead
    // of carrying it lands on the full row count, which is both different and non-zero.
    expect(typeof moved.mc?.SEX).toBe('string');
    expect(typeof moved.mc?.RACE).toBe('string');
    // A DIFFERENT combination than the clicked one — that is what "the current cell moved"
    // means. Compared as combinations, not as counts: two neighbouring cells can
    // legitimately hold the same number of rows.
    expect({SEX: moved.mc.SEX, RACE: moved.mc.RACE}).not.toEqual({SEX: 'F', RACE: 'Asian'});
    const expectedAfterArrow = await comboRowCount(page, 'SEX', moved.mc.SEX, 'RACE', moved.mc.RACE);
    expect(expectedAfterArrow).toBeGreaterThan(0);
    expect(moved.filterCount).toBe(expectedAfterArrow);
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
      tp.props.onClick = 'None';
      df.filter.setAll(true); df.rows.requestFilter();
      await new Promise((r) => setTimeout(r, 500));
    });
  });

  v.finishSpec();
});
