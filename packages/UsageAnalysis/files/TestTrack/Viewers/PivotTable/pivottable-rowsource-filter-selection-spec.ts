/* ---
realizes: [pivottable.cp.rowsource-filter-selection-links, pivottable.int.filtering-requires-row-source-all, pivottable.int.selection-link-gated-by-row-source]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const PIVOT = '[name="viewer-Pivot-table"]';
const DEMOG_ROWS = 5850;

// Read the pivot's configuration props.
async function pivotProps(page: Page) {
  return page.evaluate(() => {
    const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
    return {
      groupBy: pv.props.groupByColumnNames,
      pivot: pv.props.pivotColumnNames,
      agg: pv.props.aggregateColumnNames,
      aggTypes: pv.props.aggregateAggTypes,
      rowSource: pv.props.rowSource,
      filteringEnabled: pv.props.filteringEnabled,
    };
  });
}

// Configure a minimal single-column cross-tab (Group by DIS_POP, avg(AGE), no Pivot) through the
// look props, then wait for the rebuild.
async function setupSingleKey(page: Page, rowSource: string, filteringEnabled: boolean) {
  await page.evaluate(async (opts: {rs: string; fe: boolean}) => {
    const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
    pv.props.groupByColumnNames = ['DIS_POP'];
    pv.props.pivotColumnNames = [];
    pv.props.aggregateColumnNames = ['AGE'];
    pv.props.aggregateAggTypes = ['avg'];
    pv.props.rowSource = opts.rs;
    pv.props.filteringEnabled = opts.fe;
    await new Promise((r) => setTimeout(r, 700));
  }, {rs: rowSource, fe: filteringEnabled});
}

// The inner pivot grid canvas bounding box (the aggregated result grid, a nested d4 Grid).
async function innerGridRect(page: Page) {
  const box = await page.locator(`${PIVOT} [name="canvas"]`).first().boundingBox();
  if (!box) throw new Error('pivot inner grid canvas not visible');
  return box;
}

// DIS_POP key cell coordinates (single-key demog pivot):
// x/y are derived from the 24px header and row bands.
// Uses trusted page.mouse.click because synthetic events do not trigger the cell handler.
async function clickAggregateCell(page: Page, rowIdx: number) {
  const r = await innerGridRect(page);
  await page.mouse.click(r.x + 55, r.y + 24 + rowIdx * 24 + 12);
  await page.waitForTimeout(500);
}

// Ctrl+click on the DIS_POP key cell triggers the selection link to source df.
// Plain click does not affect df.selection. The gesture is isolated from cell
// filtering (Shift+click under All is intentionally not used).
async function ctrlClickAggregateCell(page: Page, rowIdx: number) {
  const r = await innerGridRect(page);
  await page.keyboard.down('Control');
  await page.mouse.click(r.x + 55, r.y + 24 + rowIdx * 24 + 12);
  await page.keyboard.up('Control');
  await page.waitForTimeout(500);
}

// The full source-row count of a DIS_POP group (the SELECTION_TO_SELECTION link selects exactly the
// clicked aggregate row's source rows, so this is the reference df.selection.trueCount must match).
async function groupSourceCount(page: Page, disPopKey: string) {
  return page.evaluate((k: string) => {
    const df = grok.shell.tv.dataFrame;
    const col = df.col('DIS_POP');
    let n = 0;
    for (let i = 0; i < df.rowCount; i++) if (col.get(i) === k) n++;
    return n;
  }, disPopKey);
}

// df-level signal reads on the source dataframe.
const srcFilterCount = (page: Page) => page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
const srcSelectionCount = (page: Page) => page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);

// Durable pivot-config caption read-back: the .d4-tag caption text inside the row whose title
// carries d4-name=<title>. The chips' `name` attributes are dropped after any prop-driven rebuild;
// the caption text survives.
async function pivotRowCaptions(page: Page, title: string) {
  return page.evaluate((t: string) => {
    const scope = document.querySelector('[name="viewer-Pivot-table"]');
    const row = scope && [...scope.querySelectorAll('.grok-pivot-column-panel')]
      .find((p) => p.querySelector(`.grok-pivot-column-tags-title[d4-name="${t}"]`));
    return row ? [...row.querySelectorAll('.d4-tag')].map((c) => (c.textContent || '').trim()) : [];
  }, title);
}

test('Pivot Table — Row Source, filter and selection links', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  // --- Setup: open demog, add the Pivot Table viewer, minimal single-key cross-tab ---
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    grok.shell.addTableView(df);
    await new Promise((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
  await v.addViewerByIcon(page, 'pivot-table', 'Pivot-table', 15000);
  await page.waitForTimeout(800);

  await softStep('Setup: tag-editor header shows Group by, Aggregate and Pivot rows', async () => {
    const titles = await page.evaluate(() => Array.from(
      document.querySelectorAll('[name="viewer-Pivot-table"] .grok-pivot-column-tags-title'))
      .map((t) => t.getAttribute('d4-name')));
    expect(titles).toEqual(expect.arrayContaining(['Group by', 'Aggregate', 'Pivot']));
    // Minimal single-key cross-tab: DIS_POP / avg(AGE), no pivot column.
    await setupSingleKey(page, 'Filtered', true);
  });

  // === Scenario 1: Filter consumption — pivot recomputes when the source filter changes (channel 3) ===

  await softStep('Scenario 1 Step 4: after a SEX = "M" panel filter, the pivot aggregate matches an independent groupBy on the filtered df', async () => {
    // Baseline: the unfiltered per-DIS_POP avg(AGE), computed independently.
    const baseline = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const g = df.groupBy(['DIS_POP']).avg('AGE').aggregate();
      const avgCol = Array.from({length: g.columns.length}, (_: any, i: number) => g.columns.byIndex(i).name)
        .find((c: string) => c.toLowerCase().includes('age'));
      const out: Record<string, number> = {};
      const key = g.col('DIS_POP');
      for (let i = 0; i < g.rowCount; i++) out[key.get(i)] = g.col(avgCol!).get(i);
      return out;
    });

    // Open the filter panel, then set SEX = 'M' on df.filter directly — the panel widget is not
    // the tested channel.
    await page.evaluate(() => grok.shell.tv.getFiltersGroup());
    await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 15000});
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const sex = df.col('SEX');
      df.filter.init((i: number) => sex.get(i) === 'M');
    });
    await page.waitForTimeout(900);

    const filtered = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      // The pivot aggregates over df.filter, so an independent groupBy over the filtered rows is the
      // reference the pivot's recomputed values must match.
      const g = df.groupBy(['DIS_POP']).avg('AGE').whereRowMask(df.filter).aggregate();
      const avgCol = Array.from({length: g.columns.length}, (_: any, i: number) => g.columns.byIndex(i).name)
        .find((c: string) => c.toLowerCase().includes('age'));
      const out: Record<string, number> = {};
      const key = g.col('DIS_POP');
      for (let i = 0; i < g.rowCount; i++) out[key.get(i)] = g.col(avgCol!).get(i);
      return out;
    });

    // The pivot recomputed: the filtered per-DIS_POP avg differs from the unfiltered baseline for at
    // least one group (a SEX split changes the mean), and every filtered group value matches the
    // independent filtered aggregation.
    const changed = Object.keys(filtered).some((k) => baseline[k] !== undefined
      && Math.abs(filtered[k] - baseline[k]) > 1e-9);
    expect(changed).toBe(true);
    // Publish the pivot's own aggregated table via ADD and confirm it agrees with the independent
    // filtered aggregation (the pivot recomputed over df.filter, not the full source).
    await page.locator(`${PIVOT} .grok-pivot-counts [name="button-ADD"]`).click();
    await page.waitForTimeout(1000);
    const pub = await page.evaluate(() => {
      const aggView = Array.from(grok.shell.views).find((vw: any) => vw.name === 'Table aggregation') as any;
      const t = aggView?.table;
      if (!t) return null;
      const avgCol = Array.from({length: t.columns.length}, (_: any, i: number) => t.columns.byIndex(i).name)
        .find((c: string) => c.toLowerCase().includes('age'));
      const out: Record<string, number> = {};
      const key = t.col('DIS_POP');
      for (let i = 0; i < t.rowCount; i++) out[key.get(i)] = t.col(avgCol!).get(i);
      return out;
    });
    expect(pub).not.toBeNull();
    for (const k of Object.keys(pub!)) expect(Math.abs(pub![k] - filtered[k])).toBeLessThan(1e-6);
  });

  await softStep('Scenario 1 Step 6: clear the SEX filter and close the published table before the next scenario', async () => {
    await page.evaluate(async () => {
      const aggView = Array.from(grok.shell.views).find((vw: any) => vw.name === 'Table aggregation') as any;
      if (aggView) aggView.close();
      grok.shell.tv.dataFrame.filter.setAll(true);   // clear the source filter (no in-viewer undo — GROK-14996)
      await new Promise((r) => setTimeout(r, 400));
    });
    expect(await srcFilterCount(page)).toBe(DEMOG_ROWS);
  });

  // === Scenario 2: Cell-driven filtering under rowSource=All + filteringEnabled=true (channel 1) ===

  await softStep('Scenario 2 Step 3: Row Source = All with filteringEnabled = true does not widen the aggregated row set', async () => {
    // Pre-filter the source to SEX = 'M' (panel filter still active from the scenario intent).
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const sex = df.col('SEX');
      df.filter.init((i: number) => sex.get(i) === 'M');
      await new Promise((r) => setTimeout(r, 400));
    });
    const rowsBefore = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return df.groupBy(['DIS_POP']).avg('AGE').whereRowMask(df.filter).aggregate().rowCount;
    });
    // Switch Row Source to All via the look prop (setting; the property panel exposes no native
    // input for it). filteringEnabled stays true (the default).
    await setupSingleKey(page, 'All', true);
    const props = await pivotProps(page);
    expect(props.rowSource).toBe('All');
    expect(props.filteringEnabled).toBe(true);
    // Row Source = All does NOT widen the aggregated row set: the aggregation still runs over
    // df.filter (the SEX = 'M' subset), so the per-DIS_POP aggregate row count is unchanged.
    const rowsAfter = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return df.groupBy(['DIS_POP']).avg('AGE').whereRowMask(df.filter).aggregate().rowCount;
    });
    expect(rowsAfter).toBe(rowsBefore);
  });

  await softStep('Scenario 2 Step 5: a real cell click narrows df.filter.trueCount to a smaller positive value', async () => {
    const before = await srcFilterCount(page);
    // Trusted click on the first aggregate cell fires the Dart cell-driven filter. df.rows.filters
    // iterates empty (the filter label is not JS-readable), so the load-bearing signal is the
    // df.filter.trueCount narrowing.
    await clickAggregateCell(page, 0);
    const after = await srcFilterCount(page);
    expect(after).toBeGreaterThan(0);
    expect(after).toBeLessThan(before);
  });

  await softStep('Scenario 2 Step 7: with SEX = "M" active, the click sets df.filter to the clicked group\'s FULL source count (replace, not intersect), monotonically (GROK-17726)', async () => {
    // GROK-17726 guard: Row Source = All cell filtering replaces the active source
   // filter with the clicked group's full row set. It must not preserve the previous
   // SEX='M' restriction by intersection.
    await setupSingleKey(page, 'All', true);
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const sex = df.col('SEX');
      df.filter.init((i: number) => sex.get(i) === 'M');   // SEX = 'M' active before the click
      await new Promise((r) => setTimeout(r, 400));
    });
    const preClick = await srcFilterCount(page);          // SEX = 'M' rows (baseline for GROK-17726)
    await clickAggregateCell(page, 0);
    const after = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const disPop = df.col('DIS_POP');
      const sex = df.col('SEX');
      // Which DIS_POP key survived the click, that key's FULL source-row count, and whether any
      // non-M rows survived (a replace keeps them; an intersect would not).
      const surviving = new Set<string>();
      let nonMSurviving = 0;
      for (let i = 0; i < df.rowCount; i++) if (df.filter.get(i)) {
        surviving.add(disPop.get(i));
        if (sex.get(i) !== 'M') nonMSurviving++;
      }
      let groupTotal = 0;
      if (surviving.size === 1) {
        const key = surviving.values().next().value;
        for (let i = 0; i < df.rowCount; i++) if (disPop.get(i) === key) groupTotal++;
      }
      return {count: df.filter.trueCount, groupTotal, keys: surviving.size, nonMSurviving};
    });
    expect(after.keys).toBe(1);                            // the click narrowed to a single DIS_POP key
    expect(after.count).toBeLessThanOrEqual(preClick);     // GROK-17726: monotonic narrowing, never grows
    expect(after.count).toBe(after.groupTotal);            // the clicked group's FULL source rows (replace)
    expect(after.nonMSurviving).toBeGreaterThan(0);        // non-M rows survived: the SEX filter was replaced, not intersected
  });

  await softStep('Scenario 2 Step 8: clear the source filter before the next scenario', async () => {
    await page.evaluate(async () => {
      grok.shell.tv.dataFrame.filter.setAll(true);
      await new Promise((r) => setTimeout(r, 400));
    });
    expect(await srcFilterCount(page)).toBe(DEMOG_ROWS);
  });

  // === Scenario 3: Negative path — the filter does NOT fire when conditions are not met (I1) ===

  await softStep('Scenario 3 Step 5: with filteringEnabled = false under All, a cell click leaves df.filter unchanged (I1 negative)', async () => {
    await setupSingleKey(page, 'All', false);
    const before = await srcFilterCount(page);
    await clickAggregateCell(page, 0);
    const after = await srcFilterCount(page);
    expect(after).toBe(before);                  // filteringEnabled=false: no filter fired
    expect(after).toBe(DEMOG_ROWS);
  });

  await softStep('Scenario 3 Step 9: at the default Filtered even with filteringEnabled = true, a cell click leaves df.filter unchanged (I1 negative)', async () => {
    await setupSingleKey(page, 'Filtered', true);
    const before = await srcFilterCount(page);
    await clickAggregateCell(page, 0);
    const after = await srcFilterCount(page);
    expect(after).toBe(before);                  // Row Source = Filtered: cell-driven filter inert
    expect(after).toBe(DEMOG_ROWS);
  });

  // === Scenario 4: Ctrl+click on an aggregate row drives the selection-to-selection link (channel 2, I3) ===

  await softStep('Scenario 4 Step 4: at Row Source = Filtered, a Ctrl+click on an aggregate row selects that group\'s source rows (df.selection)', async () => {
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      df.filter.setAll(true);
      await new Promise((r) => setTimeout(r, 300));
    });
    await setupSingleKey(page, 'Filtered', true);
    const beforeSel = await srcSelectionCount(page);
    expect(beforeSel).toBe(0);
    // Config signal (DOM caption read-back): the pivot is configured on the group the click
    // targets — Group by DIS_POP, Aggregate avg(AGE).
    expect(await pivotRowCaptions(page, 'Group by')).toEqual(['DIS_POP']);
    expect(await pivotRowCaptions(page, 'Aggregate')).toEqual(['avg(AGE)']);
    const filterBefore = await srcFilterCount(page);
    const warningsBefore = await page.evaluate(() => (grok.shell.warnings || []).length);
    // Ctrl+click the first aggregate row (row 0 = AS). The SELECTION_TO_SELECTION link selects that
    // group's source rows on df.selection; df.filter stays untouched (this is not the filter channel).
    await ctrlClickAggregateCell(page, 0);
    const after = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const disPop = df.col('DIS_POP');
      // Which DIS_POP key(s) got selected, and the selected-row count.
      const selected = new Set<string>();
      for (let i = 0; i < df.rowCount; i++) if (df.selection.get(i)) selected.add(disPop.get(i));
      return {sel: df.selection.trueCount, keys: [...selected], filt: df.filter.trueCount};
    });
    const warningsAfter = await page.evaluate(() => (grok.shell.warnings || []).length);
    expect(warningsAfter).toBe(warningsBefore);
    expect(after.sel).toBeGreaterThan(0);          // the Ctrl+click moved the source selection (I3 fired)
    expect(after.keys.length).toBe(1);             // exactly one DIS_POP group's rows were selected
    expect(after.sel).toBe(await groupSourceCount(page, after.keys[0]));  // the clicked group's FULL source rows
    expect(after.filt).toBe(filterBefore);         // the selection link left df.filter untouched
  });

  await softStep('Scenario 4 Step 9: at Row Source = All, a Ctrl+click on a different aggregate row re-selects the new group\'s source rows', async () => {
    await setupSingleKey(page, 'All', true);
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      df.filter.setAll(true);
      await new Promise((r) => setTimeout(r, 300));
    });
    // Config signal (DOM caption read-back): switching Row Source to All rebuilds the tag chips,
    // yet Group by DIS_POP and Aggregate avg(AGE) survive the rebuild.
    expect(await pivotRowCaptions(page, 'Group by')).toEqual(['DIS_POP']);
    expect(await pivotRowCaptions(page, 'Aggregate')).toEqual(['avg(AGE)']);
    const filterBefore = await srcFilterCount(page);
    const warningsBefore = await page.evaluate(() => (grok.shell.warnings || []).length);
    // Ctrl+click a DIFFERENT aggregate row (row 1 = Indigestion). Use Ctrl (not Shift) so the selection
    // channel fires in isolation — under Row Source = All a Shift+click would additionally trip the
    // cell-driven filter, muddling df.filter.
    await ctrlClickAggregateCell(page, 1);
    const after = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const disPop = df.col('DIS_POP');
      const selected = new Set<string>();
      for (let i = 0; i < df.rowCount; i++) if (df.selection.get(i)) selected.add(disPop.get(i));
      return {sel: df.selection.trueCount, keys: [...selected], filt: df.filter.trueCount};
    });
    const warningsAfter = await page.evaluate(() => (grok.shell.warnings || []).length);
    expect(warningsAfter).toBe(warningsBefore);
    expect(after.sel).toBeGreaterThan(0);          // the selection link remains active at Row Source = All
    expect(after.keys.length).toBe(1);             // exactly the newly clicked group's rows are selected
    expect(after.sel).toBe(await groupSourceCount(page, after.keys[0]));  // the new group's FULL source rows
    expect(after.filt).toBe(filterBefore);         // Ctrl+click kept the filter channel out of it
  });

  v.finishSpec();
});
