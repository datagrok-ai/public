/* ---
realizes: [pivottable.cp.rowsource-filter-selection-links, pivottable.int.filtering-requires-row-source-all, pivottable.int.selection-link-gated-by-row-source]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {armBalloonRecorderProved, expectNoBalloonSinceArmed} from '../../helpers/balloons';

declare const grok: any;

test.use(specTestOptions);

const PIVOT = '[name="viewer-Pivot-table"]';
const DEMOG_ROWS = 5850;

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

async function setupSingleKey(page: Page, rowSource: string, filteringEnabled: boolean) {
  await v.setViewerProps(page, 'Pivot table', [{set: {
    groupByColumnNames: ['DIS_POP'],
    pivotColumnNames: [],
    aggregateColumnNames: ['AGE'],
    aggregateAggTypes: ['avg'],
    rowSource: rowSource,
    filteringEnabled: filteringEnabled,
  }}], 700);
}

async function waitForDfEvent(page: Page, events: string[], capMs: number) {
  await page.evaluate(({evs, cap}) => new Promise<void>((resolve) => {
    const df = grok.shell.tv.dataFrame;
    const subs: any[] = [];
    const done = () => {
      for (const s of subs) { try { s.unsubscribe(); } catch (_) {} }
      resolve();
    };
    for (const ev of evs) {
      try { subs.push(df[ev].subscribe(() => done())); } catch (_) {  }
    }
    setTimeout(done, cap);
  }), {evs: events, cap: capMs});
}

async function innerGridRect(page: Page) {
  const box = await page.locator(`${PIVOT} [name="canvas"]`).first().boundingBox();
  if (!box) throw new Error('pivot inner grid canvas not visible');
  return box;
}

async function clickAggregateCell(page: Page, rowIdx: number) {
  const r = await innerGridRect(page);
  await page.mouse.click(r.x + 55, r.y + 24 + rowIdx * 24 + 12);
  await waitForDfEvent(page, ['onRowsFiltered', 'onFilterChanged'], 500);
}

async function ctrlClickAggregateCell(page: Page, rowIdx: number) {
  const r = await innerGridRect(page);
  await page.keyboard.down('Control');
  await page.mouse.click(r.x + 55, r.y + 24 + rowIdx * 24 + 12);
  await page.keyboard.up('Control');
  await waitForDfEvent(page, ['onSelectionChanged'], 500);
}

async function groupSourceCount(page: Page, disPopKey: string) {
  return page.evaluate((k: string) => {
    const df = grok.shell.tv.dataFrame;
    const col = df.col('DIS_POP');
    let n = 0;
    for (let i = 0; i < df.rowCount; i++) if (col.get(i) === k) n++;
    return n;
  }, disPopKey);
}

const srcFilterCount = (page: Page) => page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
const srcSelectionCount = (page: Page) => page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);

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
  await v.waitForViewerRendered(page, 'Pivot table', 800);

  await softStep('Setup: tag-editor header shows Group by, Aggregate and Pivot rows', async () => {
    const titles = await page.evaluate(() => Array.from(
      document.querySelectorAll('[name="viewer-Pivot-table"] .grok-pivot-column-tags-title'))
      .map((t) => t.getAttribute('d4-name')));
    expect(titles).toEqual(expect.arrayContaining(['Group by', 'Aggregate', 'Pivot']));

    await setupSingleKey(page, 'Filtered', true);
  });

  await softStep('Scenario 1 Step 4: after a SEX = "M" panel filter, the pivot aggregate matches an independent groupBy on the filtered df', async () => {

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

    await page.evaluate(() => grok.shell.tv.getFiltersGroup());
    await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 15000});
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const sex = df.col('SEX');
      df.filter.init((i: number) => sex.get(i) === 'M');
    });
    await v.waitForViewerRendered(page, 'Pivot table', 900);

    const filtered = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;

      const g = df.groupBy(['DIS_POP']).avg('AGE').whereRowMask(df.filter).aggregate();
      const avgCol = Array.from({length: g.columns.length}, (_: any, i: number) => g.columns.byIndex(i).name)
        .find((c: string) => c.toLowerCase().includes('age'));
      const out: Record<string, number> = {};
      const key = g.col('DIS_POP');
      for (let i = 0; i < g.rowCount; i++) out[key.get(i)] = g.col(avgCol!).get(i);
      return out;
    });

    const changed = Object.keys(filtered).some((k) => baseline[k] !== undefined
      && Math.abs(filtered[k] - baseline[k]) > 1e-9);
    expect(changed).toBe(true);

    await page.locator(`${PIVOT} .grok-pivot-counts [name="button-ADD"]`).click();
    await v.pollValue(() => page.evaluate(() =>
      Array.from(grok.shell.views).some((vw: any) => vw.name === 'Table aggregation')),
    (present) => present, 1000, 100);
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
    await page.evaluate(() => {
      const aggView = Array.from(grok.shell.views).find((vw: any) => vw.name === 'Table aggregation') as any;
      if (aggView) aggView.close();
      grok.shell.tv.dataFrame.filter.setAll(true);   
    });
    await v.pollValue(() => srcFilterCount(page), (n) => n === DEMOG_ROWS, 400, 100);
    expect(await srcFilterCount(page)).toBe(DEMOG_ROWS);
  });

  await softStep('Scenario 2 Step 3: Row Source = All with filteringEnabled = true does not widen the aggregated row set', async () => {

    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const sex = df.col('SEX');
      df.filter.init((i: number) => sex.get(i) === 'M');
    });
    await v.waitForViewerRendered(page, 'Pivot table', 400);
    const rowsBefore = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return df.groupBy(['DIS_POP']).avg('AGE').whereRowMask(df.filter).aggregate().rowCount;
    });

    await setupSingleKey(page, 'All', true);
    const props = await pivotProps(page);
    expect(props.rowSource).toBe('All');
    expect(props.filteringEnabled).toBe(true);

    const rowsAfter = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return df.groupBy(['DIS_POP']).avg('AGE').whereRowMask(df.filter).aggregate().rowCount;
    });
    expect(rowsAfter).toBe(rowsBefore);
  });

  await softStep('Scenario 2 Step 5: a real cell click narrows df.filter.trueCount to a smaller positive value', async () => {
    const before = await srcFilterCount(page);

    await clickAggregateCell(page, 0);

    const after = await v.pollValue(() => srcFilterCount(page), (n) => n > 0 && n < before, 3000, 150);
    expect(after).toBeGreaterThan(0);
    expect(after).toBeLessThan(before);
  });

  await softStep('Scenario 2 Step 7: with SEX = "M" active, the click sets df.filter to the clicked group\'s FULL source count (replace, not intersect), monotonically (GROK-17726)', async () => {

    await setupSingleKey(page, 'All', true);
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const sex = df.col('SEX');
      df.filter.init((i: number) => sex.get(i) === 'M');   
    });
    await v.waitForViewerRendered(page, 'Pivot table', 400);
    const preClick = await srcFilterCount(page);          
    await clickAggregateCell(page, 0);

    const after = await v.pollValue(() => page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const disPop = df.col('DIS_POP');
      const sex = df.col('SEX');

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
    }), (s) => s.keys === 1, 3000, 150);
    expect(after.keys).toBe(1);                            
    expect(after.count).toBeLessThanOrEqual(preClick);     
    expect(after.count).toBe(after.groupTotal);            
    expect(after.nonMSurviving).toBeGreaterThan(0);        
  });

  await softStep('Scenario 2 Step 8: clear the source filter before the next scenario', async () => {
    await page.evaluate(() => grok.shell.tv.dataFrame.filter.setAll(true));
    await v.pollValue(() => srcFilterCount(page), (n) => n === DEMOG_ROWS, 400, 100);
    expect(await srcFilterCount(page)).toBe(DEMOG_ROWS);
  });

  await softStep('Scenario 3 Step 5: with filteringEnabled = false under All, a cell click leaves df.filter unchanged (I1 negative)', async () => {
    await setupSingleKey(page, 'All', false);
    const before = await srcFilterCount(page);
    await clickAggregateCell(page, 0);
    const after = await srcFilterCount(page);
    expect(after).toBe(before);                  
    expect(after).toBe(DEMOG_ROWS);
  });

  await softStep('Scenario 3 Step 9: at the default Filtered even with filteringEnabled = true, a cell click leaves df.filter unchanged (I1 negative)', async () => {
    await setupSingleKey(page, 'Filtered', true);
    const before = await srcFilterCount(page);
    await clickAggregateCell(page, 0);
    const after = await srcFilterCount(page);
    expect(after).toBe(before);                  
    expect(after).toBe(DEMOG_ROWS);
  });

  await softStep('Scenario 4 Step 4: at Row Source = Filtered, a Ctrl+click on an aggregate row selects that group\'s source rows (df.selection)', async () => {
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      df.filter.setAll(true);
    });
    await v.pollValue(() => page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {sel: df.selection.trueCount, filt: df.filter.trueCount};
    }), (s) => s.sel === 0 && s.filt === DEMOG_ROWS, 300, 100);
    await setupSingleKey(page, 'Filtered', true);
    const beforeSel = await srcSelectionCount(page);
    expect(beforeSel).toBe(0);

    expect(await pivotRowCaptions(page, 'Group by')).toEqual(['DIS_POP']);
    expect(await pivotRowCaptions(page, 'Aggregate')).toEqual(['avg(AGE)']);
    const filterBefore = await srcFilterCount(page);
    // grok.shell.warnings does not exist (js-api/src/shell.ts:176 exposes only warning()), so a
    // before/after count on it is not an alternative — it compares undefined with undefined. The
    // probe is raised INSIDE the window this arm opens, so an empty reading afterwards means the
    // click raised nothing rather than that nothing was watching; the probe's own balloon carries a
    // marker the absence read filters out.
    await armBalloonRecorderProved(page, 'pivottable ctrl-click');

    await ctrlClickAggregateCell(page, 0);

    const after = await v.pollValue(() => page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const disPop = df.col('DIS_POP');

      const selected = new Set<string>();
      for (let i = 0; i < df.rowCount; i++) if (df.selection.get(i)) selected.add(disPop.get(i));
      return {sel: df.selection.trueCount, keys: [...selected], filt: df.filter.trueCount};
    }), (s) => s.sel > 0 && s.keys.length === 1, 3000, 150);
    await expectNoBalloonSinceArmed(page, 'a Ctrl+click on an aggregate row', /warning/);
    expect(after.sel).toBeGreaterThan(0);          
    expect(after.keys.length).toBe(1);             
    expect(after.sel).toBe(await groupSourceCount(page, after.keys[0]));  
    expect(after.filt).toBe(filterBefore);         
  });

  await softStep('Scenario 4 Step 9: at Row Source = All, a Ctrl+click on a different aggregate row re-selects the new group\'s source rows', async () => {
    await setupSingleKey(page, 'All', true);
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      df.filter.setAll(true);
    });
    await v.pollValue(() => page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {sel: df.selection.trueCount, filt: df.filter.trueCount};
    }), (s) => s.sel === 0 && s.filt === DEMOG_ROWS, 300, 100);

    expect(await pivotRowCaptions(page, 'Group by')).toEqual(['DIS_POP']);
    expect(await pivotRowCaptions(page, 'Aggregate')).toEqual(['avg(AGE)']);
    const filterBefore = await srcFilterCount(page);
    await armBalloonRecorderProved(page, 'pivottable ctrl-click');

    await ctrlClickAggregateCell(page, 1);
    const after = await v.pollValue(() => page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const disPop = df.col('DIS_POP');
      const selected = new Set<string>();
      for (let i = 0; i < df.rowCount; i++) if (df.selection.get(i)) selected.add(disPop.get(i));
      return {sel: df.selection.trueCount, keys: [...selected], filt: df.filter.trueCount};
    }), (s) => s.sel > 0 && s.keys.length === 1, 3000, 150);
    await expectNoBalloonSinceArmed(page, 'a Ctrl+click on an aggregate row', /warning/);
    expect(after.sel).toBeGreaterThan(0);          
    expect(after.keys.length).toBe(1);             
    expect(after.sel).toBe(await groupSourceCount(page, after.keys[0]));  
    expect(after.filt).toBe(filterBefore);         
  });

  v.finishSpec();
});
