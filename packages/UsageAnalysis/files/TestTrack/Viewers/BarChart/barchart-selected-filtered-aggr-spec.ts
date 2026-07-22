/* ---
realizes: [barchart.int.selected-filtered-rows-need-cumulative-aggr]
--- */

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const splitCol = 'race';
const valueCol = 'age';

test('Bar Chart — Selected / Filtered Rows overlays with cumulative aggregations', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });
  const errCount = () => pageErrors.length + consoleErrors.length;

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 4000});
  await v.addViewerByIcon(page, 'bar-chart', 'Bar-chart');

  await page.evaluate(() => {
    const bcEl = document.querySelector('[name="viewer-Bar-chart"]') as HTMLElement;
    const panelBase = bcEl.closest('.panel-base') as HTMLElement;
    const gear = panelBase?.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
    gear?.click();
  });
  await page.waitForTimeout(500);

  const base = await page.evaluate(async ({split, value}) => {
    const df = grok.shell.tv.dataFrame;
    const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
    bc.props.splitColumnName = split;
    bc.props.valueColumnName = value;
    await new Promise((r) => setTimeout(r, 700));
    const race = df.col(split);
    const sex = df.col('sex');
    let asian = 0, caucasian = 0, female = 0;
    for (let i = 0; i < df.rowCount; i++) {
      if (race.get(i) === 'Asian') asian++;
      if (race.get(i) === 'Caucasian') caucasian++;
      if (sex.get(i) === 'F') female++;
    }
    return {rowCount: df.rowCount, asian, caucasian, female};
  }, {split: splitCol, value: valueCol});
  expect(base.asian).toBeGreaterThan(0);
  expect(base.caucasian).toBeGreaterThan(0);
  expect(base.female).toBeGreaterThan(0);
  expect(base.female).toBeLessThan(base.rowCount);

  await softStep('Scenario 1 Step 4: selecting Asian renders the Selected Rows overlay (orange hue on canvas); filter untouched', async () => {
    const errBefore = errCount();
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.valueAggrType = 'count';
      bc.props.showSelectedRows = true;
      bc.props.showFilteredRows = true;
      df.filter.setAll(true);
      df.selection.setAll(false);
      await new Promise((r) => setTimeout(r, 600));
    });
    await page.waitForTimeout(500);
    // No selection yet — the orange selection hue must be absent from the canvas.
    const huePxBefore = await v.countSelectionHuePixels(page, 'Bar chart');
    const info = await page.evaluate(async ({split}) => {
      const df = grok.shell.tv.dataFrame;
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const race = df.col(split);
      df.selection.init((i: number) => race.get(i) === 'Asian');
      await new Promise((r) => setTimeout(r, 600));
      let selInCat = 0;
      for (let i = 0; i < df.rowCount; i++)
        if (df.selection.get(i) && race.get(i) === 'Asian') selInCat++;
      return {
        aggr: bc.props.valueAggrType,
        sel: df.selection.trueCount,
        selInCat,
        filt: df.filter.trueCount,
        rowCount: df.rowCount,
        hasCanvas: !!bc.root.querySelector('canvas'),
      };
    }, {split: splitCol});
    await page.waitForTimeout(500);
    // Overlay signal: selection paints the DarkOrange (#ff8c00) hue onto the bars.
    const huePx = await v.countSelectionHuePixels(page, 'Bar chart');
    expect(huePxBefore).toBe(0); // also rejects a -1 fault masking the delta below
    expect(huePx).toBeGreaterThan(0);
    expect(huePx).toBeGreaterThan(huePxBefore);
    // Secondary state checks (selection/filter bitsets behind the overlay).
    expect(info.aggr).toBe('count');
    expect(info.sel).toBe(base.asian);
    expect(info.selInCat).toBe(info.sel);
    expect(info.filt).toBe(info.rowCount);
    expect(info.hasCanvas).toBe(true);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 6: filtering to sex=F re-renders the bars (canvas pixel delta); selection unchanged', async () => {
    const errBefore = errCount();
    // Enter the overlay mode first and let it settle, so the color snapshot
    // isolates the FILTER effect from the rowSource/showFilteredRows flip.
    await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.rowSource = 'All';
      bc.props.showFilteredRows = true;
      await new Promise((r) => setTimeout(r, 900));
    });
    expect(await v.snapshotCanvasColors(page, 'Bar chart')).toBe(true);
    // Stability precheck: prove the mode-flip repaint has settled BEFORE the
    // filter lands, so the measured delta below is the filter's effect, not a
    // late tail of the rowSource/showFilteredRows re-render. diffCanvasColors
    // replaces the snapshot, so the stable frame becomes the baseline.
    await page.waitForTimeout(400);
    const settle = await v.diffCanvasColors(page, 'Bar chart');
    expect(settle.deltaPx).toBeGreaterThanOrEqual(0);
    expect(settle.deltaPx).toBeLessThan(50);
    const info = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const sex = df.col('sex');
      df.filter.init((i: number) => sex.get(i) === 'F');
      await new Promise((r) => setTimeout(r, 900));
      return {
        rowSource: bc.props.rowSource,
        sel: df.selection.trueCount,
        filt: df.filter.trueCount,
        rowCount: df.rowCount,
        hasCanvas: !!bc.root.querySelector('canvas'),
      };
    });
    await page.waitForTimeout(500);
    // Filtered Rows overlay signal: under rowSource 'All' the bar geometry is
    // fixed — the overlay RECOLORS bar segments in place, so assert the
    // color-histogram delta, not the non-white total (which stays equal by
    // design here).
    const {deltaPx} = await v.diffCanvasColors(page, 'Bar chart');
    expect(deltaPx).toBeGreaterThan(100);
    expect(info.rowSource).toBe('All');
    expect(info.filt).toBe(base.female);
    expect(info.filt).toBeLessThan(info.rowCount);
    expect(info.sel).toBe(base.asian);
    expect(info.hasCanvas).toBe(true);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 8: switching aggregation to min removes the overlays; selection & filter state preserved', async () => {
    const errBefore = errCount();
    const info = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.valueAggrType = 'min';
      await new Promise((r) => setTimeout(r, 600));
      return {
        aggr: bc.props.valueAggrType,
        sel: df.selection.trueCount,
        filt: df.filter.trueCount,
        hasCanvas: !!bc.root.querySelector('canvas'),
      };
    });
    await page.waitForTimeout(500);
    // Non-cumulative aggregation removes the overlays: no orange selection hue
    // may remain on the canvas (0 expected; -1 fault also fails here).
    expect(await v.countSelectionHuePixels(page, 'Bar chart')).toBe(0);
    expect(info.aggr).toBe('min');
    expect(info.sel).toBe(base.asian);
    expect(info.filt).toBe(base.female);
    expect(info.hasCanvas).toBe(true);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 10: re-select Asian + re-filter F under sum restores the Selected Rows overlay', async () => {
    const errBefore = errCount();
    const info = await page.evaluate(async ({split}) => {
      const df = grok.shell.tv.dataFrame;
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const race = df.col(split);
      const sex = df.col('sex');
      df.selection.init((i: number) => race.get(i) === 'Asian');
      df.filter.init((i: number) => sex.get(i) === 'F');
      bc.props.valueAggrType = 'sum';
      await new Promise((r) => setTimeout(r, 600));
      return {
        aggr: bc.props.valueAggrType,
        sel: df.selection.trueCount,
        filt: df.filter.trueCount,
        rowCount: df.rowCount,
        hasCanvas: !!bc.root.querySelector('canvas'),
      };
    }, {split: splitCol});
    await page.waitForTimeout(500);
    // Cumulative aggregation restores the overlays: the orange selection hue
    // reappears on the bars.
    expect(await v.countSelectionHuePixels(page, 'Bar chart')).toBeGreaterThan(0);
    expect(info.aggr).toBe('sum');
    expect(info.sel).toBe(base.asian);
    expect(info.filt).toBe(base.female);
    expect(info.filt).toBeLessThan(info.rowCount);
    expect(info.hasCanvas).toBe(true);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 4: selecting Asian+Caucasian renders the Selected Rows overlay under values aggr; filter unchanged', async () => {
    const errBefore = errCount();
    const info = await page.evaluate(async ({split}) => {
      const df = grok.shell.tv.dataFrame;
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.valueAggrType = 'values';
      bc.props.showSelectedRows = true;
      const race = df.col(split);
      df.selection.init((i: number) => race.get(i) === 'Asian' || race.get(i) === 'Caucasian');
      await new Promise((r) => setTimeout(r, 600));
      let selInCats = 0;
      for (let i = 0; i < df.rowCount; i++)
        if (df.selection.get(i) && (race.get(i) === 'Asian' || race.get(i) === 'Caucasian')) selInCats++;
      return {
        aggr: bc.props.valueAggrType,
        sel: df.selection.trueCount,
        selInCats,
        filt: df.filter.trueCount,
        hasCanvas: !!bc.root.querySelector('canvas'),
      };
    }, {split: splitCol});
    await page.waitForTimeout(500);
    // Selected Rows overlay under value-count (cumulative): orange hue present.
    expect(await v.countSelectionHuePixels(page, 'Bar chart')).toBeGreaterThan(0);
    expect(info.aggr).toBe('values');
    expect(info.sel).toBe(base.asian + base.caucasian);
    expect(info.selInCats).toBe(info.sel);
    expect(info.filt).toBe(base.female);
    expect(info.hasCanvas).toBe(true);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 6: filtering under values aggregation re-renders the Filtered Rows overlay (canvas color delta); selection unchanged', async () => {
    const errBefore = errCount();
    // Enter the overlay mode and let it settle so the color snapshot isolates
    // the FILTER effect from the rowSource/showFilteredRows repaint (mirrors
    // Scenario 1 Step 6). This asserts the Filtered Rows overlay under the
    // values (value-count) aggregation, before switching to avg.
    await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.rowSource = 'All';
      bc.props.showFilteredRows = true;
      await new Promise((r) => setTimeout(r, 900));
    });
    expect(await v.snapshotCanvasColors(page, 'Bar chart')).toBe(true);
    await page.waitForTimeout(400);
    const settle = await v.diffCanvasColors(page, 'Bar chart');
    expect(settle.deltaPx).toBeGreaterThanOrEqual(0);
    // Ceiling 40 sits well below the filter recolor signal floor 120.
    expect(settle.deltaPx).toBeLessThan(40);
    const info = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const sex = df.col('sex');
      df.filter.init((i: number) => sex.get(i) === 'M');
      await new Promise((r) => setTimeout(r, 900));
      return {
        aggr: bc.props.valueAggrType,
        rowSource: bc.props.rowSource,
        sel: df.selection.trueCount,
        filt: df.filter.trueCount,
        rowCount: df.rowCount,
        hasCanvas: !!bc.root.querySelector('canvas'),
      };
    });
    await page.waitForTimeout(500);
    // Under rowSource 'All' the filter RECOLORS the bar segments in place, so
    // assert the per-color histogram delta, not the non-white total.
    const {deltaPx} = await v.diffCanvasColors(page, 'Bar chart');
    expect(deltaPx).toBeGreaterThan(120);
    expect(info.aggr).toBe('values');
    expect(info.rowSource).toBe('All');
    expect(info.filt).toBeLessThan(info.rowCount);
    expect(info.filt).toBeGreaterThan(0);
    expect(info.sel).toBe(base.asian + base.caucasian);
    expect(info.hasCanvas).toBe(true);
    expect(errCount()).toBe(errBefore);
    // Restore the sex=F filter so the following avg step sees the documented
    // filtered-row count.
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const sex = df.col('sex');
      df.filter.init((i: number) => sex.get(i) === 'F');
      await new Promise((r) => setTimeout(r, 600));
    });
  });

  await softStep('Scenario 2 Step 7: switching aggregation to avg removes the overlays; selection & filter state preserved', async () => {
    const errBefore = errCount();
    const info = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.showFilteredRows = true;
      bc.props.valueAggrType = 'avg';
      await new Promise((r) => setTimeout(r, 600));
      return {
        aggr: bc.props.valueAggrType,
        sel: df.selection.trueCount,
        filt: df.filter.trueCount,
        hasCanvas: !!bc.root.querySelector('canvas'),
      };
    });
    await page.waitForTimeout(500);
    // Non-cumulative avg removes the overlays: no orange selection hue on the
    // canvas (0 expected; -1 fault also fails here).
    expect(await v.countSelectionHuePixels(page, 'Bar chart')).toBe(0);
    expect(info.aggr).toBe('avg');
    expect(info.sel).toBe(base.asian + base.caucasian);
    expect(info.filt).toBe(base.female);
    expect(info.hasCanvas).toBe(true);
    expect(errCount()).toBe(errBefore);
  });

  await page.evaluate(async () => {
    const tv = grok.shell.tv;
    const df = tv.dataFrame;
    const bc = Array.from(tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
    bc.props.showSelectedRows = false;
    bc.props.showFilteredRows = false;
    df.selection.setAll(false);
    df.filter.setAll(true);
    await new Promise((r) => setTimeout(r, 300));
  });

  v.finishSpec();
});
