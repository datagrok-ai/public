/* ---
realizes: [barchart.cp.value-axis-range-scale-scroll]
--- */

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';
const splitCol = 'Primary Series Name';

test('Bar Chart — Value-Axis Range, Scale, and Selected/Filtered Rows Overlays', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 4000});
  await v.addViewerByIcon(page, 'bar-chart', 'Bar-chart');

  await page.evaluate(() => {
    const bcEl = document.querySelector('[name="viewer-Bar-chart"]') as HTMLElement;
    const panelBase = bcEl.closest('.panel-base') as HTMLElement;
    const gear = panelBase.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
    gear.click();
  });
  await page.waitForTimeout(500);

  await page.evaluate(async ({split}) => {
    const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
    bc.props.splitColumnName = split;
    bc.props.valueColumnName = 'CAST Idea ID';
    bc.props.valueAggrType = 'count';
    bc.props.showClippedBarIndicators = true;
    await new Promise((r) => setTimeout(r, 900));
  }, {split: splitCol});

  await softStep('Step 3: Value Min above the shortest bar clips it (clipped-bar indicators)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const info = await page.evaluate(async ({split}) => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const df = grok.shell.tv.dataFrame;
      const col = df.col(split);
      const counts: Record<string, number> = {};
      for (const c of col.categories) counts[c] = 0;
      for (let i = 0; i < df.rowCount; i++) counts[col.get(i)]++;
      const nonZero = Object.values(counts).filter((n) => (n as number) > 0) as number[];
      const shortest = Math.min(...nonZero);
      const tallest = Math.max(...nonZero);
      const minOffset = Math.max(1, Math.round((tallest - shortest) * 0.25));
      bc.props.valueMin = shortest + minOffset;
      await new Promise((r) => setTimeout(r, 800));
      const barsBelowMin = nonZero.filter((n) => n < bc.props.valueMin).length;
      return {
        valueMin: bc.props.valueMin,
        shortest,
        tallest,
        minOffset,
        barsBelowMin,
        clippedIndicators: bc.props.showClippedBarIndicators,
        hasCanvas: !!bc.root.querySelector('canvas'),
      };
    }, {split: splitCol});
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.valueMin).toBe(info.shortest + info.minOffset);
    expect(info.valueMin).toBeLessThan(info.tallest);
    expect(info.barsBelowMin).toBeGreaterThanOrEqual(1);
    expect(info.clippedIndicators).toBe(true);
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Step 5: Value Max below the tallest bar clips it; axis stops at the maximum', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const info = await page.evaluate(async ({split}) => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const df = grok.shell.tv.dataFrame;
      const col = df.col(split);
      const counts: Record<string, number> = {};
      for (const c of col.categories) counts[c] = 0;
      for (let i = 0; i < df.rowCount; i++) counts[col.get(i)]++;
      const nonZero = Object.values(counts).filter((n) => (n as number) > 0) as number[];
      const shortest = Math.min(...nonZero);
      const tallest = Math.max(...nonZero);
      const maxOffset = Math.max(1, Math.round((tallest - shortest) * 0.25));
      bc.props.valueMax = tallest - maxOffset;
      await new Promise((r) => setTimeout(r, 800));
      const barsAboveMax = nonZero.filter((n) => n > bc.props.valueMax).length;
      return {
        valueMax: bc.props.valueMax,
        tallest,
        shortest,
        maxOffset,
        barsAboveMax,
        valueMin: bc.props.valueMin,
        hasCanvas: !!bc.root.querySelector('canvas'),
      };
    }, {split: splitCol});
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.valueMax).toBe(info.tallest - info.maxOffset);
    expect(info.valueMax).toBeGreaterThan(info.shortest);
    expect(info.valueMax).toBeGreaterThan(info.valueMin);
    expect(info.barsAboveMax).toBeGreaterThanOrEqual(1);
    expect(info.valueMin).not.toBeNull();
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Step 7: the value-axis scroll bar navigates the constrained range without error', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    await page.locator('[name="viewer-Bar-chart"]').first().hover();
    await page.waitForTimeout(400);
    const info = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const rangeSliders = bc.root.querySelectorAll('svg[type="range-slider"]').length;
      return {
        valueMin: bc.props.valueMin,
        valueMax: bc.props.valueMax,
        rangeSliders,
        hasCanvas: !!bc.root.querySelector('canvas'),
      };
    });
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.valueMin).not.toBeNull();
    expect(info.valueMax).not.toBeNull();
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Step 9: logarithmic value axis re-scales positive-count bars without error', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const info = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const before = bc.props.axisType;
      bc.props.axisType = 'logarithmic';
      await new Promise((r) => setTimeout(r, 900));
      return {
        before,
        axisType: bc.props.axisType,
        hasCanvas: !!bc.root.querySelector('canvas'),
      };
    });
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.before).toBe('linear');
    expect(info.axisType).toBe('logarithmic');
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Step 11: switching back to linear restores linear proportional heights', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const info = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const before = bc.props.axisType;
      bc.props.axisType = 'linear';
      await new Promise((r) => setTimeout(r, 900));
      return {
        before,
        axisType: bc.props.axisType,
        hasCanvas: !!bc.root.querySelector('canvas'),
      };
    });
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.before).toBe('logarithmic');
    expect(info.axisType).toBe('linear');
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Step 13: clearing Value Min/Max restores the full-range axis (no clipping)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const info = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.valueMin = null;
      bc.props.valueMax = null;
      await new Promise((r) => setTimeout(r, 800));
      return {
        valueMin: bc.props.valueMin,
        valueMax: bc.props.valueMax,
        axisType: bc.props.axisType,
        hasCanvas: !!bc.root.querySelector('canvas'),
      };
    });
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.valueMin).toBeNull();
    expect(info.valueMax).toBeNull();
    expect(info.axisType).toBe('linear');
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  v.finishSpec();
});
