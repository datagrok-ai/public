/* ---
realizes: [barchart.int.stack-needs-additive-aggr, barchart.int.datetime-split-enables-map]
--- */

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

test('Bar Chart — Stack Aggregation Precondition and DateTime Split Map', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });
  const errCount = () => pageErrors.length + consoleErrors.length;

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'bar-chart', 'Bar-chart');

  await page.evaluate(() => {
    const bcEl = document.querySelector('[name="viewer-Bar-chart"]') as HTMLElement;
    const panelBase = bcEl.closest('.panel-base') as HTMLElement;
    const gear = panelBase.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
    gear.click();
  });
  await page.waitForTimeout(500);

  await page.evaluate(async () => {
    const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
    bc.props.splitColumnName = 'RACE';
    await new Promise((r) => setTimeout(r, 400));
  });

  await softStep('Scenario 1 Step 4: Stack + avg aggregation suppresses the legend column', async () => {
    const errBefore = errCount();
    await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.stackColumnName = 'SEX';
      bc.props.valueColumnName = 'AGE';
      bc.props.valueAggrType = 'avg';
      await new Promise((r) => setTimeout(r, 700));
    });
    const legendAvg = await v.readLegend(page, 'Bar chart');
    expect(legendAvg.legendRendered).toBe(false);
    expect(legendAvg.itemCount).toBe(0);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 5: Stack + sum aggregation produces a visible legend', async () => {
    const errBefore = errCount();
    await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.valueAggrType = 'sum';
      await new Promise((r) => setTimeout(r, 700));
    });
    const legendSum = await v.readLegend(page, 'Bar chart');
    expect(legendSum.legendRendered).toBe(true);
    expect(legendSum.itemCount).toBeGreaterThan(0);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 7: count keeps the legend; removing Stack collapses it', async () => {
    const errBefore = errCount();
    const legendCount = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.valueAggrType = 'count';
      await new Promise((r) => setTimeout(r, 700));
      const root = bc.root;
      return {
        rendered: !!root.querySelector('[name="legend"]'),
        itemCount: root.querySelectorAll('[name="legend"] .d4-legend-item').length,
      };
    });
    expect(legendCount.rendered).toBe(true);
    expect(legendCount.itemCount).toBeGreaterThan(0);

    const legendNoStack = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.stackColumnName = '';
      await new Promise((r) => setTimeout(r, 700));
      const root = bc.root;
      return {
        stack: bc.props.stackColumnName,
        rendered: !!root.querySelector('[name="legend"]'),
        itemCount: root.querySelectorAll('[name="legend"] .d4-legend-item').length,
      };
    });
    expect(legendNoStack.stack).toBe('');
    expect(legendNoStack.rendered).toBe(false);
    expect(legendNoStack.itemCount).toBe(0);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 9: String Split column hides the Split Map selector', async () => {
    const errBefore = errCount();
    const strMap = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.valueColumnName = 'AGE';
      bc.props.valueAggrType = 'count';
      bc.props.splitColumnName = 'RACE';
      grok.shell.o = bc;
      await new Promise((r) => setTimeout(r, 1200));
      const splitType = bc.dataFrame.col('RACE').type;
      const host = document.querySelector('[name="input-aggr-selector-split-map"]') as HTMLElement | null;
      // A host absent from the DOM is a legitimate "selector hidden" outcome.
      if (!host) return {present: false, visible: false, splitType};
      const rect = host.getBoundingClientRect();
      const display = getComputedStyle(host).display;
      return {
        present: true,
        visible: rect.width > 0 && rect.height > 0 && display !== 'none' && host.offsetParent !== null,
        splitType,
      };
    });
    expect(strMap.splitType).toBe('string');
    expect(strMap.visible).toBe(false);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 12: DateTime Split column enables the Split Map selector', async () => {
    const errBefore = errCount();
    const dtMap = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'STARTED';
      bc.props.splitMap = 'year';
      grok.shell.o = bc;
      await new Promise((r) => setTimeout(r, 1200));
      const host = document.querySelector('[name="input-aggr-selector-split-map"]') as HTMLElement | null;
      if (!host) return {present: false, visible: false};
      const rect = host.getBoundingClientRect();
      const display = getComputedStyle(host).display;
      return {
        present: true,
        visible: rect.width > 0 && rect.height > 0 && display !== 'none' && host.offsetParent !== null,
        splitType: bc.dataFrame.col('STARTED').type,
        splitMap: bc.props.splitMap,
      };
    });
    expect(dtMap.splitType).toBe('datetime');
    expect(dtMap.visible).toBe(true);
    expect(dtMap.splitMap).toBe('year');
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 13: the in-chart Category selector reflects the STARTED split column', async () => {
    const errBefore = errCount();
    // The effective split column name is shown in the in-chart Category
    // selector. The Split Map suffix (e.g. "STARTED (Year)") is canvas-rendered
    // on the axis header and not exposed to the DOM, so asserting the "(Year)"
    // suffix is a documented reduction; the reachable signal is the selector
    // text reflecting the STARTED column.
    const sel = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const texts = Array.from(bc.root.querySelectorAll('.d4-column-selector-column'))
        .map((e) => (e.textContent || '').trim());
      return {texts, splitCol: bc.props.splitColumnName};
    });
    expect(sel.splitCol).toBe('STARTED');
    expect(sel.texts.some((t) => t.includes('STARTED'))).toBe(true);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 14-15: Split Map re-categorizes Year → Month → Year (distinct-count + canvas delta)', async () => {
    const errBefore = errCount();
    // Distinct month buckets outnumber distinct year buckets in demog.csv, so a
    // month map produces more categories than a year map — the re-categorization
    // signal. The map change also repaints the chart (canvas color delta).
    const dist = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col('STARTED');
      const years = new Set<number>();
      const months = new Set<string>();
      for (let i = 0; i < df.rowCount; i++) {
        const d = col.get(i);
        if (d == null) continue;
        const dd = new Date(d);
        years.add(dd.getFullYear());
        months.add(dd.getFullYear() + '-' + dd.getMonth());
      }
      return {yearN: years.size, monthN: months.size};
    });
    expect(dist.monthN).toBeGreaterThan(dist.yearN);
    expect(await v.snapshotCanvasColors(page, 'Bar chart')).toBe(true);
    await page.waitForTimeout(400);
    // Settle-precheck: the pre-map frame is quiescent. Ceiling 300 sits far
    // below the map-change signal floor 5000, so a residual repaint tail cannot
    // masquerade as the re-categorization.
    const settle = await v.diffCanvasColors(page, 'Bar chart');
    expect(settle.deltaPx).toBeGreaterThanOrEqual(0); // -1 = canvas fault
    expect(settle.deltaPx).toBeLessThan(300);
    const monthMap = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.splitMap = 'month';
      await new Promise((r) => setTimeout(r, 900));
      return {splitMap: bc.props.splitMap, hasCanvas: !!bc.root.querySelector('canvas')};
    });
    await page.waitForTimeout(300);
    // Year→Month re-categorization repaints the whole plot. Floor 5000 is an
    // order of magnitude above the settle ceiling.
    const monthDelta = await v.diffCanvasColors(page, 'Bar chart');
    const yearMap = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.splitMap = 'year';
      await new Promise((r) => setTimeout(r, 900));
      return {splitMap: bc.props.splitMap, hasCanvas: !!bc.root.querySelector('canvas'), splitCol: bc.props.splitColumnName};
    });
    // Month→Year revert re-categorizes back — a real render signal for the
    // year-revert leg.
    const yearDelta = await v.diffCanvasColors(page, 'Bar chart');
    expect(monthMap.splitMap).toBe('month');
    expect(monthMap.hasCanvas).toBe(true);
    expect(monthDelta.deltaPx).toBeGreaterThan(5000);
    expect(yearMap.splitMap).toBe('year');
    expect(yearMap.hasCanvas).toBe(true);
    expect(yearMap.splitCol).toBe('STARTED');
    expect(yearDelta.deltaPx).toBeGreaterThan(5000);
    expect(errCount()).toBe(errBefore);
  });

  v.finishSpec();
});
