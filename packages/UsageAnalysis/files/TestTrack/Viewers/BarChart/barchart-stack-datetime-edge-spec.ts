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

  const gearBox = await page.evaluate(() => {
    const bcEl = document.querySelector('[name="viewer-Bar-chart"]') as HTMLElement;
    const panelBase = bcEl.closest('.panel-base') as HTMLElement;
    const gear = panelBase.querySelector('.panel-titlebar [name="icon-font-icon-settings"]') as HTMLElement;
    const r = gear.getBoundingClientRect();
    return {x: r.x + r.width / 2, y: r.y + r.height / 2};
  });
  await page.mouse.click(gearBox.x, gearBox.y);
  await v.pollValue(() => page.locator('.property-grid').count(), (n) => n > 0, 500, 100);

  await v.setViewerProps(page, 'Bar chart', [{set: {splitColumnName: 'RACE'}, wait: 400}]);

  await softStep('Scenario 1 Step 4: Stack + avg aggregation suppresses the legend column', async () => {
    const errBefore = errCount();
    await v.setViewerProps(page, 'Bar chart', [{set: {
      stackColumnName: 'SEX', valueColumnName: 'AGE', valueAggrType: 'avg',
    }, wait: 700}]);

    const legendAvg = await v.readLegend(page, 'Bar chart');
    expect(legendAvg.legendRendered).toBe(false);
    expect(legendAvg.itemCount).toBe(0);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 5: Stack + sum aggregation produces a visible legend', async () => {
    const errBefore = errCount();
    await v.setViewerProps(page, 'Bar chart', [{set: {valueAggrType: 'sum'}, wait: 700}]);
    const legendSum = await v.pollValue(() => v.readLegend(page, 'Bar chart'),
      (l) => l.legendRendered && l.itemCount > 0, 3000, 150);
    expect(legendSum.legendRendered).toBe(true);
    expect(legendSum.itemCount).toBeGreaterThan(0);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 7: count keeps the legend; removing Stack collapses it', async () => {
    const errBefore = errCount();
    await v.setViewerProps(page, 'Bar chart', [{set: {valueAggrType: 'count'}, wait: 700}]);
    const legendCount = await v.pollValue(() => v.readLegend(page, 'Bar chart'),
      (l) => l.laidOut && l.itemCount > 0, 3000, 150);
    expect(legendCount.laidOut).toBe(true);
    expect(legendCount.itemCount).toBeGreaterThan(0);

    await v.setViewerProps(page, 'Bar chart', [{set: {stackColumnName: ''}, wait: 700}]);

    const legendNoStack = await v.pollValue(() => v.readLegend(page, 'Bar chart'),
      (l) => !l.laidOut && l.itemCount === 0, 3000, 150);
    const stackAfter = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      return bc.props.stackColumnName;
    });
    expect(stackAfter).toBe('');
    expect(legendNoStack.laidOut).toBe(false);
    expect(legendNoStack.itemCount).toBe(0);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 9: String Split column hides the Split Map selector', async () => {
    const errBefore = errCount();
    await v.setViewerProps(page, 'Bar chart', [{set: {
      valueColumnName: 'AGE', valueAggrType: 'count', splitColumnName: 'RACE',
    }, wait: 1200}]);
    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      grok.shell.o = bc;
    });

    const strMap = await v.pollValue(() => page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const splitType = bc.dataFrame.col('RACE').type;
      const host = document.querySelector('[name="input-aggr-selector-split-map"]') as HTMLElement | null;

      if (!host) return {present: false, visible: false, splitType};
      const rect = host.getBoundingClientRect();
      const display = getComputedStyle(host).display;
      return {
        present: true,
        visible: rect.width > 0 && rect.height > 0 && display !== 'none' && host.offsetParent !== null,
        splitType,
      };
    }), (m) => m.present, 1200, 100);
    expect(strMap.splitType).toBe('string');
    expect(strMap.visible).toBe(false);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 12: DateTime Split column enables the Split Map selector', async () => {
    const errBefore = errCount();
    await v.setViewerProps(page, 'Bar chart', [{set: {
      splitColumnName: 'STARTED', splitMap: 'year',
    }, wait: 1200}]);
    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      grok.shell.o = bc;
    });
    const dtMap = await v.pollValue(() => page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
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
    }), (m) => m.visible, 1200, 100);
    expect(dtMap.splitType).toBe('datetime');
    expect(dtMap.visible).toBe(true);
    expect(dtMap.splitMap).toBe('year');
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 13: the in-chart Category selector reflects the STARTED split column', async () => {
    const errBefore = errCount();

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

    const settle = await v.diffCanvasColors(page, 'Bar chart');
    expect(settle.deltaPx).toBeGreaterThanOrEqual(0); 
    expect(settle.deltaPx).toBeLessThan(300);
    await v.setViewerProps(page, 'Bar chart', [{set: {splitMap: 'month'}, wait: 900}]);
    const monthMap = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      return {splitMap: bc.props.splitMap, hasCanvas: !!bc.root.querySelector('canvas')};
    });
    await v.waitForViewerRendered(page, 'Bar chart', 300);

    const monthDelta = await v.diffCanvasColors(page, 'Bar chart');
    await v.setViewerProps(page, 'Bar chart', [{set: {splitMap: 'year'}, wait: 900}]);
    const yearMap = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      return {splitMap: bc.props.splitMap, hasCanvas: !!bc.root.querySelector('canvas'), splitCol: bc.props.splitColumnName};
    });

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
