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

  await softStep('Stack + avg aggregation suppresses the legend column', async () => {
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
  });

  await softStep('Stack + sum aggregation produces a visible legend', async () => {
    await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.valueAggrType = 'sum';
      await new Promise((r) => setTimeout(r, 700));
    });
    const legendSum = await v.readLegend(page, 'Bar chart');
    expect(legendSum.legendRendered).toBe(true);
    expect(legendSum.itemCount).toBeGreaterThan(0);
  });

  await softStep('count keeps the legend; removing Stack collapses it', async () => {
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
  });

  await softStep('String Split column hides the Split Map selector', async () => {
    const strMap = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.valueColumnName = 'AGE';
      bc.props.valueAggrType = 'count';
      bc.props.splitColumnName = 'RACE';
      grok.shell.o = bc;
      await new Promise((r) => setTimeout(r, 1200));
      const host = document.querySelector('[name="input-aggr-selector-split-map"]') as HTMLElement | null;
      if (!host) return {present: false, visible: false};
      const rect = host.getBoundingClientRect();
      const display = getComputedStyle(host).display;
      return {
        present: true,
        visible: rect.width > 0 && rect.height > 0 && display !== 'none' && host.offsetParent !== null,
        splitType: bc.dataFrame.col('RACE').type,
      };
    });
    expect(strMap.splitType).toBe('string');
    expect(strMap.visible).toBe(false);
  });

  await softStep('DateTime Split column enables the Split Map selector', async () => {
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
  });

  await softStep('Split Map re-categorizes Year → Month → Year without breaking', async () => {
    const roundTrip = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const out: string[] = [];

      bc.props.splitMap = 'month';
      await new Promise((r) => setTimeout(r, 700));
      out.push(bc.props.splitMap);

      bc.props.splitMap = 'year';
      await new Promise((r) => setTimeout(r, 700));
      out.push(bc.props.splitMap);

      const hasCanvas = !!bc.root.querySelector('canvas');
      return {maps: out, hasCanvas, splitCol: bc.props.splitColumnName};
    });
    expect(roundTrip.maps).toEqual(['month', 'year']);
    expect(roundTrip.hasCanvas).toBe(true);
    expect(roundTrip.splitCol).toBe('STARTED');
  });

  v.finishSpec();
});
