/* ---
realizes: [barchart.int.relative-requires-stack]
--- */

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';
const splitCol = 'Primary Series Name';
const valueCol = 'CAST Idea ID';
const stackCol = 'Stereo Category';

test('Bar Chart — Relative Values requires Stack column', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 4000});
  await v.addViewerByIcon(page, 'bar-chart', 'Bar-chart');

  await page.locator('[name="viewer-Bar-chart"]').first().hover();
  await page.waitForTimeout(300);

  await page.evaluate(() => {
    const bcEl = document.querySelector('[name="viewer-Bar-chart"]') as HTMLElement;
    const panelBase = bcEl.closest('.panel-base') as HTMLElement;
    const gear = panelBase.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
    gear?.click();
  });
  await page.waitForTimeout(500);

  await page.evaluate(async ({split, value}) => {
    const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
    bc.props.splitColumnName = split;
    bc.props.valueColumnName = value;
    bc.props.valueAggrType = 'count';
    bc.props.stackColumnName = null;
    bc.props.relativeValues = false;
    await new Promise((r) => setTimeout(r, 1000));
  }, {split: splitCol, value: valueCol});

  await softStep('Scenario 1 Step 3: Relative Values without a Stack column has no stacking effect, no error', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const info = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const legendAtBaseline = bc.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
      const stackAtEntry = bc.props.stackColumnName;
      bc.props.relativeValues = true;
      await new Promise((r) => setTimeout(r, 900));
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      const rect = cv ? cv.getBoundingClientRect() : {width: 0, height: 0};
      return {
        stackAtEntry,
        legendAtBaseline,
        rel: bc.props.relativeValues,
        stackAfter: bc.props.stackColumnName,
        hasCanvas: !!cv,
        canvasW: rect.width,
        canvasH: rect.height,
        legendRendered: !!bc.root.querySelector('[name="legend"]'),
        legendItems: bc.root.querySelectorAll('[name="legend"] .d4-legend-item').length,
      };
    });
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.stackAtEntry).toBeNull();
    expect(info.legendAtBaseline).toBe(0);
    expect(info.rel).toBe(true);
    expect(info.stackAfter).toBeNull();
    expect(info.hasCanvas).toBe(true);
    expect(info.canvasW).toBeGreaterThan(0);
    expect(info.canvasH).toBeGreaterThan(0);
    expect(info.legendRendered).toBe(false);
    expect(info.legendItems).toBe(0);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 5-6: setting a Stack column activates stacking (legend renders)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const info = await page.evaluate(async ({stack}) => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.stackColumnName = stack;
      await new Promise((r) => setTimeout(r, 1000));
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      const rect = cv ? cv.getBoundingClientRect() : {width: 0, height: 0};
      return {
        stack: bc.props.stackColumnName,
        rel: bc.props.relativeValues,
        hasCanvas: !!cv,
        canvasW: rect.width,
        legendRendered: !!bc.root.querySelector('[name="legend"]'),
        legendItems: bc.root.querySelectorAll('[name="legend"] .d4-legend-item').length,
        stackCats: grok.shell.tv.dataFrame.col(stack).categories.length,
      };
    }, {stack: stackCol});
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.stack).toBe(stackCol);
    expect(info.rel).toBe(true);
    expect(info.hasCanvas).toBe(true);
    expect(info.canvasW).toBeGreaterThan(0);
    expect(info.legendRendered).toBe(true);
    expect(info.legendItems).toBeGreaterThanOrEqual(2);
    expect(info.stackCats).toBeGreaterThanOrEqual(2);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 7-8: remove Stack reverts to inert; re-add re-activates (repeatable)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const info = await page.evaluate(async ({stack}) => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.stackColumnName = null;
      await new Promise((r) => setTimeout(r, 900));
      const afterRemove = {
        stack: bc.props.stackColumnName,
        legendRendered: !!bc.root.querySelector('[name="legend"]'),
        legendItems: bc.root.querySelectorAll('[name="legend"] .d4-legend-item').length,
        hasCanvas: !!bc.root.querySelector('canvas'),
        canvasW: (bc.root.querySelector('canvas') as HTMLCanvasElement)?.getBoundingClientRect().width ?? 0,
      };
      bc.props.stackColumnName = stack;
      await new Promise((r) => setTimeout(r, 900));
      const afterReadd = {
        stack: bc.props.stackColumnName,
        legendRendered: !!bc.root.querySelector('[name="legend"]'),
        legendItems: bc.root.querySelectorAll('[name="legend"] .d4-legend-item').length,
      };
      return {afterRemove, afterReadd, rel: bc.props.relativeValues};
    }, {stack: stackCol});
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.afterRemove.stack).toBeNull();
    expect(info.afterRemove.legendRendered).toBe(false);
    expect(info.afterRemove.legendItems).toBe(0);
    expect(info.afterRemove.hasCanvas).toBe(true);
    expect(info.afterRemove.canvasW).toBeGreaterThan(0);
    expect(info.afterReadd.stack).toBe(stackCol);
    expect(info.afterReadd.legendRendered).toBe(true);
    expect(info.afterReadd.legendItems).toBeGreaterThanOrEqual(2);
    expect(info.rel).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 3 Step 9-10: disabling Relative Values (no Stack) restores the baseline', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const info = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.stackColumnName = null;
      await new Promise((r) => setTimeout(r, 700));
      bc.props.relativeValues = false;
      await new Promise((r) => setTimeout(r, 700));
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      const rect = cv ? cv.getBoundingClientRect() : {width: 0, height: 0};
      return {
        rel: bc.props.relativeValues,
        stack: bc.props.stackColumnName,
        hasCanvas: !!cv,
        canvasW: rect.width,
        canvasH: rect.height,
        legendRendered: !!bc.root.querySelector('[name="legend"]'),
        legendItems: bc.root.querySelectorAll('[name="legend"] .d4-legend-item').length,
      };
    });
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.rel).toBe(false);
    expect(info.stack).toBeNull();
    expect(info.hasCanvas).toBe(true);
    expect(info.canvasW).toBeGreaterThan(0);
    expect(info.canvasH).toBeGreaterThan(0);
    expect(info.legendRendered).toBe(false);
    expect(info.legendItems).toBe(0);
    expect(errAfter).toBe(errBefore);
  });

  v.finishSpec();
});
