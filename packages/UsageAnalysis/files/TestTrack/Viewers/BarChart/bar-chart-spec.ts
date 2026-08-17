/* ---
realizes: []
--- */

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const T = {
  colorColumn: 800,
  invertScheme: 800,
  includeNulls: 100000,
  barBorder: 2000,
  maxBarHeight: 400,
  labels: 1000,
  aggrType: 400,
  valueSwitch: 400,
  showValues: 800,
};
const PRECHECK_CEIL = 250;

test('Bar chart tests', async ({page}) => {
  test.setTimeout(600_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await v.addViewerByIcon(page, 'bar-chart', 'Bar-chart');

  await page.evaluate(() => {
    const bcEl = document.querySelector('[name="viewer-Bar-chart"]') as HTMLElement;
    const panelBase = bcEl.closest('.panel-base') as HTMLElement;
    const gear = panelBase.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
    gear.click();
  });
  await v.pollValue(() => page.locator('.property-grid').count(), (n) => n > 0, 500, 100);

  async function canvasBaseline(): Promise<number> {
    await v.snapshotCanvasColors(page, 'Bar chart');
    await page.waitForTimeout(500);
    return (await v.diffCanvasColors(page, 'Bar chart')).deltaPx;
  }
  async function canvasDelta(): Promise<number> {
    return (await v.diffCanvasColors(page, 'Bar chart')).deltaPx;
  }

  await softStep('Color coding', async () => {
    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      bc.props.valueColumnName = 'AGE';
      bc.props.colorColumnName = '';
      bc.props.invertColorScheme = false;
    });
    await v.waitForViewerRendered(page, 'Bar chart', 600);

    const preColor = await canvasBaseline();
    const setColor = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.colorColumnName = 'HEIGHT';
      return {colorColumnName: bc.props.colorColumnName, aggr: bc.props.colorAggrType};
    });
    await v.waitForViewerRendered(page, 'Bar chart', 600);
    const colorDelta = await canvasDelta();
    console.log(`[bar-chart] colorDelta=${colorDelta}`);

    const aggrReads = await v.setViewerProps(page, 'Bar chart', [
      {set: {colorAggrType: 'min'}, wait: 200, read: 'colorAggrType'},
      {set: {colorAggrType: 'max'}, wait: 200, read: 'colorAggrType'},
      {set: {colorAggrType: 'med'}, wait: 200, read: 'colorAggrType'},
    ]);

    await page.waitForTimeout(400);

    const preInvert = await canvasBaseline();
    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.invertColorScheme = true;
    });
    await v.waitForViewerRendered(page, 'Bar chart', 600);
    const invertDelta = await canvasDelta();
    console.log(`[bar-chart] invertDelta=${invertDelta}`);

    const colColRead = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.colorColumnName = '';
      bc.props.invertColorScheme = false;
      return bc.props.colorColumnName;
    });
    await v.waitForViewerRendered(page, 'Bar chart', 300);

    expect(setColor.colorColumnName).toBe('HEIGHT');
    expect(aggrReads).toEqual(['min', 'max', 'med']);
    expect(preColor).toBeGreaterThanOrEqual(0);
    expect(preColor).toBeLessThan(PRECHECK_CEIL);
    expect(colorDelta).toBeGreaterThan(T.colorColumn);
    expect(preInvert).toBeGreaterThanOrEqual(0);
    expect(preInvert).toBeLessThan(PRECHECK_CEIL);
    expect(invertDelta).toBeGreaterThan(T.invertScheme);
    expect(colColRead).toBe('');
  });

  await softStep('Include nulls', async () => {
    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'HEIGHT';
      bc.props.includeNulls = true;
    });
    await v.waitForViewerRendered(page, 'Bar chart', 600);

    const preOff = await canvasBaseline();
    const offRead = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.includeNulls = false;
      return bc.props.includeNulls;
    });
    await v.waitForViewerRendered(page, 'Bar chart', 500);
    const offDelta = await canvasDelta();
    console.log(`[bar-chart] includeNulls offDelta=${offDelta}`);

    const preOn = await canvasBaseline();
    const onRead = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.includeNulls = true;
      return bc.props.includeNulls;
    });
    await v.waitForViewerRendered(page, 'Bar chart', 500);
    const onDelta = await canvasDelta();
    console.log(`[bar-chart] includeNulls onDelta=${onDelta}`);

    expect(offRead).toBe(false);
    expect(onRead).toBe(true);
    expect(preOff).toBeGreaterThanOrEqual(0);
    expect(preOff).toBeLessThan(PRECHECK_CEIL);
    expect(offDelta).toBeGreaterThan(T.includeNulls);
    expect(preOn).toBeGreaterThanOrEqual(0);
    expect(preOn).toBeLessThan(PRECHECK_CEIL);
    expect(onDelta).toBeGreaterThan(T.includeNulls);
  });

  await softStep('Bar style', async () => {
    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      bc.props.valueColumnName = 'AGE';
    });
    await v.waitForViewerRendered(page, 'Bar chart', 600);

    const preBorder = await canvasBaseline();
    const borderRead = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.barBorderLineWidth = 2;
      return bc.props.barBorderLineWidth;
    });
    await v.waitForViewerRendered(page, 'Bar chart', 500);
    const borderDelta = await canvasDelta();
    console.log(`[bar-chart] borderDelta=${borderDelta}`);

    const preHeight = await canvasBaseline();
    const heightRead = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.maxBarHeight = 20;
      return bc.props.maxBarHeight;
    });
    await v.waitForViewerRendered(page, 'Bar chart', 500);
    const heightDelta = await canvasDelta();
    console.log(`[bar-chart] heightDelta=${heightDelta}`);

    const styleReads = await v.setViewerProps(page, 'Bar chart', [
      {set: {barCornerRadius: 10}, wait: 200, read: 'barCornerRadius'},
      {set: {verticalAlign: 'Top'}, wait: 200, read: 'verticalAlign'},
      {set: {verticalAlign: 'Bottom'}, wait: 200, read: 'verticalAlign'},
      {set: {verticalAlign: 'Center'}, wait: 200, read: 'verticalAlign'},
      {set: {showCategoryZeroBaseline: false}, wait: 200, read: 'showCategoryZeroBaseline'},
    ]);

    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.barBorderLineWidth = 0;
      bc.props.barCornerRadius = 0;
      bc.props.maxBarHeight = 50;
      bc.props.showCategoryZeroBaseline = true;
    });
    await v.waitForViewerRendered(page, 'Bar chart', 300);

    expect(borderRead).toBe(2);
    expect(heightRead).toBe(20);
    expect(styleReads).toEqual([10, 'Top', 'Bottom', 'Center', false]);
    expect(preBorder).toBeGreaterThanOrEqual(0);
    expect(preBorder).toBeLessThan(PRECHECK_CEIL);
    expect(borderDelta).toBeGreaterThan(T.barBorder);
    expect(preHeight).toBeGreaterThanOrEqual(0);
    expect(preHeight).toBeLessThan(PRECHECK_CEIL);
    expect(heightDelta).toBeGreaterThan(T.maxBarHeight);
  });

  await softStep('Labels', async () => {
    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      bc.props.valueColumnName = 'AGE';
      bc.props.showLabels = 'inside';
    });
    await v.waitForViewerRendered(page, 'Bar chart', 600);

    const preLabels = await canvasBaseline();
    const neverRead = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.showLabels = 'never';
      return bc.props.showLabels;
    });
    await v.waitForViewerRendered(page, 'Bar chart', 500);
    const labelsDelta = await canvasDelta();
    console.log(`[bar-chart] labelsDelta=${labelsDelta}`);

    const labelReads = await v.setViewerProps(page, 'Bar chart', [
      {set: {showLabels: 'outside'}, wait: 200, read: 'showLabels'},
      {set: {showLabels: 'auto'}, wait: 200, read: 'showLabels'},
    ]);

    expect(neverRead).toBe('never');
    expect(labelReads).toEqual(['outside', 'auto']);
    expect(preLabels).toBeGreaterThanOrEqual(0);
    expect(preLabels).toBeLessThan(PRECHECK_CEIL);
    expect(labelsDelta).toBeGreaterThan(T.labels);
  });

  await softStep('Controls visibility', async () => {
    const ctrls = ['showValueSelector', 'showCategorySelector', 'showStackSelector',
      'showValueAxis', 'showCategoryValues'];
    const off = Object.fromEntries(ctrls.map((k) => [k, false]));
    const on = Object.fromEntries(ctrls.map((k) => [k, true]));

    const countVisibleSelectors = () => page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const root = bc.root as HTMLElement;
      return Array.from(root.querySelectorAll('.d4-column-selector'))
        .filter((e: any) => getComputedStyle(e).display !== 'none').length;
    });

    const offReads = await v.setViewerProps(page, 'Bar chart', [{set: off, wait: 400, read: ctrls}]);
    const offVisible = await countVisibleSelectors();
    const onReads = await v.setViewerProps(page, 'Bar chart', [{set: on, wait: 400, read: ctrls}]);
    const onVisible = await countVisibleSelectors();

    expect(offReads[0]).toEqual(off);
    expect(onReads[0]).toEqual(on);
    expect(offVisible).toBe(0);
    expect(onVisible).toBe(3);
  });

  await softStep('Aggregation types', async () => {
    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      bc.props.valueColumnName = 'AGE';
      bc.props.valueAggrType = 'avg';
    });
    await v.waitForViewerRendered(page, 'Bar chart', 600);

    const preAggr = await canvasBaseline();
    const aggrRead = (await v.setViewerProps(page, 'Bar chart',
      [{set: {valueAggrType: 'max'}, wait: 500, read: 'valueAggrType'}]))[0];
    const aggrDelta = await canvasDelta();
    console.log(`[bar-chart] aggrDelta=${aggrDelta}`);

    const preSwitch = await canvasBaseline();
    const valueRead = (await v.setViewerProps(page, 'Bar chart',
      [{set: {valueColumnName: 'WEIGHT'}, wait: 500, read: 'valueColumnName'}]))[0];
    const switchDelta = await canvasDelta();
    console.log(`[bar-chart] switchDelta=${switchDelta}`);

    expect(aggrRead).toBe('max');
    expect(valueRead).toBe('WEIGHT');
    expect(preAggr).toBeGreaterThanOrEqual(0);
    expect(preAggr).toBeLessThan(PRECHECK_CEIL);
    expect(aggrDelta).toBeGreaterThan(T.aggrType);
    expect(preSwitch).toBeGreaterThanOrEqual(0);
    expect(preSwitch).toBeLessThan(PRECHECK_CEIL);
    expect(switchDelta).toBeGreaterThan(T.valueSwitch);
  });

  await softStep('Legend position', async () => {
    await page.evaluate(async () => {
      const old = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      if (old) old.close();
      await new Promise((r) => setTimeout(r, 400));
    });
    await v.addViewerByIcon(page, 'bar-chart', 'Bar-chart');
    await page.waitForTimeout(500);

    const legend = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const root = bc.root as HTMLElement;
      const laidOut = () => {
        const el = root.querySelector('[name="legend"]') as HTMLElement | null;
        return !!el && getComputedStyle(el).display !== 'none' && el.getBoundingClientRect().height > 0;
      };
      const items = () => root.querySelectorAll('[name="legend"] .d4-legend-item').length;

      bc.props.splitColumnName = 'RACE';
      bc.props.valueColumnName = 'AGE';
      bc.props.stackColumnName = '';
      bc.props.legendVisibility = 'Always';
      await new Promise((r) => setTimeout(r, 800));
      const before = laidOut();

      bc.props.stackColumnName = 'SEX';
      await new Promise((r) => setTimeout(r, 1500));
      const after = laidOut();
      const afterItems = items();

      const positions: string[] = [];
      for (const pos of ['Left', 'Right', 'Top', 'Bottom']) {
        bc.props.legendPosition = pos;
        await new Promise((r) => setTimeout(r, 200));
        positions.push(bc.props.legendPosition);
      }

      bc.props.stackColumnName = '';
      await new Promise((r) => setTimeout(r, 1500));
      const clearedLaidOut = laidOut();
      const clearedItems = items();

      return {before, after, afterItems, positions, clearedLaidOut, clearedItems};
    });

    expect(legend.before).toBe(false);
    expect(legend.after).toBe(true);
    expect(legend.afterItems).toBeGreaterThanOrEqual(2);
    expect(legend.positions).toEqual(['Left', 'Right', 'Top', 'Bottom']);
    expect(legend.clearedLaidOut).toBe(false);
    expect(legend.clearedItems).toBe(0);
  });

  await softStep('Title and description', async () => {
    const info = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.showTitle = true;
      bc.props.title = 'Demographics';
      bc.props.description = 'By race';
      await new Promise((r) => setTimeout(r, 400));
      const root = bc.root as HTMLElement;

      const panel = (root.closest('.panel-base') as HTMLElement) ?? root;
      const panelText = panel.innerText ?? panel.textContent ?? '';
      const rootText = root.innerText ?? root.textContent ?? '';
      const r: any = {showTitle: bc.props.showTitle, titleInDom: panelText.includes('Demographics'),
        descInDom: rootText.includes('By race')};

      const positions: string[] = [];
      for (const pos of ['Top', 'Bottom', 'Left', 'Right']) {
        bc.props.descriptionPosition = pos;
        await new Promise((res) => setTimeout(res, 150));
        positions.push(bc.props.descriptionPosition);
      }
      r.positions = positions;

      bc.props.descriptionVisibilityMode = 'Never';
      await new Promise((res) => setTimeout(res, 300));
      r.hidden = bc.props.descriptionVisibilityMode;

      const rootTextHidden = root.innerText ?? root.textContent ?? '';
      r.descGoneFromDom = rootTextHidden.includes('By race');

      bc.props.showTitle = false;
      return r;
    });

    expect(info.showTitle).toBe(true);
    expect(info.titleInDom).toBe(true);
    expect(info.descInDom).toBe(true);
    expect(info.positions).toEqual(['Top', 'Bottom', 'Left', 'Right']);
    expect(info.hidden).toBe('Never');
    expect(info.descGoneFromDom).toBe(false);
  });

  await softStep('Show values instead of categories', async () => {
    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      bc.props.valueColumnName = 'AGE';
      bc.props.valueAggrType = 'avg';
      bc.props.showValuesInsteadOfCategories = false;
    });
    await v.waitForViewerRendered(page, 'Bar chart', 600);

    const preShow = await canvasBaseline();
    const onRead = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.showValuesInsteadOfCategories = true;
      return bc.props.showValuesInsteadOfCategories;
    });
    await v.waitForViewerRendered(page, 'Bar chart', 500);
    const showDelta = await canvasDelta();
    console.log(`[bar-chart] showValuesDelta=${showDelta}`);

    const offRead = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.showValuesInsteadOfCategories = false;
      return bc.props.showValuesInsteadOfCategories;
    });
    await v.waitForViewerRendered(page, 'Bar chart', 300);

    expect(onRead).toBe(true);
    expect(offRead).toBe(false);
    expect(preShow).toBeGreaterThanOrEqual(0);
    expect(preShow).toBeLessThan(PRECHECK_CEIL);
    expect(showDelta).toBeGreaterThan(T.showValues);
  });

  await softStep('Context menu', async () => {
    const result = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((view: any) => view.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      bc.props.valueColumnName = 'AGE';
      bc.props.stackColumnName = 'SEX';
      await new Promise((r) => setTimeout(r, 600));
      const canvas = bc.root.querySelector('canvas')!;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2,
      }));
      await new Promise((r) => setTimeout(r, 600));
      const items = Array.from(document.querySelectorAll('.d4-menu-item-label')).map((e) => e.textContent!.trim());
      const before = bc.props.showValueAxis;
      const sva = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find((e) => e.textContent!.trim() === 'Show Value Axis');
      if (sva) (sva.closest('.d4-menu-item') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 500));
      const after = bc.props.showValueAxis;
      bc.props.showValueAxis = before;
      return {items, before, after};
    });
    for (const label of ['Reset View', 'Orientation', 'On Click', 'Order', 'Controls',
      'Selection', 'Show Value Axis', 'Show Category Values', 'Show Selected Rows',
      'Include Nulls', 'Axis Type', 'Legend Visibility', 'Legend Position'])
      expect(result.items).toContain(label);
    expect(result.after).toBe(!result.before);
    await page.keyboard.press('Escape');
  });

  await softStep('Data panel', async () => {
    const result = await page.evaluate(async () => {
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 500));

      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      grok.shell.addTableView(df);
      await new Promise((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });

      const df2 = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      df2.name = 'SPGI';
      grok.shell.addTableView(df2);
      await new Promise((resolve) => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });

      const views = Array.from(grok.shell.views).filter((view: any) => view.type === 'TableView');
      const demogView = views.find((view: any) => view.dataFrame.name !== 'SPGI') as any;
      if (demogView) grok.shell.v = demogView;
      await new Promise((r) => setTimeout(r, 500));

      const icon = document.querySelector('[name="icon-bar-chart"]') as HTMLElement;
      icon.click();
      await new Promise((r) => setTimeout(r, 1000));

      const bc = Array.from(grok.shell.tv.viewers).find((view: any) => view.type === 'Bar chart') as any;
      const r: any[] = [];

      for (const src of ['Filtered', 'Selected', 'All']) {
        bc.props.rowSource = src;
        await new Promise((res) => setTimeout(res, 200));
        r.push(bc.props.rowSource);
      }
      bc.props.rowSource = 'All';

      const spgi = Array.from(grok.shell.tables).find((t: any) => t.name === 'SPGI') as any;
      bc.dataFrame = spgi;
      await new Promise((res) => setTimeout(res, 500));
      r.push(bc.dataFrame.name);

      bc.props.filter = '${CAST Idea ID} < 634835';
      await new Promise((res) => setTimeout(res, 500));
      r.push(bc.props.filter);

      bc.props.colorColumnName = 'Chemical Space Y';
      await new Promise((res) => setTimeout(res, 300));
      r.push(bc.props.colorColumnName);

      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise((res) => setTimeout(res, 1000));

      bc.close();
      await new Promise((res) => setTimeout(res, 500));

      const saved = await grok.dapi.layouts.find(layoutId);
      grok.shell.tv.loadLayout(saved);
      await new Promise((res) => setTimeout(res, 3000));

      const bc2 = Array.from(grok.shell.tv.viewers).find((view: any) => view.type === 'Bar chart') as any;
      r.push(bc2 ? bc2.props.colorColumnName : 'NOT_RESTORED');
      r.push(bc2 ? bc2.props.filter : 'NOT_RESTORED');

      await grok.dapi.layouts.delete(saved);

      return r;
    });
    expect(result.slice(0, 3)).toEqual(['Filtered', 'Selected', 'All']);
    expect(result[3]).toBe('SPGI');
    expect(result[4]).toBe('${CAST Idea ID} < 634835');
    expect(result[5]).toBe('Chemical Space Y');
    expect(result[6]).toBe('Chemical Space Y');
    expect(result[7]).toBe('${CAST Idea ID} < 634835');
  });

  await softStep('No page errors', async () => {
    expect(pageErrors).toEqual([]);
  });

  v.finishSpec();
});
