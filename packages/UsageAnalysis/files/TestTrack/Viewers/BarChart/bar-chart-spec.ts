import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

test('Bar chart tests', async ({page}) => {
  test.setTimeout(600_000);

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

  await softStep('Color coding', async () => {
    const result = await v.setViewerProps(page, 'Bar chart', [
      {set: {splitColumnName: 'RACE', valueColumnName: 'AGE', colorColumnName: 'HEIGHT'}, read: 'colorColumnName'},
      {set: {colorAggrType: 'min'}, read: 'colorAggrType'},
      {set: {colorAggrType: 'max'}, wait: 200, read: 'colorAggrType'},
      {set: {colorAggrType: 'med'}, wait: 200, read: 'colorAggrType'},
      {set: {invertColorScheme: true}, read: 'invertColorScheme'},
      {set: {colorColumnName: '', invertColorScheme: false}, read: 'colorColumnName'},
    ]);
    expect(result[0]).toBe('HEIGHT');
    expect(result[1]).toBe('min');
    expect(result[2]).toBe('max');
    expect(result[3]).toBe('med');
    expect(result[4]).toBe(true);
    expect(result[5]).toBe('');
  });

  await softStep('Include nulls', async () => {
    const result = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'DIS_POP';
      const r: any[] = [];

      r.push(bc.props.includeNulls);

      bc.props.includeNulls = false;
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.includeNulls);

      bc.props.includeNulls = true;
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.includeNulls);

      return r;
    });
    expect(result).toEqual([true, false, true]);
  });

  await softStep('Bar style', async () => {
    const result = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      const r: any[] = [];

      bc.props.barBorderLineWidth = 2;
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.barBorderLineWidth);

      bc.props.barCornerRadius = 10;
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.barCornerRadius);

      bc.props.maxBarHeight = 20;
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.maxBarHeight);

      bc.props.verticalAlign = 'Top';
      await new Promise(res => setTimeout(res, 200));
      r.push(bc.props.verticalAlign);
      bc.props.verticalAlign = 'Bottom';
      await new Promise(res => setTimeout(res, 200));
      r.push(bc.props.verticalAlign);
      bc.props.verticalAlign = 'Center';
      await new Promise(res => setTimeout(res, 200));
      r.push(bc.props.verticalAlign);

      bc.props.showCategoryZeroBaseline = false;
      await new Promise(res => setTimeout(res, 200));
      r.push(bc.props.showCategoryZeroBaseline);

      bc.props.barBorderLineWidth = 0;
      bc.props.barCornerRadius = 0;
      bc.props.maxBarHeight = 50;
      bc.props.showCategoryZeroBaseline = true;

      return r;
    });
    expect(result[0]).toBe(2);
    expect(result[1]).toBe(10);
    expect(result[2]).toBe(20);
    expect(result[3]).toBe('Top');
    expect(result[4]).toBe('Bottom');
    expect(result[5]).toBe('Center');
    expect(result[6]).toBe(false);
  });

  await softStep('Labels', async () => {
    const result = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      bc.props.valueColumnName = 'AGE';
      const r: string[] = [];
      for (const val of ['inside', 'outside', 'never', 'auto']) {
        bc.props.showLabels = val;
        await new Promise(res => setTimeout(res, 200));
        r.push(bc.props.showLabels);
      }
      return r;
    });
    expect(result).toEqual(['inside', 'outside', 'never', 'auto']);
  });

  await softStep('Controls visibility', async () => {
    const ctrls = ['showValueSelector', 'showCategorySelector', 'showStackSelector',
      'showValueAxis', 'showCategoryValues'];
    const off = Object.fromEntries(ctrls.map((k) => [k, false]));
    const on = Object.fromEntries(ctrls.map((k) => [k, true]));
    const result = await v.setViewerProps(page, 'Bar chart', [
      {set: off, wait: 200, read: ctrls},
      {set: on, wait: 200, read: ctrls},
    ]);
    expect(result[0]).toEqual(off);
    expect(result[1]).toEqual(on);
  });

  await softStep('Aggregation types (max + value-column switch residual)', async () => {
    const result = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      bc.props.valueColumnName = 'AGE';
      bc.props.valueAggrType = 'max';
      await new Promise(res => setTimeout(res, 200));
      const aggr = bc.props.valueAggrType;
      bc.props.valueColumnName = 'WEIGHT';
      await new Promise(res => setTimeout(res, 200));
      return {aggr, value: bc.props.valueColumnName};
    });
    expect(result.aggr).toBe('max');
    expect(result.value).toBe('WEIGHT');
  });

  await softStep('Legend position (Left/Right/Top/Bottom residual)', async () => {
    const result = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      bc.props.stackColumnName = 'SEX';
      bc.props.legendVisibility = 'Always';
      await new Promise(res => setTimeout(res, 300));
      const r: string[] = [];
      for (const pos of ['Left', 'Right', 'Top', 'Bottom']) {
        bc.props.legendPosition = pos;
        await new Promise(res => setTimeout(res, 200));
        r.push(bc.props.legendPosition);
      }
      bc.props.stackColumnName = '';
      await new Promise(res => setTimeout(res, 300));
      return r;
    });
    expect(result).toEqual(['Left', 'Right', 'Top', 'Bottom']);
  });

  await softStep('Title and description', async () => {
    const result = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      const r: any[] = [];

      bc.props.showTitle = true;
      r.push(bc.props.showTitle);

      bc.props.title = 'Demographics';
      r.push(bc.props.title);

      bc.props.description = 'By race';
      r.push(bc.props.description);

      for (const pos of ['Top', 'Bottom', 'Left', 'Right']) {
        bc.props.descriptionPosition = pos;
        await new Promise(res => setTimeout(res, 200));
        r.push(bc.props.descriptionPosition);
      }

      bc.props.descriptionVisibilityMode = 'Never';
      r.push(bc.props.descriptionVisibilityMode);

      bc.props.showTitle = false;
      return r;
    });
    expect(result[0]).toBe(true);
    expect(result[1]).toBe('Demographics');
    expect(result[2]).toBe('By race');
    expect(result.slice(3, 7)).toEqual(['Top', 'Bottom', 'Left', 'Right']);
    expect(result[7]).toBe('Never');
  });

  await softStep('Show values instead of categories', async () => {
    const result = await v.setViewerProps(page, 'Bar chart', [
      {
        set: {splitColumnName: 'RACE', valueColumnName: 'AGE', valueAggrType: 'avg',
          showValuesInsteadOfCategories: true},
        read: 'showValuesInsteadOfCategories',
      },
      {set: {showValuesInsteadOfCategories: false}, read: 'showValuesInsteadOfCategories'},
    ]);
    expect(result).toEqual([true, false]);
  });

  await softStep('Data panel (demog → SPGI table switch)', async () => {
    const result = await page.evaluate(async () => {
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));

      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      const tv = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });

      const df2 = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      df2.name = 'SPGI';
      grok.shell.addTableView(df2);
      await new Promise(resolve => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });

      const views = Array.from(grok.shell.views).filter((v: any) => v.type === 'TableView');
      const demogView = views.find((v: any) => v.dataFrame.name !== 'SPGI') as any;
      if (demogView) grok.shell.v = demogView;
      await new Promise(r => setTimeout(r, 500));

      const icon = document.querySelector('[name="icon-bar-chart"]') as HTMLElement;
      icon.click();
      await new Promise(r => setTimeout(r, 1000));

      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      const r: any[] = [];

      bc.props.rowSource = 'Filtered';
      await new Promise(res => setTimeout(res, 200));
      r.push(bc.props.rowSource);
      bc.props.rowSource = 'All';

      const spgi = Array.from(grok.shell.tables).find((t: any) => t.name === 'SPGI') as any;
      bc.dataFrame = spgi;
      await new Promise(res => setTimeout(res, 500));
      r.push(bc.dataFrame.name);

      bc.props.filter = '${CAST Idea ID} < 636500';
      await new Promise(res => setTimeout(res, 500));
      r.push(bc.props.filter);

      bc.props.colorColumnName = 'Chemical Space Y';
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.colorColumnName);

      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(res => setTimeout(res, 1000));

      bc.close();
      await new Promise(res => setTimeout(res, 500));

      const saved = await grok.dapi.layouts.find(layoutId);
      grok.shell.tv.loadLayout(saved);
      await new Promise(res => setTimeout(res, 3000));

      const bc2 = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      r.push(bc2 ? bc2.props.colorColumnName : 'NOT_RESTORED');
      r.push(bc2 ? bc2.props.filter : 'NOT_RESTORED');

      await grok.dapi.layouts.delete(saved);

      return r;
    });
    expect(result[0]).toBe('Filtered');
    expect(result[1]).toBe('SPGI');
    expect(result[2]).toBe('${CAST Idea ID} < 636500');
    expect(result[3]).toBe('Chemical Space Y');
    expect(result[4]).toBe('Chemical Space Y');
    expect(result[5]).toBe('${CAST Idea ID} < 636500');
  });

  v.finishSpec();
});
