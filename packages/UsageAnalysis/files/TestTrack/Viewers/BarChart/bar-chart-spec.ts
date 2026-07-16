import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import * as v from '../helpers/viewers';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:DemoFiles/SPGI.csv';

test('Bar chart tests', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await v.addViewerByIcon(page, 'bar-chart', 'Bar-chart');

  // Open settings
  await page.evaluate(() => {
    const bcEl = document.querySelector('[name="viewer-Bar-chart"]') as HTMLElement;
    const panelBase = bcEl.closest('.panel-base') as HTMLElement;
    const gear = panelBase.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
    gear.click();
  });
  await page.waitForTimeout(500);

  // #### Stack column
  await softStep('Stack column', async () => {
    const result = await v.setViewerProps(page, 'Bar chart', [
      {set: {stackColumnName: 'SEX'}, wait: 500, read: 'stackColumnName'},
      {set: {stackColumnName: 'RACE'}, wait: 500, read: 'stackColumnName'},
      {set: {relativeValues: true}, wait: 500, read: 'relativeValues'},
      {set: {relativeValues: false}, wait: 500, read: 'relativeValues'},
      {set: {stackColumnName: ''}, wait: 500, read: 'stackColumnName'},
    ]);
    expect(result[0]).toBe('SEX');
    expect(result[1]).toBe('RACE');
    expect(result[2]).toBe(true);
    expect(result[3]).toBe(false);
    expect(result[4]).toBe('');
  });

  // #### Sorting
  await softStep('Sorting', async () => {
    const result = await v.setViewerProps(page, 'Bar chart', [
      {set: {splitColumnName: 'RACE', barSortType: 'by value'}, read: 'barSortType'},
      {set: {barSortOrder: 'asc'}, read: 'barSortOrder'},
      {set: {barSortOrder: 'desc'}, read: 'barSortOrder'},
      {set: {barSortType: 'by category'}, read: 'barSortType'},
      {set: {barSortOrder: 'asc'}, read: 'barSortOrder'},
      {set: {barSortOrder: 'desc'}, read: 'barSortOrder'},
    ]);
    expect(result).toEqual(['by value', 'asc', 'desc', 'by category', 'asc', 'desc']);
  });

  // #### Value axis type
  await softStep('Value axis type', async () => {
    const result = await v.setViewerProps(page, 'Bar chart', [
      {set: {valueColumnName: 'HEIGHT', splitColumnName: 'RACE', axisType: 'logarithmic'}, read: 'axisType'},
      {set: {axisType: 'linear'}, read: 'axisType'},
      {set: {valueMin: 100, valueMax: 200}, read: ['valueMin', 'valueMax']},
      {set: {valueMin: null, valueMax: null}, read: ['valueMin', 'valueMax']},
    ]);
    expect(result[0]).toBe('logarithmic');
    expect(result[1]).toBe('linear');
    expect(result[2]).toEqual({valueMin: 100, valueMax: 200});
    expect(result[3]).toEqual({valueMin: null, valueMax: null});
  });

  // #### Color coding
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

  // #### Include nulls
  await softStep('Include nulls', async () => {
    const result = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'DIS_POP';
      const r: any[] = [];

      r.push(bc.props.includeNulls); // default true

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

  // #### Bar style
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

      // Reset
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

  // #### Labels
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

  // #### Controls visibility
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

  // #### Aggregation types
  await softStep('Aggregation types', async () => {
    const result = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      bc.props.valueColumnName = 'AGE';
      const r: string[] = [];
      for (const aggr of ['avg', 'min', 'max', 'sum', 'count']) {
        bc.props.valueAggrType = aggr;
        await new Promise(res => setTimeout(res, 200));
        r.push(bc.props.valueAggrType);
      }
      bc.props.valueColumnName = 'WEIGHT';
      bc.props.valueAggrType = 'avg';
      await new Promise(res => setTimeout(res, 200));
      r.push(bc.props.valueColumnName + ':' + bc.props.valueAggrType);
      return r;
    });
    expect(result).toEqual(['avg', 'min', 'max', 'sum', 'count', 'WEIGHT:avg']);
  });

  // #### Date/time split column
  await softStep('Date/time split column', async () => {
    const result = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'STARTED';
      await new Promise(res => setTimeout(res, 500));
      const r: string[] = [];
      for (const map of ['year', 'month', 'quarter']) {
        bc.props.splitMap = map;
        await new Promise(res => setTimeout(res, 300));
        r.push(bc.props.splitMap);
      }
      return r;
    });
    expect(result).toEqual(['year', 'month', 'quarter']);
  });

  // #### Legend
  await softStep('Legend', async () => {
    const result = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      bc.props.stackColumnName = 'SEX';
      await new Promise(res => setTimeout(res, 300));
      const r: any[] = [];

      bc.props.legendVisibility = 'Always';
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.legendVisibility);

      for (const pos of ['Left', 'Right', 'Top', 'Bottom']) {
        bc.props.legendPosition = pos;
        await new Promise(res => setTimeout(res, 200));
        r.push(bc.props.legendPosition);
      }

      bc.props.legendVisibility = 'Never';
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.legendVisibility);

      bc.props.stackColumnName = '';
      await new Promise(res => setTimeout(res, 300));

      return r;
    });
    expect(result[0]).toBe('Always');
    expect(result.slice(1, 5)).toEqual(['Left', 'Right', 'Top', 'Bottom']);
    expect(result[5]).toBe('Never');
  });

  // #### Title and description
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

  // #### Show values instead of categories
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

  // #### Orientation
  await softStep('Orientation', async () => {
    const result = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      const r: string[] = [];
      for (const o of ['horizontal', 'vertical', 'auto']) {
        bc.props.orientation = o;
        await new Promise(res => setTimeout(res, 300));
        r.push(bc.props.orientation);
      }
      return r;
    });
    expect(result).toEqual(['horizontal', 'vertical', 'auto']);
  });

  // #### Data panel (SPGI dataset)
  await softStep('Data panel (SPGI dataset)', async () => {
    const result = await page.evaluate(async () => {
      // Setup: close all, open demog + SPGI
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

      // Switch to demog view
      const views = Array.from(grok.shell.views).filter((v: any) => v.type === 'TableView');
      const demogView = views.find((v: any) => v.dataFrame.name !== 'SPGI') as any;
      if (demogView) grok.shell.v = demogView;
      await new Promise(r => setTimeout(r, 500));

      // Add bar chart
      const icon = document.querySelector('[name="icon-bar-chart"]') as HTMLElement;
      icon.click();
      await new Promise(r => setTimeout(r, 1000));

      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      const r: any[] = [];

      // Row Source
      bc.props.rowSource = 'Filtered';
      await new Promise(res => setTimeout(res, 200));
      r.push(bc.props.rowSource);
      bc.props.rowSource = 'All';

      // Switch table to SPGI
      const spgi = Array.from(grok.shell.tables).find((t: any) => t.name === 'SPGI') as any;
      bc.dataFrame = spgi;
      await new Promise(res => setTimeout(res, 500));
      r.push(bc.dataFrame.name);

      // Set Filter
      bc.props.filter = '${CAST Idea ID} < 636500';
      await new Promise(res => setTimeout(res, 500));
      r.push(bc.props.filter);

      // Color Column
      bc.props.colorColumnName = 'Chemical Space Y';
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.colorColumnName);

      // Save layout
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(res => setTimeout(res, 1000));

      // Close viewer
      bc.close();
      await new Promise(res => setTimeout(res, 500));

      // Apply layout
      const saved = await grok.dapi.layouts.find(layoutId);
      grok.shell.tv.loadLayout(saved);
      await new Promise(res => setTimeout(res, 3000));

      const bc2 = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      r.push(bc2 ? bc2.props.colorColumnName : 'NOT_RESTORED');
      r.push(bc2 ? bc2.props.filter : 'NOT_RESTORED');

      // Cleanup
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

  // #### Filter Panel interaction
  await softStep('Filter Panel interaction', async () => {
    const result = await page.evaluate(async () => {
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));

      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      const tv = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });

      const icon = document.querySelector('[name="icon-bar-chart"]') as HTMLElement;
      icon.click();
      await new Promise(r => setTimeout(r, 1000));

      const bc = Array.from(tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      bc.props.valueColumnName = 'AGE';

      // Open filters
      tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 1000));
      const fg = tv.getFiltersGroup();

      const r: any[] = [];

      // Categorical filter
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: ['Asian', 'Caucasian']});
      await new Promise(res => setTimeout(res, 500));
      r.push(df.filter.trueCount);

      // Numeric filter
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 30, max: 60});
      await new Promise(res => setTimeout(res, 500));
      r.push(df.filter.trueCount);

      // Reset
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: df.col('RACE').categories});
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: df.col('AGE').min, max: df.col('AGE').max});
      await new Promise(res => setTimeout(res, 500));
      r.push(df.filter.trueCount);

      return {filtered1: r[0], filtered2: r[1], reset: r[2], total: df.rowCount};
    });
    expect(result.filtered1).toBeLessThan(result.total);
    expect(result.filtered2).toBeLessThan(result.filtered1);
    expect(result.reset).toBe(result.total);
  });

  // #### Scrolling with range slider
  await softStep('Scrolling with range slider', async () => {
    const result = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'USUBJID';
      await new Promise(res => setTimeout(res, 500));
      const split1 = bc.props.splitColumnName;

      bc.props.valueColumnName = 'AGE';
      await new Promise(res => setTimeout(res, 300));
      const value1 = bc.props.valueColumnName;

      bc.props.splitColumnName = 'RACE';
      return {split1, value1};
    });
    expect(result.split1).toBe('USUBJID');
    expect(result.value1).toBe('AGE');
  });

  v.finishSpec();
});
