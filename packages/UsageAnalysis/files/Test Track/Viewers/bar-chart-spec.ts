import {test, expect} from '@playwright/test';

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:DemoFiles/SPGI.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Bar chart tests', async ({page}) => {
  // Phase 1: Navigate
  await page.goto(baseUrl, {timeout: 60000, waitUntil: 'networkidle'});
  await page.waitForFunction(() => {
    try {
      return typeof grok !== 'undefined'
        && grok.shell
        && typeof grok.shell.closeAll === 'function'
        && grok.dapi
        && grok.dapi.files;
    } catch { return false; }
  }, {timeout: 120000});
  await page.waitForTimeout(3000);

  // Phase 2: Open dataset
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch {}
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    const hasBioChem = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
      .some(c => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 5000));
    }
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Phase 3: Add Bar chart
  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-bar-chart"]') as HTMLElement;
    icon.click();
  });
  await page.locator('[name="viewer-Bar-chart"]').waitFor({timeout: 5000});

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
    const result = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      const r: any[] = [];

      bc.props.stackColumnName = 'SEX';
      await new Promise(res => setTimeout(res, 500));
      r.push(bc.props.stackColumnName);

      bc.props.stackColumnName = 'RACE';
      await new Promise(res => setTimeout(res, 500));
      r.push(bc.props.stackColumnName);

      bc.props.relativeValues = true;
      await new Promise(res => setTimeout(res, 500));
      r.push(bc.props.relativeValues);

      bc.props.relativeValues = false;
      await new Promise(res => setTimeout(res, 500));
      r.push(bc.props.relativeValues);

      bc.props.stackColumnName = '';
      await new Promise(res => setTimeout(res, 500));
      r.push(bc.props.stackColumnName);

      return r;
    });
    expect(result[0]).toBe('SEX');
    expect(result[1]).toBe('RACE');
    expect(result[2]).toBe(true);
    expect(result[3]).toBe(false);
    expect(result[4]).toBe('');
  });

  // #### Sorting
  await softStep('Sorting', async () => {
    const result = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      const r: any[] = [];

      bc.props.barSortType = 'by value';
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.barSortType);

      bc.props.barSortOrder = 'asc';
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.barSortOrder);

      bc.props.barSortOrder = 'desc';
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.barSortOrder);

      bc.props.barSortType = 'by category';
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.barSortType);

      bc.props.barSortOrder = 'asc';
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.barSortOrder);

      bc.props.barSortOrder = 'desc';
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.barSortOrder);

      return r;
    });
    expect(result).toEqual(['by value', 'asc', 'desc', 'by category', 'asc', 'desc']);
  });

  // #### Value axis type
  await softStep('Value axis type', async () => {
    const result = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      bc.props.valueColumnName = 'HEIGHT';
      bc.props.splitColumnName = 'RACE';
      const r: any[] = [];

      bc.props.axisType = 'logarithmic';
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.axisType);

      bc.props.axisType = 'linear';
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.axisType);

      bc.props.valueMin = 100;
      bc.props.valueMax = 200;
      await new Promise(res => setTimeout(res, 300));
      r.push({min: bc.props.valueMin, max: bc.props.valueMax});

      bc.props.valueMin = null;
      bc.props.valueMax = null;
      await new Promise(res => setTimeout(res, 300));
      r.push({min: bc.props.valueMin, max: bc.props.valueMax});

      return r;
    });
    expect(result[0]).toBe('logarithmic');
    expect(result[1]).toBe('linear');
    expect(result[2]).toEqual({min: 100, max: 200});
    expect(result[3]).toEqual({min: null, max: null});
  });

  // #### Color coding
  await softStep('Color coding', async () => {
    const result = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      bc.props.valueColumnName = 'AGE';
      const r: any[] = [];

      bc.props.colorColumnName = 'HEIGHT';
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.colorColumnName);

      bc.props.colorAggrType = 'min';
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.colorAggrType);

      bc.props.colorAggrType = 'max';
      await new Promise(res => setTimeout(res, 200));
      r.push(bc.props.colorAggrType);

      bc.props.colorAggrType = 'med';
      await new Promise(res => setTimeout(res, 200));
      r.push(bc.props.colorAggrType);

      bc.props.invertColorScheme = true;
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.invertColorScheme);

      bc.props.colorColumnName = '';
      bc.props.invertColorScheme = false;
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.colorColumnName);

      return r;
    });
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
    const result = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      const r: any = {};

      bc.props.showValueSelector = false;
      bc.props.showCategorySelector = false;
      bc.props.showStackSelector = false;
      bc.props.showValueAxis = false;
      bc.props.showCategoryValues = false;
      await new Promise(res => setTimeout(res, 200));
      r.allOff = [bc.props.showValueSelector, bc.props.showCategorySelector,
        bc.props.showStackSelector, bc.props.showValueAxis, bc.props.showCategoryValues];

      bc.props.showValueSelector = true;
      bc.props.showCategorySelector = true;
      bc.props.showStackSelector = true;
      bc.props.showValueAxis = true;
      bc.props.showCategoryValues = true;
      await new Promise(res => setTimeout(res, 200));
      r.allOn = [bc.props.showValueSelector, bc.props.showCategorySelector,
        bc.props.showStackSelector, bc.props.showValueAxis, bc.props.showCategoryValues];

      return r;
    });
    expect(result.allOff).toEqual([false, false, false, false, false]);
    expect(result.allOn).toEqual([true, true, true, true, true]);
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
    const result = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      bc.props.splitColumnName = 'RACE';
      bc.props.valueColumnName = 'AGE';
      bc.props.valueAggrType = 'avg';
      const r: boolean[] = [];

      bc.props.showValuesInsteadOfCategories = true;
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.showValuesInsteadOfCategories);

      bc.props.showValuesInsteadOfCategories = false;
      await new Promise(res => setTimeout(res, 300));
      r.push(bc.props.showValuesInsteadOfCategories);

      return r;
    });
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
      bc.props.splitColumnName = 'SUBJ';
      await new Promise(res => setTimeout(res, 500));
      const split1 = bc.props.splitColumnName;

      bc.props.valueColumnName = 'AGE';
      await new Promise(res => setTimeout(res, 300));
      const value1 = bc.props.valueColumnName;

      bc.props.splitColumnName = 'RACE';
      return {split1, value1};
    });
    expect(result.split1).toBe('SUBJ');
    expect(result.value1).toBe('AGE');
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
