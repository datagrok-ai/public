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

test('Pie chart tests', async ({page}) => {
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

  // Phase 3: Add Pie chart
  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-pie-chart"]') as HTMLElement;
    icon.click();
  });
  await page.locator('[name="viewer-Pie-chart"]').waitFor({timeout: 5000});

  // #### Sorting
  await softStep('Sorting', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      pie.props.categoryColumnName = 'RACE';
      const r: any[] = [];

      pie.props.pieSortType = 'by value';
      await new Promise(res => setTimeout(res, 300));
      r.push(pie.props.pieSortType);

      pie.props.pieSortOrder = 'desc';
      await new Promise(res => setTimeout(res, 300));
      r.push(pie.props.pieSortOrder);

      pie.props.pieSortOrder = 'asc';
      await new Promise(res => setTimeout(res, 300));
      r.push(pie.props.pieSortOrder);

      pie.props.pieSortType = 'by category';
      await new Promise(res => setTimeout(res, 300));
      r.push(pie.props.pieSortType);

      pie.props.pieSortOrder = 'asc';
      await new Promise(res => setTimeout(res, 300));
      r.push(pie.props.pieSortOrder);

      pie.props.pieSortOrder = 'desc';
      await new Promise(res => setTimeout(res, 300));
      r.push(pie.props.pieSortOrder);

      return r;
    });
    expect(result).toEqual(['by value', 'desc', 'asc', 'by category', 'asc', 'desc']);
  });

  // #### Segment angle and length
  await softStep('Segment angle and length', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      pie.props.categoryColumnName = 'RACE';
      const r: any[] = [];

      pie.props.segmentAngleColumnName = 'AGE';
      await new Promise(res => setTimeout(res, 300));
      r.push(pie.props.segmentAngleColumnName);

      pie.props.segmentAngleAggrType = 'sum';
      await new Promise(res => setTimeout(res, 300));
      r.push(pie.props.segmentAngleAggrType);

      pie.props.segmentAngleAggrType = 'count';
      await new Promise(res => setTimeout(res, 300));
      r.push(pie.props.segmentAngleAggrType);

      pie.props.segmentLengthColumnName = 'WEIGHT';
      await new Promise(res => setTimeout(res, 300));
      r.push(pie.props.segmentLengthColumnName);

      pie.props.segmentLengthAggrType = 'max';
      await new Promise(res => setTimeout(res, 300));
      r.push(pie.props.segmentLengthAggrType);

      pie.props.segmentAngleColumnName = '';
      pie.props.segmentLengthColumnName = '';
      await new Promise(res => setTimeout(res, 300));

      return r;
    });
    expect(result[0]).toBe('AGE');
    expect(result[1]).toBe('sum');
    expect(result[2]).toBe('count');
    expect(result[3]).toBe('WEIGHT');
    expect(result[4]).toBe('max');
  });

  // #### Appearance
  await softStep('Appearance', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const r: any[] = [];

      pie.props.startAngle = 90;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.startAngle);

      pie.props.startAngle = 180;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.startAngle);

      pie.props.startAngle = 0;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.startAngle);

      pie.props.maxRadius = 100;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.maxRadius);

      pie.props.maxRadius = 150;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.maxRadius);

      pie.props.shift = 10;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.shift);

      pie.props.shift = 0;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.shift);

      return r;
    });
    expect(result).toEqual([90, 180, 0, 100, 150, 10, 0]);
  });

  // #### Labels
  await softStep('Labels', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const r: any[] = [];

      for (const pos of ['Inside', 'Outside', 'Auto']) {
        pie.props.labelPosition = pos;
        await new Promise(res => setTimeout(res, 200));
        r.push(pie.props.labelPosition);
      }

      pie.props.showLabel = false;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.showLabel);

      pie.props.showPercentage = false;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.showPercentage);

      pie.props.showValue = true;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.showValue);

      pie.props.showLabel = true;
      pie.props.showPercentage = true;

      return r;
    });
    expect(result).toEqual(['Inside', 'Outside', 'Auto', false, false, true]);
  });

  // #### Outline
  await softStep('Outline', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const r: number[] = [];

      pie.props.outlineLineWidth = 5;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.outlineLineWidth);

      pie.props.outlineLineWidth = 0;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.outlineLineWidth);

      pie.props.outlineLineWidth = 1;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.outlineLineWidth);

      return r;
    });
    expect(result).toEqual([5, 0, 1]);
  });

  // #### Include nulls
  await softStep('Include nulls', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      pie.props.categoryColumnName = 'DIS_POP';
      const r: boolean[] = [];

      pie.props.includeNulls = true;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.includeNulls);

      pie.props.includeNulls = false;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.includeNulls);

      pie.props.includeNulls = true;
      pie.props.categoryColumnName = 'RACE';
      return r;
    });
    expect(result).toEqual([true, false]);
  });

  // #### Column selector
  await softStep('Column selector', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const r: boolean[] = [];

      pie.props.showColumnSelector = false;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.showColumnSelector);

      pie.props.showColumnSelector = true;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.showColumnSelector);

      return r;
    });
    expect(result).toEqual([false, true]);
  });

  // #### Legend
  await softStep('Legend', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const r: any[] = [];

      pie.props.legendVisibility = 'Always';
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.legendVisibility);

      pie.props.legendPosition = 'LeftTop';
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.legendPosition);

      pie.props.legendPosition = 'RightBottom';
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.legendPosition);

      pie.props.legendVisibility = 'Never';
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.legendVisibility);

      pie.props.legendVisibility = 'Auto';
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.legendVisibility);

      return r;
    });
    expect(result).toEqual(['Always', 'LeftTop', 'RightBottom', 'Never', 'Auto']);
  });

  // #### Category map (dates)
  await softStep('Category map (dates)', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      pie.props.categoryColumnName = 'STARTED';
      await new Promise(res => setTimeout(res, 500));
      const r: string[] = [];

      r.push(pie.props.categoryMap); // default year

      pie.props.categoryMap = 'month';
      await new Promise(res => setTimeout(res, 300));
      r.push(pie.props.categoryMap);

      pie.props.categoryMap = 'quarter';
      await new Promise(res => setTimeout(res, 300));
      r.push(pie.props.categoryMap);

      pie.props.categoryColumnName = 'RACE';
      return r;
    });
    expect(result).toEqual(['year', 'month', 'quarter']);
  });

  // #### Row source
  await softStep('Row source', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const df = grok.shell.tv.dataFrame;
      df.selection.init((i: number) => i < 50);
      await new Promise(res => setTimeout(res, 300));
      const r: any[] = [];

      pie.props.rowSource = 'Selected';
      await new Promise(res => setTimeout(res, 300));
      r.push(pie.props.rowSource);

      pie.props.rowSource = 'Filtered';
      await new Promise(res => setTimeout(res, 300));
      r.push(pie.props.rowSource);

      pie.props.rowSource = 'All';
      await new Promise(res => setTimeout(res, 300));
      r.push(pie.props.rowSource);

      df.selection.setAll(false);
      return r;
    });
    expect(result).toEqual(['Selected', 'Filtered', 'All']);
  });

  // #### Aggregation functions
  await softStep('Aggregation functions', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      pie.props.segmentAngleColumnName = 'AGE';
      const r: string[] = [];
      for (const aggr of ['avg', 'min', 'max', 'sum', 'med', 'stdev', 'count']) {
        pie.props.segmentAngleAggrType = aggr;
        await new Promise(res => setTimeout(res, 200));
        r.push(pie.props.segmentAngleAggrType);
      }
      pie.props.segmentAngleColumnName = '';
      return r;
    });
    expect(result).toEqual(['avg', 'min', 'max', 'sum', 'med', 'stdev', 'count']);
  });

  // #### Title and description
  await softStep('Title and description', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const r: any[] = [];

      pie.props.showTitle = true;
      r.push(pie.props.showTitle);

      pie.props.title = 'Demographics';
      r.push(pie.props.title);

      pie.props.description = 'By race';
      r.push(pie.props.description);

      pie.props.descriptionPosition = 'Top';
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.descriptionPosition);

      pie.props.descriptionVisibilityMode = 'Never';
      r.push(pie.props.descriptionVisibilityMode);

      pie.props.title = '';
      pie.props.showTitle = false;
      return r;
    });
    expect(result[0]).toBe(true);
    expect(result[1]).toBe('Demographics');
    expect(result[2]).toBe('By race');
    expect(result[3]).toBe('Top');
    expect(result[4]).toBe('Never');
  });

  // #### Layout persistence
  await softStep('Layout persistence', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;

      pie.props.categoryColumnName = 'RACE';
      pie.props.segmentAngleColumnName = 'AGE';
      pie.props.startAngle = 45;
      pie.props.shift = 5;
      await new Promise(res => setTimeout(res, 500));

      const before = {
        cat: pie.props.categoryColumnName,
        angle: pie.props.segmentAngleColumnName,
        startAngle: pie.props.startAngle,
        shift: pie.props.shift,
      };

      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(res => setTimeout(res, 1000));

      pie.close();
      await new Promise(res => setTimeout(res, 500));

      const saved = await grok.dapi.layouts.find(layoutId);
      grok.shell.tv.loadLayout(saved);
      await new Promise(res => setTimeout(res, 3000));

      const pie2 = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const after = pie2 ? {
        cat: pie2.props.categoryColumnName,
        angle: pie2.props.segmentAngleColumnName,
        startAngle: pie2.props.startAngle,
        shift: pie2.props.shift,
      } : null;

      await grok.dapi.layouts.delete(saved);
      return {before, after};
    });
    expect(result.after).toEqual(result.before);
  });

  // #### Selection and interaction
  await softStep('Selection and interaction', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      pie.props.categoryColumnName = 'RACE';
      const df = grok.shell.tv.dataFrame;

      const pieEl = document.querySelector('[name="viewer-Pie-chart"]') as HTMLElement;
      const canvas = pieEl.querySelector('canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();

      // Click slice
      const x = rect.left + rect.width * 0.65;
      const y = rect.top + rect.height * 0.4;
      canvas.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, clientX: x, clientY: y}));
      canvas.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, clientX: x, clientY: y}));
      canvas.dispatchEvent(new MouseEvent('click', {bubbles: true, clientX: x, clientY: y}));
      await new Promise(r => setTimeout(r, 500));
      const sel1 = df.selection.trueCount;

      // Toggle showSelectedRows
      pie.props.showSelectedRows = false;
      const sOff = pie.props.showSelectedRows;
      pie.props.showSelectedRows = true;
      const sOn = pie.props.showSelectedRows;

      // Toggle showMouseOverRowGroup
      pie.props.showMouseOverRowGroup = false;
      const mOff = pie.props.showMouseOverRowGroup;
      pie.props.showMouseOverRowGroup = true;
      const mOn = pie.props.showMouseOverRowGroup;

      df.selection.setAll(false);
      return {sel1, sOff, sOn, mOff, mOn};
    });
    expect(result.sel1).toBeGreaterThan(0);
    expect(result.sOff).toBe(false);
    expect(result.sOn).toBe(true);
    expect(result.mOff).toBe(false);
    expect(result.mOn).toBe(true);
  });

  // #### On Click modes
  await softStep('On Click modes', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      pie.props.categoryColumnName = 'RACE';
      const df = grok.shell.tv.dataFrame;
      const r: any[] = [];

      pie.props.onClick = 'Select';
      const pieEl = document.querySelector('[name="viewer-Pie-chart"]') as HTMLElement;
      const canvas = pieEl.querySelector('canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      const x = rect.left + rect.width * 0.65;
      const y = rect.top + rect.height * 0.4;

      canvas.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, clientX: x, clientY: y}));
      canvas.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, clientX: x, clientY: y}));
      canvas.dispatchEvent(new MouseEvent('click', {bubbles: true, clientX: x, clientY: y}));
      await new Promise(res => setTimeout(res, 500));
      r.push({mode: 'Select', selected: df.selection.trueCount});

      pie.props.onClick = 'Filter';
      canvas.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, clientX: x, clientY: y}));
      canvas.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, clientX: x, clientY: y}));
      canvas.dispatchEvent(new MouseEvent('click', {bubbles: true, clientX: x, clientY: y}));
      await new Promise(res => setTimeout(res, 500));
      r.push({mode: 'Filter', filtered: df.filter.trueCount, total: df.rowCount});

      // Click empty area to clear
      const ex = rect.left + 5;
      const ey = rect.top + 5;
      canvas.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, clientX: ex, clientY: ey}));
      canvas.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, clientX: ex, clientY: ey}));
      canvas.dispatchEvent(new MouseEvent('click', {bubbles: true, clientX: ex, clientY: ey}));
      await new Promise(res => setTimeout(res, 500));
      r.push({mode: 'clear', filtered: df.filter.trueCount});

      pie.props.onClick = 'Select';
      df.selection.setAll(false);
      return r;
    });
    expect(result[0].selected).toBeGreaterThan(0);
    expect(result[1].filtered).toBeLessThan(result[1].total);
    expect(result[2].filtered).toBe(result[1].total);
  });

  // #### Selection between grid and pie chart
  await softStep('Selection between grid and pie chart', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.selection.init((i: number) => i < 50);
      await new Promise(res => setTimeout(res, 300));
      const gridSel = df.selection.trueCount;

      const pieEl = document.querySelector('[name="viewer-Pie-chart"]') as HTMLElement;
      const canvas = pieEl.querySelector('canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      const x = rect.left + rect.width * 0.65;
      const y = rect.top + rect.height * 0.4;
      canvas.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, clientX: x, clientY: y}));
      canvas.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, clientX: x, clientY: y}));
      canvas.dispatchEvent(new MouseEvent('click', {bubbles: true, clientX: x, clientY: y}));
      await new Promise(res => setTimeout(res, 500));
      const pieSel = df.selection.trueCount;

      df.selection.setAll(false);
      return {gridSel, pieSel};
    });
    expect(result.gridSel).toBe(50);
    expect(result.pieSel).toBeGreaterThan(0);
  });

  // #### Auto layout
  await softStep('Auto layout', async () => {
    const result = await page.evaluate(async () => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const r: any[] = [];

      pie.props.autoLayout = false;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.autoLayout);

      pie.props.marginLeft = 50;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.marginLeft);

      pie.props.marginTop = 50;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.marginTop);

      pie.props.autoLayout = true;
      await new Promise(res => setTimeout(res, 200));
      r.push(pie.props.autoLayout);

      return r;
    });
    expect(result).toEqual([false, 50, 50, true]);
  });

  // #### Table switching and row source (SPGI)
  await softStep('Table switching and row source (SPGI)', async () => {
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

      const icon = document.querySelector('[name="icon-pie-chart"]') as HTMLElement;
      icon.click();
      await new Promise(r => setTimeout(r, 1000));

      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const r: any[] = [];

      // Switch to SPGI
      const spgi = Array.from(grok.shell.tables).find((t: any) => t.name === 'SPGI') as any;
      pie.dataFrame = spgi;
      await new Promise(res => setTimeout(res, 500));
      r.push(pie.dataFrame.name);

      // Switch back
      pie.dataFrame = df;
      await new Promise(res => setTimeout(res, 500));
      r.push(pie.dataFrame.name);

      // Row Source = Selected
      pie.props.rowSource = 'Selected';
      df.selection.init((i: number) => i < 100);
      await new Promise(res => setTimeout(res, 300));
      r.push({rowSource: pie.props.rowSource, selCount: df.selection.trueCount});

      // Row Source = Filtered + filter
      pie.props.rowSource = 'Filtered';
      tv.getFiltersGroup();
      await new Promise(res => setTimeout(res, 1000));
      const fg = tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: ['Asian']});
      await new Promise(res => setTimeout(res, 500));
      r.push({rowSource: pie.props.rowSource, filtered: df.filter.trueCount});

      // Reset
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: df.col('RACE').categories});
      df.selection.setAll(false);
      pie.props.rowSource = 'All';

      return r;
    });
    expect(result[0]).toBe('SPGI');
    expect(result[2].rowSource).toBe('Selected');
    expect(result[2].selCount).toBe(100);
    expect(result[3].rowSource).toBe('Filtered');
    expect(result[3].filtered).toBeGreaterThan(0);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
