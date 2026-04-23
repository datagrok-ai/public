import {test, expect} from '@playwright/test';

test.use({
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/demog.csv';
const curvesPath = 'System:DemoFiles/curves.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Trellis plot tests', async ({page}) => {
  // Phase 1: Navigate
  page.setDefaultTimeout(120000);
  await page.goto(baseUrl, {timeout: 120000, waitUntil: 'domcontentloaded'});
  await page.waitForFunction(() => {
    try {
      return typeof grok !== 'undefined'
        && grok.shell
        && typeof grok.shell.closeAll === 'function'
        && grok.dapi
        && grok.dapi.files
        && document.querySelector('.d4-root') !== null;
    } catch { return false; }
  }, {timeout: 120000});
  await page.waitForTimeout(5000);

  // Phase 2: Open dataset
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch {}
    try { grok.shell.windows.simpleMode = false; } catch {}
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

  // Phase 3: Add Trellis plot
  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-trellis-plot"]') as HTMLElement;
    icon.click();
  });
  await page.locator('[name="viewer-Trellis-plot"]').waitFor({timeout: 10000});

  // #### Inner viewer types
  await softStep('Inner viewer types', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const r: any[] = [];

      // Scatter plot
      tp.props.viewerType = 'Scatter plot';
      await new Promise(res => setTimeout(res, 800));
      tp.setOptions({innerViewerLook: {xColumnName: 'WEIGHT', yColumnName: 'HEIGHT', colorColumnName: 'RACE'}});
      await new Promise(res => setTimeout(res, 800));
      r.push(tp.props.viewerType);

      // Bar chart
      tp.props.viewerType = 'Bar chart';
      await new Promise(res => setTimeout(res, 800));
      tp.setOptions({innerViewerLook: {splitColumnName: 'RACE', valueColumnName: 'AGE', valueAggrType: 'avg'}});
      await new Promise(res => setTimeout(res, 800));
      r.push(tp.props.viewerType);

      // Histogram
      tp.props.viewerType = 'Histogram';
      await new Promise(res => setTimeout(res, 800));
      tp.setOptions({innerViewerLook: {valueColumnName: 'AGE', splitColumnName: 'RACE'}});
      await new Promise(res => setTimeout(res, 800));
      r.push(tp.props.viewerType);

      // Line chart
      tp.props.viewerType = 'Line chart';
      await new Promise(res => setTimeout(res, 800));
      tp.setOptions({innerViewerLook: {xColumnName: 'STARTED', yColumnName: 'AGE', splitColumnName: 'RACE'}});
      await new Promise(res => setTimeout(res, 800));
      r.push(tp.props.viewerType);

      // Box plot
      tp.props.viewerType = 'Box plot';
      await new Promise(res => setTimeout(res, 800));
      tp.setOptions({innerViewerLook: {categoryColumnName: 'SEX', valueColumnName: 'AGE'}});
      await new Promise(res => setTimeout(res, 800));
      r.push(tp.props.viewerType);

      // Pie chart
      tp.props.viewerType = 'Pie chart';
      await new Promise(res => setTimeout(res, 800));
      tp.setOptions({innerViewerLook: {categoryColumnName: 'RACE'}});
      await new Promise(res => setTimeout(res, 800));
      r.push(tp.props.viewerType);

      // Density plot
      tp.props.viewerType = 'Density plot';
      await new Promise(res => setTimeout(res, 800));
      tp.setOptions({innerViewerLook: {xColumnName: 'WEIGHT', yColumnName: 'HEIGHT'}});
      await new Promise(res => setTimeout(res, 800));
      r.push(tp.props.viewerType);

      // Statistics (Summary)
      tp.props.viewerType = 'Statistics';
      await new Promise(res => setTimeout(res, 800));
      tp.setOptions({innerViewerLook: {visualization: 'bars'}});
      await new Promise(res => setTimeout(res, 500));
      r.push(tp.props.viewerType);

      // Sparklines
      tp.props.viewerType = 'Sparklines';
      await new Promise(res => setTimeout(res, 800));
      tp.setOptions({innerViewerLook: {sparklineType: 'Bar Chart'}});
      await new Promise(res => setTimeout(res, 500));
      r.push(tp.props.viewerType);

      // PC Plot
      tp.props.viewerType = 'PC Plot';
      await new Promise(res => setTimeout(res, 800));
      tp.setOptions({innerViewerLook: {colorColumnName: 'SEX'}});
      await new Promise(res => setTimeout(res, 500));
      r.push(tp.props.viewerType);

      return r;
    });
    expect(result).toEqual([
      'Scatter plot', 'Bar chart', 'Histogram', 'Line chart', 'Box plot',
      'Pie chart', 'Density plot', 'Statistics', 'Sparklines', 'PC Plot',
    ]);
  });

  // #### Global scale
  await softStep('Global scale', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.viewerType = 'Scatter plot';
      await new Promise(r => setTimeout(r, 800));
      const r: boolean[] = [];

      tp.props.globalScale = true;
      await new Promise(res => setTimeout(res, 500));
      r.push(tp.props.globalScale);

      tp.props.globalScale = false;
      await new Promise(res => setTimeout(res, 500));
      r.push(tp.props.globalScale);

      tp.props.globalScale = true;
      await new Promise(res => setTimeout(res, 500));
      r.push(tp.props.globalScale);

      return r;
    });
    expect(result).toEqual([true, false, true]);
  });

  // #### Axes visibility
  await softStep('Axes visibility', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const r: any[] = [];

      for (const val of ['Always', 'Never', 'Auto']) {
        tp.props.showXAxes = val;
        await new Promise(res => setTimeout(res, 300));
        r.push(tp.props.showXAxes);
      }
      for (const val of ['Always', 'Never', 'Auto']) {
        tp.props.showYAxes = val;
        await new Promise(res => setTimeout(res, 300));
        r.push(tp.props.showYAxes);
      }

      tp.props.showRangeSliders = false;
      await new Promise(res => setTimeout(res, 300));
      r.push(tp.props.showRangeSliders);

      tp.props.showRangeSliders = true;
      await new Promise(res => setTimeout(res, 300));
      r.push(tp.props.showRangeSliders);

      return r;
    });
    expect(result).toEqual(['Always', 'Never', 'Auto', 'Always', 'Never', 'Auto', false, true]);
  });

  // #### Range sliders with global scale
  await softStep('Range sliders with global scale', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.globalScale = true;
      tp.props.showRangeSliders = true;
      tp.props.showXAxes = 'Always';
      tp.props.showYAxes = 'Always';
      await new Promise(r => setTimeout(r, 500));
      return {
        globalScale: tp.props.globalScale,
        showRangeSliders: tp.props.showRangeSliders,
        showXAxes: tp.props.showXAxes,
        showYAxes: tp.props.showYAxes,
      };
      // Note: hover/drag range sliders and Reset Inner Range Sliders are canvas interactions -- not automatable
    });
    expect(result.globalScale).toBe(true);
    expect(result.showRangeSliders).toBe(true);
    expect(result.showXAxes).toBe('Always');
    expect(result.showYAxes).toBe('Always');
  });

  // #### Gridlines
  await softStep('Gridlines', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const r: string[] = [];
      for (const val of ['always', 'never', 'auto']) {
        tp.props.showGridlines = val;
        await new Promise(res => setTimeout(res, 300));
        r.push(tp.props.showGridlines);
      }
      return r;
    });
    expect(result).toEqual(['always', 'never', 'auto']);
  });

  // #### Tiles mode
  await softStep('Tiles mode', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const r: any[] = [];

      tp.props.useTiledView = true;
      await new Promise(res => setTimeout(res, 500));
      r.push(tp.props.useTiledView);

      tp.props.tilesPerRow = 2;
      await new Promise(res => setTimeout(res, 300));
      r.push(tp.props.tilesPerRow);

      tp.props.tilesPerRow = 6;
      await new Promise(res => setTimeout(res, 300));
      r.push(tp.props.tilesPerRow);

      tp.props.useTiledView = false;
      await new Promise(res => setTimeout(res, 300));
      r.push(tp.props.useTiledView);

      return r;
    });
    expect(result).toEqual([true, 2, 6, false]);
  });

  // #### Category management
  await softStep('Category management', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const r: any[] = [];

      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      await new Promise(res => setTimeout(res, 800));
      r.push({x: [...tp.props.xColumnNames], y: [...tp.props.yColumnNames]});

      tp.props.xColumnNames = ['SEX', 'DIS_POP'];
      await new Promise(res => setTimeout(res, 500));
      r.push([...tp.props.xColumnNames]);

      tp.props.xColumnNames = ['SEX'];
      await new Promise(res => setTimeout(res, 500));
      r.push([...tp.props.xColumnNames]);

      tp.props.showXLabels = false;
      tp.props.showYLabels = false;
      await new Promise(res => setTimeout(res, 300));
      r.push({showXLabels: tp.props.showXLabels, showYLabels: tp.props.showYLabels});

      tp.props.showXLabels = true;
      tp.props.showYLabels = true;
      await new Promise(res => setTimeout(res, 300));
      r.push({showXLabels: tp.props.showXLabels, showYLabels: tp.props.showYLabels});

      return r;
    });
    expect(result[0]).toEqual({x: ['SEX'], y: ['RACE']});
    expect(result[1]).toEqual(['SEX', 'DIS_POP']);
    expect(result[2]).toEqual(['SEX']);
    expect(result[3]).toEqual({showXLabels: false, showYLabels: false});
    expect(result[4]).toEqual({showXLabels: true, showYLabels: true});
  });

  // #### Pack categories
  await softStep('Pack categories', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const df = grok.shell.tv.dataFrame;
      tp.props.xColumnNames = ['SEX'];
      tp.props.yColumnNames = ['RACE'];
      await new Promise(r => setTimeout(r, 800));

      // Open filter panel and filter out Asian
      grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 2000));
      const fg = grok.shell.tv.getFiltersGroup();
      const cats = df.col('RACE').categories;
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: cats.filter((c: string) => c !== 'Asian')});
      await new Promise(r => setTimeout(r, 1000));

      const r: any[] = [];
      r.push({filteredRows: df.filter.trueCount, packCategories: tp.props.packCategories});

      tp.props.packCategories = false;
      await new Promise(res => setTimeout(res, 500));
      r.push(tp.props.packCategories);

      tp.props.packCategories = true;
      await new Promise(res => setTimeout(res, 500));
      r.push(tp.props.packCategories);

      // Reset filter
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: cats});
      await new Promise(r => setTimeout(r, 500));

      return r;
    });
    expect(result[0].packCategories).toBe(true);
    expect(result[0].filteredRows).toBeLessThan(5850);
    expect(result[1]).toBe(false);
    expect(result[2]).toBe(true);
  });

  // #### On Click functionality
  await softStep('On Click functionality', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const r: string[] = [];

      tp.props.onClick = 'Select';
      await new Promise(res => setTimeout(res, 300));
      r.push(tp.props.onClick);

      tp.props.onClick = 'Filter';
      await new Promise(res => setTimeout(res, 300));
      r.push(tp.props.onClick);

      tp.props.onClick = 'None';
      await new Promise(res => setTimeout(res, 300));
      r.push(tp.props.onClick);

      return r;
      // Note: actual cell clicks are canvas-based and not automatable
    });
    expect(result).toEqual(['Select', 'Filter', 'None']);
  });

  // #### Selectors
  await softStep('Selectors', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const r: any[] = [];

      tp.props.showXSelectors = false;
      tp.props.showYSelectors = false;
      tp.props.showControlPanel = false;
      await new Promise(res => setTimeout(res, 300));
      r.push({x: tp.props.showXSelectors, y: tp.props.showYSelectors, cp: tp.props.showControlPanel});

      tp.props.showXSelectors = true;
      tp.props.showYSelectors = true;
      tp.props.showControlPanel = true;
      await new Promise(res => setTimeout(res, 300));
      r.push({x: tp.props.showXSelectors, y: tp.props.showYSelectors, cp: tp.props.showControlPanel});

      return r;
    });
    expect(result[0]).toEqual({x: false, y: false, cp: false});
    expect(result[1]).toEqual({x: true, y: true, cp: true});
  });

  // #### Allow viewer full screen
  await softStep('Allow viewer full screen', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      // Note: hover/click full-screen icon are canvas interactions -- not automatable
      tp.props.allowViewerFullScreen = false;
      await new Promise(r => setTimeout(r, 300));
      const off = tp.props.allowViewerFullScreen;
      tp.props.allowViewerFullScreen = true;
      await new Promise(r => setTimeout(r, 300));
      return {off, on: tp.props.allowViewerFullScreen};
    });
    expect(result.off).toBe(false);
    expect(result.on).toBe(true);
  });

  // #### Scrolling
  await softStep('Scrolling', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['SITE'];
      await new Promise(r => setTimeout(r, 1000));
      const x1 = [...tp.props.xColumnNames];

      tp.props.yColumnNames = ['RACE'];
      await new Promise(r => setTimeout(r, 500));
      const y1 = [...tp.props.yColumnNames];

      // Reset
      tp.props.xColumnNames = ['SEX'];
      await new Promise(r => setTimeout(r, 500));

      return {x1, y1};
      // Note: scroll slider drag and mouse wheel are canvas interactions
    });
    expect(result.x1).toEqual(['SITE']);
    expect(result.y1).toEqual(['RACE']);
  });

  // #### Legend
  await softStep('Legend', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.viewerType = 'Scatter plot';
      await new Promise(r => setTimeout(r, 800));
      tp.setOptions({innerViewerLook: {colorColumnName: 'SEX'}});
      await new Promise(r => setTimeout(r, 800));

      const r: any[] = [];
      tp.props.legendVisibility = 'Always';
      await new Promise(res => setTimeout(res, 300));
      r.push(tp.props.legendVisibility);

      for (const pos of ['Left', 'Right', 'Top', 'Bottom']) {
        tp.props.legendPosition = pos;
        await new Promise(res => setTimeout(res, 300));
        r.push(tp.props.legendPosition);
      }

      tp.props.legendVisibility = 'Never';
      await new Promise(res => setTimeout(res, 300));
      r.push(tp.props.legendVisibility);

      return r;
    });
    expect(result[0]).toBe('Always');
    expect(result.slice(1, 5)).toEqual(['Left', 'Right', 'Top', 'Bottom']);
    expect(result[5]).toBe('Never');
  });

  // #### Context menu
  await softStep('Context menu', async () => {
    const result = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const cells = root.querySelectorAll('.d4-trellis-plot-cell');
      if (cells.length === 0) return {error: 'No trellis cells'};
      const cell = cells[0];
      const rect = cell.getBoundingClientRect();
      cell.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2,
      }));
      await new Promise(r => setTimeout(r, 1000));
      const labels = Array.from(document.querySelectorAll('.d4-menu-item-label')).map(el => (el as HTMLElement).textContent?.trim());
      const hasScatterPlot = labels.includes('Scatter plot');
      const hasProperties = labels.includes('Properties...');
      document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
      return {hasScatterPlot, hasProperties, totalItems: labels.length};
    });
    expect(result.hasScatterPlot).toBe(true);
  });

  // #### Inner viewer properties (via JS API -- gear icon hidden in DOM)
  await softStep('Inner viewer properties', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.setOptions({innerViewerLook: {xColumnName: 'AGE', yColumnName: 'WEIGHT'}});
      await new Promise(r => setTimeout(r, 800));
      return {viewerType: tp.props.viewerType};
    });
    expect(result.viewerType).toBe('Scatter plot');
  });

  // #### Use in Trellis (Scatter plot only -- other types tested in MCP run)
  await softStep('Use in Trellis', async () => {
    const result = await page.evaluate(async () => {
      // Close current trellis
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.close();
      await new Promise(res => setTimeout(res, 500));

      // Add Scatter plot
      const sp = grok.shell.tv.addViewer('Scatter plot') as any;
      sp.setOptions({xColumnName: 'AGE', yColumnName: 'HEIGHT', colorColumnName: 'SEX'});
      await new Promise(res => setTimeout(res, 1000));

      const root = document.querySelector('[name="viewer-Scatter-plot"]') as HTMLElement;
      if (!root) return {error: 'Scatter plot root not found'};
      const canvas = root.querySelector('canvas') as HTMLElement;
      if (!canvas) return {error: 'No canvas'};
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2,
      }));
      await new Promise(res => setTimeout(res, 800));

      let labels = Array.from(document.querySelectorAll('.d4-menu-item-label'));
      const general = labels.find(el => (el as HTMLElement).textContent?.trim() === 'General');
      if (!general) {
        document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
        return {error: 'General not found'};
      }
      (general.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('mousemove', {bubbles: true}));
      await new Promise(res => setTimeout(res, 500));
      labels = Array.from(document.querySelectorAll('.d4-menu-item-label'));
      const uit = labels.find(el => (el as HTMLElement).textContent?.trim() === 'Use in Trellis');
      if (!uit) {
        document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
        return {error: 'Use in Trellis not found'};
      }
      (uit.closest('.d4-menu-item') as HTMLElement).click();
      await new Promise(res => setTimeout(res, 2000));
      const trellisFound = Array.from(grok.shell.tv.viewers).some((vw: any) => vw.type === 'Trellis plot');

      // Cleanup: close scatter plot if still open
      for (const vw of Array.from(grok.shell.tv.viewers) as any[])
        if (vw.type === 'Scatter plot') vw.close();

      // Ensure trellis exists for remaining tests
      if (!Array.from(grok.shell.tv.viewers).some((vw: any) => vw.type === 'Trellis plot'))
        grok.shell.tv.addViewer('Trellis plot');
      await new Promise(r => setTimeout(r, 1000));

      return {trellisCreated: trellisFound};
    });
    expect(result.trellisCreated).toBe(true);
  });

  // #### Auto layout
  await softStep('Auto layout', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const r: boolean[] = [];
      r.push(tp.props.autoLayout); // default true

      tp.props.autoLayout = false;
      await new Promise(res => setTimeout(res, 300));
      r.push(tp.props.autoLayout);

      tp.props.autoLayout = true;
      await new Promise(res => setTimeout(res, 300));
      r.push(tp.props.autoLayout);

      return r;
      // Note: viewer resize to trigger auto layout is not automatable
    });
    expect(result).toEqual([true, false, true]);
  });

  // #### Title and description
  await softStep('Title and description', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const r: any[] = [];

      tp.props.showTitle = true;
      r.push(tp.props.showTitle);

      tp.props.title = 'My Trellis';
      r.push(tp.props.title);

      tp.props.description = 'Test description';
      r.push(tp.props.description);

      for (const pos of ['Bottom', 'Top', 'Left', 'Right']) {
        tp.props.descriptionPosition = pos;
        await new Promise(res => setTimeout(res, 200));
        r.push(tp.props.descriptionPosition);
      }

      return r;
    });
    expect(result[0]).toBe(true);
    expect(result[1]).toBe('My Trellis');
    expect(result[2]).toBe('Test description');
    expect(result.slice(3)).toEqual(['Bottom', 'Top', 'Left', 'Right']);
  });

  // #### Label orientation
  await softStep('Label orientation', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      tp.props.xColumnNames = ['RACE'];
      await new Promise(r => setTimeout(r, 500));

      const r: string[] = [];
      for (const val of ['Horz', 'Vert', 'Auto']) {
        tp.props.xLabelsOrientation = val;
        await new Promise(res => setTimeout(res, 300));
        r.push(tp.props.xLabelsOrientation);
      }
      for (const val of ['Horz', 'Vert', 'Auto']) {
        tp.props.yLabelsOrientation = val;
        await new Promise(res => setTimeout(res, 300));
        r.push(tp.props.yLabelsOrientation);
      }
      return r;
    });
    expect(result).toEqual(['Horz', 'Vert', 'Auto', 'Horz', 'Vert', 'Auto']);
  });

  // #### Pick Up / Apply
  await softStep('Pick Up / Apply', async () => {
    const result = await page.evaluate(async () => {
      const tp1 = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const tp2 = grok.shell.tv.addViewer('Trellis plot') as any;
      await new Promise(r => setTimeout(r, 1000));

      // Configure first trellis
      tp1.props.yColumnNames = ['RACE'];
      tp1.props.viewerType = 'Bar chart';
      tp1.props.legendVisibility = 'Always';
      tp1.props.legendPosition = 'Top';
      tp1.props.showTitle = true;
      tp1.props.title = 'First Trellis';
      await new Promise(r => setTimeout(r, 800));

      // Pick up from first
      const root1 = document.querySelectorAll('[name="viewer-Trellis-plot"]')[0] as HTMLElement;
      const canvas1 = root1.querySelector('canvas') as HTMLElement;
      const rect1 = canvas1.getBoundingClientRect();
      canvas1.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect1.left + rect1.width / 2, clientY: rect1.top + rect1.height / 2,
      }));
      await new Promise(r => setTimeout(r, 800));

      let labels = Array.from(document.querySelectorAll('.d4-menu-item-label'));
      const pickUpGroup = labels.find(el => (el as HTMLElement).textContent?.trim() === 'Pick Up / Apply');
      if (!pickUpGroup) return {error: 'Pick Up / Apply not found'};
      (pickUpGroup.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('mousemove', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      labels = Array.from(document.querySelectorAll('.d4-menu-item-label'));
      const pick = labels.find(el => (el as HTMLElement).textContent?.trim() === 'Pick Up');
      if (!pick) return {error: 'Pick Up not found'};
      (pick.closest('.d4-menu-item') as HTMLElement).click();
      await new Promise(r => setTimeout(r, 800));

      // Apply to second
      const root2 = document.querySelectorAll('[name="viewer-Trellis-plot"]')[1] as HTMLElement;
      if (!root2) return {error: 'Second trellis not found'};
      const canvas2 = root2.querySelector('canvas') || root2;
      const rect2 = canvas2.getBoundingClientRect();
      canvas2.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect2.left + rect2.width / 2, clientY: rect2.top + rect2.height / 2,
      }));
      await new Promise(r => setTimeout(r, 800));
      labels = Array.from(document.querySelectorAll('.d4-menu-item-label'));
      const applyGroup = labels.find(el => (el as HTMLElement).textContent?.trim() === 'Pick Up / Apply');
      if (!applyGroup) return {error: 'Pick Up / Apply not found on second'};
      (applyGroup.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('mousemove', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      labels = Array.from(document.querySelectorAll('.d4-menu-item-label'));
      const apply = labels.find(el => (el as HTMLElement).textContent?.trim() === 'Apply');
      if (!apply) return {error: 'Apply not found'};
      (apply.closest('.d4-menu-item') as HTMLElement).click();
      await new Promise(r => setTimeout(r, 1000));

      const tps = Array.from(grok.shell.tv.viewers).filter((v: any) => v.type === 'Trellis plot') as any[];
      const result = {
        tp2Type: tps[1]?.props.viewerType,
        tp2Title: tps[1]?.props.title,
      };

      // Cleanup second trellis
      if (tps.length > 1) tps[1].close();

      return result;
    });
    expect(result.tp2Type).toBe('Bar chart');
    expect(result.tp2Title).toBe('First Trellis');
  });

  // #### Layout and Project save/restore
  await softStep('Layout save/restore', async () => {
    const result = await page.evaluate(async () => {
      const r: any = {};

      // Save layout
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(res => setTimeout(res, 1000));
      r.layoutSaved = true;

      // Add extra viewers
      grok.shell.tv.addViewer('Histogram');
      grok.shell.tv.addViewer('Bar chart');
      await new Promise(res => setTimeout(res, 1000));
      r.viewersBefore = Array.from(grok.shell.tv.viewers).map((v: any) => v.type);

      // Apply saved layout
      const saved = await grok.dapi.layouts.find(layoutId);
      grok.shell.tv.loadLayout(saved);
      await new Promise(res => setTimeout(res, 3000));
      r.viewersAfter = Array.from(grok.shell.tv.viewers).map((v: any) => v.type);

      await grok.dapi.layouts.delete(saved);
      return r;
    });
    expect(result.layoutSaved).toBe(true);
    expect(result.viewersBefore.length).toBeGreaterThan(result.viewersAfter.length);
    expect(result.viewersAfter).toContain('Trellis plot');
  });

  // #### Viewer filter formula
  await softStep('Viewer filter formula', async () => {
    const result = await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      const r: string[] = [];

      tp.props.filter = '${AGE} > 40';
      await new Promise(res => setTimeout(res, 500));
      r.push(tp.props.filter);

      tp.props.filter = '';
      await new Promise(res => setTimeout(res, 500));
      r.push(tp.props.filter);

      return r;
    });
    expect(result[0]).toBe('${AGE} > 40');
    expect(result[1]).toBe('');
  });

  // #### Multi Curve inner viewer (and table switching)
  await softStep('Multi Curve inner viewer', async () => {
    const result = await page.evaluate(async (cPath) => {
      const dfCurves = await grok.dapi.files.readCsv(cPath);
      grok.shell.addTableView(dfCurves);
      await new Promise(r => setTimeout(r, 2000));

      // Go back to demog
      const views = Array.from(grok.shell.views) as any[];
      const demog = views.find(v => v.name && v.name.includes('Table'));
      if (demog) grok.shell.v = demog;
      await new Promise(r => setTimeout(r, 500));

      const tp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Trellis plot') as any;
      if (!tp) return {error: 'Trellis plot not found'};

      // Switch table
      try { tp.props.table = dfCurves.name; } catch {}
      await new Promise(r => setTimeout(r, 1000));

      // Set Multi curve viewer
      tp.props.viewerType = 'Multi curve viewer';
      await new Promise(r => setTimeout(r, 1000));
      const vt = tp.props.viewerType;

      // Reset
      tp.props.viewerType = 'Scatter plot';
      try { tp.props.table = 'demog'; } catch {}

      return {viewerType: vt, curvesRows: dfCurves.rowCount};
    }, curvesPath);
    expect(result.viewerType).toBe('Multi curve viewer');
  });

  // #### To Script
  await softStep('To Script', async () => {
    const result = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
      const canvas = root.querySelector('canvas') as HTMLElement;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2,
      }));
      await new Promise(r => setTimeout(r, 800));

      let labels = Array.from(document.querySelectorAll('.d4-menu-item-label'));
      const toScript = labels.find(el => (el as HTMLElement).textContent?.trim() === 'To Script');
      if (!toScript) {
        document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
        return {error: 'To Script not found'};
      }
      (toScript.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('mousemove', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      labels = Array.from(document.querySelectorAll('.d4-menu-item-label'));
      const toJs = labels.find(el => (el as HTMLElement).textContent?.trim() === 'To JavaScript');
      if (toJs) {
        (toJs.closest('.d4-menu-item') as HTMLElement).click();
        await new Promise(r => setTimeout(r, 1500));
        const balloon = document.querySelector('.d4-balloon');
        if (balloon) {
          const close = balloon.querySelector('.close') || balloon.querySelector('[name="icon-times"]');
          if (close) (close as HTMLElement).click();
        }
        return {scriptGenerated: true};
      }
      document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
      return {error: 'To JavaScript not found'};
    });
    expect(result.scriptGenerated).toBe(true);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
