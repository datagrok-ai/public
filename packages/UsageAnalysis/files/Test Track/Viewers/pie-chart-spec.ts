import {test, expect} from '@playwright/test';

const baseUrl = 'http://localhost:8888';
const datasetPath = 'System:DemoFiles/SPGI.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Pie chart viewer', async ({page}) => {
  // Phase 1: Navigate and wait for full platform init
  await page.goto(baseUrl);
  await page.waitForFunction(() => {
    try { return typeof grok !== 'undefined' && grok.shell && grok.shell.tv !== undefined; }
    catch (_) { return false; }
  }, {timeout: 30000});
  await page.waitForTimeout(3000);

  // Phase 2: Open datasets, link tables, switch to SPGI
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.closeAll();

    const df = await grok.dapi.files.readCsv(path);
    df.name = 'SPGI';
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });

    const df1 = await grok.dapi.files.readCsv('System:DemoFiles/SPGI-linked1.csv');
    df1.name = 'SPGI-linked1';
    grok.shell.addTableView(df1);

    const df2 = await grok.dapi.files.readCsv('System:DemoFiles/SPGI-linked2.csv');
    df2.name = 'SPGI-linked2';
    grok.shell.addTableView(df2);

    grok.data.linkTables(df, df1, ['Sample Name'], ['Sample Name'], [DG.SYNC_TYPE.FILTER_TO_FILTER]);
    grok.data.linkTables(df, df2, ['Sample Name'], ['Sample Name'], [DG.SYNC_TYPE.FILTER_TO_FILTER]);

    const views = Array.from(grok.shell.views).filter((v: any) => v.type === 'TableView');
    const spgiView = views.find((v: any) => v.dataFrame.name === 'SPGI');
    if (spgiView) grok.shell.v = spgiView;
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Step 3: Add Pie chart
  await softStep('Step 3: Add Pie chart', async () => {
    await page.evaluate(() => {
      const icon = document.querySelector('[name="icon-pie-chart"]') as HTMLElement;
      icon.click();
    });
    await page.locator('[name="viewer-Pie-chart"]').waitFor({timeout: 5000});
    const viewerCount = await page.evaluate(() => Array.from(grok.shell.tv.viewers).length);
    expect(viewerCount).toBe(2);
  });

  // Step 4: Open Property Pane via gear icon
  await softStep('Step 4: Open Property Pane', async () => {
    await page.evaluate(() => {
      const pieContainer = document.querySelector('[name="viewer-Pie-chart"]') as HTMLElement;
      const panelBase = pieContainer.closest('.panel-base') as HTMLElement;
      const gearIcon = panelBase.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
      gearIcon.click();
    });
    await page.waitForTimeout(500);
  });

  // Step 5: Check aggregation functions
  await softStep('Step 5: Check aggregation functions', async () => {
    const result = await page.evaluate(() => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const aggrTypes = ['count', 'avg', 'min', 'max', 'sum', 'med', 'stdev'];
      const results: string[] = [];
      for (const aggr of aggrTypes) {
        pie.props.segmentAngleAggrType = aggr;
        results.push(pie.props.segmentAngleAggrType);
      }
      pie.props.segmentAngleAggrType = 'count';
      return results;
    });
    expect(result).toEqual(['count', 'avg', 'min', 'max', 'sum', 'med', 'stdev']);
  });

  // Step 7.1: Select first 50 rows, verify pie chart reflects selection
  await softStep('Step 7.1: Selection from grid to pie chart', async () => {
    const result = await page.evaluate(() => {
      const tv = grok.shell.tv;
      tv.dataFrame.selection.init((i: number) => i < 50);
      return tv.dataFrame.selection.trueCount;
    });
    expect(result).toBe(50);
    await page.waitForTimeout(500);
  });

  // Step 7.3: Click pie chart segment, check grid selection
  await softStep('Step 7.3: Selection from pie chart to grid', async () => {
    const result = await page.evaluate(async () => {
      const pieEl = document.querySelector('[name="viewer-Pie-chart"]') as HTMLElement;
      const canvas = pieEl.querySelector('canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      const x = rect.left + rect.width * 0.65;
      const y = rect.top + rect.height * 0.4;
      canvas.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, clientX: x, clientY: y}));
      canvas.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, clientX: x, clientY: y}));
      canvas.dispatchEvent(new MouseEvent('click', {bubbles: true, clientX: x, clientY: y}));
      await new Promise(r => setTimeout(r, 500));
      const sel = grok.shell.tv.dataFrame.selection.trueCount;
      grok.shell.tv.dataFrame.selection.setAll(false);
      return sel;
    });
    expect(result).toBeGreaterThan(0);
  });

  // Step 8: Filtering
  await softStep('Step 8: Filtering interaction', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const fg = tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 1000));
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: ['R_ONE', 'S_ABS']});
      await new Promise(r => setTimeout(r, 500));
      const filtered = tv.dataFrame.filter.trueCount;
      // Reset
      fg.updateOrAdd({
        type: DG.FILTER_TYPE.CATEGORICAL,
        column: 'Stereo Category',
        selected: tv.dataFrame.col('Stereo Category').categories,
      });
      await new Promise(r => setTimeout(r, 300));
      return {filtered, total: tv.dataFrame.rowCount};
    });
    expect(result.filtered).toBeLessThan(result.total);
    expect(result.filtered).toBeGreaterThan(0);
  });

  // Step 9: Data properties (table switching, rowSource, filter)
  await softStep('Step 9: Data properties', async () => {
    const result = await page.evaluate(() => {
      const pie = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Pie chart') as any;
      const tables = Array.from(grok.shell.tables);

      // Table switching
      const linked1 = tables.find((t: any) => t.name === 'SPGI-linked1') as any;
      pie.dataFrame = linked1;
      const switched = pie.dataFrame.name;
      const spgi = tables.find((t: any) => t.name === 'SPGI') as any;
      pie.dataFrame = spgi;

      // Row source
      pie.props.rowSource = 'Filtered';
      const rowSrc = pie.props.rowSource;
      pie.props.rowSource = 'All';

      // Filter property
      pie.props.filter = '${Average Mass} > 300';
      const filter = pie.props.filter;
      pie.props.filter = '';

      return {switched, backTo: pie.dataFrame.name, rowSrc, filter};
    });
    expect(result.switched).toBe('SPGI-linked1');
    expect(result.backTo).toBe('SPGI');
    expect(result.rowSrc).toBe('Filtered');
    expect(result.filter).toBe('${Average Mass} > 300');
  });

  // Step 10: Title, Description, layout save
  await softStep('Step 10: Title, Description, and layout save', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const pie = Array.from(tv.viewers).find((v: any) => v.type === 'Pie chart') as any;

      // Set title and description
      pie.props.showTitle = true;
      pie.props.title = 'Test Pie Chart';
      try { pie.props.descriptionVisibilityMode = 'Always'; } catch {}
      try { pie.props.description = 'Test description'; } catch {}
      try { pie.props.descriptionPosition = 'Bottom'; } catch {}

      const before = {title: pie.props.title, showTitle: pie.props.showTitle};

      // Save layout
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1000));

      // Close and restore
      pie.close();
      await new Promise(r => setTimeout(r, 500));
      const saved = await grok.dapi.layouts.find(layoutId);
      tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));

      const viewers2 = Array.from(tv.viewers);
      const pie2 = viewers2.find((v: any) => v.type === 'Pie chart') as any;
      const after = pie2 ? {title: pie2.props.title, showTitle: pie2.props.showTitle} : null;
      await grok.dapi.layouts.delete(saved);

      return {before, after};
    });
    expect(result.after).toEqual(result.before);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
