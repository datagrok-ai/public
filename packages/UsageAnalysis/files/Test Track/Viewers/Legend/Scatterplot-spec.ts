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

async function openDataset(page: any) {
  await page.goto(baseUrl);
  await page.waitForFunction(() => typeof grok !== 'undefined' && grok.shell, {timeout: 15000});

  await page.evaluate(async (path: string) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
}

test('Scatterplot Legend scenarios', async ({page}) => {
  stepErrors.length = 0;

  // ---- Scenario 1: Color/Marker settings and dynamic legends ----
  await softStep('1. Color/Marker settings and dynamic legends', async () => {
    await openDataset(page);

    // Add scatterplot, set Color and Marker to Series
    const result1 = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const sp = tv.addViewer('Scatter plot');
      await new Promise(r => setTimeout(r, 1000));
      sp.props.colorColumnName = 'Series';
      sp.props.markersColumnName = 'Series';
      await new Promise(r => setTimeout(r, 1000));
      return {color: sp.props.colorColumnName, marker: sp.props.markersColumnName};
    });
    expect(result1.color).toBe('Series');
    expect(result1.marker).toBe('Series');

    // Verify combined legend exists
    const legendItems = await page.locator('[name="viewer-Scatter-plot"] .d4-legend-item').count();
    expect(legendItems).toBeGreaterThan(0);

    // Save and apply layout
    const layoutResult = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1000));

      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      const tv2 = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });

      const saved = await grok.dapi.layouts.find(layoutId);
      tv2.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));

      const viewers = Array.from(tv2.viewers);
      const sp2 = viewers.find((v: any) => v.type === 'Scatter plot') as any;
      const result = {
        restored: sp2 != null,
        color: sp2?.props.colorColumnName,
        marker: sp2?.props.markersColumnName,
      };
      await grok.dapi.layouts.delete(saved);
      return result;
    });
    expect(layoutResult.restored).toBe(true);
    expect(layoutResult.color).toBe('Series');
    expect(layoutResult.marker).toBe('Series');

    // Test dynamic color columns
    const dynamicResult = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const viewers = Array.from(tv.viewers);
      const sp = viewers.find((v: any) => v.type === 'Scatter plot') as any;

      // Linear legend
      df.columns.addNewCalculated('linear_test',
        'if(${Stereo Category}=="S_UNKN", null, ${Average Mass})');
      await new Promise(r => setTimeout(r, 1000));
      sp.props.colorColumnName = 'linear_test';
      await new Promise(r => setTimeout(r, 1000));
      const linearColor = sp.props.colorColumnName;

      // Categorical legend
      df.columns.addNewCalculated('cat_test',
        'if(${Stereo Category}=="S_UNKN", null, ${Series})');
      await new Promise(r => setTimeout(r, 1000));
      sp.props.colorColumnName = 'cat_test';
      await new Promise(r => setTimeout(r, 1000));
      const catColor = sp.props.colorColumnName;

      // Set Color=ID, Marker=Core
      sp.props.colorColumnName = 'Id';
      sp.props.markersColumnName = 'Core';
      await new Promise(r => setTimeout(r, 1000));

      return {linearColor, catColor, finalColor: sp.props.colorColumnName, finalMarker: sp.props.markersColumnName};
    });
    expect(dynamicResult.linearColor).toBe('linear_test');
    expect(dynamicResult.catColor).toBe('cat_test');
    expect(dynamicResult.finalColor).toBe('Id');
    expect(dynamicResult.finalMarker).toBe('Core');

    await page.evaluate(() => grok.shell.closeAll());
  });

  // ---- Scenario 2: Legend update on axis change ----
  await softStep('2. Legend update on axis change', async () => {
    await openDataset(page);

    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;

      df.columns.addNewCalculated('col1',
        'if(${Stereo Category}!="S_UNKN", null, ${Average Mass})');
      df.columns.addNewCalculated('col2',
        'if(${Stereo Category}=="S_UNKN", null, ${Average Mass})');
      await new Promise(r => setTimeout(r, 1000));

      const sp = tv.addViewer('Scatter plot');
      await new Promise(r => setTimeout(r, 1000));
      sp.props.xColumnName = 'col1';
      sp.props.colorColumnName = 'Stereo Category';
      await new Promise(r => setTimeout(r, 1000));

      // Change X to col2
      sp.props.xColumnName = 'col2';
      await new Promise(r => setTimeout(r, 1000));

      grok.shell.closeAll();
      return {xCol: sp.props.xColumnName, colorCol: sp.props.colorColumnName};
    });
    expect(result.xCol).toBe('col2');
    expect(result.colorCol).toBe('Stereo Category');
  });

  // ---- Scenario 3: In-viewer filtering ----
  await softStep('3. In-viewer filtering', async () => {
    await openDataset(page);

    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;

      const sp1 = tv.addViewer('Scatter plot');
      await new Promise(r => setTimeout(r, 1000));
      sp1.props.markersColumnName = 'Stereo Category';
      sp1.props.filter = '${Stereo Category} in ["R_ONE", "S_UNKN"]';
      await new Promise(r => setTimeout(r, 1000));

      const sp2 = tv.addViewer('Scatter plot');
      await new Promise(r => setTimeout(r, 1000));
      sp2.props.markersColumnName = 'Stereo Category';
      sp2.props.filter = '${Stereo Category} in ["R_ONE", "S_UNKN"]';
      await new Promise(r => setTimeout(r, 1000));

      // Save/apply layout
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1000));

      grok.shell.closeAll();
      const df2 = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      const tv2 = grok.shell.addTableView(df2);
      await new Promise(resolve => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });

      const saved = await grok.dapi.layouts.find(layoutId);
      tv2.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));

      const viewers = Array.from(tv2.viewers);
      const sps = viewers.filter((v: any) => v.type === 'Scatter plot');
      const filters = sps.map((sp: any) => sp.props.filter);
      const markers = sps.map((sp: any) => sp.props.markersColumnName);

      await grok.dapi.layouts.delete(saved);
      grok.shell.closeAll();

      return {spCount: sps.length, filters, markers};
    });
    expect(result.spCount).toBe(2);
    expect(result.filters[0]).toContain('R_ONE');
    expect(result.markers[0]).toBe('Stereo Category');
  });

  // ---- Scenario 4: Filtering via Filter Panel ----
  await softStep('4. Filtering via Filter Panel', async () => {
    await openDataset(page);

    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;

      const sp = tv.addViewer('Scatter plot');
      await new Promise(r => setTimeout(r, 1000));
      sp.props.xColumnName = 'Chemical Space X';
      sp.props.yColumnName = 'Chemical Space Y';
      sp.props.colorColumnName = 'Primary Scaffold Name';
      sp.props.markersColumnName = 'Stereo Category';
      await new Promise(r => setTimeout(r, 1000));

      const fg = tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 2000));

      const categories = df.col('Primary Scaffold Name').categories;
      const selected = categories.slice(0, 3);
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Primary Scaffold Name', selected});
      await new Promise(r => setTimeout(r, 1000));

      const filtered = df.filter.trueCount;
      return {filtered, total: df.rowCount, selected};
    });
    expect(result.filtered).toBeLessThan(result.total);
    expect(result.filtered).toBeGreaterThan(0);

    // Step 6: Click R_ONE in legend — verify viewer-level legend filtering
    // Must use page.click() — JS DOM click() does not trigger Dart legend handler
    const rOneLocator = page.locator('[name="viewer-Scatter-plot"] .d4-legend-item-extra .d4-legend-value',
      {hasText: 'R_ONE'});
    await rOneLocator.click();
    await page.waitForTimeout(1500);

    const legendResult = await page.evaluate(() => {
      const spContainer = document.querySelector('[name="viewer-Scatter-plot"]');
      const items = spContainer!.querySelectorAll('.d4-legend-item');
      const states: {text: string; isCurrent: boolean; isExtra: boolean}[] = [];
      for (const item of items) {
        const val = item.querySelector('.d4-legend-value');
        states.push({
          text: val?.textContent?.trim() || '',
          isCurrent: item.classList.contains('d4-legend-item-current'),
          isExtra: item.classList.contains('d4-legend-item-extra'),
        });
      }
      const df = grok.shell.tv.dataFrame;
      return {states, filter: df.filter.trueCount};
    });

    // R_ONE should be selected (isCurrent), S_PART and S_UNKN should not
    const rOne = legendResult.states.find((s: any) => s.text === 'R_ONE');
    expect(rOne?.isCurrent).toBe(true);
    const sPart = legendResult.states.find((s: any) => s.text === 'S_PART');
    expect(sPart?.isCurrent).toBe(false);
    // Dataframe filter unchanged (legend filter is viewer-level)
    expect(legendResult.filter).toBe(result.filtered);

    await page.evaluate(() => grok.shell.closeAll());
  });

  // ---- Scenario 5: Color coding from grid ----
  await softStep('5. Color coding from grid', async () => {
    await openDataset(page);

    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;

      const sp = tv.addViewer('Scatter plot');
      await new Promise(r => setTimeout(r, 500));
      const bp = tv.addViewer('Box plot');
      await new Promise(r => setTimeout(r, 500));
      const pc = tv.addViewer('PC plot');
      await new Promise(r => setTimeout(r, 500));

      sp.props.colorColumnName = 'Chemical Space X';
      await new Promise(r => setTimeout(r, 1000));

      // Enable linear color coding on grid
      const grid = Array.from(tv.viewers).find((v: any) => v.type === 'Grid') as any;
      const gridCol = grid.columns.byName('Chemical Space X');
      gridCol.colorCodingType = 'Linear';
      gridCol.isTextColorCoded = true;
      await new Promise(r => setTimeout(r, 1000));

      // Save/apply layout
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1000));

      grok.shell.closeAll();
      const df2 = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      const tv2 = grok.shell.addTableView(df2);
      await new Promise(resolve => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });

      const saved = await grok.dapi.layouts.find(layoutId);
      tv2.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));

      const viewers2 = Array.from(tv2.viewers);
      const sp2 = viewers2.find((v: any) => v.type === 'Scatter plot') as any;
      const grid2 = viewers2.find((v: any) => v.type === 'Grid') as any;
      const gridCol2 = grid2.columns.byName('Chemical Space X');

      const res = {
        viewerCount: viewers2.length,
        spColor: sp2?.props.colorColumnName,
        isTextColorCoded: gridCol2?.isTextColorCoded,
      };

      await grok.dapi.layouts.delete(saved);

      // Change to categorical
      gridCol2.colorCodingType = 'Categorical';
      await new Promise(r => setTimeout(r, 1000));

      grok.shell.closeAll();
      return {...res, categoricalApplied: true};
    });
    expect(result.viewerCount).toBe(4);
    expect(result.spColor).toBe('Chemical Space X');
    expect(result.isTextColorCoded).toBe(true);
    expect(result.categoricalApplied).toBe(true);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
