import {test, expect, chromium} from '@playwright/test';

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
  await page.evaluate(async (path: string) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(() => resolve(undefined), 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
}

test('Scatterplot Legend scenarios', async () => {
  const browser = await chromium.connectOverCDP('http://localhost:9222');
  const context = browser.contexts()[0];
  const page = context.pages()[0] || await context.newPage();

  await page.goto(baseUrl);
  await page.waitForFunction(() => typeof grok !== 'undefined' && grok.shell, {timeout: 15000});

  // ---- Scenario 1: Color/Marker settings and dynamic legends ----
  await softStep('1. Color/Marker settings and dynamic legends', async () => {
    await openDataset(page);

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

    const legendItems = await page.locator('[name="viewer-Scatter-plot"] .d4-legend-item').count();
    expect(legendItems).toBeGreaterThan(0);

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
        setTimeout(() => resolve(undefined), 3000);
      });
      const saved = await grok.dapi.layouts.find(layoutId);
      tv2.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));

      let sp2 = null;
      for (let i = 0; i < tv2.viewers.length; i++) {
        if (tv2.viewers[i].type === 'Scatter plot') { sp2 = tv2.viewers[i]; break; }
      }
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

    const dynamicResult = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      let sp = null;
      for (let i = 0; i < tv.viewers.length; i++) {
        if (tv.viewers[i].type === 'Scatter plot') { sp = tv.viewers[i]; break; }
      }

      df.columns.addNewCalculated('linear_test',
        'if(${Stereo Category}=="S_UNKN", null, ${Average Mass})');
      await new Promise(r => setTimeout(r, 1000));
      sp.props.colorColumnName = 'linear_test';
      await new Promise(r => setTimeout(r, 1000));
      const linearColor = sp.props.colorColumnName;

      df.columns.addNewCalculated('cat_test',
        'if(${Stereo Category}=="S_UNKN", null, ${Series})');
      await new Promise(r => setTimeout(r, 1000));
      sp.props.colorColumnName = 'cat_test';
      await new Promise(r => setTimeout(r, 1000));
      const catColor = sp.props.colorColumnName;

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

      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1000));

      grok.shell.closeAll();
      const df2 = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      const tv2 = grok.shell.addTableView(df2);
      await new Promise(resolve => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(() => resolve(undefined), 3000);
      });
      const saved = await grok.dapi.layouts.find(layoutId);
      tv2.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));

      const sps = [];
      for (let i = 0; i < tv2.viewers.length; i++) {
        if (tv2.viewers[i].type === 'Scatter plot')
          sps.push({filter: tv2.viewers[i].props.filter, marker: tv2.viewers[i].props.markersColumnName});
      }

      await grok.dapi.layouts.delete(saved);
      grok.shell.closeAll();
      return {spCount: sps.length, filters: sps.map(s => s.filter), markers: sps.map(s => s.marker)};
    });
    expect(result.spCount).toBe(2);
    expect(result.filters[0]).toContain('R_ONE');
    expect(result.markers[0]).toBe('Stereo Category');
  });

  // ---- Scenario 4: Filtering via Filter Panel ----
  await softStep('4. Filtering via Filter Panel', async () => {
    await openDataset(page);

    const filterResult = await page.evaluate(async () => {
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

      return {filtered: df.filter.trueCount, total: df.rowCount};
    });
    expect(filterResult.filtered).toBeLessThan(filterResult.total);
    expect(filterResult.filtered).toBeGreaterThan(0);

    // Click R_ONE in legend via pointer events
    const legendResult = await page.evaluate(async () => {
      const spContainer = document.querySelector('[name="viewer-Scatter-plot"]');
      const items = spContainer?.querySelectorAll('.d4-legend-item-extra') || [];
      let target = null;
      for (const item of Array.from(items)) {
        const val = item.querySelector('.d4-legend-value');
        if (val?.textContent?.trim() === 'R_ONE') { target = item; break; }
      }
      if (!target) return {error: 'R_ONE not found'};

      const r = target.getBoundingClientRect();
      const x = r.x + r.width / 2;
      const y = r.y + r.height / 2;
      target.dispatchEvent(new PointerEvent('pointerdown', {clientX: x, clientY: y, bubbles: true, button: 0}));
      target.dispatchEvent(new PointerEvent('pointerup', {clientX: x, clientY: y, bubbles: true, button: 0}));
      target.dispatchEvent(new PointerEvent('click', {clientX: x, clientY: y, bubbles: true, button: 0}));
      await new Promise(r => setTimeout(r, 1500));

      const items2 = spContainer?.querySelectorAll('.d4-legend-item') || [];
      const states: {text: string; isCurrent: boolean}[] = [];
      for (const item of Array.from(items2)) {
        const val = item.querySelector('.d4-legend-value');
        states.push({
          text: val?.textContent?.trim() || '',
          isCurrent: item.classList.contains('d4-legend-item-current'),
        });
      }
      const df = grok.shell.tv.dataFrame;
      return {states, filterCount: df.filter.trueCount};
    });

    const rOne = legendResult.states?.find((s: any) => s.text === 'R_ONE');
    expect(rOne?.isCurrent).toBe(true);
    const sPart = legendResult.states?.find((s: any) => s.text === 'S_PART');
    expect(sPart?.isCurrent).toBe(false);
    expect(legendResult.filterCount).toBe(filterResult.filtered);

    await page.evaluate(() => grok.shell.closeAll());
  });

  // ---- Scenario 5: Color coding from grid ----
  await softStep('5. Color coding from grid', async () => {
    await openDataset(page);

    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;

      const sp = tv.addViewer('Scatter plot');
      await new Promise(r => setTimeout(r, 500));
      tv.addViewer('Box plot');
      await new Promise(r => setTimeout(r, 500));
      tv.addViewer('PC plot');
      await new Promise(r => setTimeout(r, 500));

      sp.props.colorColumnName = 'Chemical Space X';
      await new Promise(r => setTimeout(r, 1000));

      const grid = tv.viewers[0];
      const gridCol = grid.columns.byName('Chemical Space X');
      gridCol.colorCodingType = 'Linear';
      gridCol.isTextColorCoded = true;
      await new Promise(r => setTimeout(r, 1000));

      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1000));

      grok.shell.closeAll();
      const df2 = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      const tv2 = grok.shell.addTableView(df2);
      await new Promise(resolve => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(() => resolve(undefined), 3000);
      });
      const saved = await grok.dapi.layouts.find(layoutId);
      tv2.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));

      let sp2 = null;
      for (let i = 0; i < tv2.viewers.length; i++) {
        if (tv2.viewers[i].type === 'Scatter plot') { sp2 = tv2.viewers[i]; break; }
      }
      const grid2 = tv2.viewers[0];
      const gridCol2 = grid2.columns.byName('Chemical Space X');

      const res = {
        viewerCount: tv2.viewers.length,
        spColor: sp2?.props.colorColumnName,
        isTextColorCoded: gridCol2?.isTextColorCoded,
      };

      await grok.dapi.layouts.delete(saved);

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
