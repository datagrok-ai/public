import {test, expect, type Page} from '@playwright/test';

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/demog.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

// -- UI Helpers --

/** Right-click center of line chart canvas to open context menu */
async function openLineChartContextMenu(page: Page) {
  await page.evaluate(() => {
    document.querySelectorAll('.d4-menu-popup').forEach(m => m.remove());
  });
  await page.waitForTimeout(200);
  const box = await page.evaluate(() => {
    const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart');
    const canvases = lc!.root.querySelectorAll('canvas');
    let mc: HTMLCanvasElement | null = null, ma = 0;
    for (const c of canvases) {
      const r = c.getBoundingClientRect();
      if (r.width * r.height > ma) { ma = r.width * r.height; mc = c; }
    }
    const rect = mc!.getBoundingClientRect();
    return {cx: rect.left + rect.width * 0.4, cy: rect.top + rect.height * 0.4};
  });
  await page.mouse.click(box.cx, box.cy, {button: 'right'});
  await page.waitForTimeout(500);
}

/** Click a context menu item by its name= attribute, making parent container visible */
async function clickNamedMenuItem(page: Page, nameAttr: string) {
  await page.evaluate((name) => {
    const item = document.querySelector(`[name="${name}"]`);
    if (!item) throw new Error(`Menu item not found: ${name}`);
    const container = item.closest('.d4-menu-item-container') as HTMLElement;
    if (container) container.style.display = '';
    item.dispatchEvent(new MouseEvent('click', {bubbles: true}));
  }, nameAttr);
  await page.waitForTimeout(500);
}

/** Open context menu on line chart and click a named submenu item */
async function contextMenuClick(page: Page, nameAttr: string) {
  await openLineChartContextMenu(page);
  await clickNamedMenuItem(page, nameAttr);
}

/** Get the line chart viewer props via evaluate */
async function lcProps(page: Page, ...propNames: string[]): Promise<Record<string, any>> {
  return page.evaluate((names) => {
    const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
    const result: Record<string, any> = {};
    for (const n of names) result[n] = (lc.props as any)[n];
    return result;
  }, propNames);
}

/** Set line chart viewer props via evaluate */
async function lcSetProps(page: Page, props: Record<string, any>) {
  await page.evaluate((p) => {
    const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
    for (const [k, v] of Object.entries(p)) (lc.props as any)[k] = v;
  }, props);
  await page.waitForTimeout(300);
}

/** Open a dataset with full Bio/Chem wait, closeAll first */
async function openDataset(page: Page, path: string) {
  await page.evaluate(async (p) => {
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(p);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
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
  }, path);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 60000});
}

/** Open dataset and add a line chart */
async function openDatasetWithLineChart(page: Page, path: string) {
  await openDataset(page, path);
  await page.evaluate(() => {
    document.querySelector('[name="icon-line-chart"]')!.dispatchEvent(new MouseEvent('click', {bubbles: true}));
  });
  await page.locator('[name="viewer-Line-chart"]').waitFor({timeout: 10000});
}

// -- Setup --

async function setupDemogLineChart(page: Page) {
  // Phase 1: Navigate and wait for full Dart interop
  await page.goto(baseUrl);
  // Wait for the platform UI to be fully loaded (sidebar visible = Dart app initialized)
  await page.locator('[name="Toolbox"], [name="Browse"], .d4-sidebar').first().waitFor({timeout: 60000});
  // Extra wait for Dart interop functions to register
  await page.waitForFunction(() => {
    try {
      if (typeof grok === 'undefined' || !grok.shell) return false;
      grok.shell.settings.showFiltersIconsConstantly;
      return true;
    } catch { return false; }
  }, {timeout: 30000});

  // Phase 2: Open dataset
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
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

  // Phase 3: Add line chart via Toolbox icon
  await page.evaluate(() => {
    document.querySelector('[name="icon-line-chart"]')!.dispatchEvent(new MouseEvent('click', {bubbles: true}));
  });
  await page.locator('[name="viewer-Line-chart"]').waitFor({timeout: 10000});
}

// -- Test --

test('Line chart tests (Playwright) — UI-first', async ({page}) => {
  test.setTimeout(600000);
  stepErrors.length = 0;
  await setupDemogLineChart(page);

  // ── Chart types ────────────────────────────────────────────────────────
  await softStep('Chart types', async () => {
    await contextMenuClick(page, 'div-Chart-Type---Area-Chart');
    expect((await lcProps(page, 'chartTypes')).chartTypes).toContain('Area Chart');

    await contextMenuClick(page, 'div-Chart-Type---Stacked-Area-Chart');
    expect((await lcProps(page, 'chartTypes')).chartTypes).toContain('Stacked Area Chart');

    await contextMenuClick(page, 'div-Chart-Type---Stacked-Bar-Chart');
    expect((await lcProps(page, 'chartTypes')).chartTypes).toContain('Stacked Bar Chart');

    await contextMenuClick(page, 'div-Chart-Type---Line-Chart');
    expect((await lcProps(page, 'chartTypes')).chartTypes).toContain('Line Chart');
  });

  // ── Whiskers (error bars) ──────────────────────────────────────────────
  await softStep('Whiskers (error bars)', async () => {
    await lcSetProps(page, {xColumnName: 'RACE', yColumnNames: ['AGE']});

    await lcSetProps(page, {whiskersType: 'Avg | Min, Max'});
    expect((await lcProps(page, 'whiskersType')).whiskersType).toBe('Avg | Min, Max');

    await lcSetProps(page, {whiskersType: 'Med | Q1, Q3'});
    expect((await lcProps(page, 'whiskersType')).whiskersType).toBe('Med | Q1, Q3');

    await lcSetProps(page, {whiskersType: 'Avg | +/-StDev'});
    expect((await lcProps(page, 'whiskersType')).whiskersType).toBe('Avg | +/-StDev');

    await lcSetProps(page, {whiskersType: 'Avg | +/-StError'});
    expect((await lcProps(page, 'whiskersType')).whiskersType).toBe('Avg | +/-StError');

    await lcSetProps(page, {whiskersType: 'None'});
    expect((await lcProps(page, 'whiskersType')).whiskersType).toBe('None');
  });

  // ── Markers ────────────────────────────────────────────────────────────
  await softStep('Markers', async () => {
    // Step 1: Visibility = Always via context menu
    await contextMenuClick(page, 'div-Markers---Always');
    expect((await lcProps(page, 'showMarkers')).showMarkers).toBe('Always');

    // Steps 2-9: Marker properties via Context Panel (JS API)
    await lcSetProps(page, {markerType: 'square'});
    expect((await lcProps(page, 'markerType')).markerType).toBe('square');

    await lcSetProps(page, {markerType: 'triangle up'});
    expect((await lcProps(page, 'markerType')).markerType).toBe('triangle up');

    await lcSetProps(page, {markerSize: 10});
    expect((await lcProps(page, 'markerSize')).markerSize).toBe(10);

    await lcSetProps(page, {markerOpacity: 50});
    expect((await lcProps(page, 'markerOpacity')).markerOpacity).toBe(50);

    await lcSetProps(page, {markersColumnName: 'SEX'});
    expect((await lcProps(page, 'markersColumnName')).markersColumnName).toBe('SEX');

    await lcSetProps(page, {markersSizeColumnName: 'WEIGHT'});
    expect((await lcProps(page, 'markersSizeColumnName')).markersSizeColumnName).toBe('WEIGHT');

    // Clear marker column
    await lcSetProps(page, {markersColumnName: null});
    expect((await lcProps(page, 'markersColumnName')).markersColumnName).toBeNull();

    // Visibility = Auto via context menu
    await contextMenuClick(page, 'div-Markers---Auto');
    expect((await lcProps(page, 'showMarkers')).showMarkers).toBe('Auto');

    // Reset
    await lcSetProps(page, {markerType: 'circle', markerSize: null, markerOpacity: 100, markersSizeColumnName: null});
  });

  // ── Axis configuration ─────────────────────────────────────────────────
  await softStep('Axis configuration', async () => {
    await lcSetProps(page, {xColumnName: 'AGE', yColumnNames: ['WEIGHT']});

    await lcSetProps(page, {xAxisType: 'logarithmic'});
    expect((await lcProps(page, 'xAxisType')).xAxisType).toBe('logarithmic');

    await lcSetProps(page, {invertXAxis: true});
    expect((await lcProps(page, 'invertXAxis')).invertXAxis).toBe(true);

    await lcSetProps(page, {yAxisType: 'logarithmic'});
    expect((await lcProps(page, 'yAxisType')).yAxisType).toBe('logarithmic');

    await lcSetProps(page, {showVerticalGridLines: false});
    expect((await lcProps(page, 'showVerticalGridLines')).showVerticalGridLines).toBe(false);

    await lcSetProps(page, {showHorizontalGridLines: false});
    expect((await lcProps(page, 'showHorizontalGridLines')).showHorizontalGridLines).toBe(false);

    await lcSetProps(page, {xAxisLabelOrientation: 'Vert'});
    expect((await lcProps(page, 'xAxisLabelOrientation')).xAxisLabelOrientation).toBe('Vert');

    // Reset
    await lcSetProps(page, {
      xAxisType: 'linear', yAxisType: 'linear', invertXAxis: false,
      showVerticalGridLines: true, showHorizontalGridLines: true, xAxisLabelOrientation: 'Auto'
    });
  });

  // ── Interpolation ──────────────────────────────────────────────────────
  await softStep('Interpolation', async () => {
    await lcSetProps(page, {interpolation: 'Spline'});
    expect((await lcProps(page, 'interpolation')).interpolation).toBe('Spline');

    await lcSetProps(page, {splineTension: 1.0});
    expect((await lcProps(page, 'splineTension')).splineTension).toBe(1);

    await lcSetProps(page, {interpolation: 'None'});
    expect((await lcProps(page, 'interpolation')).interpolation).toBe('None');
  });

  // ── Aggregation types ──────────────────────────────────────────────────
  await softStep('Aggregation types', async () => {
    await lcSetProps(page, {xColumnName: 'RACE', yColumnNames: ['AGE']});

    for (const aggr of ['avg', 'min', 'max', 'med', 'sum', 'stdev', 'avg']) {
      await lcSetProps(page, {aggrType: aggr});
      expect((await lcProps(page, 'aggrType')).aggrType).toBe(aggr);
    }
  });

  // ── Left panel histogram ───────────────────────────────────────────────
  await softStep('Left panel histogram', async () => {
    await lcSetProps(page, {leftPanel: 'Histogram'});
    expect((await lcProps(page, 'leftPanel')).leftPanel).toBe('Histogram');

    await lcSetProps(page, {leftPanel: 'None'});
    expect((await lcProps(page, 'leftPanel')).leftPanel).toBe('None');
  });

  // ── Controls visibility ────────────────────────────────────────────────
  await softStep('Controls visibility', async () => {
    const controls = ['showXSelector', 'showYSelectors', 'showAggrTypeSelector', 'showSplitSelector', 'showXAxis', 'showYAxis'];
    for (const ctrl of controls) {
      await lcSetProps(page, {[ctrl]: false});
      expect((await lcProps(page, ctrl))[ctrl]).toBe(false);
    }
    // Recheck all
    const restoreProps: Record<string, boolean> = {};
    for (const ctrl of controls) restoreProps[ctrl] = true;
    await lcSetProps(page, restoreProps);
    for (const ctrl of controls)
      expect((await lcProps(page, ctrl))[ctrl]).toBe(true);
  });

  // ── Y global scale ────────────────────────────────────────────────────
  await softStep('Y global scale', async () => {
    await lcSetProps(page, {yColumnNames: ['AGE', 'HEIGHT']});

    await lcSetProps(page, {multiAxis: true});
    expect((await lcProps(page, 'multiAxis')).multiAxis).toBe(true);

    await lcSetProps(page, {yGlobalScale: true});
    expect((await lcProps(page, 'yGlobalScale')).yGlobalScale).toBe(true);

    await lcSetProps(page, {yGlobalScale: false});
    expect((await lcProps(page, 'yGlobalScale')).yGlobalScale).toBe(false);

    await lcSetProps(page, {multiAxis: false});
  });

  // ── Split by column ────────────────────────────────────────────────────
  await softStep('Split by column', async () => {
    await lcSetProps(page, {splitColumnName: 'SEX'});
    expect((await lcProps(page, 'splitColumnName')).splitColumnName).toBe('SEX');

    await lcSetProps(page, {splitColumnName: 'RACE'});
    expect((await lcProps(page, 'splitColumnName')).splitColumnName).toBe('RACE');

    await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
      lc.props.splitColumnNames = ['SEX', 'RACE'];
    });
    await page.waitForTimeout(300);
    const splitCols = await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
      return lc.props.splitColumnNames;
    });
    expect(splitCols).toEqual(['SEX', 'RACE']);

    await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
      lc.props.splitColumnNames = [];
      lc.props.splitColumnName = '';
    });
    await page.waitForTimeout(300);
  });

  // ── Multi-axis mode ────────────────────────────────────────────────────
  await softStep('Multi-axis mode', async () => {
    await lcSetProps(page, {yColumnNames: ['AGE', 'HEIGHT', 'WEIGHT']});

    await lcSetProps(page, {multiAxis: true});
    expect((await lcProps(page, 'multiAxis')).multiAxis).toBe(true);

    await lcSetProps(page, {multiAxis: false});
    expect((await lcProps(page, 'multiAxis')).multiAxis).toBe(false);
  });

  // ── Title and description ──────────────────────────────────────────────
  await softStep('Title and description', async () => {
    await lcSetProps(page, {showTitle: true});
    expect((await lcProps(page, 'showTitle')).showTitle).toBe(true);

    await lcSetProps(page, {title: 'My Line Chart'});
    expect((await lcProps(page, 'title')).title).toBe('My Line Chart');

    await lcSetProps(page, {description: 'Test description'});
    expect((await lcProps(page, 'description')).description).toBe('Test description');

    await lcSetProps(page, {descriptionPosition: 'Top'});
    expect((await lcProps(page, 'descriptionPosition')).descriptionPosition).toBe('Top');

    await lcSetProps(page, {descriptionPosition: 'Bottom'});
    expect((await lcProps(page, 'descriptionPosition')).descriptionPosition).toBe('Bottom');

    await lcSetProps(page, {descriptionVisibilityMode: 'Never'});
    expect((await lcProps(page, 'descriptionVisibilityMode')).descriptionVisibilityMode).toBe('Never');

    // Reset
    await lcSetProps(page, {showTitle: false, title: '', description: ''});
  });

  // ── Custom axis min/max ────────────────────────────────────────────────
  await softStep('Custom axis min/max', async () => {
    await lcSetProps(page, {xColumnName: 'AGE', yColumnNames: ['HEIGHT']});

    await lcSetProps(page, {xMin: 30});
    expect((await lcProps(page, 'xMin')).xMin).toBe(30);

    await lcSetProps(page, {xMax: 60});
    expect((await lcProps(page, 'xMax')).xMax).toBe(60);

    await lcSetProps(page, {yMin: 150});
    expect((await lcProps(page, 'yMin')).yMin).toBe(150);

    await lcSetProps(page, {yMax: 190});
    expect((await lcProps(page, 'yMax')).yMax).toBe(190);

    // Clear
    await lcSetProps(page, {xMin: null, xMax: null, yMin: null, yMax: null});
    const cleared = await lcProps(page, 'xMin', 'xMax', 'yMin', 'yMax');
    expect(cleared.xMin).toBeNull();
    expect(cleared.yMax).toBeNull();
  });

  // ── Date/time X axis ───────────────────────────────────────────────────
  await softStep('Date/time X axis', async () => {
    await lcSetProps(page, {xColumnName: 'STARTED'});
    expect((await lcProps(page, 'xColumnName')).xColumnName).toBe('STARTED');

    await lcSetProps(page, {xMap: 'Year'});
    expect((await lcProps(page, 'xMap')).xMap).toBe('Year');

    await lcSetProps(page, {xMap: 'Month'});
    expect((await lcProps(page, 'xMap')).xMap).toBe('Month');

    await lcSetProps(page, {xMap: 'Day of week'});
    expect((await lcProps(page, 'xMap')).xMap).toBe('Day of week');

    await lcSetProps(page, {xMap: 'None'});
    expect((await lcProps(page, 'xMap')).xMap).toBe('None');
  });

  // ── Line styling ───────────────────────────────────────────────────────
  await softStep('Line styling', async () => {
    await lcSetProps(page, {xColumnName: 'AGE', lineWidth: 3});
    expect((await lcProps(page, 'lineWidth')).lineWidth).toBe(3);

    await lcSetProps(page, {lineTransparency: 0.5});
    expect((await lcProps(page, 'lineTransparency')).lineTransparency).toBe(0.5);

    await lcSetProps(page, {lineColoringType: 'Custom'});
    expect((await lcProps(page, 'lineColoringType')).lineColoringType).toBe('Custom');

    await lcSetProps(page, {lineWidth: 1, lineTransparency: 0, lineColoringType: 'Auto'});
  });

  // ── Axis tickmarks modes ───────────────────────────────────────────────
  await softStep('Axis tickmarks modes', async () => {
    await lcSetProps(page, {xColumnName: 'AGE', yColumnNames: ['HEIGHT']});

    await lcSetProps(page, {xAxisTickmarksMode: 'MinMax'});
    expect((await lcProps(page, 'xAxisTickmarksMode')).xAxisTickmarksMode).toBe('MinMax');

    await lcSetProps(page, {xAxisTickmarksMode: 'Auto'});
    expect((await lcProps(page, 'xAxisTickmarksMode')).xAxisTickmarksMode).toBe('Auto');

    await lcSetProps(page, {yAxisTickmarksMode: 'MinMax'});
    expect((await lcProps(page, 'yAxisTickmarksMode')).yAxisTickmarksMode).toBe('MinMax');

    await lcSetProps(page, {yAxisTickmarksMode: 'Auto'});
  });

  // ── Overview chart ─────────────────────────────────────────────────────
  await softStep('Overview chart', async () => {
    await contextMenuClick(page, 'div-Overview---Line-Chart');
    expect((await lcProps(page, 'overviewType')).overviewType).toBe('Line Chart');

    await contextMenuClick(page, 'div-Overview---Area-Chart');
    expect((await lcProps(page, 'overviewType')).overviewType).toBe('Area Chart');

    await contextMenuClick(page, 'div-Overview---Stacked-Bar-Chart');
    expect((await lcProps(page, 'overviewType')).overviewType).toBe('Stacked Bar Chart');

    await contextMenuClick(page, 'div-Overview---None');
    expect((await lcProps(page, 'overviewType')).overviewType).toBe('None');
  });

  // ── Legend ──────────────────────────────────────────────────────────────
  await softStep('Legend', async () => {
    await lcSetProps(page, {splitColumnName: 'SEX'});

    await lcSetProps(page, {legendVisibility: 'Always'});
    expect((await lcProps(page, 'legendVisibility')).legendVisibility).toBe('Always');

    for (const pos of ['Left', 'Top', 'Bottom', 'Right']) {
      await lcSetProps(page, {legendPosition: pos});
      expect((await lcProps(page, 'legendPosition')).legendPosition).toBe(pos);
    }

    await lcSetProps(page, {legendVisibility: 'Never'});
    expect((await lcProps(page, 'legendVisibility')).legendVisibility).toBe('Never');

    await lcSetProps(page, {legendVisibility: 'Auto'});
    expect((await lcProps(page, 'legendVisibility')).legendVisibility).toBe('Auto');

    await lcSetProps(page, {splitColumnName: ''});
  });

  // ── Axes follow filter ─────────────────────────────────────────────────
  await softStep('Axes follow filter', async () => {
    await lcSetProps(page, {xColumnName: 'AGE'});

    const defaultVal = (await lcProps(page, 'axesFollowFilter')).axesFollowFilter;
    expect(defaultVal).toBe(true);

    await lcSetProps(page, {axesFollowFilter: false});
    expect((await lcProps(page, 'axesFollowFilter')).axesFollowFilter).toBe(false);

    // Open filter panel and narrow AGE range
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 30, max: 60});
    });
    await page.waitForTimeout(500);

    const filteredCount = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filteredCount).toBeLessThan(5850);

    await lcSetProps(page, {axesFollowFilter: true});
    expect((await lcProps(page, 'axesFollowFilter')).axesFollowFilter).toBe(true);

    // Reset filter
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 0, max: 100});
    });
    await page.waitForTimeout(300);
  });

  // ── Context menu — chart area ──────────────────────────────────────────
  await softStep('Context menu — chart area', async () => {
    await openLineChartContextMenu(page);
    const items = await page.evaluate(() => {
      const popup = document.querySelector('.d4-menu-popup');
      if (!popup) return [];
      return Array.from(popup.querySelectorAll(':scope > .d4-menu-item > .d4-menu-item-label'))
        .map(el => el.textContent!.trim());
    });

    for (const expected of ['Reset View', 'Tools', 'Data', 'Markers', 'Chart Type', 'Overview', 'Selection', 'Controls'])
      expect(items).toContain(expected);

    await page.keyboard.press('Escape');
  });

  // ── Layout save and restore ────────────────────────────────────────────
  await softStep('Layout save and restore', async () => {
    // Configure
    await lcSetProps(page, {
      xColumnName: 'STARTED', yColumnNames: ['AGE', 'HEIGHT'],
      splitColumnName: 'SEX', multiAxis: true, lineWidth: 3, interpolation: 'Spline'
    });

    // Save layout
    const layoutId = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      return layout.id;
    });

    // Close line chart
    await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart');
      lc!.close();
    });
    await page.waitForTimeout(500);

    // Apply saved layout
    await page.evaluate(async (id) => {
      const saved = await grok.dapi.layouts.find(id);
      grok.shell.tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));
    }, layoutId);

    // Verify
    const restored = await lcProps(page, 'xColumnName', 'yColumnNames', 'splitColumnName', 'multiAxis', 'lineWidth', 'interpolation');
    expect(restored.xColumnName).toBe('STARTED');
    expect(restored.yColumnNames).toEqual(['AGE', 'HEIGHT']);
    expect(restored.splitColumnName).toBe('SEX');
    expect(restored.multiAxis).toBe(true);
    expect(restored.lineWidth).toBe(3);
    expect(restored.interpolation).toBe('Spline');

    // Delete layout
    await page.evaluate(async (id) => {
      const saved = await grok.dapi.layouts.find(id);
      await grok.dapi.layouts.delete(saved);
    }, layoutId);
  });

  // ── Table switching and row source ─────────────────────────────────────
  await softStep('Table switching and row source', async () => {
    // Close existing line charts from previous steps, then open fresh
    await page.evaluate(() => {
      for (const v of Array.from(grok.shell.tv.viewers))
        if (v.type === 'Line chart') v.close();
    });
    await page.waitForTimeout(300);

    // Open SPGI (with Bio/Chem wait)
    await page.evaluate(async () => {
      const spgi = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      grok.shell.addTableView(spgi);
      await new Promise(resolve => {
        const sub = spgi.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });
      const hasBioChem = Array.from({length: spgi.columns.length}, (_, i) => spgi.columns.byIndex(i))
        .some(c => c.semType === 'Molecule' || c.semType === 'Macromolecule');
      if (hasBioChem) {
        for (let i = 0; i < 50; i++) {
          if (document.querySelector('[name="viewer-Grid"] canvas')) break;
          await new Promise(r => setTimeout(r, 200));
        }
        await new Promise(r => setTimeout(r, 5000));
      }
    });

    // Switch back to demog view
    await page.evaluate(() => {
      const demogTv = Array.from(grok.shell.tableViews).find(v => v.dataFrame.rowCount === 5850);
      if (demogTv) grok.shell.v = demogTv;
    });
    await page.waitForTimeout(500);

    // Add line chart on demog
    await page.evaluate(() => {
      document.querySelector('[name="icon-line-chart"]')!.dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.locator('[name="viewer-Line-chart"]').first().waitFor({timeout: 10000});

    // Switch table to SPGI
    const spgiTableName = await page.evaluate(() => {
      const spgiTv = Array.from(grok.shell.tableViews).find(v => v.dataFrame.rowCount === 3624);
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
      lc.props.table = spgiTv!.dataFrame.name;
      return lc.props.table;
    });
    expect(spgiTableName).toBeTruthy();

    // Switch back to demog
    await page.evaluate(() => {
      const demogTv = Array.from(grok.shell.tableViews).find(v => v.dataFrame.rowCount === 5850);
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
      lc.props.table = demogTv!.dataFrame.name;
    });
    await page.waitForTimeout(300);

    // Row Source = Selected
    await lcSetProps(page, {rowSource: 'Selected'});
    expect((await lcProps(page, 'rowSource')).rowSource).toBe('Selected');

    // Select some rows
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      for (let i = 0; i < 100; i++) df.selection.set(i, true);
    });
    await page.waitForTimeout(300);

    // Row Source = Filtered
    await lcSetProps(page, {rowSource: 'Filtered'});

    // Apply filter
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 20, max: 40});
    });
    await page.waitForTimeout(500);

    const filtered = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filtered).toBeLessThan(5850);

    // Reset
    await page.evaluate(() => {
      grok.shell.tv.dataFrame.selection.setAll(false);
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 0, max: 100});
    });
    await page.waitForTimeout(300);

    // Close SPGI view
    await page.evaluate(() => {
      const spgiTv = Array.from(grok.shell.tableViews).find(v => v.dataFrame.rowCount === 3624);
      spgiTv?.close();
    });
  });

  // ── Filter expression and collaborative filtering (SPGI) ─────────────
  await softStep('Filter expression and collaborative filtering (SPGI)', async () => {
    await openDatasetWithLineChart(page, 'System:DemoFiles/SPGI.csv');

    await lcSetProps(page, {filter: '${CAST Idea ID} <636500'});
    const filterVal = (await lcProps(page, 'filter')).filter;
    expect(filterVal).toBe('${CAST Idea ID} <636500');

    await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 500));
      fg.updateOrAdd({type: 'histogram', column: 'Average Mass', min: 200, max: 400});
    });
    await page.waitForTimeout(500);

    const filtered = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filtered).toBeLessThan(3624);

    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'Average Mass', min: -Infinity, max: Infinity});
    });
    await page.waitForTimeout(300);
    await lcSetProps(page, {filter: ''});
  });

  // ── Split and Y-columns sync with Context Panel ──────────────────────
  await softStep('Split and Y-columns sync with Context Panel', async () => {
    await openDatasetWithLineChart(page, 'System:DemoFiles/demog.csv');

    await lcSetProps(page, {splitColumnName: 'SEX'});
    expect((await lcProps(page, 'splitColumnName')).splitColumnName).toBe('SEX');

    await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
      lc.props.splitColumnNames = ['SEX', 'RACE'];
    });
    await page.waitForTimeout(300);
    const splitCols = await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
      return lc.props.splitColumnNames;
    });
    expect(splitCols).toEqual(['SEX', 'RACE']);

    await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
      lc.props.splitColumnNames = ['SEX'];
    });
    await page.waitForTimeout(300);

    await lcSetProps(page, {yColumnNames: ['AGE', 'HEIGHT', 'WEIGHT']});
    expect((await lcProps(page, 'yColumnNames')).yColumnNames).toHaveLength(3);

    await lcSetProps(page, {yColumnNames: ['AGE', 'WEIGHT']});
    expect((await lcProps(page, 'yColumnNames')).yColumnNames).toHaveLength(2);

    await lcSetProps(page, {yColumnNames: ['AGE']});
    expect((await lcProps(page, 'yColumnNames')).yColumnNames).toHaveLength(1);

    await lcSetProps(page, {yColumnNames: ['AGE', 'HEIGHT'], multiAxis: true});
    expect((await lcProps(page, 'multiAxis')).multiAxis).toBe(true);

    await lcSetProps(page, {multiAxis: false});
    await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
      lc.props.splitColumnNames = [];
      lc.props.splitColumnName = '';
    });
  });

  // ── Selection checkboxes ──────────────────────────────────────────────
  await softStep('Selection checkboxes', async () => {
    await lcSetProps(page, {xColumnName: 'AGE', yColumnNames: ['HEIGHT']});

    const box = await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
      const canvases = lc.root.querySelectorAll('canvas');
      let mc: HTMLCanvasElement | null = null, ma = 0;
      for (const c of canvases) {
        const r = c.getBoundingClientRect();
        if (r.width * r.height > ma) { ma = r.width * r.height; mc = c; }
      }
      const rect = mc!.getBoundingClientRect();
      return {x1: rect.left + rect.width * 0.3, y1: rect.top + rect.height * 0.3,
        x2: rect.left + rect.width * 0.6, y2: rect.top + rect.height * 0.6};
    });

    await page.keyboard.down('Shift');
    await page.mouse.move(box.x1, box.y1);
    await page.mouse.down();
    await page.mouse.move(box.x2, box.y2, {steps: 10});
    await page.mouse.up();
    await page.keyboard.up('Shift');
    await page.waitForTimeout(500);

    const sel = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    console.log(`Selection checkboxes: ${sel} rows selected`);

    // Layout save/restore with selection
    const layoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      return layout.id;
    });

    await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart');
      lc!.close();
    });
    await page.waitForTimeout(500);

    await page.evaluate(async (id) => {
      const saved = await grok.dapi.layouts.find(id);
      grok.shell.tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));
    }, layoutId);

    const lcRestored = await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).some(v => v.type === 'Line chart'));
    expect(lcRestored).toBe(true);

    await page.evaluate(async (id) => {
      const saved = await grok.dapi.layouts.find(id);
      await grok.dapi.layouts.delete(saved);
    }, layoutId);
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));
  });

  // ── Data panel checkboxes ─────────────────────────────────────────────
  await softStep('Data panel checkboxes', async () => {
    const layoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      return layout.id;
    });

    await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart');
      lc!.close();
    });
    await page.waitForTimeout(500);

    await page.evaluate(async (id) => {
      const saved = await grok.dapi.layouts.find(id);
      grok.shell.tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));
    }, layoutId);

    const lcRestored = await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).some(v => v.type === 'Line chart'));
    expect(lcRestored).toBe(true);

    await page.evaluate(async (id) => {
      const saved = await grok.dapi.layouts.find(id);
      await grok.dapi.layouts.delete(saved);
    }, layoutId);
  });

  // ── GROK-17835 regression ──────────────────────────────────────────────
  await softStep('GROK-17835 regression (SPGI)', async () => {
    await openDatasetWithLineChart(page, 'System:DemoFiles/SPGI.csv');

    // Multi Axis
    await lcSetProps(page, {multiAxis: true});

    // Split by Series and Scaffold Names
    await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
      lc.props.splitColumnNames = ['Series', 'Scaffold Names'];
    });
    await page.waitForTimeout(500);

    // Set X to Chemist 521
    await lcSetProps(page, {xColumnName: 'Chemist 521'});

    // Hover over viewer — should not cause errors
    const warningsBefore = await page.evaluate(() => grok.shell.warnings?.length ?? 0);

    const hoverBox = await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
      const canvases = lc.root.querySelectorAll('canvas');
      let mc: HTMLCanvasElement | null = null, ma = 0;
      for (const c of canvases) {
        const r = c.getBoundingClientRect();
        if (r.width * r.height > ma) { ma = r.width * r.height; mc = c; }
      }
      const rect = mc!.getBoundingClientRect();
      return {cx: rect.left + rect.width * 0.5, cy: rect.top + rect.height * 0.5};
    });

    await page.mouse.move(hoverBox.cx, hoverBox.cy);
    await page.waitForTimeout(1000);

    const warningsAfter = await page.evaluate(() => grok.shell.warnings?.length ?? 0);
    expect(warningsAfter).toBe(warningsBefore);
  });

  // ── Final summary ──────────────────────────────────────────────────────
  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
