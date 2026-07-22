import {test, expect, type Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

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

async function contextMenuClick(page: Page, nameAttr: string) {
  await openLineChartContextMenu(page);
  await clickNamedMenuItem(page, nameAttr);
}

async function lcProps(page: Page, ...propNames: string[]): Promise<Record<string, any>> {
  return page.evaluate((names) => {
    const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
    const result: Record<string, any> = {};
    for (const n of names) result[n] = (lc.props as any)[n];
    return result;
  }, propNames);
}

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

async function openDatasetWithLineChart(page: Page, path: string) {
  await openDataset(page, path);
  await page.evaluate(() => {
    document.querySelector('[name="icon-line-chart"]')!.dispatchEvent(new MouseEvent('click', {bubbles: true}));
  });
  await page.locator('[name="viewer-Line-chart"]').waitFor({timeout: 10000});
}

async function setupDemogLineChart(page: Page) {
  await loginToDatagrok(page);

  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await page.evaluate(() => {
    document.querySelector('[name="icon-line-chart"]')!.dispatchEvent(new MouseEvent('click', {bubbles: true}));
  });
  await page.locator('[name="viewer-Line-chart"]').waitFor({timeout: 10000});
}

test('Line chart tests (Playwright) — UI-first', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  // Settings that only change how the chart is painted have no DOM counterpart, so
  // the check is that driving them raises nothing and leaves the chart rendering.
  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });
  const errorCount = () => pageErrors.length + consoleErrors.length;
  const chartAlive = () => page.evaluate(() => {
    const lc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Line chart');
    const root = document.querySelector('[name="viewer-Line-chart"]');
    return !!lc && !!root && root.querySelectorAll('canvas').length > 0;
  });

  await setupDemogLineChart(page);

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

  await softStep('Axis configuration', async () => {
    const errBefore = errorCount();
    await lcSetProps(page, {xColumnName: 'AGE', yColumnNames: ['WEIGHT']});

    await lcSetProps(page, {xAxisType: 'logarithmic'});

    await lcSetProps(page, {invertXAxis: true});

    await lcSetProps(page, {yAxisType: 'logarithmic'});

    await lcSetProps(page, {showVerticalGridLines: false});

    await lcSetProps(page, {showHorizontalGridLines: false});

    await lcSetProps(page, {xAxisLabelOrientation: 'Vert'});

    await lcSetProps(page, {
      xAxisType: 'linear', yAxisType: 'linear', invertXAxis: false,
      showVerticalGridLines: true, showHorizontalGridLines: true, xAxisLabelOrientation: 'Auto'
    });
    expect(await chartAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);
  });

  await softStep('Interpolation', async () => {
    const errBefore = errorCount();
    await lcSetProps(page, {interpolation: 'Spline'});

    await lcSetProps(page, {splineTension: 1.0});

    await lcSetProps(page, {interpolation: 'None'});
    expect(await chartAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);
  });

  await softStep('Left panel histogram', async () => {
    const errBefore = errorCount();
    await lcSetProps(page, {leftPanel: 'Histogram'});

    await lcSetProps(page, {leftPanel: 'None'});
    expect(await chartAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);
  });

  await softStep('Controls visibility', async () => {
    const errBefore = errorCount();
    const controls = ['showXSelector', 'showYSelectors', 'showAggrTypeSelector', 'showSplitSelector', 'showXAxis', 'showYAxis'];
    for (const ctrl of controls) {
      await lcSetProps(page, {[ctrl]: false});
    }
    const restoreProps: Record<string, boolean> = {};
    for (const ctrl of controls) restoreProps[ctrl] = true;
    await lcSetProps(page, restoreProps);
    expect(await chartAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);
  });

  await softStep('Y global scale', async () => {
    const errBefore = errorCount();
    await lcSetProps(page, {yColumnNames: ['AGE', 'HEIGHT']});

    await lcSetProps(page, {multiAxis: true});

    await lcSetProps(page, {yGlobalScale: true});

    await lcSetProps(page, {yGlobalScale: false});

    await lcSetProps(page, {multiAxis: false});
    expect(await chartAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);
  });

  await softStep('Title and description', async () => {
    const errBefore = errorCount();
    await lcSetProps(page, {showTitle: true});

    await lcSetProps(page, {title: 'My Line Chart'});

    await lcSetProps(page, {description: 'Test description'});

    await lcSetProps(page, {descriptionPosition: 'Top'});

    await lcSetProps(page, {descriptionPosition: 'Bottom'});

    await lcSetProps(page, {descriptionVisibilityMode: 'Never'});

    await lcSetProps(page, {showTitle: false, title: '', description: ''});
    expect(await chartAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);
  });

  await softStep('Date/time X axis', async () => {
    const errBefore = errorCount();
    await lcSetProps(page, {xColumnName: 'STARTED'});

    await lcSetProps(page, {xMap: 'year'});

    await lcSetProps(page, {xMap: 'month'});

    await lcSetProps(page, {xMap: 'day'});

    await lcSetProps(page, {xMap: ''});
    expect(await chartAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);
  });

  await softStep('Line styling', async () => {
    const errBefore = errorCount();
    await lcSetProps(page, {xColumnName: 'AGE', lineWidth: 3});

    await lcSetProps(page, {lineTransparency: 0.5});

    await lcSetProps(page, {lineColoringType: 'Custom'});

    await lcSetProps(page, {lineWidth: 1, lineTransparency: 0, lineColoringType: 'Auto'});
    expect(await chartAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);
  });

  await softStep('Axis tickmarks modes', async () => {
    const errBefore = errorCount();
    await lcSetProps(page, {xColumnName: 'AGE', yColumnNames: ['HEIGHT']});

    await lcSetProps(page, {xAxisTickmarksMode: 'MinMax'});

    await lcSetProps(page, {xAxisTickmarksMode: 'Auto'});

    await lcSetProps(page, {yAxisTickmarksMode: 'MinMax'});

    await lcSetProps(page, {yAxisTickmarksMode: 'Auto'});
    expect(await chartAlive()).toBe(true);
    expect(errorCount()).toBe(errBefore);
  });

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

  // The legend is a real DOM element whose entries are the split column's
  // categories, so visibility is checked by looking for it. Position is a
  // canvas-layout outcome and is only driven here.
  await softStep('Legend', async () => {
    const legend = () => page.evaluate(() => {
      const el = document.querySelector('[name="viewer-Line-chart"] .d4-legend') as HTMLElement | null;
      return {
        present: !!el,
        labels: el ? (el.innerText || '').split('\n').map(s => s.trim()).filter(Boolean) : [],
      };
    });

    // The legend is re-rendered asynchronously, so wait for the state rather than
    // for a fixed delay.
    const legendSettles = async (present: boolean) =>
      expect.poll(async () => (await legend()).present, {timeout: 10_000}).toBe(present);

    // Earlier steps leave colouring and axis settings behind, several of which put
    // a legend on the chart, so start from a known baseline.
    await lcSetProps(page, {
      splitColumnName: '', legendVisibility: 'Auto', legendPosition: 'Auto',
      lineColoringType: 'Auto', multiAxis: false, yColumnNames: ['AGE'],
    });
    await legendSettles(false);

    await lcSetProps(page, {splitColumnName: 'SEX', legendVisibility: 'Always'});
    await legendSettles(true);
    // demog's SEX column has exactly two values, and they are what the legend lists.
    expect((await legend()).labels.sort()).toEqual(['F', 'M']);

    for (const pos of ['Left', 'Top', 'Bottom', 'Right'])
      await lcSetProps(page, {legendPosition: pos});
    await legendSettles(true);

    await lcSetProps(page, {legendVisibility: 'Never'});
    await legendSettles(false);

    await lcSetProps(page, {legendVisibility: 'Auto'});
    await legendSettles(true);
    expect((await legend()).labels.sort()).toEqual(['F', 'M']);

    // Clearing splitColumnName alone leaves the legend on screen until something
    // triggers a repaint, so the second assignment is what retires it.
    await lcSetProps(page, {splitColumnName: ''});
    await lcSetProps(page, {splitColumnNames: []});
    await legendSettles(false);
  });

  await softStep('Axes follow filter', async () => {
    await lcSetProps(page, {xColumnName: 'AGE'});

    const defaultVal = (await lcProps(page, 'axesFollowFilter')).axesFollowFilter;
    expect(defaultVal).toBe(true);

    await lcSetProps(page, {axesFollowFilter: false});
    expect((await lcProps(page, 'axesFollowFilter')).axesFollowFilter).toBe(false);

    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 30, max: 60});
    });
    await page.waitForTimeout(500);

    const filteredCount = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filteredCount).toBeLessThan(5850);

    await lcSetProps(page, {axesFollowFilter: true});
    expect((await lcProps(page, 'axesFollowFilter')).axesFollowFilter).toBe(true);

    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 0, max: 100});
    });
    await page.waitForTimeout(300);
  });

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

  await softStep('Layout save and restore', async () => {
    await lcSetProps(page, {
      xColumnName: 'STARTED', yColumnNames: ['AGE', 'HEIGHT'],
      splitColumnName: 'SEX', multiAxis: true, lineWidth: 3, interpolation: 'Spline'
    });

    const layoutId = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const layout = tv.saveLayout();
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

    const restored = await lcProps(page, 'xColumnName', 'yColumnNames', 'splitColumnName', 'multiAxis', 'lineWidth', 'interpolation');
    expect(restored.xColumnName).toBe('STARTED');
    expect(restored.yColumnNames).toEqual(['AGE', 'HEIGHT']);
    expect(restored.splitColumnName).toBe('SEX');
    expect(restored.multiAxis).toBe(true);
    expect(restored.lineWidth).toBe(3);
    expect(restored.interpolation).toBe('Spline');

    await page.evaluate(async (id) => {
      const saved = await grok.dapi.layouts.find(id);
      await grok.dapi.layouts.delete(saved);
    }, layoutId);
  });

  await softStep('Table switching and row source', async () => {
    await page.evaluate(() => {
      for (const v of Array.from(grok.shell.tv.viewers))
        if (v.type === 'Line chart') v.close();
    });
    await page.waitForTimeout(300);

    // SPGI has Bio/Chem columns, so the grid needs the extra render wait.
    await page.evaluate(async () => {
      const spgi = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
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

    await page.evaluate(() => {
      const demogTv = Array.from(grok.shell.tableViews).find(v => v.dataFrame.rowCount === 5850);
      if (demogTv) grok.shell.v = demogTv;
    });
    await page.waitForTimeout(500);

    await page.evaluate(() => {
      document.querySelector('[name="icon-line-chart"]')!.dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page.locator('[name="viewer-Line-chart"]').first().waitFor({timeout: 10000});

    // Bind the chart to SPGI: props.table echoes the SPGI name and the viewer's
    // bound dataFrame switches to SPGI's, so the row source moved between tables,
    // not just the stored string.
    const toSpgi = await page.evaluate(() => {
      const spgiTv = Array.from(grok.shell.tableViews).find(v => v.dataFrame.rowCount === 100);
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
      lc.props.table = spgiTv!.dataFrame.name;
      return {name: lc.props.table, wanted: spgiTv!.dataFrame.name};
    });
    await page.waitForTimeout(500);
    expect(toSpgi.name).toBe(toSpgi.wanted);
    expect(await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
      return (lc as any).dataFrame.rowCount;
    })).toBe(100);

    // Switch back to demog: the bound dataFrame returns to demog's row count.
    await page.evaluate(() => {
      const demogTv = Array.from(grok.shell.tableViews).find(v => v.dataFrame.rowCount === 5850);
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
      lc.props.table = demogTv!.dataFrame.name;
    });
    await page.waitForTimeout(500);
    expect(await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
      return (lc as any).dataFrame.rowCount;
    })).toBe(5850);

    // Row source = Selected: selecting rows re-renders the chart to the selected
    // subset. Under rowSource=Selected an empty selection still paints, so the
    // signal is that selecting CHANGES the canvas (per-color diff), not a
    // blank-vs-drawn threshold. Drain the render tail with a settle-precheck diff
    // before the measured one.
    await lcSetProps(page, {xColumnName: 'AGE', yColumnNames: ['HEIGHT'], rowSource: 'Selected'});
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));
    await page.waitForTimeout(800);
    await v.snapshotCanvasColors(page, 'Line chart');
    await page.waitForTimeout(700);
    // Settle-precheck: prove the rowSource-flip repaint has drained before the
    // baseline snapshot, so the measured delta is the selection's effect only.
    const preSel = (await v.diffCanvasColors(page, 'Line chart')).deltaPx;
    expect(preSel).toBeGreaterThanOrEqual(0);
    expect(preSel).toBeLessThan(200);
    await v.snapshotCanvasColors(page, 'Line chart');
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      for (let i = 0; i < 100; i++) df.selection.set(i, true);
    });
    await page.waitForTimeout(1000);
    const selDelta = (await v.diffCanvasColors(page, 'Line chart')).deltaPx;
    expect(selDelta).toBeGreaterThan(1000);

    await lcSetProps(page, {rowSource: 'Filtered'});

    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 20, max: 40});
    });
    await page.waitForTimeout(500);

    const filtered = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filtered).toBeLessThan(5850);

    await page.evaluate(() => {
      grok.shell.tv.dataFrame.selection.setAll(false);
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 0, max: 100});
    });
    await page.waitForTimeout(300);
    await lcSetProps(page, {rowSource: 'All'});

    await page.evaluate(() => {
      const spgiTv = Array.from(grok.shell.tableViews).find(v => v.dataFrame.rowCount === 100);
      spgiTv?.close();
    });
  });

  await softStep('Filter expression and collaborative filtering', async () => {
    await openDatasetWithLineChart(page, 'System:AppData/Chem/tests/spgi-100.csv');

    const full = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);

    // The expression filter is a viewer-local predicate: it narrows the rows the
    // line chart plots (lc.filter) without touching the shared df.filter. The
    // real signal is that in-viewer count dropping below the full row set — the
    // points on the chart are filtered — not the string echoing back.
    await lcSetProps(page, {filter: '${CAST Idea ID} <634834'});
    const lcFiltered = await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
      return (lc as any).filter.trueCount;
    });
    expect(lcFiltered).toBeGreaterThan(0);
    expect(lcFiltered).toBeLessThan(full);

    // Collaborative filtering: a Filter Panel numeric filter narrows the shared
    // df.filter on top of the in-viewer expression filter.
    await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 500));
      fg.updateOrAdd({type: 'histogram', column: 'Average Mass', min: 200, max: 400});
    });
    await page.waitForTimeout(500);

    const filtered = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filtered).toBeLessThan(full);

    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'Average Mass', min: -Infinity, max: Infinity});
    });
    await page.waitForTimeout(300);
    await lcSetProps(page, {filter: ''});
  });

  await softStep('Selection checkboxes', async () => {
    // The previous step left SPGI open; this step needs demog's AGE/HEIGHT
    // columns, so reopen demog with a fresh line chart.
    await openDatasetWithLineChart(page, datasetPath);
    const errBefore = errorCount();
    // Single-line chart with no split column: no categorical-palette orange is
    // painted, so the only orange on the canvas is the selection overlay
    // (SELECTION_HUE_RANGE), which makes the hue count a clean selection signal.
    await lcSetProps(page, {xColumnName: 'AGE', yColumnNames: ['HEIGHT'], splitColumnName: '',
      rowSource: 'All', showSelectedRows: true});

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
    await page.waitForTimeout(700);

    // Shift-drag selects rows in the dataframe...
    const sel = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(sel).toBeGreaterThan(0);
    // ...and the line chart paints the selection overlay for them.
    const hueSelected = await v.countSelectionHuePixels(page, 'Line chart');
    expect(hueSelected).toBeGreaterThan(0);

    // Show Selected Rows checkbox: unchecking it removes the overlay, rechecking
    // restores it — the checkbox's real visual effect, not just its stored value.
    await lcSetProps(page, {showSelectedRows: false});
    const hueOff = await v.countSelectionHuePixels(page, 'Line chart');
    // Near-zero: unchecking must remove the overlay, not just halve it.
    expect(hueOff).toBeLessThan(hueSelected * 0.15);
    await lcSetProps(page, {showSelectedRows: true});
    expect(await v.countSelectionHuePixels(page, 'Line chart')).toBeGreaterThan(0);

    // The remaining Selection checkboxes drive mouse-over highlights that need a
    // live pointer path (no headless canvas signal), so exercise each through its
    // stored value and guard that toggling raises nothing.
    for (const p of ['showCurrentRowLine', 'showMouseOverCategory', 'showMouseOverRowLine']) {
      const cur = (await lcProps(page, p))[p];
      await lcSetProps(page, {[p]: !cur});
      expect((await lcProps(page, p))[p]).toBe(!cur);
    }
    expect(errorCount()).toBe(errBefore);

    // Persist a non-default Selection combo and confirm it survives a layout
    // save → close → restore round-trip.
    const combo = {showCurrentRowLine: true, showMouseOverCategory: false,
      showSelectedRows: false, showMouseOverRowLine: false};
    await lcSetProps(page, combo);

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

    const restored = await lcProps(page, 'showCurrentRowLine', 'showMouseOverCategory',
      'showSelectedRows', 'showMouseOverRowLine');
    expect(restored).toEqual(combo);

    await page.evaluate(async (id) => {
      const saved = await grok.dapi.layouts.find(id);
      await grok.dapi.layouts.delete(saved);
    }, layoutId);
    await lcSetProps(page, {showCurrentRowLine: false, showMouseOverCategory: true,
      showSelectedRows: true, showMouseOverRowLine: true});
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));
  });

  await softStep('Data panel checkboxes', async () => {
    const errBefore = errorCount();
    await lcSetProps(page, {xColumnName: 'AGE', yColumnNames: ['AGE', 'HEIGHT']});

    // A demog line chart exposes no Show Null / Show Missing checkboxes (those
    // are grid props); the Data-section booleans that exist are Pack Categories
    // and Multi Axis. Toggle each off its default and read the stored value back.
    for (const p of ['packCategories', 'multiAxis']) {
      const cur = (await lcProps(page, p))[p];
      await lcSetProps(page, {[p]: !cur});
      expect((await lcProps(page, p))[p]).toBe(!cur);
    }
    expect(errorCount()).toBe(errBefore);

    // The toggled Data state must survive a layout save → close → restore.
    const combo = {packCategories: false, multiAxis: true};
    await lcSetProps(page, combo);

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

    const restored = await lcProps(page, 'packCategories', 'multiAxis');
    expect(restored).toEqual(combo);

    await page.evaluate(async (id) => {
      const saved = await grok.dapi.layouts.find(id);
      await grok.dapi.layouts.delete(saved);
    }, layoutId);
    await lcSetProps(page, {packCategories: true, multiAxis: false});
  });

  await softStep('GROK-17835 regression (SPGI)', async () => {
    await openDatasetWithLineChart(page, 'System:AppData/Chem/tests/spgi-100.csv');

    await lcSetProps(page, {multiAxis: true});

    await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Line chart')!;
      lc.props.splitColumnNames = ['Series', 'Scaffold Names'];
    });
    await page.waitForTimeout(500);

    await lcSetProps(page, {xColumnName: 'Chemist 521'});

    // Hovering the chart must not raise errors. grok.shell.warnings is
    // undefined on this build, so page and console errors are the channel.
    const errBefore = errorCount();

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

    expect(errorCount()).toBe(errBefore);
  });

  v.finishSpec();
});
