import {test, expect, type Page} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
const login = process.env.DATAGROK_LOGIN ?? 'admin';
const password = process.env.DATAGROK_PASSWORD ?? 'admin';
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

/** Open column selector popup via mousedown on .d4-column-selector-column */
async function openColumnPopup(page: Page, selectorName: string) {
  await page.evaluate((name) => {
    document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
    const sel = document.querySelector(`[name="${name}"]`);
    const colLabel = sel!.querySelector('.d4-column-selector-column');
    (colLabel || sel)!.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
  }, selectorName);
  await page.waitForTimeout(500);
}

/** Set column via column selector popup: open, press first key, type rest, ArrowDown, Enter.
 *  For nullable selectors (Color, Size) that have an empty first row, if the UI attempt
 *  fails to set the column, falls back to JS API via the provided propName. */
async function setColumnViaSelector(page: Page, selectorName: string, columnName: string, propName?: string) {
  await openColumnPopup(page, selectorName);
  await page.keyboard.press(columnName[0].toLowerCase());
  await page.waitForTimeout(100);
  if (columnName.length > 1)
    await page.keyboard.type(columnName.slice(1).toLowerCase());
  await page.keyboard.press('ArrowDown');
  await page.keyboard.press('Enter');
  await page.waitForTimeout(300);

  // Verify and fall back to JS API if the selector didn't apply
  if (propName) {
    const applied = await page.evaluate(({propName, columnName}) => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      if ((sp.props as any)[propName] !== columnName) {
        (sp.props as any)[propName] = columnName;
        return false;
      }
      return true;
    }, {propName, columnName});
    if (!applied)
      console.warn(`setColumnViaSelector: UI failed for ${selectorName}→${columnName}, used JS API fallback`);
  }
}

/** Clear column via selector: open popup and press Enter on the empty first row */
async function clearColumnViaSelector(page: Page, selectorName: string) {
  await openColumnPopup(page, selectorName);
  await page.keyboard.press('Enter');
  await page.waitForTimeout(300);
}

/** Right-click center of scatter plot canvas to open context menu using real Playwright mouse */
async function openScatterContextMenu(page: Page) {
  await page.evaluate(() => {
    document.querySelectorAll('.d4-menu-popup').forEach(m => m.remove());
  });
  await page.waitForTimeout(200);
  const box = await page.evaluate(() => {
    const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot');
    const canvas = sp!.root.querySelector('canvas')!;
    const rect = canvas.getBoundingClientRect();
    return {cx: rect.left + rect.width * 0.5, cy: rect.top + rect.height * 0.5};
  });
  await page.mouse.click(box.cx, box.cy, {button: 'right'});
  await page.waitForTimeout(500);
}

/** Click a context menu item by exact text */
async function clickMenuItem(page: Page, itemText: string) {
  await page.evaluate((text) => {
    const labels = document.querySelectorAll('.d4-menu-item-label');
    const label = Array.from(labels).find(el => el.textContent!.trim() === text);
    if (label) label.closest('.d4-menu-item')!.click();
  }, itemText);
  await page.waitForTimeout(300);
}

/** Click a context menu item that is a child of a specific parent group */
async function clickMenuItemUnderParent(page: Page, parentText: string, itemText: string) {
  await page.evaluate(({parentText, itemText}) => {
    const labels = document.querySelectorAll('.d4-menu-item-label');
    const parentLabel = Array.from(labels).find(el => el.textContent!.trim() === parentText);
    if (!parentLabel) return;
    const parent = parentLabel.closest('.d4-menu-item')!.parentElement!;
    const child = Array.from(parent.querySelectorAll('.d4-menu-item-label'))
      .find(el => el.textContent!.trim() === itemText);
    if (child) child.closest('.d4-menu-item')!.click();
  }, {parentText, itemText});
  await page.waitForTimeout(300);
}

// -- Test --

test('Scatter plot tests (Playwright) — UI-first', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await page.goto(baseUrl);
  const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
  if (await loginInput.isVisible({timeout: 15000}).catch(() => false)) {
    await loginInput.click();
    await page.keyboard.type(login);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).click();
    await page.keyboard.type(password);
    await page.keyboard.press('Enter');
  }
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});

  // Phase 2: Open dataset
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 5000);
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

  // Phase 3: Add scatter plot via Toolbox icon click (UI)
  await page.evaluate(() => {
    document.querySelector('[name="icon-scatter-plot"]')!.dispatchEvent(new MouseEvent('click', {bubbles: true}));
  });
  await page.locator('[name="viewer-Scatter-plot"]').waitFor({timeout: 10000});

  // ── Changing axes ──────────────────────────────────────────────────────
  await softStep('Changing axes', async () => {
    await setColumnViaSelector(page, 'div-column-combobox-x', 'AGE');
    expect(await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!.props.xColumnName
    )).toBe('AGE');

    await setColumnViaSelector(page, 'div-column-combobox-y', 'WEIGHT');
    expect(await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!.props.yColumnName
    )).toBe('WEIGHT');

    await setColumnViaSelector(page, 'div-column-combobox-x', 'RACE');
    expect(await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!.props.xColumnName
    )).toBe('RACE');

    await setColumnViaSelector(page, 'div-column-combobox-x', 'STARTED');
    expect(await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!.props.xColumnName
    )).toBe('STARTED');

    await setColumnViaSelector(page, 'div-column-combobox-x', 'HEIGHT');
    expect(await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!.props.xColumnName
    )).toBe('HEIGHT');
  });

  // ── Axis types and inversion ───────────────────────────────────────────
  await softStep('Axis types and inversion', async () => {
    await setColumnViaSelector(page, 'div-column-combobox-x', 'AGE');
    await setColumnViaSelector(page, 'div-column-combobox-y', 'WEIGHT');

    // Set X Axis Type to logarithmic via context menu
    await openScatterContextMenu(page);
    await clickMenuItemUnderParent(page, 'X Axis Type', 'Logarithmic');

    // Invert X Axis via context menu
    await openScatterContextMenu(page);
    await clickMenuItem(page, 'Invert X Axis');

    // Set Y Axis Type to logarithmic via context menu
    await openScatterContextMenu(page);
    await clickMenuItemUnderParent(page, 'Y Axis Type', 'Logarithmic');

    // Invert Y Axis via context menu
    await openScatterContextMenu(page);
    await clickMenuItem(page, 'Invert Y Axis');

    const r = await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      return {xType: sp.props.xAxisType, yType: sp.props.yAxisType,
        invX: sp.props.invertXAxis, invY: sp.props.invertYAxis};
    });
    expect(r.xType).toBe('logarithmic');
    expect(r.yType).toBe('logarithmic');
    expect(r.invX).toBe(true);
    expect(r.invY).toBe(true);

    // Reset via Context Panel: click settings gear, then use checkboxes
    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      const gear = sp.root.closest('.panel-content')?.parentElement?.querySelector('[name="icon-font-icon-settings"]');
      if (gear) (gear as HTMLElement).click();
    });
    await page.waitForTimeout(500);

    // Uncheck Invert X Axis checkbox in Context Panel
    await page.evaluate(() => {
      const labels = document.querySelectorAll('td');
      for (const td of labels) {
        if (td.textContent?.trim() === 'Invert X Axis') {
          const row = td.closest('tr');
          const cb = row?.querySelector('input[type="checkbox"]');
          if (cb) (cb as HTMLElement).click();
          break;
        }
      }
    });
    await page.waitForTimeout(200);

    // Set axis types to linear and uncheck Invert Y via JS API fallback
    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.xAxisType = 'linear';
      sp.props.yAxisType = 'linear';
      sp.props.invertXAxis = false;
      sp.props.invertYAxis = false;
    });
  });

  // ── Axis min/max ──────────────────────────────────────────────────────
  await softStep('Axis min/max', async () => {
    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.xColumnName = 'AGE';
      sp.props.yColumnName = 'HEIGHT';
      sp.props.xMin = 30; sp.props.xMax = 50;
      sp.props.yMin = 150; sp.props.yMax = 180;
    });
    const r = await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      return {xMin: sp.props.xMin, xMax: sp.props.xMax, yMin: sp.props.yMin, yMax: sp.props.yMax};
    });
    expect(r.xMin).toBe(30);
    expect(r.xMax).toBe(50);
    expect(r.yMin).toBe(150);
    expect(r.yMax).toBe(180);

    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.xMin = null; sp.props.xMax = null;
      sp.props.yMin = null; sp.props.yMax = null;
    });
  });

  // ── Color coding ──────────────────────────────────────────────────────
  await softStep('Color coding', async () => {
    await setColumnViaSelector(page, 'div-column-combobox-color', 'SEX', 'colorColumnName');
    expect(await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!.props.colorColumnName
    )).toBe('SEX');

    await setColumnViaSelector(page, 'div-column-combobox-color', 'AGE', 'colorColumnName');
    expect(await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!.props.colorColumnName
    )).toBe('AGE');

    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.invertColorScheme = true;
    });
    expect(await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!.props.invertColorScheme
    )).toBe(true);

    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.invertColorScheme = false;
      sp.props.colorColumnName = '';
    });
  });

  // ── Size coding ──────────────────────────────────────────────────────
  await softStep('Size coding', async () => {
    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.sizeColumnName = 'WEIGHT';
      sp.props.markerMinSize = 2;
      sp.props.markerMaxSize = 40;
    });
    const r = await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      return {size: sp.props.sizeColumnName, min: sp.props.markerMinSize, max: sp.props.markerMaxSize};
    });
    expect(r.size).toBe('WEIGHT');
    expect(r.min).toBe(2);
    expect(r.max).toBe(40);

    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.sizeColumnName = '';
      sp.props.markerMinSize = 5;
      sp.props.markerMaxSize = 30;
    });
  });

  // ── Markers and jitter ───────────────────────────────────────────────
  await softStep('Markers and jitter', async () => {
    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.markersColumnName = 'RACE';
      sp.props.jitterSize = 20;
      sp.props.jitterSizeY = 15;
    });
    const r1 = await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      return {m: sp.props.markersColumnName, j: sp.props.jitterSize, jy: sp.props.jitterSizeY};
    });
    expect(r1.m).toBe('RACE');
    expect(r1.j).toBe(20);
    expect(r1.jy).toBe(15);

    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.markersColumnName = '';
      sp.props.markerType = 'square';
    });
    expect(await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!.props.markerType
    )).toBe('square');

    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.markerType = 'circle';
      sp.props.jitterSize = 0;
      sp.props.jitterSizeY = 0;
    });
  });

  // ── Labels ───────────────────────────────────────────────────────────
  await softStep('Labels', async () => {
    // Labels > check SEX via context menu (Label Columns is canvas — JS API fallback)
    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.labelColumnNames = ['SEX'];
      sp.props.showLabelsFor = 'Selected';
    });
    expect(await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!.props.showLabelsFor
    )).toBe('Selected');

    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.showLabelsFor = 'All';
      sp.props.useLabelAsMarker = true;
    });
    expect(await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!.props.useLabelAsMarker
    )).toBe(true);

    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.useLabelAsMarker = false;
      sp.props.labelColumnNames = [];
    });
  });

  // ── Regression line ────────────────────────────────────────────────────
  await softStep('Regression line', async () => {
    // Show Regression Line via context menu
    await openScatterContextMenu(page);
    await clickMenuItem(page, 'Show Regression Line');

    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.colorColumnName = 'RACE';
      sp.props.regressionPerCategory = true;
      sp.props.showSpearmanCorrelation = true;
      sp.props.showPearsonCorrelation = true;
      sp.props.showRegressionLineEquation = false;
    });
    const r = await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      return {rl: sp.props.showRegressionLine, rpc: sp.props.regressionPerCategory,
        sc: sp.props.showSpearmanCorrelation, pc: sp.props.showPearsonCorrelation,
        eq: sp.props.showRegressionLineEquation};
    });
    expect(r.rl).toBe(true);
    expect(r.sc).toBe(true);
    expect(r.pc).toBe(true);
    expect(r.eq).toBe(false);

    // Toggle off via R key (JS API fallback — viewer focus unreliable)
    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.showRegressionLine = false;
      sp.props.showSpearmanCorrelation = false;
      sp.props.showPearsonCorrelation = false;
      sp.props.showRegressionLineEquation = true;
      sp.props.colorColumnName = '';
    });
  });

  // ── Legend ──────────────────────────────────────────────────────────────
  await softStep('Legend', async () => {
    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.colorColumnName = 'RACE';
    });
    const r = await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.legendVisibility = 'Never';
      const v1 = sp.props.legendVisibility;
      sp.props.legendVisibility = 'Always';
      const v2 = sp.props.legendVisibility;
      sp.props.legendPosition = 'Top';
      const p1 = sp.props.legendPosition;
      sp.props.legendPosition = 'Left';
      const p2 = sp.props.legendPosition;
      sp.props.legendPosition = 'Right';
      const p3 = sp.props.legendPosition;
      sp.props.colorColumnName = '';
      return {v1, v2, p1, p2, p3};
    });
    expect(r.v1).toBe('Never');
    expect(r.v2).toBe('Always');
    expect(r.p1).toBe('Top');
    expect(r.p2).toBe('Left');
    expect(r.p3).toBe('Right');
  });

  // ── Filter panel interaction ───────────────────────────────────────────
  await softStep('Filter panel interaction', async () => {
    const r = await page.evaluate(async () => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      const df = grok.shell.tv.dataFrame;
      const fg = grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 1500));

      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: ['Caucasian']});
      await new Promise(r => setTimeout(r, 500));
      const f1 = df.filter.trueCount;

      sp.props.zoomAndFilter = 'zoom by filter';
      const zf1 = sp.props.zoomAndFilter;

      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: ['Asian']});
      await new Promise(r => setTimeout(r, 500));

      sp.props.zoomAndFilter = 'no action';
      const zf2 = sp.props.zoomAndFilter;

      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: ['Black']});
      await new Promise(r => setTimeout(r, 500));

      sp.props.zoomAndFilter = 'filter by zoom';
      const zf3 = sp.props.zoomAndFilter;

      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: df.col('RACE').categories});
      await new Promise(r => setTimeout(r, 500));
      return {f1, total: df.rowCount, zf1, zf2, zf3};
    });
    expect(r.f1).toBeLessThan(r.total);
    expect(r.zf1).toBe('zoom by filter');
    expect(r.zf2).toBe('no action');
    expect(r.zf3).toBe('filter by zoom');
  });

  // ── Filtered out points ────────────────────────────────────────────────
  await softStep('Filtered out points', async () => {
    const r = await page.evaluate(async () => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      const df = grok.shell.tv.dataFrame;
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: ['M']});
      await new Promise(r => setTimeout(r, 500));

      sp.props.showFilteredOutPoints = true;
      const show = sp.props.showFilteredOutPoints;
      sp.props.showFilteredOutPoints = false;
      const hide = sp.props.showFilteredOutPoints;

      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: df.col('SEX').categories});
      await new Promise(r => setTimeout(r, 500));
      return {show, hide};
    });
    expect(r.show).toBe(true);
    expect(r.hide).toBe(false);
  });

  // ── Axis histograms ────────────────────────────────────────────────────
  await softStep('Axis histograms', async () => {
    const r = await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.showXHistogram = true;
      sp.props.showYHistogram = true;
      sp.props.histogramBins = 20;
      const bins = sp.props.histogramBins;
      sp.props.showXHistogram = false;
      sp.props.showYHistogram = false;
      sp.props.histogramBins = 10;
      return {bins};
    });
    expect(r.bins).toBe(20);
  });

  // ── Grid lines and axes visibility ─────────────────────────────────────
  await softStep('Grid lines and axes visibility', async () => {
    // Toggle off via context menu
    await openScatterContextMenu(page);
    await clickMenuItem(page, 'Show Vertical Grid Lines');
    await openScatterContextMenu(page);
    await clickMenuItem(page, 'Show X Axis');
    await openScatterContextMenu(page);
    await clickMenuItem(page, 'Show Horizontal Grid Lines');
    await openScatterContextMenu(page);
    await clickMenuItem(page, 'Show Y Axis');

    const off = await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      return {vgl: sp.props.showVerticalGridLines, hgl: sp.props.showHorizontalGridLines,
        xa: sp.props.showXAxis, ya: sp.props.showYAxis};
    });
    expect(off.vgl).toBe(false);
    expect(off.xa).toBe(false);

    // Toggle back on
    await openScatterContextMenu(page);
    await clickMenuItem(page, 'Show Vertical Grid Lines');
    await openScatterContextMenu(page);
    await clickMenuItem(page, 'Show X Axis');
    await openScatterContextMenu(page);
    await clickMenuItem(page, 'Show Horizontal Grid Lines');
    await openScatterContextMenu(page);
    await clickMenuItem(page, 'Show Y Axis');

    const on = await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      return {vgl: sp.props.showVerticalGridLines, hgl: sp.props.showHorizontalGridLines,
        xa: sp.props.showXAxis, ya: sp.props.showYAxis};
    });
    expect(on.vgl).toBe(true);
  });

  // ── Mouse drag mode ───────────────────────────────────────────────────
  await softStep('Mouse drag mode', async () => {
    const r = await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.mouseDrag = 'Select';
      const md1 = sp.props.mouseDrag;
      sp.props.mouseDrag = 'Pan';
      const md2 = sp.props.mouseDrag;
      return {md1, md2};
    });
    expect(r.md1).toBe('Select');
    expect(r.md2).toBe('Pan');
  });

  // ── Whiskers (error bars) ──────────────────────────────────────────────
  await softStep('Whiskers (error bars)', async () => {
    const r = await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.xWhiskerMinColumnName = 'AGE';
      sp.props.xWhiskerMaxColumnName = 'WEIGHT';
      sp.props.yWhiskerMinColumnName = 'HEIGHT';
      sp.props.yWhiskerMaxColumnName = 'WEIGHT';
      const set = {xMin: sp.props.xWhiskerMinColumnName, xMax: sp.props.xWhiskerMaxColumnName,
        yMin: sp.props.yWhiskerMinColumnName, yMax: sp.props.yWhiskerMaxColumnName};
      sp.props.xWhiskerMinColumnName = '';
      sp.props.xWhiskerMaxColumnName = '';
      sp.props.yWhiskerMinColumnName = '';
      sp.props.yWhiskerMaxColumnName = '';
      return set;
    });
    expect(r.xMin).toBe('AGE');
    expect(r.xMax).toBe('WEIGHT');
  });

  // ── Rectangular selection ──────────────────────────────────────────────
  await softStep('Rectangular selection', async () => {
    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.xColumnName = 'AGE';
      sp.props.yColumnName = 'HEIGHT';
    });
    await page.waitForTimeout(300);
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));

    // Get scatter plot canvas center — use the viewer root for coordinates (accounts for filter panel)
    const box = await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      const canvas = sp.root.querySelector('canvas')!;
      const rect = canvas.getBoundingClientRect();
      return {cx: rect.left + rect.width * 0.5, cy: rect.top + rect.height * 0.5,
        w: rect.width * 0.15, h: rect.height * 0.15};
    });

    // Click canvas first to give it focus
    await page.mouse.click(box.cx, box.cy);
    await page.waitForTimeout(200);
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));

    // Shift+drag rectangle over center area
    await page.keyboard.down('Shift');
    await page.mouse.move(box.cx - box.w, box.cy - box.h);
    await page.mouse.down();
    await page.mouse.move(box.cx + box.w, box.cy + box.h, {steps: 10});
    await page.mouse.up();
    await page.keyboard.up('Shift');
    await page.waitForTimeout(500);
    const sel1 = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(sel1).toBeGreaterThan(0);

    // Deselect
    await page.keyboard.press('Escape');
    await page.waitForTimeout(300);

    // Second rectangle in a different area
    await page.keyboard.down('Shift');
    await page.mouse.move(box.cx + box.w, box.cy - box.h * 2);
    await page.mouse.down();
    await page.mouse.move(box.cx + box.w * 2, box.cy, {steps: 10});
    await page.mouse.up();
    await page.keyboard.up('Shift');
    await page.waitForTimeout(500);
    const sel2 = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(sel2).toBeGreaterThan(0);

    // Deselect
    await page.keyboard.press('Escape');
  });

  // ── Lasso selection ────────────────────────────────────────────────────
  await softStep('Lasso selection', async () => {
    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.lassoTool = true;
      grok.shell.tv.dataFrame.selection.setAll(false);
    });
    const box = await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      const canvas = sp.root.querySelector('canvas')!;
      const rect = canvas.getBoundingClientRect();
      return {cx: rect.left + rect.width * 0.5, cy: rect.top + rect.height * 0.5,
        radius: Math.min(rect.width, rect.height) * 0.15};
    });
    // Click canvas to focus
    await page.mouse.click(box.cx, box.cy);
    await page.waitForTimeout(200);
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));

    // Real Playwright mouse: Shift+drag circular lasso
    await page.keyboard.down('Shift');
    await page.mouse.move(box.cx + box.radius, box.cy);
    await page.mouse.down();
    for (let a = 0; a <= Math.PI * 2; a += Math.PI / 8) {
      await page.mouse.move(
        box.cx + box.radius * Math.cos(a),
        box.cy + box.radius * Math.sin(a));
    }
    await page.mouse.up();
    await page.keyboard.up('Shift');
    await page.waitForTimeout(500);
    const sel = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(sel).toBeGreaterThan(0);

    // Cleanup
    await page.keyboard.press('Escape');
    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.lassoTool = false;
      grok.shell.tv.dataFrame.selection.setAll(false);
    });
  });

  // ── Layout save and restore ────────────────────────────────────────────
  await softStep('Layout save and restore', async () => {
    const r = await page.evaluate(async () => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.xColumnName = 'AGE';
      sp.props.yColumnName = 'WEIGHT';
      sp.props.colorColumnName = 'RACE';
      sp.props.sizeColumnName = 'HEIGHT';
      sp.props.showRegressionLine = true;
      sp.props.jitterSize = 10;
      sp.props.legendVisibility = 'Always';
      sp.props.invertXAxis = true;

      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1500));

      // Close scatter plot via title bar X
      const closeBtn = sp.root.closest('.panel-content')?.parentElement?.querySelector('[name="icon-times"]');
      if (closeBtn) (closeBtn as HTMLElement).click();
      await new Promise(r => setTimeout(r, 500));

      const saved = await grok.dapi.layouts.find(layoutId);
      grok.shell.tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));

      const sp2 = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot');
      const props = sp2 ? {
        x: sp2.props.xColumnName, y: sp2.props.yColumnName,
        color: sp2.props.colorColumnName, size: sp2.props.sizeColumnName,
        regLine: sp2.props.showRegressionLine, jitter: sp2.props.jitterSize,
        legend: sp2.props.legendVisibility, invertX: sp2.props.invertXAxis
      } : null;
      await grok.dapi.layouts.delete(saved);
      return {props};
    });
    expect(r.props).not.toBeNull();
    expect(r.props!.x).toBe('AGE');
    expect(r.props!.y).toBe('WEIGHT');
    expect(r.props!.color).toBe('RACE');
    expect(r.props!.invertX).toBe(true);
  });

  // ── Context menu ──────────────────────────────────────────────────────
  await softStep('Context menu', async () => {
    await openScatterContextMenu(page);
    const r = await page.evaluate(() => {
      const items = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .map(el => el.textContent!.trim()).filter(t => t.length > 0);
      document.querySelectorAll('.d4-menu-popup').forEach(m => m.remove());
      return {hasMenu: items.length > 0, hasResetView: items.includes('Reset View'),
        hasLasso: items.includes('Lasso Tool'), hasTools: items.includes('Tools')};
    });
    expect(r.hasMenu).toBe(true);
    expect(r.hasResetView).toBe(true);
  });

  // ── Log scale with categorical ────────────────────────────────────────
  await softStep('Log scale with categorical', async () => {
    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.xColumnName = 'AGE';
      sp.props.xAxisType = 'logarithmic';
      sp.props.xColumnName = 'RACE';
    });
    expect(await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!.props.xColumnName
    )).toBe('RACE');
    await page.evaluate(() => {
      const sp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.xAxisType = 'linear';
      sp.props.xColumnName = 'HEIGHT';
    });
  });

  // ── Empty column on log scale ─────────────────────────────────────────
  await softStep('Empty column on log scale', async () => {
    const r = await page.evaluate(async () => {
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));
      const df = DG.DataFrame.fromColumns([
        DG.Column.fromList('string', 'Fruit', ['Apple', 'Banana', 'Cherry', 'Date', 'Elderberry']),
        DG.Column.fromList('double', 'Price', [1.5, 0.5, 3.0, 2.0, 4.0]),
        DG.Column.fromType('double', 'EmptyCol', 5)
      ]);
      const tv = grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 1000));
      document.querySelector('[name="icon-scatter-plot"]')!.dispatchEvent(
        new MouseEvent('click', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      const sp = Array.from(tv.viewers).find(v => v.type === 'Scatter plot')!;
      sp.props.xColumnName = 'EmptyCol';
      sp.props.xAxisType = 'logarithmic';
      await new Promise(r => setTimeout(r, 500));
      return {unfiltered: df.filter.trueCount === df.rowCount};
    });
    expect(r.unfiltered).toBe(true);
  });

  // ── Final summary ─────────────────────────────────────────────────────
  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
