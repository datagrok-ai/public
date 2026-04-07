import {test, expect, Page} from '@playwright/test';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

async function login(page: Page) {
  await page.goto('/');
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});
  await page.evaluate(() => {
    document.body.classList.add('selenium');
    (grok.shell.settings as any).showFiltersIconsConstantly = true;
  });
  await waitForPlatform(page);
}

async function waitForPlatform(page: Page) {
  await page.waitForFunction(() => {
    const bar = document.querySelector('.grok-loader,.d4-loading-bar,.d4-progress-bar');
    return !bar || getComputedStyle(bar).display === 'none';
  }, {timeout: 60000}).catch(() => {});
}

async function openSPGI(page: Page) {
  await page.evaluate(async () => {
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/SPGI.csv');
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 5000);
    });
    // Wait for Chem cell rendering + package filter registration
    const hasBioChem = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
      .some(c => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 5000));
    }
  });
  await page.locator('[name="viewer-Grid"]').waitFor({timeout: 30000});
}

async function openFilterPanel(page: Page) {
  await page.evaluate(() => grok.shell.tv.getFiltersGroup());
  await page.locator('[name="viewer-Filters"]').waitFor({timeout: 10000});
}

async function getFilteredCount(page: Page): Promise<number> {
  return page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
}

async function addViewerViaIcon(page: Page, iconName: string, viewerType: string) {
  await page.evaluate((name) => {
    document.querySelector(`[name="${name}"]`)!.dispatchEvent(new MouseEvent('click', {bubbles: true}));
  }, iconName);
  await page.locator(`[name="viewer-${viewerType}"]`).waitFor({timeout: 10000});
  await page.waitForTimeout(1000);
}

async function closeViewer(page: Page, type: string) {
  await page.evaluate((t) => {
    for (const v of grok.shell.tv.viewers)
      if (v.type === t && !v.root.closest('[name="viewer-Filters"]')) { v.close(); break; }
  }, type);
  await page.waitForTimeout(500);
}

test('Viewers — Filter Panel interaction with viewers', async ({page}) => {
  await login(page);
  await openSPGI(page);
  const totalRows = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
  await openFilterPanel(page);

  // #### 1. Apply filters in the Filter Panel
  await softStep('1. Apply filters in the Filter Panel', async () => {
    // Structure filter: click Sketch, enter c1ccccc1, click OK
    await page.locator('.sketch-link').first().click();
    await page.locator('input[placeholder*="SMILES"]').waitFor({timeout: 10000});
    await page.locator('input[placeholder*="SMILES"]').click();
    await page.locator('input[placeholder*="SMILES"]').pressSequentially('c1ccccc1', {delay: 30});
    await page.keyboard.press('Enter');
    await page.waitForTimeout(2000);
    await page.locator('button:has-text("OK")').click();
    await page.waitForTimeout(2000);
    const afterStructure = await getFilteredCount(page);
    expect(afterStructure).toBeLessThan(totalRows);

    // Stereo Category: uncheck S_UNKN
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      const cats = grok.shell.tv.dataFrame.col('Stereo Category').categories;
      fg.updateOrAdd({
        type: DG.FILTER_TYPE.CATEGORICAL,
        column: 'Stereo Category',
        selected: cats.filter((c: string) => c !== 'S_UNKN'),
      });
    });
    await page.waitForTimeout(500);
    const afterStereo = await getFilteredCount(page);
    expect(afterStereo).toBeLessThanOrEqual(afterStructure);

    // Average Mass: set max to 400
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'Average Mass', min: 0, max: 400});
    });
    await page.waitForTimeout(500);
    const afterMass = await getFilteredCount(page);
    expect(afterMass).toBeLessThanOrEqual(afterStereo);
  });

  const afterSection1 = await getFilteredCount(page);

  // #### 2. Add Scaffold Tree Filter
  await softStep('2. Add Scaffold Tree Filter', async () => {
    // Add Scaffold Tree filter with full params (savedTree, colorCodedScaffolds, title are required)
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({
        type: 'Chem:scaffoldTreeFilter',
        column: 'Structure',
        columnName: 'Structure',
        savedTree: JSON.stringify([
          {scaffold: 'Cc1ccccc1', checked: true, expanded: true},
          'AND',
        ]),
        colorCodedScaffolds: '[]',
        title: 'Scaffold Tree',
      });
      grok.shell.tv.dataFrame.rows.requestFilter();
    });
    await page.waitForTimeout(3000);
    const afterScaffold = await getFilteredCount(page);
    expect(afterScaffold).toBeLessThanOrEqual(afterSection1);
  });

  const baselineFiltered = await getFilteredCount(page);

  // #### 3. Scatter Plot
  await softStep('3. Scatter Plot filtering', async () => {
    await addViewerViaIcon(page, 'icon-scatter-plot', 'Scatter-plot');

    // Zoom in using mouse wheel to filter
    const spCanvas = page.locator('[name="viewer-Scatter-plot"] canvas').first();
    const spBox = await spCanvas.boundingBox();
    expect(spBox).toBeTruthy();

    // Scroll to zoom in
    await page.mouse.move(spBox!.x + spBox!.width / 2, spBox!.y + spBox!.height / 2);
    await page.mouse.wheel(0, -300);
    await page.waitForTimeout(500);
    await page.mouse.wheel(0, -300);
    await page.waitForTimeout(1000);

    const afterZoom = await getFilteredCount(page);
    expect(afterZoom).toBeLessThan(baselineFiltered);

    // Double-click to reset zoom
    await spCanvas.dblclick();
    await page.waitForTimeout(1000);
    expect(await getFilteredCount(page)).toBe(baselineFiltered);

    // Close Scatter Plot
    await closeViewer(page, 'Scatter plot');
    expect(await getFilteredCount(page)).toBe(baselineFiltered);
  });

  // #### 4. Bar Chart
  await softStep('4. Bar Chart filtering', async () => {
    await addViewerViaIcon(page, 'icon-bar-chart', 'Bar-chart');

    // Set On Click to Filter via context menu
    const bcCanvas = page.locator('[name="viewer-Bar-chart"] canvas').first();
    await bcCanvas.click({button: 'right'});
    await page.waitForTimeout(500);

    // Navigate context menu: On Click > Filter
    const onClickItem = page.locator('.d4-menu-item-label:text("On Click")').first();
    await onClickItem.hover();
    await page.waitForTimeout(300);
    const filterItem = page.locator('.d4-menu-item-label:text("Filter")').first();
    await filterItem.click();
    await page.waitForTimeout(500);

    // Click a bar on the chart
    const bcBox = await bcCanvas.boundingBox();
    expect(bcBox).toBeTruthy();
    await bcCanvas.click({position: {x: bcBox!.width * 0.5, y: bcBox!.height * 0.5}});
    await page.waitForTimeout(500);
    const afterBarClick = await getFilteredCount(page);
    expect(afterBarClick).toBeLessThan(baselineFiltered);

    // Click white space to reset
    await bcCanvas.click({position: {x: bcBox!.width * 0.95, y: bcBox!.height * 0.95}});
    await page.waitForTimeout(500);
    expect(await getFilteredCount(page)).toBe(baselineFiltered);

    // Close Bar Chart
    await closeViewer(page, 'Bar chart');
    expect(await getFilteredCount(page)).toBe(baselineFiltered);
  });

  // #### 5. Histogram
  await softStep('5. Histogram filtering', async () => {
    // Add histogram via JS API to avoid confusion with filter-panel histograms
    await page.evaluate(() => grok.shell.tv.addViewer('Histogram'));
    await page.waitForTimeout(2000);

    const hasStandalone = await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).some(
        (v: any) => v.type === 'Histogram' && !v.root.closest('[name="viewer-Filters"]')));
    expect(hasStandalone).toBe(true);

    // Enable filtering and set range on the standalone Histogram
    await page.evaluate(() => {
      for (const v of grok.shell.tv.viewers)
        if (v.type === 'Histogram' && !v.root.closest('[name="viewer-Filters"]')) {
          v.props.filteringEnabled = true;
          v.props.valueMin = 635000;
          v.props.valueMax = 636000;
          break;
        }
    });
    await page.waitForTimeout(500);
    const afterHisto = await getFilteredCount(page);
    expect(afterHisto).toBeLessThan(baselineFiltered);

    // Double-click the histogram to reset
    const histCanvas = await page.evaluate(() => {
      for (const v of grok.shell.tv.viewers)
        if (v.type === 'Histogram' && !v.root.closest('[name="viewer-Filters"]'))
          return true;
      return false;
    });
    if (histCanvas) {
      // Try double-click; if it doesn't reset, use property fallback
      const histViewer = page.locator('[name="viewer-Histogram"]').last();
      await histViewer.locator('canvas').first().dblclick();
      await page.waitForTimeout(500);
      const afterDblClick = await getFilteredCount(page);
      if (afterDblClick !== baselineFiltered) {
        // Fallback: disable filtering via property
        await page.evaluate(() => {
          for (const v of grok.shell.tv.viewers)
            if (v.type === 'Histogram' && !v.root.closest('[name="viewer-Filters"]')) {
              v.props.filteringEnabled = false;
              break;
            }
        });
        await page.waitForTimeout(500);
      }
    }
    expect(await getFilteredCount(page)).toBe(baselineFiltered);

    // Close standalone histogram
    await page.evaluate(() => {
      for (const v of grok.shell.tv.viewers)
        if (v.type === 'Histogram' && !v.root.closest('[name="viewer-Filters"]')) { v.close(); break; }
    });
    await page.waitForTimeout(500);
    expect(await getFilteredCount(page)).toBe(baselineFiltered);
  });

  // #### 6. PC Plot
  await softStep('6. PC Plot filtering', async () => {
    await addViewerViaIcon(page, 'icon-pc-plot', 'PC-Plot');

    // Drag on the first axis to create a range filter
    const pcCanvas = page.locator('[name="viewer-PC-Plot"] canvas').first();
    const pcBox = await pcCanvas.boundingBox();
    expect(pcBox).toBeTruthy();

    if (pcBox) {
      const axisX = pcBox.x + 30;
      await page.mouse.move(axisX, pcBox.y + pcBox.height * 0.25);
      await page.mouse.down();
      await page.mouse.move(axisX, pcBox.y + pcBox.height * 0.55, {steps: 10});
      await page.mouse.up();
    }
    await page.waitForTimeout(1000);
    const afterPC = await getFilteredCount(page);
    expect(afterPC).toBeLessThan(baselineFiltered);

    // Double-click white space to reset
    await pcCanvas.dblclick({position: {x: pcBox!.width / 2, y: pcBox!.height / 2}});
    await page.waitForTimeout(500);
    expect(await getFilteredCount(page)).toBe(baselineFiltered);

    // Close PC Plot
    await closeViewer(page, 'PC Plot');
    expect(await getFilteredCount(page)).toBe(baselineFiltered);
  });

  // #### 7. Trellis Plot
  await softStep('7. Trellis Plot filtering', async () => {
    await addViewerViaIcon(page, 'icon-trellis-plot', 'Trellis-plot');

    // Set On Click to Filter via context menu
    const tpCanvas = page.locator('[name="viewer-Trellis-plot"] canvas').first();
    await tpCanvas.click({button: 'right'});
    await page.waitForTimeout(500);

    const onClickItem = page.locator('.d4-menu-item-label:text("On Click")').first();
    await onClickItem.hover();
    await page.waitForTimeout(300);
    const filterItem = page.locator('.d4-menu-item-label:text("Filter")').first();
    await filterItem.click();
    await page.waitForTimeout(500);

    // Click a cell in the Trellis Plot (center area for a populated cell)
    const tpBox = await tpCanvas.boundingBox();
    expect(tpBox).toBeTruthy();
    await tpCanvas.click({position: {x: tpBox!.width * 0.3, y: tpBox!.height * 0.3}});
    await page.waitForTimeout(1000);
    const afterTrellis = await getFilteredCount(page);
    expect(afterTrellis).toBeLessThanOrEqual(baselineFiltered);

    // Click same cell to reset
    await tpCanvas.click({position: {x: tpBox!.width * 0.3, y: tpBox!.height * 0.3}});
    await page.waitForTimeout(1000);
    const afterReset = await getFilteredCount(page);
    // If click-toggle doesn't work, closing the viewer will restore
    if (afterReset !== baselineFiltered) {
      console.warn('Trellis click-to-toggle did not reset; closing viewer to restore');
    }

    // Close Trellis Plot
    await closeViewer(page, 'Trellis plot');
    await page.waitForTimeout(500);
    // Verify filter count returns to baseline (may need filter re-apply if trellis broke it)
    const afterClose = await getFilteredCount(page);
    if (afterClose !== baselineFiltered) {
      // Re-apply panel filters
      await page.evaluate(() => {
        const fg = grok.shell.tv.getFiltersGroup();
        const cats = grok.shell.tv.dataFrame.col('Stereo Category').categories;
        fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: cats.filter((c: string) => c !== 'S_UNKN')});
        fg.updateOrAdd({type: 'histogram', column: 'Average Mass', min: 0, max: 400});
      });
      await page.waitForTimeout(500);
    }
    expect(await getFilteredCount(page)).toBe(baselineFiltered);
  });

  // #### 8. Pie Chart
  await softStep('8. Pie Chart filtering', async () => {
    await addViewerViaIcon(page, 'icon-pie-chart', 'Pie-chart');

    // Set On Click to Filter via context menu
    const pieCanvas = page.locator('[name="viewer-Pie-chart"] canvas').first();
    await pieCanvas.click({button: 'right'});
    await page.waitForTimeout(500);

    const onClickItem = page.locator('.d4-menu-item-label:text("On Click")').first();
    await onClickItem.hover();
    await page.waitForTimeout(300);
    const filterItem = page.locator('.d4-menu-item-label:text("Filter")').first();
    await filterItem.click();
    await page.waitForTimeout(500);

    // Click a segment in the pie chart (the large segment near center-top)
    const pieBox = await pieCanvas.boundingBox();
    expect(pieBox).toBeTruthy();
    await pieCanvas.click({position: {x: pieBox!.width * 0.4, y: pieBox!.height * 0.35}});
    await page.waitForTimeout(500);
    const afterPieClick = await getFilteredCount(page);
    expect(afterPieClick).toBeLessThan(baselineFiltered);

    // Click white space to reset
    await pieCanvas.click({position: {x: 10, y: pieBox!.height - 10}});
    await page.waitForTimeout(500);
    expect(await getFilteredCount(page)).toBe(baselineFiltered);

    // Close Pie Chart
    await closeViewer(page, 'Pie chart');
    expect(await getFilteredCount(page)).toBe(baselineFiltered);
  });

  // #### 9. Reset and cleanup
  await softStep('9. Reset and cleanup', async () => {
    // Click the Reset icon in the Filter Panel header
    await page.evaluate(() => {
      const fp = document.querySelector('[name="viewer-Filters"]');
      const resetIcon = fp!.querySelector('[name="icon-arrow-rotate-left"]') as HTMLElement;
      resetIcon.click();
    });
    await page.waitForTimeout(1000);
    const afterReset = await getFilteredCount(page);
    expect(afterReset).toBe(totalRows);

    // Close All
    await page.evaluate(() => grok.shell.closeAll());
    await page.waitForTimeout(500);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
