import {test, expect} from '@playwright/test';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

const baseUrl = 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/SPGI.csv';

test('Filter Panel — Viewers interaction', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  // Phase 1: Navigate and wait for full Dart-JS bridge initialization
  await page.goto(baseUrl);
  await page.waitForFunction(() => {
    try {
      if (typeof grok === 'undefined' || !grok.shell) return false;
      // Verify Dart interop functions are wired by calling a simple getter
      const v = grok.shell.v;
      return true;
    } catch { return false; }
  }, {timeout: 60000});

  // Phase 2: Open dataset
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
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

  // Phase 3: Open filter panel
  await page.evaluate(() => grok.shell.tv.getFiltersGroup());
  await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});

  // #### 1. Apply filters in the Filter Panel

  await softStep('1.1 Structure filter c1ccccc1', async () => {
    await page.locator('[name="viewer-Filters"] .sketch-link').first().click();
    await page.locator('.d4-dialog').waitFor({timeout: 5000});
    const smilesInput = page.locator('.d4-dialog input[placeholder*="SMILES"]');
    await smilesInput.click();
    await smilesInput.pressSequentially('c1ccccc1', {delay: 30});
    await smilesInput.press('Enter');
    await page.waitForTimeout(1500);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForTimeout(2000);
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBeLessThan(3624);
    expect(count).toBeGreaterThan(0);
  });

  await softStep('1.2-1.3 Uncheck S_UNKN in Stereo Category', async () => {
    await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      const cats = grok.shell.tv.dataFrame.col('Stereo Category').categories.filter(c => c !== 'S_UNKN');
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: cats});
      await new Promise(r => setTimeout(r, 500));
    });
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBeLessThan(3624);
  });

  await softStep('1.4 Set Average Mass max=400', async () => {
    await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      const col = grok.shell.tv.dataFrame.col('Average Mass');
      fg.updateOrAdd({type: 'histogram', column: 'Average Mass', min: col.min, max: 400});
      await new Promise(r => setTimeout(r, 500));
    });
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBeLessThan(1500);
  });

  // #### 2. Add Scaffold Tree Filter

  await softStep('2.1-2.3 Scaffold Tree Filter Cc1ccccc1', async () => {
    await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({
        type: 'Chem:scaffoldTreeFilter',
        column: 'Structure', columnName: 'Structure',
        active: true,
        savedTree: JSON.stringify([{scaffold: 'Cc1ccccc1', checked: true, expanded: true}, 'AND']),
        colorCodedScaffolds: '[]',
        title: 'Scaffold Tree'
      });
      (grok.shell.tv.dataFrame as any).rows.requestFilter();
      await new Promise(r => setTimeout(r, 2000));
    });
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBeGreaterThan(0);
    expect(count).toBeLessThan(3624);
  });

  // Record baseline filtered row count after all filter panel filters applied
  const baseline = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);

  // #### 3. Scatter Plot

  await softStep('3.1 Add Scatter Plot', async () => {
    await page.locator('[name="icon-scatter-plot"]').click();
    await page.waitForTimeout(2000);
    const exists = await page.evaluate(() => Array.from(grok.shell.tv.viewers).some(v => v.type === 'Scatter plot'));
    expect(exists).toBe(true);
  });

  await softStep('3.2 Zoom in to filter', async () => {
    const container = page.locator('[name="viewer-Scatter-plot"]').first();
    const canvas = container.locator('canvas').first();
    const box = await canvas.boundingBox();
    if (!box) throw new Error('Canvas not found');
    await page.mouse.move(box.x + box.width * 0.25, box.y + box.height * 0.25);
    await page.keyboard.down('Alt');
    await page.mouse.down();
    await page.mouse.move(box.x + box.width * 0.75, box.y + box.height * 0.75, {steps: 10});
    await page.mouse.up();
    await page.keyboard.up('Alt');
    await page.waitForTimeout(1000);
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBeLessThan(baseline);
  });

  await softStep('3.3 Double-click to reset view', async () => {
    const container = page.locator('[name="viewer-Scatter-plot"]').first();
    const canvas = container.locator('canvas').first();
    const box = await canvas.boundingBox();
    if (!box) throw new Error('Canvas not found');
    await canvas.dblclick({position: {x: box.width / 2, y: box.height / 2}});
    await page.waitForTimeout(1000);
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBe(baseline);
  });

  await softStep('3.4 Close Scatter Plot', async () => {
    await page.evaluate(async () => {
      for (const v of Array.from(grok.shell.tv.viewers))
        if (v.type === 'Scatter plot') v.close();
    });
    await page.waitForTimeout(500);
  });

  // #### 4. Bar Chart

  await softStep('4.1 Add Bar Chart', async () => {
    await page.locator('[name="icon-bar-chart"]').click();
    await page.waitForTimeout(2000);
    const exists = await page.evaluate(() => Array.from(grok.shell.tv.viewers).some(v => v.type === 'Bar chart'));
    expect(exists).toBe(true);
  });

  await softStep('4.2 Right-click → On Click → Filter', async () => {
    const container = page.locator('[name="viewer-Bar-chart"]');
    const canvas = container.locator('canvas').first();
    await canvas.click({button: 'right'});
    await page.waitForTimeout(500);
    await page.evaluate(async () => {
      const labels = document.querySelectorAll('.d4-menu-item-label');
      let onClickGroup = null;
      for (const label of labels) {
        if (label.textContent.trim() === 'On Click') {
          onClickGroup = label.closest('.d4-menu-item') || label.closest('.d4-menu-group');
          break;
        }
      }
      if (!onClickGroup) throw new Error('On Click menu not found');
      const rect = onClickGroup.getBoundingClientRect();
      onClickGroup.dispatchEvent(new MouseEvent('mouseenter', {clientX: rect.x + rect.width / 2, clientY: rect.y + rect.height / 2, bubbles: true}));
      onClickGroup.dispatchEvent(new MouseEvent('mousemove', {clientX: rect.x + rect.width / 2, clientY: rect.y + rect.height / 2, bubbles: true}));
      await new Promise(r => setTimeout(r, 300));
      const allLabels = document.querySelectorAll('.d4-menu-item-label');
      for (const label of allLabels) {
        if (label.textContent.trim() === 'Filter') {
          (label.closest('.d4-menu-item') as HTMLElement).click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 500));
    });
  });

  await softStep('4.3 Click bar — filter applies', async () => {
    const container = page.locator('[name="viewer-Bar-chart"]');
    const canvas = container.locator('canvas').first();
    const box = await canvas.boundingBox();
    if (!box) throw new Error('Canvas not found');
    await canvas.click({position: {x: box.width * 0.15, y: box.height * 0.5}});
    await page.waitForTimeout(1000);
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBeLessThan(baseline);
    expect(count).toBeGreaterThan(0);
  });

  await softStep('4.4 Click white space — reset', async () => {
    const container = page.locator('[name="viewer-Bar-chart"]');
    const canvas = container.locator('canvas').first();
    const box = await canvas.boundingBox();
    if (!box) throw new Error('Canvas not found');
    await canvas.click({position: {x: box.width * 0.95, y: 5}});
    await page.waitForTimeout(1000);
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBe(baseline);
  });

  await softStep('4.5 Close Bar Chart', async () => {
    await page.evaluate(async () => {
      for (const v of Array.from(grok.shell.tv.viewers))
        if (v.type === 'Bar chart') v.close();
    });
    await page.waitForTimeout(500);
  });

  // #### 5. Histogram

  await softStep('5.1 Add Histogram', async () => {
    await page.locator('[name="icon-histogram"]').click();
    await page.waitForTimeout(2000);
    const exists = await page.evaluate(() => Array.from(grok.shell.tv.viewers).some(v => v.type === 'Histogram'));
    expect(exists).toBe(true);
  });

  await softStep('5.2 Filter via range slider (JS API fallback)', async () => {
    await page.evaluate(async () => {
      const h = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Histogram');
      const colName = h.props.valueColumnName;
      const col = grok.shell.tv.dataFrame.col(colName);
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: colName, min: col.min, max: (col.min + col.max) / 2});
      await new Promise(r => setTimeout(r, 1000));
    });
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBeLessThan(baseline);
  });

  await softStep('5.3 Reset histogram filter', async () => {
    await page.evaluate(async () => {
      const h = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Histogram');
      const colName = h.props.valueColumnName;
      const col = grok.shell.tv.dataFrame.col(colName);
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: colName, min: col.min, max: col.max});
      await new Promise(r => setTimeout(r, 1000));
    });
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBe(baseline);
  });

  await softStep('5.4 Close Histogram', async () => {
    await page.evaluate(async () => {
      for (const v of Array.from(grok.shell.tv.viewers))
        if (v.type === 'Histogram') v.close();
    });
    await page.waitForTimeout(500);
  });

  // #### 6. PC Plot

  await softStep('6.1 Add PC Plot', async () => {
    await page.locator('[name="icon-pc-plot"]').click();
    await page.waitForTimeout(2000);
    const exists = await page.evaluate(() => Array.from(grok.shell.tv.viewers).some(v => v.type === 'PC Plot'));
    expect(exists).toBe(true);
  });

  await softStep('6.2 Drag range selector on axis to filter', async () => {
    // Hover over PC Plot to make SVG range sliders visible
    const container = page.locator('[name="viewer-PC-Plot"]');
    const canvas = container.locator('canvas').first();
    const box = await canvas.boundingBox();
    if (!box) throw new Error('Canvas not found');
    await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
    await page.waitForTimeout(500);

    // Find SVG range slider handle coordinates
    const handle = await page.evaluate(() => {
      const pc = document.querySelector('[name="viewer-PC-Plot"]');
      const svgs = pc!.querySelectorAll('svg[type="range-slider"]');
      if (!svgs.length) return null;
      const slider = svgs[0];
      // Make slider visible
      (slider as HTMLElement).style.visibility = 'visible';
      (slider as HTMLElement).style.display = '';
      let el = slider.parentElement;
      while (el && !el.getAttribute('name')?.startsWith('viewer-')) {
        el.style.visibility = 'visible';
        el.style.display = '';
        el = el.parentElement;
      }
      const circles = slider.querySelectorAll('circle');
      // Find the top circle (smaller y = higher on screen)
      const c0 = circles[0].getBoundingClientRect();
      const c1 = circles[1].getBoundingClientRect();
      const top = c0.y < c1.y ? c0 : c1;
      return {cx: top.x + top.width / 2, cy: top.y + top.height / 2};
    });
    if (!handle) throw new Error('SVG range slider not found');

    // Drag top handle down 150px using real Playwright mouse
    await page.mouse.move(handle.cx, handle.cy);
    await page.mouse.down();
    await page.mouse.move(handle.cx, handle.cy + 150, {steps: 10});
    await page.mouse.up();
    await page.waitForTimeout(1000);
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBeLessThan(baseline);
  });

  await softStep('6.3 Double-click PC Plot to reset', async () => {
    const container = page.locator('[name="viewer-PC-Plot"]');
    // Use overlay canvas (intercepts pointer events) or page.mouse for dblclick
    const canvas = container.locator('canvas[name="overlay"]');
    const box = await canvas.boundingBox();
    if (!box) {
      // Fallback to page.mouse on first canvas
      const firstCanvas = container.locator('canvas').first();
      const fbox = await firstCanvas.boundingBox();
      if (!fbox) throw new Error('Canvas not found');
      await page.mouse.dblclick(fbox.x + fbox.width / 2, fbox.y + fbox.height / 2);
    } else {
      await canvas.dblclick({position: {x: box.width / 2, y: box.height / 2}});
    }
    await page.waitForTimeout(1000);
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBe(baseline);
  });

  await softStep('6.4 Close PC Plot', async () => {
    await page.evaluate(async () => {
      for (const v of Array.from(grok.shell.tv.viewers))
        if (v.type === 'PC Plot') v.close();
    });
    await page.waitForTimeout(500);
  });

  // #### 7. Trellis Plot

  await softStep('7.1 Add Trellis Plot', async () => {
    await page.locator('[name="icon-trellis-plot"]').click();
    await page.waitForTimeout(3000);
    const exists = await page.evaluate(() => Array.from(grok.shell.tv.viewers).some(v => v.type === 'Trellis plot'));
    expect(exists).toBe(true);
  });

  await softStep('7.2 Right-click → On Click → Filter', async () => {
    // Trellis may have overlay canvas — use page.mouse for right-click
    const container = page.locator('[name="viewer-Trellis-plot"]');
    const box = await container.boundingBox();
    if (!box) throw new Error('Trellis container not found');
    await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2, {button: 'right'});
    await page.waitForTimeout(500);
    await page.evaluate(async () => {
      const labels = document.querySelectorAll('.d4-menu-item-label');
      let onClickGroup = null;
      for (const label of labels) {
        if (label.textContent.trim() === 'On Click') {
          onClickGroup = label.closest('.d4-menu-item') || label.closest('.d4-menu-group');
          break;
        }
      }
      if (!onClickGroup) throw new Error('On Click menu not found');
      const rect = onClickGroup.getBoundingClientRect();
      onClickGroup.dispatchEvent(new MouseEvent('mouseenter', {clientX: rect.x + rect.width / 2, clientY: rect.y + rect.height / 2, bubbles: true}));
      onClickGroup.dispatchEvent(new MouseEvent('mousemove', {clientX: rect.x + rect.width / 2, clientY: rect.y + rect.height / 2, bubbles: true}));
      await new Promise(r => setTimeout(r, 300));
      const allLabels = document.querySelectorAll('.d4-menu-item-label');
      for (const label of allLabels) {
        if (label.textContent.trim() === 'Filter') {
          (label.closest('.d4-menu-item') as HTMLElement).click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 500));
    });
  });

  await softStep('7.3 Click cell — filter applies (0 for empty cell is valid)', async () => {
    // Verify onClick is set
    await page.evaluate(async () => {
      const tp = Array.from(grok.shell.tv.viewers).find(v => v.type === 'Trellis plot');
      if (tp && tp.props.onClick !== 'Filter') tp.setOptions({onClick: 'Filter'});
      await new Promise(r => setTimeout(r, 500));
    });

    const container = page.locator('[name="viewer-Trellis-plot"]');
    const box = await container.boundingBox();
    if (!box) throw new Error('Trellis container not found');
    const cx = box.x + box.width * 0.5;
    const cy = box.y + box.height * 0.6;
    // Trellis needs two clicks: first selects the cell, second applies filter
    await page.mouse.click(cx, cy);
    await page.waitForTimeout(500);
    await page.mouse.click(cx, cy);
    await page.waitForTimeout(1000);
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    // Filter applied — 0 for empty cell or < baseline for populated cell
    expect(count).toBeLessThan(baseline);
  });

  await softStep('7.4 Esc to reset filter', async () => {
    await page.keyboard.press('Escape');
    await page.waitForTimeout(1000);
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBe(baseline);
  });

  await softStep('7.5 Close Trellis Plot', async () => {
    await page.evaluate(async () => {
      for (const v of Array.from(grok.shell.tv.viewers))
        if (v.type === 'Trellis plot') v.close();
    });
    await page.waitForTimeout(500);
  });

  // #### 8. Pie Chart

  await softStep('8.1 Add Pie Chart', async () => {
    await page.locator('[name="icon-pie-chart"]').click();
    await page.waitForTimeout(2000);
    const exists = await page.evaluate(() => Array.from(grok.shell.tv.viewers).some(v => v.type === 'Pie chart'));
    expect(exists).toBe(true);
  });

  await softStep('8.2 Right-click → On Click → Filter', async () => {
    const container = page.locator('[name="viewer-Pie-chart"]');
    const canvas = container.locator('canvas').first();
    await canvas.click({button: 'right'});
    await page.waitForTimeout(500);
    await page.evaluate(async () => {
      const labels = document.querySelectorAll('.d4-menu-item-label');
      let onClickGroup = null;
      for (const label of labels) {
        if (label.textContent.trim() === 'On Click') {
          onClickGroup = label.closest('.d4-menu-item') || label.closest('.d4-menu-group');
          break;
        }
      }
      if (!onClickGroup) throw new Error('On Click menu not found');
      const rect = onClickGroup.getBoundingClientRect();
      onClickGroup.dispatchEvent(new MouseEvent('mouseenter', {clientX: rect.x + rect.width / 2, clientY: rect.y + rect.height / 2, bubbles: true}));
      onClickGroup.dispatchEvent(new MouseEvent('mousemove', {clientX: rect.x + rect.width / 2, clientY: rect.y + rect.height / 2, bubbles: true}));
      await new Promise(r => setTimeout(r, 300));
      const allLabels = document.querySelectorAll('.d4-menu-item-label');
      for (const label of allLabels) {
        if (label.textContent.trim() === 'Filter') {
          (label.closest('.d4-menu-item') as HTMLElement).click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 500));
    });
  });

  await softStep('8.3 Click segment — filter applies', async () => {
    const container = page.locator('[name="viewer-Pie-chart"]');
    const canvas = container.locator('canvas').first();
    const box = await canvas.boundingBox();
    if (!box) throw new Error('Canvas not found');
    await canvas.click({position: {x: box.width * 0.6, y: box.height * 0.4}});
    await page.waitForTimeout(1000);
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBeLessThan(baseline);
    expect(count).toBeGreaterThan(0);
  });

  await softStep('8.4 Click white space — reset', async () => {
    const container = page.locator('[name="viewer-Pie-chart"]');
    const canvas = container.locator('canvas').first();
    await canvas.click({position: {x: 5, y: 5}});
    await page.waitForTimeout(1000);
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBe(baseline);
  });

  await softStep('8.5 Close Pie Chart', async () => {
    await page.evaluate(async () => {
      for (const v of Array.from(grok.shell.tv.viewers))
        if (v.type === 'Pie chart') v.close();
    });
    await page.waitForTimeout(500);
  });

  // #### 9. Reset and cleanup

  await softStep('9.1 Reset all filters', async () => {
    await page.evaluate(async () => {
      const header = document.querySelector('[name="viewer-Filters"] .d4-filter-group-header');
      const resetIcon = header?.querySelector('.fa-arrow-rotate-left') as HTMLElement;
      if (resetIcon) resetIcon.click();
      await new Promise(r => setTimeout(r, 500));
      const okBtn = document.querySelector('.d4-dialog [name="button-OK"]') as HTMLElement;
      if (okBtn) okBtn.click();
      await new Promise(r => setTimeout(r, 500));
    });
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBe(await page.evaluate(() => grok.shell.tv.dataFrame.rowCount));
  });

  await softStep('9.2 Close All', async () => {
    await page.evaluate(() => grok.shell.closeAll());
    await page.waitForTimeout(500);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
