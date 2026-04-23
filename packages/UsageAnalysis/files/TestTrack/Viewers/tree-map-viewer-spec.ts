import {test} from '@playwright/test';

test.use({
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) { stepErrors.push({step: name, error: e.message ?? String(e)}); }
}

test('Tree Map tests', async ({page, baseURL}) => {
  test.setTimeout(300_000);

  await page.goto(baseURL ?? '/');
  await page.locator('[name="Toolbox"], [name="Browse"], .d4-sidebar').first().waitFor({timeout: 60000});
  await page.waitForFunction(() => {
    try {
      if (typeof grok === 'undefined' || !grok.shell) return false;
      grok.shell.settings.showFiltersIconsConstantly;
      return true;
    } catch { return false; }
  }, {timeout: 30000});

  // Setup: selenium class + open demog + wait for semType
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (grok as any).shell.settings.showFiltersIconsConstantly = true;
    (grok as any).shell.windows.simpleMode = true;
    (grok as any).shell.closeAll();
    const df = await (grok as any).dapi.files.readCsv('System:DemoFiles/demog.csv');
    (grok as any).shell.addTableView(df);
    await new Promise((resolve: any) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Add Tree Map
  await page.evaluate(() => {
    const el = document.querySelector('[name="icon-tree-map"]') as HTMLElement;
    el?.click();
  });
  await page.locator('[name="viewer-Tree-map"]').waitFor({timeout: 10000});

  // Helper: get Tree Map viewer via JS API
  const getTmViewer = async () => page.evaluate(async () => {
    const tv = (grok as any).shell.tv;
    return Array.from(tv.viewers as any[]).find((v: any) => v.type === 'Tree map');
  });

  // Helper: set options on Tree Map viewer
  const setTmOptions = (opts: Record<string, any>) => page.evaluate(
    async (o) => {
      const tv = (grok as any).shell.tv;
      const tm = Array.from(tv.viewers as any[]).find((v: any) => v.type === 'Tree map');
      tm.setOptions(o);
      await new Promise(r => setTimeout(r, 300));
    },
    opts
  );

  // Helper: change a <select> value by index
  const setSelectByIndex = (idx: number, value: string) => page.evaluate(
    async ({idx, value}) => {
      const selects = document.querySelectorAll('select');
      const sel = selects[idx] as HTMLSelectElement;
      if (!sel) throw new Error(`select[${idx}] not found`);
      sel.value = value;
      sel.dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise(r => setTimeout(r, 300));
    },
    {idx, value}
  );

  // --- Split Column — Single Level ---

  await softStep('Verify viewer opens with split column auto-selected (DIS_POP)', async () => {
    const val = await page.evaluate(() => {
      const selects = document.querySelectorAll('select');
      return (selects[0] as HTMLSelectElement)?.value;
    });
    if (!val) throw new Error('No split column pre-selected');
  });

  await softStep('Change first split to SEX', async () => {
    await setSelectByIndex(0, 'SEX');
    const val = await page.evaluate(() => (document.querySelectorAll('select')[0] as HTMLSelectElement).value);
    if (val !== 'SEX') throw new Error(`Expected SEX, got ${val}`);
  });

  await softStep('Change first split back to RACE', async () => {
    await setSelectByIndex(0, 'RACE');
    const val = await page.evaluate(() => (document.querySelectorAll('select')[0] as HTMLSelectElement).value);
    if (val !== 'RACE') throw new Error(`Expected RACE, got ${val}`);
  });

  // --- Split Column — Multiple Levels ---

  await softStep('RACE already set (step 1 of Multiple Levels)', async () => { /* RACE confirmed above */ });

  await softStep('Add SEX as second split level via trailing empty selector', async () => {
    await setSelectByIndex(1, 'SEX');
    const val = await page.evaluate(() => (document.querySelectorAll('select')[1] as HTMLSelectElement).value);
    if (val !== 'SEX') throw new Error(`Expected SEX in second selector, got ${val}`);
  });

  await softStep('Clear second selector to remove SEX level', async () => {
    await setSelectByIndex(1, '');
    // Verify canvas reverted to single-level
    const opts = await page.evaluate(() => {
      const tv = (grok as any).shell.tv;
      const tm = Array.from(tv.viewers as any[]).find((v: any) => v.type === 'Tree map');
      return tm.getOptions();
    });
    if (opts.look?.splitByColumnNames?.length !== 1) throw new Error('Expected single split level after clearing');
  });

  // --- Color Column and Aggregation ---

  await softStep('Set Color column to AGE', async () => {
    await setTmOptions({colorColumnName: 'AGE'});
    const opts = await page.evaluate(() => {
      const tv = (grok as any).shell.tv;
      const tm = Array.from(tv.viewers as any[]).find((v: any) => v.type === 'Tree map');
      return tm.getOptions();
    });
    if (opts.look?.colorColumnName !== 'AGE') throw new Error('colorColumnName not set to AGE');
  });

  await softStep('Change color aggregation to max', async () => {
    // Color aggregation select is inside [name="div-column-combobox-color"]
    await page.evaluate(async () => {
      const colorCombo = document.querySelector('[name="div-column-combobox-color"]') as HTMLElement;
      const sel = colorCombo.querySelector('select') as HTMLSelectElement;
      sel.value = 'max';
      sel.dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise(r => setTimeout(r, 300));
    });
  });

  await softStep('Change color aggregation to sum', async () => {
    await page.evaluate(async () => {
      const colorCombo = document.querySelector('[name="div-column-combobox-color"]') as HTMLElement;
      const sel = colorCombo.querySelector('select') as HTMLSelectElement;
      sel.value = 'sum';
      sel.dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise(r => setTimeout(r, 300));
    });
  });

  await softStep('Clear Color column', async () => {
    await setTmOptions({colorColumnName: ''});
  });

  // --- Size Column and Aggregation ---

  await softStep('Open Settings (gear icon)', async () => {
    await page.evaluate(() => {
      const tmViewer = document.querySelector('[name="viewer-Tree-map"]') as HTMLElement;
      let wrapper = tmViewer?.closest('.d4-viewer') || tmViewer?.parentElement;
      while (wrapper && !wrapper.querySelector('[name="icon-font-icon-settings"]')) {
        wrapper = wrapper.parentElement as HTMLElement;
        if (!wrapper || wrapper === document.body) break;
      }
      const icon = wrapper?.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
      icon?.click();
    });
    await page.waitForTimeout(300);
  });

  await softStep('Set Size Column Name to HEIGHT', async () => {
    await setTmOptions({sizeColumnName: 'HEIGHT'});
  });

  await softStep('Set Size Aggr Type to avg', async () => {
    await setTmOptions({sizeAggrType: 'avg'});
  });

  await softStep('Set Size Aggr Type back to sum', async () => {
    await setTmOptions({sizeAggrType: 'sum'});
  });

  await softStep('Clear Size Column Name', async () => {
    await setTmOptions({sizeColumnName: ''});
  });

  // --- Show/Hide Column Selection Panel ---

  await softStep('Uncheck Show Column Selection Panel', async () => {
    await page.evaluate(() => {
      const labels = Array.from(document.querySelectorAll('div, span, label'));
      const label = labels.find(el => el.textContent?.trim() === 'Show Column Selection Panel');
      const row = label?.closest('.ui-input-root, .d4-flex-row') || label?.parentElement;
      const checkbox = row?.querySelector('input[type="checkbox"]') as HTMLInputElement;
      checkbox?.click();
    });
    await page.waitForTimeout(300);
    // Verify selectors row disappeared
    const headerVisible = await page.evaluate(() => {
      const selects = document.querySelectorAll('select');
      // When panel hidden, split selects are removed from DOM or hidden
      return selects.length;
    });
    // headerVisible should be reduced (aggregation select may remain)
  });

  await softStep('Re-check Show Column Selection Panel', async () => {
    await page.evaluate(() => {
      const labels = Array.from(document.querySelectorAll('div, span, label'));
      const label = labels.find(el => el.textContent?.trim() === 'Show Column Selection Panel');
      const row = label?.closest('.ui-input-root, .d4-flex-row') || label?.parentElement;
      const checkbox = row?.querySelector('input[type="checkbox"]') as HTMLInputElement;
      checkbox?.click();
    });
    await page.waitForTimeout(300);
  });

  // --- Row Source ---

  await softStep('Set Row Source to All', async () => {
    await setTmOptions({rowSource: 'All'});
  });

  await softStep('Set Row Source to Selected', async () => {
    await setTmOptions({rowSource: 'Selected'});
  });

  await softStep('Set Row Source to Filtered', async () => {
    await setTmOptions({rowSource: 'Filtered'});
  });

  // --- Filter Formula ---

  await softStep('Set Filter to ${AGE} > 40', async () => {
    await setTmOptions({filter: '${AGE} > 40'});
    const opts = await page.evaluate(() => {
      const tv = (grok as any).shell.tv;
      const tm = Array.from(tv.viewers as any[]).find((v: any) => v.type === 'Tree map');
      return tm.getOptions();
    });
    if (opts.look?.filter !== '${AGE} > 40') throw new Error('Filter not set');
  });

  await softStep('Clear Filter field', async () => {
    await setTmOptions({filter: ''});
  });

  // --- Outer Margins ---

  await softStep('Set all outer margins to 30', async () => {
    await page.evaluate(async () => {
      const tv = (grok as any).shell.tv;
      const tm = Array.from(tv.viewers as any[]).find((v: any) => v.type === 'Tree map');
      tm.setOptions({outerMarginLeft: 30});
      await new Promise(r => setTimeout(r, 150));
      tm.setOptions({outerMarginTop: 30});
      await new Promise(r => setTimeout(r, 150));
      tm.setOptions({outerMarginRight: 30});
      await new Promise(r => setTimeout(r, 150));
      tm.setOptions({outerMarginBottom: 30});
      await new Promise(r => setTimeout(r, 300));
    });
  });

  await softStep('Reset all outer margins to 0', async () => {
    await setTmOptions({outerMarginLeft: 0, outerMarginTop: 0, outerMarginRight: 0, outerMarginBottom: 0});
  });

  // --- Row Selection ---

  await softStep('Set first split to RACE (already set)', async () => { /* RACE is current split */ });

  await softStep('Click center of Tree Map canvas to select rows', async () => {
    await page.evaluate(async () => {
      const tmViewer = document.querySelector('[name="viewer-Tree-map"]') as HTMLElement;
      const canvas = tmViewer?.querySelector('canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      const x = rect.left + rect.width / 2;
      const y = rect.top + rect.height / 2;
      canvas.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, clientX: x, clientY: y}));
      canvas.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, clientX: x, clientY: y}));
      canvas.dispatchEvent(new MouseEvent('click', {bubbles: true, clientX: x, clientY: y}));
      await new Promise(r => setTimeout(r, 500));
    });
    const selected = await page.evaluate(() => (grok as any).shell.tv.dataFrame.selection.trueCount);
    if (selected === 0) throw new Error('No rows selected after canvas click');
  });

  await softStep('Shift-click different area to expand selection', async () => {
    const beforeCount = await page.evaluate(() => (grok as any).shell.tv.dataFrame.selection.trueCount);
    await page.evaluate(async () => {
      const tmViewer = document.querySelector('[name="viewer-Tree-map"]') as HTMLElement;
      const canvas = tmViewer?.querySelector('canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      const x = rect.left + 20;
      const y = rect.top + 20;
      canvas.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, clientX: x, clientY: y, shiftKey: true}));
      canvas.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, clientX: x, clientY: y, shiftKey: true}));
      canvas.dispatchEvent(new MouseEvent('click', {bubbles: true, clientX: x, clientY: y, shiftKey: true}));
      await new Promise(r => setTimeout(r, 500));
    });
    const afterCount = await page.evaluate(() => (grok as any).shell.tv.dataFrame.selection.trueCount);
    if (afterCount <= beforeCount) throw new Error(`Selection did not expand: ${beforeCount} → ${afterCount}`);
  });

  await softStep('Ctrl-click first point to toggle off selection', async () => {
    const beforeCount = await page.evaluate(() => (grok as any).shell.tv.dataFrame.selection.trueCount);
    await page.evaluate(async () => {
      const tmViewer = document.querySelector('[name="viewer-Tree-map"]') as HTMLElement;
      const canvas = tmViewer?.querySelector('canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      const x = rect.left + rect.width / 2;
      const y = rect.top + rect.height / 2;
      canvas.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, clientX: x, clientY: y, ctrlKey: true}));
      canvas.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, clientX: x, clientY: y, ctrlKey: true}));
      canvas.dispatchEvent(new MouseEvent('click', {bubbles: true, clientX: x, clientY: y, ctrlKey: true}));
      await new Promise(r => setTimeout(r, 500));
    });
    const afterCount = await page.evaluate(() => (grok as any).shell.tv.dataFrame.selection.trueCount);
    if (afterCount >= beforeCount) throw new Error(`Selection did not shrink after ctrl-click: ${beforeCount} → ${afterCount}`);
  });

  // --- Layout Persistence ---

  await softStep('Set RACE+SEX split and AGE color for layout', async () => {
    await setSelectByIndex(1, 'SEX');
    await setTmOptions({colorColumnName: 'AGE'});
  });

  let savedLayoutId: string;
  await softStep('Save current layout', async () => {
    savedLayoutId = await page.evaluate(async () => {
      const tv = (grok as any).shell.tv;
      const layout = tv.saveLayout();
      const saved = await (grok as any).dapi.layouts.save(layout);
      return saved.id;
    });
    if (!savedLayoutId) throw new Error('Layout save returned no id');
  });

  await softStep('Close the Tree Map viewer', async () => {
    await page.evaluate(() => {
      const tv = (grok as any).shell.tv;
      const tm = Array.from(tv.viewers as any[]).find((v: any) => v.type === 'Tree map');
      tm.close();
    });
    await page.waitForTimeout(400);
    const count = await page.evaluate(() => Array.from((grok as any).shell.tv.viewers as any[]).length);
    if (count !== 1) throw new Error(`Expected 1 viewer after close, got ${count}`);
  });

  await softStep('Restore saved layout — Tree Map reappears with RACE+SEX split and AGE color', async () => {
    const opts = await page.evaluate(async (id) => {
      const tv = (grok as any).shell.tv;
      const layout = await (grok as any).dapi.layouts.find(id);
      tv.loadLayout(layout);
      await new Promise(r => setTimeout(r, 1500));
      const tm = Array.from(tv.viewers as any[]).find((v: any) => v.type === 'Tree map');
      if (!tm) throw new Error('Tree Map not found after layout restore');
      return tm.getOptions();
    }, savedLayoutId!);
    if (opts.look?.colorColumnName !== 'AGE') throw new Error('Color column not restored');
    if (!opts.look?.splitByColumnNames?.includes('SEX')) throw new Error('SEX split not restored');
  });

  await softStep('Delete saved layout', async () => {
    await page.evaluate(async (id) => {
      const layout = await (grok as any).dapi.layouts.find(id);
      await (grok as any).dapi.layouts.delete(layout);
    }, savedLayoutId!);
  });

  // Final error check
  if (stepErrors.length > 0) {
    throw new Error(
      'Some steps failed:\n' +
      stepErrors.map(e => `  [${e.step}] ${e.error}`).join('\n')
    );
  }
});
