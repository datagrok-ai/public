import {test, expect, Page} from '@playwright/test';

const baseUrl = process.env.BASE_URL ?? 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

/** Click a property-grid view label to open the inline editor, set a <select> value, press Enter */
async function setPropSelect(page: Page, rowName: string, value: string) {
  await page.evaluate((args: {rowName: string; value: string}) => {
    const row = document.querySelector(`[name="${args.rowName}"]`) as HTMLElement;
    (row?.querySelector('[name^="prop-view-"]') as HTMLElement)?.click();
    const sel = row?.querySelector('select.property-grid-item-editor-spinner') as HTMLSelectElement;
    if (sel) { sel.value = args.value; sel.dispatchEvent(new Event('change', {bubbles: true})); }
  }, {rowName, value});
  await page.keyboard.press('Enter');
}

/** Toggle a property-grid boolean checkbox row */
async function togglePropCheckbox(page: Page, rowName: string) {
  await page.evaluate((rn: string) => {
    (document.querySelector(`[name="${rn}"] input[type="checkbox"]`) as HTMLInputElement)?.click();
  }, rowName);
}

/** Open the density plot column popup via mousedown, type column name, confirm with Enter */
async function setColumnViaPopup(page: Page, axis: 'x' | 'y', colName: string) {
  await page.evaluate((ax: string) => {
    const viewer = document.querySelector('[name="viewer-Density-plot"]') as HTMLElement;
    const combo = viewer.querySelector(`[name="div-column-combobox-${ax}"]`) as HTMLElement;
    combo?.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
  }, axis);
  await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'), {timeout: 3000});
  await page.keyboard.type(colName);
  await page.keyboard.press('Enter');
}

test('Density plot tests', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  // Phase 1: Navigate
  await page.goto(baseUrl);
  await page.waitForFunction(() => {
    try { return typeof grok !== 'undefined' && grok.shell &&
      typeof grok.shell.settings?.showFiltersIconsConstantly === 'boolean'; }
    catch (e) { return false; }
  }, {timeout: 30000});

  // Phase 2: Open dataset — combine setup + open in one evaluate
  await page.evaluate(async (path: string) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = false;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Phase 3: Add Density Plot by clicking toolbox icon
  await page.evaluate(() => {
    (document.querySelector('[name="icon-density-plot"]') as HTMLElement)?.click();
  });
  await page.locator('[name="viewer-Density-plot"]').waitFor({timeout: 10000});

  // ── Zoom and Reset View ──────────────────────────────────────────────────
  await softStep('Zoom in with mouse wheel', async () => {
    const box = await page.locator('[name="viewer-Density-plot"] canvas').boundingBox();
    expect(box).toBeTruthy();
    await page.evaluate((b: {x: number; y: number; width: number; height: number}) => {
      const canvas = document.querySelector('[name="viewer-Density-plot"] canvas') as HTMLElement;
      canvas.dispatchEvent(new WheelEvent('wheel', {
        bubbles: true, cancelable: true,
        clientX: b.x + b.width / 2, clientY: b.y + b.height / 2, deltaY: -300,
      }));
    }, box!);
  });

  await softStep('Right-click → Reset View', async () => {
    const box = await page.locator('[name="viewer-Density-plot"] canvas').boundingBox();
    expect(box).toBeTruthy();
    await page.evaluate((b: {x: number; y: number; width: number; height: number}) => {
      const canvas = document.querySelector('[name="viewer-Density-plot"] canvas') as HTMLElement;
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true,
        clientX: b.x + b.width / 2, clientY: b.y + b.height / 2,
      }));
    }, box!);
    // Find Reset View menu item by text (context menu items have no name= attribute)
    await page.waitForFunction(() =>
      Array.from(document.querySelectorAll('[role="menuitem"]'))
        .some(el => el.textContent?.trim() === 'Reset View')
    , {timeout: 3000});
    await page.evaluate(() => {
      const item = Array.from(document.querySelectorAll('[role="menuitem"]'))
        .find(el => el.textContent?.trim() === 'Reset View') as HTMLElement;
      item?.click();
    });
  });

  // ── Axis Column Assignment ───────────────────────────────────────────────
  // Column selectors: [name="div-column-combobox-x"] / [name="div-column-combobox-y"] (lowercase)
  // Popup opens on mousedown; search by typing; confirm with Enter
  await softStep('Set X column to WEIGHT (JS API fallback — click on combo did not open popup)', async () => {
    // mousedown doesn't work before settings are opened; fall back to JS API
    await page.evaluate(() => {
      const dp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any;
      dp.props.xColumnName = 'WEIGHT';
    });
    const x = await page.evaluate(() =>
      (Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any).props.xColumnName
    );
    expect(x).toBe('WEIGHT');
  });

  await softStep('Set Y column to HEIGHT via UI popup (mousedown + type + Enter)', async () => {
    await setColumnViaPopup(page, 'y', 'HEIGHT');
    const y = await page.evaluate(() =>
      (Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any).props.yColumnName
    );
    expect(y).toBe('HEIGHT');
  });

  await softStep('Set X column to AGE via UI popup (mousedown + type + Enter)', async () => {
    await setColumnViaPopup(page, 'x', 'AGE');
    const x = await page.evaluate(() =>
      (Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any).props.xColumnName
    );
    expect(x).toBe('AGE');
  });

  await softStep('Set Y column to WEIGHT via UI popup (mousedown + type + Enter)', async () => {
    await setColumnViaPopup(page, 'y', 'WEIGHT');
    const y = await page.evaluate(() =>
      (Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any).props.yColumnName
    );
    expect(y).toBe('WEIGHT');
  });

  // ── Bins ─────────────────────────────────────────────────────────────────
  // Open settings: gear is 2 levels up from viewer-Density-plot container
  await softStep('Open settings panel', async () => {
    await page.evaluate(() => {
      const viewer = document.querySelector('[name="viewer-Density-plot"]') as HTMLElement;
      const gear = viewer.parentElement?.parentElement?.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
      gear?.click();
    });
  });

  // Expand Misc section (icon-plus → click row → icon-minus)
  await softStep('Expand Misc section', async () => {
    await page.evaluate(() => {
      const all = Array.from(document.querySelectorAll('*')) as HTMLElement[];
      const miscEl = all.find(el =>
        el.childNodes.length === 1 && el.childNodes[0].nodeType === 3 &&
        (el.childNodes[0] as Text).textContent?.trim() === 'Misc'
      );
      const row = miscEl?.closest('tr') as HTMLElement;
      const icon = row?.querySelector('.property-grid-category-icon');
      if (icon?.classList.contains('property-grid-icon-plus')) row?.click();
    });
    // Scroll Bins into view
    await page.evaluate(() => {
      document.querySelector('[name="prop-bins"]')?.scrollIntoView({behavior: 'instant', block: 'center'});
    });
  });

  // Bin Shape: click view label → <select class="property-grid-item-editor-spinner"> + Enter
  for (const shape of ['rectangle', 'hexagon'] as const) {
    await softStep(`Set Bin Shape to ${shape}`, async () => {
      await setPropSelect(page, 'prop-bin-shape', shape);
      const val = await page.evaluate(() =>
        (Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any).props.binShape
      );
      expect(val).toBe(shape);
    });
  }

  // Bins: focus input.property-grid-slider-textbox → Ctrl+A → type → Enter
  for (const bins of [5, 100, 200, 50]) {
    await softStep(`Set Bins to ${bins}`, async () => {
      await page.evaluate(() => {
        (document.querySelector('[name="prop-bins"] input.property-grid-slider-textbox') as HTMLInputElement)?.focus();
      });
      await page.keyboard.press('Control+A');
      await page.keyboard.type(String(bins));
      await page.keyboard.press('Enter');
      const val = await page.evaluate(() =>
        (Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any).props.bins
      );
      expect(val).toBe(bins);
    });
  }

  // ── Show/Hide Color Scale ─────────────────────────────────────────────────
  await softStep('Set Show Color Scale to false', async () => {
    await togglePropCheckbox(page, 'prop-show-color-scale');
    const val = await page.evaluate(() =>
      (Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any).props.showColorScale
    );
    expect(val).toBe(false);
  });
  await softStep('Set Show Color Scale to true', async () => {
    await togglePropCheckbox(page, 'prop-show-color-scale');
    const val = await page.evaluate(() =>
      (Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any).props.showColorScale
    );
    expect(val).toBe(true);
  });

  // ── Axis Visibility ───────────────────────────────────────────────────────
  await softStep('Set Show X Axis to false', async () => { await togglePropCheckbox(page, 'prop-show-x-axis'); });
  await softStep('Set Show Y Axis to false', async () => { await togglePropCheckbox(page, 'prop-show-y-axis'); });
  await softStep('Set Show X Axis to true', async () => { await togglePropCheckbox(page, 'prop-show-x-axis'); });
  await softStep('Set Show Y Axis to true', async () => {
    await togglePropCheckbox(page, 'prop-show-y-axis');
    const vals = await page.evaluate(() => {
      const dp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any;
      return {x: dp.props.showXAxis, y: dp.props.showYAxis};
    });
    expect(vals.x).toBe(true);
    expect(vals.y).toBe(true);
  });

  // ── Axis Inversion and Logarithmic Axis ───────────────────────────────────
  // Invert checkboxes: [name="prop-invert-x-axis"] / [name="prop-invert-y-axis"]
  // Axis type selects: [name="prop-x-axis-type"] / [name="prop-y-axis-type"]
  await softStep('Set Invert X Axis to true', async () => { await togglePropCheckbox(page, 'prop-invert-x-axis'); });
  await softStep('Set X Axis Type to logarithmic', async () => { await setPropSelect(page, 'prop-x-axis-type', 'logarithmic'); });
  await softStep('Set Invert Y Axis to true', async () => { await togglePropCheckbox(page, 'prop-invert-y-axis'); });
  await softStep('Set Y Axis Type to logarithmic', async () => { await setPropSelect(page, 'prop-y-axis-type', 'logarithmic'); });
  await softStep('Set Invert X Axis to false', async () => { await togglePropCheckbox(page, 'prop-invert-x-axis'); });
  await softStep('Set X Axis Type to linear', async () => { await setPropSelect(page, 'prop-x-axis-type', 'linear'); });
  await softStep('Set Invert Y Axis to false', async () => { await togglePropCheckbox(page, 'prop-invert-y-axis'); });
  await softStep('Set Y Axis Type to linear', async () => {
    await setPropSelect(page, 'prop-y-axis-type', 'linear');
    const vals = await page.evaluate(() => {
      const dp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any;
      return {invertX: dp.props.invertXAxis, invertY: dp.props.invertYAxis, xType: dp.props.xAxisType, yType: dp.props.yAxisType};
    });
    expect(vals.invertX).toBe(false);
    expect(vals.invertY).toBe(false);
    expect(vals.xType).toBe('linear');
    expect(vals.yType).toBe('linear');
  });

  // ── Show/Hide Selectors and Bin Slider ────────────────────────────────────
  await softStep('Set Show X Selector to false', async () => { await togglePropCheckbox(page, 'prop-show-x-selector'); });
  await softStep('Set Show Y Selector to false', async () => { await togglePropCheckbox(page, 'prop-show-y-selector'); });
  await softStep('Set Show Bin Selector to false', async () => { await togglePropCheckbox(page, 'prop-show-bin-selector'); });
  await softStep('Set Show X Selector to true', async () => { await togglePropCheckbox(page, 'prop-show-x-selector'); });
  await softStep('Set Show Y Selector to true', async () => { await togglePropCheckbox(page, 'prop-show-y-selector'); });
  await softStep('Set Show Bin Selector to true', async () => {
    await togglePropCheckbox(page, 'prop-show-bin-selector');
    const vals = await page.evaluate(() => {
      const dp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any;
      return {x: dp.props.showXSelector, y: dp.props.showYSelector, bin: dp.props.showBinSelector};
    });
    expect(vals.x).toBe(true);
    expect(vals.y).toBe(true);
    expect(vals.bin).toBe(true);
  });

  // ── Min/Max Axis Bounds ───────────────────────────────────────────────────
  // Min/Max inputs are nullable; property grid editor complex — use JS API
  await softStep('Set X Min=20, X Max=50, Y Min=100, Y Max=200 via JS API', async () => {
    const result = await page.evaluate(() => {
      const dp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any;
      dp.props['xMin'] = 20; dp.props['xMax'] = 50;
      dp.props['yMin'] = 100; dp.props['yMax'] = 200;
      return {xMin: dp.props['xMin'], xMax: dp.props['xMax'], yMin: dp.props['yMin'], yMax: dp.props['yMax']};
    });
    expect(result.xMin).toBe(20);
    expect(result.xMax).toBe(50);
    expect(result.yMin).toBe(100);
    expect(result.yMax).toBe(200);
  });
  await softStep('Clear X Min, X Max, Y Min, Y Max', async () => {
    const result = await page.evaluate(() => {
      const dp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any;
      dp.props['xMin'] = null; dp.props['xMax'] = null;
      dp.props['yMin'] = null; dp.props['yMax'] = null;
      return {xMin: dp.props['xMin'], xMax: dp.props['xMax'], yMin: dp.props['yMin'], yMax: dp.props['yMax']};
    });
    expect(result.xMin).toBeNull();
    expect(result.xMax).toBeNull();
    expect(result.yMin).toBeNull();
    expect(result.yMax).toBeNull();
  });

  // ── Color Transform Type ──────────────────────────────────────────────────
  // Row found by aria-label (not by name) since property is in Style category
  await softStep('Set Color Transform Type to logarithmic', async () => {
    await page.evaluate(() => {
      const row = Array.from(document.querySelectorAll('.property-grid-item'))
        .find(r => r.getAttribute('aria-label') === 'Color Transform Type') as HTMLElement;
      (row?.querySelector('[name^="prop-view-"]') as HTMLElement)?.click();
      const sel = row?.querySelector('select.property-grid-item-editor-spinner') as HTMLSelectElement;
      if (sel) { sel.value = 'logarithmic'; sel.dispatchEvent(new Event('change', {bubbles: true})); }
    });
    await page.keyboard.press('Enter');
  });
  await softStep('Set Color Transform Type to linear', async () => {
    await page.evaluate(() => {
      const row = Array.from(document.querySelectorAll('.property-grid-item'))
        .find(r => r.getAttribute('aria-label') === 'Color Transform Type') as HTMLElement;
      (row?.querySelector('[name^="prop-view-"]') as HTMLElement)?.click();
      const sel = row?.querySelector('select.property-grid-item-editor-spinner') as HTMLSelectElement;
      if (sel) { sel.value = 'linear'; sel.dispatchEvent(new Event('change', {bubbles: true})); }
    });
    await page.keyboard.press('Enter');
    const val = await page.evaluate(() =>
      (Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any).props.colorTransformType
    );
    expect(val).toBe('linear');
  });

  // ── Title and Description ─────────────────────────────────────────────────
  // Show Title: checkbox [name="prop-show-title"]
  // Title / Description: property grid text editor doesn't fire Dart events reliably → JS API
  // Visibility mode / position: [name="prop-description-visibility-mode"] / [name="prop-description-position"]
  await softStep('Enable Show Title', async () => {
    await page.evaluate(() => {
      const cb = document.querySelector('[name="prop-show-title"] input[type="checkbox"]') as HTMLInputElement;
      if (cb && !cb.checked) cb.click();
    });
  });
  await softStep('Set Title to "Density Distribution" via JS API', async () => {
    await page.evaluate(() => {
      (Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any).props.title = 'Density Distribution';
    });
  });
  await softStep('Set Description to "AGE vs HEIGHT density" via JS API', async () => {
    await page.evaluate(() => {
      (Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any).props.description = 'AGE vs HEIGHT density';
    });
    const vals = await page.evaluate(() => {
      const dp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any;
      return {title: dp.props.title, desc: dp.props.description};
    });
    expect(vals.title).toBe('Density Distribution');
    expect(vals.desc).toBe('AGE vs HEIGHT density');
  });
  await softStep('Set Description Visibility Mode to Always', async () => {
    await setPropSelect(page, 'prop-description-visibility-mode', 'Always');
    const val = await page.evaluate(() =>
      (Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any).props.descriptionVisibilityMode
    );
    expect(val).toBe('Always');
  });
  await softStep('Set Description Position to Bottom', async () => {
    await setPropSelect(page, 'prop-description-position', 'Bottom');
    const val = await page.evaluate(() =>
      (Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any).props.descriptionPosition
    );
    expect(val).toBe('Bottom');
  });

  // ── Selection ─────────────────────────────────────────────────────────────
  await softStep('Click a bin — verify rows selected', async () => {
    const box = await page.locator('[name="viewer-Density-plot"] canvas').boundingBox();
    expect(box).toBeTruthy();
    await page.evaluate((b: {x: number; y: number; width: number; height: number}) => {
      const canvas = document.querySelector('[name="viewer-Density-plot"] canvas') as HTMLElement;
      const cx = b.x + b.width * 0.4;
      const cy = b.y + b.height * 0.5;
      canvas.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, clientX: cx, clientY: cy, button: 0}));
      canvas.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, clientX: cx, clientY: cy, button: 0}));
      canvas.dispatchEvent(new MouseEvent('click', {bubbles: true, clientX: cx, clientY: cy, button: 0}));
    }, box!);
    const selected = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(selected).toBeGreaterThan(0);
  });
  await softStep('Esc — verify selection cleared', async () => {
    await page.keyboard.press('Escape');
    const selected = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(selected).toBe(0);
  });

  // ── Layout Persistence ────────────────────────────────────────────────────
  let savedLayoutId = '';
  await softStep('Set WEIGHT/HEIGHT/25 bins/rectangle/invertColor=true then save layout', async () => {
    await page.evaluate(() => {
      const dp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any;
      dp.props.xColumnName = 'WEIGHT';
      dp.props.yColumnName = 'HEIGHT';
      dp.props.bins = 25;
      dp.props.binShape = 'rectangle';
      dp.props.invertColorScheme = true;
    });
    savedLayoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      return layout.id;
    });
    expect(savedLayoutId).toBeTruthy();
  });

  await softStep('Close Density Plot viewer via JS API', async () => {
    // Close button [name="icon-times"] not present in title bar during testing; use dp.close()
    await page.evaluate(() => {
      (Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any)?.close();
    });
    await page.waitForFunction(
      () => !Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot'),
      {timeout: 5000}
    );
  });

  await softStep('Apply saved layout and verify properties restored', async () => {
    await page.evaluate(async (id: string) => {
      const saved = await grok.dapi.layouts.find(id);
      grok.shell.tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));
    }, savedLayoutId);
    const restored = await page.evaluate(() => {
      const dp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any;
      return {
        found: !!dp,
        x: dp?.props.xColumnName, y: dp?.props.yColumnName,
        bins: dp?.props.bins, binShape: dp?.props.binShape,
        invertColor: dp?.props.invertColorScheme,
      };
    });
    expect(restored.found).toBe(true);
    expect(restored.x).toBe('WEIGHT');
    expect(restored.y).toBe('HEIGHT');
    expect(restored.bins).toBe(25);
    expect(restored.binShape).toBe('rectangle');
    expect(restored.invertColor).toBe(true);
    // Cleanup
    await page.evaluate(async (id: string) => {
      const saved = await grok.dapi.layouts.find(id);
      if (saved) await grok.dapi.layouts.delete(saved);
    }, savedLayoutId);
  });

  // ── Row Source Filtering ──────────────────────────────────────────────────
  await softStep('Set Filter to ${AGE} > 30', async () => {
    await page.evaluate(() => {
      (Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any).props.filter = '${AGE} > 30';
    });
    const val = await page.evaluate(() =>
      (Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any).props.filter
    );
    expect(val).toBe('${AGE} > 30');
  });

  await softStep('Clear Filter', async () => {
    await page.evaluate(() => {
      (Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any).props.filter = '';
    });
  });

  await softStep('Open spgi-100', async () => {
    await page.evaluate(async (path: string) => {
      const df2 = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df2);
      await new Promise(resolve => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
        setTimeout(resolve, 5000);
      });
    }, spgiPath);
  });

  await softStep('Go back to demog view', async () => {
    // Switch to the view that contains the Density plot (demog view)
    await page.evaluate(() => {
      const demogView = Array.from(grok.shell.tableViews)
        .find((tv: any) => tv.dataFrame?.rowCount === 5850) as any;
      if (demogView) grok.shell.v = demogView;
    });
  });

  await softStep('Set Table to spgi-100 on density plot — no errors', async () => {
    const result = await page.evaluate(() => {
      const dp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any;
      // spgi-100 table is named "Table (2)" (100 rows) — identify by row count
      const spgi = Array.from(grok.shell.tables).find((t: any) => t.rowCount === 100) as any;
      if (!dp || !spgi) return {error: `dp=${!!dp}, spgi=${!!spgi}`};
      dp.props.table = spgi.name;
      return {tableSet: dp.props.table, spgiName: spgi.name};
    });
    expect(result.error).toBeUndefined();
  });

  await softStep('Set table back to demog', async () => {
    await page.evaluate(() => {
      const dp = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Density plot') as any;
      // demog table is named "Table" (5850 rows)
      const demog = Array.from(grok.shell.tables).find((t: any) => t.rowCount === 5850) as any;
      if (dp && demog) dp.props.table = demog.name;
    });
  });

  await softStep('Close All', async () => {
    await page.evaluate(() => grok.shell.closeAll());
  });

  // ── Final summary ─────────────────────────────────────────────────────────
  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
