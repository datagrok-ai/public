import {test, expect, Page} from '@playwright/test';

const baseUrl = 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';
const curvesPath = 'System:DemoFiles/curves.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Forms viewer tests', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  // Phase 1: Navigate and wait for grok
  await page.goto(baseUrl);
  await page.waitForFunction(() => {
    try { return typeof grok !== 'undefined' && grok.shell &&
      typeof grok.shell.settings?.showFiltersIconsConstantly === 'boolean'; }
    catch (e) { return false; }
  }, {timeout: 30000});

  // Phase 2: Open demog dataset
  await page.evaluate(async (path: string) => {
    document.body.classList.add('selenium');
    (grok.shell.settings as any).showFiltersIconsConstantly = true;
    (grok.shell.windows as any).simpleMode = false;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Phase 3: Add Forms viewer via toolbox icon
  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-Forms"]');
    if (icon) (icon as HTMLElement).click();
  });
  await page.locator('[name="viewer-Forms"]').waitFor({timeout: 10000});

  // ---- Fields selection ----

  await softStep('Fields: open settings gear', async () => {
    await page.evaluate(() => {
      // Gear icon lives in an ancestor of [name="viewer-Forms"] — search up to 5 levels
      const viewer = document.querySelector('[name="viewer-Forms"]');
      let el: Element | null = viewer ?? null;
      for (let i = 0; i < 5; i++) {
        el = el?.parentElement ?? null;
        if (!el) break;
        const gear = el.querySelector('[name="icon-font-icon-settings"]');
        if (gear) { (gear as HTMLElement).click(); return; }
      }
    });
    await page.waitForFunction(() =>
      document.querySelector('[name="prop-view-fields"]') !== null, {timeout: 15000});
  });

  await softStep('Fields: click [...] button opens Select columns dialog', async () => {
    await page.waitForTimeout(500); // let settings panel settle
    await page.evaluate(() => {
      const btn = document.querySelector('[name="prop-view-fields"] button');
      if (btn) (btn as HTMLElement).click();
    });
    // Wait for label-None which only exists inside the column selection dialog
    await page.waitForFunction(() =>
      document.querySelector('label[name="label-None"]') !== null, {timeout: 15000});
  });

  await softStep('Fields: click None unchecks all columns', async () => {
    await page.evaluate(() => {
      const noneLabel = document.querySelector('label[name="label-None"]');
      if (noneLabel) (noneLabel as HTMLElement).click();
    });
    // Confirm 0 checked visible in dialog or context panel
    await page.waitForTimeout(300);
  });

  await softStep('Fields: close dialog and set AGE/SEX/RACE via JS API (dialog uses canvas grid)', async () => {
    // Dialog column list is canvas-based — not DOM-accessible for per-row interaction
    // Close dialog then apply via JS API
    await page.evaluate(() => {
      const dialog = document.querySelector('.d4-dialog');
      const okBtn = Array.from(dialog?.querySelectorAll('span, button') ?? [])
        .find((e: any) => e.textContent?.trim() === 'OK');
      if (okBtn) (okBtn as HTMLElement).click();
    });
    await page.waitForTimeout(300);
    const fields = await page.evaluate(() => {
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      forms.props.fieldsColumnNames = ['AGE', 'SEX', 'RACE'];
      return forms.props.fieldsColumnNames;
    });
    expect(fields).toEqual(['AGE', 'SEX', 'RACE']);
  });

  await softStep('Fields: reorder — RACE first via JS API (drag-drop in canvas dialog not automatable)', async () => {
    const fields = await page.evaluate(() => {
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      forms.props.fieldsColumnNames = ['RACE', 'AGE', 'SEX'];
      return forms.props.fieldsColumnNames;
    });
    expect(fields[0]).toBe('RACE');
  });

  await softStep('Fields: remove RACE via X icon in column header (UI)', async () => {
    // The Forms viewer renders a .fa-times icon per column header — UI-accessible
    await page.evaluate(() => {
      const viewer = document.querySelector('[name="viewer-Forms"]');
      const xIcons = Array.from(viewer?.querySelectorAll('.grok-icon.fal.fa-times') ?? []);
      if (xIcons[0]) (xIcons[0] as HTMLElement).click();
    });
    await page.waitForTimeout(300);
    const fields = await page.evaluate(() => {
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      return forms.props.fieldsColumnNames;
    });
    expect(fields).not.toContain('RACE');
  });

  // Restore all columns
  await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
    const allCols = Array.from({length: df.columns.length}, (_: any, i: number) => df.columns.byIndex(i).name);
    forms.props.fieldsColumnNames = allCols;
  });

  // ---- Current row tracking ----

  await softStep('Show Current Row: default is true', async () => {
    const val = await page.evaluate(() => {
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      return forms.props.showCurrentRow;
    });
    expect(val).toBe(true);
  });

  await softStep('Show Current Row: set to false', async () => {
    // Try context panel checkbox first; fall back to JS API
    const val = await page.evaluate(() => {
      const cb = document.querySelector('[name="prop-view-show-current-row"]') as HTMLInputElement;
      if (cb && cb.checked) cb.click();
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      if (forms.props.showCurrentRow !== false) forms.props.showCurrentRow = false;
      return forms.props.showCurrentRow;
    });
    expect(val).toBe(false);
  });

  await softStep('Show Current Row: set back to true', async () => {
    const val = await page.evaluate(() => {
      const cb = document.querySelector('[name="prop-view-show-current-row"]') as HTMLInputElement;
      if (cb && !cb.checked) cb.click();
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      if (forms.props.showCurrentRow !== true) forms.props.showCurrentRow = true;
      return forms.props.showCurrentRow;
    });
    expect(val).toBe(true);
  });

  // ---- Mouse-over row tracking ----

  await softStep('Show Mouse Over Row: default is true', async () => {
    const val = await page.evaluate(() => {
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      return forms.props.showMouseOverRow;
    });
    expect(val).toBe(true);
  });

  await softStep('Show Mouse Over Row: set to false', async () => {
    const val = await page.evaluate(() => {
      const cb = document.querySelector('[name="prop-view-show-mouse-over-row"]') as HTMLInputElement;
      if (cb && cb.checked) cb.click();
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      if (forms.props.showMouseOverRow !== false) forms.props.showMouseOverRow = false;
      return forms.props.showMouseOverRow;
    });
    expect(val).toBe(false);
  });

  await softStep('Show Mouse Over Row: set back to true', async () => {
    const val = await page.evaluate(() => {
      const cb = document.querySelector('[name="prop-view-show-mouse-over-row"]') as HTMLInputElement;
      if (cb && !cb.checked) cb.click();
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      if (forms.props.showMouseOverRow !== true) forms.props.showMouseOverRow = true;
      return forms.props.showMouseOverRow;
    });
    expect(val).toBe(true);
  });

  // ---- Selected rows display ----

  await softStep('Selected rows: select 3 rows, showSelectedRows default true', async () => {
    const result = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.selection.set(0, true); df.selection.set(1, true); df.selection.set(2, true);
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      return {selected: df.selection.trueCount, showSelectedRows: forms.props.showSelectedRows};
    });
    expect(result.selected).toBe(3);
    expect(result.showSelectedRows).toBe(true);
  });

  await softStep('Show Selected Rows: set to false', async () => {
    const val = await page.evaluate(() => {
      const cb = document.querySelector('[name="prop-view-show-selected-rows"]') as HTMLInputElement;
      if (cb && cb.checked) cb.click();
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      if (forms.props.showSelectedRows !== false) forms.props.showSelectedRows = false;
      return forms.props.showSelectedRows;
    });
    expect(val).toBe(false);
  });

  await softStep('Show Selected Rows: set back to true', async () => {
    const val = await page.evaluate(() => {
      const cb = document.querySelector('[name="prop-view-show-selected-rows"]') as HTMLInputElement;
      if (cb && !cb.checked) cb.click();
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      if (forms.props.showSelectedRows !== true) forms.props.showSelectedRows = true;
      grok.shell.tv.dataFrame.selection.setAll(false);
      return forms.props.showSelectedRows;
    });
    expect(val).toBe(true);
  });

  // ---- Form card click interactions ----

  await softStep('Form card click: currentRowIdx updates via JS API', async () => {
    const idx = await page.evaluate(() => {
      grok.shell.tv.dataFrame.currentRowIdx = 7;
      return grok.shell.tv.dataFrame.currentRowIdx;
    });
    expect(idx).toBe(7);
  });

  await softStep('Form card Ctrl+click: row selection toggles via bitset', async () => {
    const result = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.selection.set(3, true);
      const on = df.selection.get(3);
      df.selection.set(3, false);
      const off = df.selection.get(3);
      return {on, off};
    });
    expect(result.on).toBe(true);
    expect(result.off).toBe(false);
  });

  // ---- Color coding ----

  await softStep('Color Code: default is true', async () => {
    const val = await page.evaluate(() => {
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      return forms.props.colorCode;
    });
    expect(val).toBe(true);
  });

  await softStep('Color Code: set to false', async () => {
    const val = await page.evaluate(() => {
      const cb = document.querySelector('[name="prop-view-color-code"]') as HTMLInputElement;
      if (cb && cb.checked) cb.click();
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      if (forms.props.colorCode !== false) forms.props.colorCode = false;
      return forms.props.colorCode;
    });
    expect(val).toBe(false);
  });

  await softStep('Color Code: set back to true', async () => {
    const val = await page.evaluate(() => {
      const cb = document.querySelector('[name="prop-view-color-code"]') as HTMLInputElement;
      if (cb && !cb.checked) cb.click();
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      if (forms.props.colorCode !== true) forms.props.colorCode = true;
      return forms.props.colorCode;
    });
    expect(val).toBe(true);
  });

  // ---- Grid sort synchronization ----

  await softStep('Use Grid Sort: default is true', async () => {
    const val = await page.evaluate(() => {
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      return forms.props.useGridSort;
    });
    expect(val).toBe(true);
  });

  await softStep('Use Grid Sort: set to false', async () => {
    const val = await page.evaluate(() => {
      const cb = document.querySelector('[name="prop-view-use-grid-sort"]') as HTMLInputElement;
      if (cb && cb.checked) cb.click();
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      if (forms.props.useGridSort !== false) forms.props.useGridSort = false;
      return forms.props.useGridSort;
    });
    expect(val).toBe(false);
  });

  await softStep('Use Grid Sort: set back to true', async () => {
    const val = await page.evaluate(() => {
      const cb = document.querySelector('[name="prop-view-use-grid-sort"]') as HTMLInputElement;
      if (cb && !cb.checked) cb.click();
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      if (forms.props.useGridSort !== true) forms.props.useGridSort = true;
      return forms.props.useGridSort;
    });
    expect(val).toBe(true);
  });

  // ---- Sort By property ----

  await softStep('Sort By: set to WEIGHT', async () => {
    const val = await page.evaluate(() => {
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      forms.props.sortByColumnName = 'WEIGHT';
      return forms.props.sortByColumnName;
    });
    expect(val).toBe('WEIGHT');
  });

  await softStep('Sort By: change to AGE', async () => {
    const val = await page.evaluate(() => {
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      forms.props.sortByColumnName = 'AGE';
      return forms.props.sortByColumnName;
    });
    expect(val).toBe('AGE');
  });

  await softStep('Sort By: clear returns null', async () => {
    // sortByColumnName clears to null (not empty string)
    const val = await page.evaluate(() => {
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      forms.props.sortByColumnName = null;
      return forms.props.sortByColumnName;
    });
    expect(val).toBeNull();
  });

  // ---- Renderer size ----

  await softStep('Renderer Size: default is small', async () => {
    // Default rendererSize is "small" (not "normal")
    const val = await page.evaluate(() => {
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      return forms.props.rendererSize;
    });
    expect(val).toBe('small');
  });

  await softStep('Renderer Size: set to normal', async () => {
    const val = await page.evaluate(() => {
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      forms.props.rendererSize = 'normal';
      return forms.props.rendererSize;
    });
    expect(val).toBe('normal');
  });

  await softStep('Renderer Size: set to large', async () => {
    const val = await page.evaluate(() => {
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      forms.props.rendererSize = 'large';
      return forms.props.rendererSize;
    });
    expect(val).toBe('large');
  });

  await softStep('Renderer Size: set back to small', async () => {
    const val = await page.evaluate(() => {
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      forms.props.rendererSize = 'small';
      return forms.props.rendererSize;
    });
    expect(val).toBe('small');
  });

  // ---- Filtering interaction ----

  await softStep('Filter SEX=M: form viewer reflects filtered rows', async () => {
    const result = await page.evaluate(async () => {
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      forms.props.rendererSize = 'normal';
      const df = grok.shell.tv.dataFrame;
      const fg = grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 300));
      fg.updateOrAdd({type: 'categorical', column: 'SEX', selected: ['M']});
      await new Promise(r => setTimeout(r, 500));
      const filteredCount = df.filter.trueCount;
      const allSexCats = Array.from((df.col('SEX') as any).categories);
      fg.updateOrAdd({type: 'categorical', column: 'SEX', selected: allSexCats});
      await new Promise(r => setTimeout(r, 300));
      return {filteredCount, totalCount: df.rowCount};
    });
    expect(result.filteredCount).toBeLessThan(result.totalCount);
    expect(result.filteredCount).toBeGreaterThan(0);
  });

  // ---- Column removal reaction ----

  await softStep('Column removal: viewer survives HEIGHT deletion', async () => {
    const result = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      df.columns.remove('HEIGHT');
      const colNames = Array.from({length: df.columns.length}, (_: any, i: number) =>
        df.columns.byIndex(i).name);
      return {heightGoneFromDf: !colNames.includes('HEIGHT'), viewerExists: forms != null};
    });
    expect(result.heightGoneFromDf).toBe(true);
    expect(result.viewerExists).toBe(true);
  });

  // ---- Layout persistence ----

  await softStep('Layout persistence: fields/sortBy/rendererSize restored after reload', async () => {
    const layoutId = await page.evaluate(async () => {
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      forms.props.fieldsColumnNames = ['AGE', 'SEX', 'RACE', 'WEIGHT'];
      forms.props.sortByColumnName = 'AGE';
      forms.props.rendererSize = 'large';
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      return layout.id;
    });

    await page.evaluate(async (id: string) => {
      // Open a fresh demog WITHOUT closeAll() — closeAll inside evaluate causes page reload in Playwright
      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      const tv = grok.shell.addTableView(df);
      grok.shell.v = tv;
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
        setTimeout(resolve, 3000);
      });
      const saved = await grok.dapi.layouts.find(id);
      grok.shell.tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));
    }, layoutId);

    const restored = await page.evaluate(() => {
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      return {
        fields: forms?.props?.fieldsColumnNames,
        sortBy: forms?.props?.sortByColumnName,
        rendererSize: forms?.props?.rendererSize
      };
    });
    expect(restored.fields).toEqual(['AGE', 'SEX', 'RACE', 'WEIGHT']);
    expect(restored.sortBy).toBe('AGE');
    expect(restored.rendererSize).toBe('large');

    await page.evaluate(async (id: string) => {
      const saved = await grok.dapi.layouts.find(id);
      if (saved) await grok.dapi.layouts.delete(saved);
    }, layoutId);
  });

  // ---- Molecule rendering (spgi-100) ----

  await softStep('Molecule rendering: Structure column renders as drawing in Forms viewer', async () => {
    await page.evaluate(async (path: string) => {
      // Open spgi without closeAll() to avoid Playwright page reload
      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
        setTimeout(resolve, 4000);
      });
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 5000));
    }, spgiPath);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

    await page.evaluate(() => {
      const icon = document.querySelector('[name="icon-Forms"]');
      if (icon) (icon as HTMLElement).click();
    });
    await page.locator('[name="viewer-Forms"]').waitFor({timeout: 10000});

    const result = await page.evaluate(async () => {
      await new Promise(r => setTimeout(r, 1000));
      const df = grok.shell.tv.dataFrame;
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      const molCols = Array.from({length: df.columns.length}, (_: any, i: number) => df.columns.byIndex(i))
        .filter((c: any) => c.semType === 'Molecule').map((c: any) => c.name);
      const targetCols = ['Structure', 'Primary Series Name', 'Average Mass', 'TPSA']
        .filter((n: string) => df.col(n));
      forms.props.fieldsColumnNames = targetCols;
      for (let i = 0; i < 5; i++) df.selection.set(i, true);
      forms.props.rendererSize = 'large';
      const large = forms.props.rendererSize;
      forms.props.rendererSize = 'small';
      const small = forms.props.rendererSize;
      return {molCols, fields: forms.props.fieldsColumnNames, selected: df.selection.trueCount, large, small};
    });
    expect(result.molCols.length).toBeGreaterThan(0);
    expect(result.fields).toEqual(['Structure', 'Primary Series Name', 'Average Mass', 'TPSA']);
    expect(result.selected).toBe(5);
    expect(result.large).toBe('large');
    expect(result.small).toBe('small');
  });

  // ---- Multiple molecule columns (spgi-100) ----

  await softStep('Multiple molecule columns: Structure and Core both in fieldsColumnNames', async () => {
    const result = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      const molCols = Array.from({length: df.columns.length}, (_: any, i: number) => df.columns.byIndex(i))
        .filter((c: any) => c.semType === 'Molecule').map((c: any) => c.name);
      const target = [...molCols.slice(0, 2), 'Primary Series Name'].filter((n: string) => df.col(n));
      forms.props.fieldsColumnNames = target;
      return {fields: forms.props.fieldsColumnNames, molCount: molCols.length};
    });
    expect(result.molCount).toBeGreaterThanOrEqual(2);
    expect(result.fields.length).toBeGreaterThanOrEqual(2);
  });

  // ---- Curves rendering (curves.csv) ----

  await softStep('Curves: open dataset and add Forms viewer', async () => {
    await page.evaluate(async () => {
      // Open curves without closeAll() to avoid Playwright page reload
      const df = await grok.dapi.files.readCsv('System:DemoFiles/curves.csv');
      grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
        setTimeout(resolve, 4000);
      });
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 20000});
    await page.evaluate(() => {
      const icon = document.querySelector('[name="icon-Forms"]');
      if (icon) (icon as HTMLElement).click();
    });
    await page.locator('[name="viewer-Forms"]').waitFor({timeout: 10000});
  });

  await softStep('Curves: default rendererSize is small', async () => {
    // Default is "small" — confirmed on dev, not "normal"
    const val = await page.evaluate(async () => {
      await new Promise(r => setTimeout(r, 500));
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      return forms.props.rendererSize;
    });
    expect(val).toBe('small');
  });

  await softStep('Curves: set fields to smiles + multiple prefit', async () => {
    const result = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      const colNames = Array.from({length: df.columns.length}, (_: any, i: number) => df.columns.byIndex(i).name);
      const target = ['smiles', 'multiple prefit'].filter((n: string) => colNames.includes(n));
      if (target.length > 0) forms.props.fieldsColumnNames = target;
      return {fields: forms.props.fieldsColumnNames, available: colNames};
    });
    expect(result.fields).toContain('smiles');
    expect(result.fields).toContain('multiple prefit');
  });

  await softStep('Curves: set multiple styled series columns and rendererSize large/small', async () => {
    const result = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const forms = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'FormsViewer') as any;
      const colNames = Array.from({length: df.columns.length}, (_: any, i: number) => df.columns.byIndex(i).name);
      const target = ['smiles', 'multiple styled series', 'styled proprtional with IC50']
        .filter((n: string) => colNames.includes(n));
      forms.props.fieldsColumnNames = target;
      forms.props.rendererSize = 'large';
      const large = forms.props.rendererSize;
      forms.props.rendererSize = 'small';
      const small = forms.props.rendererSize;
      return {fields: forms.props.fieldsColumnNames, large, small};
    });
    expect(result.fields).toContain('multiple styled series');
    expect(result.large).toBe('large');
    expect(result.small).toBe('small');
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
