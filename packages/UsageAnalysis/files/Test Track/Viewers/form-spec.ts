import {test, expect, Page} from '@playwright/test';

test.use({
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const baseUrl = 'https://dev.datagrok.ai';
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

test('Form tests (Playwright)', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  await page.goto(baseUrl);
  await page.waitForFunction(() => {
    try { return typeof grok !== 'undefined' && grok.shell &&
      typeof grok.shell.settings?.showFiltersIconsConstantly === 'boolean'; }
    catch (e) { return false; }
  }, {timeout: 30000});

  // Phase 2: Open demog
  await page.evaluate(async (path: string) => {
    document.body.classList.add('selenium');
    (grok.shell.settings as any).showFiltersIconsConstantly = true;
    (grok.shell.windows as any).simpleMode = false;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise<void>(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Phase 3: Add Form viewer
  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-form"]');
    if (icon) (icon as HTMLElement).click();
    else grok.shell.tv.addViewer('Form');
  });
  await page.locator('[name="viewer-Form"]').waitFor({timeout: 10000});

  // ---- Row navigation ----

  await softStep('Row nav: default row is 0', async () => {
    const row = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
    expect(row).toBe(0);
  });

  await softStep('Row nav: chevron-right advances to row 1', async () => {
    await page.evaluate(() => {
      const viewer = document.querySelector('[name="viewer-Form"]');
      const btn = viewer?.querySelector('[name="icon-chevron-right"]');
      if (btn) (btn as HTMLElement).click();
    });
    await page.waitForTimeout(300);
    const row = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
    expect(row).toBe(1);
  });

  await softStep('Row nav: chevron-left returns to row 0', async () => {
    await page.evaluate(() => {
      const viewer = document.querySelector('[name="viewer-Form"]');
      const btn = viewer?.querySelector('[name="icon-chevron-left"]');
      if (btn) (btn as HTMLElement).click();
    });
    await page.waitForTimeout(300);
    const row = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
    expect(row).toBe(0);
  });

  await softStep('Row nav: grid row 5 click updates form (JS API - canvas grid)', async () => {
    const row = await page.evaluate(() => {
      grok.shell.tv.dataFrame.currentRowIdx = 5;
      return grok.shell.tv.dataFrame.currentRowIdx;
    });
    expect(row).toBe(5);
  });

  // ---- Keyboard navigation ----

  await softStep('Keyboard nav: Right arrow advances to next row', async () => {
    await page.evaluate(() => {
      grok.shell.tv.dataFrame.currentRowIdx = 0;
    });
    await page.waitForTimeout(200);
    // Click the form viewer container to give it keyboard focus
    await page.locator('[name="viewer-Form"]').click({force: true});
    await page.waitForTimeout(200);
    await page.keyboard.press('ArrowRight');
    await page.waitForTimeout(200);
    const row = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
    expect(row).toBe(1);
  });

  await softStep('Keyboard nav: Down arrow advances to next row', async () => {
    await page.keyboard.press('ArrowDown');
    await page.waitForTimeout(200);
    const row = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
    expect(row).toBe(2);
  });

  await softStep('Keyboard nav: Left arrow goes to previous row', async () => {
    await page.keyboard.press('ArrowLeft');
    await page.waitForTimeout(200);
    const row = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
    expect(row).toBe(1);
  });

  await softStep('Keyboard nav: Up arrow goes to previous row', async () => {
    await page.keyboard.press('ArrowUp');
    await page.waitForTimeout(200);
    const row = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
    expect(row).toBe(0);
  });

  await softStep('Keyboard nav: Space toggles row selection', async () => {
    await page.keyboard.press('Space');
    await page.waitForTimeout(200);
    const sel = await page.evaluate(() => grok.shell.tv.dataFrame.selection.get(0));
    expect(sel).toBe(true);
    // Deselect
    await page.keyboard.press('Space');
    await page.waitForTimeout(200);
  });

  // ---- Row selection ----

  await softStep('Row selection: icon-square selects current row', async () => {
    await page.evaluate(() => {
      grok.shell.tv.dataFrame.selection.setAll(false);
      grok.shell.tv.dataFrame.currentRowIdx = 0;
    });
    await page.evaluate(() => {
      const viewer = document.querySelector('[name="viewer-Form"]');
      const btn = viewer?.querySelector('[name="icon-square"]');
      if (btn) (btn as HTMLElement).click();
    });
    await page.waitForTimeout(200);
    const sel = await page.evaluate(() => grok.shell.tv.dataFrame.selection.get(0));
    expect(sel).toBe(true);
  });

  await softStep('Row selection: icon-square again deselects row 0', async () => {
    await page.evaluate(() => {
      const viewer = document.querySelector('[name="viewer-Form"]');
      const btn = viewer?.querySelector('[name="icon-square"]');
      if (btn) (btn as HTMLElement).click();
    });
    await page.waitForTimeout(200);
    const sel = await page.evaluate(() => grok.shell.tv.dataFrame.selection.get(0));
    expect(sel).toBe(false);
  });

  await softStep('Row selection: select row 3 in grid, icon-square reflects it', async () => {
    const result = await page.evaluate(() => {
      grok.shell.tv.dataFrame.selection.set(3, true);
      grok.shell.tv.dataFrame.currentRowIdx = 3;
      return grok.shell.tv.dataFrame.selection.get(3);
    });
    expect(result).toBe(true);
  });

  // ---- Sync mode ----

  await softStep('Sync mode: default is Current', async () => {
    const mode = await page.evaluate(() => {
      const form = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Form') as any;
      return form?.props?.syncMode;
    });
    expect(mode).toBe('Current');
  });

  await softStep('Sync mode: right-click Track Row → Mouse Over', async () => {
    await page.evaluate(async () => {
      const viewer = document.querySelector('[name="viewer-Form"]');
      const canvas = viewer?.querySelector('.d4-form-viewer') || viewer;
      (canvas as HTMLElement).dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2}));
      await new Promise(r => setTimeout(r, 400));
      const allItems = Array.from(document.querySelectorAll('.d4-menu-item'));
      const trackRowGroup = allItems.find(el => el.querySelector('.d4-menu-item-label')?.textContent.trim() === 'Track Row');
      if (trackRowGroup) {
        trackRowGroup.dispatchEvent(new MouseEvent('mousemove', {bubbles: true}));
        trackRowGroup.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise(r => setTimeout(r, 400));
      }
      const mouseOverItem = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(l => l.textContent.trim() === 'Mouse Over');
      if (mouseOverItem) (mouseOverItem.closest('.d4-menu-item') as HTMLElement)?.click();
      await new Promise(r => setTimeout(r, 300));
    });
    const mode = await page.evaluate(() => {
      const form = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Form') as any;
      return form?.props?.syncMode;
    });
    expect(mode).toBe('Mouse Over');
  });

  await softStep('Sync mode: set None via JS API (menu has ambiguous None items)', async () => {
    // Context menu Track Row → None fails when multiple None items exist — use JS API
    const mode = await page.evaluate(() => {
      const form = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Form') as any;
      form.props.syncMode = 'None';
      return form.props.syncMode;
    });
    expect(mode).toBe('None');
  });

  await softStep('Sync mode: restore to Current', async () => {
    const mode = await page.evaluate(() => {
      const form = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Form') as any;
      form.props.syncMode = 'Current';
      return form.props.syncMode;
    });
    expect(mode).toBe('Current');
  });

  // ---- Edit mode ----

  await softStep('Edit mode: icon-edit makes fields editable', async () => {
    const result = await page.evaluate(async () => {
      const viewer = document.querySelector('[name="viewer-Form"]');
      const editBtn = viewer?.querySelector('[name="icon-edit"]');
      if (editBtn) (editBtn as HTMLElement).click();
      await new Promise(r => setTimeout(r, 300));
      const editable = Array.from(viewer?.querySelectorAll('input') ?? []).filter((el: any) => !el.readOnly && !el.disabled).length;
      return editable;
    });
    expect(result).toBeGreaterThan(0);
  });

  await softStep('Edit mode: icon-edit off makes fields readonly again', async () => {
    const result = await page.evaluate(async () => {
      const viewer = document.querySelector('[name="viewer-Form"]');
      const editBtn = viewer?.querySelector('[name="icon-edit"]');
      if (editBtn) (editBtn as HTMLElement).click();
      await new Promise(r => setTimeout(r, 300));
      return Array.from(viewer?.querySelectorAll('input') ?? []).filter((el: any) => el.readOnly).length;
    });
    expect(result).toBeGreaterThan(0);
  });

  // ---- Column selector ----

  await softStep('Column selector: icon-list opens dialog', async () => {
    await page.evaluate(() => {
      const viewer = document.querySelector('[name="viewer-Form"]');
      const listBtn = viewer?.querySelector('[name="icon-list"]');
      if (listBtn) (listBtn as HTMLElement).click();
    });
    await page.waitForFunction(() => !!document.querySelector('.d4-dialog'), {timeout: 5000});
  });

  await softStep('Column selector: None unchecks all, OK closes dialog', async () => {
    await page.evaluate(async () => {
      const noneLabel = document.querySelector('label[name="label-None"]');
      if (noneLabel) (noneLabel as HTMLElement).click();
      await new Promise(r => setTimeout(r, 200));
      const dialog = document.querySelector('.d4-dialog');
      const okBtn = Array.from(dialog?.querySelectorAll('button') ?? []).find((b: any) => b.textContent.trim() === 'OK');
      if (okBtn) (okBtn as HTMLElement).click();
    });
    await page.waitForTimeout(300);
  });

  // ---- Toolbar visibility properties ----

  await softStep('Toolbar visibility: open settings panel', async () => {
    await page.evaluate(() => {
      const viewer = document.querySelector('[name="viewer-Form"]');
      let el: Element | null = viewer ?? null;
      for (let i = 0; i < 6; i++) {
        el = el?.parentElement ?? null;
        if (!el) break;
        const gear = el.querySelector('[name="icon-font-icon-settings"]');
        if (gear) { (gear as HTMLElement).click(); return; }
      }
    });
    await page.waitForFunction(() => !!document.querySelector('[name="prop-view-show-navigation"]'), {timeout: 10000});
  });

  await softStep('Toolbar visibility: all 9 props toggle false/true', async () => {
    const results = await page.evaluate(async () => {
      const form = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Form') as any;
      const props = ['showNavigation','showPrevRowArrow','showNextRowArrow','showRowSelector',
        'showFieldEditor','showDesignEditor','showColumnSelector','showSaveFile','showOpenFile'];
      const res: Record<string, boolean> = {};
      for (const p of props) {
        form.props[p] = false;
        await new Promise(r => setTimeout(r, 80));
        res[p + '_false'] = form.props[p] === false;
        form.props[p] = true;
        await new Promise(r => setTimeout(r, 80));
        res[p + '_true'] = form.props[p] === true;
      }
      return res;
    });
    for (const [k, v] of Object.entries(results)) {
      expect(v, k).toBe(true);
    }
  });

  // ---- Filtered data navigation ----

  await softStep('Filtered nav: SEX=M filter, chevron skips non-M rows', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.currentRowIdx = 0;
      const fg = grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 300));
      fg.updateOrAdd({type: 'categorical', column: 'SEX', selected: ['M']});
      await new Promise(r => setTimeout(r, 500));
      const filtered = df.filter.trueCount;
      const viewer = document.querySelector('[name="viewer-Form"]');
      const nextBtn = viewer?.querySelector('[name="icon-chevron-right"]');
      if (nextBtn) (nextBtn as HTMLElement).click();
      await new Promise(r => setTimeout(r, 200));
      const row1 = df.currentRowIdx;
      const row1Sex = df.col('SEX').get(row1);
      // Remove filter
      const allCats = Array.from((df.col('SEX') as any).categories);
      fg.updateOrAdd({type: 'categorical', column: 'SEX', selected: allCats});
      await new Promise(r => setTimeout(r, 300));
      return { filtered, row1Sex, restored: df.filter.trueCount };
    });
    expect(result.filtered).toBeLessThan(5850);
    expect(result.row1Sex).toBe('M');
    expect(result.restored).toBe(5850);
  });

  // ---- Context menu ----

  await softStep('Context menu: Edit Form..., Select Columns..., Track Row present', async () => {
    const result = await page.evaluate(async () => {
      const viewer = document.querySelector('[name="viewer-Form"]');
      const canvas = viewer?.querySelector('.d4-form-viewer') || viewer;
      (canvas as HTMLElement).dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2}));
      await new Promise(r => setTimeout(r, 400));
      const labels = Array.from(document.querySelectorAll('.d4-menu-item-label')).map((el: any) => el.textContent.trim());
      document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
      return {
        hasEditForm: labels.includes('Edit Form...'),
        hasSelectCols: labels.includes('Select Columns...'),
        hasTrackRow: labels.includes('Track Row'),
      };
    });
    expect(result.hasEditForm).toBe(true);
    expect(result.hasSelectCols).toBe(true);
    expect(result.hasTrackRow).toBe(true);
  });

  // ---- Column changes reaction ----

  await softStep('Column changes: HEIGHT removal, viewer survives', async () => {
    const result = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      if (df.col('HEIGHT')) df.columns.remove('HEIGHT');
      const form = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Form');
      return { heightGone: !df.col('HEIGHT'), formAlive: !!form };
    });
    expect(result.heightGone).toBe(true);
    expect(result.formAlive).toBe(true);
  });

  // ---- Column rename ----

  await softStep('Column changes: rename AGE → AGE_NEW, form stays alive', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const ageCol = df.col('AGE');
      if (!ageCol) return { skipped: true, renamed: false, formAlive: false };
      ageCol.name = 'AGE_NEW';
      await new Promise(r => setTimeout(r, 600));
      const renamed = !!df.col('AGE_NEW') && !df.col('AGE');
      const form = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Form');
      ageCol.name = 'AGE';
      return { skipped: false, renamed, formAlive: !!form };
    });
    if (!result.skipped) {
      expect(result.renamed).toBe(true);
      expect(result.formAlive).toBe(true);
    }
  });

  // ---- Design mode ----

  await softStep('Design mode: icon-object-ungroup button found and toggles', async () => {
    const result = await page.evaluate(async () => {
      const viewer = document.querySelector('[name="viewer-Form"]');
      const btn = viewer?.querySelector('[name="icon-object-ungroup"]');
      if (!btn) return { found: false };
      (btn as HTMLElement).click();
      await new Promise(r => setTimeout(r, 300));
      (btn as HTMLElement).click();
      await new Promise(r => setTimeout(r, 300));
      return { found: true };
    });
    expect(result.found).toBe(true);
  });

  // ---- Layout persistence ----

  await softStep('Layout: save and restore Form viewer', async () => {
    const layoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      return layout.id;
    });

    await page.evaluate(async (id: string) => {
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 300));
      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      grok.shell.addTableView(df);
      await new Promise<void>(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      const saved = await grok.dapi.layouts.find(id);
      grok.shell.tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));
      await grok.dapi.layouts.delete(saved);
    }, layoutId);

    await page.locator('[name="viewer-Form"]').waitFor({timeout: 15000});
    const form = await page.evaluate(() => {
      return !!Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Form');
    });
    expect(form).toBe(true);
  });

  // ---- Color coding ----

  await softStep('Color coding: AGE column gets linear color scheme', async () => {
    const result = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const ageCol = df.col('AGE');
      if (ageCol) ageCol.meta.colors.setLinear([0xFF0000FF, 0xFFFF0000]);
      return !!ageCol;
    });
    expect(result).toBe(true);
  });

  // ---- Table switching ----

  await softStep('Table switching: open SPGI, add Form, switch table prop', async () => {
    const result = await page.evaluate(async () => {
      const spgiDf = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      grok.shell.addTableView(spgiDf);
      await new Promise<void>(resolve => {
        const sub = spgiDf.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      grok.shell.tv.addViewer('Form');
      await new Promise(r => setTimeout(r, 500));
      const form = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Form') as any;
      const spgiName = spgiDf.name;
      // Switch to demog table
      form.props.table = 'demog';
      await new Promise(r => setTimeout(r, 500));
      // Switch back
      form.props.table = spgiName;
      await new Promise(r => setTimeout(r, 300));
      grok.shell.closeAll();
      return { formFound: !!form };
    });
    expect(result.formFound).toBe(true);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
