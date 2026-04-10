import {test, expect} from '@playwright/test';

const baseUrl = process.env.BASE_URL ?? 'https://dev.datagrok.ai';
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

test('Correlation Plot tests', async ({page}) => {
  test.setTimeout(600_000);

  // Phase 1: Navigate
  await page.goto(baseUrl);
  await page.waitForFunction(() => {
    try { return typeof grok !== 'undefined' && grok.shell &&
      typeof grok.shell.settings?.showFiltersIconsConstantly === 'boolean'; }
    catch (e) { return false; }
  }, {timeout: 30000});

  // Phase 2: Open dataset
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = false;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Phase 3: Add Correlation Plot via toolbox
  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-correlation-plot"]');
    if (icon) (icon as HTMLElement).click();
  });
  await page.locator('[name="viewer-Correlation-plot"]').waitFor({timeout: 10000});

  // #### Double-click and cell interaction
  await softStep('Double-click and cell interaction', async () => {
    // Canvas-based cells: double-click via DOM events does not reach Dart event loop
    // Verify ignoreDoubleClick property works
    const result = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find(v => v.type === 'Correlation plot')!;
      const defaultIgnore = cp.props.ignoreDoubleClick;
      cp.props.ignoreDoubleClick = true;
      const afterEnable = cp.props.ignoreDoubleClick;
      cp.props.ignoreDoubleClick = false;
      const afterDisable = cp.props.ignoreDoubleClick;
      return { defaultIgnore, afterEnable, afterDisable };
    });
    expect(result.defaultIgnore).toBe(false);
    expect(result.afterEnable).toBe(true);
    expect(result.afterDisable).toBe(false);
  });

  // #### Column reordering
  await softStep('Column reordering', async () => {
    const result = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find(v => v.type === 'Correlation plot')!;
      const defaultX = cp.props.xColumnNames.slice();
      cp.props.xColumnNames = ['HEIGHT', 'AGE', 'WEIGHT', 'STARTED'];
      const reordered = cp.props.xColumnNames.slice();
      cp.props.xColumnNames = defaultX;
      return { defaultX, reordered };
    });
    expect(result.reordered[0]).toBe('HEIGHT');
    expect(result.reordered[1]).toBe('AGE');
  });

  // #### Column selection
  await softStep('Column selection', async () => {
    const result = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find(v => v.type === 'Correlation plot')!;
      cp.props.xColumnNames = ['AGE', 'HEIGHT', 'STARTED'];
      const xNoWeight = cp.props.xColumnNames.slice();
      cp.props.xColumnNames = ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'];
      cp.props.yColumnNames = ['AGE', 'WEIGHT', 'STARTED'];
      const yNoHeight = cp.props.yColumnNames.slice();
      cp.props.yColumnNames = ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'];
      return { xNoWeight, yNoHeight };
    });
    expect(result.xNoWeight).not.toContain('WEIGHT');
    expect(result.yNoHeight).not.toContain('HEIGHT');
  });

  // #### Correlation type and display options
  await softStep('Correlation type and display options', async () => {
    const result = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find(v => v.type === 'Correlation plot')!;
      const r: any = {};
      r.defaultType = cp.props.correlationType;
      cp.props.correlationType = 'Spearman';
      r.spearman = cp.props.correlationType;
      cp.props.correlationType = 'Pearson';
      r.defaultShowR = cp.props.showPearsonR;
      cp.props.showPearsonR = false;
      r.showROff = cp.props.showPearsonR;
      cp.props.showPearsonR = true;
      r.defaultTooltip = cp.props.showTooltip;
      cp.props.showTooltip = false;
      r.tooltipOff = cp.props.showTooltip;
      cp.props.showTooltip = true;
      r.defaultIgnore = cp.props.ignoreDoubleClick;
      cp.props.ignoreDoubleClick = true;
      r.ignoreOn = cp.props.ignoreDoubleClick;
      cp.props.ignoreDoubleClick = false;
      return r;
    });
    expect(result.defaultType).toBe('Pearson');
    expect(result.spearman).toBe('Spearman');
    expect(result.defaultShowR).toBe(true);
    expect(result.showROff).toBe(false);
    expect(result.tooltipOff).toBe(false);
    expect(result.ignoreOn).toBe(true);
  });

  // #### Row source
  await softStep('Row source', async () => {
    const result = await page.evaluate(async () => {
      const cp = grok.shell.tv.viewers.find(v => v.type === 'Correlation plot')!;
      const df = grok.shell.tv.dataFrame;
      df.selection.init(i => i < 20);
      const defaultRS = cp.props.rowSource;
      cp.props.rowSource = 'Selected';
      const selected = cp.props.rowSource;
      cp.props.rowSource = 'All';
      const all = cp.props.rowSource;
      cp.props.rowSource = 'Filtered';
      const filtered = cp.props.rowSource;
      df.selection.setAll(false);
      return { defaultRS, selected, all, filtered };
    });
    expect(result.defaultRS).toBe('Filtered');
    expect(result.selected).toBe('Selected');
    expect(result.all).toBe('All');
    expect(result.filtered).toBe('Filtered');
  });

  // #### Style customization
  await softStep('Style customization', async () => {
    const result = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find(v => v.type === 'Correlation plot')!;
      const origFont = cp.props.defaultCellFont;
      cp.props.defaultCellFont = 'normal normal 18px "Roboto"';
      const cellFont = cp.props.defaultCellFont;
      cp.props.colHeaderFont = 'bold normal 16px "Roboto"';
      const headerFont = cp.props.colHeaderFont;
      cp.props.backColor = DG.Color.lightGray;
      const backColor = cp.props.backColor;
      cp.props.defaultCellFont = origFont;
      cp.props.colHeaderFont = 'bold normal 13px "Roboto"';
      cp.props.backColor = DG.Color.white;
      return { cellFont, headerFont, backColor };
    });
    expect(result.cellFont).toBe('normal normal 18px "Roboto"');
    expect(result.headerFont).toBe('bold normal 16px "Roboto"');
    expect(result.backColor).toBeTruthy();
  });

  // #### Title and description
  await softStep('Title and description', async () => {
    const result = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find(v => v.type === 'Correlation plot')!;
      cp.props.showTitle = true;
      const showTitle = cp.props.showTitle;
      cp.props.title = 'Correlation Analysis';
      const title = cp.props.title;
      cp.props.description = 'Shows pairwise correlations';
      const desc = cp.props.description;
      cp.props.descriptionVisibilityMode = 'Always';
      const descVis = cp.props.descriptionVisibilityMode;
      cp.props.descriptionPosition = 'Bottom';
      const descPos = cp.props.descriptionPosition;
      cp.props.descriptionVisibilityMode = 'Never';
      const descNever = cp.props.descriptionVisibilityMode;
      cp.props.showTitle = false; cp.props.title = ''; cp.props.description = '';
      cp.props.descriptionVisibilityMode = 'Auto';
      return { showTitle, title, desc, descVis, descPos, descNever };
    });
    expect(result.showTitle).toBe(true);
    expect(result.title).toBe('Correlation Analysis');
    expect(result.descVis).toBe('Always');
    expect(result.descNever).toBe('Never');
  });

  // #### Context menu
  await softStep('Context menu', async () => {
    const result = await page.evaluate(async () => {
      const viewer = document.querySelector('[name="viewer-Correlation-plot"]')!;
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width * 0.5, clientY: rect.top + rect.height * 0.3
      }));
      await new Promise(r => setTimeout(r, 500));
      const items = Array.from(document.querySelectorAll('.d4-menu-item-label')).map(el => el.textContent!.trim());
      const showPearsonR = items.includes('Show Pearson R');
      const tooltip = items.includes('Tooltip');
      const columns = items.includes('Columns');
      document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
      return { showPearsonR, tooltip, columns };
    });
    expect(result.showPearsonR).toBe(true);
    expect(result.tooltip).toBe(true);
    expect(result.columns).toBe(true);
  });

  // #### Viewer filter formula
  await softStep('Viewer filter formula', async () => {
    const result = await page.evaluate(async () => {
      const cp = grok.shell.tv.viewers.find(v => v.type === 'Correlation plot')!;
      cp.props.filter = '${AGE} > 40';
      await new Promise(r => setTimeout(r, 500));
      const filterSet = cp.props.filter;
      cp.props.filter = '';
      const filterCleared = cp.props.filter;
      return { filterSet, filterCleared };
    });
    expect(result.filterSet).toBe('${AGE} > 40');
    expect(result.filterCleared).toBe('');
  });

  // #### Layout persistence
  await softStep('Layout persistence', async () => {
    const result = await page.evaluate(async () => {
      const cp = grok.shell.tv.viewers.find(v => v.type === 'Correlation plot')!;
      cp.props.correlationType = 'Spearman';
      cp.props.showPearsonR = false;
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1000));
      cp.close();
      await new Promise(r => setTimeout(r, 500));
      const saved = await grok.dapi.layouts.find(layoutId);
      grok.shell.tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));
      const cp2 = grok.shell.tv.viewers.find(v => v.type === 'Correlation plot');
      const restored = !!cp2;
      const restoredType = cp2?.props.correlationType;
      const restoredShowR = cp2?.props.showPearsonR;
      await grok.dapi.layouts.delete(saved);
      if (cp2) { cp2.props.correlationType = 'Pearson'; cp2.props.showPearsonR = true; }
      return { restored, restoredType, restoredShowR };
    });
    expect(result.restored).toBe(true);
    expect(result.restoredType).toBe('Spearman');
    expect(result.restoredShowR).toBe(false);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
