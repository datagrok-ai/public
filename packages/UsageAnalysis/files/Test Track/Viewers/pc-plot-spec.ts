import {test, expect} from '@playwright/test';

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

test('PC Plot tests', async ({page}) => {
  test.setTimeout(600_000);

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
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Phase 3: Add PC Plot via toolbox
  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-pc-plot"]');
    if (icon) (icon as HTMLElement).click();
  });
  await page.locator('[name="viewer-PC-Plot"]').waitFor({timeout: 10000});

  // #### Menu Ribbon and To Script
  await softStep('Menu Ribbon and To Script', async () => {
    const viewerExists = await page.evaluate(() => !!grok.shell.tv.viewers.find(v => v.type === 'PC Plot'));
    expect(viewerExists).toBe(true);

    await page.evaluate(async () => {
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2
      }));
      await new Promise(r => setTimeout(r, 500));
      const items = Array.from(document.querySelectorAll('.d4-menu-item-label'));
      const toScript = items.find(el => el.textContent!.trim() === 'To Script');
      if (toScript) {
        const parent = toScript.closest('.d4-menu-item')!;
        parent.dispatchEvent(new MouseEvent('mousemove', { bubbles: true }));
        parent.dispatchEvent(new MouseEvent('mouseenter', { bubbles: true }));
        await new Promise(r => setTimeout(r, 300));
        const sub = Array.from(document.querySelectorAll('.d4-menu-item-label'));
        const js = sub.find(el => el.textContent!.trim() === 'To JavaScript');
        if (js) js.closest('.d4-menu-item')!.click();
      }
      await new Promise(r => setTimeout(r, 500));
    });

    await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot');
      if (pc) pc.close();
      await new Promise(r => setTimeout(r, 500));
      const icon = document.querySelector('[name="icon-pc-plot"]');
      if (icon) (icon as HTMLElement).click();
      await new Promise(r => setTimeout(r, 1000));
    });
    const reopened = await page.evaluate(() => !!grok.shell.tv.viewers.find(v => v.type === 'PC Plot'));
    expect(reopened).toBe(true);
  });

  // #### Axis scale & normalization
  await softStep('Axis scale & normalization', async () => {
    const result = await page.evaluate(() => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const defaultNorm = pc.props.normalizeEachColumn;
      pc.props.normalizeEachColumn = false;
      const afterDisable = pc.props.normalizeEachColumn;
      pc.props.normalizeEachColumn = true;
      const afterReenable = pc.props.normalizeEachColumn;
      pc.props.logColumnsColumnNames = ['AGE'];
      const logAge = pc.props.logColumnsColumnNames.slice();
      pc.props.logColumnsColumnNames = ['AGE', 'WEIGHT'];
      const logBoth = pc.props.logColumnsColumnNames.slice();
      pc.props.logColumnsColumnNames = [];
      return { defaultNorm, afterDisable, afterReenable, logAge, logBoth };
    });
    expect(result.defaultNorm).toBe(true);
    expect(result.afterDisable).toBe(false);
    expect(result.afterReenable).toBe(true);
    expect(result.logAge).toEqual(['AGE']);
    expect(result.logBoth).toEqual(['AGE', 'WEIGHT']);

    const menuResult = await page.evaluate(async () => {
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const rect = canvas.getBoundingClientRect();
      const openMenu = async () => {
        canvas.dispatchEvent(new MouseEvent('contextmenu', {
          bubbles: true, cancelable: true, button: 2,
          clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2
        }));
        await new Promise(r => setTimeout(r, 500));
      };
      const clickSub = async (parent: string, child: string) => {
        const items = Array.from(document.querySelectorAll('.d4-menu-item-label'));
        const p = items.find(el => el.textContent!.trim() === parent);
        if (!p) return false;
        const pm = p.closest('.d4-menu-item')!;
        pm.dispatchEvent(new MouseEvent('mousemove', { bubbles: true }));
        pm.dispatchEvent(new MouseEvent('mouseenter', { bubbles: true }));
        await new Promise(r => setTimeout(r, 300));
        const sub = Array.from(document.querySelectorAll('.d4-menu-item-label'));
        const c = sub.find(el => el.textContent!.trim() === child);
        if (c) c.closest('.d4-menu-item')!.click();
        await new Promise(r => setTimeout(r, 300));
        return !!c;
      };
      await openMenu();
      await clickSub('Y Axis', 'Global');
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const afterGlobal = pc.props.normalizeEachColumn;
      await openMenu();
      await clickSub('Y Axis', 'Normalized');
      const afterNorm = pc.props.normalizeEachColumn;
      return { afterGlobal, afterNorm };
    });
    expect(menuResult.afterGlobal).toBe(false);
    expect(menuResult.afterNorm).toBe(true);
  });

  // #### Selection & line display
  await softStep('Selection & line display', async () => {
    const result = await page.evaluate(() => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const r: any = {};
      r.defaultCurrent = pc.props.showCurrentLine;
      r.defaultMouseOver = pc.props.showMouseOverLine;
      r.defaultAll = pc.props.showAllLines;
      pc.props.showCurrentLine = false;
      r.currentOff = pc.props.showCurrentLine;
      pc.props.showCurrentLine = true;
      pc.props.showMouseOverLine = false;
      r.mouseOverOff = pc.props.showMouseOverLine;
      pc.props.showMouseOverLine = true;
      pc.props.showMouseOverRowGroup = true;
      pc.props.showAllLines = false;
      r.allOff = pc.props.showAllLines;
      pc.props.showAllLines = true;
      return r;
    });
    expect(result.defaultCurrent).toBe(true);
    expect(result.currentOff).toBe(false);
    expect(result.mouseOverOff).toBe(false);
    expect(result.allOff).toBe(false);
  });

  // #### Style & layout
  await softStep('Style & layout', async () => {
    const result = await page.evaluate(() => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      pc.props.lineWidth = 3;
      const lw = pc.props.lineWidth;
      pc.props.currentLineWidth = 5;
      const clw = pc.props.currentLineWidth;
      pc.props.mouseOverLineWidth = 5;
      pc.props.labelsOrientation = 'Vert';
      const lo = pc.props.labelsOrientation;
      pc.props.minMaxOrientation = 'Vert';
      pc.props.horzMargin = 60;
      const hm = pc.props.horzMargin;
      pc.props.autoLayout = false;
      const al = pc.props.autoLayout;
      pc.props.lineWidth = 0.5; pc.props.currentLineWidth = 2; pc.props.mouseOverLineWidth = 2;
      pc.props.labelsOrientation = 'Auto'; pc.props.minMaxOrientation = 'Auto';
      pc.props.horzMargin = 40; pc.props.autoLayout = true;
      return { lw, clw, lo, hm, al };
    });
    expect(result.lw).toBe(3);
    expect(result.clw).toBe(5);
    expect(result.lo).toBe('Vert');
    expect(result.al).toBe(false);
  });

  // #### In-chart filtering & reset
  await softStep('In-chart filtering & reset', async () => {
    const result = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      pc.props.showFilteredOutLines = true;
      const sfo = pc.props.showFilteredOutLines;
      const viewer = document.querySelector('[name="viewer-PC-Plot"]')!;
      const canvas = viewer.querySelector('canvas[name="canvas"]')!;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('dblclick', {
        bubbles: true, cancelable: true,
        clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2
      }));
      await new Promise(r => setTimeout(r, 500));
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2
      }));
      await new Promise(r => setTimeout(r, 500));
      const items = Array.from(document.querySelectorAll('.d4-menu-item-label'));
      const rv = items.find(el => el.textContent!.trim() === 'Reset View');
      if (rv) rv.closest('.d4-menu-item')!.click();
      await new Promise(r => setTimeout(r, 300));
      pc.props.showFilteredOutLines = false;
      return { sfo, resetFound: !!rv };
    });
    expect(result.sfo).toBe(true);
    expect(result.resetFound).toBe(true);
  });

  // #### Filter panel interaction
  await softStep('Filter panel interaction', async () => {
    const result = await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 500));
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 30, max: 50});
      await new Promise(r => setTimeout(r, 500));
      const afterFilter = grok.shell.tv.dataFrame.filter.trueCount;
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 18, max: 89});
      await new Promise(r => setTimeout(r, 500));
      const afterReset = grok.shell.tv.dataFrame.filter.trueCount;
      return { afterFilter, afterReset };
    });
    expect(result.afterFilter).toBeLessThan(5850);
    expect(result.afterReset).toBe(5850);
  });

  // #### Column management & reordering
  await softStep('Column management & reordering', async () => {
    const result = await page.evaluate(() => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const defaultCols = pc.props.columnNames.slice();
      pc.props.columnNames = defaultCols.filter((c: string) => c !== 'HEIGHT');
      const removed = pc.props.columnNames.slice();
      pc.props.columnNames = [...pc.props.columnNames, 'HEIGHT'];
      const added = pc.props.columnNames.slice();
      pc.props.columnNames = ['WEIGHT', 'AGE', 'HEIGHT', 'STARTED'];
      const reordered = pc.props.columnNames.slice();
      pc.props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'];
      return { removed, added, reordered };
    });
    expect(result.removed).not.toContain('HEIGHT');
    expect(result.added).toContain('HEIGHT');
    expect(result.reordered[0]).toBe('WEIGHT');
  });

  // #### Density styles
  await softStep('Density styles', async () => {
    const result = await page.evaluate(() => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      const def = pc.props.densityStyle;
      pc.props.densityStyle = 'box plot';
      const box = pc.props.densityStyle;
      pc.props.showInterquartileRange = false;
      const iqrOff = pc.props.showInterquartileRange;
      pc.props.showInterquartileRange = true;
      pc.props.showUpperDash = false; pc.props.showUpperDash = true;
      pc.props.showLowerDash = false; pc.props.showLowerDash = true;
      pc.props.showMeanCross = false; pc.props.showMeanCross = true;
      pc.props.showMedian = false; pc.props.showMedian = true;
      pc.props.showCircles = true;
      pc.props.densityStyle = 'violin plot';
      const violin = pc.props.densityStyle;
      pc.props.bins = 200;
      const bins = pc.props.bins;
      pc.props.whiskerLineWidth = 5;
      pc.props.densityStyle = 'circles';
      pc.props.bins = 100; pc.props.whiskerLineWidth = 2;
      return { def, box, iqrOff, violin, bins };
    });
    expect(result.def).toBe('circles');
    expect(result.box).toBe('box plot');
    expect(result.iqrOff).toBe(false);
    expect(result.violin).toBe('violin plot');
    expect(result.bins).toBe(200);
  });

  // #### Color coding, legend & grid coloring
  await softStep('Color coding, legend & grid coloring', async () => {
    const result = await page.evaluate(async () => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      pc.props.colorColumnName = 'AGE';
      const colorAge = pc.props.colorColumnName;
      pc.props.colorAxisType = 'logarithmic';
      const logColor = pc.props.colorAxisType;
      pc.props.invertColorScheme = true;
      const inverted = pc.props.invertColorScheme;
      pc.props.invertColorScheme = false;
      pc.props.colorMin = 30; pc.props.colorMax = 60;
      const cMin = pc.props.colorMin;
      pc.props.colorMin = null; pc.props.colorMax = null; pc.props.colorAxisType = 'linear';
      pc.props.colorColumnName = 'RACE';
      const colorRace = pc.props.colorColumnName;
      pc.props.legendPosition = 'Left';
      const lLeft = pc.props.legendPosition;
      pc.props.legendPosition = 'Right'; pc.props.legendPosition = 'Top'; pc.props.legendPosition = 'Bottom';
      pc.props.legendVisibility = 'Never';
      const lNever = pc.props.legendVisibility;
      pc.props.legendVisibility = 'Auto';
      pc.props.colorColumnName = 'HEIGHT';
      const df = grok.shell.tv.dataFrame;
      df.col('HEIGHT').meta.colors.setLinear([DG.Color.blue, DG.Color.red]);
      await new Promise(r => setTimeout(r, 500));
      df.col('HEIGHT').meta.colors.setConditional({'20-150': DG.Color.green, '150-250': DG.Color.orange});
      await new Promise(r => setTimeout(r, 500));
      df.col('HEIGHT').meta.colors.setLinear();
      pc.props.colorColumnName = '';
      return { colorAge, logColor, inverted, cMin, colorRace, lLeft, lNever };
    });
    expect(result.colorAge).toBe('AGE');
    expect(result.logColor).toBe('logarithmic');
    expect(result.colorRace).toBe('RACE');
    expect(result.lNever).toBe('Never');
  });

  // #### Title and description
  await softStep('Title and description', async () => {
    const result = await page.evaluate(() => {
      const pc = grok.shell.tv.viewers.find(v => v.type === 'PC Plot')!;
      pc.props.title = 'My PC Plot';
      const title = pc.props.title;
      pc.props.description = 'Test description';
      const desc = pc.props.description;
      pc.props.descriptionPosition = 'Bottom';
      const pos = pc.props.descriptionPosition;
      pc.props.title = ''; pc.props.description = '';
      return { title, desc, pos };
    });
    expect(result.title).toBe('My PC Plot');
    expect(result.desc).toBe('Test description');
  });

  // #### Pick Up / Apply
  await softStep('Pick Up / Apply', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const pc1 = tv.viewers.find(v => v.type === 'PC Plot')!;
      pc1.props.columnNames = ['AGE', 'WEIGHT', 'STARTED'];
      pc1.props.logColumnsColumnNames = ['AGE'];
      pc1.props.colorColumnName = 'RACE';
      pc1.props.legendPosition = 'Left';
      pc1.props.title = 'Source Plot';
      tv.addViewer('PC Plot');
      await new Promise(r => setTimeout(r, 500));
      const clickSub = async (idx: number, parent: string, child: string) => {
        const viewers = document.querySelectorAll('[name="viewer-PC-Plot"]');
        const canvas = viewers[idx].querySelector('canvas[name="canvas"]')!;
        const rect = canvas.getBoundingClientRect();
        canvas.dispatchEvent(new MouseEvent('contextmenu', {
          bubbles: true, cancelable: true, button: 2,
          clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2
        }));
        await new Promise(r => setTimeout(r, 500));
        const items = Array.from(document.querySelectorAll('.d4-menu-item-label'));
        const p = items.find(el => el.textContent!.trim() === parent);
        if (!p) return;
        const pm = p.closest('.d4-menu-item')!;
        pm.dispatchEvent(new MouseEvent('mousemove', { bubbles: true }));
        pm.dispatchEvent(new MouseEvent('mouseenter', { bubbles: true }));
        await new Promise(r => setTimeout(r, 300));
        const sub = Array.from(document.querySelectorAll('.d4-menu-item-label'));
        const c = sub.find(el => el.textContent!.trim() === child);
        if (c) c.closest('.d4-menu-item')!.click();
        await new Promise(r => setTimeout(r, 500));
      };
      await clickSub(0, 'Pick Up / Apply', 'Pick Up');
      await clickSub(1, 'Pick Up / Apply', 'Apply');
      const pcs = tv.viewers.filter(v => v.type === 'PC Plot');
      const r = { pc2Cols: pcs[1]?.props.columnNames?.slice(), pc2Color: pcs[1]?.props.colorColumnName, pc2Title: pcs[1]?.props.title };
      if (pcs[1]) pcs[1].close();
      const pcF = tv.viewers.find(v => v.type === 'PC Plot')!;
      pcF.props.columnNames = ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'];
      pcF.props.logColumnsColumnNames = []; pcF.props.colorColumnName = '';
      pcF.props.legendPosition = 'Auto'; pcF.props.title = '';
      return r;
    });
    expect(result.pc2Cols).toEqual(['AGE', 'WEIGHT', 'STARTED']);
    expect(result.pc2Color).toBe('RACE');
    expect(result.pc2Title).toBe('Source Plot');
  });

  // #### Layout and project save/restore
  await softStep('Layout and project save/restore', async () => {
    const result = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1000));
      tv.addViewer('Scatter plot');
      await new Promise(r => setTimeout(r, 500));
      const saved = await grok.dapi.layouts.find(layoutId);
      tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));
      const hasScatter = tv.viewers.some(v => v.type === 'Scatter plot');
      const hasPc = tv.viewers.some(v => v.type === 'PC Plot');
      await grok.dapi.layouts.delete(saved);
      return { hasScatter, hasPc };
    });
    expect(result.hasScatter).toBe(false);
    expect(result.hasPc).toBe(true);
  });

  // #### Table switching and transformation
  await softStep('Table switching and transformation', async () => {
    const result = await page.evaluate(async (path) => {
      const df2 = await grok.dapi.files.readCsv(path);
      const tv2 = grok.shell.addTableView(df2);
      await new Promise(resolve => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
        setTimeout(resolve, 3000);
      });
      const hasBioChem = Array.from({length: df2.columns.length}, (_, i) => df2.columns.byIndex(i))
        .some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
      if (hasBioChem) {
        for (let i = 0; i < 50; i++) {
          if (document.querySelector('[name="viewer-Grid"] canvas')) break;
          await new Promise(r => setTimeout(r, 200));
        }
        await new Promise(r => setTimeout(r, 5000));
      }
      const views = Array.from(grok.shell.views);
      const demogView = views.find((v: any) => v.name === 'Table');
      if (demogView) grok.shell.v = demogView;
      await new Promise(r => setTimeout(r, 500));
      const pc = grok.shell.tv.addViewer('PC Plot');
      await new Promise(r => setTimeout(r, 500));
      pc.props.table = df2.name;
      await new Promise(r => setTimeout(r, 500));
      const tableSet = pc.dataFrame?.name;
      pc.props.transformation = '[{"#type":"GroupAggregation","aggType":"key","colName":"Chemist 521"},{"#type":"GroupAggregation","aggType":"pivot","colName":"Series"},{"#type":"GroupAggregation","aggType":"count","colName":"Id"}]';
      await new Promise(r => setTimeout(r, 1000));
      pc.close();
      return { spgiRows: df2.rowCount, tableSet };
    }, spgiPath);
    expect(result.spgiRows).toBe(100);
    expect(result.tableSet).toBeTruthy();
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
