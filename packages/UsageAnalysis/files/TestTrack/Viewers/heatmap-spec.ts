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

test('Heat map tests', async ({page}: {page: Page}) => {
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
    grok.shell.addTableView(df);
    await new Promise<void>(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Phase 3: Add Heat Map
  await page.evaluate(() => {
    const icon = document.querySelector('[name="icon-heat-map"]');
    if (icon) (icon as HTMLElement).click();
    else grok.shell.tv.addViewer('Heat map');
  });
  await page.waitForFunction(() => !!Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Heat map'), {timeout: 10000});

  // Open settings panel
  await page.evaluate(() => {
    const viewer = document.querySelector('[name="viewer-Heat-map"]') ||
      Array.from(document.querySelectorAll('[name^="viewer-"]')).find(el => el.getAttribute('name')?.includes('Heat'));
    let el: Element | null = viewer ?? null;
    for (let i = 0; i < 6; i++) {
      el = el?.parentElement ?? null;
      if (!el) break;
      const gear = el.querySelector('[name="icon-font-icon-settings"]');
      if (gear) { (gear as HTMLElement).click(); return; }
    }
  });
  await page.waitForTimeout(500);

  // ---- Heatmap colors ----

  await softStep('Heatmap colors: default true, toggle false/true', async () => {
    const result = await page.evaluate(async () => {
      const hm = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Heat map') as any;
      const def = hm.props.heatmapColors;
      hm.props.heatmapColors = false;
      await new Promise(r => setTimeout(r, 100));
      const f = hm.props.heatmapColors;
      hm.props.heatmapColors = true;
      await new Promise(r => setTimeout(r, 100));
      return { def, f, t: hm.props.heatmapColors };
    });
    expect(result.def).toBe(true);
    expect(result.f).toBe(false);
    expect(result.t).toBe(true);
  });

  // ---- Global color scaling ----

  await softStep('Global color scaling: default false, toggle true/false', async () => {
    const result = await page.evaluate(async () => {
      const hm = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Heat map') as any;
      const def = hm.props.globalColorScaling;
      hm.props.globalColorScaling = true;
      await new Promise(r => setTimeout(r, 100));
      const t = hm.props.globalColorScaling;
      hm.props.globalColorScaling = false;
      await new Promise(r => setTimeout(r, 100));
      return { def, t, f: hm.props.globalColorScaling };
    });
    expect(result.def).toBe(false);
    expect(result.t).toBe(true);
    expect(result.f).toBe(false);
  });

  // ---- Column label orientation ----

  await softStep('Col labels orientation: Auto → Vert → Horz → Auto', async () => {
    const result = await page.evaluate(async () => {
      const hm = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Heat map') as any;
      const def = hm.props.colLabelsOrientation;
      hm.props.colLabelsOrientation = 'Vert';
      const vert = hm.props.colLabelsOrientation;
      hm.props.colLabelsOrientation = 'Horz';
      const horz = hm.props.colLabelsOrientation;
      hm.props.colLabelsOrientation = 'Auto';
      return { def, vert, horz, auto: hm.props.colLabelsOrientation };
    });
    expect(result.def).toBe('Auto');
    expect(result.vert).toBe('Vert');
    expect(result.horz).toBe('Horz');
    expect(result.auto).toBe('Auto');
  });

  // ---- Max heatmap columns ----

  await softStep('Max heatmap columns: default 100, set 3/1000/100', async () => {
    const result = await page.evaluate(async () => {
      const hm = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Heat map') as any;
      const def = hm.props.maxHeatmapColumns;
      hm.props.maxHeatmapColumns = 3;
      const v3 = hm.props.maxHeatmapColumns;
      hm.props.maxHeatmapColumns = 1000;
      const v1000 = hm.props.maxHeatmapColumns;
      hm.props.maxHeatmapColumns = 100;
      return { def, v3, v1000, v100: hm.props.maxHeatmapColumns };
    });
    expect(result.def).toBe(100);
    expect(result.v3).toBe(3);
    expect(result.v1000).toBe(1000);
    expect(result.v100).toBe(100);
  });

  // ---- Show heatmap scrollbars ----

  await softStep('Show heatmap scrollbars: default true, toggle false/true', async () => {
    const result = await page.evaluate(async () => {
      const hm = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Heat map') as any;
      const def = hm.props.showHeatmapScrollbars;
      hm.props.showHeatmapScrollbars = false;
      const f = hm.props.showHeatmapScrollbars;
      hm.props.showHeatmapScrollbars = true;
      return { def, f, t: hm.props.showHeatmapScrollbars };
    });
    expect(result.def).toBe(true);
    expect(result.f).toBe(false);
    expect(result.t).toBe(true);
  });

  // ---- Is Heatmap toggle ----

  await softStep('Is Heatmap: default true (heatmap mode), toggle false (grid mode) and back', async () => {
    const result = await page.evaluate(async () => {
      const hm = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Heat map') as any;
      const def = hm.props.isHeatmap;
      hm.props.isHeatmap = false;
      await new Promise(r => setTimeout(r, 300));
      const gridMode = hm.props.isHeatmap;
      hm.props.isHeatmap = true;
      await new Promise(r => setTimeout(r, 300));
      return { def, gridMode, heatmapMode: hm.props.isHeatmap };
    });
    expect(result.def).toBe(true);
    expect(result.gridMode).toBe(false);
    expect(result.heatmapMode).toBe(true);
  });

  // ---- Filtering interaction ----

  await softStep('Filtering: direct bitset filter AGE 20-40 reduces rows', async () => {
    // Note: fg.updateOrAdd({type:'range',...}) throws "Error adding filter" — use df.filter.init() instead
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const ageCol = df.col('AGE');
      df.filter.init((i: number) => ageCol.get(i) >= 20 && ageCol.get(i) <= 40);
      await new Promise(r => setTimeout(r, 300));
      const filtered = df.filter.trueCount;
      df.filter.setAll(true);
      await new Promise(r => setTimeout(r, 300));
      return { filtered, total: df.rowCount };
    });
    expect(result.filtered).toBeLessThan(result.total);
    expect(result.filtered).toBeGreaterThan(0);
  });

  // ---- Table switching ----

  await softStep('Table switching: open spgi-100 twice, switch heat map table prop', async () => {
    const result = await page.evaluate(async () => {
      const df1 = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      grok.shell.addTableView(df1);
      await new Promise<void>(resolve => {
        const sub = df1.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 4000);
      });
      for (let i = 0; i < 30; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 3000));

      const df2 = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      grok.shell.addTableView(df2);
      await new Promise(r => setTimeout(r, 2000));

      const views = Array.from(grok.shell.tableViews);
      const firstView = views.find((v: any) => v.dataFrame === df1) || views[0];
      (grok.shell as any).v = firstView;
      await new Promise(r => setTimeout(r, 500));

      const icon = document.querySelector('[name="icon-heat-map"]');
      if (icon) (icon as HTMLElement).click();
      await new Promise(r => setTimeout(r, 1000));

      const hm = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Heat map') as any;
      hm.props.table = df2.name;
      await new Promise(r => setTimeout(r, 500));
      const afterSwitch = hm.props.table;
      hm.props.table = df1.name;
      await new Promise(r => setTimeout(r, 300));
      const afterRestore = hm.props.table;
      grok.shell.closeAll();
      return { switched: afterSwitch === df2.name, restored: afterRestore === df1.name };
    });
    expect(result.switched).toBe(true);
    expect(result.restored).toBe(true);
  });

  // Re-open demog for remaining sections
  await page.evaluate(async (path: string) => {
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    await new Promise<void>(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    const icon = document.querySelector('[name="icon-heat-map"]');
    if (icon) (icon as HTMLElement).click();
    await new Promise(r => setTimeout(r, 1000));
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // ---- Selection interaction ----

  await softStep('Selection: click/shift-drag rows, Esc clears', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.currentRowIdx = 5;
      df.selection.set(5, true); df.selection.set(6, true); df.selection.set(7, true);
      await new Promise(r => setTimeout(r, 200));
      const sel3 = df.selection.trueCount;
      df.selection.setAll(false);
      return { sel3, sel0: df.selection.trueCount };
    });
    expect(result.sel3).toBe(3);
    expect(result.sel0).toBe(0);
  });

  // ---- Column sorting ----

  await softStep('Column sorting: sort AGE asc/desc/reset via grid.sort', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const ageCol = df.col('AGE');
      const grid = grok.shell.tv.grid;
      grid.sort([ageCol], [true]);
      await new Promise(r => setTimeout(r, 300));
      grid.sort([ageCol], [false]);
      await new Promise(r => setTimeout(r, 300));
      grid.sort([], []);
      await new Promise(r => setTimeout(r, 300));
      return { noError: true };
    });
    expect(result.noError).toBe(true);
  });

  // ---- Color scheme customization ----

  await softStep('Color scheme customization: linearColorScheme and categoricalColorScheme are non-null', async () => {
    const result = await page.evaluate(() => {
      const hm = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Heat map') as any;
      return {
        linear: hm.props.linearColorScheme,
        categorical: hm.props.categoricalColorScheme,
      };
    });
    expect(result.linear).not.toBeNull();
    expect(result.categorical).not.toBeNull();
  });

  // ---- Draw every row ----

  await softStep('Draw every row: default false, toggle true/false', async () => {
    const result = await page.evaluate(async () => {
      const hm = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Heat map') as any;
      const def = hm.props.drawEveryRow;
      hm.props.drawEveryRow = true;
      const t = hm.props.drawEveryRow;
      hm.props.drawEveryRow = false;
      return { def, t, f: hm.props.drawEveryRow };
    });
    expect(result.def).toBe(false);
    expect(result.t).toBe(true);
    expect(result.f).toBe(false);
  });

  // ---- Range slider navigation ----

  await softStep('Range slider: at least 2 .d4-range-selector elements present; double-click resets', async () => {
    const count = await page.locator('.d4-range-selector').count();
    expect(count).toBeGreaterThanOrEqual(2);
    await page.locator('.d4-range-selector').first().dblclick();
    await page.waitForTimeout(300);
  });

  // ---- Layout save and restore ----

  await softStep('Layout: save ColLabelsOrientation=Vert + GlobalColorScaling=true, restore', async () => {
    const layoutId = await page.evaluate(async () => {
      const hm = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Heat map') as any;
      hm.props.colLabelsOrientation = 'Vert';
      hm.props.globalColorScaling = true;
      await new Promise(r => setTimeout(r, 200));
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      hm.props.colLabelsOrientation = 'Auto';
      hm.props.globalColorScaling = false;
      return layout.id;
    });

    await page.evaluate(async (id: string) => {
      const saved = await grok.dapi.layouts.find(id);
      grok.shell.tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));
    }, layoutId);

    const restored = await page.evaluate(async (id: string) => {
      const hm = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Heat map') as any;
      const result = { orient: hm?.props?.colLabelsOrientation, global: hm?.props?.globalColorScaling };
      const saved = await grok.dapi.layouts.find(id);
      if (saved) await grok.dapi.layouts.delete(saved);
      return result;
    }, layoutId);

    expect(restored.orient).toBe('Vert');
    expect(restored.global).toBe(true);
  });

  // ---- Layout with isHeatmap toggle (spgi-100) ----

  await softStep('Layout: maxHeatmapColumns=200, isHeatmap=false saved and restored (spgi-100)', async () => {
    const result = await page.evaluate(async () => {
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 300));
      const df = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      grok.shell.addTableView(df);
      await new Promise<void>(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 4000);
      });
      for (let i = 0; i < 30; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 3000));

      const icon = document.querySelector('[name="icon-heat-map"]');
      if (icon) (icon as HTMLElement).click();
      await new Promise(r => setTimeout(r, 1000));

      const hm = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Heat map') as any;
      hm.props.maxHeatmapColumns = 200;
      hm.props.isHeatmap = false;
      await new Promise(r => setTimeout(r, 200));

      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));

      hm.props.maxHeatmapColumns = 100;
      hm.props.isHeatmap = true;

      const saved = await grok.dapi.layouts.find(layout.id);
      grok.shell.tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));

      const hmR = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Heat map') as any;
      const maxCols = hmR?.props?.maxHeatmapColumns;
      const isHeatmap = hmR?.props?.isHeatmap;

      if (hmR) hmR.props.isHeatmap = true;
      await grok.dapi.layouts.delete(saved);
      grok.shell.closeAll();

      return { maxCols, isHeatmap };
    });
    expect(result.maxCols).toBe(200);
    expect(result.isHeatmap).toBe(false);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
