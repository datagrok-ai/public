import {test, expect, chromium, Page} from '@playwright/test';

declare const grok: any;
declare const DG: any;

const baseUrl = process.env.DATAGROK_URL ?? process.env.BASE_URL ?? 'https://dev.datagrok.ai';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Grid: comprehensive scenario (SPGI)', async () => {
  test.setTimeout(600_000);

  const cdp = await chromium.connectOverCDP('http://127.0.0.1:9222');
  const ctx = cdp.contexts()[0];
  const pages = ctx.pages();
  let page: Page = pages.find((p) => p.url().includes('datagrok')) ?? pages[0];
  if (!page) page = await ctx.newPage();
  await page.bringToFront();

  await page.goto(baseUrl, {timeout: 60000, waitUntil: 'networkidle'});
  await page.waitForFunction(() => {
    try {
      if (typeof grok === 'undefined' || !grok.shell || !grok.dapi || !grok.dapi.files) return false;
      grok.shell.closeAll();
      return true;
    } catch { return false; }
  }, {timeout: 180000, polling: 1000});
  await page.waitForTimeout(2000);

  // Setup: open SPGI with chem wait
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch {}
    try { grok.shell.windows.simpleMode = false; } catch {}
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
    df.name = 'SPGI';
    const tv = grok.shell.addTableView(df);
    tv.name = 'SPGI';
    await new Promise((r) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); r(undefined); });
      setTimeout(r, 4000);
    });
    for (let i = 0; i < 60; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise((r) => setTimeout(r, 200));
    }
    await new Promise((r) => setTimeout(r, 4000));
  });
  await page.locator('[name="viewer-Grid"]').first().waitFor({timeout: 30000});

  // #### 1. Viewer basics
  await softStep('1. Viewer basics — sort, resize, select, edit', async () => {
    const r = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const grid = tv.grid;
      const df = tv.dataFrame;
      grid.sort(['Average Mass'], [true]);
      grid.sort([]);
      grid.col('Series').width = 200;
      const newW = grid.col('Series').width;
      const origRH = grid.props.rowHeight;
      grid.props.rowHeight = 60;
      const newRH = grid.props.rowHeight;
      grid.props.rowHeight = origRH;
      df.selection.init((i: number) => i < 5);
      const selCnt = df.selection.trueCount;
      df.selection.setAll(false);
      const c = df.col('Series');
      const orig = c.get(0);
      c.set(0, 'TEST_VALUE');
      const after = c.get(0);
      c.set(0, orig);
      return {newW, newRH, selCnt, after};
    });
    expect(r.newW).toBe(200);
    expect(r.newRH).toBe(60);
    expect(r.selCnt).toBe(5);
    expect(r.after).toBe('TEST_VALUE');
  });

  // #### 2. Columns menu — filter panel, color coding
  await softStep('2. Columns — filter panel, color coding', async () => {
    const r = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      tv.getFiltersGroup({createDefaultFilters: true});
      await new Promise((r) => setTimeout(r, 600));
      const df = tv.dataFrame;
      const am = df.col('Average Mass');
      am.meta.colors.setLinear();
      const lin = am.meta.colors.getType();
      df.col('Series').meta.colors.setCategorical();
      const cat = df.col('Series').meta.colors.getType();
      am.meta.colors.setConditional({'<400': '#00ff00', '>=400': '#ff0000'});
      const cond = am.meta.colors.getType();
      am.meta.colors.setLinear();
      return {lin, cat, cond};
    });
    expect(r.lin).toBe('Linear');
    expect(r.cat).toBe('Categorical');
    expect(r.cond).toBe('Conditional');
    await page.locator('[name="viewer-Filters"]').waitFor({timeout: 10000});
  });

  // #### 3. Column context menu — renderer / sort / format / hide
  await softStep('3. Structure/Series/Chemical Space X', async () => {
    const r = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const grid = tv.grid;
      const s = df.col('Structure');
      s.tags['cell.renderer'] = 'Text';
      grid.invalidate();
      await new Promise((r) => setTimeout(r, 300));
      s.tags['cell.renderer'] = 'Molecule';
      grid.invalidate();
      grid.sort(['Series'], [true]);
      grid.sort([]);
      const gc = grid.col('Series');
      gc.visible = false;
      const v1 = gc.visible;
      gc.visible = true;
      df.col('Chemical Space X').tags['format'] = '0.000';
      const fmt = df.col('Chemical Space X').tags['format'];
      return {v1, fmt};
    });
    expect(r.v1).toBe(false);
    expect(r.fmt).toBe('0.000');
  });

  // #### 4. Context Panel
  await softStep('4. Context panel — current col, filter panel, permissions', async () => {
    const r = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      grok.shell.windows.showContextPanel = true;
      df.currentCol = df.col('Series');
      tv.getFiltersGroup({createDefaultFilters: true});
      await new Promise((r) => setTimeout(r, 500));
      df.col('Chemical Space X').tags['.editedBy'] = grok.shell.user.login;
      return {ctxShown: grok.shell.windows.showContextPanel, curName: df.currentCol?.name};
    });
    expect(r.ctxShown).toBe(true);
    expect(r.curName).toBe('Series');
  });

  // #### 5. Grid context menu — add column, summaries, link tables
  await softStep('5. Grid context menu — add / summaries / link tables', async () => {
    const r = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const grid = tv.grid;
      await df.columns.addNewCalculated('TestCol', '${Average Mass} * 2');
      const hasNew = !!df.col('TestCol');
      const cols = ['Average Mass', 'TPSA', 'NIBR logP'];
      const added: string[] = [];
      const summaries: [string, string][] = [['Summary', 'sparkline'], ['Bars', 'barchart'], ['PieSum', 'piechart'], ['SmartForm', 'Smart Form']];
      for (const [name, cellType] of summaries) {
        const c = grid.columns.add({gridColumnName: name, cellType});
        c.settings = {columnNames: cols};
        added.push(name);
      }
      const hv = tv.addViewer('Histogram', {valueColumnName: 'Average Mass'});
      hv.close();
      const origName = df.col('Chemical Space X').name;
      df.col('Chemical Space X').name = 'Chemical Space X Renamed';
      df.col('Chemical Space X Renamed').name = origName;
      const df2 = await grok.dapi.files.readCsv('System:DemoFiles/SPGI-linked1.csv');
      df2.name = 'SPGI-linked1';
      const tv2 = grok.shell.addTableView(df2);
      tv2.name = 'SPGI-linked1';
      await new Promise((r) => setTimeout(r, 1500));
      grok.data.linkTables(df, df2, ['Id'], ['Id'], [DG.SYNC_TYPE.CURRENT_ROW_TO_SELECTION]);
      let spgiView: any = null;
      for (const v of grok.shell.views) if (v.type === 'TableView' && v.dataFrame === df) spgiView = v;
      grok.shell.v = spgiView;
      await new Promise((r) => setTimeout(r, 500));
      grid.col('Lab Notebook').visible = false;
      const hidden = grid.col('Lab Notebook').visible;
      grid.col('Lab Notebook').visible = true;
      return {hasNew, added, hidden};
    });
    expect(r.hasNew).toBe(true);
    expect(r.added.length).toBe(4);
    expect(r.hidden).toBe(false);
  });

  // #### 6. Pick Up / Apply coloring
  await softStep('6. Pick Up / Apply coloring', async () => {
    await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const grid = tv.grid;
      const src = df.col('Average Mass');
      src.meta.colors.setLinear(['#ff0000', '#00ff00']);
      const dst = df.col('TPSA');
      for (const k of Object.keys(src.tags)) if (k.startsWith('.color-coding')) dst.tags[k] = src.tags[k];
      grid.invalidate();
      grid.props.colorCoding = 'All';
      src.meta.colors.setLinear(['#0000ff', '#ffff00']);
      const df2 = await grok.dapi.files.readCsv('System:DemoFiles/SPGI-linked2.csv');
      df2.name = 'SPGI-linked2';
      grok.shell.addTableView(df2);
      await new Promise((r) => setTimeout(r, 3000));
    });
  });

  // #### 7. Column groups — no public JS API; documented, not asserted
  await softStep('7. Column groups (SKIP — no public JS API)', async () => {
    const hasApi = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      return typeof grid.columns.addColumnGroup === 'function';
    });
    test.info().annotations.push({type: 'skip', description:
      `grid.columns.addColumnGroup exposed? ${hasApi} — column groups have no public JS API; canvas-only feature`});
  });

  // #### 8. Filtering
  await softStep('8. Filtering — categorical and numeric', async () => {
    const r = await page.evaluate(async () => {
      const df = Array.from(grok.shell.tables).find((t: any) => t.name === 'SPGI') as any;
      let tv: any = null;
      for (const v of grok.shell.views) if (v.type === 'TableView' && v.dataFrame === df) tv = v;
      grok.shell.v = tv;
      tv.getFiltersGroup({createDefaultFilters: true});
      await new Promise((r) => setTimeout(r, 1000));
      const first = df.col('Series').categories[0];
      df.rows.match({'Series': first}).filter();
      const catCount = df.filter.trueCount;
      df.filter.setAll(true);
      const am = df.col('Average Mass');
      df.filter.init((i: number) => am.get(i) > 400);
      const numCount = df.filter.trueCount;
      df.filter.setAll(true);
      return {catCount, numCount};
    });
    expect(r.catCount).toBeGreaterThan(0);
    expect(r.numCount).toBeGreaterThan(0);
    await page.locator('[name="viewer-Filters"]').waitFor({timeout: 10000});
  });

  // #### 9. Layout save / modify / restore
  await softStep('9. Layout save / modify / restore', async () => {
    const layoutId = await page.evaluate(async () => {
      const df = Array.from(grok.shell.tables).find((t: any) => t.name === 'SPGI') as any;
      let tv: any = null;
      for (const v of grok.shell.views) if (v.type === 'TableView' && v.dataFrame === df) tv = v;
      grok.shell.v = tv;
      const layout = tv.saveLayout();
      layout.name = 'grid-debug-layout-' + Date.now();
      await grok.dapi.layouts.save(layout);
      await new Promise((r) => setTimeout(r, 1200));
      tv.addViewer('Scatter plot', {xColumnName: 'Chemical Space X', yColumnName: 'Chemical Space Y'});
      await new Promise((r) => setTimeout(r, 1000));
      return layout.id;
    });
    expect(layoutId).toBeTruthy();
    await page.evaluate(async (id: string) => {
      const saved = await grok.dapi.layouts.find(id);
      const df = Array.from(grok.shell.tables).find((t: any) => t.name === 'SPGI') as any;
      let tv: any = null;
      for (const v of grok.shell.views) if (v.type === 'TableView' && v.dataFrame === df) tv = v;
      tv.loadLayout(saved);
      await new Promise((r) => setTimeout(r, 3000));
      await grok.dapi.layouts.delete(saved);
    }, layoutId);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
