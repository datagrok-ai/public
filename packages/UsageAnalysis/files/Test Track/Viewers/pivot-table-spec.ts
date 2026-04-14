import {test, expect, chromium} from '@playwright/test';

declare const grok: any;

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:DemoFiles/SPGI.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Pivot table viewer', async () => {
  test.setTimeout(600_000);

  const browser = await chromium.connectOverCDP('http://127.0.0.1:9222');
  const context = browser.contexts()[0];
  let page = context.pages().find(p => p.url().includes('datagrok'));
  if (!page) {
    page = await context.newPage();
    await page.goto(baseUrl, {waitUntil: 'networkidle', timeout: 60000});
  }
  await page.waitForFunction(() => {
    try {
      return typeof grok !== 'undefined' && grok.shell
        && typeof grok.shell.closeAll === 'function'
        && grok.dapi && grok.dapi.files;
    } catch { return false; }
  }, {timeout: 60000});

  // Phase 2: Open demog
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch {}
    try { grok.shell.windows.simpleMode = false; } catch {}
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
  await page.evaluate(() => grok.shell.tv.getFiltersGroup());

  // #### Step 2: Add Pivot table via toolbox icon
  await softStep('Add Pivot table — defaults DIS_POP/SEVERITY/AGE(avg)', async () => {
    const result = await page.evaluate(async () => {
      (document.querySelector('[name="icon-pivot-table"]') as HTMLElement)?.click();
      await new Promise(r => setTimeout(r, 1500));
      const pv = grok.shell.tv.viewers.find((v: any) => v.type === 'Pivot table');
      return {
        present: !!document.querySelector('[name="viewer-Pivot-table"]'),
        gb: pv?.props.groupByColumnNames,
        p: pv?.props.pivotColumnNames,
        a: pv?.props.aggregateColumnNames,
        at: pv?.props.aggregateAggTypes,
      };
    });
    expect(result.present).toBe(true);
    expect(result.gb).toEqual(['DIS_POP']);
    expect(result.p).toEqual(['SEVERITY']);
    expect(result.a).toEqual(['AGE']);
    expect(result.at).toEqual(['avg']);
  });

  // #### Step 3: Close via panel titlebar and re-add — defaults persist
  await softStep('Close via titlebar then re-add', async () => {
    const result = await page.evaluate(async () => {
      const tbs = Array.from(document.querySelectorAll('.panel-titlebar'));
      const pivotTb = tbs.find(tb => tb.textContent?.trim() === 'Pivot table');
      const closeBtn = pivotTb?.querySelector('[name="Close"], .panel-titlebar-button-close') as HTMLElement;
      closeBtn?.click();
      await new Promise(r => setTimeout(r, 800));
      const closed = !document.querySelector('[name="viewer-Pivot-table"]');
      (document.querySelector('[name="icon-pivot-table"]') as HTMLElement)?.click();
      await new Promise(r => setTimeout(r, 1500));
      const pv = grok.shell.tv.viewers.find((v: any) => v.type === 'Pivot table');
      return {
        closed,
        reopened: !!document.querySelector('[name="viewer-Pivot-table"]'),
        gb: pv?.props.groupByColumnNames,
        p: pv?.props.pivotColumnNames,
        a: pv?.props.aggregateColumnNames,
      };
    });
    expect(result.closed).toBe(true);
    expect(result.reopened).toBe(true);
    expect(result.gb).toEqual(['DIS_POP']);
    expect(result.p).toEqual(['SEVERITY']);
    expect(result.a).toEqual(['AGE']);
  });

  // #### Step 4: Modify tag editors via property API (UI combo popup is brittle headlessly)
  await softStep('Modify tag editors (Group by/Pivot/Aggregate)', async () => {
    const result = await page.evaluate(async () => {
      const pv = grok.shell.tv.viewers.find((v: any) => v.type === 'Pivot table');
      pv.props.groupByColumnNames = ['SEX'];
      await new Promise(r => setTimeout(r, 300));
      pv.props.pivotColumnNames = ['RACE'];
      await new Promise(r => setTimeout(r, 300));
      pv.props.aggregateColumnNames = ['HEIGHT', 'WEIGHT'];
      pv.props.aggregateAggTypes = ['avg', 'sum'];
      await new Promise(r => setTimeout(r, 400));
      pv.props.aggregateAggTypes = ['avg', 'min'];
      await new Promise(r => setTimeout(r, 400));
      return {
        gb: pv.props.groupByColumnNames,
        p: pv.props.pivotColumnNames,
        a: pv.props.aggregateColumnNames,
        at: pv.props.aggregateAggTypes,
      };
    });
    expect(result.gb).toEqual(['SEX']);
    expect(result.p).toEqual(['RACE']);
    expect(result.a).toEqual(['HEIGHT', 'WEIGHT']);
    expect(result.at).toEqual(['avg', 'min']);
  });

  // #### Step 5-6: Property pane toggles
  await softStep('Toggle Show Header / Command Bar / Filtering / Row Source', async () => {
    const result = await page.evaluate(async () => {
      const pv = grok.shell.tv.viewers.find((v: any) => v.type === 'Pivot table');
      (document.querySelector('[name="viewer-Pivot-table"] [name="icon-font-icon-settings"]') as HTMLElement)?.click();
      await new Promise(r => setTimeout(r, 600));
      pv.props.showHeader = false;
      await new Promise(r => setTimeout(r, 200));
      const top = document.querySelector('[name="viewer-Pivot-table"] .grok-pivot-top') as HTMLElement;
      const headerHidden = !top || top.offsetHeight === 0;
      pv.props.showHeader = true;
      pv.props.showCommandBar = false;
      await new Promise(r => setTimeout(r, 200));
      pv.props.showCommandBar = true;
      for (const rs of ['Filtered', 'All', 'Selected', 'Filtered']) {
        pv.props.rowSource = rs;
        await new Promise(r => setTimeout(r, 200));
      }
      return { headerHidden, finalRowSource: pv.props.rowSource };
    });
    expect(result.headerHidden).toBe(true);
    expect(result.finalRowSource).toBe('Filtered');
  });

  // #### Step 7: Title and layout persistence (SPGI)
  await softStep('SPGI: title and coloring persist after layout reapply', async () => {
    const result = await page.evaluate(async (path) => {
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });
      await new Promise(r => setTimeout(r, 1000));
      (document.querySelector('[name="icon-pivot-table"]') as HTMLElement)?.click();
      await new Promise(r => setTimeout(r, 1500));
      const pv = tv.viewers.find((v: any) => v.type === 'Pivot table');
      pv.props.title = 'My Pivot SPGI';
      const aggColName = pv.props.aggregateColumnNames[0];
      const col = tv.dataFrame.col(aggColName);
      col.meta.colors.setLinear([0xFF0000FF, 0xFFFF0000]);
      await new Promise(r => setTimeout(r, 400));
      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1500));
      pv.props.title = 'CHANGED';
      const saved = await grok.dapi.layouts.find(layout.id);
      tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));
      const pv2 = tv.viewers.find((v: any) => v.type === 'Pivot table');
      const titleAfter = pv2?.props.title;
      const colorAfter = col.meta.colors.getType();
      try { await grok.dapi.layouts.delete(saved); } catch {}
      return { titleAfter, colorAfter };
    }, spgiPath);
    expect(result.titleAfter).toBe('My Pivot SPGI');
    expect(result.colorAfter).toBe('Linear');
  });

  // #### Step 8: Coloring preserved across rowSource changes (demog)
  await softStep('Coloring preserved across rowSource changes', async () => {
    const result = await page.evaluate(async () => {
      const tvs = Array.from(grok.shell.tableViews) as any[];
      const demogTv = tvs.find(t => t.dataFrame.columns.contains('DIS_POP'));
      grok.shell.v = demogTv;
      await new Promise(r => setTimeout(r, 500));
      const pv = demogTv.viewers.find((v: any) => v.type === 'Pivot table');
      pv.props.aggregateColumnNames = ['HEIGHT'];
      pv.props.aggregateAggTypes = ['avg'];
      pv.props.rowSource = 'Filtered';
      await new Promise(r => setTimeout(r, 400));
      const col = demogTv.dataFrame.col('HEIGHT');
      col.meta.colors.setLinear([0xFF00FF00, 0xFFFF00FF]);
      await new Promise(r => setTimeout(r, 300));
      const before = col.meta.colors.getType();
      pv.props.rowSource = 'Selected';
      await new Promise(r => setTimeout(r, 400));
      pv.props.rowSource = 'Filtered';
      await new Promise(r => setTimeout(r, 400));
      const after = col.meta.colors.getType();
      return { before, after };
    });
    expect(result.before).toBe('Linear');
    expect(result.after).toBe('Linear');
  });

  // #### Step 9: Two-way property sync (props -> tag DOM)
  await softStep('Property changes reflected in tag editors', async () => {
    const result = await page.evaluate(async () => {
      const pv = grok.shell.tv.viewers.find((v: any) => v.type === 'Pivot table');
      pv.props.groupByColumnNames = ['DIS_POP'];
      pv.props.aggregateColumnNames = ['AGE'];
      pv.props.aggregateAggTypes = ['avg'];
      await new Promise(r => setTimeout(r, 400));
      const v = document.querySelector('[name="viewer-Pivot-table"]')!;
      const panels = Array.from(v.querySelectorAll('.grok-pivot-column-panel'));
      const findRow = (rowName: string) => {
        for (const p of panels) {
          const t = p.querySelector('.grok-pivot-column-tags-title');
          if (t && t.textContent?.trim() === rowName)
            return Array.from(p.querySelectorAll('.d4-tag')).map(x => x.textContent?.trim());
        }
        return null;
      };
      return { gbTags: findRow('Group by'), aggTags: findRow('Aggregate') };
    });
    expect(result.gbTags).toContain('DIS_POP');
    expect(result.aggTags?.some(t => t?.includes('AGE'))).toBe(true);
  });

  // #### Step 10: ADD button pushes pivot result to workspace
  await softStep('ADD button creates new aggregated table', async () => {
    const result = await page.evaluate(async () => {
      const before = grok.shell.tables.length;
      const addBtn = document.querySelector('[name="viewer-Pivot-table"] .grok-pivot-counts [name="button-ADD"]') as HTMLElement;
      addBtn?.click();
      await new Promise(r => setTimeout(r, 1500));
      const after = grok.shell.tables.length;
      const newDf = grok.shell.tables[grok.shell.tables.length - 1];
      return { added: after - before, name: newDf?.name, rows: newDf?.rowCount };
    });
    expect(result.added).toBe(1);
    expect(result.name).toMatch(/aggregat/i);
    expect(result.rows).toBeGreaterThan(0);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
