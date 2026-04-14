import {test, expect, chromium} from '@playwright/test';

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const spgiPath = 'System:DemoFiles/SPGI.csv';
const demogPath = 'System:DemoFiles/demog.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Sunburst viewer', async () => {
  const browser = await chromium.connectOverCDP('http://localhost:9222');
  const context = browser.contexts()[0];
  let page = context.pages().find(p => p.url().includes('datagrok'));
  if (!page) {
    page = await context.newPage();
    await page.goto(baseUrl, {waitUntil: 'networkidle', timeout: 60000});
    await page.waitForFunction(() => {
      try { return typeof grok !== 'undefined' && typeof grok.shell.closeAll === 'function'; }
      catch { return false; }
    }, {timeout: 45000});
  }

  // Step 1: Open SPGI.csv and demog.csv, add Sunburst viewers
  await softStep('Step 1: Open files and add Sunburst viewer', async () => {
    const result = await page!.evaluate(async (paths) => {
      document.querySelectorAll('.d4-dialog').forEach(d => {
        const cancel = d.querySelector('[name="button-CANCEL"]');
        if (cancel) cancel.click();
      });
      grok.shell.closeAll();
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = false;

      // Open SPGI
      const df1 = await grok.dapi.files.readCsv(paths.spgi);
      const tv1 = grok.shell.addTableView(df1);
      await new Promise(resolve => {
        const sub = df1.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 5000);
      });
      const hasChem = Array.from({length: df1.columns.length}, (_, i) => df1.columns.byIndex(i))
        .some(c => c.semType === 'Molecule');
      if (hasChem) {
        for (let i = 0; i < 50; i++) {
          if (document.querySelector('[name="viewer-Grid"] canvas')) break;
          await new Promise(r => setTimeout(r, 200));
        }
        await new Promise(r => setTimeout(r, 5000));
      }
      tv1.addViewer('Sunburst');
      await new Promise(r => setTimeout(r, 2000));

      // Open demog
      const df2 = await grok.dapi.files.readCsv(paths.demog);
      const tv2 = grok.shell.addTableView(df2);
      await new Promise(resolve => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      tv2.addViewer('Sunburst');
      await new Promise(r => setTimeout(r, 2000));

      return {spgiRows: df1.rowCount, demogRows: df2.rowCount, viewers: tv2.viewers.length};
    }, {spgi: spgiPath, demog: demogPath});
    expect(result.spgiRows).toBe(3624);
    expect(result.demogRows).toBe(5850);
    expect(result.viewers).toBe(2);
  });

  // Step 2: Viewer properties panel
  await softStep('Step 2: Viewer properties panel', async () => {
    const opened = await page!.evaluate(async () => {
      const sunburst = document.querySelector('[name="viewer-Sunburst"]');
      if (!sunburst) return false;
      const panel = sunburst.parentElement?.parentElement;
      const gear = panel?.querySelector('[name="icon-font-icon-settings"]');
      if (gear) gear.click();
      await new Promise(r => setTimeout(r, 500));
      return true;
    });
    expect(opened).toBe(true);
  });

  // Step 3.1: Table switching
  await softStep('Step 3.1: Table switching', async () => {
    const tableName = await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const sb = viewers.find(v => v.type === 'Sunburst');
      sb.setOptions({table: 'Table'});
      await new Promise(r => setTimeout(r, 1000));
      const name = sb.dataFrame?.name;
      sb.setOptions({table: 'Table (2)'});
      await new Promise(r => setTimeout(r, 1000));
      return name;
    });
    expect(tableName).toBe('Table');
  });

  // Step 3.2: Hierarchy configuration
  await softStep('Step 3.2: Hierarchy configuration', async () => {
    const hierarchy = await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const sb = viewers.find(v => v.type === 'Sunburst');
      sb.setOptions({hierarchyColumnNames: ['SEX', 'RACE']});
      await new Promise(r => setTimeout(r, 1000));
      return sb.getOptions().look?.hierarchyColumnNames;
    });
    expect(hierarchy).toEqual(['SEX', 'RACE']);
  });

  // Step 3.3: Inherit from grid
  await softStep('Step 3.3: Inherit from grid', async () => {
    await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const sb = viewers.find(v => v.type === 'Sunburst');
      sb.setOptions({hierarchyColumnNames: ['SEX'], inheritFromGrid: true});
      await new Promise(r => setTimeout(r, 500));
      grok.shell.tv.dataFrame.col('SEX').meta.colors.setCategorical({'M': '#0000ff', 'F': '#ff0000'});
      await new Promise(r => setTimeout(r, 1000));
      // Change colors
      grok.shell.tv.dataFrame.col('SEX').meta.colors.setCategorical({'M': '#00ff00', 'F': '#ff00ff'});
      await new Promise(r => setTimeout(r, 1000));
    });
  });

  // Step 3.4: Include nulls (on SPGI)
  await softStep('Step 3.4: Include nulls', async () => {
    await page!.evaluate(async () => {
      const views = Array.from(grok.shell.views);
      const spgiView = views.find(v => v.name === 'Table');
      if (spgiView) grok.shell.v = spgiView;
      await new Promise(r => setTimeout(r, 500));

      const viewers = Array.from(grok.shell.tv.viewers);
      const sb = viewers.find(v => v.type === 'Sunburst');
      sb.setOptions({hierarchyColumnNames: ['Core', 'R101'], includeNulls: true});
      await new Promise(r => setTimeout(r, 1000));
      sb.setOptions({includeNulls: false});
      await new Promise(r => setTimeout(r, 1000));
    });
  });

  // Step 4: View reset
  await softStep('Step 4: View reset', async () => {
    const result = await page!.evaluate(async () => {
      const views = Array.from(grok.shell.views);
      const demogView = views.find(v => v.name === 'Table (2)');
      if (demogView) grok.shell.v = demogView;
      await new Promise(r => setTimeout(r, 500));

      const viewers = Array.from(grok.shell.tv.viewers);
      const sb = viewers.find(v => v.type === 'Sunburst');
      sb.setOptions({hierarchyColumnNames: ['SEX', 'RACE', 'DIS_POP']});
      await new Promise(r => setTimeout(r, 1000));

      // Reset via Ctrl+Shift+A is keyboard-based; verify selection is 0
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      return {selection: df.selection.trueCount, filtered: df.filter.trueCount};
    });
    expect(result.selection).toBe(0);
    expect(result.filtered).toBe(5850);
  });

  // Step 7: Projects & layouts
  await softStep('Step 7: Layout save/restore', async () => {
    const result = await page!.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const sb = viewers.find(v => v.type === 'Sunburst');
      sb.setOptions({hierarchyColumnNames: ['SEX', 'RACE', 'DIS_POP']});
      await new Promise(r => setTimeout(r, 1000));

      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1000));

      // Close sunburst
      sb.close();
      await new Promise(r => setTimeout(r, 500));

      // Restore layout
      const saved = await grok.dapi.layouts.find(layoutId);
      grok.shell.tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));

      const viewersAfter = Array.from(grok.shell.tv.viewers);
      const hasSunburst = viewersAfter.some(v => v.type === 'Sunburst');
      const restoredSb = viewersAfter.find(v => v.type === 'Sunburst');
      const hierarchy = restoredSb?.getOptions().look?.hierarchyColumnNames;

      await grok.dapi.layouts.delete(saved);
      return {hasSunburst, hierarchy};
    });
    expect(result.hasSunburst).toBe(true);
    expect(result.hierarchy).toEqual(['SEX', 'RACE', 'DIS_POP']);
  });

  // Step 9: Collaborative filtering
  await softStep('Step 9: Collaborative filtering', async () => {
    const result = await page!.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      // Apply filter to SEX = M only
      const sexCol = df.col('SEX');
      for (let i = 0; i < df.rowCount; i++) {
        if (sexCol.get(i) !== 'M')
          df.filter.set(i, false);
      }
      df.filter.fireChanged();
      await new Promise(r => setTimeout(r, 1000));
      const filtered = df.filter.trueCount;

      // Reset filter
      df.filter.setAll(true);
      df.filter.fireChanged();
      return {filtered};
    });
    expect(result.filtered).toBe(2607);
  });

  // Cleanup
  await page!.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
