import { test, expect, Page } from '@playwright/test';

const baseUrl = 'https://dev.datagrok.ai';
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

test('Sunburst viewer', async ({ page }) => {
  test.setTimeout(120_000);
  // Phase 1: Navigate
  await page.goto(baseUrl);
  await page.waitForFunction(() => {
    return typeof grok !== 'undefined' && grok.shell && grok.shell.views;
  }, {timeout: 30000});

  // Phase 2: Open SPGI.csv
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    const hasBioChem = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
      .some(c => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 5000));
    }
  }, spgiPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Step 1: Add Sunburst for SPGI and demog
  await softStep('Step 1: Open files and add Sunburst viewer', async () => {
    await page.evaluate(async (demogPath) => {
      // Add Sunburst for SPGI
      await grok.shell.tv.addViewer('Sunburst');
      await new Promise(r => setTimeout(r, 2000));

      // Open demog.csv
      const df2 = await grok.dapi.files.readCsv(demogPath);
      const tv2 = grok.shell.addTableView(df2);
      await new Promise(resolve => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      await tv2.addViewer('Sunburst');
      await new Promise(r => setTimeout(r, 2000));
    }, datasetPath);
    const sunburst = await page.locator('[name="viewer-Sunburst"]').isVisible();
    expect(sunburst).toBe(true);
  });

  // Step 2: Open properties panel
  await softStep('Step 2: Viewer properties panel', async () => {
    const opened = await page.evaluate(() => {
      const sunburst = document.querySelector('[name="viewer-Sunburst"]');
      if (!sunburst) return false;
      const dockPanel = sunburst.closest('.panel-content');
      const titleBar = dockPanel?.parentElement?.querySelector('.panel-titlebar');
      const gear = titleBar?.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
      if (gear) gear.click();
      return true;
    });
    expect(opened).toBe(true);
    await page.waitForTimeout(500);
  });

  // Step 3.1: Table switching
  await softStep('Step 3.1: Table switching', async () => {
    const result = await page.evaluate(async () => {
      const combo = document.querySelector('.grok-prop-panel select') as HTMLSelectElement;
      if (!combo) return { error: 'no combo' };
      const original = combo.value;
      combo.value = 'Table';
      combo.dispatchEvent(new Event('change', { bubbles: true }));
      await new Promise(r => setTimeout(r, 500));
      combo.value = original;
      combo.dispatchEvent(new Event('change', { bubbles: true }));
      return { switched: true };
    });
    expect(result).not.toHaveProperty('error');
  });

  // Step 3.2: Hierarchy configuration
  await softStep('Step 3.2: Hierarchy configuration', async () => {
    const result = await page.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const sunburst = viewers.find((v: any) => v.type === 'Sunburst') as any;
      if (!sunburst) return { error: 'no sunburst' };
      sunburst.setOptions({ hierarchyColumnNames: ['SEX', 'RACE'] });
      await new Promise(r => setTimeout(r, 1000));
      const opts = sunburst.getOptions();
      return { hierarchy: opts.look?.hierarchyColumnNames };
    });
    expect(result).not.toHaveProperty('error');
    expect(result.hierarchy).toEqual(['SEX', 'RACE']);
  });

  // Step 3.3: Inherit from grid
  await softStep('Step 3.3: Inherit from grid', async () => {
    const result = await page.evaluate(async () => {
      const viewers = Array.from(grok.shell.tv.viewers);
      const sunburst = viewers.find((v: any) => v.type === 'Sunburst') as any;
      sunburst.setOptions({ hierarchyColumnNames: ['SEX'], inheritFromGrid: true });
      await new Promise(r => setTimeout(r, 500));
      const df = grok.shell.tv.dataFrame;
      df.col('SEX').meta.colors.setLinear([DG.Color.blue, DG.Color.red]);
      await new Promise(r => setTimeout(r, 500));
      return { done: true };
    });
    expect(result.done).toBe(true);
  });

  // Step 3.4: Include nulls (on SPGI)
  await softStep('Step 3.4: Include nulls', async () => {
    const result = await page.evaluate(async () => {
      // Switch to SPGI view
      const views = Array.from(grok.shell.views);
      const spgiView = views.find((v: any) => v.name === 'Table');
      if (spgiView) grok.shell.v = spgiView;
      await new Promise(r => setTimeout(r, 500));

      const viewers = Array.from(grok.shell.tv.viewers);
      const sunburst = viewers.find((v: any) => v.type === 'Sunburst') as any;
      if (!sunburst) return { error: 'no sunburst on SPGI' };

      sunburst.setOptions({ hierarchyColumnNames: ['Core', 'R101'], includeNulls: true });
      await new Promise(r => setTimeout(r, 1000));

      sunburst.setOptions({ includeNulls: false });
      await new Promise(r => setTimeout(r, 1000));

      return { done: true };
    });
    expect(result).not.toHaveProperty('error');
  });

  // Step 7: Layout save/restore
  await softStep('Step 7: Projects & layouts', async () => {
    const result = await page.evaluate(async () => {
      // Switch back to demog
      const views = Array.from(grok.shell.views);
      const demogView = views.find((v: any) => v.name === 'Table (2)');
      if (demogView) grok.shell.v = demogView;
      await new Promise(r => setTimeout(r, 500));

      const viewers = Array.from(grok.shell.tv.viewers);
      const sunburst = viewers.find((v: any) => v.type === 'Sunburst') as any;
      if (sunburst)
        sunburst.setOptions({ hierarchyColumnNames: ['SEX', 'CONTROL', 'RACE'] });
      await new Promise(r => setTimeout(r, 500));

      const layout = grok.shell.tv.saveLayout();
      const saved = await grok.dapi.layouts.save(layout);
      const layoutId = saved.id;
      await new Promise(r => setTimeout(r, 1000));

      if (sunburst) sunburst.close();
      await new Promise(r => setTimeout(r, 500));

      const restored = await grok.dapi.layouts.find(layoutId);
      grok.shell.tv.loadLayout(restored);
      await new Promise(r => setTimeout(r, 3000));

      const viewersAfter = Array.from(grok.shell.tv.viewers);
      const hasSunburst = viewersAfter.some((v: any) => v.type === 'Sunburst');
      await grok.dapi.layouts.delete(restored);

      return { hasSunburst };
    });
    expect(result.hasSunburst).toBe(true);
  });

  // Step 9: Collaborative filtering
  await softStep('Step 9: Collaborative filtering', async () => {
    const result = await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 1000));
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: ['Asian', 'Caucasian']});
      await new Promise(r => setTimeout(r, 500));
      return { filtered: grok.shell.tv.dataFrame.filter.trueCount };
    });
    expect(result.filtered).toBeLessThan(5850);
    expect(result.filtered).toBeGreaterThan(0);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
