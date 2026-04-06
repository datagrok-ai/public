import {test, expect} from '@playwright/test';

const baseUrl = 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/SPGI.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Filter panel — Cloned views', async ({page}) => {
  test.setTimeout(600_000);

  // Phase 1: Navigate
  await page.goto(baseUrl);
  await page.waitForFunction(() => typeof grok !== 'undefined' && grok.shell, {timeout: 15000});

  // Phase 2: Open dataset
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 5000);
    });
    // Wait for Chem cell rendering + package filter registration
    const hasBioChem = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
      .some(c => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 5000));
    }
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Phase 3: Close Context Panel and Context Help, open filters
  await page.evaluate(() => {
    grok.shell.windows.showContextPanel = false;
    grok.shell.windows.showHelp = false;
    grok.shell.tv.getFiltersGroup();
  });
  await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});

  await softStep('Set Competition assay filter + Stereo Category + Structure', async () => {
    await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.HISTOGRAM, column: 'Competition assay', filterOutMissingValues: true});
      await new Promise(r => setTimeout(r, 500));
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: ['S_ACHIR']});
      await new Promise(r => setTimeout(r, 500));
    });

    // Open sketcher for structure
    await page.locator('[name="viewer-Filters"] .sketch-link').first().click();
    await page.locator('.d4-dialog').waitFor({timeout: 5000});
    const smilesInput = page.locator('.d4-dialog input[placeholder*="SMILES"]');
    await smilesInput.fill('c1ccccc1');
    await smilesInput.press('Enter');
    await page.waitForTimeout(1000);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForTimeout(2000);

    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBe(391);
  });

  await softStep('Clone View', async () => {
    await page.locator('[name="div-View"]').click();
    await page.waitForTimeout(500);
    await page.locator('[name="div-View---Layout"]').hover();
    await page.waitForTimeout(300);
    await page.locator('[name="div-View---Layout---Clone-View"]').click();
    await page.waitForTimeout(3000);

    const viewName = await page.evaluate(() => grok.shell.tv.name);
    expect(viewName).toContain('copy');
  });

  await softStep('Verify cloned view filter state', async () => {
    const result = await page.evaluate(() => ({
      count: grok.shell.tv.dataFrame.filter.trueCount,
      panelOpen: !!document.querySelector('[name="viewer-Filters"]'),
    }));
    expect(result.count).toBe(391);
    expect(result.panelOpen).toBe(true);
  });

  await softStep('Toggle all filters off/on', async () => {
    const result = await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      const filters = fg.filters;
      for (let i = 0; i < filters.length; i++) fg.setEnabled(filters[i], false);
      await new Promise(r => setTimeout(r, 1000));
      const countOff = grok.shell.tv.dataFrame.filter.trueCount;
      for (let i = 0; i < filters.length; i++) fg.setEnabled(filters[i], true);
      await new Promise(r => setTimeout(r, 1000));
      const countOn = grok.shell.tv.dataFrame.filter.trueCount;
      return {countOff, countOn};
    });
    expect(result.countOff).toBe(3624);
    expect(result.countOn).toBe(391);
  });

  await softStep('Clear Structure, set C1CCCCC1', async () => {
    await page.evaluate(async () => {
      const fv = document.querySelector('[name="viewer-Filters"]');
      const clearBtn = fv?.querySelector('.chem-clear-sketcher-button') as HTMLElement;
      clearBtn?.click();
      await new Promise(r => setTimeout(r, 1000));
    });

    await page.locator('[name="viewer-Filters"] .sketch-link').first().click();
    await page.locator('.d4-dialog').waitFor({timeout: 5000});
    const smilesInput = page.locator('.d4-dialog input[placeholder*="SMILES"]');
    await smilesInput.fill('C1CCCCC1');
    await smilesInput.press('Enter');
    await page.waitForTimeout(1000);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForTimeout(2000);

    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBe(10);
  });

  await softStep('Remove Structure filter from clone', async () => {
    await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      const filters = fg.filters;
      for (let i = 0; i < filters.length; i++) {
        if (filters[i].column?.name === 'Structure' || filters[i].columnName === 'Structure') {
          fg.remove(filters[i]);
          break;
        }
      }
      await new Promise(r => setTimeout(r, 1000));
    });
  });

  await softStep('Save layout, close, apply — verify no Structure filter', async () => {
    const result = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(r => setTimeout(r, 1000));
      const fg = grok.shell.tv.getFiltersGroup();
      fg.close();
      await new Promise(r => setTimeout(r, 1000));
      const saved = await grok.dapi.layouts.find(layoutId);
      grok.shell.tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 3000));
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      let structurePresent = false;
      cards.forEach(c => {
        const t = c.querySelector('.d4-filter-header')?.textContent?.trim();
        if (t?.includes('Structure')) structurePresent = true;
      });
      try { await grok.dapi.layouts.delete(saved); } catch (e) {}
      return {panelOpen: !!document.querySelector('[name="viewer-Filters"]'), structurePresent};
    });
    expect(result.panelOpen).toBe(true);
    expect(result.structurePresent).toBe(false);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
