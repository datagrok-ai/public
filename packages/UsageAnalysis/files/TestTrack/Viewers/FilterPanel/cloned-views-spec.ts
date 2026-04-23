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
const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) { stepErrors.push({step: name, error: e?.message ?? String(e)}); }
}

test('FilterPanel — Cloned Views', async ({page}) => {
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
      setTimeout(resolve, 5000);
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
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Phase 3: Open filter panel
  await page.evaluate(() => grok.shell.tv.getFiltersGroup());
  await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});

  let originalFilteredCount: number;

  await softStep('Set up filters: Competition assay missing values, Stereo Category, Structure', async () => {
    const result = await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      const df = grok.shell.tv.dataFrame;

      // Filter out missing values for Competition assay
      fg.updateOrAdd({type: DG.FILTER_TYPE.HISTOGRAM, column: 'Competition assay', filterOutMissingValues: true});
      // Set Stereo Category to S_ACHIR
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: ['S_ACHIR']});

      await new Promise(r => setTimeout(r, 2000));
      return {filteredCount: df.filter.trueCount, totalRows: df.rowCount};
    });
    expect(result.filteredCount).toBeLessThan(result.totalRows);

    // Set Structure filter via sketch link → SMILES input → OK
    await page.evaluate(() => {
      (document.querySelectorAll('.sketch-link')[0] as HTMLElement).click();
    });
    await page.locator('.d4-dialog input[placeholder*="SMILES"]').waitFor({timeout: 5000});
    const smilesInput = page.locator('.d4-dialog input[placeholder*="SMILES"]');
    await smilesInput.focus();
    await smilesInput.press('Control+A');
    await smilesInput.type('c1ccccc1', {delay: 30});
    await smilesInput.press('Enter');
    await page.waitForTimeout(2000);
    await page.locator('[name="button-OK"]').click();
    await page.waitForTimeout(5000);

    originalFilteredCount = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(originalFilteredCount).toBeLessThan(result.filteredCount);
  });

  await softStep('Clone View and verify filter state matches', async () => {
    // Use evaluate for menu interaction — Playwright hover doesn't reliably trigger Dart submenus
    await page.evaluate(async () => {
      document.querySelector('[name="div-View"]')?.click();
      await new Promise(r => setTimeout(r, 500));
      const layout = document.querySelector('[name="div-View---Layout"]');
      if (layout) {
        layout.dispatchEvent(new MouseEvent('mousemove', {bubbles: true}));
        layout.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      }
      await new Promise(r => setTimeout(r, 500));
      document.querySelector('[name="div-View---Layout---Clone-View"]')?.click();
    });
    await page.waitForTimeout(5000);

    const state = await page.evaluate(() => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      return {
        currentView: tv.name,
        filterPanelOpen: document.querySelector('[name="viewer-Filters"]') !== null,
        filteredCount: df.filter.trueCount,
        hasStructureCanvas: document.querySelector('[name="viewer-Filters"] .chem-external-sketcher-canvas') !== null,
      };
    });

    expect(state.currentView).toContain('copy');
    expect(state.filterPanelOpen).toBe(true);
    expect(state.filteredCount).toBe(originalFilteredCount);
    expect(state.hasStructureCanvas).toBe(true);
  });

  await softStep('Toggle all filters off — all rows visible', async () => {
    const result = await page.evaluate(async () => {
      const cb = document.querySelector('[name="viewer-Filters"] .d4-filter-group-header input[type="checkbox"]') as HTMLInputElement;
      cb.click();
      await new Promise(r => setTimeout(r, 2000));
      const df = grok.shell.tv.dataFrame;
      return {filteredCount: df.filter.trueCount, totalRows: df.rowCount};
    });
    expect(result.filteredCount).toBe(result.totalRows);
  });

  await softStep('Toggle all filters on — rows filtered again', async () => {
    const result = await page.evaluate(async () => {
      const cb = document.querySelector('[name="viewer-Filters"] .d4-filter-group-header input[type="checkbox"]') as HTMLInputElement;
      cb.click();
      await new Promise(r => setTimeout(r, 2000));
      const df = grok.shell.tv.dataFrame;
      return {filteredCount: df.filter.trueCount, totalRows: df.rowCount};
    });
    expect(result.filteredCount).toBe(originalFilteredCount);
  });

  await softStep('Clear Structure filter and set to C1CCCCC1', async () => {
    await page.evaluate(async () => {
      const clearBtn = document.querySelector('[name="viewer-Filters"] .chem-clear-sketcher-button') as HTMLElement;
      clearBtn.click();
      await new Promise(r => setTimeout(r, 2000));
    });

    await page.evaluate(() => {
      (document.querySelectorAll('.sketch-link')[0] as HTMLElement).click();
    });
    await page.locator('.d4-dialog input[placeholder*="SMILES"]').waitFor({timeout: 5000});
    const smilesInput2 = page.locator('.d4-dialog input[placeholder*="SMILES"]');
    await smilesInput2.focus();
    await smilesInput2.press('Control+A');
    await smilesInput2.type('C1CCCCC1', {delay: 30});
    await smilesInput2.press('Enter');
    await page.waitForTimeout(2000);
    await page.locator('[name="button-OK"]').click();
    await page.waitForTimeout(5000);

    const newCount = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(newCount).not.toBe(originalFilteredCount);
  });

  await softStep('Remove Structure filter — filtered state should not change', async () => {
    const beforeRemove = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);

    await page.evaluate(async () => {
      const filterPanel = document.querySelector('[name="viewer-Filters"]')!;
      const headers = filterPanel.querySelectorAll('.d4-filter-header');
      for (const header of headers) {
        if (header.textContent!.trim().startsWith('Structure')) {
          const card = header.closest('.d4-filter');
          const closeIcon = card?.querySelector('[name="icon-times"]') as HTMLElement;
          closeIcon?.click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 2000));
    });

    const afterRemove = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(afterRemove).toBe(beforeRemove);
  });

  let layoutId: string;

  await softStep('Save layout, close filter panel, apply layout', async () => {
    layoutId = await page.evaluate(async () => {
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1000));
      return layout.id;
    });

    const beforeClose = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);

    await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.close();
      await new Promise(r => setTimeout(r, 1000));
    });
    expect(await page.locator('[name="viewer-Filters"]').count()).toBe(0);

    const state = await page.evaluate(async (id) => {
      const saved = await grok.dapi.layouts.find(id);
      grok.shell.tv.loadLayout(saved);
      await new Promise(r => setTimeout(r, 5000));

      const filterPanelOpen = document.querySelector('[name="viewer-Filters"]') !== null;
      const headers: string[] = [];
      document.querySelectorAll('[name="viewer-Filters"] .d4-filter-header').forEach(h => {
        headers.push(h.textContent!.trim());
      });
      const hasStructure = headers.some(h => h.startsWith('Structure'));
      const df = grok.shell.tv.dataFrame;
      return {filterPanelOpen, hasStructure, filteredCount: df.filter.trueCount};
    }, layoutId);

    expect(state.filterPanelOpen).toBe(true);
    expect(state.hasStructure).toBe(false);
    expect(state.filteredCount).toBe(beforeClose);
  });

  // Cleanup
  await page.evaluate(async (id) => {
    try {
      const saved = await grok.dapi.layouts.find(id);
      await grok.dapi.layouts.delete(saved);
    } catch (_) {}
    grok.shell.closeAll();
  }, layoutId!);

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
