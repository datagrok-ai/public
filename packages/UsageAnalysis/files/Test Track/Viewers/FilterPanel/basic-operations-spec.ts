import {test, expect} from '@playwright/test';

const baseUrl = 'https://dev.datagrok.ai';
const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Filter panel — Basic operations', async ({page}) => {
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

  // Phase 3: Open filters
  await page.evaluate(() => grok.shell.tv.getFiltersGroup());
  await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});

  // ==== Section 1: Filtering and resetting ====

  await softStep('1.1-1.3 Structure filter c1ccccc1 → 32 rows', async () => {
    await page.locator('[name="viewer-Filters"] .sketch-link').first().click();
    await page.locator('.d4-dialog').waitFor({timeout: 5000});
    const smilesInput = page.locator('.d4-dialog input[placeholder*="SMILES"]');
    await smilesInput.click();
    await smilesInput.pressSequentially('c1ccccc1', {delay: 30});
    await smilesInput.press('Enter');
    await page.waitForTimeout(1500);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForTimeout(2000);
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBe(32);
  });

  await softStep('1.4-1.5 Set Stereo Category to R_ONE → 15 rows', async () => {
    await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: ['R_ONE']});
      await new Promise(r => setTimeout(r, 500));
    });
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBe(15);
  });

  await softStep('1.6-1.7 Set Average Mass max=400 → ~4 rows', async () => {
    await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      const col = grok.shell.tv.dataFrame.col('Average Mass');
      fg.updateOrAdd({type: 'histogram', column: 'Average Mass', min: col.min, max: 400});
      await new Promise(r => setTimeout(r, 500));
    });
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBe(4);
  });

  await softStep('1.8 Disable all filters globally → 100 rows', async () => {
    await page.evaluate(async () => {
      const checkbox = document.querySelector('[name="viewer-Filters"] .d4-filter-group-header input[type="checkbox"]') as HTMLInputElement;
      checkbox?.click();
      await new Promise(r => setTimeout(r, 500));
    });
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBe(100);
  });

  await softStep('1.9 Re-enable all filters → 4 rows restored', async () => {
    await page.evaluate(async () => {
      const checkbox = document.querySelector('[name="viewer-Filters"] .d4-filter-group-header input[type="checkbox"]') as HTMLInputElement;
      checkbox?.click();
      await new Promise(r => setTimeout(r, 500));
    });
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBe(4);
  });

  await softStep('1.10 Disable Stereo Category → count increases to 9', async () => {
    await page.evaluate(async () => {
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      for (const card of cards) {
        const header = card.querySelector('.d4-filter-header');
        if (header?.textContent?.trim()?.startsWith('Stereo Category')) {
          const checkbox = card.querySelector('input[type="checkbox"]') as HTMLInputElement;
          checkbox?.click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 500));
    });
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBe(9);
  });

  await softStep('1.11 Close Filter Panel → 100 rows', async () => {
    await page.evaluate(async () => {
      grok.shell.tv.getFiltersGroup().close();
      await new Promise(r => setTimeout(r, 500));
    });
    const result = await page.evaluate(() => ({
      open: !!document.querySelector('[name="viewer-Filters"]'),
      rows: grok.shell.tv.dataFrame.filter.trueCount,
    }));
    expect(result.open).toBe(false);
    expect(result.rows).toBe(100);
  });

  await softStep('1.12 Reopen — Stereo Category still disabled', async () => {
    await page.evaluate(async () => {
      grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 2000));
    });
    const checked = await page.evaluate(() => {
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      for (const card of cards) {
        const header = card.querySelector('.d4-filter-header');
        if (header?.textContent?.trim()?.startsWith('Stereo Category')) {
          const cb = card.querySelector('input[type="checkbox"]') as HTMLInputElement;
          return cb?.checked;
        }
      }
      return undefined;
    });
    expect(checked).toBe(false);
  });

  await softStep('1.13 Reset all filters → 100 rows', async () => {
    await page.evaluate(async () => {
      const header = document.querySelector('[name="viewer-Filters"] .d4-filter-group-header');
      const resetIcon = header?.querySelector('.fa-arrow-rotate-left') as HTMLElement;
      resetIcon?.click();
      await new Promise(r => setTimeout(r, 500));
      const okBtn = document.querySelector('.d4-dialog [name="button-OK"]') as HTMLElement;
      okBtn?.click();
      await new Promise(r => setTimeout(r, 500));
    });
    const count = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(count).toBe(100);
  });

  await softStep('1.14-1.15 Close/reopen — default state', async () => {
    const result = await page.evaluate(async () => {
      grok.shell.tv.getFiltersGroup().close();
      await new Promise(r => setTimeout(r, 500));
      grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 2000));
      return {
        rows: grok.shell.tv.dataFrame.filter.trueCount,
        cards: document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length,
      };
    });
    expect(result.rows).toBe(100);
    expect(result.cards).toBeGreaterThan(0);
  });

  await softStep('1.16-1.17 Remove Structure and Core filters', async () => {
    await page.evaluate(async () => {
      const toRemove = ['Structure', 'Core'];
      for (const name of toRemove) {
        const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
        for (const card of cards) {
          const h = card.querySelector('.d4-filter-header');
          if (h?.textContent?.trim() === name) {
            const closeBtn = card.querySelector('[name="icon-times"]') as HTMLElement;
            closeBtn?.click();
            await new Promise(r => setTimeout(r, 500));
            break;
          }
        }
      }
    });
    const names = await page.evaluate(() => {
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      return Array.from(cards).map(c => c.querySelector('.d4-filter-header')?.textContent?.trim());
    });
    expect(names).not.toContain('Structure');
    expect(names).not.toContain('Core');
  });

  await softStep('1.18 Close/reopen — removed filters absent', async () => {
    const result = await page.evaluate(async () => {
      grok.shell.tv.getFiltersGroup().close();
      await new Promise(r => setTimeout(r, 500));
      grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 2000));
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      const names = Array.from(cards).map(c => c.querySelector('.d4-filter-header')?.textContent?.trim());
      return {names, count: cards.length};
    });
    expect(result.names).not.toContain('Structure');
    expect(result.names).not.toContain('Core');
  });

  // ==== Section 2: Adding and Reordering Filters ====

  await softStep('2.1 Hamburger → Remove All', async () => {
    const count = await page.evaluate(async () => {
      const allH = document.querySelectorAll('[name="icon-font-icon-menu"]');
      for (const h of allH) {
        if (!h.closest('[name="viewer-Filters"]') && h.parentElement?.className?.includes('panel-titlebar')) {
          (h as HTMLElement).click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 500));
      const items = document.querySelectorAll('.d4-menu-item-label');
      for (const item of items) {
        if (item.textContent?.trim() === 'Remove All') {
          (item.closest('.d4-menu-item') as HTMLElement).click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 500));
      (document.querySelector('.d4-dialog [name="button-OK"]') as HTMLElement)?.click();
      await new Promise(r => setTimeout(r, 500));
      return document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length;
    });
    expect(count).toBe(0);
  });

  await softStep('2.2-2.5 Add 4 filters via different methods', async () => {
    const names = await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      // Method 1: ID (drag column header equivalent)
      fg.updateOrAdd({type: 'histogram', column: 'Id'});
      await new Promise(r => setTimeout(r, 500));
      // Method 2: CAST Idea ID (column header menu equivalent)
      fg.updateOrAdd({type: 'histogram', column: 'CAST Idea ID'});
      await new Promise(r => setTimeout(r, 500));
      // Method 3: Structure (right-click cell → Use as filter equivalent)
      fg.updateOrAdd({type: 'Chem:substructureFilter', column: 'Structure', columnName: 'Structure'});
      await new Promise(r => setTimeout(r, 1000));
      // Method 4: Scaffold Tree Filter (context menu → Add filter equivalent)
      fg.updateOrAdd({
        type: 'Chem:scaffoldTreeFilter', column: 'Structure', columnName: 'Structure',
        active: true, savedTree: JSON.stringify([]),
        colorCodedScaffolds: '[]', title: 'Scaffold Tree',
      });
      await new Promise(r => setTimeout(r, 1000));
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      return Array.from(cards).map(c => c.querySelector('.d4-filter-header')?.textContent?.trim());
    });
    expect(names.length).toBe(4);
  });

  await softStep('2.6 Verify filter order', async () => {
    const types = await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      return Array.from(fg.filters).map(f => ({
        col: f.columnName,
        cls: f.root?.className?.substring(0, 50),
      }));
    });
    // Expected: Scaffold Tree, Structure (substructure), CAST Idea ID, Id
    expect(types.length).toBe(4);
    expect(types[0].col).toBe('Structure'); // Scaffold Tree
    expect(types[1].col).toBe('Structure'); // Substructure
  });

  await softStep('2.7 Hamburger → Remove All', async () => {
    const count = await page.evaluate(async () => {
      const allH = document.querySelectorAll('[name="icon-font-icon-menu"]');
      for (const h of allH) {
        if (!h.closest('[name="viewer-Filters"]') && h.parentElement?.className?.includes('panel-titlebar')) {
          (h as HTMLElement).click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 500));
      const items = document.querySelectorAll('.d4-menu-item-label');
      for (const item of items) {
        if (item.textContent?.trim() === 'Remove All') {
          (item.closest('.d4-menu-item') as HTMLElement).click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 500));
      (document.querySelector('.d4-dialog [name="button-OK"]') as HTMLElement)?.click();
      await new Promise(r => setTimeout(r, 500));
      return document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length;
    });
    expect(count).toBe(0);
  });

  await softStep('2.8 Close Filter Panel', async () => {
    await page.evaluate(async () => {
      grok.shell.tv.getFiltersGroup().close();
      await new Promise(r => setTimeout(r, 500));
    });
    const open = await page.evaluate(() => !!document.querySelector('[name="viewer-Filters"]'));
    expect(open).toBe(false);
  });

  // ==== Section 3: Hidden Columns ====

  await softStep('3.1-3.3 Hide Structure/Core/R1 → verify grid', async () => {
    const result = await page.evaluate(async () => {
      const grid = grok.shell.tv.grid;
      grid.columns.byName('Structure').visible = false;
      grid.columns.byName('Core').visible = false;
      grid.columns.byName('R1').visible = false;
      grid.invalidate();
      await new Promise(r => setTimeout(r, 500));
      const visible: string[] = [];
      for (let i = 0; i < grid.columns.length; i++) {
        const gc = grid.columns.byIndex(i);
        if (gc.visible) visible.push(gc.name);
      }
      return {
        hasStructure: visible.includes('Structure'),
        hasCore: visible.includes('Core'),
        hasR1: visible.includes('R1'),
      };
    });
    expect(result.hasStructure).toBe(false);
    expect(result.hasCore).toBe(false);
    expect(result.hasR1).toBe(false);
  });

  await softStep('3.4 Open Filter Panel → no hidden column cards', async () => {
    const result = await page.evaluate(async () => {
      grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 2000));
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      const names = Array.from(cards).map(c => c.querySelector('.d4-filter-header')?.textContent?.trim());
      return {
        hasStructure: names.some(n => n === 'Structure'),
        hasCore: names.some(n => n === 'Core'),
        hasR1: names.some(n => n === 'R1'),
      };
    });
    expect(result.hasStructure).toBe(false);
    expect(result.hasCore).toBe(false);
    expect(result.hasR1).toBe(false);
  });

  await softStep('3.5 Restore columns', async () => {
    await page.evaluate(async () => {
      grok.shell.tv.getFiltersGroup().close();
      await new Promise(r => setTimeout(r, 500));
      const grid = grok.shell.tv.grid;
      grid.columns.byName('Structure').visible = true;
      grid.columns.byName('Core').visible = true;
      grid.columns.byName('R1').visible = true;
      grid.invalidate();
      await new Promise(r => setTimeout(r, 500));
    });
    const result = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      return {
        hasStructure: grid.columns.byName('Structure')?.visible,
        hasCore: grid.columns.byName('Core')?.visible,
        hasR1: grid.columns.byName('R1')?.visible,
      };
    });
    expect(result.hasStructure).toBe(true);
    expect(result.hasCore).toBe(true);
    expect(result.hasR1).toBe(true);
  });

  await softStep('3.6-3.7 Remove All, close Filter Panel', async () => {
    await page.evaluate(async () => {
      grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 2000));
      const allH = document.querySelectorAll('[name="icon-font-icon-menu"]');
      for (const h of allH) {
        if (!h.closest('[name="viewer-Filters"]') && h.parentElement?.className?.includes('panel-titlebar')) {
          (h as HTMLElement).click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 500));
      const items = document.querySelectorAll('.d4-menu-item-label');
      for (const item of items) {
        if (item.textContent?.trim() === 'Remove All') {
          (item.closest('.d4-menu-item') as HTMLElement).click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 500));
      (document.querySelector('.d4-dialog [name="button-OK"]') as HTMLElement)?.click();
      await new Promise(r => setTimeout(r, 500));
      grok.shell.tv.getFiltersGroup().close();
      await new Promise(r => setTimeout(r, 500));
    });
    const open = await page.evaluate(() => !!document.querySelector('[name="viewer-Filters"]'));
    expect(open).toBe(false);
  });

  await softStep('3.8 Reopen → Structure/Core/R1 cards present', async () => {
    const result = await page.evaluate(async () => {
      grok.shell.tv.getFiltersGroup();
      await new Promise(r => setTimeout(r, 2000));
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      const names = Array.from(cards).map(c => c.querySelector('.d4-filter-header')?.textContent?.trim());
      return {
        hasStructure: names.some(n => n === 'Structure'),
        hasCore: names.some(n => n === 'Core'),
        hasR1: names.some(n => n === 'R1'),
        total: cards.length,
      };
    });
    expect(result.hasStructure).toBe(true);
    expect(result.hasCore).toBe(true);
    expect(result.hasR1).toBe(true);
    expect(result.total).toBeGreaterThanOrEqual(42);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
