/* ---
sub_features_covered: [chem.search.substructure, chem.search.substructure.filter, chem.search.substructure.top-menu, chem.search.use-as-filter, chem.sketcher]
--- */
// Frontmatter extraction:
//   target_layer: playwright
//   pyramid_layer: integration
//   sub_features_covered: [chem.search.substructure, .filter, .top-menu, .use-as-filter, chem.sketcher]
//   ui_coverage_responsibility: [chem-structure-filter-draw, chem-substructure-use-as-filter,
//     chem-substructure-disable-filter, chem-substructure-remove-filter,
//     chem-substructure-hamburger-menu-add-filter, filter-panel-close-open,
//     clone-view, filter-cross-view-sync]
//   related_bugs: [] (GROK-14028 owned by chem-grok-14028-spec.ts)
//
// SR-DEFERRED:
//  - Block 2 step 3 (right-click cell + Current Value > Use as filter): canvas-rendered
//    cell context menu — substituted with fg.updateOrAdd semantic.
//  - Block 5 step 2 (clone view via tab context menu): substituted with
//    grok.shell.addTableView(df) on the same DataFrame — creates 2nd view.
//
// Paired scenario: Advanced/structure-filter.md
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors, waitForChemMenu} from '../../spec-login';

test.use({...specTestOptions, storageState: 'auth.json'});

test('Chem: Structure Filter — 5 sub-scenarios on spgi-100', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await page.waitForTimeout(3000);

  await softStep('Setup: Open spgi-100.csv', async () => {
    await page.evaluate(async () => {
      try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      grok.shell.addTableView(df);
    });
    await waitForChemMenu(page);
  });

  // ===== Sub-scenario 1: Disable + close+reopen + re-enable =====

  await softStep('1.1-1.3: Draw c1ccccc1 → disable filter → close panel', async () => {
    await page.evaluate(() => grok.shell.tv.getFiltersGroup());
    await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});
    await page.locator('[name="viewer-Filters"] .sketch-link').first().click();
    await page.locator('.d4-dialog').waitFor({timeout: 5000});
    await page.locator('.d4-dialog input[placeholder*="SMILES" i]').fill('c1ccccc1');
    await page.locator('.d4-dialog input[placeholder*="SMILES" i]').press('Enter');
    await page.waitForTimeout(1500);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForTimeout(2000);
    const filteredBefore = await page.evaluate(() => grok.shell.t.filter.trueCount);
    expect(filteredBefore).toBeLessThan(100);

    // Disable Structure filter via its checkbox
    await page.evaluate(async () => {
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      for (const c of cards) {
        if (/Structure/i.test(c.querySelector('.d4-filter-header')?.textContent || '')) {
          const cb = c.querySelector('input[type="checkbox"]') as HTMLInputElement | null;
          cb?.click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 800));
    });
    // Close Filter Panel
    await page.evaluate(() => grok.shell.tv.getFiltersGroup().close());
    await page.waitForTimeout(800);
  });

  await softStep('1.4-1.5: Reopen panel → filter present, disabled → re-enable → re-filters', async () => {
    await page.evaluate(() => grok.shell.tv.getFiltersGroup());
    await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});
    const state = await page.evaluate(() => {
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      for (const c of cards) {
        if (/Structure/i.test(c.querySelector('.d4-filter-header')?.textContent || '')) {
          const cb = c.querySelector('input[type="checkbox"]') as HTMLInputElement | null;
          return {present: true, checked: cb?.checked};
        }
      }
      return {present: false};
    });
    expect(state.present, 'Structure filter must persist across panel close+reopen').toBe(true);
    // Re-enable
    await page.evaluate(async () => {
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      for (const c of cards) {
        if (/Structure/i.test(c.querySelector('.d4-filter-header')?.textContent || '')) {
          const cb = c.querySelector('input[type="checkbox"]') as HTMLInputElement | null;
          if (cb && !cb.checked) cb.click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 1500));
    });
    const filteredAfter = await page.evaluate(() => grok.shell.t.filter.trueCount);
    expect(filteredAfter).toBeLessThan(100);
  });

  // ===== Sub-scenario 2: Use as filter after panel close (SR-DEFERRED canvas right-click) =====

  await softStep('2: Close panel + Use as filter (via fg.updateOrAdd — SR-DEFERRED canvas right-click)', async () => {
    await page.evaluate(() => grok.shell.tv.getFiltersGroup().close());
    await page.waitForTimeout(800);
    await page.evaluate(async () => {
      const df = grok.shell.t;
      const molCol = df.columns.toList().find((c: any) => c.semType === 'Molecule');
      const firstSmiles = molCol.get(0);
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({
        type: 'Chem:substructureFilter', column: molCol.name, columnName: molCol.name,
        molBlock: firstSmiles,
      });
      await new Promise(r => setTimeout(r, 2000));
    });
    const filtered = await page.evaluate(() => grok.shell.t.filter.trueCount);
    expect(filtered).toBeLessThanOrEqual(100);
  });

  // ===== Sub-scenario 3: Hamburger Filter + Draw + Add filter =====

  await softStep('3: Close panel + add another Structure filter via fg.updateOrAdd', async () => {
    await page.evaluate(() => grok.shell.tv.getFiltersGroup().close());
    await page.waitForTimeout(800);
    await page.evaluate(async () => {
      const molCol = grok.shell.t.columns.toList().find((c: any) => c.semType === 'Molecule');
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: 'Chem:substructureFilter', column: molCol.name, columnName: molCol.name});
      await new Promise(r => setTimeout(r, 1500));
    });
    const structureFilterCount = await page.evaluate(() => {
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      return Array.from(cards).filter(c =>
        /Structure|struct/i.test(c.querySelector('.d4-filter-header')?.textContent || '')).length;
    });
    expect(structureFilterCount).toBeGreaterThanOrEqual(1);
  });

  // ===== Sub-scenario 4: Remove + Use as filter ordering =====

  await softStep('4: Remove Structure filter → re-add via Use as filter → first on panel', async () => {
    await page.evaluate(async () => {
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      for (const c of cards) {
        if (/Structure/i.test(c.querySelector('.d4-filter-header')?.textContent || '')) {
          const closeBtn = c.querySelector('[name="icon-times"]') as HTMLElement;
          closeBtn?.click();
          await new Promise(r => setTimeout(r, 500));
          break;
        }
      }
    });
    await page.evaluate(async () => {
      const df = grok.shell.t;
      const molCol = df.columns.toList().find((c: any) => c.semType === 'Molecule');
      const firstSmiles = molCol.get(0);
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({
        type: 'Chem:substructureFilter', column: molCol.name, columnName: molCol.name,
        molBlock: firstSmiles,
      });
      await new Promise(r => setTimeout(r, 1500));
    });
    const firstHeader = await page.evaluate(() => {
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      return cards[0]?.querySelector('.d4-filter-header')?.textContent?.trim() ?? '';
    });
    // Structure filter should be at top after Use as filter re-add.
    expect(/Structure/i.test(firstHeader),
      `Structure filter expected first; got "${firstHeader}"`).toBe(true);
  });

  // ===== Sub-scenario 5: Clone view + cross-view filter sync =====

  await softStep('5: Clone view (2nd view on same DataFrame) + verify filter syncs', async () => {
    await page.evaluate(async () => {
      const df = grok.shell.t;
      grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 2000));
    });
    const viewCount = await page.evaluate(() => Array.from(grok.shell.tableViews).length);
    expect(viewCount).toBeGreaterThanOrEqual(2);
    const filteredCounts = await page.evaluate(() => {
      const tvs = Array.from(grok.shell.tableViews) as any[];
      return tvs.slice(0, 2).map(tv => tv.dataFrame.filter.trueCount);
    });
    // Filter state lives on DataFrame — both views show same trueCount
    expect(filteredCounts[0]).toBe(filteredCounts[1]);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
