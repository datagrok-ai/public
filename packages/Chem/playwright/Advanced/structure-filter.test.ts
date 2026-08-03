/* ---
sub_features_covered: [chem.search.substructure, chem.search.substructure.filter, chem.search.substructure.top-menu, chem.search.use-as-filter, chem.sketcher]
--- */
// SR-DEFERRED:
//  - Block 2 step 3 (right-click cell + Current Value > Use as filter): canvas-rendered
//    cell context menu — substituted with fg.updateOrAdd semantic.
//  - Block 5 step 2 (clone view via tab context menu): substituted with
//    grok.shell.addTableView(df) on the same DataFrame — creates 2nd view.
//
// Paired scenario: Advanced/structure-filter.md
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

test('Chem: Structure Filter — 5 sub-scenarios on spgi-100', async ({page}) => {
  test.setTimeout(120_000);

  let filteredBefore = 0;

  await loginToDatagrok(page);
  await page.waitForFunction(() => (window as any).grok?.shell != null, null, {timeout: 15000});

  await softStep('Setup: Open spgi-100.csv', async () => {
    await page.evaluate(async () => {
      try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      grok.shell.addTableView(df);
    });
    await waitForChemMenu(page);
    await waitForMolecule(page);
    // Load the Chem package bundle so `Chem:substructureFilter` resolves for the synchronous
    // filter operations below (fg.updateOrAdd) — the semType detector alone does not load it.
    await page.evaluate(async () => { await grok.functions.call('Chem:getRdKitModule', {}); });
  });

  // ===== Sub-scenario 1: Disable + close+reopen + re-enable =====

  await softStep('1.1-1.3: Draw c1ccccc1 → disable filter → close panel', async () => {
    // fg.updateOrAdd substitutes the hidden sketch-link → sketcher UI (SR-DEFERRED header).
    await page.evaluate(() => {
      const molCol = grok.shell.t.columns.toList().find((c: any) => c.semType === 'Molecule');
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({
        type: 'Chem:substructureFilter', column: molCol.name, columnName: molCol.name,
        molBlock: 'c1ccccc1',
      });
    });
    await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});
    await page.waitForFunction(() => grok.shell.t.filter.trueCount < grok.shell.t.rowCount,
      null, {timeout: 15000});
    filteredBefore = await page.evaluate(() => grok.shell.t.filter.trueCount);
    expect(filteredBefore, 'benzene substructure must match some but not all rows').toBeGreaterThan(0);
    expect(filteredBefore).toBeLessThan(100);

    // Disable Structure filter via its checkbox
    await page.evaluate(() => {
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      for (const c of cards) {
        if (/Structure/i.test(c.querySelector('.d4-filter-header')?.textContent || '')) {
          const cb = c.querySelector('input[type="checkbox"]') as HTMLInputElement | null;
          cb?.click();
          break;
        }
      }
    });
    // Disabling the filter must un-filter the table (all rows pass again)
    await page.waitForFunction(() => grok.shell.t.filter.trueCount === grok.shell.t.rowCount,
      null, {timeout: 10000});
    // Close Filter Panel
    await page.evaluate(() => grok.shell.tv.getFiltersGroup().close());
    await page.locator('[name="viewer-Filters"]').first().waitFor({state: 'detached', timeout: 10000});
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
    expect(state.checked, 'Structure filter must remain DISABLED after panel reopen').toBe(false);
    // Re-enable
    await page.evaluate(() => {
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      for (const c of cards) {
        if (/Structure/i.test(c.querySelector('.d4-filter-header')?.textContent || '')) {
          const cb = c.querySelector('input[type="checkbox"]') as HTMLInputElement | null;
          if (cb && !cb.checked) cb.click();
          break;
        }
      }
    });
    await page.waitForFunction(() => grok.shell.t.filter.trueCount < grok.shell.t.rowCount,
      null, {timeout: 15000});
    const filteredAfter = await page.evaluate(() => grok.shell.t.filter.trueCount);
    // Re-filtered set must be identical to the pre-disable set (.md step 5 round-trip)
    expect(filteredAfter, 'Re-enabled filter must restore the pre-disable row set').toBe(filteredBefore);
  });

  // ===== Sub-scenario 2: Use as filter after panel close (SR-DEFERRED canvas right-click) =====

  await softStep('2: Close panel + Use as filter (via fg.updateOrAdd — SR-DEFERRED canvas right-click)', async () => {
    await page.evaluate(() => grok.shell.tv.getFiltersGroup().close());
    await page.locator('[name="viewer-Filters"]').first().waitFor({state: 'detached', timeout: 10000});
    await page.evaluate(() => {
      const df = grok.shell.t;
      const molCol = df.columns.toList().find((c: any) => c.semType === 'Molecule');
      const firstSmiles = molCol.get(0);
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({
        type: 'Chem:substructureFilter', column: molCol.name, columnName: molCol.name,
        molBlock: firstSmiles,
      });
    });
    // Row 0's own molecule matches itself (>0) yet excludes non-matches (<100)
    await page.waitForFunction(() =>
      grok.shell.t.filter.trueCount > 0 && grok.shell.t.filter.trueCount < grok.shell.t.rowCount,
    null, {timeout: 15000});
    const filtered = await page.evaluate(() => grok.shell.t.filter.trueCount);
    expect(filtered).toBeGreaterThan(0);
    expect(filtered).toBeLessThan(100);
  });

  // ===== Sub-scenario 3: Hamburger Filter + Draw + Add filter =====

  await softStep('3: Close panel + draw a distinct structure → table re-filters to the new predicate', async () => {
    await page.evaluate(() => grok.shell.tv.getFiltersGroup().close());
    await page.locator('[name="viewer-Filters"]').first().waitFor({state: 'detached', timeout: 10000});
    // Distinct non-empty substructure (pyridine) so this is a real re-filter, not a no-op.
    // NOTE: fg.updateOrAdd updates the existing Structure card in place; the .md's
    // "Add filter" 1→2 additional-card combination lives on the canvas hamburger
    // "Add filter" UI which has no stable selector (run.md) — deferred for human review.
    await page.evaluate(() => {
      const molCol = grok.shell.t.columns.toList().find((c: any) => c.semType === 'Molecule');
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({
        type: 'Chem:substructureFilter', column: molCol.name, columnName: molCol.name,
        molBlock: 'c1ccncc1',
      });
    });
    await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});
    await page.waitForFunction(() =>
      grok.shell.t.filter.trueCount > 0 && grok.shell.t.filter.trueCount < grok.shell.t.rowCount,
    null, {timeout: 15000});
    const structureFilterCount = await page.evaluate(() => {
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      return Array.from(cards).filter(c =>
        /Structure|struct/i.test(c.querySelector('.d4-filter-header')?.textContent || '')).length;
    });
    expect(structureFilterCount, 'a Structure filter card must be present').toBeGreaterThanOrEqual(1);
    const filtered = await page.evaluate(() => grok.shell.t.filter.trueCount);
    expect(filtered).toBeGreaterThan(0);
    expect(filtered).toBeLessThan(100);
  });

  // ===== Sub-scenario 4: Remove + Use as filter ordering =====

  await softStep('4: Remove Structure filter → re-add via Use as filter → first on panel', async () => {
    await page.evaluate(() => {
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      for (const c of cards) {
        if (/Structure/i.test(c.querySelector('.d4-filter-header')?.textContent || '')) {
          const closeBtn = c.querySelector('[name="icon-times"]') as HTMLElement;
          closeBtn?.click();
          break;
        }
      }
    });
    // Removal: no Structure card should remain on the panel
    await page.waitForFunction(() => {
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      return !Array.from(cards).some(c =>
        /Structure/i.test(c.querySelector('.d4-filter-header')?.textContent || ''));
    }, null, {timeout: 10000});
    await page.evaluate(() => {
      const df = grok.shell.t;
      const molCol = df.columns.toList().find((c: any) => c.semType === 'Molecule');
      const firstSmiles = molCol.get(0);
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({
        type: 'Chem:substructureFilter', column: molCol.name, columnName: molCol.name,
        molBlock: firstSmiles,
      });
    });
    // Re-added Structure filter must land first on the panel (.md block 4 ordering)
    await page.waitForFunction(() => {
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      return /Structure/i.test(cards[0]?.querySelector('.d4-filter-header')?.textContent || '');
    }, null, {timeout: 15000});
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
    await page.evaluate(() => grok.shell.addTableView(grok.shell.t));
    await page.waitForFunction(() => Array.from(grok.shell.tableViews).length >= 2,
      null, {timeout: 15000});
    const viewCount = await page.evaluate(() => Array.from(grok.shell.tableViews).length);
    expect(viewCount).toBeGreaterThanOrEqual(2);
    // Apply a FRESH substructure filter through the cloned (2nd) view; the other view
    // must witness the same filtered subset (filter state lives on the shared DataFrame).
    await page.evaluate(() => {
      const tvs = Array.from(grok.shell.tableViews) as any[];
      const molCol = tvs[1].dataFrame.columns.toList().find((c: any) => c.semType === 'Molecule');
      tvs[1].getFiltersGroup().updateOrAdd({
        type: 'Chem:substructureFilter', column: molCol.name, columnName: molCol.name,
        molBlock: 'c1ccccc1',
      });
    });
    await page.waitForFunction(() => {
      const tv = (Array.from(grok.shell.tableViews) as any[])[0];
      return tv.dataFrame.filter.trueCount > 0 && tv.dataFrame.filter.trueCount < tv.dataFrame.rowCount;
    }, null, {timeout: 15000});
    const sync = await page.evaluate(() => {
      const tvs = (Array.from(grok.shell.tableViews) as any[]).slice(0, 2);
      return {
        counts: tvs.map(tv => tv.dataFrame.filter.trueCount),
        sameDf: tvs[0].dataFrame === tvs[1].dataFrame,
        rowCount: tvs[0].dataFrame.rowCount,
      };
    });
    expect(sync.sameDf, 'both views must reference the same DataFrame instance').toBe(true);
    expect(sync.counts[0], 'cross-view filtered count must be non-trivial').toBeGreaterThan(0);
    expect(sync.counts[0]).toBeLessThan(sync.rowCount);
    expect(sync.counts[0], 'both views must witness the same filtered subset').toBe(sync.counts[1]);
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
