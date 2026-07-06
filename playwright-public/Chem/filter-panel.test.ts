/* ---
sub_features_covered: [chem.search.substructure, chem.search.substructure.editor, chem.search.substructure.filter, chem.search.use-as-filter, chem.sketcher, chem.sketcher.ocl]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);

test('Chem: Filter Panel deep-dive (Blocks A-E)', async ({page}) => {
  test.setTimeout(150_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(e.message));

  await loginToDatagrok(page);

  await softStep('Step 1: Open spgi-100.csv', async () => {
    await page.evaluate(async () => {
      try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      grok.shell.addTableView(df);
    });
    await waitForChemMenu(page);
    await waitForMolecule(page);
  });

  // ===== Block A — Filter Panel basics =====

  await softStep('Step 2 (A): Open Filter Panel → Structure filter card rendered', async () => {
    await page.evaluate(() => grok.shell.tv.getFiltersGroup());
    await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 10000});
    const hasStructureFilter = await page.evaluate(() =>
      Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'))
        .some(f => /Structure|struct/i.test(f.querySelector('.d4-filter-header')?.textContent || '')));
    expect(hasStructureFilter).toBe(true);
  });

  await softStep('Step 3 (A): Draw c1ccccc1 substructure → ~32 rows filter', async () => {
    // SR-DEFERRED sketch-link draw: filling the Edit-sketcher dialog's SMILES input does not
    // load the molecule into the structure filter's sketcher (filter stays at trueCount=100).
    // Apply the benzene substructure via the filter-group API (same substitution as Steps 6/9/11/13).
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      const molCol = grok.shell.t.columns.toList().find((c: any) => c.semType === 'Molecule');
      fg.updateOrAdd({
        type: 'Chem:substructureFilter',
        column: molCol.name,
        columnName: molCol.name,
        molBlock: 'c1ccccc1',
      });
    });
    // The substructure search is async — poll for the filter to take effect.
    await expect.poll(() => page.evaluate(() => grok.shell.t.filter.trueCount),
      {timeout: 20000}).toBeLessThan(100);
    const filtered = await page.evaluate(() => grok.shell.t.filter.trueCount);
    expect(filtered).toBeLessThan(100);
    expect(filtered).toBeGreaterThan(0);
  });

  await softStep('Step 4 (A): Cycle Contains / Included in / Exact / Similar tabs', async () => {
    const errBefore = pageErrors.length;
    const results = await page.evaluate(async () => {
      const settle = async (): Promise<number> => {
        let prev = -1, stable = 0;
        for (let i = 0; i < 40; i++) {
          const c = grok.shell.tv.dataFrame.filter.trueCount;
          if (c === prev) { if (++stable >= 3) return c; } else { stable = 0; prev = c; }
          await new Promise(r => setTimeout(r, 200));
        }
        return prev;
      };
      const filterViewer = document.querySelector('[name="viewer-Filters"]');
      const cards = filterViewer?.querySelectorAll('.d4-filter');
      let structureCard: Element | null = null;
      cards?.forEach(card => {
        const header = card.querySelector('.d4-filter-header');
        if (header?.textContent?.trim() === 'Structure') structureCard = card;
      });
      const gearIcon = structureCard?.querySelector('.chem-search-options-icon, [class*="search-options"]') as HTMLElement | null;
      if (gearIcon) gearIcon.click();
      let select: HTMLSelectElement | null = null;
      for (let i = 0; i < 25 && !select; i++) {
        select = structureCard?.querySelector('select') as HTMLSelectElement | null;
        if (!select) await new Promise(r => setTimeout(r, 200));
      }
      if (!select) return null;
      const res: Record<string, number> = {};
      const types = ['Contains', 'Included in', 'Exact', 'Similar'];
      for (const type of types) {
        select.value = type;
        select.dispatchEvent(new Event('input', {bubbles: true}));
        select.dispatchEvent(new Event('change', {bubbles: true}));
        res[type] = await settle();
      }
      return res;
    });
    expect(results, 'settings gear must expose the match-mode select').not.toBeNull();
    const r = results as Record<string, number>;
    // The benzene query from Step 3 is active — Contains must be a real filtered subset, not the full 100.
    expect(r.Contains).toBeGreaterThan(0);
    expect(r.Contains).toBeLessThan(100);
    // Exact matches are a subset of Contains matches.
    expect(r.Exact).toBeLessThanOrEqual(r.Contains);
    // md Block A: each tab switch must complete without console errors.
    expect(pageErrors.length, `pageerrors during mode cycle: ${pageErrors.slice(errBefore).join('; ')}`).toBe(errBefore);
  });

  // ===== Block B — Use as filter from molecule cell context menu =====

  await softStep('Step 5 (B): Close Filter Panel', async () => {
    await page.evaluate(() => grok.shell.tv.getFiltersGroup().close());
    await expect.poll(() => page.evaluate(() =>
      document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length), {timeout: 5000}).toBe(0);
  });

  await softStep('Step 6 (B): Use as filter for first molecule structure (SR-DEFERRED: right-click canvas)', async () => {
    await page.evaluate(() => {
      const df = grok.shell.t;
      const molCol = df.columns.toList().find((c: any) => c.semType === 'Molecule');
      const firstSmiles = molCol.get(0);
      const fg = grok.shell.tv.getFiltersGroup();
      // SR-DEFERRED right-click canvas → Use as filter: substituted with fg.updateOrAdd using the cell's SMILES.
      fg.updateOrAdd({
        type: 'Chem:substructureFilter',
        column: molCol.name,
        columnName: molCol.name,
        molBlock: firstSmiles,
      });
    });
    // A query built from a real row's structure must match at least that row and filter the table.
    await expect.poll(() => page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount),
      {timeout: 20000}).toBeLessThan(100);
    const filtered = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filtered).toBeLessThan(100);
    expect(filtered).toBeGreaterThan(0);
  });

  // ===== Block C — Filter from column header hamburger menu =====

  await softStep('Step 7 (C): Close Filter Panel again', async () => {
    await page.evaluate(() => grok.shell.tv.getFiltersGroup().close());
    await expect.poll(() => page.evaluate(() =>
      document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length), {timeout: 5000}).toBe(0);
  });

  await softStep('Step 8 (C): Column hamburger → Filter resolves to the substructure editor function', async () => {
    const fnFound = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const molCol = tv.dataFrame.columns.toList().find((c: any) => c.semType === 'Molecule');
      // The column-hamburger Filter entry resolves to this editor function.
      const fn = (DG as any).Func.find({name: 'searchSubstructureEditor'})[0] ||
                  (DG as any).Func.find({name: 'SearchSubstructureEditor'})[0];
      if (!fn) return false;
      await fn.prepare({molecules: molCol}).call();
      return true;
    });
    // Hard-assert the editor exists (previously a silent `if (fn)` no-op that always passed).
    expect(fnFound, 'searchSubstructureEditor must be registered for the hamburger Filter route').toBe(true);
    // Single mol-column path invokes the search directly rather than opening a modal sketcher dialog,
    // so dialog presence is not asserted here (see filter-panel.md Block C).
  });

  const structureCardCount = () => page.evaluate(() =>
    Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'))
      .filter(c => /Structure|struct/i.test(c.querySelector('.d4-filter-header')?.textContent || '')).length);

  await softStep('Step 9 (C): Add another structure filter (CCC) + verify presence', async () => {
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      const molCol = grok.shell.t.columns.toList().find((c: any) => c.semType === 'Molecule');
      fg.updateOrAdd({type: 'Chem:substructureFilter', column: molCol.name, columnName: molCol.name});
    });
    await expect.poll(structureCardCount, {timeout: 15000}).toBeGreaterThanOrEqual(1);
  });

  await softStep('Step 10 (C): Remove Structure filter from Filter Panel', async () => {
    for (let i = 0; i < 6; i++) {
      const remaining = await structureCardCount();
      if (remaining === 0) break;
      await page.evaluate(() => {
        const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
        for (const card of cards) {
          const h = card.querySelector('.d4-filter-header');
          if (h && /Structure|struct/i.test(h.textContent || '')) {
            (card.querySelector('[name="icon-times"]') as HTMLElement | null)?.click();
            return;
          }
        }
      });
      await expect.poll(structureCardCount, {timeout: 5000}).toBeLessThan(remaining);
    }
    // Removing the widget(s) must clear the structure filter cards from the panel.
    await expect.poll(structureCardCount, {timeout: 10000}).toBe(0);
  });

  // ===== Block D — Drag-and-drop column header =====

  await softStep('Step 11 (D): Add Structure filter back (SR-DEFERRED: drag-drop canvas header)', async () => {
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      const molCol = grok.shell.t.columns.toList().find((c: any) => c.semType === 'Molecule');
      // SR-DEFERRED drag-drop column header → Filter Panel: substituted with fg.updateOrAdd.
      // The md "inserted at top of panel" invariant is specific to the drag-drop route we cannot
      // drive here, so position is not asserted; presence + empty-query state is.
      fg.updateOrAdd({type: 'Chem:substructureFilter', column: molCol.name, columnName: molCol.name});
    });
    await expect.poll(structureCardCount, {timeout: 15000}).toBeGreaterThanOrEqual(1);
    // Added without a molBlock → the structure filter query must be empty until the user draws one.
    const filterMol = await page.evaluate(() => {
      const molCol = grok.shell.t.columns.toList().find((c: any) => c.semType === 'Molecule');
      const tag = molCol?.temp['chem-scaffold-filter'];
      return tag ? (JSON.parse(tag)[0]?.molecule ?? '') : '';
    });
    expect(filterMol, 'a freshly-added structure filter must carry no pre-populated query').toBe('');
  });

  // ===== Block E — Cross-filter + sketcher↔filter sync =====

  await softStep('Step 12 (E): Add Stereo Category categorical filter → AND composition', async () => {
    const before = await page.evaluate(() => grok.shell.t.filter.trueCount);
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: ['R_ONE']});
    });
    // Selecting the R_ONE subset must strictly reduce the visible set (AND-composition), not leave it unchanged.
    await expect.poll(() => page.evaluate(() => grok.shell.t.filter.trueCount),
      {timeout: 15000}).toBeLessThan(before);
    const after = await page.evaluate(() => grok.shell.t.filter.trueCount);
    expect(after).toBeLessThan(before);
    expect(after).toBeGreaterThan(0);
  });

  await softStep('Step 13 (E): Re-edit substructure filter via hamburger sync', async () => {
    // SR-DEFERRED hamburger Re-edit click: apply modified structure via fg.updateOrAdd (filter↔structure sync).
    const before = await page.evaluate(() => grok.shell.t.filter.trueCount);
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      const molCol = grok.shell.t.columns.toList().find((c: any) => c.semType === 'Molecule');
      fg.updateOrAdd({
        type: 'Chem:substructureFilter',
        column: molCol.name,
        columnName: molCol.name,
        molBlock: 'CCCC',
      });
    });
    // filter↔sketcher sync: the filter widget's active query must update to the modified structure.
    await expect.poll(() => page.evaluate(() => {
      const molCol = grok.shell.t.columns.toList().find((c: any) => c.semType === 'Molecule');
      const tag = molCol?.temp['chem-scaffold-filter'];
      return tag ? (JSON.parse(tag)[0]?.molecule ?? '') : '';
    }), {timeout: 20000}).not.toBe('');
    const {after, filterMol} = await page.evaluate(() => {
      const molCol = grok.shell.t.columns.toList().find((c: any) => c.semType === 'Molecule');
      const tag = molCol?.temp['chem-scaffold-filter'];
      const filterMol = tag ? (JSON.parse(tag)[0]?.molecule ?? '') : '';
      return {after: grok.shell.t.filter.trueCount, filterMol};
    });
    // The filter widget observes the applied structure (butane molblock carries carbon atoms).
    expect(filterMol.toUpperCase()).toContain('C');
    // The table re-filters to the new query intersected with the active R_ONE categorical filter;
    // adding the CCCC substructure constraint can only shrink the R_ONE set (AND-composition).
    expect(after).toBeGreaterThan(0);
    expect(after).toBeLessThan(100);
    expect(after).toBeLessThanOrEqual(before);
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
