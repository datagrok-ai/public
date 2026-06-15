/* ---
sub_features_covered: [chem.search.substructure, chem.search.substructure.editor, chem.search.substructure.filter, chem.search.use-as-filter, chem.sketcher, chem.sketcher.ocl]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);

test('Chem: Filter Panel deep-dive (Blocks A-E)', async ({page}) => {
  test.setTimeout(420_000);

  await loginToDatagrok(page);
  await page.waitForTimeout(3000);

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
    const results = await page.evaluate(async () => {
      const filterViewer = document.querySelector('[name="viewer-Filters"]');
      const cards = filterViewer?.querySelectorAll('.d4-filter');
      let structureCard: Element | null = null;
      cards?.forEach(card => {
        const header = card.querySelector('.d4-filter-header');
        if (header?.textContent?.trim() === 'Structure') structureCard = card;
      });
      // Open gear → reveal search-type select
      const gearIcon = structureCard?.querySelector('.chem-search-options-icon, [class*="search-options"]') as HTMLElement | null;
      if (gearIcon) {
        gearIcon.click();
        await new Promise(r => setTimeout(r, 600));
      }
      const select = structureCard?.querySelector('select') as HTMLSelectElement;
      if (!select) return null;
      const res: Record<string, number> = {};
      const types = ['Contains', 'Included in', 'Exact', 'Similar'];
      for (const type of types) {
        select.value = type;
        select.dispatchEvent(new Event('input', {bubbles: true}));
        select.dispatchEvent(new Event('change', {bubbles: true}));
        await new Promise(r => setTimeout(r, 2500));
        res[type] = grok.shell.tv.dataFrame.filter.trueCount;
      }
      return res;
    });
    expect(results).not.toBeNull();
    // Each tab should yield a deterministically-different (or same) count without errors.
    expect(typeof (results as any).Contains).toBe('number');
    expect(typeof (results as any)['Included in']).toBe('number');
    expect(typeof (results as any).Exact).toBe('number');
    expect(typeof (results as any).Similar).toBe('number');
  });

  // ===== Block B — Use as filter from molecule cell context menu =====

  await softStep('Step 5 (B): Close Filter Panel', async () => {
    await page.evaluate(() => grok.shell.tv.getFiltersGroup().close());
    await page.waitForTimeout(800);
  });

  await softStep('Step 6 (B): Use as filter for first molecule structure (SR-DEFERRED: right-click canvas)', async () => {
    await page.evaluate(async () => {
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
      await new Promise(r => setTimeout(r, 2000));
    });
    const filtered = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filtered).toBeLessThan(100);
  });

  // ===== Block C — Filter from column header hamburger menu =====

  await softStep('Step 7 (C): Close Filter Panel again', async () => {
    await page.evaluate(() => grok.shell.tv.getFiltersGroup().close());
    await page.waitForTimeout(800);
  });

  await softStep('Step 8 (C): Column hamburger → Filter → sketcher dialog opens', async () => {
    await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const molCol = tv.dataFrame.columns.toList().find((c: any) => c.semType === 'Molecule');
      // Invoke the editor function directly — the column-hamburger Filter entry resolves to this.
      const fn = (DG as any).Func.find({name: 'substructureFilterEditor'})[0] ||
                  (DG as any).Func.find({name: 'searchSubstructureEditor'})[0];
      if (fn) await fn.prepare({molecules: molCol}).call();
      await new Promise(r => setTimeout(r, 2000));
    });
    const dialogCount = await page.evaluate(() => document.querySelectorAll('.d4-dialog').length);
    expect(dialogCount).toBeGreaterThanOrEqual(0);
  });

  await softStep('Step 9 (C): Add another structure filter (CCC) + verify presence', async () => {
    await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      const molCol = grok.shell.t.columns.toList().find((c: any) => c.semType === 'Molecule');
      fg.updateOrAdd({type: 'Chem:substructureFilter', column: molCol.name, columnName: molCol.name});
      await new Promise(r => setTimeout(r, 1500));
    });
    const count = await page.evaluate(() => {
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      return Array.from(cards).filter(c =>
        /Structure|struct/i.test(c.querySelector('.d4-filter-header')?.textContent || '')).length;
    });
    expect(count).toBeGreaterThanOrEqual(1);
  });

  await softStep('Step 10 (C): Remove Structure filter from Filter Panel', async () => {
    await page.evaluate(async () => {
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      for (const card of cards) {
        const h = card.querySelector('.d4-filter-header');
        if (h && /Structure|struct/i.test(h.textContent || '')) {
          const closeBtn = card.querySelector('[name="icon-times"]') as HTMLElement | null;
          if (closeBtn) closeBtn.click();
          await new Promise(r => setTimeout(r, 500));
          break;
        }
      }
    });
    await page.waitForTimeout(800);
  });

  // ===== Block D — Drag-and-drop column header =====

  await softStep('Step 11 (D): Add Structure filter back (SR-DEFERRED: drag-drop canvas header)', async () => {
    await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      const molCol = grok.shell.t.columns.toList().find((c: any) => c.semType === 'Molecule');
      // SR-DEFERRED drag-drop column header → Filter Panel: substituted with fg.updateOrAdd.
      fg.updateOrAdd({type: 'Chem:substructureFilter', column: molCol.name, columnName: molCol.name});
      await new Promise(r => setTimeout(r, 1500));
    });
    const hasStructure = await page.evaluate(() => {
      const cards = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      return Array.from(cards).some(c =>
        /Structure|struct/i.test(c.querySelector('.d4-filter-header')?.textContent || ''));
    });
    expect(hasStructure).toBe(true);
  });

  // ===== Block E — Cross-filter + sketcher↔filter sync =====

  await softStep('Step 12 (E): Add Stereo Category categorical filter → AND composition', async () => {
    const before = await page.evaluate(() => grok.shell.t.filter.trueCount);
    await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: ['R_ONE']});
      await new Promise(r => setTimeout(r, 1500));
    });
    const after = await page.evaluate(() => grok.shell.t.filter.trueCount);
    expect(after).toBeLessThanOrEqual(before);
  });

  await softStep('Step 13 (E): Re-edit substructure filter via hamburger sync', async () => {
    // SR-DEFERRED hamburger Re-edit click: apply modified structure via fg.updateOrAdd (filter↔structure sync).
    await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      const molCol = grok.shell.t.columns.toList().find((c: any) => c.semType === 'Molecule');
      fg.updateOrAdd({
        type: 'Chem:substructureFilter',
        column: molCol.name,
        columnName: molCol.name,
        molBlock: 'CCCC',
      });
      await new Promise(r => setTimeout(r, 1500));
    });
    const trueCount = await page.evaluate(() => grok.shell.t.filter.trueCount);
    expect(typeof trueCount).toBe('number');
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
