import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';

test('FilterPanel — Cloned Views', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, withFilterPanel: true});

  let originalFilteredCount: number;

  await softStep('Set up filters: Competition assay missing values, Stereo Category, Structure', async () => {
    const result = await page.evaluate(async () => {
      const fg = grok.shell.tv.getFiltersGroup();
      const df = grok.shell.tv.dataFrame;

      fg.updateOrAdd({type: DG.FILTER_TYPE.HISTOGRAM, column: 'Competition assay', filterOutMissingValues: true});
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: ['S_ACHIR']});

      await new Promise(r => setTimeout(r, 2000));
      return {filteredCount: df.filter.trueCount, totalRows: df.rowCount};
    });
    expect(result.filteredCount).toBeLessThan(result.totalRows);

    // Add the Structure substructure filter via the canonical FilterPanel API
    // (same fg.updateOrAdd path the histogram/categorical filters use above).
    // The `.sketch-link` SMILES-dialog UI commits the molecule to the sketcher
    // but does NOT engage substructure row-filtering (verified on dev), so
    // updateOrAdd with a SUBSTRUCTURE molBlock is the reliable way to drive it.
    // 'c1ccncc1' (pyridine) is a discriminating substructure on spgi-100 — it
    // narrows the already-filtered subset (the original benzene 'c1ccccc1' is
    // present in every survivor, so it could never reduce the count).
    originalFilteredCount = await page.evaluate(async (baseline) => {
      const fg = grok.shell.tv.getFiltersGroup();
      const df = grok.shell.tv.dataFrame;
      fg.updateOrAdd({type: DG.FILTER_TYPE.SUBSTRUCTURE, column: 'Structure',
        columnName: 'Structure', molBlock: 'c1ccncc1'});
      for (let i = 0; i < 25; i++) {
        await new Promise(r => setTimeout(r, 1000));
        if (df.filter.trueCount !== baseline) break;
      }
      return df.filter.trueCount;
    }, result.filteredCount);
    // The Structure filter's sketcher canvas should now be rendered.
    expect(await page.locator('[name="viewer-Filters"] .chem-external-sketcher-canvas').count())
      .toBeGreaterThan(0);
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
    // Clear via the sketcher's clear button (UI gesture on the panel), then
    // re-arm the substructure filter with cyclohexane through the canonical
    // updateOrAdd path (the SMILES-dialog UI does not engage row-filtering).
    const newCount = await page.evaluate(async (orig) => {
      const clearBtn = document.querySelector('[name="viewer-Filters"] .chem-clear-sketcher-button') as HTMLElement;
      clearBtn?.click();
      await new Promise(r => setTimeout(r, 2000));
      const fg = grok.shell.tv.getFiltersGroup();
      const df = grok.shell.tv.dataFrame;
      fg.updateOrAdd({type: DG.FILTER_TYPE.SUBSTRUCTURE, column: 'Structure',
        columnName: 'Structure', molBlock: DG.WHITE_MOLBLOCK});
      await new Promise(r => setTimeout(r, 1500));
      fg.updateOrAdd({type: DG.FILTER_TYPE.SUBSTRUCTURE, column: 'Structure',
        columnName: 'Structure', molBlock: 'C1CCCCC1'});
      for (let i = 0; i < 25; i++) {
        await new Promise(r => setTimeout(r, 1000));
        if (df.filter.trueCount !== orig) break;
      }
      return df.filter.trueCount;
    }, originalFilteredCount);
    expect(newCount).not.toBe(originalFilteredCount);
  });

  await softStep('Remove Structure filter — filtered state should not change', async () => {
    const beforeRemove = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);

    const removed = await page.evaluate(async () => {
      const filterPanel = document.querySelector('[name="viewer-Filters"]')!;
      const headers = filterPanel.querySelectorAll('.d4-filter-header');
      let clicked = false;
      for (const header of headers) {
        if (header.textContent!.trim().startsWith('Structure')) {
          const card = header.closest('.d4-filter');
          const closeIcon = card?.querySelector('[name="icon-times"]') as HTMLElement;
          closeIcon?.click();
          clicked = true;
          break;
        }
      }
      await new Promise(r => setTimeout(r, 2500));
      // The Structure filter card is removed from this (cloned) view's group.
      const fg = grok.shell.tv.getFiltersGroup();
      const structFiltersLeft = (fg.filters as any[])
        .filter((f: any) => f.columnName === 'Structure').length;
      return {clicked, structFiltersLeft};
    });

    const afterRemove = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    // Per the scenario: removing the Structure filter from the cloned view does
    // not change the filtered state — the substructure constraint was applied to
    // the shared DataFrame and the cloned view's removal leaves that row state in
    // place (verified on dev). Assert the card was actually removed and the
    // filtered count is unchanged.
    expect(removed.clicked).toBe(true);
    expect(removed.structFiltersLeft).toBe(0);
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

  v.finishSpec();
});
