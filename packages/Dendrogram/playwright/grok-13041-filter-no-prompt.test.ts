/* ---
sub_features_covered: [dendrogram.clustering.inject-tree-for-grid, dendrogram.event.selection-changed]
--- */
// GROK-13041: a filter alone MUST NOT raise the dendrogram sort-revert overlay (regression
// guard via filterChangeCounter in inject-tree-for-grid2.ts:313-352); a real sort with no
// preceding filter MUST raise it (positive contrast). The prompt surface is .dendrogram-overlay
// (not a .d4-dialog) with a [name="button-Revert-sort"]. Filter is driven via df.rows.filter and
// sort via grid.sort() (grid column headers are canvas, not DOM); both fire the same df events as
// the UI. Resource-load 404s on every neighbor mount are non-fatal noise (see isFatalConsoleError).
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

// True for fatal app errors; false for non-fatal noise (404 resource loads, ResizeObserver loop).
function isFatalConsoleError(text: string): boolean {
  if (/Failed to load resource[\s\S]*404/i.test(text)) return false;
  if (/ResizeObserver loop/i.test(text)) return false;
  return true;
}

async function openHierarchicalClusteringDialog(page: Page): Promise<void> {
  await page.evaluate(async () => {
    const chem = document.querySelector('[name="div-Chem"]') as HTMLElement | null;
    if (!chem) throw new Error('Top-menu Chem entry not found');
    chem.dispatchEvent(new MouseEvent('click', {bubbles: true}));
    const findLabel = async (pred: (t: string) => boolean): Promise<HTMLElement | undefined> => {
      for (let i = 0; i < 60; i++) {
        const el = Array.from(document.querySelectorAll('.d4-menu-item-label'))
          .find(m => pred((m.textContent || '').trim())) as HTMLElement | undefined;
        if (el) return el;
        await new Promise(r => setTimeout(r, 100));
      }
      return undefined;
    };
    const analyze = await findLabel(t => t === 'Analyze');
    if (!analyze) throw new Error('"Analyze" sub-menu item not found');
    (analyze.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
    const hc = await findLabel(t => /Hierarchical\s+Clustering/i.test(t));
    if (!hc) throw new Error('"Hierarchical Clustering..." sub-menu item not found');
    (hc.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
  });
  await page.locator('[name="dialog-Hierarchical-Clustering"]').waitFor({timeout: 15_000});
}

async function clickOkAndWaitForNeighbor(page: Page): Promise<number> {
  await page.locator('[name="dialog-Hierarchical-Clustering"] [name="button-OK"]').click();
  // Neighbor mounts when .dendrogram-assign-clusters-bttn (magic wand) appears.
  const foundAtMs: number = await page.evaluate(async () => {
    const start = Date.now();
    for (let i = 0; i < 120; i++) {
      if (document.querySelector('.dendrogram-assign-clusters-bttn'))
        return Date.now() - start;
      await new Promise(r => setTimeout(r, 500));
    }
    return -1;
  });
  return foundAtMs;
}

// Snapshot UI state relevant to the bug invariant: prompt dialog / sort-revert overlay / neighbor.
async function readPromptState(page: Page): Promise<{
  dialogCount: number;
  overlayCount: number;
  revertOverlayTextPresent: boolean;
  neighborMounted: boolean;
  tempKeySet: boolean;
  filterTrueCount: number;
  rowCount: number;
}> {
  return await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    const overlayEls = document.querySelectorAll('.dendrogram-overlay');
    const overlayText = overlayEls.length > 0 ? (overlayEls[0].textContent || '') : '';
    return {
      dialogCount: document.querySelectorAll('.d4-dialog').length,
      overlayCount: overlayEls.length,
      revertOverlayTextPresent: overlayText.includes('Revert columns sort order to see Dendrogram Tree'),
      neighborMounted: !!document.querySelector('.dendrogram-assign-clusters-bttn'),
      tempKeySet: !!grok.shell.tv.grid.temp['__dendrogram_neighbor_temp__'],
      filterTrueCount: df.filter.trueCount,
      rowCount: df.rowCount,
    };
  });
}

test('Dendrogram / GROK-13041 — filter does NOT trigger remove/revert prompt; sort DOES (positive contrast)', async ({page}) => {
  test.setTimeout(240_000);

  await loginToDatagrok(page);

  // Setup — open mol1K and wait for the Molecule semType + Chem package.
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    await new Promise(r => setTimeout(r, 1000));
    const df = await grok.dapi.files.readCsv('System:AppData/Chem/mol1K.csv');
    grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 4000);
    });
    // Chem dataset: wait for Grid canvas + extra settle for Chem package warmup.
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise(r => setTimeout(r, 200));
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  // Gate on the real Chem-ready signal (top-menu Chem attaches once semType is detected).
  await waitForChemMenu(page);
  // The Chem menu attaches before the package finishes warming up its feature
  // auto-detection; without this settle the Hierarchical Clustering dialog opens with
  // Features(0) and clustering produces no dendrogram (HEAD kept a 5s settle here).
  await page.waitForTimeout(5_000);

  // Console-error listener for the full test body (Step 5 asserts no console error on filter).
  const consoleErrors: string[] = [];
  const consoleListener = (msg: any) => { if (msg.type() === 'error') consoleErrors.push(msg.text()); };
  page.on('console', consoleListener);

  try {
    await softStep('1. Open mol1K and verify molecule + pIC50_HIV_Integrase columns', async () => {
      const info = await page.evaluate(() => {
        const df = grok.shell.tv?.dataFrame;
        const molCol = df?.col('molecule');
        const pIC50 = df?.col('pIC50_HIV_Integrase');
        return {
          rows: df?.rowCount,
          molSemType: molCol?.semType,
          pIC50Type: pIC50?.type,
          pIC50Min: pIC50?.min,
          pIC50Max: pIC50?.max,
        };
      });
      expect(info.rows, 'mol1K row count').toBe(1000);
      expect(info.molSemType, 'molecule semType').toBe('Molecule');
      expect(info.pIC50Type, 'pIC50_HIV_Integrase type').toBe('double');
      // Step 5 filters < 6 — verify the range is non-degenerate so the filter excludes a subset.
      expect(info.pIC50Min!, 'pIC50 min < 6 (rows survive filter)').toBeLessThan(6);
      expect(info.pIC50Max!, 'pIC50 max > 6 (rows excluded by filter)').toBeGreaterThan(6);
    });

    await softStep('2. Run Chem | Analyze | Hierarchical Clustering... — dialog opens with euclidean+ward+molecule', async () => {
      await openHierarchicalClusteringDialog(page);
      const defaults = await page.evaluate(() => ({
        dialogPresent: !!document.querySelector('[name="dialog-Hierarchical-Clustering"]'),
        distance: (document.querySelector('[name="input-Distance"]') as HTMLSelectElement)?.value,
        linkage: (document.querySelector('[name="input-Linkage"]') as HTMLSelectElement)?.value,
        features: document.querySelector('[name="input-host-Features"]')?.textContent?.trim(),
      }));
      expect(defaults.dialogPresent, 'Hierarchical Clustering dialog opened').toBe(true);
      expect(defaults.distance, 'Distance default').toBe('euclidean');
      expect(defaults.linkage, 'Linkage default').toBe('ward');
      expect(defaults.features, 'Features defaults to molecule').toContain('molecule');
    });

    await softStep('3. Click OK → dendrogram neighbor injected; magic wand + close icon + temp key all set', async () => {
      const foundAtMs = await clickOkAndWaitForNeighbor(page);
      expect(foundAtMs, 'Magic-wand mount time (ms; -1 = timeout)').toBeGreaterThan(0);
      const state = await readPromptState(page);
      expect(state.neighborMounted, 'magic-wand icon present (neighbor mounted)').toBe(true);
      expect(state.tempKeySet, '__dendrogram_neighbor_temp__ set on grid').toBe(true);
      // Pre-filter / pre-sort state: no prompt dialog, no overlay yet.
      expect(state.dialogCount, 'no prompt dialog before filter/sort').toBe(0);
      expect(state.overlayCount, 'no .dendrogram-overlay before filter/sort').toBe(0);
      expect(state.filterTrueCount, 'filter.trueCount == rowCount (1000) initially').toBe(state.rowCount);
    });

    // Scenario 1 — GROK-13041: filter MUST NOT trigger the remove/revert prompt.
    await softStep('4. Open Filter Panel for the active table view', async () => {
      await page.evaluate(() => grok.shell.tv.getFiltersGroup());
      await page.locator('[name="viewer-Filters"]').waitFor({timeout: 10_000});
      const state = await page.evaluate(() => ({
        filterPanelMounted: !!document.querySelector('[name="viewer-Filters"]'),
        filtersListed: document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length,
        pIC50FilterPresent: Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter'))
          .some(f => (f.textContent || '').includes('pIC50_HIV_Integrase')),
        promptState: {
          dialogCount: document.querySelectorAll('.d4-dialog').length,
          overlayCount: document.querySelectorAll('.dendrogram-overlay').length,
        },
      }));
      expect(state.filterPanelMounted, 'Filter Panel mounted').toBe(true);
      expect(state.filtersListed, 'at least one filter card listed').toBeGreaterThan(0);
      expect(state.pIC50FilterPresent, 'pIC50_HIV_Integrase filter card present (the filter step 5 drives)').toBe(true);
      // Opening the panel itself MUST not raise the prompt.
      expect(state.promptState.dialogCount, 'no prompt dialog from opening Filter Panel').toBe(0);
      expect(state.promptState.overlayCount, 'no .dendrogram-overlay from opening Filter Panel').toBe(0);
    });

    await softStep('5. Apply filter pIC50_HIV_Integrase < 6 → NO remove/revert prompt; neighbor intact; no console error', async () => {
      // Snapshot row order before the filter — it must be unchanged after filtering.
      const beforeOrder = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const order = [];
        for (let i = 0; i < 10; i++) order.push(df.get('prID', i));
        return order;
      });
      // df.rows.filter fires df.onFilterChanged — the bug-relevant entry (inject-tree-for-grid2.ts:316).
      await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        df.rows.filter((row: any) => row.pIC50_HIV_Integrase < 6);
      });
      // Deterministically wait for the filter to have taken effect; alignGridWithTree is
      // dispatched synchronously off the same onFilterChanged event, so a short yield suffices.
      await page.waitForFunction(() =>
        grok.shell.tv.dataFrame.filter.trueCount < grok.shell.tv.dataFrame.rowCount, null, {timeout: 15_000});
      await page.waitForTimeout(300);

      const state = await readPromptState(page);
      const afterOrder = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const order = [];
        for (let i = 0; i < 10; i++) order.push(df.get('prID', i));
        return order;
      });

      // GROK-13041: the filter trigger MUST NOT surface either prompt-shaped UI.
      expect(state.dialogCount, 'GROK-13041 invariant: no .d4-dialog opens on filter').toBe(0);
      expect(state.overlayCount, 'GROK-13041 invariant: no .dendrogram-overlay appears on filter').toBe(0);
      expect(state.revertOverlayTextPresent, 'GROK-13041 invariant: no "Revert columns sort order" text on filter').toBe(false);
      // The neighbor stays attached after the filter (no auto-detach).
      expect(state.neighborMounted, 'magic-wand still present after filter').toBe(true);
      expect(state.tempKeySet, '__dendrogram_neighbor_temp__ still set after filter').toBe(true);
      // The filter actually restricted a non-empty subset of rows.
      expect(state.filterTrueCount, 'filter.trueCount reduced').toBeLessThan(state.rowCount);
      expect(state.filterTrueCount, 'filter.trueCount > 0 (non-empty subset)').toBeGreaterThan(0);
      // Row order unchanged — filtering does not reorder the underlying df.
      expect(afterOrder, 'underlying row order unchanged by filter').toEqual(beforeOrder);
      // No fatal console errors (non-fatal 404 / ResizeObserver noise filtered).
      const fatalErrors = consoleErrors.filter(isFatalConsoleError);
      expect(fatalErrors, 'no fatal console errors during filter step').toEqual([]);
    });

    await softStep('6. Clear the filter → neighbor still attached; no late-arriving prompt', async () => {
      await page.evaluate(() => grok.shell.tv.dataFrame.filter.setAll(true));
      await page.waitForFunction(() =>
        grok.shell.tv.dataFrame.filter.trueCount === grok.shell.tv.dataFrame.rowCount, null, {timeout: 15_000});
      const state = await readPromptState(page);
      expect(state.filterTrueCount, 'filter cleared (trueCount == rowCount)').toBe(state.rowCount);
      expect(state.dialogCount, 'no .d4-dialog after clearing filter').toBe(0);
      expect(state.overlayCount, 'no .dendrogram-overlay after clearing filter').toBe(0);
      expect(state.neighborMounted, 'neighbor still attached after clearing filter').toBe(true);
      expect(state.tempKeySet, '__dendrogram_neighbor_temp__ still set after clearing filter').toBe(true);
    });

    // Scenario 2 — Sort DOES trigger the overlay (positive contrast); prompt surface is the overlay.
    await softStep('7. Sort pIC50_HIV_Integrase ascending → .dendrogram-overlay raised with "Revert columns sort order"', async () => {
      // grid.sort is the canvas-fallback path — grid column headers are canvas, not DOM.
      await page.evaluate(() => {
        grok.shell.tv.grid.sort(['pIC50_HIV_Integrase'], [true]);
      });
      await page.locator('.dendrogram-overlay').waitFor({state: 'visible', timeout: 15_000});

      const state = await readPromptState(page);
      const overlayDetail = await page.evaluate(() => {
        const o = document.querySelector('.dendrogram-overlay');
        return {
          revertBtnPresent: !!document.querySelector('[name="button-Revert-sort"]'),
          overlayClassList: o ? Array.from(o.classList) : [],
          overlayParentIsNeighborRoot: !!o?.parentElement?.querySelector('.dendrogram-assign-clusters-bttn'),
        };
      });

      // The overlay IS the prompt surface — assert overlay present + text present.
      expect(state.overlayCount, 'positive contrast: .dendrogram-overlay raised on sort').toBe(1);
      expect(state.revertOverlayTextPresent, 'positive contrast: overlay carries "Revert columns sort order to see Dendrogram Tree"').toBe(true);
      expect(overlayDetail.revertBtnPresent, '[name="button-Revert-sort"] present inside overlay').toBe(true);
      expect(overlayDetail.overlayClassList, '.dendrogram-overlay class present').toContain('dendrogram-overlay');
      // Overlay is attached to the neighbor root, not a stray div.
      expect(overlayDetail.overlayParentIsNeighborRoot, 'overlay is attached to the neighbor root').toBe(true);
      // No dialog on the sort path — surface any future regression that swaps the prompt surface.
      expect(state.dialogCount, 'sort path: no .d4-dialog opens (prompt is overlay, not dialog)').toBe(0);
      // The neighbor stays beneath the overlay (the overlay is an extra div, not a removal).
      expect(state.neighborMounted, 'neighbor magic-wand still in DOM beneath overlay').toBe(true);
    });

    await softStep('8. Click [name="button-Revert-sort"] in overlay → overlay dismissed; neighbor reattached cleanly', async () => {
      await page.locator('[name="button-Revert-sort"]').click();
      await page.locator('.dendrogram-overlay').waitFor({state: 'detached', timeout: 10_000});
      const state = await readPromptState(page);
      expect(state.overlayCount, 'overlay dismissed after Revert sort').toBe(0);
      expect(state.neighborMounted, 'neighbor magic-wand still present after revert').toBe(true);
      expect(state.dialogCount, 'no stray dialog after revert').toBe(0);
    });

    // Cleanup
    await page.evaluate(() => grok.shell.closeAll());
  } finally {
    page.off('console', consoleListener);
  }

  finishSpec();
});
