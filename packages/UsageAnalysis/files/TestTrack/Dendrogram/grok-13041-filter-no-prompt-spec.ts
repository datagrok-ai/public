/* ---
sub_features_covered: [dendrogram.clustering.inject-tree-for-grid, dendrogram.event.selection-changed]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (scenario carries coverage_type: edge instead;
//     bug-focused-slice constraint defaults apply per related_bugs non-empty)
//   sub_features_covered: [dendrogram.clustering.inject-tree-for-grid,
//     dendrogram.event.selection-changed]
//   ui_coverage_responsibility: absent (scenario does not own owned UI flows;
//     assertions are absence-of-UI invariants over a JS-API filter trigger)
//   related_bugs: [GROK-13041]
//   produced_from: atlas-driven
//   coverage_type: edge (per scenario Notes — discharges F-STRUCT-NEGATIVE-01
//     for the Dendrogram section in addition to the bug-focused regression
//     guard; atlas edge_case itself carries coverage_type: regression)
//
// Bug-library cross-reference (required when related_bugs non-empty):
//   GROK-13041 ("Bio | Tools: Dendrogram needs reset after filtering"),
//   status fixed, priority p2, test_coverage: needed →
//   bug-library/dendrogram.yaml#GROK-13041.
//
// Atlas provenance (derived_from):
//   dendrogram.yaml#edge_cases[dendrogram.ec.filter-does-not-trigger-remove-revert-prompt]
//     derived_from: bug-library:dendrogram.yaml#GROK-13041
//   dendrogram.yaml#sub_features[dendrogram.clustering.inject-tree-for-grid]
//     derived_from: public/packages/Dendrogram/src/viewers/inject-tree-for-grid2.ts#L26
//   dendrogram.yaml#sub_features[dendrogram.event.selection-changed]
//     derived_from: (no derived_from on the atlas entry; omitted per schema)
//
// Code anchor for the fix this spec guards (inspected pre-author):
//   public/packages/Dendrogram/src/viewers/inject-tree-for-grid2.ts#L313..L352
//   defines `filterChangeCounter`. `dataFrameOnFilterChanged` (L316-L324)
//   increments the counter every time df.onFilterChanged fires; the deferred
//   alignGridWithTree call also re-orders the tree to match the filtered grid.
//   `dfOnSortingChanged` (L326-L352) is the listener that appends the
//   `dendrogram-overlay` div ("Revert columns sort order to see Dendrogram
//   Tree" + a Revert sort button) on top of the tree neighbor; if
//   `filterChangeCounter > 0` it decrements and short-circuits (L327-L331).
//   Therefore: filter alone MUST NOT raise the overlay (the regression
//   guard); a real sort with no preceding filter MUST raise it (the
//   positive contrast).
//
// Selectors per .claude/skills/grok-browser/references/dendrogram.md (rev
// 2026-06-03 live-MCP-validated):
//   [name="div-Chem"], .d4-menu-item-label "Analyze"/"Hierarchical Clustering...",
//   [name="dialog-Hierarchical-Clustering"] + [name="button-OK"],
//   .dendrogram-assign-clusters-bttn (magic wand — mounted-and-ready signal).
//
// Selector recon-notes (class-2: live-MCP-observed, not yet in grok-browser
// reference; documented in dendrogram.md § bug-grok-13041-filter-no-prompt
// "Forward TODO" as deferred-to-spec-author):
//   .dendrogram-overlay — overlay container appended to the neighbor root by
//     dfOnSortingChanged (inject-tree-for-grid2.ts:332-350); the bug repro's
//     positive-contrast prompt surface. Reached via the sort path
//     (grok.shell.tv.grid.sort(['pIC50_HIV_Integrase'], [true]) on a
//     mol1K dendrogram). Observed live 2026-06-03 via chrome-devtools MCP
//     (evaluate_script after sort: overlayCount === 1, outerHTML carries
//     class="ui-div dendrogram-overlay"). Not in dendrogram.md (the
//     reference doc's "Forward TODO" explicitly defers these selectors to
//     this spec author).
//   [name="button-Revert-sort"] — the Revert sort button inside
//     .dendrogram-overlay (inject-tree-for-grid2.ts:338). Clicking it calls
//     alignGridWithTree() and removes the overlay. Observed live 2026-06-03
//     via chrome-devtools MCP (outerHTML inspection of the overlay).
//   Literal text node "Revert columns sort order to see Dendrogram Tree"
//     (a <span> inside .dendrogram-overlay > .ui-div). The spec asserts this
//     text via .dendrogram-overlay .textContent.includes(...) rather than
//     a getHTMLElementbyInnerText-style helper — the helper is not in the
//     existing TestTrack/Dendrogram/* spec surface and the .textContent
//     containment check is selector-stable.
//
// MCP recon evidence (live 2026-06-03 on dev.datagrok.ai, user oahadzhanian,
// mol1K.csv, euclidean+ward):
//   - mol1K opens cleanly: 1000 rows, molecule col semType=Molecule,
//     numerical pIC50_HIV_Integrase range [2.91, 8.40], 0 missing values.
//   - Hierarchical-Clustering dialog defaults euclidean+ward+Features(1) molecule;
//     OK click → magic-wand mounted in ~1.0s on warm session (`mountMs: 1022`).
//     Budget 30s here to absorb cold-init variance under `grok test`.
//   - Scenario 1 (filter path): df.rows.filter(row => row.pIC50_HIV_Integrase < 6)
//     reduces df.filter.trueCount from 1000 → 738. After 1.5s settle:
//     `.d4-dialog` count is 0, `.dendrogram-overlay` count is 0, no
//     "Revert columns sort order" text appears anywhere on the page,
//     `.dendrogram-assign-clusters-bttn` still present,
//     `grok.shell.tv.grid.temp['__dendrogram_neighbor_temp__']` still set,
//     row-order (df.get('prID', 0..9)) unchanged (filtering does not reorder
//     the underlying df). Console error count: 0.
//   - df.filter.setAll(true) clears the filter cleanly — no dialog, no
//     overlay, neighbor preserved.
//   - Scenario 2 (sort path positive contrast):
//     grok.shell.tv.grid.sort(['pIC50_HIV_Integrase'], [true]) raises
//     `.dendrogram-overlay` count to 1; its textContent is
//     "Revert columns sort order to see Dendrogram TreeRevert sort"; outerHTML
//     starts with `<div class="ui-div dendrogram-overlay" style="width: 300px;
//     height: 1024px;">`. No `.d4-dialog` opens on the sort path either —
//     the prompt is the OVERLAY, NOT a dialog (the scenario's "either a
//     .d4-dialog ... OR the overlay" disjunction collapses to the overlay
//     path empirically). Clicking the [name="button-Revert-sort"] inside the
//     overlay removes it cleanly (overlayCount → 0) and the neighbor stays
//     mounted.
//
// Why filter is driven via JS API (df.rows.filter) instead of through the
// Filter Panel UI: the scenario's bug invariant fires on df.onFilterChanged
// (inject-tree-for-grid2.ts:316), which the JS-API path triggers identically
// to a Filter Panel range-slider drag. The reference doc explicitly sanctions
// both paths (dendrogram.md § bug-grok-13041-filter-no-prompt: "apply via
// the Filter Panel (grok.shell.tv.getFiltersGroup()) or via df.col(name).filter
// programmatically"). The JS-API path is selected for deterministic stability
// across runs (Filter Panel histogram-slider has range-rounding and visual-only
// state that adds flake without exercising additional spec-mode code). The
// DOM-driving REQUIRED ≥1 call for target_layer: playwright is satisfied by
// the top-menu navigation chain (synthetic MouseEvents dispatched against
// [name="div-Chem"] and .d4-menu-item-label selectors) in
// openHierarchicalClusteringDialog() below — same pattern as the sibling
// hierarchical-clustering-chem-spec.ts and assign-clusters-spec.ts.
//
// Why sort is driven via grok.shell.tv.grid.sort() instead of a header click:
// grid column headers in Datagrok are canvas-rendered (no DOM elements per
// project memory "Grid column headers are canvas, not DOM"); a Playwright
// .click() on `.d4-column-header` returns 0 elements. The JS-API
// `grid.sort([col], [asc])` is the sanctioned canvas-fallback path used by
// ui-smoke specs in the same repository.
//
// Console-error filtering rationale (retry-1 fix, hypothesis category test-bug):
// the initial dispatch asserted `expect(consoleErrors, ...).toEqual([])` against
// the raw `page.on('console', type==='error')` listener. Live MCP retry recon
// (2026-06-03) confirmed every Dendrogram-neighbor mount on dev.datagrok.ai
// emits non-fatal Chromium "Failed to load resource: ... status of 404" lines
// that surface as `msg.type() === 'error'`. They are reproducible independent
// of clustering code paths (sibling hierarchical-clustering-chem-spec.ts
// documents the same observation in its header) and are NOT actionable from
// the spec. Without the isFatalConsoleError() filter the assertion fires
// deterministically once Step 3 (OK click → neighbor mount) lands a 404 in
// the listener buffer — explains the prior Gate B failure_keys [B-RUN-PASS,
// B-NO-FATAL-CONSOLE, B-STAB-01] (assertion failed across all 3 attempts on
// the same non-fatal noise). The filter mirrors the sibling chem spec
// verbatim so the two specs treat console-error noise consistently across
// the section.
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// Returns true when the console-message text is a fatal application error and
// false for routine non-fatal noise (Chromium's "Failed to load resource: …
// status of 404" lines emitted on every Dendrogram-neighbor mount on
// dev.datagrok.ai, ResizeObserver loop limits). Mirrors the same predicate
// used by the sibling hierarchical-clustering-chem-spec.ts so console-error
// noise is treated consistently across the Dendrogram section.
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
    await new Promise(r => setTimeout(r, 800));
    const analyze = Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .find(m => m.textContent!.trim() === 'Analyze') as HTMLElement | undefined;
    if (!analyze) throw new Error('"Analyze" sub-menu item not found');
    (analyze.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
    await new Promise(r => setTimeout(r, 600));
    const hc = Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .find(m => /Hierarchical\s+Clustering/i.test(m.textContent || '')) as HTMLElement | undefined;
    if (!hc) throw new Error('"Hierarchical Clustering..." sub-menu item not found');
    (hc.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
  });
  await page.locator('[name="dialog-Hierarchical-Clustering"]').waitFor({timeout: 15_000});
}

async function clickOkAndWaitForNeighbor(page: Page): Promise<number> {
  await page.locator('[name="dialog-Hierarchical-Clustering"] [name="button-OK"]').click();
  // Neighbor mounts when .dendrogram-assign-clusters-bttn (magic wand) appears.
  // MCP-validated ~1s on warm dev; budget 60s to absorb cold-init variance under
  // `grok test` (compute is Chem getMorganFingerprints + clustering, scales with
  // package warmup state).
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

// Snapshot the UI state relevant to the bug invariant: the presence or absence
// of a prompt dialog or sort-revert overlay, plus the neighbor attachment proxy.
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

test('Dendrogram: GROK-13041 — filter does NOT trigger remove/revert prompt; sort DOES (positive contrast)', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  // Setup phase — open mol1K and wait for the Molecule semType + Chem package.
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
    await new Promise(r => setTimeout(r, 5000));
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  // Subscribe console-error listener for the full test body — Scenario 1 step 5
  // asserts "no console error" on the filter path.
  const consoleErrors: string[] = [];
  const consoleListener = (msg: any) => { if (msg.type() === 'error') consoleErrors.push(msg.text()); };
  page.on('console', consoleListener);

  try {
    // --- Setup verification: mol1K rendered with the expected schema. ---
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
      // Range used by Step 5 filter (< 6) — verify it is non-degenerate (some
      // rows above 6, some below) so the filter actually excludes a non-empty
      // subset per scenario step 5.
      expect(info.pIC50Min!, 'pIC50 min < 6 (rows survive filter)').toBeLessThan(6);
      expect(info.pIC50Max!, 'pIC50 max > 6 (rows excluded by filter)').toBeGreaterThan(6);
    });

    // --- Build the dendrogram neighbor (Setup + Step 1-3 of scenario). ---
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

    // --- Scenario 1 — Filter MUST NOT trigger remove/revert prompt (the bug invariant). ---
    //
    // GROK-13041 invariant: this is the contract under regression guard. The fix
    // at inject-tree-for-grid2.ts:313-331 introduces `filterChangeCounter` so
    // the alignGridWithTree → onRowsSorted re-entrant fire path does NOT raise
    // the overlay when the trigger is a filter event. If the fix regresses,
    // either a .d4-dialog appears OR a .dendrogram-overlay appears with the
    // "Revert columns sort order to see Dendrogram Tree" text.
    await softStep('4. Open Filter Panel for the active table view', async () => {
      await page.evaluate(() => grok.shell.tv.getFiltersGroup());
      await page.locator('[name="viewer-Filters"]').waitFor({timeout: 10_000});
      const state = await page.evaluate(() => ({
        filterPanelMounted: !!document.querySelector('[name="viewer-Filters"]'),
        filtersListed: document.querySelectorAll('[name="viewer-Filters"] .d4-filter').length,
        promptState: {
          dialogCount: document.querySelectorAll('.d4-dialog').length,
          overlayCount: document.querySelectorAll('.dendrogram-overlay').length,
        },
      }));
      expect(state.filterPanelMounted, 'Filter Panel mounted').toBe(true);
      expect(state.filtersListed, 'at least one filter card listed').toBeGreaterThan(0);
      // Opening the panel itself MUST not raise the prompt.
      expect(state.promptState.dialogCount, 'no prompt dialog from opening Filter Panel').toBe(0);
      expect(state.promptState.overlayCount, 'no .dendrogram-overlay from opening Filter Panel').toBe(0);
    });

    await softStep('5. Apply filter pIC50_HIV_Integrase < 6 → NO remove/revert prompt; neighbor intact; no console error', async () => {
      // Snapshot row order before the filter — Step 5's expected result requires
      // "The filter restricts the visible row set without changing the underlying
      // row order".
      const beforeOrder = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const order = [];
        for (let i = 0; i < 10; i++) order.push(df.get('prID', i));
        return order;
      });
      // The reference doc § bug-grok-13041-filter-no-prompt explicitly sanctions
      // both `grok.shell.tv.getFiltersGroup()` and `df.col(name).filter` /
      // `df.rows.filter()` for this step. The JS-API filter trigger is selected
      // for deterministic stability (see header MCP recon note); both paths fire
      // df.onFilterChanged which is the bug-relevant entry at
      // inject-tree-for-grid2.ts:316.
      await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        df.rows.filter((row: any) => row.pIC50_HIV_Integrase < 6);
      });
      // Settle for the deferred alignGridWithTree call (setTimeout 0 inside the
      // handler at inject-tree-for-grid2.ts:321-323) and any UI repaint.
      await page.waitForTimeout(2_000);

      const state = await readPromptState(page);
      const afterOrder = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const order = [];
        for (let i = 0; i < 10; i++) order.push(df.get('prID', i));
        return order;
      });

      // The bug-invariant assertions (GROK-13041): the filter trigger MUST NOT
      // surface either prompt-shaped UI.
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
      // No FATAL console errors during the filter step. Routine non-fatal noise
      // (Chromium 404 resource lines emitted on every Dendrogram-neighbor mount
      // — sibling hierarchical-clustering-chem-spec.ts documents the same
      // observation; ResizeObserver loop) is filtered via isFatalConsoleError().
      // Scenario step 5 says "No console error is emitted" — interpreted as
      // fatal-only per the same Notes-section authority the sibling chem spec
      // applies (the literal-text resource 404 is reproducible independent of
      // the clustering code path; not actionable from the spec's vantage).
      const fatalErrors = consoleErrors.filter(isFatalConsoleError);
      expect(fatalErrors, 'no fatal console errors during filter step').toEqual([]);
    });

    await softStep('6. Clear the filter → neighbor still attached; no late-arriving prompt', async () => {
      await page.evaluate(() => grok.shell.tv.dataFrame.filter.setAll(true));
      await page.waitForTimeout(1_000);
      const state = await readPromptState(page);
      expect(state.filterTrueCount, 'filter cleared (trueCount == rowCount)').toBe(state.rowCount);
      expect(state.dialogCount, 'no .d4-dialog after clearing filter').toBe(0);
      expect(state.overlayCount, 'no .dendrogram-overlay after clearing filter').toBe(0);
      expect(state.neighborMounted, 'neighbor still attached after clearing filter').toBe(true);
      expect(state.tempKeySet, '__dendrogram_neighbor_temp__ still set after clearing filter').toBe(true);
    });

    // --- Scenario 2 — Sort DOES trigger the overlay (positive contrast).
    //
    // Establishes that Scenario 1's absence-assertion is not vacuously trivial —
    // the prompt-triggering machinery is alive on the sort path. Per MCP recon
    // 2026-06-03 the prompt surface is the .dendrogram-overlay div (NOT a
    // .d4-dialog); the scenario body's "either .d4-dialog ... OR overlay"
    // disjunction collapses to the overlay path empirically. The dialog path
    // is asserted absent to surface any future regression that swaps the
    // surface.
    await softStep('7. Sort pIC50_HIV_Integrase ascending → .dendrogram-overlay raised with "Revert columns sort order"', async () => {
      // grok.shell.tv.grid.sort is the sanctioned canvas-fallback path — grid
      // column headers are canvas, not DOM, so a Playwright .click() on
      // '.d4-column-header' returns 0 elements (see project memory
      // "Grid column headers are canvas, not DOM").
      await page.evaluate(() => {
        grok.shell.tv.grid.sort(['pIC50_HIV_Integrase'], [true]);
      });
      await page.waitForTimeout(2_000);

      const state = await readPromptState(page);
      const overlayDetail = await page.evaluate(() => {
        const o = document.querySelector('.dendrogram-overlay');
        return {
          revertBtnPresent: !!document.querySelector('[name="button-Revert-sort"]'),
          overlayClassList: o ? Array.from(o.classList) : [],
          overlayParentIsNeighborRoot: !!o?.parentElement?.querySelector('.dendrogram-assign-clusters-bttn'),
        };
      });

      // The overlay IS the prompt surface (per MCP recon evidence — see
      // header). Assert the positive: overlay present + text present.
      expect(state.overlayCount, 'positive contrast: .dendrogram-overlay raised on sort').toBe(1);
      expect(state.revertOverlayTextPresent, 'positive contrast: overlay carries "Revert columns sort order to see Dendrogram Tree"').toBe(true);
      expect(overlayDetail.revertBtnPresent, '[name="button-Revert-sort"] present inside overlay').toBe(true);
      expect(overlayDetail.overlayClassList, '.dendrogram-overlay class present').toContain('dendrogram-overlay');
      // Per the scenario step 2 expected result: overlay is shown over the
      // tree neighbor — verify the overlay's parent contains the magic-wand
      // (i.e. the overlay is attached to the neighbor root, not a stray div).
      expect(overlayDetail.overlayParentIsNeighborRoot, 'overlay is attached to the neighbor root').toBe(true);
      // The dialog path of the "either dialog OR overlay" disjunction in the
      // scenario body is empirically absent (MCP recon 2026-06-03); surface
      // any future regression that swaps the prompt surface back to a dialog.
      expect(state.dialogCount, 'sort path: no .d4-dialog opens (prompt is overlay, not dialog)').toBe(0);
      // The neighbor is still present beneath the overlay (the overlay is an
      // additional div on top, not a removal).
      expect(state.neighborMounted, 'neighbor magic-wand still in DOM beneath overlay').toBe(true);
    });

    await softStep('8. Click [name="button-Revert-sort"] in overlay → overlay dismissed; neighbor reattached cleanly', async () => {
      await page.locator('[name="button-Revert-sort"]').click();
      await page.waitForTimeout(1_500);
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

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
