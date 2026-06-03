/* ---
sub_features_covered: [dendrogram.clustering.menu.chem, dendrogram.clustering.dialog, dendrogram.clustering.inject-tree-for-grid]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: ui-smoke
//   sub_features_covered: [dendrogram.clustering.menu.chem, dendrogram.clustering.dialog,
//     dendrogram.clustering.inject-tree-for-grid]
//   ui_coverage_responsibility: [hierarchical-clustering-dialog,
//     hierarchical-clustering-dialog-distance-dropdown,
//     hierarchical-clustering-dialog-linkage-dropdown,
//     hierarchical-clustering-dialog-features-column-input,
//     hierarchical-clustering-dialog-ok-button,
//     chem-analyze-hierarchical-clustering-menu-entry] (delegated_to: null)
//   related_bugs: []
//   produced_from: migrated
//
// Atlas provenance (derived_from):
//   dendrogram.yaml#sub_features[dendrogram.clustering.menu.chem].interactions[0]
//     derived_from: [SRC Dendrogram:hierarchicalClusteringMolecules
//                    public/packages/Dendrogram/src/package.g.ts#L106]
//   dendrogram.yaml#sub_features[dendrogram.clustering.dialog]
//     derived_from: public/packages/Dendrogram/src/utils/hierarchical-clustering.ts#L16
//   dendrogram.yaml#sub_features[dendrogram.clustering.inject-tree-for-grid]
//     derived_from: public/packages/Dendrogram/src/viewers/inject-tree-for-grid2.ts#L26
//
// Scope reductions (carried from scenario frontmatter source_text_fixes
// and discovered during MCP recon 2026-06-03):
//   SR-01: scenario Step 7 "non-monotonic centroid/median tree" visual assertion
//     is NOT exercised — the non-monotonic property is a canvas-pixel /
//     merge-height-sequence property and is not addressable via Playwright DOM.
//     The atlas matrix coverage of all 14 distance × linkage combinations lives
//     in the sibling apitest (hierarchical-clustering-chem-api.md / -api.ts).
//     The UI smoke here asserts only what the scenario actually claims at the
//     UI layer: OK click → dendrogram neighbor mounts + no fatal console error
//     surfaces. The scenario's source_text_fixes carry the matching slug
//     `downgrade-step-7-non-monotonic-visual-assertion-to-no-throw` and the
//     `unresolved_ambiguities:
//      step-7-non-monotonic-tree-not-assertable-via-playwright-canvas` entry.
//   SR-02: scenario Step 7 "Features = numeric columns
//     (pIC50_HIV_Integrase, Q)" substitution is NOT exercised via the
//     Hierarchical Clustering dialog's Features picker. MCP recon
//     (2026-06-03) established that clicking the [name="input-Features"]
//     editor opens a [name="dialog-Select-columns..."] popup whose
//     column-grid is `.d4-column-grid` — a canvas-rendered DG.ColumnGrid
//     with NO DOM checkboxes (verified: 0 input[type="checkbox"] elements
//     in the popup; only the .d4-search-input, the `[name="label-All"]` /
//     `[name="label-None"]` toggles and a type-filter SELECT are
//     DOM-addressable). Individual-column toggling is canvas-pixel only.
//     DG.Dialog.getOpenDialogs()[0].input("Features") also throws
//     "Input \"Features\" not found." — the dialog's inputs[] bag is empty
//     for this column-list widget.
//     `grok.functions.call('Dendrogram:hierarchicalClustering', {df,
//     colNameList, distance, linkage})` bypasses the dialog entirely
//     (atlas `dendrogram.clustering.api`) but is FORBIDDEN as a
//     substitution for the ui-smoke owned flow
//     `hierarchical-clustering-dialog-features-column-input`.
//     Canvas-click on column-grid rows at bounding-box coordinates was
//     NOT exhaustively refuted this session, so the manual-only split
//     trigger (per §"ui-affordance manual-only split") does not yet fire.
//     Step 7's numeric-features assertion is therefore deferred to the
//     sibling apitest (hierarchical-clustering-chem-api.md), which can
//     pass the column-list to the API directly. Step 7 here exercises
//     the median linkage on the DEFAULT (molecule) features (see SR-03
//     below for why centroid was swapped for median).
//   SR-03: scenario Step 7 "centroid (or median)" — the OR is load-bearing.
//     MCP recon (2026-06-03 retry) established that the centroid linkage
//     on the MOLECULE features path triggers a fatal platform crash:
//     `Error: memory access out of bounds` (WASM memory access fault),
//     surfaced through `hierarchical-clustering.ts:217` (the catch-and-rethrow
//     in `hierarchicalClustering()`). The dendrogram never mounts; the
//     magic-wand button never appears. This is a platform-level core-bug
//     in the centroid-linkage + Morgan-fingerprint compute chain
//     (`getTreeHelper().calcDistanceMatrix(df, ['molecule'], 'euclidean')` →
//     `getClusterMatrixWorker(matrix.data, rowCount, 4 /* centroid */)`),
//     NOT a test-bug. The scenario's "(or median)" alternative is the
//     correct path for the UI smoke — median on molecule mounts in ~3.1s
//     with no console error. The centroid-molecule bug is reported via the
//     dispatch-yaml mcp_observations[] (decision-log category
//     core-bug-candidate) for triage; full distance×linkage coverage stays
//     owned by the sibling apitest, which exercises every combination via
//     the JS API and surfaces the same crash structurally if it persists.
//
// Selectors per .claude/skills/grok-browser/references/dendrogram.md
// (rev 2026-06-03 live-MCP-validated). All selectors used here are class-1
// (in-reference) — no Selector recon-notes block needed.
//
// MCP recon evidence (live 2026-06-03 on dev.datagrok.ai, user oahadzhanian,
// mol1K.csv, dialog defaults verified; centroid-molecule core-bug + non-fatal
// 404 filter added on retry):
//   - dialog [name="dialog-Hierarchical-Clustering"] opens via Chem | Analyze |
//     Hierarchical Clustering... (synthetic-event nav chain validated).
//   - [name="input-Distance"] is a SELECT with options ["euclidean", "manhattan"]
//     in that order; default "euclidean". Matches scenario Block A step 3.
//   - [name="input-Linkage"] is a SELECT with options ["single", "complete",
//     "average", "weighted", "centroid", "median", "ward"] in that order;
//     default "ward". Matches scenario Block A step 4 verbatim.
//   - [name="input-host-Features"].textContent === "Features(1) molecule".
//     Matches scenario step 2 "Features defaulting to molecule".
//   - [name="input-Table"] SELECT defaults to "Table" (the DataFrame name
//     assigned by grok.dapi.files.readCsv — NOT "mol1K"). The scenario's
//     literal "Table = mol1K" is the dataset filename, not the in-memory
//     DataFrame name; we assert "Table" SELECT has exactly one option
//     and the dialog's input-host-Table is populated.
//   - Block B mount path is identical to the assign-clusters-spec sibling
//     (poll for .dendrogram-assign-clusters-bttn — magic-wand is the
//     mounted-and-ready signal). euclidean+ward mounted in ~3.7s, median
//     in ~3.1s on warm session. The retained 30s budget absorbs cold-init.
//   - CENTROID + MOLECULE features triggers a FATAL platform crash
//     "Error: memory access out of bounds" at hierarchical-clustering.ts:217
//     (catch-and-console.error of a WASM OOB inside the centroid-linkage
//     compute). Magic-wand never mounts; tested 90s budget on warm session.
//     This is a core-bug, surfaced via mcp_observations[] in the retry
//     dispatch yaml. SR-03 swaps Step 7's centroid to median per the
//     scenario's explicit "centroid (or median)" alternative.
//   - Resource-load 404s ("Failed to load resource: the server responded
//     with a status of 404") are surfaced by Chromium on every Dendrogram
//     neighbor mount on dev.datagrok.ai (consistent / reproducible). They
//     are NOT application-level errors — the platform serves them on every
//     visit independent of Hierarchical-Clustering. Scenario Step 6/7 says
//     "no fatal console errors" — these 404s are non-fatal by definition
//     (page is fully functional). Step 5 says "no console errors" but the
//     scenario authority's Notes section is broader than the literal step
//     wording; the same non-fatal filter applies symmetrically. See the
//     isFatalConsoleError() helper for the exclusion predicate.
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

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
  // MCP-validated ~2s on mol1K; budget 30s to absorb cold-init variance.
  const foundAtMs: number = await page.evaluate(async () => {
    const start = Date.now();
    for (let i = 0; i < 60; i++) {
      if (document.querySelector('.dendrogram-assign-clusters-bttn'))
        return Date.now() - start;
      await new Promise(r => setTimeout(r, 500));
    }
    return -1;
  });
  return foundAtMs;
}

async function closeDendrogramNeighbor(page: Page): Promise<void> {
  await page.evaluate(async () => {
    const close = document.querySelector('.dendrogram-close-bttn') as HTMLElement | null;
    if (!close) throw new Error('.dendrogram-close-bttn not found (neighbor not mounted)');
    close.click();
    for (let i = 0; i < 20; i++) {
      if (!document.querySelector('.dendrogram-assign-clusters-bttn')) return;
      await new Promise(r => setTimeout(r, 200));
    }
    throw new Error('Neighbor did not detach within 4s');
  });
}

async function setDialogSelect(page: Page, name: 'Distance' | 'Linkage', value: string): Promise<void> {
  await page.evaluate(([n, v]: [string, string]) => {
    const sel = document.querySelector(`[name="dialog-Hierarchical-Clustering"] [name="input-${n}"]`) as HTMLSelectElement | null;
    if (!sel) throw new Error(`input-${n} SELECT not found`);
    sel.value = v;
    sel.dispatchEvent(new Event('input', {bubbles: true}));
    sel.dispatchEvent(new Event('change', {bubbles: true}));
  }, [name, value]);
}

// Returns true when the console-message text is a fatal application error and
// false for routine non-fatal noise (Chromium's "Failed to load resource: …
// status of 404 ()" lines emitted on every Dendrogram-neighbor mount on
// dev.datagrok.ai, ResizeObserver loop limits, etc.). Scenario Notes treat
// "no console errors" / "no fatal console errors" as the same fatal-only
// predicate — the literal-text resource 404 is not actionable from the spec's
// vantage and is reproducible independent of the clustering code path.
function isFatalConsoleError(text: string): boolean {
  if (/Failed to load resource[\s\S]*404/i.test(text)) return false;
  if (/ResizeObserver loop/i.test(text)) return false;
  return true;
}

test('Dendrogram: Hierarchical Clustering (Chem) — dialog gateway + representative OK runs', async ({page}) => {
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

  // --- Block A — Dialog exposes all Distance and Linkage values ---

  await softStep('1. Open mol1K and verify molecule column rendered', async () => {
    const info = await page.evaluate(() => {
      const df = grok.shell.tv?.dataFrame;
      const molCol = df?.col('molecule');
      return {rows: df?.rowCount, molSemType: molCol?.semType};
    });
    expect(info.rows, 'mol1K row count').toBe(1000);
    expect(info.molSemType, 'molecule semType').toBe('Molecule');
  });

  await softStep('2. Run Chem | Analyze | Hierarchical Clustering... — dialog opens with Features=molecule', async () => {
    await openHierarchicalClusteringDialog(page);
    const defaults = await page.evaluate(() => ({
      dialogPresent: !!document.querySelector('[name="dialog-Hierarchical-Clustering"]'),
      table: (document.querySelector('[name="input-Table"]') as HTMLSelectElement)?.value,
      tableOptionCount: (document.querySelector('[name="input-Table"]') as HTMLSelectElement)?.options.length,
      features: document.querySelector('[name="input-host-Features"]')?.textContent?.trim(),
      distancePresent: !!document.querySelector('[name="input-Distance"]'),
      linkagePresent: !!document.querySelector('[name="input-Linkage"]'),
      okBtn: !!document.querySelector('[name="dialog-Hierarchical-Clustering"] [name="button-OK"]'),
    }));
    expect(defaults.dialogPresent, 'Hierarchical Clustering dialog opened').toBe(true);
    // Table SELECT: exactly one TableView is open. The DataFrame name assigned
    // by readCsv is "Table" (not "mol1K") — see header MCP recon note.
    expect(defaults.tableOptionCount, 'Table SELECT has one option (one TableView open)').toBe(1);
    expect(defaults.features, 'Features defaults to molecule').toContain('molecule');
    expect(defaults.distancePresent, 'Distance input visible').toBe(true);
    expect(defaults.linkagePresent, 'Linkage input visible').toBe(true);
    expect(defaults.okBtn, 'OK button present').toBe(true);
  });

  await softStep('3. Distance dropdown enumerates exactly [euclidean, manhattan] (default euclidean)', async () => {
    const distance = await page.evaluate(() => {
      const sel = document.querySelector('[name="input-Distance"]') as HTMLSelectElement | null;
      return {
        default: sel?.value,
        options: sel ? Array.from(sel.options).map(o => o.value) : null,
      };
    });
    expect(distance.options, 'Distance options exact list+order').toEqual(['euclidean', 'manhattan']);
    expect(distance.default, 'Distance default').toBe('euclidean');
  });

  await softStep('4. Linkage dropdown enumerates exactly 7 values in spec order (default ward)', async () => {
    const linkage = await page.evaluate(() => {
      const sel = document.querySelector('[name="input-Linkage"]') as HTMLSelectElement | null;
      return {
        default: sel?.value,
        options: sel ? Array.from(sel.options).map(o => o.value) : null,
      };
    });
    expect(linkage.options, 'Linkage options exact list+order').toEqual([
      'single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward',
    ]);
    expect(linkage.default, 'Linkage default').toBe('ward');
  });

  // --- Block B — Representative end-to-end runs (spot-check) ---

  await softStep('5. OK with euclidean+ward (Features=molecule) → dendrogram neighbor injected, no fatal console error', async () => {
    // Defaults are already euclidean+ward+molecule; click OK and wait for the
    // magic-wand mount (the deterministic ready signal per
    // grok-browser/dendrogram.md § common observability).
    const consoleErrors: string[] = [];
    const listener = (msg: any) => { if (msg.type() === 'error') consoleErrors.push(msg.text()); };
    page.on('console', listener);
    try {
      const foundAtMs = await clickOkAndWaitForNeighbor(page);
      expect(foundAtMs, 'Magic-wand mount time (ms; -1 = timeout)').toBeGreaterThan(0);
      const state = await page.evaluate(() => ({
        magicWand: !!document.querySelector('.dendrogram-assign-clusters-bttn'),
        closeBtn: !!document.querySelector('.dendrogram-close-bttn'),
        neighborHasCanvas: !!document.querySelector('.dendrogram-assign-clusters-bttn')
          ?.parentElement?.querySelector('canvas'),
        viewerTypes: Array.from(grok.shell.tv.viewers).map((v: any) => v.type),
      }));
      expect(state.magicWand, 'magic-wand icon present').toBe(true);
      expect(state.closeBtn, 'close icon present').toBe(true);
      expect(state.neighborHasCanvas, 'neighbor canvas mounted').toBe(true);
      // GridNeighbor is NOT a DG.Viewer per grok-browser/dendrogram.md.
      expect(state.viewerTypes, 'viewer types list (neighbor is NOT a Viewer)').toEqual(['Grid']);
      // Scenario step 5 says "No console errors" — interpreted as fatal-only
      // per Notes section authority (non-fatal resource 404s are emitted on
      // every dev.datagrok.ai mount independent of the clustering code path;
      // see isFatalConsoleError() helper above).
      const fatalErrors = consoleErrors.filter(isFatalConsoleError);
      expect(fatalErrors, 'no fatal console errors on euclidean+ward run').toEqual([]);
    } finally {
      page.off('console', listener);
    }
  });

  await softStep('6. Close, re-open dialog, set manhattan+single (Features=molecule), OK → second dendrogram builds', async () => {
    await closeDendrogramNeighbor(page);
    await openHierarchicalClusteringDialog(page);
    await setDialogSelect(page, 'Distance', 'manhattan');
    await setDialogSelect(page, 'Linkage', 'single');
    // Features defaults to molecule on the Chem leaf; no change needed.
    const verify = await page.evaluate(() => ({
      distance: (document.querySelector('[name="input-Distance"]') as HTMLSelectElement)?.value,
      linkage: (document.querySelector('[name="input-Linkage"]') as HTMLSelectElement)?.value,
      features: document.querySelector('[name="input-host-Features"]')?.textContent?.trim(),
    }));
    expect(verify.distance, 'Distance set to manhattan').toBe('manhattan');
    expect(verify.linkage, 'Linkage set to single').toBe('single');
    expect(verify.features, 'Features still includes molecule').toContain('molecule');

    const consoleErrors: string[] = [];
    const unsupportedTypeErrors: string[] = [];
    const listener = (msg: any) => {
      if (msg.type() === 'error') {
        consoleErrors.push(msg.text());
        if (/Unsupported\s+column\s+type/i.test(msg.text())) unsupportedTypeErrors.push(msg.text());
      }
    };
    page.on('console', listener);
    try {
      const foundAtMs = await clickOkAndWaitForNeighbor(page);
      expect(foundAtMs, 'second dendrogram mount time (ms; -1 = timeout)').toBeGreaterThan(0);
      const state = await page.evaluate(() => ({
        magicWand: !!document.querySelector('.dendrogram-assign-clusters-bttn'),
        // Exactly one neighbor — the replace-on-rerun path closes the previous
        // before injecting the new one (per atlas edge_case + assign-clusters-spec
        // Scenario 5; not the focus here, but a useful invariant).
        neighborCount: document.querySelectorAll('.dendrogram-assign-clusters-bttn').length,
      }));
      expect(state.magicWand, 'magic-wand icon present after manhattan+single run').toBe(true);
      expect(state.neighborCount, 'exactly one neighbor attached').toBe(1);
      // Scenario step 6 explicit assertions:
      expect(unsupportedTypeErrors, 'no "Unsupported column type" error').toEqual([]);
      const fatalErrors = consoleErrors.filter(isFatalConsoleError);
      expect(fatalErrors, 'no fatal console errors on manhattan+single run').toEqual([]);
    } finally {
      page.off('console', listener);
    }
  });

  await softStep('7. Close, re-open dialog, set euclidean+median (centroid (or median) per scenario; default Features=molecule per SR-02), OK → dendrogram builds', async () => {
    await closeDendrogramNeighbor(page);
    await openHierarchicalClusteringDialog(page);
    await setDialogSelect(page, 'Distance', 'euclidean');
    // SR-03 (see header): scenario says "centroid (or median)". Live MCP recon
    // 2026-06-03 established that centroid + molecule features triggers a fatal
    // WASM "memory access out of bounds" platform crash at
    // hierarchical-clustering.ts:217; the dendrogram never mounts. Median (the
    // scenario's explicit alternative) mounts cleanly in ~3s on warm session.
    // The full distance × linkage matrix coverage lives in the sibling apitest.
    await setDialogSelect(page, 'Linkage', 'median');
    // SR-02 (see header): Features stays at the default `molecule` selection —
    // the numeric-features substitution from the scenario step 7 is deferred
    // to the sibling apitest because the column-list picker is canvas and
    // has no verified DOM toggle path. Median linkage on MOLECULE features
    // still exercises the non-ward path the scenario step targets.
    const verify = await page.evaluate(() => ({
      distance: (document.querySelector('[name="input-Distance"]') as HTMLSelectElement)?.value,
      linkage: (document.querySelector('[name="input-Linkage"]') as HTMLSelectElement)?.value,
      features: document.querySelector('[name="input-host-Features"]')?.textContent?.trim(),
    }));
    expect(verify.distance, 'Distance set to euclidean').toBe('euclidean');
    expect(verify.linkage, 'Linkage set to median (centroid (or median) per scenario; SR-03)').toBe('median');
    expect(verify.features, 'Features remains at default molecule (SR-02)').toContain('molecule');

    const consoleErrors: string[] = [];
    const unsupportedTypeErrors: string[] = [];
    const listener = (msg: any) => {
      if (msg.type() === 'error') {
        consoleErrors.push(msg.text());
        if (/Unsupported\s+column\s+type/i.test(msg.text())) unsupportedTypeErrors.push(msg.text());
      }
    };
    page.on('console', listener);
    try {
      const foundAtMs = await clickOkAndWaitForNeighbor(page);
      // SR-01: the non-monotonic-tree visual property is NOT asserted — see
      // header. The scenario downgrades the Step 7 assertion to "dendrogram
      // builds without error", which is what we exercise here.
      expect(foundAtMs, 'median mount time (ms; -1 = timeout)').toBeGreaterThan(0);
      const state = await page.evaluate(() => ({
        magicWand: !!document.querySelector('.dendrogram-assign-clusters-bttn'),
        neighborCount: document.querySelectorAll('.dendrogram-assign-clusters-bttn').length,
      }));
      expect(state.magicWand, 'magic-wand icon present after median run').toBe(true);
      expect(state.neighborCount, 'exactly one neighbor attached').toBe(1);
      // Scenario step 7 explicit assertions:
      expect(unsupportedTypeErrors, 'no "Unsupported column type" error on median path').toEqual([]);
      const fatalErrors = consoleErrors.filter(isFatalConsoleError);
      expect(fatalErrors, 'no fatal console errors on euclidean+median run').toEqual([]);
    } finally {
      page.off('console', listener);
    }
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
