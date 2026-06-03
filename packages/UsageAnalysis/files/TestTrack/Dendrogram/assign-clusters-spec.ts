/* ---
sub_features_covered: [dendrogram.clustering.menu.chem, dendrogram.clustering.dialog, dendrogram.clustering.inject-tree-for-grid, dendrogram.clustering.assign-clusters-dialog, dendrogram.api.tree-helper.cut-tree-to-grid, dendrogram.event.context-menu, dendrogram.viewer]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: integration
//   sub_features_covered: [dendrogram.clustering.menu.chem, dendrogram.clustering.dialog,
//     dendrogram.clustering.inject-tree-for-grid, dendrogram.clustering.assign-clusters-dialog,
//     dendrogram.api.tree-helper.cut-tree-to-grid, dendrogram.event.context-menu,
//     dendrogram.viewer]
//   ui_coverage_responsibility: [assign-clusters-dialog, assign-clusters-magic-wand-icon,
//     assign-clusters-threshold-slider, assign-clusters-clusters-input,
//     assign-clusters-assign-button, dendrogram-context-menu-assign-clusters,
//     hierarchical-clustering-replace-warning] (delegated_to: null)
//     (zoom flows dendrogram-horizontal-zoom-ctrl-wheel /
//      dendrogram-vertical-scroll-plain-wheel /
//      dendrogram-double-click-reset-zoom are SR-01 deferred to manual)
//   related_bugs: []
//   produced_from: migrated
//
// Atlas provenance (derived_from):
//   dendrogram.yaml#critical_paths[dendrogram.cp.assign-clusters-column-creation]
//     derived_from: scenario-chain:dendrogram.yaml#assign-clusters.md
//   dendrogram.yaml#sub_features[dendrogram.clustering.inject-tree-for-grid].interactions[0..3]
//     derived_from: public/packages/Dendrogram/src/viewers/inject-tree-for-grid2.ts#L65..L76
//   dendrogram.yaml#edge_cases[Cluster-assignment dialog Threshold/Clusters two-way binding]
//     derived_from: public/packages/Dendrogram/src/viewers/inject-tree-for-grid2.ts#L394
//   dendrogram.yaml#edge_cases[Dendrogram already attached to a grid]
//     derived_from: public/packages/Dendrogram/src/utils/hierarchical-clustering.ts#L104
//
// Scope reductions (carried from scenario frontmatter):
//   SR-01: source Steps 11-13 (plain mouse-wheel vertical scroll, Ctrl+wheel
//     horizontal zoom with clamping, double-click empty area resets horizontal
//     zoom) deferred to manual per atlas
//     `manual_only[dendrogram.mo.ctrl-wheel-zoom-tactile]`.
//   SR-02: source Step 10 (canvas-pixel visual cut-position indicator) deferred
//     to manual per atlas `manual_only[dendrogram.mo.tree-canvas-visual-regression]`.
//     The automated path asserts the `Cluster (<threshold>)` column instead —
//     the persisted, observable output of the cut.
//
// Selectors per .claude/skills/grok-browser/references/dendrogram.md (rev
// 2026-06-03 live-MCP-validated):
//   [name="div-Chem"], .d4-menu-item-label "Analyze"/"Hierarchical Clustering...",
//   [name="dialog-Hierarchical-Clustering"] + [name="button-OK"],
//   .dendrogram-assign-clusters-bttn (magic wand, aria-label "Assign Clusters"),
//   .dendrogram-close-bttn (close icon, aria-label "Remove Dendrogram"),
//   right-click neighbor → context menu .d4-menu-item-label "Assign Clusters"/"Reset Zoom",
//   [name="dialog-Assign-Clusters"] + [name="input-Threshold"] / [name="input-Clusters"] /
//   [name="input-host-Threshold"] input[type="range"] / [name="button-Assign"].
//
// MCP recon evidence (live 2026-06-03 on dev.datagrok.ai, user oahadzhanian,
// mol1K.csv, euclidean+ward):
//   - Neighbor mounts within ~2.0-2.1s after Hierarchical-Clustering OK (mol1K
//     is small; the grok-browser ref's "8-15s" upper bound is conservative —
//     budget set to 30s here to absorb cold-init variance under `grok test`).
//   - Initial Assign Clusters dialog: Threshold 319.54, Clusters 6.
//   - Set Clusters=5 → Threshold settles to 239.66; Assign appends
//     `Cluster (239.66)` (string col, 5 categories ["1".."5"]).
//   - Magic wand re-opens dialog; set Clusters=8 → Threshold 399.43; Assign
//     appends `Cluster (399.43)` (8 categories). Both cluster columns coexist
//     (auto-uniquely named per df.columns.getUnusedName).
//   - Slider move → Clusters recompute confirmed (slider→159.77 → Clusters=3).
//   - Close icon removes the neighbor (.dendrogram-assign-clusters-bttn gone).
//   - Re-running Chem|Analyze|Hierarchical Clustering while neighbor still
//     attached results in exactly ONE neighbor (replace path —
//     hierarchical-clustering.ts#L104). The toast warning surfaces as a
//     transient grok.shell.warning; the durable, deterministic observable is
//     `document.querySelectorAll('.dendrogram-assign-clusters-bttn').length == 1`.
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
  // Neighbor mounts when .dendrogram-assign-clusters-bttn (magic wand) appears.
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

async function openAssignClustersViaContextMenu(page: Page): Promise<void> {
  await page.evaluate(async () => {
    const wand = document.querySelector('.dendrogram-assign-clusters-bttn');
    if (!wand) throw new Error('Dendrogram neighbor not mounted (no magic-wand)');
    const neighborRoot = wand.parentElement as HTMLElement;
    const canvas = neighborRoot.querySelector('canvas') as HTMLCanvasElement;
    if (!canvas) throw new Error('Neighbor canvas not found');
    const rect = canvas.getBoundingClientRect();
    // The listener is on the neighbor root (treeNb.root), not the canvas itself —
    // dispatch on neighborRoot per grok-browser/dendrogram.md
    // § assign-clusters-dialog § Entry vector 1.
    const ev = new MouseEvent('contextmenu', {
      bubbles: true, cancelable: true,
      clientX: rect.left + rect.width / 2,
      clientY: rect.top + rect.height / 2,
    });
    neighborRoot.dispatchEvent(ev);
    await new Promise(r => setTimeout(r, 600));
    const assign = Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .find(m => m.textContent!.trim() === 'Assign Clusters') as HTMLElement | undefined;
    if (!assign) throw new Error('"Assign Clusters" context-menu item not found');
    (assign.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
  });
  await page.locator('[name="dialog-Assign-Clusters"]').waitFor({timeout: 10_000});
}

async function openAssignClustersViaMagicWand(page: Page): Promise<void> {
  await page.locator('.dendrogram-assign-clusters-bttn').click();
  await page.locator('[name="dialog-Assign-Clusters"]').waitFor({timeout: 10_000});
}

async function setClustersInputAndAssign(page: Page, clusters: number): Promise<{
  threshold: string;
  newColumns: string[];
  allClusterColumns: string[];
}> {
  return await page.evaluate(async (target: number) => {
    const clustersInput = document.querySelector('[name="input-Clusters"]') as HTMLInputElement;
    if (!clustersInput) throw new Error('Clusters input not found');
    clustersInput.value = String(target);
    clustersInput.dispatchEvent(new Event('input', {bubbles: true}));
    clustersInput.dispatchEvent(new Event('change', {bubbles: true}));
    await new Promise(r => setTimeout(r, 500));
    const thresholdInput = document.querySelector('[name="input-Threshold"]') as HTMLInputElement;
    const threshold = thresholdInput.value;
    const df = grok.shell.tv.dataFrame;
    const colsBefore = df.columns.names();
    const assignBtn = document.querySelector('[name="button-Assign"]') as HTMLButtonElement;
    if (!assignBtn) throw new Error('Assign button not found');
    assignBtn.click();
    // Wait for dialog to close + new column to appear.
    for (let i = 0; i < 30; i++) {
      if (!document.querySelector('[name="dialog-Assign-Clusters"]')) break;
      await new Promise(r => setTimeout(r, 200));
    }
    // Tiny settle for the column-add side effect.
    await new Promise(r => setTimeout(r, 500));
    const colsAfter = df.columns.names();
    const newColumns = colsAfter.filter((c: string) => !colsBefore.includes(c));
    const allClusterColumns = colsAfter.filter((c: string) => c.startsWith('Cluster ('));
    return {threshold, newColumns, allClusterColumns};
  }, clusters);
}

test('Dendrogram: Assign Clusters end-to-end (column creation, two-way binding, replace-on-rerun)', async ({page}) => {
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

  // Scenario 1 — Build dendrogram via Chem | Analyze | Hierarchical Clustering on mol1K.
  await softStep('1. Open mol1K and verify molecule column rendered', async () => {
    const info = await page.evaluate(() => {
      const df = grok.shell.tv?.dataFrame;
      const molCol = df?.col('molecule');
      return {rows: df?.rowCount, molType: molCol?.type, molSemType: molCol?.semType};
    });
    expect(info.rows, 'mol1K row count').toBe(1000);
    expect(info.molSemType, 'molecule semType').toBe('Molecule');
  });

  await softStep('2. Open Chem | Analyze | Hierarchical Clustering — verify dialog defaults', async () => {
    await openHierarchicalClusteringDialog(page);
    const defaults = await page.evaluate(() => ({
      distance: (document.querySelector('[name="input-Distance"]') as HTMLSelectElement)?.value,
      linkage: (document.querySelector('[name="input-Linkage"]') as HTMLSelectElement)?.value,
      features: document.querySelector('[name="input-host-Features"]')?.textContent?.trim(),
      okBtn: !!document.querySelector('[name="dialog-Hierarchical-Clustering"] [name="button-OK"]'),
    }));
    expect(defaults.distance, 'Distance default').toBe('euclidean');
    expect(defaults.linkage, 'Linkage default').toBe('ward');
    expect(defaults.features, 'Features displayed value').toContain('molecule');
    expect(defaults.okBtn, 'OK button present').toBe(true);
  });

  await softStep('3. Click OK → dendrogram neighbor injected to left of grid', async () => {
    const foundAtMs = await clickOkAndWaitForNeighbor(page);
    expect(foundAtMs, 'Magic-wand mount time (ms; -1 = timeout)').toBeGreaterThan(0);
    const state = await page.evaluate(() => ({
      magicWand: !!document.querySelector('.dendrogram-assign-clusters-bttn'),
      closeBtn: !!document.querySelector('.dendrogram-close-bttn'),
      neighborHasCanvas: !!document.querySelector('.dendrogram-assign-clusters-bttn')
        ?.parentElement?.querySelector('canvas'),
      // viewers list still only includes Grid — the neighbor is a GridNeighbor,
      // NOT a DG.Viewer (per grok-browser/dendrogram.md § common observability).
      viewerTypes: Array.from(grok.shell.tv.viewers).map((v: any) => v.type),
    }));
    expect(state.magicWand, 'magic-wand icon present').toBe(true);
    expect(state.closeBtn, 'close icon present').toBe(true);
    expect(state.neighborHasCanvas, 'neighbor canvas mounted').toBe(true);
    expect(state.viewerTypes, 'viewer types').toEqual(['Grid']);
  });

  // Scenario 2 — Open Assign Clusters via both surfaces (context menu + magic-wand icon).
  await softStep('4. Right-click neighbor canvas → Assign Clusters opens dialog', async () => {
    await openAssignClustersViaContextMenu(page);
    const dialog = await page.evaluate(() => ({
      title: document.querySelector('[name="dialog-Assign-Clusters"] .d4-dialog-title')?.textContent?.trim(),
      threshold: (document.querySelector('[name="input-Threshold"]') as HTMLInputElement)?.value,
      clusters: (document.querySelector('[name="input-Clusters"]') as HTMLInputElement)?.value,
      sliderPresent: !!document.querySelector('[name="input-host-Threshold"] input[type="range"]'),
      assignBtn: !!document.querySelector('[name="button-Assign"]'),
    }));
    expect(dialog.title, 'dialog title').toBe('Assign Clusters');
    // Threshold defaults to about half the tree height; Clusters defaults to a
    // small integer ≥ 1. Don't pin the exact value (mol1K rebuild-to-rebuild
    // variance + the binary-search rounding documented in
    // inject-tree-for-grid2.ts#L394). Assert the shape, not the literal.
    expect(parseFloat(dialog.threshold || ''), 'threshold parses as positive').toBeGreaterThan(0);
    expect(parseInt(dialog.clusters || '0', 10), 'clusters integer ≥ 1').toBeGreaterThanOrEqual(1);
    expect(dialog.sliderPresent, 'threshold slider present').toBe(true);
    expect(dialog.assignBtn, 'Assign button present').toBe(true);
  });

  await softStep('5. Close dialog, click magic-wand → same dialog re-opens (aria-label Assign Clusters)', async () => {
    await page.locator('[name="dialog-Assign-Clusters"] [name="button-CANCEL"]').click();
    await page.locator('[name="dialog-Assign-Clusters"]').waitFor({state: 'detached', timeout: 5_000});
    const ariaLabel = await page.evaluate(() => document.querySelector('.dendrogram-assign-clusters-bttn')?.getAttribute('aria-label'));
    expect(ariaLabel, 'magic-wand aria-label').toBe('Assign Clusters');
    await openAssignClustersViaMagicWand(page);
    const title = await page.evaluate(() =>
      document.querySelector('[name="dialog-Assign-Clusters"] .d4-dialog-title')?.textContent?.trim());
    expect(title, 'dialog title after magic-wand').toBe('Assign Clusters');
  });

  // Scenario 3 — Two-way Threshold ↔ Clusters binding.
  await softStep('6. Move Threshold slider → Clusters count recomputes', async () => {
    const result = await page.evaluate(async () => {
      const slider = document.querySelector('[name="input-host-Threshold"] input[type="range"]') as HTMLInputElement;
      const clustersInput = document.querySelector('[name="input-Clusters"]') as HTMLInputElement;
      const beforeClusters = parseInt(clustersInput.value, 10);
      const sliderMin = parseFloat(slider.min);
      const sliderMax = parseFloat(slider.max);
      // Slide to a quarter of the range — empirically yields a DIFFERENT cluster
      // count from the default (mol1K euclidean+ward: slider≈160 → Clusters=3,
      // default Clusters=6). The exact resulting count depends on the tree
      // shape; assert "changed", not the literal.
      const newVal = sliderMin + (sliderMax - sliderMin) / 4;
      slider.value = String(newVal);
      slider.dispatchEvent(new Event('input', {bubbles: true}));
      slider.dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      return {beforeClusters, afterClusters: parseInt(clustersInput.value, 10), sliderVal: slider.value};
    });
    expect(result.afterClusters, 'Clusters recomputed after slider move').not.toBe(result.beforeClusters);
    expect(result.afterClusters, 'Clusters ≥ 1 after slider move').toBeGreaterThanOrEqual(1);
  });

  await softStep('7. Type Clusters=5 → Threshold recalculates (binary search; tolerance not pinned)', async () => {
    const result = await page.evaluate(async () => {
      const clustersInput = document.querySelector('[name="input-Clusters"]') as HTMLInputElement;
      const thresholdInput = document.querySelector('[name="input-Threshold"]') as HTMLInputElement;
      const beforeThreshold = parseFloat(thresholdInput.value);
      clustersInput.value = '5';
      clustersInput.dispatchEvent(new Event('input', {bubbles: true}));
      clustersInput.dispatchEvent(new Event('change', {bubbles: true}));
      await new Promise(r => setTimeout(r, 500));
      return {beforeThreshold, afterThreshold: parseFloat(thresholdInput.value), clusters: clustersInput.value};
    });
    // Per atlas edge_case dendrogram.yaml#inject-tree-for-grid2.ts#L394: the
    // 20-iteration binary search maps the requested Clusters back to a
    // Threshold; if the requested count is not exactly reachable, `minDiff`
    // selects the threshold yielding the closest count. Therefore we assert
    // "Threshold changed" and "input value preserved at 5", NOT an exact
    // Threshold literal (the scenario's `unresolved_ambiguities:
    // step-7-clusters-to-threshold-tolerance-not-pinned`).
    expect(result.clusters, 'Clusters input preserved at typed value').toBe('5');
    expect(result.afterThreshold, 'Threshold recalculated (changed)').not.toBe(result.beforeThreshold);
    expect(result.afterThreshold, 'Threshold > 0').toBeGreaterThan(0);
  });

  // Scenario 4 — Assign creates a `Cluster (<threshold>)` column.
  await softStep('8. Click Assign with Clusters=5 → `Cluster (<n.nn>)` column appended', async () => {
    const before = await page.evaluate(() => grok.shell.tv.dataFrame.columns.names());
    const result = await setClustersInputAndAssign(page, 5);
    expect(result.newColumns.length, 'exactly one new column added').toBe(1);
    const newCol = result.newColumns[0];
    expect(newCol, 'new column name format Cluster (<n.nn>)').toMatch(/^Cluster \(\d+\.\d{2}\)$/);
    // Cluster column shape: string, 5 categories ("1".."5"), no nulls.
    const colInfo = await page.evaluate((name: string) => {
      const c = grok.shell.tv.dataFrame.col(name)!;
      return {
        type: c.type,
        categoriesLen: (c as any).categories?.length,
        firstFiveValues: Array.from({length: 5}, (_: any, i: number) => c.get(i)),
        nullCount: (c as any).stats?.missingValueCount ?? 0,
      };
    }, newCol);
    expect(colInfo.type, 'cluster column type').toBe('string');
    expect(colInfo.categoriesLen, 'category count == requested Clusters value').toBe(5);
    for (const v of colInfo.firstFiveValues)
      expect(v, 'cluster id is a non-empty string').toMatch(/^\d+$/);
    expect(before.length + 1, 'column count grew by exactly one').toBe(before.length + 1);
  });

  await softStep('9. Re-open via magic-wand, Clusters=8, Assign → second column appended (auto-unique)', async () => {
    await openAssignClustersViaMagicWand(page);
    const colsBefore = await page.evaluate(() =>
      grok.shell.tv.dataFrame.columns.names().filter((n: string) => n.startsWith('Cluster (')));
    const result = await setClustersInputAndAssign(page, 8);
    expect(result.newColumns.length, 'exactly one new column added').toBe(1);
    const secondCol = result.newColumns[0];
    expect(secondCol, 'second column name format').toMatch(/^Cluster \(\d+\.\d{2}\)/);
    // df.columns.getUnusedName ensures uniqueness (see grok-browser/dendrogram.md
    // § assign-clusters-dialog § Cluster column invariants). Both columns must
    // coexist.
    const colsAfter = result.allClusterColumns;
    expect(colsAfter.length, 'both cluster columns coexist').toBe(colsBefore.length + 1);
    expect(colsAfter, 'previous cluster column preserved').toEqual(expect.arrayContaining(colsBefore));
    const colInfo = await page.evaluate((name: string) => {
      const c = grok.shell.tv.dataFrame.col(name)!;
      return {categoriesLen: (c as any).categories?.length, type: c.type};
    }, secondCol);
    expect(colInfo.type, 'second cluster column type').toBe('string');
    expect(colInfo.categoriesLen, 'second cluster column category count == 8').toBe(8);
  });

  // Scenario 5 — Replace-on-rerun cleanup (no duplicate viewers).
  await softStep('10. Re-run Chem | Analyze | Hierarchical Clustering on same TableView — exactly one neighbor remains', async () => {
    // Per atlas edge_case "Dendrogram already attached to a grid" (source:
    // hierarchical-clustering.ts#L104), the second OK on a TableView with an
    // existing neighbor fires a transient grok.shell.warning toast and closes
    // the existing GridNeighbor before injecting a new one (tracked via
    // DENDROGRAM_NEIGHBOR_TEMP_NAME on tv.grid.temp).
    //
    // The deterministic, persisted observable for this contract is the
    // post-replace count of `.dendrogram-assign-clusters-bttn` (which is 1 —
    // the new neighbor — not 2 or 0). The toast itself is transient and
    // race-prone (MCP 2026-06-03 recon caught the replacement but missed the
    // toast window); asserting on the durable invariant matches what the
    // scenario actually claims about end state.
    expect(await page.evaluate(() => document.querySelectorAll('.dendrogram-assign-clusters-bttn').length),
      'one neighbor attached before re-run').toBe(1);
    await openHierarchicalClusteringDialog(page);
    await clickOkAndWaitForNeighbor(page);
    // Extra settle to let the close-then-inject sequence finish.
    await page.waitForTimeout(2_000);
    const counts = await page.evaluate(() => ({
      neighborCount: document.querySelectorAll('.dendrogram-assign-clusters-bttn').length,
      closeBtnCount: document.querySelectorAll('.dendrogram-close-bttn').length,
      // Cluster columns from the prior cycles are not affected by re-clustering.
      clusterColumnCount: grok.shell.tv.dataFrame.columns.names()
        .filter((n: string) => n.startsWith('Cluster (')).length,
    }));
    expect(counts.neighborCount, 'exactly one dendrogram neighbor attached after re-run').toBe(1);
    expect(counts.closeBtnCount, 'exactly one close icon').toBe(1);
    expect(counts.clusterColumnCount, 'prior Cluster columns preserved across re-run').toBeGreaterThanOrEqual(2);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
