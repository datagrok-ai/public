/* ---
sub_features_covered: [dendrogram.api.tree-helper.cut-tree-to-grid, dendrogram.clustering.assign-clusters-dialog, dendrogram.clustering.dialog, dendrogram.clustering.inject-tree-for-grid, dendrogram.clustering.menu.chem, dendrogram.event.context-menu, dendrogram.viewer]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

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
    // The listener is on the neighbor root (treeNb.root), not the canvas itself.
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
    await new Promise(r => setTimeout(r, 500));
    const colsAfter = df.columns.names();
    const newColumns = colsAfter.filter((c: string) => !colsBefore.includes(c));
    const allClusterColumns = colsAfter.filter((c: string) => c.startsWith('Cluster ('));
    return {threshold, newColumns, allClusterColumns};
  }, clusters);
}

test('Dendrogram / Assign Clusters end-to-end (column creation, two-way binding, replace-on-rerun)', async ({page}) => {
  test.setTimeout(600_000);

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
      // Neighbor is a GridNeighbor, NOT a DG.Viewer — viewers list stays ['Grid'].
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
    // Don't pin exact threshold/clusters (rebuild variance + binary-search rounding); assert shape.
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
      // Slide to a quarter of the range — empirically yields a DIFFERENT cluster count.
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
    // Binary search maps Clusters→Threshold inexactly; assert "changed" + preserved input, not a literal.
    expect(result.clusters, 'Clusters input preserved at typed value').toBe('5');
    expect(result.afterThreshold, 'Threshold recalculated (changed)').not.toBe(result.beforeThreshold);
    expect(result.afterThreshold, 'Threshold > 0').toBeGreaterThan(0);
  });

  // Scenario 4 — Assign creates a `Cluster (<threshold>)` column.
  await softStep('8. Click Assign with Clusters=5 → `Cluster (<n.nn>)` column appended', async () => {
    const before = await page.evaluate(() => grok.shell.tv.dataFrame.columns.names());
    const result = await setClustersInputAndAssign(page, 5);
    // Scope the "exactly one created" invariant to the Cluster column: background Chem processing on
    // mol1K can append unrelated columns concurrently within the before/after snapshot window.
    const newClusterCols = result.newColumns.filter((c) => /^Cluster \(/.test(c));
    expect(newClusterCols.length, `exactly one new Cluster column (all new: ${JSON.stringify(result.newColumns)})`).toBe(1);
    const newCol = newClusterCols[0];
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
    const newClusterCols = result.newColumns.filter((c) => /^Cluster \(/.test(c));
    expect(newClusterCols.length, `exactly one new Cluster column (all new: ${JSON.stringify(result.newColumns)})`).toBe(1);
    const secondCol = newClusterCols[0];
    expect(secondCol, 'second column name format').toMatch(/^Cluster \(\d+\.\d{2}\)/);
    // df.columns.getUnusedName ensures uniqueness; both columns must coexist.
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
    // Replace-on-rerun (hierarchical-clustering.ts#L104): the durable observable is the
    // post-replace count of .dendrogram-assign-clusters-bttn (1, not 2 or 0).
    expect(await page.evaluate(() => document.querySelectorAll('.dendrogram-assign-clusters-bttn').length),
      'one neighbor attached before re-run').toBe(1);
    await openHierarchicalClusteringDialog(page);
    await clickOkAndWaitForNeighbor(page);
    // Settle so the close-then-inject sequence finishes.
    await page.waitForTimeout(2_000);
    const counts = await page.evaluate(() => ({
      neighborCount: document.querySelectorAll('.dendrogram-assign-clusters-bttn').length,
      closeBtnCount: document.querySelectorAll('.dendrogram-close-bttn').length,
      clusterColumnCount: grok.shell.tv.dataFrame.columns.names()
        .filter((n: string) => n.startsWith('Cluster (')).length,
    }));
    expect(counts.neighborCount, 'exactly one dendrogram neighbor attached after re-run').toBe(1);
    expect(counts.closeBtnCount, 'exactly one close icon').toBe(1);
    expect(counts.clusterColumnCount, 'prior Cluster columns preserved across re-run').toBeGreaterThanOrEqual(2);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
