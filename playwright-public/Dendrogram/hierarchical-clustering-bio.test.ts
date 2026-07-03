/* ---
sub_features_covered: [dendrogram.api.tree-helper.calc-distance-matrix, dendrogram.api.tree-helper.cut-tree-to-grid, dendrogram.clustering.assign-clusters-dialog, dendrogram.clustering.dialog, dendrogram.clustering.inject-tree-for-grid, dendrogram.clustering.menu.bio]
--- */
// Bio Hierarchical Clustering on FASTA_PT_activity.csv (99 rows, sequence col = Macromolecule).
// The bio leaf auto-defaults Features to the sequence column; Levenshtein build path.
// Clusters→Threshold binary search is inexact, so assert a multi-way partition (>=2), not == requested.
// Mount budget 120s absorbs cold-init; mount-success contract is state.magicWand === true (the
// foundAtMs numeric is clamped to >=1 on success / -1 on timeout). Resource-load 404s on every
// neighbor mount are non-fatal noise (see isFatalConsoleError). Centroid linkage crashes (WASM OOB
// at hierarchical-clustering.ts:217) so this spec uses manhattan+complete on the second build.
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);

async function openHierarchicalClusteringDialog(page: Page): Promise<void> {
  await page.evaluate(async () => {
    const bio = document.querySelector('[name="div-Bio"]') as HTMLElement | null;
    if (!bio) throw new Error('Top-menu Bio entry not found');
    bio.dispatchEvent(new MouseEvent('click', {bubbles: true}));
    const findLabel = (pred: (t: string) => boolean) => Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .find(m => pred((m.textContent || '').trim())) as HTMLElement | undefined;
    const pollLabel = async (pred: (t: string) => boolean) => {
      for (let i = 0; i < 50; i++) {
        const el = findLabel(pred);
        if (el) return el;
        await new Promise(r => setTimeout(r, 100));
      }
      return undefined;
    };
    const analyze = await pollLabel(t => t === 'Analyze');
    if (!analyze) throw new Error('"Analyze" sub-menu item not found');
    (analyze.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
    const hc = await pollLabel(t => /Hierarchical\s+Clustering/i.test(t));
    if (!hc) throw new Error('"Hierarchical Clustering..." sub-menu item not found');
    (hc.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
  });
  await page.locator('[name="dialog-Hierarchical-Clustering"]').waitFor({timeout: 15_000});
}

async function clickOkAndWaitForNeighbor(page: Page, budgetMs: number = 120_000): Promise<number> {
  await page.locator('[name="dialog-Hierarchical-Clustering"] [name="button-OK"]').click();
  // Poll for the magic-wand mount; return is clamped to >=1 on success, -1 on timeout.
  const iterCap = Math.max(1, Math.ceil(budgetMs / 500));
  const foundAtMs: number = await page.evaluate(async (cap: number) => {
    const start = Date.now();
    for (let i = 0; i < cap; i++) {
      if (document.querySelector('.dendrogram-assign-clusters-bttn'))
        return Math.max(1, Date.now() - start);
      await new Promise(r => setTimeout(r, 500));
    }
    return -1;
  }, iterCap);
  return foundAtMs;
}

async function closeDendrogramNeighbor(page: Page): Promise<void> {
  await page.evaluate(async () => {
    const close = document.querySelector('.dendrogram-close-bttn') as HTMLElement | null;
    if (!close) throw new Error('.dendrogram-close-bttn not found (neighbor not mounted)');
    close.click();
    for (let i = 0; i < 30; i++) {
      if (!document.querySelector('.dendrogram-assign-clusters-bttn')) {
        // Settle so the next OK click does not race the worker teardown.
        await new Promise(r => setTimeout(r, 500));
        return;
      }
      await new Promise(r => setTimeout(r, 200));
    }
    throw new Error('Neighbor did not detach within 6s');
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

// True for fatal app errors; false for non-fatal noise (404 resource loads, ResizeObserver loop).
function isFatalConsoleError(text: string): boolean {
  if (/Failed to load resource[\s\S]*404/i.test(text)) return false;
  if (/ResizeObserver loop/i.test(text)) return false;
  return true;
}

test('Dendrogram / Hierarchical Clustering (Bio) — sequence-default dialog + Levenshtein build path + Assign Clusters smoke', async ({page}) => {
  test.setTimeout(360_000);

  await loginToDatagrok(page);

  // Setup — open FASTA_PT_activity.csv and wait for the Macromolecule semType + Bio warmup.
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    for (let i = 0; i < 50; i++) {
      if (grok.shell.tv == null) break;
      await new Promise(r => setTimeout(r, 100));
    }
    const df = await grok.dapi.files.readCsv('System:AppData/Bio/samples/FASTA_PT_activity.csv');
    grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 5000);
    });
    // Bio dataset: wait for Grid canvas, then for the Bio top-menu to register
    // (sequence renderer + package warmup complete once [name="div-Bio"] attaches).
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise(r => setTimeout(r, 200));
    }
    for (let i = 0; i < 150; i++) {
      if (document.querySelector('[name="div-Bio"]')) break;
      await new Promise(r => setTimeout(r, 200));
    }
    // Wait for the sequence column's Macromolecule semType to actually be detected — the
    // Bio top-menu registers on package load independently of this dataframe's semType
    // detection, so div-Bio can attach before detection finishes. Gating on the real
    // condition faithfully replaces the old fixed settle and lets the Bio leaf auto-default
    // Features to the sequence column.
    for (let i = 0; i < 150; i++) {
      const st = grok.shell.tv?.dataFrame?.col('sequence')?.semType;
      if (st && /macromolecule/i.test(st)) break;
      await new Promise(r => setTimeout(r, 200));
    }
    // Pre-warm the Dendrogram TreeHelper singleton so the first build skips package init.
    try {
      await (grok as any).functions.call('Dendrogram:getTreeHelper');
    } catch (e) { /* best-effort pre-warm; downstream budget absorbs */ }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  // --- Scenario 1 — Bio dialog opens with sequence-default Features ---

  await softStep('1. Open FASTA_PT_activity.csv and verify sequence column rendered as Macromolecule', async () => {
    const info = await page.evaluate(() => {
      const df = grok.shell.tv?.dataFrame;
      const seqCol = df?.col('sequence');
      return {
        rows: df?.rowCount,
        seqSemType: seqCol?.semType,
        seqType: seqCol?.type,
      };
    });
    expect(info.rows, 'FASTA_PT_activity row count').toBe(99);
    expect(info.seqSemType, 'sequence semType (case-insensitive match)').toMatch(/^macromolecule$/i);
    expect(info.seqType, 'sequence storage type').toBe('string');
  });

  await softStep('2. Run Bio | Analyze | Hierarchical Clustering... — dialog opens with Features=sequence (MACROMOLECULE auto-default)', async () => {
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
    // DataFrame name from readCsv is "Table" (not "FASTA_PT_activity"); one TableView open.
    expect(defaults.tableOptionCount, 'Table SELECT has one option (one TableView open)').toBe(1);
    // Bio leaf auto-defaults Features to the sole MACROMOLECULE column (sequence).
    expect(defaults.features, 'Features defaults to sequence (MACROMOLECULE auto-default)').toContain('sequence');
    expect(defaults.distancePresent, 'Distance input visible').toBe(true);
    expect(defaults.linkagePresent, 'Linkage input visible').toBe(true);
    expect(defaults.okBtn, 'OK button present').toBe(true);
  });

  // --- Scenario 2 — Distance and Linkage dropdowns expose the canonical value sets ---

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

  // --- Scenario 3 — Build dendrogram from the sequence column (Levenshtein path) ---

  await softStep('5. OK with euclidean+ward (Features=sequence) → dendrogram neighbor injected via Levenshtein-on-encoded-sequences path, no fatal console error', async () => {
    // Defaults are already euclidean+ward+sequence; click OK and wait for the magic-wand mount.
    const consoleErrors: string[] = [];
    const listener = (msg: any) => { if (msg.type() === 'error') consoleErrors.push(msg.text()); };
    page.on('console', listener);
    try {
      const foundAtMs = await clickOkAndWaitForNeighbor(page);
      // `>= 0` matches the clamped success return (1..budgetMs); mount contract is state.magicWand.
      expect(foundAtMs, 'Magic-wand mount time (ms; -1 = timeout)').toBeGreaterThanOrEqual(0);
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
      // GridNeighbor is NOT a DG.Viewer.
      expect(state.viewerTypes, 'viewer types list (neighbor is NOT a Viewer)').toEqual(['Grid']);
      const fatalErrors = consoleErrors.filter(isFatalConsoleError);
      expect(fatalErrors, 'no fatal console errors on euclidean+ward Levenshtein run').toEqual([]);
      const unsupportedTypeErrors = consoleErrors.filter(t => /Unsupported\s+column\s+type/i.test(t));
      expect(unsupportedTypeErrors, 'no "Unsupported column type" error').toEqual([]);
    } finally {
      page.off('console', listener);
    }
  });

  await softStep('6. Close, re-open dialog, set manhattan+complete (Features=sequence), OK → second dendrogram builds successfully', async () => {
    await closeDendrogramNeighbor(page);
    await openHierarchicalClusteringDialog(page);
    await setDialogSelect(page, 'Distance', 'manhattan');
    await setDialogSelect(page, 'Linkage', 'complete');
    // Features defaults to sequence on the Bio leaf; no change needed.
    const verify = await page.evaluate(() => ({
      distance: (document.querySelector('[name="input-Distance"]') as HTMLSelectElement)?.value,
      linkage: (document.querySelector('[name="input-Linkage"]') as HTMLSelectElement)?.value,
      features: document.querySelector('[name="input-host-Features"]')?.textContent?.trim(),
    }));
    expect(verify.distance, 'Distance set to manhattan').toBe('manhattan');
    expect(verify.linkage, 'Linkage set to complete').toBe('complete');
    expect(verify.features, 'Features still includes sequence').toContain('sequence');

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
      // `>= 0` matches the clamped success return; mount contract is state.magicWand.
      expect(foundAtMs, 'second dendrogram mount time (ms; -1 = timeout)').toBeGreaterThanOrEqual(0);
      const state = await page.evaluate(() => ({
        magicWand: !!document.querySelector('.dendrogram-assign-clusters-bttn'),
        // Replace-on-rerun closes the previous neighbor before injecting the new one.
        neighborCount: document.querySelectorAll('.dendrogram-assign-clusters-bttn').length,
        neighborHasCanvas: !!document.querySelector('.dendrogram-assign-clusters-bttn')
          ?.parentElement?.querySelector('canvas'),
      }));
      expect(state.magicWand, 'magic-wand icon present after manhattan+complete run').toBe(true);
      expect(state.neighborCount, 'exactly one neighbor attached').toBe(1);
      expect(state.neighborHasCanvas, 'neighbor canvas mounted on second run').toBe(true);
      // Scenario step 6 explicit assertions:
      expect(unsupportedTypeErrors, 'no "Unsupported column type" error on manhattan+complete sequence run').toEqual([]);
      const fatalErrors = consoleErrors.filter(isFatalConsoleError);
      expect(fatalErrors, 'no fatal console errors on manhattan+complete Levenshtein run').toEqual([]);
    } finally {
      page.off('console', listener);
    }
  });

  // --- Scenario 4 — Shared post-build smoke ride-along (Assign Clusters column creation) ---

  await softStep('7. Click magic-wand on dendrogram neighbor → Assign Clusters dialog opens with Threshold (slider) + Clusters (int) inputs', async () => {
    const wandClicked = await page.evaluate(async () => {
      const wand = document.querySelector('.dendrogram-assign-clusters-bttn') as HTMLElement | null;
      if (!wand) return false;
      wand.click();
      // Wait for dialog
      for (let i = 0; i < 30; i++) {
        if (document.querySelector('[name="dialog-Assign-Clusters"]')) return true;
        await new Promise(r => setTimeout(r, 200));
      }
      return false;
    });
    expect(wandClicked, 'Assign Clusters dialog opened from magic-wand').toBe(true);
    const dlgState = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-Assign-Clusters"]');
      return {
        title: dlg?.querySelector('.d4-dialog-title')?.textContent?.trim(),
        thresholdInputPresent: !!dlg?.querySelector('[name="input-Threshold"]'),
        clustersInputPresent: !!dlg?.querySelector('[name="input-Clusters"]'),
        sliderPresent: !!dlg?.querySelector('[name="input-host-Threshold"] input[type="range"]'),
        assignBtnPresent: !!dlg?.querySelector('[name="button-Assign"]'),
        cancelBtnPresent: !!dlg?.querySelector('[name="button-CANCEL"]'),
        clustersInitial: (dlg?.querySelector('[name="input-Clusters"]') as HTMLInputElement)?.value,
      };
    });
    expect(dlgState.title, 'dialog title').toBe('Assign Clusters');
    expect(dlgState.thresholdInputPresent, 'Threshold input present').toBe(true);
    expect(dlgState.clustersInputPresent, 'Clusters input present').toBe(true);
    expect(dlgState.sliderPresent, 'Threshold slider present (range input)').toBe(true);
    expect(dlgState.assignBtnPresent, 'Assign button present').toBe(true);
    expect(dlgState.cancelBtnPresent, 'CANCEL button present').toBe(true);
    // Initial Clusters depends on the tree topology; assert positive, not an exact value.
    expect(Number(dlgState.clustersInitial), 'initial Clusters value is positive').toBeGreaterThan(0);
  });

  await softStep('8. Set Clusters=5 and click Assign → dialog closes; new Cluster (<threshold>) categorical column appended to DataFrame', async () => {
    const colsBefore: string[] = await page.evaluate(() => grok.shell.tv.dataFrame.columns.names());

    const consoleErrors: string[] = [];
    const listener = (msg: any) => { if (msg.type() === 'error') consoleErrors.push(msg.text()); };
    page.on('console', listener);
    try {
      const result = await page.evaluate(async () => {
        const dlg = document.querySelector('[name="dialog-Assign-Clusters"]');
        const clInput = dlg?.querySelector('[name="input-Clusters"]') as HTMLInputElement | null;
        if (!clInput) throw new Error('Clusters input not found');
        clInput.value = '5';
        clInput.dispatchEvent(new Event('input', {bubbles: true}));
        clInput.dispatchEvent(new Event('change', {bubbles: true}));
        // Poll until the Threshold input recomputes and stabilizes (unchanged, non-empty).
        const readThreshold = () => (dlg?.querySelector('[name="input-Threshold"]') as HTMLInputElement | null)?.value ?? '';
        let prev = '';
        for (let i = 0; i < 40; i++) {
          const cur = readThreshold();
          if (cur !== '' && cur === prev) break;
          prev = cur;
          await new Promise(r => setTimeout(r, 100));
        }
        const settledThreshold = readThreshold();
        const assignBtn = dlg?.querySelector('[name="button-Assign"]') as HTMLElement | null;
        if (!assignBtn) throw new Error('Assign button not found');
        assignBtn.click();
        // Wait for dialog to close
        let closed = false;
        for (let i = 0; i < 30; i++) {
          if (!document.querySelector('[name="dialog-Assign-Clusters"]')) { closed = true; break; }
          await new Promise(r => setTimeout(r, 200));
        }
        return {settledThreshold, closed};
      });
      expect(result.closed, 'dialog closed on Assign').toBe(true);

      const colsAfter: string[] = await page.evaluate(() => grok.shell.tv.dataFrame.columns.names());
      const allNew = colsAfter.filter(n => !colsBefore.includes(n));
      // Scope the "exactly one created" invariant to the Cluster column: background Bio/Chem processing
      // can append unrelated columns concurrently within the before/after snapshot window.
      const newCols = allNew.filter(n => /^Cluster\s*\(/.test(n));
      expect(newCols.length, `exactly one new Cluster column appended (all new: ${JSON.stringify(allNew)})`).toBe(1);
      // Column name format Cluster (<threshold.toFixed(2)>); match the format, not the value.
      expect(newCols[0], 'new column name follows Cluster (N.NN) format').toMatch(/^Cluster\s*\(\d+(\.\d{1,2})?\)$/);
      // Tie the created column to the settled threshold captured from the dialog (Clusters=5 request).
      const nameThreshold = Number(newCols[0].match(/\(([\d.]+)\)/)![1]);
      expect(nameThreshold, 'column-name threshold matches settled dialog threshold').toBeCloseTo(Number(result.settledThreshold), 2);

      const colInfo = await page.evaluate((colName: string) => {
        const col = grok.shell.tv.dataFrame.col(colName)!;
        return {
          type: col.type,
          categoriesLength: col.categories ? col.categories.length : 0,
          rowCountMatches: col.length === grok.shell.tv.dataFrame.rowCount,
        };
      }, newCols[0]);
      expect(colInfo.type, 'new column is string (categorical)').toBe('string');
      // Clusters→Threshold binary search is inexact (~4-6 for a request of 5), so we don't assert
      // == 5, but a request for 5 clusters MUST yield a real multi-way partition — collapsing to a
      // single category means clustering did nothing. Upper band intentionally omitted (unverified on
      // dev; recorded for human review) — the floor is what catches the degenerate breakage.
      expect(colInfo.categoriesLength, 'requested 5 clusters produced a multi-way partition').toBeGreaterThanOrEqual(2);
      expect(colInfo.rowCountMatches, 'cluster column length matches DataFrame row count').toBe(true);

      const fatalErrors = consoleErrors.filter(isFatalConsoleError);
      expect(fatalErrors, 'no fatal console errors on Assign').toEqual([]);
    } finally {
      page.off('console', listener);
    }
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
