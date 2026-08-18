import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);

async function openHierarchicalClusteringDialog(page: Page): Promise<void> {
  await page.evaluate(async () => {
    const bio = document.querySelector('[name="div-Bio"]') as HTMLElement | null;
    if (!bio) throw new Error('Top-menu Bio entry not found');
    bio.dispatchEvent(new MouseEvent('click', {bubbles: true}));
    await new Promise(r => setTimeout(r, 800));
    const analyze = Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .find(m => m.textContent!.trim() === 'Analyze') as HTMLElement | undefined;
    if (!analyze) throw new Error('"Analyze" sub-menu item not found');
    (analyze.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
    await new Promise(r => setTimeout(r, 700));
    const hc = Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .find(m => /Hierarchical\s+Clustering/i.test(m.textContent || '')) as HTMLElement | undefined;
    if (!hc) throw new Error('"Hierarchical Clustering..." sub-menu item not found');
    (hc.closest('.d4-menu-item') as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
  });
  await page.locator('[name="dialog-Hierarchical-Clustering"]').waitFor({timeout: 15_000});
}

async function clickOkAndWaitForNeighbor(page: Page, budgetMs: number = 120_000): Promise<number> {
  await page.locator('[name="dialog-Hierarchical-Clustering"] [name="button-OK"]').click();

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

function isFatalConsoleError(text: string): boolean {
  if (/Failed to load resource[\s\S]*404/i.test(text)) return false;
  if (/ResizeObserver loop/i.test(text)) return false;
  return true;
}

test('Dendrogram: Hierarchical Clustering (Bio) — sequence-default dialog + Levenshtein build path + Assign Clusters smoke', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    await new Promise(r => setTimeout(r, 1000));
    const df = await grok.dapi.files.readCsv('System:AppData/Bio/samples/FASTA_PT_activity.csv');
    grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 5000);
    });

    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise(r => setTimeout(r, 200));
    }
    await new Promise(r => setTimeout(r, 6000));

    try {
      await (grok as any).functions.call('Dendrogram:getTreeHelper');
    } catch (e) {  }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

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

    expect(defaults.tableOptionCount, 'Table SELECT has one option (one TableView open)').toBe(1);

    expect(defaults.features, 'Features defaults to sequence (MACROMOLECULE auto-default)').toContain('sequence');
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

  await softStep('5. OK with euclidean+ward (Features=sequence) → dendrogram neighbor injected via Levenshtein-on-encoded-sequences path, no fatal console error', async () => {

    const consoleErrors: string[] = [];
    const listener = (msg: any) => { if (msg.type() === 'error') consoleErrors.push(msg.text()); };
    page.on('console', listener);
    try {
      const foundAtMs = await clickOkAndWaitForNeighbor(page);

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

      expect(foundAtMs, 'second dendrogram mount time (ms; -1 = timeout)').toBeGreaterThanOrEqual(0);
      const state = await page.evaluate(() => ({
        magicWand: !!document.querySelector('.dendrogram-assign-clusters-bttn'),

        neighborCount: document.querySelectorAll('.dendrogram-assign-clusters-bttn').length,
        neighborHasCanvas: !!document.querySelector('.dendrogram-assign-clusters-bttn')
          ?.parentElement?.querySelector('canvas'),
      }));
      expect(state.magicWand, 'magic-wand icon present after manhattan+complete run').toBe(true);
      expect(state.neighborCount, 'exactly one neighbor attached').toBe(1);
      expect(state.neighborHasCanvas, 'neighbor canvas mounted on second run').toBe(true);

      expect(unsupportedTypeErrors, 'no "Unsupported column type" error on manhattan+complete sequence run').toEqual([]);
      const fatalErrors = consoleErrors.filter(isFatalConsoleError);
      expect(fatalErrors, 'no fatal console errors on manhattan+complete Levenshtein run').toEqual([]);
    } finally {
      page.off('console', listener);
    }
  });

  await softStep('7. Click magic-wand on dendrogram neighbor → Assign Clusters dialog opens with Threshold (slider) + Clusters (int) inputs', async () => {
    const wandClicked = await page.evaluate(async () => {
      const wand = document.querySelector('.dendrogram-assign-clusters-bttn') as HTMLElement | null;
      if (!wand) return false;
      wand.click();

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
        await new Promise(r => setTimeout(r, 800));
        const settledThreshold = (dlg?.querySelector('[name="input-Threshold"]') as HTMLInputElement | null)?.value;
        const assignBtn = dlg?.querySelector('[name="button-Assign"]') as HTMLElement | null;
        if (!assignBtn) throw new Error('Assign button not found');
        assignBtn.click();

        let closed = false;
        for (let i = 0; i < 30; i++) {
          if (!document.querySelector('[name="dialog-Assign-Clusters"]')) { closed = true; break; }
          await new Promise(r => setTimeout(r, 200));
        }
        return {settledThreshold, closed};
      });
      expect(result.closed, 'dialog closed on Assign').toBe(true);

      const colsAfter: string[] = await page.evaluate(() => grok.shell.tv.dataFrame.columns.names());
      const newCols = colsAfter.filter(n => !colsBefore.includes(n));
      expect(newCols.length, 'exactly one new column appended').toBe(1);

      expect(newCols[0], 'new column name follows Cluster (N.NN) format').toMatch(/^Cluster\s*\(\d+(\.\d{1,2})?\)$/);

      const colInfo = await page.evaluate((colName: string) => {
        const col = grok.shell.tv.dataFrame.col(colName)!;
        return {
          type: col.type,
          categoriesLength: col.categories ? col.categories.length : 0,
          rowCountMatches: col.length === grok.shell.tv.dataFrame.rowCount,
        };
      }, newCols[0]);
      expect(colInfo.type, 'new column is string (categorical)').toBe('string');

      expect(colInfo.categoriesLength, 'cluster column has at least one category').toBeGreaterThan(0);
      expect(colInfo.rowCountMatches, 'cluster column length matches DataFrame row count').toBe(true);

      const fatalErrors = consoleErrors.filter(isFatalConsoleError);
      expect(fatalErrors, 'no fatal console errors on Assign').toEqual([]);
    } finally {
      page.off('console', listener);
    }
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
