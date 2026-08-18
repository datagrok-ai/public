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

function isFatalConsoleError(text: string): boolean {
  if (/Failed to load resource[\s\S]*404/i.test(text)) return false;
  if (/ResizeObserver loop/i.test(text)) return false;
  return true;
}

test('Dendrogram: Hierarchical Clustering (Chem) — dialog gateway + representative OK runs', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

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

    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise(r => setTimeout(r, 200));
    }
    await new Promise(r => setTimeout(r, 5000));
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

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

  await softStep('5. OK with euclidean+ward (Features=molecule) → dendrogram neighbor injected, no fatal console error', async () => {

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

      expect(state.viewerTypes, 'viewer types list (neighbor is NOT a Viewer)').toEqual(['Grid']);
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

        neighborCount: document.querySelectorAll('.dendrogram-assign-clusters-bttn').length,
      }));
      expect(state.magicWand, 'magic-wand icon present after manhattan+single run').toBe(true);
      expect(state.neighborCount, 'exactly one neighbor attached').toBe(1);

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

    await setDialogSelect(page, 'Linkage', 'median');

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

      expect(foundAtMs, 'median mount time (ms; -1 = timeout)').toBeGreaterThan(0);
      const state = await page.evaluate(() => ({
        magicWand: !!document.querySelector('.dendrogram-assign-clusters-bttn'),
        neighborCount: document.querySelectorAll('.dendrogram-assign-clusters-bttn').length,
      }));
      expect(state.magicWand, 'magic-wand icon present after median run').toBe(true);
      expect(state.neighborCount, 'exactly one neighbor attached').toBe(1);

      expect(unsupportedTypeErrors, 'no "Unsupported column type" error on median path').toEqual([]);
      const fatalErrors = consoleErrors.filter(isFatalConsoleError);
      expect(fatalErrors, 'no fatal console errors on euclidean+median run').toEqual([]);
    } finally {
      page.off('console', listener);
    }
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
