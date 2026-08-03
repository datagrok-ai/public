import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
const DATASET_PATH = 'System:AppData/Bio/tests/filter_FASTA.csv';
test('Bio Similarity Search docks KNN viewer + row-click re-queries neighbours', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(() => resolve(), 4000);
    });
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const hasMacromolecule = cols.some((c: any) => c.semType === 'Macromolecule');
    if (hasMacromolecule) {
      for (let i = 0; i < 60; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
    }
  }, DATASET_PATH);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  await page.evaluate(async () => {
    const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
    for (const fn of probes) {
      try { await (grok as any).functions.call(fn, {}); return; } catch {  }
    }
    await new Promise((r) => setTimeout(r, 3000));
  });
  await page.evaluate(async () => {
    const candidates = ['Bio:similaritySearchTopMenu', 'Bio:similaritySearch'];
    const findAny = (names: string[]): boolean => {
      for (const n of names) {
        try {
          if ((grok as any).functions.find && (grok as any).functions.find(n)) return true;
        } catch {  }
      }
      return false;
    };
    const deadline = Date.now() + 15_000;
    while (Date.now() < deadline) {
      if (findAny(candidates)) return;
      await new Promise((r) => setTimeout(r, 300));
    }
    await new Promise((r) => setTimeout(r, 1500));
  });
  await page.waitForTimeout(2000);
  const setupProbe = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const macro = cols.find((c: any) => c.semType === 'Macromolecule') as any;
    return {
      hasMacromoleculeCol: !!macro,
      macroName: macro?.name ?? null,
      macroSemType: macro?.semType ?? null,
      rowCount: df.rowCount,
      currentRowIdx: df.currentRowIdx,
    };
  });
  expect(setupProbe.hasMacromoleculeCol,
    'atlas bio.detector contract: readCsv path MUST classify a Macromolecule column synchronously').toBe(true);
  expect(setupProbe.macroSemType).toBe('Macromolecule');
  expect(setupProbe.rowCount,
    'scenario .md Setup: ≥ 5 rows so KNN K=3..5 is meaningful and a non-current row exists').toBeGreaterThanOrEqual(5);
  await softStep('Scenario 1.1: click Bio > Search > Similarity Search (no dialog — viewer docks directly)', async () => {
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      const group = document.querySelector('[name="div-Bio---Search"]');
      if (!group) throw new Error('[name="div-Bio---Search"] group anchor not found');
      group.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      group.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      const leaf = document.querySelector('[name="div-Bio---Search---Similarity-Search"]');
      if (!leaf) throw new Error('[name="div-Bio---Search---Similarity-Search"] leaf not found under Bio > Search');
      (leaf as HTMLElement).click();
    });
  });
  await softStep('Scenario 1.2-3: viewer Sequence Similarity Search docks; initial KNN populated', async () => {
    await page.locator('[name="viewer-Sequence-Similarity-Search"]').waitFor({timeout: 60_000});
    await page.waitForFunction(() => {
      const viewers = Array.from((grok.shell.tv as any).viewers) as any[];
      const v = viewers.find((vw) => vw.type === 'Sequence Similarity Search');
      if (!v) return false;
      return (v as any).idxs != null && (v as any).scores != null;
    }, null, {timeout: 180_000});
    const dockProbe = await page.evaluate(() => {
      const viewers = Array.from((grok.shell.tv as any).viewers) as any[];
      const v = viewers.find((vw) => vw.type === 'Sequence Similarity Search');
      const viewerTypes = viewers.map((vw) => vw.type);
      let knnRowCount: number | null = null;
      try {
        const idxsCol = (v as any)?.idxs;
        if (idxsCol && typeof idxsCol.length === 'number') knnRowCount = idxsCol.length;
      } catch {  }
      const targetMoleculeIdx = (v as any)?.targetMoleculeIdx ?? null;
      return {
        viewerPresent: !!v,
        viewerType: v?.type ?? null,
        viewerTypes,
        knnRowCount,
        targetMoleculeIdx,
        currentRowIdx: grok.shell.tv.dataFrame.currentRowIdx,
        viewerDockedDOM: !!document.querySelector('[name="viewer-Sequence-Similarity-Search"]'),
      };
    });
    expect(dockProbe.viewerPresent,
      `Expected 'Sequence Similarity Search' viewer in tv.viewers; actual types: ${JSON.stringify(dockProbe.viewerTypes)}`).toBe(true);
    expect(dockProbe.viewerType).toBe('Sequence Similarity Search');
    expect(dockProbe.viewerDockedDOM,
      'viewer dock container [name="viewer-Sequence-Similarity-Search"] MUST be present in DOM').toBe(true);
    expect(dockProbe.knnRowCount,
      'Viewer KNN compute MUST populate idxs with >=2 rows (target + at least one neighbour).' +
      ` Observed: ${dockProbe.knnRowCount}`).not.toBeNull();
    expect(dockProbe.knnRowCount!).toBeGreaterThanOrEqual(2);
  });
  // Settle after Scenario 1's compute so Scenario 2's render isn't queued behind a pending renderPromise.
  await page.waitForTimeout(1500);
  // Scenario 2 — Clicking a different row re-queries the KNN viewer.
  const scenario1Baseline = await page.evaluate(() => {
    const viewers = Array.from((grok.shell.tv as any).viewers) as any[];
    const v = viewers.find((vw) => vw.type === 'Sequence Similarity Search');
    const targetMoleculeIdx = (v as any)?.targetMoleculeIdx ?? null;
    return {
      targetMoleculeIdx,
      currentRowIdx: grok.shell.tv.dataFrame.currentRowIdx,
    };
  });
  await softStep('Scenario 2.1: click a non-current row in the grid (sets df.currentRowIdx → fires onCurrentRowChanged)', async () => {
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const newIdx = df.rowCount - 1;
      if (newIdx === df.currentRowIdx) throw new Error(
        `Cannot pick a distinct non-current row: rowCount=${df.rowCount}, currentRowIdx=${df.currentRowIdx}`);
      df.currentRowIdx = newIdx;
    });
  });
  await softStep('Scenario 2.2: viewer re-queries — targetMoleculeIdx changes from Scenario 1 baseline', async () => {
    await page.waitForFunction((baseline) => {
      const viewers = Array.from((grok.shell.tv as any).viewers) as any[];
      const v = viewers.find((vw) => vw.type === 'Sequence Similarity Search');
      if (!v) return false;
      const targetMoleculeIdx = (v as any)?.targetMoleculeIdx ?? null;
      return targetMoleculeIdx !== null
        && targetMoleculeIdx !== (baseline as any).targetMoleculeIdx;
    }, scenario1Baseline, {timeout: 180_000});
    await page.waitForFunction(() => {
      const viewers = Array.from((grok.shell.tv as any).viewers) as any[];
      const v = viewers.find((vw) => vw.type === 'Sequence Similarity Search');
      if (!v) return false;
      return (v as any).idxs != null && (v as any).scores != null;
    }, null, {timeout: 180_000});
    const scenario2Probe = await page.evaluate(() => {
      const viewers = Array.from((grok.shell.tv as any).viewers) as any[];
      const v = viewers.find((vw) => vw.type === 'Sequence Similarity Search');
      const targetMoleculeIdx = (v as any)?.targetMoleculeIdx ?? null;
      let knnRowCount: number | null = null;
      try {
        const idxsCol = (v as any)?.idxs;
        if (idxsCol && typeof idxsCol.length === 'number') knnRowCount = idxsCol.length;
      } catch { /* keep null */ }
      return {
        targetMoleculeIdx,
        knnRowCount,
        currentRowIdx: grok.shell.tv.dataFrame.currentRowIdx,
        viewerPresent: !!v,
      };
    });
    expect(scenario2Probe.viewerPresent,
      'Viewer MUST remain docked after row-click re-query').toBe(true);
    expect(scenario2Probe.targetMoleculeIdx,
      `Re-query MUST surface in viewer state. Scenario 1 targetMoleculeIdx: ` +
      `${scenario1Baseline.targetMoleculeIdx}. Scenario 2: ` +
      `${scenario2Probe.targetMoleculeIdx}. ` +
      `Viewer did not react to df.currentRowIdx change.`)
      .not.toBe(scenario1Baseline.targetMoleculeIdx);
    expect(scenario2Probe.knnRowCount,
      `Re-queried KNN MUST populate >=2 rows. Observed: ${scenario2Probe.knnRowCount}`)
      .not.toBeNull();
    expect(scenario2Probe.knnRowCount!).toBeGreaterThanOrEqual(2);
    expect(scenario2Probe.currentRowIdx).toBe(setupProbe.rowCount - 1);
    const dfStable = await page.evaluate(() => !!grok.shell.tv.dataFrame);
    expect(dfStable).toBe(true);
  });
  finishSpec();
});
