/* ---
sub_features_covered:
  - bio.viewers.similarity-search
  - bio.search.similarity
  - bio.search.similarity.top-menu
  - bio.detector
  - bio.rendering
--- */
//   related_bugs: [] (positive-path; GROK-16111 empty-input edge contract is
//     "Related-bug context")
// open via the readCsv entry path (bio.md "GROK-18616 entry-path-class invariant");
// Retry-1 hypothesis (round 1 of three-touchpoint loop) — test-bug category,
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
test.use(specTestOptions);
const DATASET_PATH = 'System:AppData/Bio/tests/filter_FASTA.csv';
test('Bio Similarity Search docks KNN viewer + row-click re-queries neighbours', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  // bio.md "GROK-18616 entry-path-class invariant" — the readCsv path triggers
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
    // Scenario 1 Expected — "The viewer displays a non-empty K-nearest-
    // neighbours panel — i.e. at least one row card / list entry beyond
    // the current row reference is shown". For a 14-row dataset with K=10
    // default the contract is K+1 = 11 (target + neighbours per
    // sequence-similarity-viewer.ts L162); the assertion >= 2 covers the
    // "at least one row beyond the reference" floor and is independent of
    // the K=10 default being preserved across Bio versions.
    expect(dockProbe.knnRowCount,
      'Viewer KNN compute MUST populate idxs with >=2 rows (target + at least one neighbour).' +
      ` Observed: ${dockProbe.knnRowCount}`).not.toBeNull();
    expect(dockProbe.knnRowCount!).toBeGreaterThanOrEqual(2);
  });
  // Retry-1 tactical fix: explicit settle after Scenario 1's compute so
  // Scenario 2's render is not queued behind a still-pending
  // renderPromise (search-base-viewer.ts L78-85 chains renderInt onto
  // renderPromise). Without this settle, the Scenario 2 row-click
  // dispatch can queue BEHIND the Scenario 1 computeByMM tail and
  // targetMoleculeIdx will lag the click by the residual Scenario 1
  // compute time, racing the Scenario 2 waitForFunction window.
  await page.waitForTimeout(1500);
  // Scenario 2 — Clicking a different row re-queries the KNN viewer.
  // Capture the Scenario 1 baseline: the index the viewer last computed
  // against (targetMoleculeIdx — sequence-similarity-viewer.ts L29, L75).
  // Retry-1 tactical fix: dropped the prior topResultIdx baseline read
  // (depended on the unreliable embedded-grid-DOM-to-viewer.dataFrame
  // path; the baseline was effectively null, so the Scenario 2 polling
  // degraded to targetMoleculeIdx-only anyway). targetMoleculeIdx is the
  // canonical per-source signal.
  const scenario1Baseline = await page.evaluate(() => {
    const viewers = Array.from((grok.shell.tv as any).viewers) as any[];
    const v = viewers.find((vw) => vw.type === 'Sequence Similarity Search');
    const targetMoleculeIdx = (v as any)?.targetMoleculeIdx ?? null;
    return {
      targetMoleculeIdx,
      currentRowIdx: grok.shell.tv.dataFrame.currentRowIdx,
    };
  });
  // Step 2.1: click a different row in the underlying table grid.
  // Scenario .md Scenario 2 step 1: "click a different (non-current) row
  // in the underlying table grid". Grid cells are canvas-painted and not
  // DOM-addressable; the canonical class-1 path per bio.md and the
  // user-facing grid contract is setting df.currentRowIdx (the same
  // mutation a canvas grid click performs internally). Per scenario .md
  // Scenario 2 Expected: "the currentRow indicator in the table tracks
  // the click (atlas bio.search.similarity viewer-current-row bridge —
  // the scenario level assertion is on the visible viewer reaction, not
  // on grid.dataFrame.currentRowIdx itself; that is apitest-layer)" —
  // i.e. the spec-layer assertion is on the viewer reaction, and the
  // currentRow set is the trigger. Pick a non-current row deterministically
  // (rowCount-1) so the new index is guaranteed distinct from the
  // Scenario 1 default (0).
  await softStep('Scenario 2.1: click a non-current row in the grid (sets df.currentRowIdx → fires onCurrentRowChanged)', async () => {
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const newIdx = df.rowCount - 1;
      if (newIdx === df.currentRowIdx) throw new Error(
        `Cannot pick a distinct non-current row: rowCount=${df.rowCount}, currentRowIdx=${df.currentRowIdx}`);
      df.currentRowIdx = newIdx;
    });
  });
  // Step 2.2: viewer re-queries KNN against the newly-current sequence.
  // Scenario .md Scenario 2 Expected: "The viewer's row-card / list panel
  // updates to show neighbours of the newly-clicked row — visibly
  // different cards or a changed top-result identifier compared to the
  // Scenario 1 state".
  //
  // Source-of-truth flow (sequence-similarity-viewer.ts L74-75 +
  // search-base-viewer.ts L39-43):
  //   df.onCurrentRowChanged → DG.debounce 50ms → render(true) →
  //   renderPromise.then(renderInt) → `this.targetMoleculeIdx = df.currentRowIdx`
  //   (SYNCHRONOUS) → await this.computeByMM().
  // The targetMoleculeIdx mutation happens BEFORE the computeByMM await
  // resolves, so polling on targetMoleculeIdx is more responsive than
  // polling on idxs reassignment. 180s tolerance covers the cold-compute
  // ceiling plus the 50ms debounce.
  await softStep('Scenario 2.2: viewer re-queries — targetMoleculeIdx changes from Scenario 1 baseline', async () => {
    await page.waitForFunction((baseline) => {
      const viewers = Array.from((grok.shell.tv as any).viewers) as any[];
      const v = viewers.find((vw) => vw.type === 'Sequence Similarity Search');
      if (!v) return false;
      const targetMoleculeIdx = (v as any)?.targetMoleculeIdx ?? null;
      // Re-query happened if targetMoleculeIdx moved off the Scenario 1
      // baseline. The setter runs synchronously inside renderInt before
      // awaiting computeByMM (sequence-similarity-viewer.ts L75).
      return targetMoleculeIdx !== null
        && targetMoleculeIdx !== (baseline as any).targetMoleculeIdx;
    }, scenario1Baseline, {timeout: 180_000});
    // After the re-query trigger has been observed, wait for the new
    // compute to finish so the viewer state is stable when assertions
    // read it. idxs/scores are reassigned at the end of computeByMM
    // (sequence-similarity-viewer.ts L163-164); we await a fresh KNN
    // population by polling for both columns to remain non-null at the
    // new targetMoleculeIdx.
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
    // Invariant 1: viewer is still present after row-click (re-query MUST
    // NOT crash the viewer — positive-path complement to the GROK-16111
    // empty-input contract in empty-input-row-viewers-spec.ts).
    expect(scenario2Probe.viewerPresent,
      'Viewer MUST remain docked after row-click re-query').toBe(true);
    // Invariant 2: re-query evidence — targetMoleculeIdx moved off the
    // Scenario 1 baseline. The waitForFunction above already gates on
    // this; the assertion here documents the contract for operator
    // review.
    expect(scenario2Probe.targetMoleculeIdx,
      `Re-query MUST surface in viewer state. Scenario 1 targetMoleculeIdx: ` +
      `${scenario1Baseline.targetMoleculeIdx}. Scenario 2: ` +
      `${scenario2Probe.targetMoleculeIdx}. ` +
      `Viewer did not react to df.currentRowIdx change.`)
      .not.toBe(scenario1Baseline.targetMoleculeIdx);
    // Invariant 3: re-queried KNN result still has >= 2 rows (target +
    // at least one neighbour). Mirrors Scenario 1.2-3's contract floor
    // and guards against the silent-empty-result failure mode.
    expect(scenario2Probe.knnRowCount,
      `Re-queried KNN MUST populate >=2 rows. Observed: ${scenario2Probe.knnRowCount}`)
      .not.toBeNull();
    expect(scenario2Probe.knnRowCount!).toBeGreaterThanOrEqual(2);
    // Invariant 4 (scenario .md Scenario 2 Expected): "currentRow
    // indicator in the table tracks the click". The setter mutation has
    // visible spec-layer effect: df.currentRowIdx equals the value we
    // set. apitest-layer ownership of the indicator/canvas-paint is
    // explicitly noted in scenario .md and out of spec scope.
    expect(scenario2Probe.currentRowIdx).toBe(setupProbe.rowCount - 1);
    // df handle is the same TableView dataframe — guards against silent
    // table replacement (parallel to empty-input-row-viewers-spec.ts
    // Invariant 3).
    const dfStable = await page.evaluate(() => !!grok.shell.tv.dataFrame);
    expect(dfStable).toBe(true);
  });
  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
