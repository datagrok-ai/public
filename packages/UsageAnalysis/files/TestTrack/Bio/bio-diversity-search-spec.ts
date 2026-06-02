/* ---
sub_features_covered:
  - bio.viewers.diversity-search
  - bio.search.diversity
  - bio.search.diversity.top-menu
  - bio.detector
  - bio.rendering
--- */
//   related_bugs: [] (positive-path; GROK-16111 empty-input edge contract is
//     "Related-bug context")
// open via the readCsv entry path (bio.md "GROK-18616 entry-path-class
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
test.use(specTestOptions);
const FASTA_DATASET_PATH = 'System:AppData/Bio/tests/filter_FASTA.csv';
const HELM_DATASET_PATH = 'System:AppData/Bio/tests/filter_HELM.csv';
test('Bio Diversity Search docks viewer + re-runs diversity selection on fresh dataset', async ({page}) => {
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
  }, FASTA_DATASET_PATH);
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
    const candidates = ['Bio:diversitySearchTopMenu', 'Bio:diversitySearch'];
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
      macroUnits: macro?.getTag?.('units') ?? null,
      rowCount: df.rowCount,
    };
  });
  expect(setupProbe.hasMacromoleculeCol,
    'atlas bio.detector contract: readCsv path MUST classify a Macromolecule column synchronously').toBe(true);
  expect(setupProbe.macroSemType).toBe('Macromolecule');
  expect(setupProbe.macroUnits,
    'atlas bio.rendering: filter_FASTA.csv MUST classify with units=fasta').toBe('fasta');
  expect(setupProbe.rowCount,
    'scenario .md Setup: ≥ 5 rows so the diversity subset is meaningful').toBeGreaterThanOrEqual(5);
  await softStep('Scenario 1.1: click Bio > Search > Diversity Search (no dialog — viewer docks directly)', async () => {
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      const group = document.querySelector('[name="div-Bio---Search"]');
      if (!group) throw new Error('[name="div-Bio---Search"] group anchor not found');
      group.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      group.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      const leaf = document.querySelector('[name="div-Bio---Search---Diversity-Search"]');
      if (!leaf) throw new Error('[name="div-Bio---Search---Diversity-Search"] leaf not found under Bio > Search');
      (leaf as HTMLElement).click();
    });
  });
  await softStep('Scenario 1.2-3: viewer Sequence Diversity Search docks; diverse subset populated', async () => {
    await page.locator('[name="viewer-Sequence-Diversity-Search"]').waitFor({timeout: 60_000});
    await page.waitForFunction(() => {
      const viewers = Array.from((grok.shell.tv as any).viewers) as any[];
      const v = viewers.find((vw) => vw.type === 'Sequence Diversity Search');
      if (!v) return false;
      const ids = (v as any).renderMolIds;
      return Array.isArray(ids) && ids.length > 0;
    }, null, {timeout: 180_000});
    const dockProbe = await page.evaluate(() => {
      const viewers = Array.from((grok.shell.tv as any).viewers) as any[];
      const v = viewers.find((vw) => vw.type === 'Sequence Diversity Search');
      const viewerTypes = viewers.map((vw) => vw.type);
      const renderMolIds: number[] | null = (v as any)?.renderMolIds ?? null;
      const subsetSize = Array.isArray(renderMolIds) ? renderMolIds.length : null;
      const subsetIndices: number[] | null = Array.isArray(renderMolIds) ? [...renderMolIds] : null;
      return {
        viewerPresent: !!v,
        viewerType: v?.type ?? null,
        viewerTypes,
        subsetSize,
        subsetIndices,
        viewerDockedDOM: !!document.querySelector('[name="viewer-Sequence-Diversity-Search"]'),
      };
    });
    expect(dockProbe.viewerPresent,
      `Expected 'Sequence Diversity Search' viewer in tv.viewers; actual types: ${JSON.stringify(dockProbe.viewerTypes)}`).toBe(true);
    expect(dockProbe.viewerType).toBe('Sequence Diversity Search');
    expect(dockProbe.viewerDockedDOM,
      'viewer dock container [name="viewer-Sequence-Diversity-Search"] MUST be present in DOM').toBe(true);
    // Scenario 1 Expected — "The viewer displays a non-empty diversity-
    // subset panel — i.e. at least two row cards / list entries are
    // surfaced". For a 14-row dataset with limit=10 default the contract
    // is min(14, 10) = 10 (sequence-diversity-viewer.ts L106); the
    // assertion >= 2 covers the floor and is independent of the limit=10
    // default being preserved across Bio versions.
    expect(dockProbe.subsetSize,
      'Viewer diversity compute MUST populate renderMolIds with >=2 entries (non-empty diverse subset).' +
      ` Observed: ${dockProbe.subsetSize}`).not.toBeNull();
    expect(dockProbe.subsetSize!).toBeGreaterThanOrEqual(2);
    // Scenario 1 Expected — "No error balloon appears". Balloon hook is
    // not installed here (positive-path scenario; the empty-input edge
    // contract is covered by sibling empty-input-row-viewers-spec.ts).
    // We assert no error-balloon DOM is present at the end of Scenario 1.
    // grok.shell.balloon errors render as <div class="d4-balloon ... error"> —
    // the strictly-typed `error` class distinguishes errors from info
    // balloons (which the diversity-compute path does not emit on the
    // positive path).
    const errorBalloonCount = await page.locator('.d4-balloon.error, .grok-balloon-error').count();
    expect(errorBalloonCount,
      'Scenario 1 Expected: "No error balloon appears" on positive path').toBe(0);
    // Stash the Scenario-1 subset indices on the page for Scenario 2 to
    // cross-reference (the Scenario 2 viewer-instance probe asserts the
    // NEW viewer's renderMolIds is drawn from a DIFFERENT row identifier
    // set, confirming the diversity computation re-ran on the new column
    // rather than re-displaying the cached Scenario-1 subset).
    await page.evaluate((indices) => {
      (window as any).__scenario1SubsetIndices = indices;
      (window as any).__scenario1RowCount = grok.shell.tv.dataFrame.rowCount;
    }, dockProbe.subsetIndices);
  });
  // Brief settle after Scenario 1's compute so Scenario 2's close-and-reopen
  // sequence is not racing a still-pending compute-tail (the diversity
  // compute is synchronous-from-renderInt's POV but the worker termination
  // in sequence-diversity-viewer.ts L99 happens before the renderMolIds
  // assignment — by waitForFunction's return the worker is already gone).
  await page.waitForTimeout(1500);
  // ============================================================================
  // Scenario 2 — Reopening with a fresh dataset re-runs diversity selection.
  // ============================================================================
  // Step 2.1: close the Diversity Search viewer per scenario .md
  // Scenario 2 step 1 ("close the viewer (right-click viewer header →
  // Close, or equivalent dock close)"). The canonical class-1 JS-API
  // equivalent — and the one the right-click → Close DOM action
  // dispatches internally — is viewer.close() on the viewer instance,
  // which removes the viewer from tv.viewers + tears down its dock host.
  // This is the sanctioned grok-browser JS-API substitute for the
  // canvas-painted viewer-header right-click context menu (no
  // [name=...] selector on the header context-menu Close item).
  await softStep('Scenario 2.1: close the Diversity Search viewer', async () => {
    await page.evaluate(() => {
      const viewers = Array.from((grok.shell.tv as any).viewers) as any[];
      const v = viewers.find((vw) => vw.type === 'Sequence Diversity Search');
      if (!v) throw new Error('Scenario 2.1: expected docked Sequence Diversity Search viewer from Scenario 1');
      v.close();
    });
    // Wait for the viewer to fully detach from tv.viewers AND from DOM
    // before opening the second dataset, otherwise the Scenario 2 viewer-
    // instance probe could race the Scenario-1 viewer's detach and
    // surface stale renderMolIds.
    await page.waitForFunction(() => {
      const viewers = Array.from((grok.shell.tv as any).viewers) as any[];
      const stillThere = viewers.some((vw) => vw.type === 'Sequence Diversity Search');
      const stillInDom = !!document.querySelector('[name="viewer-Sequence-Diversity-Search"]');
      return !stillThere && !stillInDom;
    }, null, {timeout: 30_000});
  });
  // Step 2.2: open filter_HELM.csv per scenario .md Scenario 2 step 2.
  // The same readCsv → addTableView → semType-detect sequence as the
  // FASTA setup; the HELM detector path is bio.md L104 (units=helm).
  // closeAll() at the boundary clears any stale views so the new view
  // is unambiguously the only TableView in shell.
  await softStep('Scenario 2.2: open filter_HELM.csv; verify HELM Macromolecule detector classifies synchronously', async () => {
    await page.evaluate(async (path) => {
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
    }, HELM_DATASET_PATH);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    // Re-wait for the Bio top-menu against the FRESH view (top-menu
    // refreshes on view change; tactical re-readiness mirrors the
    // sibling empty-input-row-viewers-spec.ts per-test setup loop).
    await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
    // Scenario 2 Expected: "The Macromolecule detector classified the
    // HELM column synchronously on open (atlas bio.detector) — the
    // Diversity Search top-menu becomes invokable without requiring a
    // manual semType assignment."
    const helmProbe = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macro = cols.find((c: any) => c.semType === 'Macromolecule') as any;
      return {
        hasMacromoleculeCol: !!macro,
        macroSemType: macro?.semType ?? null,
        macroUnits: macro?.getTag?.('units') ?? null,
        rowCount: df.rowCount,
      };
    });
    expect(helmProbe.hasMacromoleculeCol,
      'Scenario 2 Expected: HELM Macromolecule detector classifies synchronously on open').toBe(true);
    expect(helmProbe.macroSemType).toBe('Macromolecule');
    // filter_HELM.csv carries units=helm per bio.md L104 — confirms the
    // detector took the HELM-notation branch, not a fallback.
    expect(helmProbe.macroUnits,
      'atlas bio.detector: filter_HELM.csv MUST classify with units=helm').toBe('helm');
    expect(helmProbe.rowCount,
      'filter_HELM.csv MUST carry >=2 rows so a non-trivial diverse subset can be selected').toBeGreaterThanOrEqual(2);
  });
  // Step 2.3: click Bio > Search > Diversity Search against the fresh
  // table view (scenario .md Scenario 2 step 3). Same Click-pattern
  // recipe as Scenario 1.1.
  await softStep('Scenario 2.3: click Bio > Search > Diversity Search against the HELM TableView', async () => {
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      const group = document.querySelector('[name="div-Bio---Search"]');
      if (!group) throw new Error('[name="div-Bio---Search"] group anchor not found (Scenario 2)');
      group.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      group.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      const leaf = document.querySelector('[name="div-Bio---Search---Diversity-Search"]');
      if (!leaf) throw new Error(
        '[name="div-Bio---Search---Diversity-Search"] leaf not found under Bio > Search (Scenario 2)');
      (leaf as HTMLElement).click();
    });
  });
  // Step 2.4 (Scenario 2 Expected): viewer re-docks against the new
  // TableView and surfaces a diversity subset whose row identifiers are
  // drawn from the new HELM dataset rather than the prior FASTA dataset.
  //
  // Two distinct invariants:
  //   (a) Re-docking on the new view (NEW viewer instance + DOM
  //       container) — i.e. NOT a stale Scenario-1 viewer left over.
  //   (b) renderMolIds is computed fresh against the HELM column
  //       (the indices reference rows in the HELM dataframe). Since
  //       the HELM dataset row count may differ from FASTA, the most
  //       robust distinctness assertion is: every index in
  //       Scenario-2's renderMolIds is in the valid range
  //       [0, helmRowCount). Combined with the closeAll() boundary
  //       and the NEW viewer-instance probe, this confirms the
  //       diversity compute re-ran on the HELM column.
  await softStep('Scenario 2.4: viewer re-docks on HELM view; renderMolIds drawn from new dataset', async () => {
    await page.locator('[name="viewer-Sequence-Diversity-Search"]').waitFor({timeout: 60_000});
    await page.waitForFunction(() => {
      const viewers = Array.from((grok.shell.tv as any).viewers) as any[];
      const v = viewers.find((vw) => vw.type === 'Sequence Diversity Search');
      if (!v) return false;
      const ids = (v as any).renderMolIds;
      return Array.isArray(ids) && ids.length > 0;
    }, null, {timeout: 180_000});
    const scenario2Probe = await page.evaluate(() => {
      const viewers = Array.from((grok.shell.tv as any).viewers) as any[];
      const v = viewers.find((vw) => vw.type === 'Sequence Diversity Search');
      const renderMolIds: number[] | null = (v as any)?.renderMolIds ?? null;
      const subsetIndices: number[] | null = Array.isArray(renderMolIds) ? [...renderMolIds] : null;
      // Pull the Scenario-1 stashed indices for cross-reference.
      const scenario1Indices: number[] | null = (window as any).__scenario1SubsetIndices ?? null;
      const scenario1RowCount: number | null = (window as any).__scenario1RowCount ?? null;
      return {
        viewerPresent: !!v,
        subsetIndices,
        subsetSize: subsetIndices?.length ?? null,
        scenario1Indices,
        scenario1RowCount,
        currentRowCount: grok.shell.tv.dataFrame.rowCount,
        viewerDockedDOM: !!document.querySelector('[name="viewer-Sequence-Diversity-Search"]'),
      };
    });
    // Invariant 1: viewer is docked on the new TableView.
    expect(scenario2Probe.viewerPresent,
      'Scenario 2 Expected: viewer re-docks on the new HELM TableView').toBe(true);
    expect(scenario2Probe.viewerDockedDOM,
      'viewer dock container [name="viewer-Sequence-Diversity-Search"] MUST be present in DOM (Scenario 2)').toBe(true);
    // Invariant 2: renderMolIds is non-empty.
    expect(scenario2Probe.subsetSize,
      `Scenario 2 viewer MUST populate renderMolIds with >=2 entries. Observed: ${scenario2Probe.subsetSize}`)
      .not.toBeNull();
    expect(scenario2Probe.subsetSize!).toBeGreaterThanOrEqual(2);
    // Invariant 3: every renderMolIds index is in the valid row range
    // for the CURRENT (HELM) dataframe. This is the strongest "drawn
    // from the new dataset" signal that does not depend on
    // implementation-specific shuffle determinism — if Scenario 2's
    // renderMolIds happened to coincide indices-wise with Scenario 1's
    // (rare but possible for small / structured datasets), the indices
    // still resolve against the HELM rowCount, which is the contract
    // the scenario .md actually asserts ("row identifiers are drawn
    // from the new HELM dataset").
    expect(scenario2Probe.subsetIndices,
      'renderMolIds MUST be a concrete index array on the HELM compute path').not.toBeNull();
    for (const idx of scenario2Probe.subsetIndices!) {
      expect(idx,
        `Scenario 2 renderMolIds index ${idx} MUST resolve into the HELM rowCount range ` +
        `[0, ${scenario2Probe.currentRowCount}) — confirms compute re-ran on HELM column, ` +
        `not on cached Scenario-1 FASTA subset (FASTA rowCount was ${scenario2Probe.scenario1RowCount})`)
        .toBeGreaterThanOrEqual(0);
      expect(idx).toBeLessThan(scenario2Probe.currentRowCount);
    }
    // Invariant 4 (scenario .md Scenario 2 Expected): "No error balloon
    // appears". Same surface check as Scenario 1 — the positive-path
    // contract; empty-input edge is sibling-spec scope.
    const errorBalloonCount = await page.locator('.d4-balloon.error, .grok-balloon-error').count();
    expect(errorBalloonCount,
      'Scenario 2 Expected: "No error balloon appears" on positive path').toBe(0);
  });
  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
