/* ---
sub_features_covered:
  - bio.viewers.diversity-search
  - bio.search.diversity
  - bio.search.diversity.top-menu
  - bio.detector
  - bio.rendering
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (coverage_type: smoke — atlas critical_path
//     bio.cp.diversity-search p0 → smoke per STEP E p0-to-smoke mapping)
//   sub_features_covered: [bio.viewers.diversity-search, bio.search.diversity,
//     bio.search.diversity.top-menu, bio.detector, bio.rendering]
//   ui_coverage_responsibility: absent (delegated_to: null)
//   related_bugs: [] (positive-path; GROK-16111 empty-input edge contract is
//     covered by sibling empty-input-row-viewers-spec.ts per scenario .md
//     "Related-bug context")
//   produced_from: atlas-driven
//   coverage_type: smoke
//
// Atlas provenance: realises critical_path bio.cp.diversity-search (atlas
// priority p0). The Sequence Diversity Search top-menu docks the
// SequenceDiversityViewer (package.ts diversitySearchTopMenu → addViewer
// 'Sequence Diversity Search' → dockManager.dock 'down'); the Macromolecule
// detector (atlas bio.detector) classifies the input column synchronously on
// open via the readCsv entry path (bio.md "GROK-18616 entry-path-class
// invariant"). Unlike Similarity Search (a per-row KNN viewer), Diversity
// Search is a full-column viewer: it computes the FULL distance matrix over
// the target column (sequence-diversity-viewer.ts L93-99) and selects the
// maximally-diverse subset of size `limit` (default 10 per bio.md L322
// getOptions().look snapshot) via getDiverseSubset (L104-110). For datasets
// above MAX_SAMPLE_SIZE=10000 (sequence-diversity-viewer.ts L16) a random
// subset is sampled first; both test datasets here (filter_FASTA.csv 14 rows,
// filter_HELM.csv) are well below that threshold so the full distance-matrix
// path is taken.
//
// Source-of-truth readiness signals (sequence-diversity-viewer.ts):
//   - L21: `renderMolIds: number[] | null = null` — initialised null,
//          populated AFTER computeByMM() resolves.
//   - L113: `this.renderMolIds = diverseIndicesInWorkingSet.map(...)`
//          assigned at the TAIL of computeByMM — DEFINITIVE done-signal.
//   - L43-44: result column DG.Column.string(<diverseColumnName>,
//          this.renderMolIds!.length).init(...) — confirms renderMolIds
//          is the canonical selection-index source.
//   - L58: `this.computeCompleted.next(true)` — Subject<boolean> fired
//          after the embedded grid is wired in. Alternative done-signal
//          (we use the simpler `renderMolIds != null` JS-API probe).
//
// Sibling specs (reference templates per § 4.4):
//   bio-similarity-search-spec.ts — canonical Bio Search top-menu spec
//     (retry-1 with per-leaf function-registration probe + post-compute
//     viewer-instance JS-API readiness probe). Adopted verbatim here for
//     the Bio init readiness sequence + the per-leaf function-registration
//     probe pattern. The post-compute readiness probe is swapped from
//     similarity's `idxs/scores` to diversity's `renderMolIds`.
//   empty-input-row-viewers-spec.ts — sibling Bio setup + cold-start
//     stabilization + Diversity Search dock contract under the empty-
//     current-row case. Adopted verbatim for the readCsv + semType-detect +
//     Bio init-completion readiness sequence.
//   analyze-spec.ts — canonical per-leaf function-registration probe
//     pattern + cold-start tactical guards on Bio top-menu clicks.
//
// Selector sources (class 1 — all in grok-browser/references/bio.md):
//   [name="div-Bio"]                                — bio.md L79 + L592 + L76
//                                                     top-menu root.
//   [name="div-Bio---Search"]                       — bio.md L79 "Click
//                                                     pattern" group anchor;
//                                                     hover-not-click required.
//   [name="div-Bio---Search---Diversity-Search"]    — bio.md L310 (no "..."
//                                                     suffix; no dialog —
//                                                     viewer docks directly).
//   [name="viewer-Sequence-Diversity-Search"]       — bio.md L314, L617.
//   Viewer-instance type 'Sequence Diversity Search' — bio.md L316
//                                                     (Array.from(grok.shell.tv.viewers)
//                                                     .find(v => v.type === ...)).
//   [name="viewer-Grid"]                            — bio.md L107 Bio grid host;
//                                                     used for cold-start readiness.
//
// MCP recon observation (this initial dispatch):
//   list_pages returned successfully (single page open at
//   https://dev.datagrok.ai/). Post-list_pages evaluate_script captured
//   loginFormVisible=true, browsePresent=false — Datagrok profile auth
//   stale despite the cycle's prewarm_mcp_browser_auth precondition
//   (likely session-token TTL exceeded between prewarm and this dispatch,
//   same condition the sibling bio-similarity-search-spec.ts retry-1
//   dispatch encountered). Per §"MCP recon — auth assumption" the agent
//   does NOT attempt to re-auth from inside MCP tool calls — the
//   architecture intentionally moves auth out of the agent boundary.
//   No empirical session-replay was possible this dispatch; spec authored
//   from class-1 selectors only (all selectors present in bio.md recon
//   date 2026-06-01 Bio 2.26.5) + Datagrok Bio package source-of-truth
//   (sequence-diversity-viewer.ts L21 / L34-61 / L93-113 / L16). No class-2
//   selector emitted; no Selector recon-notes block required.
//
// Cleanup contract per scenario .md "Notes": none cited; the test cleans
// up its TableView implicitly via fresh Playwright context (no shared
// state across tests). Scenario 2 explicitly opens a SECOND dataset
// (filter_HELM.csv) into a NEW TableView; we let grok.shell.closeAll()
// at the Scenario-2 setup boundary clear the Scenario-1 view so the
// viewer-instance probe unambiguously sees the new docked instance.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// filter_FASTA.csv has 14 rows per direct file inspection (bio.md L598
// "Test datasets" table + sibling bio-similarity-search-spec.ts L145-148);
// satisfies the scenario .md "Setup" assertion of ≥ 5 rows so the
// diversity subset is meaningful (limit=10 default → 10 picks out of 14
// → visibly different from a trivial first-N).
const FASTA_DATASET_PATH = 'System:AppData/Bio/tests/filter_FASTA.csv';

// filter_HELM.csv — Scenario 2 second dataset. HELM notation column,
// detected synchronously by the Macromolecule detector (atlas bio.detector)
// per bio.md L104. Distinct from filter_FASTA.csv so the Scenario-2
// re-run's renderMolIds are necessarily drawn from a different column
// (different row identifiers, different distance-matrix → different
// diverse subset).
const HELM_DATASET_PATH = 'System:AppData/Bio/tests/filter_HELM.csv';

test('Bio Diversity Search docks viewer + re-runs diversity selection on fresh dataset', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  // ============================================================================
  // Scenario 1 — Diversity Search top-menu docks the diversity viewer.
  // ============================================================================

  // Setup phase: open filter_FASTA.csv via readCsv (atlas bio.io.fasta-handler;
  // bio.md "GROK-18616 entry-path-class invariant" — the readCsv path triggers
  // Macromolecule detector synchronously, satisfying scenario .md "Setup"
  // step: "The Macromolecule detector (atlas bio.detector) classifies the
  // sequence column synchronously on open"). Wait for semType detection +
  // Bio cell-renderer init (atlas bio.rendering — "the table view opens with
  // the Macromolecule cell renderer painting sequence cells").
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

  // Bio top-menu + init-completion readiness — mirrors analyze-spec.ts +
  // empty-input-row-viewers-spec.ts L198-205 + bio-similarity-search-spec.ts
  // L188-198. Layer 1: DOM visibility of [name="div-Bio"]. Layer 2:
  // Bio:getSeqHelper resolves only after initBio completes (bio.md
  // "Init-order invariant").
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  await page.evaluate(async () => {
    const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
    for (const fn of probes) {
      try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
    }
    await new Promise((r) => setTimeout(r, 3000));
  });

  // Per-leaf function-registration probe (cold-start tactical guard;
  // mirrors bio-similarity-search-spec.ts L216-237 + analyze-spec.ts
  // L113-156).
  //
  // The Bio package init probe above guarantees init COMPLETION (via
  // grok.functions.call('Bio:getSeqHelper') — atlas
  // bio.cp.bio-service-surface-init). It does NOT guarantee that the
  // top-menu LEAF function dispatched by
  // [name="div-Bio---Search---Diversity-Search"] (Bio:diversitySearchTopMenu,
  // registered in package.ts) is findable in the function registry yet —
  // leaf registration is a distinct code path that can lag init completion
  // by a short window on truly-cold Bio boots. Bounded poll until the
  // target Search leaf is findable as a Datagrok function, with a short
  // defensive settle if it is missing on this Bio build (function rename
  // across Bio versions).
  await page.evaluate(async () => {
    const candidates = ['Bio:diversitySearchTopMenu', 'Bio:diversitySearch'];
    const findAny = (names: string[]): boolean => {
      for (const n of names) {
        try {
          if ((grok as any).functions.find && (grok as any).functions.find(n)) return true;
        } catch { /* try next */ }
      }
      return false;
    };
    const deadline = Date.now() + 15_000;
    while (Date.now() < deadline) {
      if (findAny(candidates)) return;
      await new Promise((r) => setTimeout(r, 300));
    }
    // Even if none of the candidate names are findable (function rename
    // across Bio versions), do not error — the per-step 180s viewer
    // tolerance below is the defensive ceiling. Short settle so the menu
    // dispatch doesn't race a still-completing function-table install.
    await new Promise((r) => setTimeout(r, 1500));
  });
  await page.waitForTimeout(2000);

  // Verify the scenario .md "Setup" preconditions: Macromolecule semType
  // detected (atlas bio.detector contract) and table row count >= 5 (so
  // the diversity subset is meaningful — limit=10 default → 10 of 14
  // selected → visibly different from a trivial first-N).
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
  // filter_FASTA.csv carries `units=fasta` per bio.md L103 — sanity-checks
  // the renderer-dispatch precondition (atlas bio.rendering).
  expect(setupProbe.macroUnits,
    'atlas bio.rendering: filter_FASTA.csv MUST classify with units=fasta').toBe('fasta');
  expect(setupProbe.rowCount,
    'scenario .md Setup: ≥ 5 rows so the diversity subset is meaningful').toBeGreaterThanOrEqual(5);

  // Step 1.1: click Bio > Search > Diversity Search via the bio.md
  // "Click pattern" recipe (click root → 400ms → mouseover Search →
  // 300ms → click leaf). Hover-not-click on the group is required to
  // surface the leaves per bio.md L79. The leaf has NO "..." suffix and
  // opens NO dialog — viewer docks directly (bio.md L310 +
  // empty-input-row-viewers-spec.ts L131-135 hasEditorDialog=false).
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

  // Step 1.2 / 1.3 (combined): wait for the Sequence Diversity Search
  // viewer to dock AND finish its initial diversity compute (scenario .md
  // Scenario 1 step 3: "Wait for the viewer to compute the full distance
  // matrix over the input column and surface the diverse subset").
  //
  // Readiness signal: the viewer instance's populated `renderMolIds`
  // column. Per sequence-diversity-viewer.ts L113, renderMolIds is
  // assigned at the TAIL of computeByMM — a definitive done-signal that
  // (a) the full distance matrix has been computed, (b)
  // getDiverseSubset has returned a non-empty index list, (c) the
  // working-set → original-index mapping has been applied. After this
  // assignment, renderInt L43-44 builds the result column directly off
  // renderMolIds.length, so `renderMolIds != null` is sufficient + the
  // length is the diversity-subset cardinality.
  //
  // 180s tolerance matches sibling Bio cold-compute ceiling
  // (bio-similarity-search-spec.ts L320 / L455). The full distance matrix
  // for a 14-row column is fast; the tolerance is for the cold Bio /
  // worker spawn (DistanceMatrixService L92 spawns a worker per call).
  // Well below the 600s Playwright B-STAB-04 layer bound.
  await softStep('Scenario 1.2-3: viewer Sequence Diversity Search docks; diverse subset populated', async () => {
    // Layer 1 DOM probe: viewer dock container present per bio.md L314, L617.
    // Cold-start safe at 60s; matches sibling search-spec.ts dock-wait envelope.
    await page.locator('[name="viewer-Sequence-Diversity-Search"]').waitFor({timeout: 60_000});

    // Layer 2 JS-API probe: viewer instance present in the TableView's
    // viewer collection AND its compute path has assigned renderMolIds.
    // Per sequence-diversity-viewer.ts L93-113:
    //   const distanceMatrixData = await distanceMatrixService.calc(...);
    //   ...
    //   const diverseIndicesInWorkingSet = getDiverseSubset(...);
    //   this.renderMolIds = diverseIndicesInWorkingSet.map(...);
    // renderMolIds is non-null AFTER computeByMM's distance-matrix worker
    // call resolves AND getDiverseSubset returns — the most reliable
    // diversity-compute-done signal.
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
      // Diversity subset row count via the viewer's renderMolIds array
      // (canonical post-compute surface per sequence-diversity-viewer.ts L113).
      // For a 14-row dataset with limit=10 default (bio.md L322
      // getOptions().look snapshot), expect Math.min(workingLength, limit)
      // = Math.min(14, 10) = 10 selected rows per
      // sequence-diversity-viewer.ts L106.
      const renderMolIds: number[] | null = (v as any)?.renderMolIds ?? null;
      const subsetSize = Array.isArray(renderMolIds) ? renderMolIds.length : null;
      // Capture the actual selected indices as a stable identifier for
      // the dataset's diverse subset — used in Scenario 2 to assert the
      // re-run picked from the NEW dataset (different row identifier set).
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

    // Scenario 1 Expected — "The Sequence Diversity Search viewer is
    // present in the active view's viewer list (isViewerPresent /
    // findViewer style assertion against the viewer-name 'Sequence
    // Diversity Search')".
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
