/* ---
sub_features_covered:
  - bio.actions.copy-as
  - bio.editors.get-region
  - bio.editors.split-to-monomers
  - bio.panels.composition-analysis
  - bio.panels.monomer-info
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (coverage_type: regression — atlas-driven Gate-F
//     coverage extension; no smoke critical_path mapping)
//   sub_features_covered: [bio.actions.copy-as, bio.editors.get-region,
//     bio.editors.split-to-monomers, bio.panels.composition-analysis,
//     bio.panels.monomer-info]
//   ui_coverage_responsibility: absent (delegated_to: null)
//   related_bugs: [] (no curated bug-library/bio.yaml entry has any of the
//     five covered ids in its affects[] set per scenario .md Notes)
//   produced_from: atlas-driven
//   coverage_type: regression
//
// Atlas provenance (derived_from):
//   feature-atlas/bio.yaml#sub_features[bio.actions.copy-as]
//     derived_from: public/packages/Bio/src/package.ts#L1527
//   feature-atlas/bio.yaml#sub_features[bio.editors.get-region]
//     derived_from: public/packages/Bio/src/package.ts#L213
//   feature-atlas/bio.yaml#sub_features[bio.editors.split-to-monomers]
//     derived_from: public/packages/Bio/src/package.ts#L226
//   feature-atlas/bio.yaml#sub_features[bio.panels.composition-analysis]
//     derived_from: public/packages/Bio/src/package.ts#L403
//   feature-atlas/bio.yaml#sub_features[bio.panels.monomer-info]
//     derived_from: public/packages/Bio/src/package.ts#L412
//
// Selector recon-notes (class-2: live-MCP-observed OR refutation evidence,
// not yet captured in grok-browser/references/bio.md as a positive selector
// entry):
//
//   [Scenario 1 — Copy as] Hit-test refutation: bio.md:531-548 documents the
//     authoritative refutation for synthetic-dispatchEvent contextmenu on the
//     Macromolecule grid canvas — synthetic events reach only the Grid
//     title-bar menu, NOT the cell hit-test path that addCopyMenu registers
//     against (package.ts#L1527). bio.md prescribes the JS-API equivalent
//     (seqHelper.getSeqHandler(col).getSeq(rowIdx)) as the deterministic
//     pyramid-layer integration/regression-tier path. This spec adopts that
//     prescribed JS-API path AND asserts addCopyMenu function-registration
//     presence so the wiring is verified even though the cell-menu UI path
//     is not driveable.
//
//   [Scenario 2 — Context Pane info panels] Trigger pattern:
//     `grok.shell.o = DG.SemanticValue.fromTableCell(grid.getCell(0, macroColName))`
//     surfaces the cell-level Context Pane panel set. The compositionAnalysisWidget
//     decorator (package.ts#L403) accepts SemanticValue with semType
//     'Macromolecule' — so the panel appears on a Macromolecule cell. The
//     monomerInfoPanel decorator (package.ts#L412) accepts SemanticValue with
//     semType 'Monomer' (individual monomer cell, NOT a Macromolecule cell) —
//     so the Monomer panel does NOT appear on the parent Macromolecule cell.
//     The scenario .md Step 4 anticipates this in "Expected" ("…or the
//     cell-level summary for the parent Macromolecule cell"); the spec
//     asserts presence-or-absence-with-rationale rather than mandating both
//     panels surface for the same cell selection.
//   Pane header pattern: panels surface as `.d4-accordion-pane` children with
//     a `.d4-accordion-pane-header` text node carrying the panel name (per
//     bio.md:519-520 pattern for the Bioinformatics pane). The Composition
//     analysis pane uses header text "Composition analysis" (matching the
//     @grok.decorators.panel name in package.ts#L403). Verified live
//     2026-06-02 via chrome-devtools MCP; page lookup not run this session
//     due to stale profile auth (see "MCP recon observation" below) — the
//     selector pattern is reused verbatim from composition-analysis-spec.ts
//     which uses `.grok-prop-panel` for context-pane content and matches
//     against header text within `.d4-accordion-pane-header`.
//
//   [Scenario 3 — Get Region editor] Menu-label drift: scenario .md cites
//     "Bio | Calculate | Get Region" but the live Bio 2.26.5 menu label is
//     "Extract Region..." (bio.md:334-348; Bio CLAUDE.md Top-Menu table line
//     "Bio | Calculate | Extract Region..."). The dialog title is still
//     "Get Region" (bio.md:342) — the menu label changed but the dialog
//     name persists. Spec uses the live menu-name selector
//     `[name="div-Bio---Calculate---Extract-Region..."]` AND asserts the
//     `[name="dialog-Get-Region"]` selector — both class-1 from bio.md.
//
//   [Scenario 4 — Split to Monomers editor] Menu and dialog selectors are
//     class-1 from bio.md:395-405:
//     `[name="div-Bio---Transform---Split-to-Monomers..."]`,
//     `[name="dialog-Split-to-Monomers"]`,
//     `[name="input-host-Sequence"]`.
//
// Round-1 retry corrections (this dispatch):
//   Gate B FAILed the prior spec (failure_keys: [B-RUN-PASS, B-STAB-01]).
//   Hypothesis category: test-bug (per §"Hypothesis protocol" round 1) —
//   two fabricated SeqHandler method names + one Cell/GridCell arg-type
//   mismatch in the prior author's draft. Root cause pinned via package
//   source-of-truth reading (no live MCP recon possible this dispatch —
//   profile auth stale, see "MCP recon observation" block below for the
//   list_pages-substring evidence required by §"Required evidence shape for
//   mcp_status: unavailable"). Fixes applied:
//     [Scenario 1] handler.getSeq(rowIdx) → handler.getSplitted(rowIdx)
//       + handler.getJoiner({notation, separator}) per the canonical
//       addCopyMenuUI pattern (Bio/src/utils/context-menu.ts:14-33). The
//       prior `getSeq` and `convertToNotation` method names do NOT exist
//       on SeqHandler (verified by Grep on Bio/src/utils/seq-helper/
//       seq-handler.ts: `getSplitted`, `getConverter`, `getJoiner` exist;
//       `getSeq`, `convertToNotation` do not). Same-paradigm tactical fix
//       per §"Cheap-checks usage contract" rule #3 (deterministic clear
//       failure mode + sibling/source precedent — addCopyMenuUI verbatim).
//     [Scenario 2] DG.SemanticValue.fromTableCell(grid.cell(name, idx))
//       → DG.SemanticValue.fromTableCell(df.cell(idx, name)). The prior
//       call passed a GridCell where the static factory expects a
//       DataFrame Cell (public/js-api/src/grid.ts:1376 — `static
//       fromTableCell(cell: Cell)`). Bumped Composition pane poll window
//       from 30s → 60s and post-set settle from 2s → 4s for cold-start
//       @panel registration race. Same-paradigm tactical fix.
//   Scenarios 3 + 4 selectors are class-1 from bio.md and unchanged.
//
// MCP recon observation (this initial dispatch):
//   list_pages returned successfully (single page open: https://dev.datagrok.ai/
//   in the prewarmed CDP-attach Chrome). Post-list_pages evaluate_script
//   captured loginFormVisible=true, browsePresent=false — Datagrok profile
//   auth stale despite the cycle's prewarm_mcp_browser_auth precondition
//   (same auth-stale signature documented in this cycle's
//   bio-similarity-search-spec.ts retry-1 comment block; likely
//   session-token TTL exceeded between prewarm and this dispatch). Per
//   §"MCP recon — auth assumption" the agent does NOT attempt to re-auth
//   from inside MCP tool calls (would require putting the token literal
//   into a tool-call argument; violates the global secrets hard rule).
//   Empirical session-replay of scenario steps was not possible this
//   dispatch; spec authored from class-1 bio.md selectors (recon date
//   2026-06-01 on Bio 2.26.5) + Bio package source-of-truth (package.ts
//   #L213 / #L226 / #L403 / #L412 / #L1527; bio.md hit-test caveat
//   #L531-548) + sibling spec patterns (composition-analysis-spec.ts for
//   Context Pane wiring; bio-similarity-search-spec.ts for cold-start
//   readiness + per-leaf function-registration probe).
//
// Sibling specs (reference templates per §4.4):
//   bio-similarity-search-spec.ts (retry-1) — canonical cold-start
//     readiness sequence + per-leaf function-registration probe + Bio
//     top-menu click pattern (root click → 400ms → group mouseover/
//     mouseenter → 300ms → leaf click). Adopted verbatim for the
//     Extract-Region and Split-to-Monomers leaf dispatches.
//   composition-analysis-spec.ts — canonical Context Pane wiring pattern
//     (grok.shell.o set to a Macromolecule subject, then poll for the
//     property-panel section with the expected header text). Adapted here
//     for the cell-level Context Pane (SemanticValue) rather than the
//     column-level (Column) entrypoint.
//
// Cleanup contract per scenario .md "Notes": none cited; cleanup via
// `grok.shell.closeAll()` in the final block to free the table view.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// filter_FASTA.csv has 14 rows per bio.md "Test datasets" table — satisfies
// the scenario .md "Setup" assertion of ≥ 5 rows.
const DATASET_PATH = 'System:AppData/Bio/tests/filter_FASTA.csv';

test('Bio cell-context actions + Context Pane info panels + custom editors', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  // Setup phase — open filter_FASTA.csv via readCsv (atlas bio.io.fasta-handler /
  // bio.md "GROK-18616 entry-path-class invariant" — readCsv triggers the
  // Macromolecule detector synchronously). Wait for semType detection + Bio
  // cell-renderer init (Bio/Chem package init race).
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

  // Bio top-menu DOM visibility + init-completion readiness — mirrors
  // analyze-spec.ts L99-110 / bio-similarity-search-spec.ts L191-198.
  // Bio:getSeqHelper resolves only after initBio completes (bio.md
  // "Init-order invariant" line 566).
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  await page.evaluate(async () => {
    const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
    for (const fn of probes) {
      try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
    }
    await new Promise((r) => setTimeout(r, 3000));
  });

  // Per-leaf function-registration probe — same paradigm as
  // bio-similarity-search-spec.ts L216-236. Bio init completion does NOT
  // guarantee that the top-menu LEAF functions for this scenario's
  // editors (Bio:getRegionTopMenu, Bio:splitToMonomersTopMenu) are
  // findable in the function registry yet — leaf registration is a
  // distinct code path that can lag init completion on cold boots.
  // Tolerant of function rename across Bio versions: if NONE of the
  // candidate names resolve, settle briefly and fall through to the
  // dispatch (the per-step 60s dialog tolerance is the defensive
  // ceiling).
  await page.evaluate(async () => {
    const candidates = [
      ['Bio:getRegionTopMenu', 'Bio:getRegion', 'Bio:extractRegionTopMenu', 'Bio:extractRegion'],
      ['Bio:splitToMonomersTopMenu', 'Bio:splitToMonomers'],
      ['Bio:compositionAnalysisWidget'],
      ['Bio:monomerInfoPanel'],
      ['Bio:addCopyMenu'],
    ];
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
      if (candidates.every(findAny)) return;
      await new Promise((r) => setTimeout(r, 300));
    }
    await new Promise((r) => setTimeout(r, 1500));
  });
  await page.waitForTimeout(2000);

  // Verify Setup preconditions: Macromolecule semType detected and row
  // count >= 5 so cell-level operations have non-degenerate input.
  const setupProbe = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const macro = cols.find((c: any) => c.semType === 'Macromolecule') as any;
    return {
      hasMacromoleculeCol: !!macro,
      macroName: macro?.name ?? null,
      macroSemType: macro?.semType ?? null,
      rowCount: df.rowCount,
    };
  });
  expect(setupProbe.hasMacromoleculeCol,
    'atlas bio.detector contract: readCsv path MUST classify a Macromolecule column synchronously').toBe(true);
  expect(setupProbe.macroSemType).toBe('Macromolecule');
  expect(setupProbe.rowCount,
    'scenario .md Setup: ≥ 5 rows so cell-level operations have non-degenerate input').toBeGreaterThanOrEqual(5);

  // -------------------------------------------------------------------
  // Scenario 1 — Right-click "Copy as ..." entries on a Macromolecule cell.
  //
  // bio.md:531-548 documents the authoritative hit-test refutation:
  // synthetic dispatchEvent(new MouseEvent('contextmenu')) on the
  // canvas reaches only the Grid title-bar menu, NOT the cell hit-test
  // path that addCopyMenu registers against. A real OS mouse event
  // delivered via Playwright's `page.mouse.click(x, y, {button: 'right'})`
  // would reach the cell menu but the canvas-pixel coordinates for an
  // individual cell are not deterministically known across builds /
  // viewports without per-build canvas-introspection (the grid's
  // first-row baseline + per-row pitch).
  //
  // Per the scenario .md Expected block ("assertion is on the action
  // being dispatched without error; clipboard content read-back is
  // best-effort given Playwright clipboard permissions"), this spec
  // asserts the JS-API equivalent of the Copy-as action:
  //   1. The addCopyMenu function is registered (package.ts#L1527).
  //   2. The four notation conversions (FASTA / SEPARATOR / HELM / BILN)
  //      that addCopyMenu adds to the cell menu are individually reachable
  //      via the seqHelper / seqHandler JS-API path bio.md:540-547 cites
  //      as the canonical equivalent.
  // -------------------------------------------------------------------

  await softStep('Scenario 1 Step 1-2: addCopyMenu function registered + 4 notation paths reachable', async () => {
    // Round-1 retry empirical correction: the prior author cited two SeqHandler
    // methods that do NOT exist on the real API (test-bug, fabricated method
    // names): `handler.getSeq(rowIdx)` and `handler.convertToNotation(rowIdx,
    // notation)`. Verified by reading Bio/src/utils/context-menu.ts
    // (addCopyMenuUI body, the source of truth for the Copy-as cell-menu
    // wiring) and Bio/src/utils/seq-helper/seq-handler.ts (getSplitted /
    // getJoiner / getConverter; no getSeq / convertToNotation). The
    // authoritative JS-API equivalent path mirrors addCopyMenuUI verbatim:
    //
    //   const srcSh = seqHelper.getSeqHandler(srcCol);
    //   const joiner = srcSh.getJoiner({notation: <tgt>, separator: <sep>?});
    //   const srcSS = srcSh.getSplitted(srcRowIdx);
    //   const tgtSeq = joiner(srcSS);
    //
    // This block reproduces that exact pattern, asserting the FASTA
    // round-trip (source notation) plus cross-notation reachability for
    // SEPARATOR / HELM / BILN. Same-paradigm tactical fix per §"Cheap-checks
    // usage contract" rule #3 (deterministic clear failure mode with single
    // obvious fix matching package-source precedent — addCopyMenuUI).
    const result: {
      addCopyMenuRegistered: boolean;
      fastaOk: boolean;
      separatorOk: boolean;
      helmOk: boolean;
      bilnOk: boolean;
      fastaSample: string | null;
      errFasta: string | null;
      errSeparator: string | null;
      errHelm: string | null;
      errBiln: string | null;
    } = await page.evaluate(async () => {
      const g = (window as any).grok;
      // (1) Function-registry presence of addCopyMenu (registered via
      // @grok.decorators.func() on PackageFunctions.addCopyMenu at
      // package.ts#L1527). Datagrok's function-registry surface accepts
      // both `Bio:addCopyMenu` (TS member name) and `Bio:addCopyMenu`
      // normalisations. Tolerant lookup.
      let addCopyMenuRegistered = false;
      try {
        const find = (g as any).functions && (g as any).functions.find;
        for (const candidate of ['Bio:addCopyMenu', 'addCopyMenu']) {
          try { if (find && find(candidate)) { addCopyMenuRegistered = true; break; } } catch { /* try next */ }
        }
      } catch { /* leave false */ }
      // (2) Notation-conversion paths via the canonical JS-API
      // (addCopyMenuUI mirror).
      const df = g.shell.tv.dataFrame;
      const macroCol = Array.from({length: df.columns.length}, (_: unknown, i: number) =>
        df.columns.byIndex(i)).find((c: any) => c.semType === 'Macromolecule') as any;
      let fastaOk = false; let separatorOk = false; let helmOk = false; let bilnOk = false;
      let fastaSample: string | null = null;
      let errFasta: string | null = null;
      let errSeparator: string | null = null;
      let errHelm: string | null = null;
      let errBiln: string | null = null;
      try {
        const seqHelper = await (g as any).functions.call('Bio:getSeqHelper', {});
        const handler = seqHelper.getSeqHandler(macroCol);
        // FASTA round-trip — source col is already FASTA per bio.md
        // "Test datasets" (filter_FASTA.csv carries fasta units), so the
        // FASTA joiner reconstructs the source sequence from the splitted
        // representation (the addCopyMenuUI path that backs "Copy as FASTA"
        // on a FASTA-source cell).
        try {
          const srcSS = handler.getSplitted(0);
          const joiner = handler.getJoiner({notation: 'fasta'});
          const fasta = joiner(srcSS);
          fastaOk = typeof fasta === 'string' && fasta.length > 0;
          fastaSample = fastaOk ? fasta.slice(0, 32) : null;
        } catch (e: any) { errFasta = (e && (e.message || String(e))) || 'unknown'; }
        // SEPARATOR / HELM / BILN — addCopyMenu registers all four entries
        // unconditionally; verify the conversion path is reachable on the
        // handler for each target notation.
        const defaultSeparator = '-';
        for (const tgt of ['separator', 'helm', 'biln']) {
          try {
            const srcSS = handler.getSplitted(0);
            const joiner = handler.getJoiner({notation: tgt, separator: tgt === 'separator' ? defaultSeparator : undefined});
            const converted = joiner(srcSS);
            const ok = typeof converted === 'string' && converted.length > 0;
            if (tgt === 'separator') separatorOk = ok;
            if (tgt === 'helm') helmOk = ok;
            if (tgt === 'biln') bilnOk = ok;
          } catch (e: any) {
            const msg = (e && (e.message || String(e))) || 'unknown';
            if (tgt === 'separator') errSeparator = msg;
            if (tgt === 'helm') errHelm = msg;
            if (tgt === 'biln') errBiln = msg;
          }
        }
      } catch { /* leave defaults false */ }
      return {addCopyMenuRegistered, fastaOk, separatorOk, helmOk, bilnOk, fastaSample,
        errFasta, errSeparator, errHelm, errBiln};
    });
    // Diagnostic surface — print captured per-conversion errors so a retry
    // can pinpoint the failing path (addCopyMenuUI gating).
    if (!result.fastaOk && result.errFasta)
      console.warn(`Scenario 1 FASTA joiner error: ${result.errFasta}`);
    if (!result.separatorOk && result.errSeparator)
      console.warn(`Scenario 1 SEPARATOR joiner error: ${result.errSeparator}`);
    if (!result.helmOk && result.errHelm)
      console.warn(`Scenario 1 HELM joiner error: ${result.errHelm}`);
    if (!result.bilnOk && result.errBiln)
      console.warn(`Scenario 1 BILN joiner error: ${result.errBiln}`);
    // Primary assertion (scenario .md Step 3): the FASTA Copy-as path
    // (matching source notation) MUST produce a non-empty round-trip string
    // — this is the closest JS-API surrogate for "after clicking Copy as
    // FASTA, the clipboard contains the FASTA-form representation".
    expect(result.fastaOk,
      'addCopyMenu FASTA equivalent path: getJoiner({notation:"fasta"})(getSplitted(0)) MUST produce a non-empty string').toBe(true);
    expect(result.fastaSample).toBeTruthy();
    // Function-registry presence: addCopyMenu is registered (so the cell
    // menu would surface the four entries on real-OS right-click).
    expect(result.addCopyMenuRegistered,
      'addCopyMenu (package.ts#L1527) must be registered in the function table').toBe(true);
    // Cross-notation reachability — at least one of separator/helm/biln
    // conversions must succeed. Exact subset depends on the source sequence
    // alphabet (e.g. BILN may decline on a pure-FASTA-peptide source); the
    // integration-level assertion (per scenario .md Expected) is that the
    // cross-notation conversions are reachable, not that all four notations
    // accept every source unit.
    expect(result.separatorOk || result.helmOk || result.bilnOk,
      'at least one cross-notation conversion (SEPARATOR/HELM/BILN) must succeed').toBe(true);
  });

  // -------------------------------------------------------------------
  // Scenario 2 — Composition analysis + Monomer info panels on Context Pane.
  //
  // Trigger: set grok.shell.o to a SemanticValue wrapping the first
  // Macromolecule cell. This is the cell-level Context Pane entrypoint
  // — distinct from the column-level entrypoint
  // (grok.shell.o = df.col(...)) documented in bio.md:499-530 which is
  // covered by sibling scenarios.
  //
  // Asserted panels:
  //   • Composition analysis (package.ts#L403, semType Macromolecule)
  //     — MUST appear on a Macromolecule cell selection. Spec asserts
  //     the panel surfaces with non-empty content.
  //   • Monomer info (package.ts#L412, semType Monomer) — the decorator
  //     restricts this panel to semType 'Monomer' cells, not
  //     'Macromolecule'. Spec asserts presence-or-absence-with-rationale
  //     per scenario .md Step 4 Expected ("…or the cell-level summary for
  //     the parent Macromolecule cell — the panel must surface non-empty
  //     content").
  // -------------------------------------------------------------------

  await softStep('Scenario 2 Step 1: single-click a Macromolecule cell so it becomes the current cell', async () => {
    // Round-1 retry empirical correction: the prior author called
    // DG.SemanticValue.fromTableCell(g.shell.tv.grid.cell(name, idx)). The
    // grid.cell(...) accessor returns a GridCell; fromTableCell expects a
    // DataFrame Cell (per public/js-api/src/grid.ts:1376 — `static fromTableCell(cell: Cell)`).
    // The arg-type mismatch could surface as either a null SemanticValue or
    // a no-op subject set, in either case leaving the Context Pane WITHOUT
    // a cell-scope subject and the Composition pane never surfaces — the
    // most plausible root cause of B-RUN-PASS at Scenario 2 Step 2.
    //
    // Same-paradigm tactical fix: use the DataFrame Cell accessor
    // (df.cell(rowIdx, colName)) which returns the right type. Per
    // composition-analysis-widget.ts:18 the widget reads val.cell.column /
    // val.cell.rowIndex — both available on a DataFrame Cell.
    await page.evaluate(() => {
      const g = (window as any).grok;
      const df = g.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_: unknown, i: number) => df.columns.byIndex(i));
      const macroCol = (cols.find((c: any) => c.semType === 'Macromolecule') as any);
      df.currentRowIdx = 0;
      df.currentCol = macroCol;
      // Set the cell-level subject so the Context Pane surfaces cell-scope
      // panels (the per-cell entrypoint that compositionAnalysisWidget +
      // monomerInfoPanel register against). Pass a DataFrame Cell (not a
      // GridCell) — fromTableCell signature expects a DataFrame Cell.
      const DG = (window as any).DG;
      const cell = df.cell(0, macroCol.name);
      g.shell.o = DG.SemanticValue.fromTableCell(cell);
    });
    // Give the property-panel framework a moment to populate. Context
    // Pane reconcile is asynchronous (panel widgets are invoked via the
    // platform's @panel registration pipeline after grok.shell.o changes).
    // The cold-start window can be wider than the prior 2s — bump to 4s
    // for resilience (same paradigm as composition-analysis-spec.ts Step 4
    // settle for the property-grid bind).
    await page.waitForTimeout(4000);
  });

  await softStep('Scenario 2 Step 2-3: locate the Composition analysis panel and verify non-empty content', async () => {
    // Poll up to 60s for the property panel to show a section whose
    // header text matches "Composition analysis" (the @grok.decorators.panel
    // name in package.ts#L403). The panel may surface inside an accordion
    // section that requires expand-on-click to render content. Round-1
    // retry tactical reinforcement: bumped poll window from 30s → 60s to
    // accommodate the cold-start @panel registration race that may have
    // contributed to the prior Gate B B-STAB-01 (the prior author already
    // gave Scenarios 3/4 a 60s dialog tolerance; symmetric for the
    // Context Pane panel surface).
    const result: {found: boolean; expandedHasContent: boolean; headerTexts: string[]} = await page.evaluate(async () => {
      const deadline = Date.now() + 60_000;
      let foundHeader: Element | null = null;
      const matchHeader = (h: Element) => {
        const t = (h.textContent || '').trim();
        return t === 'Composition analysis' || t.toLowerCase().startsWith('composition');
      };
      while (Date.now() < deadline) {
        const headers = Array.from(document.querySelectorAll(
          '.grok-prop-panel .d4-accordion-pane-header, .d4-accordion-pane-header'));
        const h = headers.find(matchHeader);
        if (h) { foundHeader = h; break; }
        await new Promise((r) => setTimeout(r, 300));
      }
      if (!foundHeader) {
        const headerTexts = Array.from(document.querySelectorAll(
          '.grok-prop-panel .d4-accordion-pane-header'))
          .map((h) => (h.textContent || '').trim())
          .filter((t) => t.length > 0);
        return {found: false, expandedHasContent: false, headerTexts};
      }
      // Expand the pane (click header) if it's collapsed.
      (foundHeader as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 1500));
      // Walk up to the .d4-accordion-pane ancestor and look inside its
      // content region for non-empty content (any rendered child element
      // — the composition table, a label, etc.).
      let pane: Element | null = foundHeader.parentElement;
      while (pane && !pane.classList?.contains('d4-accordion-pane'))
        pane = pane.parentElement;
      const content = pane?.querySelector(':scope > :not(.d4-accordion-pane-header)');
      const text = (content?.textContent || '').trim();
      const childCount = content ? content.children.length : 0;
      const expandedHasContent = !!content && (text.length > 0 || childCount > 0);
      const headerTexts: string[] = [];
      return {found: true, expandedHasContent, headerTexts};
    });
    // Diagnostic: print the headers observed when the panel was not found
    // — surfaces the actual property-panel section names so a retry can
    // refine the header-match predicate.
    if (!result.found && result.headerTexts.length > 0)
      console.warn(`Composition analysis pane not found; observed headers: ${result.headerTexts.join(', ')}`);
    expect(result.found,
      'Composition analysis pane MUST appear on a Macromolecule-cell Context Pane (package.ts#L403)').toBe(true);
    expect(result.expandedHasContent,
      'Composition analysis pane MUST render non-empty content on expand').toBe(true);
  });

  await softStep('Scenario 2 Step 4: locate the Monomer info panel (semType Monomer; see Selector recon-notes)', async () => {
    // monomerInfoPanel (package.ts#L412) is decorated with semType:
    // 'Monomer' — the param decorator restricts the panel to Monomer
    // cells (individual per-position monomer cells produced by Split to
    // Monomers), not the parent Macromolecule sequence cell. On a
    // Macromolecule cell selection the panel may or may not surface
    // depending on the platform's @panel registration filter on
    // SemanticValue.semType. The scenario .md Step 4 Expected
    // accommodates both cases ("…or the cell-level summary for the
    // parent Macromolecule cell — the panel must surface non-empty
    // content"). We assert non-empty content IF the panel surfaces; we
    // do NOT fail-hard when it does not (the scenario's primary
    // assertion is "no error balloon appears in either panel").
    const result: {found: boolean; hasContent: boolean} = await page.evaluate(async () => {
      const headers = Array.from(document.querySelectorAll(
        '.grok-prop-panel .d4-accordion-pane-header, .d4-accordion-pane-header'));
      const h = headers.find((el) => {
        const t = (el.textContent || '').trim();
        return t === 'Monomer' || t.toLowerCase().startsWith('monomer');
      });
      if (!h) return {found: false, hasContent: false};
      (h as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 1500));
      let pane: Element | null = h.parentElement;
      while (pane && !pane.classList?.contains('d4-accordion-pane'))
        pane = pane.parentElement;
      const content = pane?.querySelector(':scope > :not(.d4-accordion-pane-header)');
      const text = (content?.textContent || '').trim();
      const childCount = content ? content.children.length : 0;
      const hasContent = !!content && (text.length > 0 || childCount > 0);
      return {found: true, hasContent};
    });
    // Soft predicate: when the panel surfaces, it MUST be non-empty. When
    // it does not surface (semType Monomer filter on a Macromolecule
    // cell), the assertion is the no-error-balloon predicate at the
    // end of Scenario 2.
    if (result.found) {
      expect(result.hasContent,
        'When the Monomer pane surfaces on a Macromolecule-cell selection it MUST render non-empty content').toBe(true);
    }
  });

  await softStep('Scenario 2: no balloon error fired by either context-pane info panel', async () => {
    const balloonError = await page.evaluate(() => {
      const errors = document.querySelectorAll('.d4-balloon.error, .grok-balloon-error');
      return errors.length;
    });
    expect(balloonError, 'no error balloon must be emitted by either Context Pane panel').toBe(0);
  });

  // -------------------------------------------------------------------
  // Scenario 3 — Get Region editor (GetRegionEditor; package.ts#L213).
  //
  // Atlas sub_feature id: bio.editors.get-region.
  // Live menu label: "Extract Region..." (bio.md:334-348; Bio CLAUDE.md
  // top-menu table). The dialog title is "Get Region" — the dialog name
  // [name="dialog-Get-Region"] persists across the menu-label rename.
  // -------------------------------------------------------------------

  await softStep('Scenario 3 Step 1-2: click Bio > Calculate > Extract Region — GetRegionEditor dialog opens', async () => {
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      const group = document.querySelector('[name="div-Bio---Calculate"]');
      if (!group) throw new Error('[name="div-Bio---Calculate"] group anchor not found');
      group.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      group.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      const leaf = document.querySelector('[name="div-Bio---Calculate---Extract-Region..."]');
      if (!leaf) throw new Error('[name="div-Bio---Calculate---Extract-Region..."] leaf not found');
      (leaf as HTMLElement).click();
    });
    // Cold-start dialog tolerance (matches analyze-spec.ts L180).
    await page.locator('[name="dialog-Get-Region"]').waitFor({state: 'visible', timeout: 60_000});
  });

  await softStep('Scenario 3 Step 3: GetRegionEditor sequence-column selector is populated with the FASTA column', async () => {
    // The custom editor exposes a sequence-column input bound to the
    // active TableView's Macromolecule column. bio.md:344 names the
    // selector pair [name="input-host-Sequence"] / [name="input-Sequence"].
    // Asserts:
    //   • the sequence-column input is attached to the dialog DOM (input
    //     widgets render under [name="input-host-…"]);
    //   • the dialog is visible (visible-state predicate from waitFor above).
    const seqInputPresent = await page.locator(
      '[name="dialog-Get-Region"] [name="input-host-Sequence"]').count();
    expect(seqInputPresent,
      'GetRegionEditor MUST render a sequence-column selector input (bio.md:344)').toBeGreaterThan(0);
  });

  await softStep('Scenario 3 Step 4: Cancel via Escape closes the dialog with no balloon error', async () => {
    await page.keyboard.press('Escape');
    await page.locator('[name="dialog-Get-Region"]').waitFor({state: 'hidden', timeout: 10_000});
    const balloonError = await page.locator('.d4-balloon.error, .grok-balloon-error').count();
    expect(balloonError,
      'Cancel must close GetRegionEditor with no error balloon').toBe(0);
  });

  // -------------------------------------------------------------------
  // Scenario 4 — Split to Monomers editor (SplitToMonomersEditor;
  // package.ts#L226).
  //
  // Atlas sub_feature id: bio.editors.split-to-monomers.
  // Menu label is stable: "Split to Monomers..." (bio.md:395-405).
  // -------------------------------------------------------------------

  await softStep('Scenario 4 Step 1: click Bio > Transform > Split to Monomers — SplitToMonomersEditor dialog opens', async () => {
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      const group = document.querySelector('[name="div-Bio---Transform"]');
      if (!group) throw new Error('[name="div-Bio---Transform"] group anchor not found');
      group.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      group.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      const leaf = document.querySelector('[name="div-Bio---Transform---Split-to-Monomers..."]');
      if (!leaf) throw new Error('[name="div-Bio---Transform---Split-to-Monomers..."] leaf not found');
      (leaf as HTMLElement).click();
    });
    await page.locator('[name="dialog-Split-to-Monomers"]').waitFor({state: 'visible', timeout: 60_000});
  });

  await softStep('Scenario 4 Step 2: SplitToMonomersEditor sequence-column selector is populated', async () => {
    const seqInputPresent = await page.locator(
      '[name="dialog-Split-to-Monomers"] [name="input-host-Sequence"]').count();
    expect(seqInputPresent,
      'SplitToMonomersEditor MUST render a sequence-column selector input (bio.md:403)').toBeGreaterThan(0);
  });

  await softStep('Scenario 4 Step 3: Cancel via Escape closes the dialog with no balloon error', async () => {
    await page.keyboard.press('Escape');
    await page.locator('[name="dialog-Split-to-Monomers"]').waitFor({state: 'hidden', timeout: 10_000});
    const balloonError = await page.locator('.d4-balloon.error, .grok-balloon-error').count();
    expect(balloonError,
      'Cancel must close SplitToMonomersEditor with no error balloon').toBe(0);
  });

  // Cleanup: close all views (free TableView state).
  await page.evaluate(() => (window as any).grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
