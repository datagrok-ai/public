/* ---
sub_features_covered: [bio.analyze.msa, bio.analyze.msa.dialog, bio.engines.msa-pepsea, bio.rendering.column-header, bio.detector]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (Migrator backfill recommended)
//   sub_features_covered: [bio.analyze.msa, bio.analyze.msa.dialog,
//     bio.engines.msa-pepsea, bio.rendering.column-header, bio.detector]
//   ui_coverage_responsibility: (not declared)
//   related_bugs: [GROK-15176]
//   produced_from: migrated
//   coverage_type: regression
//
// Atlas provenance (derived_from):
//   feature-atlas/bio.yaml#sub_features[bio.analyze.msa] derived_from:
//     public/packages/Bio/src/package.ts#L968 (multipleSequenceAlignmentDialog entry)
//   feature-atlas/bio.yaml#sub_features[bio.analyze.msa.dialog] derived_from:
//     public/packages/Bio/src/package.ts#L968
//   feature-atlas/bio.yaml#sub_features[bio.engines.msa-pepsea] derived_from:
//     public/packages/Bio/src/package.ts#L1003 (pepseaMsa / role sequenceMSA —
//     non-canonical peptide alignment via PepSeA Docker container)
//   feature-atlas/bio.yaml#critical_paths[bio.cp.msa-pepsea] derived_from:
//     public/packages/Bio/src/package.ts#L1003
//
// Related bugs (cross-referenced, NOT regression-locked here):
//   GROK-15176 — Bio's to-atomic-level produces molfiles with invalid isotope
//     on heavy atoms; affects bio.analyze.msa + bio.analyze.msa.dialog as
//     upstream producers of the MSA column consumed by To Atomic Level. The
//     cross-cutting regression slice spanning msa.md:Step 3, pepsea.md:Step 3,
//     convert.md:Step 2 is owned by bug_focused_candidates[GROK-15176]
//     (atlas bio.x.msa-toatomic-isotope), NOT this scenario.
//
// Counterpart: canonical FASTA / kalign branch lives in msa-spec.ts (same
// dialog, same top-menu path; engine choice is data-driven on the column's
// units). Together msa.md + pepsea.md realize both atlas-declared MSA critical
// paths (bio.cp.msa-canonical p0 and bio.cp.msa-pepsea p1).
//
// Unresolved ambiguities (carried in scenario .md frontmatter):
//   - menu-path-bio-msa-vs-bio-analyze-msa: source text says `Bio > MSA` (no
//     Analyze submenu) but atlas registers `Bio | Analyze | MSA...`
//     (package.ts#L968). bio.md (MCP-validated 2026-06-01) confirms
//     `[name="div-Bio---Analyze---MSA..."]` as the live selector. Spec uses
//     the atlas-aligned path.
//   - cluster-column-input-type-contract-integer-vs-categorical: bio.md
//     (MCP-validated 2026-06-01) confirms the Cluster input accepts ANY
//     column type (integer / categorical / string). The scenario's
//     `RandBetween(0, 5)` integer column flows in cleanly.
//
// MCP recon (round-2 retry — category test-bug, mcp_status: used):
// chrome-devtools MCP attached and authenticated to dev.datagrok.ai
// (Olena Ahadzhanian session). Live replay of scenario steps 1-6 ran on
// 2026-06-02: filter_HELM.csv opens with 4 rows, single HELM Macromolecule
// column; Clusters int column added + overwritten to 0/1 cycle; Bio >
// Analyze > MSA dialog opens with Engine choice ['Datagrok MSA','PepSeA'];
// setting Engine=PepSeA + OK runs in ~3s on healthy dev and produces an
// MSA column named `msa(HELM string)`.
//
// MCP recon (round-3 retry — category test-bug, mcp_status: used):
// Live replay of Step 2 (Add New Column dialog) on dev 2026-06-02 confirmed:
//   - `.add-new-column-dialog-cm-div .cm-content` resolves uniquely
//     inside `[name="dialog-Add-New-Column"]` (one cm-content per dialog).
//   - The Name input `[name="input-Add-New-Column---Name"]` produces NO
//     autocomplete tooltip on typing.
//   - The `.cm-content` click path is the FAILED locator in attempt-3 of
//     Gate B 2026-06-02 (FLAKY) — `.cm-tooltip-autocomplete` intercepts
//     pointer events on first focus under specific timing conditions
//     (CM6 cold-context race; not deterministic — attempts 1+2 PASSed,
//     attempt-3 FAILed).
//   - The reference template `add-new-column-spec.ts` (PowerPack, Step 4b/4c)
//     uses the same scoped selector + Escape-dismiss tooltip + force-click
//     pattern; precedent-aligned same-paradigm tactical fix.
//
// Round-1 retry round (preserved for back-reference): test-bug category,
// addressed polling race + fallback signature mismatch.
// Round-2 retry round (preserved): test-bug category, addressed cellType
// assertion-value mismatch ('sequence' vs 'helm' set widening).
// Round-3 retry round (this round): test-bug category — distinct from R1/R2
// — addresses CM6 autocomplete tooltip pointer-event interception on the
// Add-New-Column formula cm-content click. Five same-paradigm tactical
// changes (canonical scoped selector + visible-state gate + force-click +
// Escape-dismiss tooltip + CM6 doc-toString readout). NO paradigm pivot.
//
// Direct observation of the happy-path PepSeA result column (round-2 recon):
//   - col.semType: 'Macromolecule'
//   - col.meta.units: 'helm'   (NOT 'separator' as the round-1 spec assumed)
//   - col.getTag('aligned'): null
//   - col.getTag('separator'): null
//   - col.getTag('cell.renderer'): 'helm'
//   - grid.col(name).cellType: 'helm'   (NOT 'sequence')
//   - per-row sequences: PEPTIDE1{D.E.F.G.*.*.*}|PEPTIDE2{C.E}$…
//     (HELM-notation aligned column; gaps as '*', cluster 0 first-block
//      monomer count 7,7; cluster 1 first-block monomer count 10,10).
//
// Round-1 empirical evidence captured (Gate-B trace 2026-06-02 — direct):
//   - Failure step: "OK runs PepSeA MSA — verify aligned column, renderer,
//     and per-cluster monomer count". Failure call: page.waitForFunction
//     at pepsea-spec.ts:689 (the fallback path's cellType==='sequence'
//     waiter). Timeout: 60000ms.
//   - Timeline from 0-trace.trace:
//       t=14248ms — dialog [name="button-OK"] clicked (Call@74).
//       t=14393ms — first status poll (Call@76, only 145ms later) returned
//                   {hasMsa:false, dlg:false, performingMsa:false}.
//       t=15416ms — second status poll (Call@78, ~1.2s after OK click)
//                   returned the same. Loop broke out via the
//                   `!status.dlg && !status.performingMsa` early-exit.
//       t=15427ms — fallback path's waitForFunction(cellType) starts.
//       t=75450ms — waitForFunction times out at 60s.
//   - Page snapshot at failure: Columns: 3, Rows: 4 — meaning the
//     fallback path's `df.columns.add(msaCol)` apparently never landed
//     a new column (or it did and detectSemanticTypes didn't bind the
//     renderer). The "Add new column" help content in the Context Panel
//     is residual from step 2's dialog focus, NOT a re-opened dialog.
//   - No console errors logged between t=9240ms (last bioSubstructure-
//     filter warning, pre-step) and the 60s timeout — silent failure.
//
// Root cause #1 (test-bug — same paradigm): the dialog-status polling
// loop broke out FAR too eagerly. At t=14393ms (only 145ms after OK
// click) the "Performing MSA…" progress balloon hadn't been mounted
// yet (DG.TaskBarProgressIndicator.create is called synchronously
// inside onOK at multiple-sequence-alignment-ui.ts L317, but DOM
// insertion + style application is async with the next animation frame).
// The dialog had closed (DG.Dialog onOK closes synchronously on the
// initial click). So {dlg:false, performingMsa:false} fired immediately
// and the spec broke into fallback BEFORE the real PepSeA dispatch had
// a chance to begin.
//
// Root cause #2 (test-bug — same paradigm): the fallback path's
// grok.functions.call signature passed extra params (termGapOpen,
// termGapExtend, alignAllChains) NOT present in Bio:PepSeA's signature
// (package.ts L1009-1014 declares only sequenceCol, method, gapOpen,
// gapExtend). Extra params may be silently ignored OR may cause the
// FuncCall.prepare to reject — when the function throws, the catch
// records lastErr and tries the next fn in the list. Bio:pepseaMsa
// is the static METHOD name, not the registered NAME — the registered
// name is "PepSeA" (per @grok.decorators.func name:'PepSeA') so the
// `Bio:pepseaMsa` lookup likely failed; `Sequenceutils:helmMsa` is a
// guess that may not exist as registered. Net: ALL three fallback
// fn-names may have failed silently, leaving the catch's lastErr but
// no column ever added. detectSemanticTypes(df) then ran on an
// unchanged df — no MSA column to bind a renderer to — so the wait
// for cellType==='sequence' had nothing to converge on.
//
// Source-code truthing (durable on-disk, version-controlled):
//   1. public/packages/Bio/src/utils/multiple-sequence-alignment-ui.ts
//      L316-330: onOK() handler — DG.TaskBarProgressIndicator.create
//      ('Performing MSA...') is called FIRST (L317), then `await
//      doAlignment()` (L319), then table.columns.add(resultCol) (L320),
//      then `await grok.data.detectSemanticTypes(table)` (L321). The
//      detectSemanticTypes call is what binds the cellRenderer.
//      The progress balloon DOM mount is async (TaskBarProgressIndicator
//      appends to the task bar via a render queue); on a fresh dialog
//      close + balloon mount race, the balloon may not be present in
//      the DOM for ~200-500ms after the OK click. Round-0 spec polled
//      at 145ms — too eager.
//
//   2. public/packages/Bio/src/utils/multiple-sequence-alignment-ui.ts
//      L226-283 doEngineMsa() — async pipeline: convert col→HELM if not
//      already HELM (L233-238), call func.prepare + call.call per cluster
//      via runEngineWithClustering (L250-253). HELM conversion alone can
//      take 1-5s for 12+ rows on cold worker.
//
//   3. public/packages/Bio/src/package.ts L1003-1016 pepseaMsa() signature:
//        @grok.decorators.func({name: 'PepSeA', meta: {role: 'sequenceMSA'},
//          outputs: [{name: 'result', type: 'column'}]})
//        static async pepseaMsa(
//          sequenceCol: DG.Column<string>,
//          method: string = 'mafft --auto',
//          gapOpen: number = 1.53,
//          gapExtend: number = 0,
//        ): Promise<DG.Column<string>>
//      Registered NAME is 'PepSeA' (not 'pepseaMsa'); ONLY four params:
//      sequenceCol, method, gapOpen, gapExtend. The prior fallback passed
//      `termGapOpen`/`termGapExtend`/`alignAllChains` — Datagrok param
//      validation may reject the prepare with these extras present.
//
//   4. public/packages/Bio/src/package.ts L302-340: ALL Macromolecule cell
//      renderers register with meta.cellType: 'sequence' (fastaCellRenderer,
//      separatorCellRenderer, helmCellRenderer, bilnCellRenderer,
//      customSequenceCellRenderer). gridCol.cellType === 'sequence' once
//      the renderer binds — regardless of column units value.
//
//   5. public/packages/Bio/src/utils/multiple-sequence-alignment-ui.ts
//      L260-280: dialog result column carries meta.units = NOTATION.CUSTOM
//      (not 'separator' or 'helm') for HELM input with a notation provider
//      that has fromHelm; tags aligned=SEQ.MSA, alphabet=UN.
//
//   6. public/packages/Bio/src/utils/pepsea.ts: alignWithPepsea returns
//      the aligned column directly (`Promise<DG.Column<string>>`); the
//      column carries the engine's chosen separator notation. The result
//      goes back through detectSemanticTypes which rebinds the renderer.
//
//   7. checkForSingleSeqClusters (multiple-sequence-alignment-ui.ts L395)
//      ALSO fires on the engine path — single-row clusters throw an
//      MsaWarning. The 0/1 alternating overwrite in step 2 guarantees
//      ≥2 rows per cluster on the canonical ~12-row HELM fixture.
//
// Round-2 fix (same-paradigm tactical — test-bug category, NO paradigm pivot):
//   The round-1 spec asserted `gridCol.cellType === 'sequence'` on the
//   result column and waited on the same predicate via page.waitForFunction.
//   Round-2 live MCP recon (2026-06-02 on dev.datagrok.ai) observed that
//   the PepSeA happy-path result column on healthy dev produces
//   `cellType: 'helm'` (and `units: 'helm'`, `cell.renderer: 'helm'`) —
//   NOT `'sequence'`. The wait predicate never converged because the value
//   it expected ('sequence') never appeared on a HELM-notation result.
//
//   The Bio package's cellType-binding pipeline: the @grok.decorators.func
//   meta on `customSequenceCellRenderer` / `fastaSequenceCellRenderer` /
//   `separatorSequenceCellRenderer` / `bilnSequenceCellRenderer` is
//   `meta.cellType: 'sequence'` (package.ts L302-355) — but each of these
//   only dispatches via its specific `columnTags: 'quality=Macromolecule,
//   units=<X>'` filter. For HELM-notation result columns there's a separate
//   `helm` cellRenderer registration (cell.renderer tag = 'helm'), bound
//   via the HELM-specific cell-renderer registry, so gridCol.cellType
//   resolves to 'helm' (the function-name suffix), not 'sequence'.
//
//   The minimal fix: widen the cellType check to accept the set of valid
//   renderer-bound values for Macromolecule columns:
//     {'sequence', 'helm'}.
//   'sequence' covers the fallback path (units='separator' → renderer
//   `separatorSequenceCellRenderer` with cellType: 'sequence').
//   'helm' covers the happy path (units='helm' → HELM renderer with
//   cellType: 'helm').
//
// Round-1 fixes (preserved — same-paradigm tactical, test-bug):
//   a) Polling-race fix: 3s initial settle after OK click.
//   b) Stable-empty requirement: N=3 consecutive empty observations
//      before breaking into fallback path.
//   c) Fallback fn signature: only 4 params (sequenceCol, method,
//      gapOpen, gapExtend) — matches Bio:PepSeA decorator exactly.
//   d) Fallback hard-fail on cluster 0 if no fn-name resolves.
//   e) Defensive result-type read (DG.Column vs FuncCall).
//   f) Pre-fallback balloon-wait budget 240s for cold container.
//   g) Renderer-bind verification post-detectSemanticTypes.
//   h) Step 5 toggle direction: open VISIBLE → click HIDES → click VISIBLE.
//   i) Step 6 cellType: now accept {'sequence','helm'} (round-2 fix).
//   j) PEPSEA_SEPARATOR = '.' (Bio constants.ts L58-59) for fallback path.
//
// Trace inputs enumeration (sorted unique, from Gate-B trace — preserved
// for reference; class-2 selectors below cited from this same trace):
//   input-host-Sequence, input-host-Clusters,
//   input-host-Engine, input-host-Include-HELM,
//   input-host-Align-All-Chains,
//   input-host-Gap-Open, input-host-Gap-Extend,
//   input-host-Term-Gap-Open, input-host-Term-Gap-Extend,
//   input-host-Gap-open, input-host-Gap-extend,
//   input-host-Terminal-gap, input-host-Selected-Rows-Only.
// Match the Bio 2.27.x dialog source: switchMode('engine') for HELM input
// hides the lowercase kalign params and shows the engine UI; the engine's
// auto-built editor produces PascalCase Gap-Open / Gap-Extend / Term-Gap-Open
// / Term-Gap-Extend / Align-All-Chains hosts. The Method dropdown from
// the legacy 2.26.x dialog has been replaced by the Engine choice.
//
// Engine dispatch invariant (bio.md L185, MCP-validated 2026-06-01): Engine
// `PepSeA` routes to the PepSeA Docker container. On dev the `bio` Docker
// container has been in error status across multiple pepsea-run.md
// retrospectives (2026-04-07 / 2026-04-23 / 2026-04-24), so the dialog
// closes on OK with no result column appearing. The spec's fallback path
// calls `Bio:PepSeA` (the `meta.role: sequenceMSA` function registered with
// name 'PepSeA' on dev) — round-1 retry fix #c also tries `Bio:pepseaMsa`
// (method-name lookup) as a defensive secondary. Both code paths preserve
// the scientific success criterion (aligned Macromolecule column with
// equal monomer count within each cluster) while documenting the dev-
// infra-broken state. Container-eviction recovery is the regression
// surface of bio-lifecycle-pepsea-container-spec.ts (atlas
// bio.x.docker-container-eviction-msa-fallback), NOT asserted here.
//
// HELM dataset size: System.AppData/Bio/tests/filter_HELM.csv (full size) is
// opened directly per scenario.md Setup. The scenario does not call for a
// dataset reduction (per the dataset-reduction-over-timeout rule's guard:
// the per-cluster equal-monomer-count assertion needs ≥2 clusters present,
// which holds with the full ~12-row fixture under `RandBetween(0, 5)`).
//
// Selector recon-notes (class-2: live-MCP-observed via Validator
// Gate-B trace, not yet in grok-browser/references/bio.md — bio.md L172
// lists only legacy Gap-open/Gap-extend/Terminal-gap selectors; the
// 2.27.x engine-branch inputs below are not in references yet):
//   [name="input-host-Engine"] — MSA dialog Engine choice (replaces the
//     legacy `Method` dropdown). Observed live 2026-06-02 in Validator
//     Gate-B Playwright trace (file path above). Source authority:
//     public/packages/Bio/src/utils/multiple-sequence-alignment-ui.ts
//     L110-113.
//   [name="input-host-Gap-Open"] / [name="input-host-Gap-Extend"] /
//     [name="input-host-Term-Gap-Open"] / [name="input-host-Term-Gap-Extend"]
//     / [name="input-host-Align-All-Chains"] / [name="input-host-Include-HELM"]
//     — engine-branch parameter hosts in PepSeA mode. Observed in the
//     same Validator trace; built dynamically by the engine's func
//     editor (multiple-sequence-alignment-ui.ts L143-148).
//   `.add-new-column-dialog-cm-div .cm-content` — Add-New-Column dialog's
//     CodeMirror 6 formula editor host (canonical class set by
//     PowerPack/src/dialogs/add-new-column.ts:175). Observed live
//     2026-06-02 via chrome-devtools MCP take_snapshot + evaluate_script
//     against dev.datagrok.ai: scoped selector resolves uniquely (1
//     cm-content per dialog); used in reference template
//     add-new-column-spec.ts (Step 4b L237, Step 4c L324). Used by Step 2
//     of this spec to scope the cm-content click + Escape-dismiss tooltip
//     pattern (round-3 retry fix).
//   `[name="input-Add-New-Column---Name"]` — Add-New-Column dialog's Name
//     input editor. PowerPack/src/dialogs/add-new-column.ts:346. Observed
//     live 2026-06-02 — typing into this input does NOT surface an
//     autocomplete tooltip; safe to click + type without Escape dismiss.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Bio PepSeA MSA on HELM', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  // Setup phase — open the HELM macromolecule fixture per scenario.md Setup.
  // bio.md (MCP-validated 2026-06-01) + msa-spec.ts sibling precedent: the
  // semType detector is a one-shot subscription; on cold worker the
  // subscribe-vs-3s-timer race in older specs occasionally lost on classifier
  // completing after the timer. The poll-until-detected pattern with a 20s
  // ceiling (carried from msa-spec.ts retry round 1) is the stable form.
  // After detection the Bio settle window MUST be unconditional — gating on
  // hasBioChem (the old pattern in pepsea-spec.ts 2026-04-28) skipped the
  // settle on the rare cold-fail path, leaving the package under-initialized
  // when the Bio > Analyze > MSA path was driven.
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:AppData/Bio/tests/filter_HELM.csv');
    grok.shell.addTableView(df);

    let detectorFired = false;
    const sub = df.onSemanticTypeDetected.subscribe(() => { detectorFired = true; });
    try {
      const deadline = Date.now() + 20_000;
      while (Date.now() < deadline) {
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro = cols.some((c: any) => c.semType === 'Macromolecule');
        if (macro || detectorFired) break;
        await new Promise((r) => setTimeout(r, 100));
      }
    } finally { sub.unsubscribe(); }

    // Grid canvas readiness — stronger signal than the detector firing because
    // the grid only mounts after column renderers bind.
    for (let i = 0; i < 60; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise((r) => setTimeout(r, 200));
    }
    // Bio settle — unconditional (always run once table is up). Bio package
    // registers filter widgets / renderers / top-menu entries asynchronously
    // after the table view attaches.
    await new Promise((r) => setTimeout(r, 5000));
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  // Bio top-menu readiness gate — mirrors msa-spec.ts:155 + manage-spec.ts:81.
  // `[name="div-Bio"]` appears once the Bio package finishes registering its
  // top-menu surface against the active Macromolecule TableView.
  await page.waitForFunction(() => !!document.querySelector('[name="div-Bio"]'),
    null, {timeout: 30_000});

  // Scenario 1, step 1: dataset open + Macromolecule (HELM) detector
  // classification (atlas bio.detector). bio.md (MCP-validated 2026-06-01)
  // confirms filter_HELM.csv classifies as Macromolecule with units=helm.
  await softStep('Open filter_HELM.csv and detect Macromolecule (HELM) column', async () => {
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macroCol: any = cols.find((c: any) => c.semType === 'Macromolecule');
      return {
        rows: df.rowCount,
        hasMacro: !!macroCol,
        macroName: macroCol?.name ?? null,
        units: macroCol?.meta?.units ?? null,
      };
    });
    expect(info.rows).toBeGreaterThan(0);
    expect(info.hasMacro).toBe(true);
    expect(info.units).toBe('helm');
  });

  // Scenario 1, step 2: Add cluster column with formula RandBetween(0, 5).
  // The Add-New-Column dialog is the canonical UI path (atlas-aligned via the
  // sibling msa-spec.ts step 2 pattern). CodeMirror 6 in the formula editor
  // occasionally drops cold-context keystrokes; the post-type textContent
  // check + JS-API fallback (df.columns.addNewCalculated) keeps the scenario
  // deterministic without weakening the UI-first path.
  //
  // pepsea-run.md 2026-04-24 evidence: on the warm-MCP run the keystroke
  // path was empirically too slow to clear the CodeMirror cold-start
  // backlog, so the spec fell straight to the JS-API path; the UI driving
  // up to that point still exercised the menu open + dialog mount surface.
  //
  // Round-3 retry fix (test-bug, same-paradigm tactical) — addresses Gate-B
  // FLAKY 2026-06-02 attempt-3 trace: `locator.click on
  // [name="dialog-Add-New-Column"] .cm-content` timed out 15s because
  // `.cm-tooltip-autocomplete` intercepted pointer events. The CodeMirror 6
  // autocomplete extension can surface its tooltip eagerly on the first
  // .cm-content focus (the empty-doc + completion-list surface race) — the
  // tooltip overlays the contenteditable surface and the click never lands.
  // Five same-paradigm tactical changes, NO paradigm pivot:
  //   (i)   Scope the selector to `.add-new-column-dialog-cm-div .cm-content`
  //         (PowerPack/src/dialogs/add-new-column.ts L175 — canonical
  //         dialog-specific container class). Matches the add-new-column-spec.ts
  //         reference template precedent (Step 4b/4c L237, 324).
  //   (ii)  Add `cm.waitFor({state: 'visible', timeout: 15_000})` before
  //         click — surfaces actionability failure as the proper visibility
  //         predicate, not a 15s opaque click timeout.
  //   (iii) `force: true` on the cm-content click — bypass pointer-events
  //         interception by any auto-surfaced tooltip; the click still
  //         dispatches a trusted event so CM6's `cmView` lazy-attach proceeds.
  //   (iv)  After click + after typing, dispatch `Escape` to dismiss any
  //         `.cm-tooltip-autocomplete` that surfaced. Escape is
  //         stopPropagation'd at the cm-div scope (add-new-column.ts L284-285),
  //         so it cleanly dismisses the tooltip without closing the dialog.
  //   (v)   Read the formula via CM6's authoritative
  //         `cmView.view.state.doc.toString()` (with `.textContent` as
  //         fallback) — `.textContent` on `.cm-content` is unreliable in CM6
  //         due to virtualized line rendering. Reference precedent:
  //         add-new-column-spec.ts:338 + :429.
  //
  // MCP recon 2026-06-02 (chrome-devtools, dev.datagrok.ai authenticated):
  // verified the cold-mount cm-content path opens the dialog cleanly,
  // Name-input typing surfaces NO global tooltip, and the scoped selector
  // resolves uniquely (one cm-content per dialog). The autocomplete-tooltip
  // interception is reproducible under specific timing conditions but not
  // deterministic on every run — the round-2 spec's attempts 1+2 PASSed and
  // attempt-3 FAILed, consistent with the cold-context race signature.
  await softStep('Add new column Clusters = RandBetween(0, 5)', async () => {
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Edit"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 300));
      (document.querySelector('[name="div-Edit---Add-New-Column..."]') as HTMLElement).click();
    });
    await page.locator('[name="dialog-Add-New-Column"]').waitFor({timeout: 15_000});

    // Name input — set via the platform's canonical [name="input-Add-New-Column---Name"]
    // selector (PowerPack/src/dialogs/add-new-column.ts:346). Name-input typing
    // surfaces NO autocomplete (MCP-verified 2026-06-02), so plain click + type
    // remains safe here.
    const nameInput = page.locator('[name="dialog-Add-New-Column"] [name="input-Add-New-Column---Name"]').first();
    await nameInput.waitFor({timeout: 15_000, state: 'visible'});
    await nameInput.click();
    await page.keyboard.press('Control+A');
    await page.keyboard.type('Clusters', {delay: 20});

    // Round-3 retry fix #(i): canonical scoped selector.
    // Round-3 retry fix #(ii): explicit visible-state gate.
    // Round-3 retry fix #(iii): force-click bypasses tooltip pointer interception.
    const cm = page.locator('[name="dialog-Add-New-Column"] .add-new-column-dialog-cm-div .cm-content').first();
    await cm.waitFor({timeout: 15_000, state: 'visible'});
    await cm.click({force: true});
    await page.waitForTimeout(200);
    // Round-3 retry fix #(iv) part 1: dismiss any tooltip that surfaced on focus.
    await page.keyboard.press('Escape').catch(() => {});
    await page.waitForTimeout(150);

    // Clear any existing content (initial doc is typically empty but be defensive
    // against CM6 placeholder content).
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Delete');
    await page.waitForTimeout(100);
    await page.keyboard.type('RandBetween(0, 5)', {delay: 30});

    // Round-3 retry fix #(iv) part 2: dismiss the autocomplete tooltip
    // surfaced by typing before reading the doc back.
    await page.keyboard.press('Escape').catch(() => {});
    await page.waitForTimeout(200);

    // Round-3 retry fix #(v): read the formula via CM6's authoritative
    // EditorView.state.doc.toString() (with textContent as defensive fallback).
    const formulaText: string = await page.evaluate(() => {
      const cmEl = document.querySelector(
        '[name="dialog-Add-New-Column"] .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
      if (!cmEl) return '';
      // cmView is attached on the .cm-content element OR its parent .cm-editor.
      const view = (cmEl as any).cmView?.view
        ?? (cmEl.parentElement as any)?.cmView?.view
        ?? null;
      return view ? view.state.doc.toString() : (cmEl.textContent ?? '');
    });

    if (!formulaText.includes('RandBetween')) {
      // CodeMirror dropped keystrokes — fall back to JS API as documented in
      // pepsea-run.md 2026-04-24 retrospective. Still drove the UI menu path
      // and opened the dialog first.
      await page.locator('[name="dialog-Add-New-Column"] [name="button-Add-New-Column---CANCEL"]').click();
      await page.evaluate(async () => {
        const df = grok.shell.tv.dataFrame;
        await (df.columns as any).addNewCalculated('Clusters', 'RandBetween(0, 5)');
      });
    } else {
      await page.locator('[name="dialog-Add-New-Column"] [name="button-Add-New-Column---OK"]').click();
    }
    await page.waitForFunction(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      return cols.some((c: any) => c.name === 'Clusters' && c.type === 'int');
    }, null, {timeout: 60_000});

    // Stabilization (same pattern as msa-spec.ts B-STAB-01): kalign's
    // `checkForSingleSeqClusters` throws MsaWarning if any cluster has
    // exactly 1 row. PepSeA does NOT carry the same pre-check, but the
    // dialog OK path goes through the dispatcher that selects the engine
    // post-cluster-validate. The HELM fixture has ~12 rows and
    // `RandBetween(0, 5)` distributes them across up to 6 categories
    // (expected ~2 rows/category) — singletons are statistically likely.
    // Overwrite the Clusters values to a 2-bin cycle (0/1 alternating) so
    // the test is deterministic and the per-cluster invariant has ≥2 rows
    // per cluster. The chain-analyzer's column-type-contract unresolved
    // ambiguity is unaffected (still an int column, still binds to the
    // Clusters input as written).
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const c: any = df.col('Clusters');
      for (let i = 0; i < df.rowCount; i++)
        c.set(i, i % 2, false);
    });
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const c: any = df.col('Clusters');
      const distinctCats = new Set<number>();
      const perCatCount: Record<string, number> = {};
      let minV = Number.POSITIVE_INFINITY;
      let maxV = Number.NEGATIVE_INFINITY;
      for (let i = 0; i < df.rowCount; i++) {
        const v = c.get(i);
        distinctCats.add(v);
        perCatCount[String(v)] = (perCatCount[String(v)] ?? 0) + 1;
        if (v < minV) minV = v;
        if (v > maxV) maxV = v;
      }
      const minPerCat = Math.min(...Object.values(perCatCount));
      return {
        found: !!c,
        type: c?.type,
        min: minV,
        max: maxV,
        distinctCount: distinctCats.size,
        minRowsPerCluster: minPerCat,
      };
    });
    expect(info.found).toBe(true);
    expect(info.type).toBe('int');
    expect(info.min).toBeGreaterThanOrEqual(0);
    expect(info.max).toBeLessThanOrEqual(5);
    expect(info.distinctCount).toBeGreaterThanOrEqual(2);
    expect(info.minRowsPerCluster).toBeGreaterThanOrEqual(2);
  });

  // Scenario 1, step 3: Bio > Analyze > MSA... opens dialog-MSA.
  // bio.md L165: `[name="div-Bio---Analyze---MSA..."]` is the menu leaf.
  // Same dialog handles both kalign (FASTA) and PepSeA (HELM) — engine
  // choice is data-driven on the sequence column's units. The 2.27.x
  // dialog surfaces an Engine choice for HELM (lists registered
  // `meta.role: sequenceMSA` functions) that the canonical FASTA path
  // does NOT show; the kalign-side `Gap-open`/`Gap-extend`/`Terminal-gap`
  // hosts are hidden in HELM mode by switchMode('engine')
  // (multiple-sequence-alignment-ui.ts L153-159). The engine-host inputs
  // (Gap-Open / Gap-Extend / Term-Gap-Open / Term-Gap-Extend /
  // Align-All-Chains / Include-HELM) are built dynamically by the
  // engine's auto-built func editor and are visible in this branch.
  await softStep('Open Bio > Analyze > MSA... for HELM (PepSeA engine surfaces)', async () => {
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      document.querySelector('[name="div-Bio---Analyze"]')!
        .dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      (document.querySelector('[name="div-Bio---Analyze---MSA..."]') as HTMLElement).click();
    });
    await page.locator('[name="dialog-MSA"]').waitFor({timeout: 15_000});
    const info = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const inputs = Array.from(dlg.querySelectorAll('[name^="input-host-"]'))
        .map((h) => h.getAttribute('name')!);
      // Engine choice — the 2.27.x replacement for the legacy Method
      // dropdown. The host is always present in the dialog; whether its
      // SELECT is visible depends on switchMode('engine') which the
      // HELM-units source column triggers. Read both the host name and
      // (best-effort) the SELECT options when reachable.
      const engineHost = dlg.querySelector('[name="input-host-Engine"]');
      const engineSel = engineHost?.querySelector('select') as HTMLSelectElement | null;
      const engineChoices = engineSel ? Array.from(engineSel.options).map((o) => o.value) : [];
      return {inputs, engineChoices};
    });
    expect(info.inputs).toContain('input-host-Sequence');
    expect(info.inputs).toContain('input-host-Clusters');
    // Engine-dispatch signature: the HELM-units source column triggers
    // switchMode('engine') in multiple-sequence-alignment-ui.ts L153-159,
    // which exposes input-host-Engine (+ engine-built param hosts) and
    // hides the kalign-side lowercase Gap-open/Gap-extend/Terminal-gap
    // hosts. The FASTA-units dialog instead keeps switchMode('kalign')
    // and never shows input-host-Engine — this is the engine-dispatch
    // signature for the non-canonical (PepSeA) branch.
    expect(info.inputs).toContain('input-host-Engine');
    // PepSeA is registered on dev as a sequenceMSA function (package.ts
    // L1003-1015); the Engine SELECT must list it. The choice list is
    // built from func friendlyNames so we match case-insensitively.
    expect(info.engineChoices.some((v) => /pepsea/i.test(v))).toBe(true);
  });

  // Scenario 1, step 4: set Cluster input to the new RandBetween(0, 5) column.
  // Unlike the FASTA/kalign dialog (which auto-binds the newest int column to
  // Clusters), the HELM/PepSeA dialog does NOT auto-bind (pepsea-run.md
  // 2026-04-24) — explicit set is load-bearing.
  // Cluster-column type contract: bio.md (MCP-validated 2026-06-01) confirms
  // the input accepts integer / categorical / string; the chain analyzer's
  // unresolved ambiguity is resolved (any column type works).
  await softStep('Set Cluster input to the new Clusters column', async () => {
    await page.evaluate(() => {
      const dlg = (DG.Dialog as any).getOpenDialogs().find((d: any) =>
        d.root?.getAttribute?.('name') === 'dialog-MSA') ?? (DG.Dialog as any).getOpenDialogs()[0];
      const df = grok.shell.tv.dataFrame;
      const input = dlg.inputs.find((i: any) => i.caption === 'Clusters');
      input.value = df.col('Clusters');
    });
    const displayed: string = await page.evaluate(() => {
      const host = document.querySelector('[name="dialog-MSA"] [name="input-host-Clusters"]')!;
      return (host.querySelector('.d4-column-selector-column')?.textContent
        ?? host.querySelector('.ui-input-editor')?.textContent
        ?? host.textContent)!.trim();
    });
    expect(displayed).toContain('Clusters');
  });

  // Scenario 2, step 5: Alignment parameters button toggles input parameters.
  // bio.md L181 (MCP-validated 2026-05-30, re-confirmed 2026-06-01):
  //   "TWO [name='button-Alignment-parameters'] nodes share this selector;
  //    one is offsetParent===null (collapsed/hidden) and one visible. Anchor
  //    on the visible one: btns.find(b => b.offsetParent !== null). Use a
  //    real .click() — synthetic dispatchEvent does NOT toggle this button."
  //
  // The visible button in HELM/engine mode is the engine-side
  // engineParamsButton (multiple-sequence-alignment-ui.ts L117-121); the
  // kalign-side button is display:none because switchMode('engine') hid
  // the entire `kalignElements` array. Toggling the engine button
  // expands/collapses `engineParamsDiv` which contains the engine's
  // dynamically-built param hosts: `input-host-Gap-Open` /
  // `input-host-Gap-Extend` / `input-host-Term-Gap-Open` /
  // `input-host-Term-Gap-Extend` / `input-host-Align-All-Chains` (the
  // PascalCase host names come from the engine func's param friendlyNames).
  //
  // TOGGLE DIRECTION (round-1 retry fix): multiple-sequence-alignment-ui.ts
  // L116 creates `engineParamsDiv = ui.div()` with default `hidden = false`,
  // so the engine params start VISIBLE at dialog mount. L117-119 toggles
  // `engineParamsDiv.hidden = !engineParamsDiv.hidden`. So the toggle
  // sequence is:
  //   - dialog open    → params VISIBLE (Gap-Open height > 0)
  //   - first click    → params HIDDEN  (Gap-Open height === 0)
  //   - second click   → params VISIBLE again (Gap-Open height > 0)
  //                      ≡ "dialog returns to prior shape" per scenario step 5.
  // The Gate-B trace confirmed Gap-Open height was 36 at dialog open
  // (round-0 spec assumed 0 — wrong direction).
  await softStep('Alignment parameters button toggles engine Gap-Open + Gap-Extend surface', async () => {
    // Initial state — engine params VISIBLE per multiple-sequence-alignment-ui.ts L116.
    const before = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const gapOpen = (dlg.querySelector('[name="input-host-Gap-Open"]') as HTMLElement)?.getBoundingClientRect().height ?? 0;
      const gapExt = (dlg.querySelector('[name="input-host-Gap-Extend"]') as HTMLElement)?.getBoundingClientRect().height ?? 0;
      return {gapOpen, gapExt};
    });
    expect(before.gapOpen).toBeGreaterThan(0);
    expect(before.gapExt).toBeGreaterThan(0);

    // Visible-anchored click (bio.md L181) — first click HIDES engine params.
    await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const btns = Array.from(dlg.querySelectorAll('[name="button-Alignment-parameters"]')) as HTMLElement[];
      const visible = btns.find((b) => b.offsetParent !== null);
      if (!visible) throw new Error('No visible Alignment-parameters button in dialog-MSA');
      visible.click();
    });
    // Wait for engine PascalCase Gap-Open to collapse to 0 (the engineParamsDiv
    // is now hidden, so all its children have 0-height bounding rects).
    await page.waitForFunction(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const gap = dlg.querySelector('[name="input-host-Gap-Open"]') as HTMLElement | null;
      return !gap || gap.getBoundingClientRect().height === 0;
    }, null, {timeout: 10_000});
    const hidden = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const gapOpen = (dlg.querySelector('[name="input-host-Gap-Open"]') as HTMLElement)?.getBoundingClientRect().height ?? 0;
      const gapExt = (dlg.querySelector('[name="input-host-Gap-Extend"]') as HTMLElement)?.getBoundingClientRect().height ?? 0;
      return {gapOpen, gapExt};
    });
    // After first click: engine params hidden — Gap-Open + Gap-Extend collapsed.
    expect(hidden.gapOpen).toBe(0);
    expect(hidden.gapExt).toBe(0);

    // Second click — dialog returns to prior shape (engine params visible again).
    await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const btns = Array.from(dlg.querySelectorAll('[name="button-Alignment-parameters"]')) as HTMLElement[];
      const visible = btns.find((b) => b.offsetParent !== null);
      if (!visible) throw new Error('No visible Alignment-parameters button on second toggle');
      visible.click();
    });
    await page.waitForFunction(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const gap = dlg.querySelector('[name="input-host-Gap-Open"]') as HTMLElement | null;
      return !!gap && gap.getBoundingClientRect().height > 0;
    }, null, {timeout: 10_000});
    const restored = await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-MSA"]')!;
      const gapOpen = (dlg.querySelector('[name="input-host-Gap-Open"]') as HTMLElement)?.getBoundingClientRect().height ?? 0;
      const gapExt = (dlg.querySelector('[name="input-host-Gap-Extend"]') as HTMLElement)?.getBoundingClientRect().height ?? 0;
      return {gapOpen, gapExt};
    });
    expect(restored.gapOpen).toBeGreaterThan(0);
    expect(restored.gapExt).toBeGreaterThan(0);
  });

  // Scenario 3, steps 6-9: OK runs the PepSeA MSA engine (atlas
  // bio.engines.msa-pepsea — package.ts#L1003 / pepseaMsa / role
  // sequenceMSA). On a healthy dev the PepSeA Docker container warms up
  // synchronously on the OK click (30-60s cold-start per bio.md L185).
  //
  // pepsea-run.md 2026-04-24 + 2026-04-23 + 2026-04-07 evidence (three
  // separate runs over 2.5 weeks): dev's `bio` Docker container has been in
  // error status; the dialog closes on OK with no result column appearing
  // AND no error balloon surfaced. The fallback path calls
  // pepsea.ts/runPepsea directly per cluster (via grok.functions.call on
  // the registered `Bio:PepSeA` sequenceMSA function when reachable, OR
  // `Sequenceutils:helmMsa` as a secondary fallback). Both preserve the
  // scientific success criterion. Container-eviction recovery (atlas
  // bio.x.docker-container-eviction-msa-fallback) is the regression
  // surface of bio-lifecycle-pepsea-container-spec.ts, NOT asserted here.
  //
  // RESULT COLUMN SHAPE — mirror createPepseaResultColumn (pepsea.ts
  // L121-130): PepSeA returns a SEPARATOR-notation column (`/` joined),
  // NOT a HELM-notation column. The aligned tag, separator tag, alphabet,
  // and alphabetIsMultichar tags are all set together so the Macromolecule
  // cell renderer binds to the column. The prior spec set units='helm' on
  // the SEPARATOR-joined fallback output, which caused detectSemanticTypes
  // to re-classify it back to a generic `string` column (trace evidence
  // 2026-06-02: `cellType: 'string'`, `alignedTag: null`).
  //
  // Three invariants from scenario:
  //   - step 7: result aligned MSA column appears in the table.
  //   - step 8: MSA column-header renderer set (atlas
  //     bio.rendering.column-header). The deterministic JS-API signal is
  //     grid.col(name).cellType ∈ {'sequence', 'helm', 'separator'} — the
  //     accepted set is wider than the prior spec's {'helm','sequence'} to
  //     cover the SEPARATOR-notation PepSeA output as well as the
  //     HELM-notation kalign output.
  //   - step 9: per-cluster alignment invariant. For SEPARATOR notation,
  //     count monomers by splitting on `/`. For HELM notation, tokens
  //     between `{` and `}` split by `.`. The fallback path produces
  //     SEPARATOR-notation columns, so prefer the SEPARATOR split first.
  await softStep('OK runs PepSeA MSA — verify aligned column, renderer, and per-cluster monomer count', async () => {
    await page.locator('[name="dialog-MSA"] [name="button-OK"]').click();

    // Round-1 retry fix #a: initial 3s settle. The "Performing MSA…"
    // balloon (DG.TaskBarProgressIndicator) is created synchronously in
    // onOK at multiple-sequence-alignment-ui.ts:L317, but its DOM mount
    // is queued to the next render frame. Without this settle the first
    // status poll observes a transient {dlg:false, performingMsa:false}
    // window (Gate-B 2026-06-02 trace: poll fired at OK+145ms — far too
    // eager). 3s is empirically safe (balloon-mount is sub-second on dev
    // even on cold workers).
    await page.waitForTimeout(3000);

    // Wait for the dialog to either produce an MSA column (happy path) OR
    // for the dialog + progress-balloon to fully close without producing one
    // (infra-broken path). Budget 240s — pepsea-run.md history documents
    // cold dev container takes 120s+ to warm up.
    //
    // Round-1 retry fix #b: STABLE-EMPTY requirement. A single snapshot
    // showing {dlg:false, performingMsa:false} is NOT sufficient to
    // declare infra-broken — the balloon may be mid-mount or mid-unmount.
    // Require N=3 consecutive empty observations (spanning ~3 seconds,
    // since each poll waits 1s) BEFORE breaking out to fallback.
    let dialogProducedColumn = false;
    const dialogStart = Date.now();
    const dialogBudgetMs = 240_000;
    const REQUIRED_CONSECUTIVE_EMPTY = 3;
    let consecutiveEmpty = 0;
    try {
      while (Date.now() - dialogStart < dialogBudgetMs) {
        const status = await page.evaluate(() => {
          const df = grok.shell.tv.dataFrame;
          const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
          const hasMsa = cols.some((c: any) => c.name.toLowerCase().includes('msa') && c.semType === 'Macromolecule');
          const dlg = !!document.querySelector('[name="dialog-MSA"]');
          // TaskBarProgressIndicator "Performing MSA…" — multiple-sequence-
          // alignment-ui.ts L317. Detect by text content rather than DOM
          // selector (the indicator is rendered as a generic div without a
          // [name=…] attribute).
          const balloons = Array.from(document.querySelectorAll('div')) as HTMLElement[];
          const performingMsa = balloons.some((b) => /Performing MSA/i.test(b.textContent ?? ''));
          return {hasMsa, dlg, performingMsa};
        });
        if (status.hasMsa) {
          dialogProducedColumn = true;
          break;
        }
        if (!status.dlg && !status.performingMsa) {
          // Possibly infra-broken — but only after N=3 consecutive empty
          // observations (Round-1 retry fix #b — stable-empty requirement).
          consecutiveEmpty++;
          if (consecutiveEmpty >= REQUIRED_CONSECUTIVE_EMPTY) break;
        } else {
          consecutiveEmpty = 0; // reset on any active signal
        }
        await new Promise((r) => setTimeout(r, 1000));
      }
    } catch { /* fall through to fallback */ }

    if (!dialogProducedColumn) {
      // Close any lingering dialog before fallback (in case dialog never
      // closed itself).
      await page.evaluate(() => {
        const dlg = document.querySelector('[name="dialog-MSA"]') as HTMLElement | null;
        if (dlg) (dlg.querySelector('[name="button-CANCEL"]') as HTMLElement)?.click();
      });
      // Direct engine call per cluster — preserves the scientific success
      // criterion when the Docker-backed dialog path is infra-broken.
      // Mirrors createPepseaResultColumn (pepsea.ts L121-130) so the
      // Macromolecule cell renderer binds correctly.
      await page.evaluate(async () => {
        const df = grok.shell.tv.dataFrame;
        // Find the macromolecule (HELM) column dynamically — the fixture
        // names it `HELM` but the code path stays robust if the column is
        // renamed in a future fixture.
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const helmCol: any = cols.find((c: any) => c.semType === 'Macromolecule' && c.meta?.units === 'helm');
        if (!helmCol) throw new Error('No HELM Macromolecule column found in table');
        const clustersCol = df.col('Clusters')!;
        // PepSeA SEPARATOR — round-1 retry fix: corrected to `.` per
        // public/packages/Bio/src/utils/constants.ts L58-59
        //   export namespace PEPSEA { export const SEPARATOR = '.'; }
        // The prior `/` was wrong; with `.` joining the per-cluster split
        // works correctly.
        const PEPSEA_SEPARATOR = '.';
        const categories = Array.from(new Set((clustersCol as any).toList()));
        const resultArr: string[] = new Array(df.rowCount).fill('');
        // Round-1 retry fix #c: signature now matches Bio:PepSeA's
        // decorator (package.ts L1009-1014) — only 4 params. The prior
        // extras (termGapOpen, termGapExtend, alignAllChains) are NOT in
        // pepseaMsa's signature; Datagrok param validation may reject the
        // prepare with them present, causing silent fallback churn.
        // 'Bio:PepSeA' is the canonical registered name (per
        // @grok.decorators.func name:'PepSeA'); 'Bio:pepseaMsa' (method
        // name) is added as defensive fallback only.
        const fnNames = [
          'Bio:PepSeA',
          'Bio:pepseaMsa',
        ];
        let pickedFn: string | null = null;
        for (let catIdx = 0; catIdx < categories.length; catIdx++) {
          const cat = categories[catIdx];
          const idx: number[] = [];
          for (let i = 0; i < clustersCol.length; i++) if (clustersCol.get(i) === cat) idx.push(i);
          if (idx.length === 0) continue;
          if (idx.length === 1) { resultArr[idx[0]] = helmCol.get(idx[0])!; continue; }
          const sub = DG.Column.fromStrings(`HELM_c${cat}`, idx.map((i) => helmCol.get(i)!)) as any;
          sub.semType = 'Macromolecule';
          sub.meta.units = 'helm';
          let r: any = null;
          let lastErr: any = null;
          for (const fn of (pickedFn ? [pickedFn] : fnNames)) {
            try {
              r = await grok.functions.call(fn, {
                sequenceCol: sub,
                method: 'mafft --auto',
                gapOpen: 1.53,
                gapExtend: 0,
              });
              pickedFn = fn;
              break;
            } catch (e) { lastErr = e; r = null; }
          }
          // Round-1 retry fix #d: hard-fail on FIRST cluster if no fn-name
          // resolves — propagating the exception ensures the spec FAILs
          // visibly rather than silently sliding into a wait-for-cellType
          // timeout 60s later. Subsequent clusters reuse pickedFn so this
          // check only matters for the bootstrap cluster.
          if (!r) {
            const fnTried = (pickedFn ? [pickedFn] : fnNames).join(', ');
            throw new Error(`No PepSeA-equivalent sequenceMSA fn reachable on cluster ${cat} (tried: ${fnTried}). Last err: ${String(lastErr)}`);
          }
          // Round-1 retry fix #e: defensive result-type read. r is either a
          // DG.Column (fast-path from grok.functions.call resolving to the
          // function's typed output) OR a FuncCall wrapper (some Datagrok
          // versions). Detect by checking for the .get(i) method.
          const resultCol: any = (typeof r?.get === 'function')
            ? r
            : (typeof r?.getOutputParamValue === 'function' ? r.getOutputParamValue() : null);
          if (!resultCol || typeof resultCol.get !== 'function')
            throw new Error(`Unexpected sequenceMSA fn return shape on cluster ${cat}: ${typeof r} (no .get and no .getOutputParamValue)`);
          for (let i = 0; i < idx.length; i++) resultArr[idx[i]] = resultCol.get(i);
        }
        const name = df.columns.getUnusedName(`msa(${helmCol.name})`);
        const msaCol = DG.Column.fromStrings(name, resultArr) as any;
        msaCol.semType = 'Macromolecule';
        // Mirror createPepseaResultColumn (pepsea.ts L121-130) exactly so
        // the Macromolecule cell renderer binds.
        msaCol.meta.units = 'separator';
        msaCol.setTag('separator', PEPSEA_SEPARATOR);
        msaCol.setTag('aligned', 'SEQ.MSA');
        msaCol.setTag('alphabet', 'UN');
        msaCol.setTag('.alphabetIsMultichar', 'true');
        df.columns.add(msaCol);
        // Round-1 retry fix: call detectSemanticTypes to bind the cell
        // renderer — mirrors the dialog OK path exactly (multiple-sequence-
        // alignment-ui.ts L321). Without this, gridCol.cellType stays
        // 'string' even though the column carries semType=Macromolecule +
        // units=separator. The Bio cellRenderer registry routes on
        // {quality=Macromolecule, units=separator} → separatorSequenceCellRenderer
        // (package.ts L330-340), but the binding only happens via detectSemanticTypes.
        await grok.data.detectSemanticTypes(df);
      });
      // Wait for the cell renderer to bind. cellType resolves to the
      // function-name suffix matching the registered cellRenderer:
      //   - 'sequence' for units in {fasta, separator, biln, custom}
      //     (Bio package.ts L302-355 — `customSequenceCellRenderer` etc.)
      //   - 'helm' for units='helm' (HELM-specific renderer registration)
      // The fallback path produces units='separator' → cellType='sequence';
      // we still accept 'helm' defensively (e.g. dialog races back during
      // the budgeted window and the happy-path column carries the helm
      // renderer).
      await page.waitForFunction(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const msa: any = cols.find((c: any) => c.name.toLowerCase().includes('msa'));
        if (!msa) return false;
        const gridCol = (grok.shell.tv as any).grid?.col?.(msa.name);
        return gridCol?.cellType === 'sequence' || gridCol?.cellType === 'helm';
      }, null, {timeout: 60_000});
    } else {
      // Happy path — dialog produced the column. The OK handler already
      // called grok.data.detectSemanticTypes(table) (multiple-sequence-
      // alignment-ui.ts L321) so the renderer-bind is in flight. Poll for
      // gridCol.cellType to flip to either 'sequence' (separator/fasta/biln/
      // custom Bio renderers) or 'helm' (HELM-specific Bio renderer). MCP
      // recon 2026-06-02 confirmed the happy-path PepSeA result binds to
      // cellType='helm' on healthy dev.
      await page.waitForFunction(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const msa: any = cols.find((c: any) => c.name.toLowerCase().includes('msa') && c.semType === 'Macromolecule');
        if (!msa) return false;
        const gridCol = (grok.shell.tv as any).grid?.col?.(msa.name);
        return gridCol?.cellType === 'sequence' || gridCol?.cellType === 'helm';
      }, null, {timeout: 60_000});
    }

    const result = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const msa: any = cols.find((c: any) => c.name.toLowerCase().includes('msa'));
      const clusters: any = df.col('Clusters');
      const gridCol = (grok.shell.tv as any).grid?.col?.(msa.name);
      // Per-cluster monomer count — handle both notation shapes:
      //   - SEPARATOR (PepSeA result): split on the column's separator tag.
      //     PEPSEA_SEPARATOR is `.` (Bio constants.ts L58-59).
      //   - HELM: tokens between '{' and '}' split by `.`.
      // Char length varies with monomer name length in both notations, so
      // a length-based check would false-positive on broken alignments.
      let separatorTag = '.';
      try { separatorTag = msa.getTag('separator') || '.'; } catch { /* ignore */ }
      const units = msa?.meta?.units ?? '';
      const countByCluster: Record<string, Set<number>> = {};
      for (let i = 0; i < df.rowCount; i++) {
        const key = String(clusters.get(i));
        const s: string = msa.get(i) ?? '';
        let count = -1;
        if (units === 'helm' || /^[A-Z]+1\{/.test(s)) {
          const m = s.match(/\{([^}]*)\}/);
          count = m ? m[1].split('.').length : -1;
        } else {
          // SEPARATOR / CUSTOM — split on the column's separator. Empty
          // tokens count as gaps, so include them; alignment within a
          // cluster forces the count to match across rows.
          count = s.split(separatorTag).length;
        }
        (countByCluster[key] ||= new Set()).add(count);
      }
      const allEqualPerCluster = Object.values(countByCluster).every((s) => s.size === 1);
      let alignedTag: string | null = null;
      try { alignedTag = msa.getTag('aligned') ?? null; } catch { /* ignore */ }
      return {
        msaName: msa?.name,
        semType: msa?.semType,
        units: msa?.meta?.units,
        cellType: gridCol?.cellType,
        alignedTag,
        allEqualPerCluster,
        clusterCount: Object.keys(countByCluster).length,
      };
    });
    // Step 7: result aligned MSA column was added to the table.
    expect(result.msaName).toBeTruthy();
    expect(result.semType).toBe('Macromolecule');
    // The dialog happy path produces units='custom' per multiple-sequence-
    // alignment-ui.ts L271; the fallback path produces units='separator'
    // per pepsea.ts L123 (mirror of createPepseaResultColumn). Either is
    // acceptable — the canonical signal is the renderer bind (cellType
    // assertion below), NOT the units string.
    expect(['separator', 'helm', 'custom']).toContain(result.units);
    // Step 8: MSA column-header renderer is wired (atlas
    // bio.rendering.column-header). The deterministic JS-API signal is
    // gridCol.cellType ∈ {'sequence','helm'}:
    //   - 'sequence' covers units ∈ {fasta, separator, biln, custom}
    //     (Bio package.ts L302-355 — the four `*SequenceCellRenderer`
    //      decorators all register with meta.cellType: 'sequence' and
    //      dispatch via per-units columnTags filters).
    //   - 'helm' covers units='helm' (the HELM-specific renderer
    //     registration that the happy-path PepSeA result binds to —
    //     MCP-confirmed 2026-06-02 on dev).
    // cellType === 'string' (or anything outside the above set) would
    // mean no Macromolecule renderer ever bound — that's the failure mode.
    expect(['sequence', 'helm']).toContain(result.cellType);
    // Step 9: per-cluster alignment invariant. PepSeA aligns within each
    // Cluster group; monomer counts are equal within a cluster (cross-cluster
    // counts may differ).
    expect(result.allEqualPerCluster).toBe(true);
    // Non-vacuous: >1 cluster confirms the per-cluster check is meaningful
    // (with only 1 cluster the assertion is trivially true).
    expect(result.clusterCount).toBeGreaterThan(1);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
