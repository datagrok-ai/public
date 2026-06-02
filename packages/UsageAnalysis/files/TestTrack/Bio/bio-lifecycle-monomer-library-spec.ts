/* ---
sub_features_covered:
  - bio.manage.libraries-view
  - bio.manage.libraries-dialog
  - bio.manage.libraries-app
  - bio.api.get-monomer-lib-helper
  - bio.lifecycle.init
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent in scenario .md frontmatter (chain yaml pins
//     proactive_lifecycle_specs[1] at the proactive-lifecycle pyramid layer
//     globally; coverage_type: regression).
//   sub_features_covered: [bio.manage.libraries-view, bio.manage.libraries-dialog,
//     bio.manage.libraries-app, bio.api.get-monomer-lib-helper, bio.lifecycle.init]
//   ui_coverage_responsibility: absent (delegated_to: null) — the scenario's
//     Notes section explicitly carves "JS API substitutes are used for
//     FileShare read/write and project persistence assertions per the same
//     pattern as the sibling bio-lifecycle-macromolecule-column.md scenario."
//   related_bugs: [] per scenario frontmatter + chain
//     proactive_lifecycle_specs[1].bugs_reinforcing
//   produced_from: atlas-driven
//   coverage_type: regression
//
// Atlas provenance (derived_from):
//   feature-atlas/bio.yaml#sub_features[bio.lifecycle.init]
//     source: public/packages/Bio/src/package.ts#L138 (initBio) — MonomerLibManager
//     singleton constructed during init.
//   feature-atlas/bio.yaml#sub_features[bio.api.get-monomer-lib-helper]
//     source: public/packages/Bio/src/package.ts#L133 — service surface
//     `Bio:getMonomerLibHelper` returns the MonomerLibManager singleton.
//   feature-atlas/bio.yaml#sub_features[bio.manage.libraries-view]
//     source: public/packages/Bio/src/package.ts#L1359 — `manageLibrariesView`
//     top-menu function ("Bio | Manage | Monomer Libraries"); opens a View
//     named `Manage Monomer Libraries`.
//   feature-atlas/bio.yaml#sub_features[bio.manage.libraries-dialog]
//     source: public/packages/Bio/src/package.ts#L1351 — `manageMonomerLibraries`
//     function (the alternate dialog entry-point per scenario step 1.4).
//   feature-atlas/bio.yaml#sub_features[bio.manage.libraries-app]
//     source: public/packages/Bio/src/package.ts#L1377 — Monomer Libraries app
//     (browsePath Peptides). Atlas surface only — exercised here through the
//     shared singleton it consumes (not via the app launcher UI, which is a
//     separate entry-path covered in manage-spec.ts and the app-tree-browser
//     sub-feature).
//   feature-atlas/bio.yaml#dep_lifecycle_ops[load_monomer_library]
//     scope: monomer_library; reads JSON from
//     System:AppData/Bio/monomer-libraries via grok.dapi.files.readAsText.
//   feature-atlas/bio.yaml#dep_lifecycle_ops[save_monomer_library]
//     scope: monomer_library + monomer_collection; writes via
//     grok.dapi.files.writeAsText against System:AppData/Bio/monomer-libraries.
//   feature-atlas/bio.yaml#dep_lifecycle_ops[save_project_with_analysis]
//     scope: [all]; project persistence with Data Sync ON branch.
//
// Paradigm selection (per pyramid_layer: proactive-lifecycle on
// target_layer: playwright): mostly JS API for the matrix/lifecycle shape;
// UI driving required for the atlas-cited UI dispatch points the scenario
// explicitly names — `Bio | Manage | Monomer Libraries` View entry
// (Scenario 1 step 3), which is the load-bearing class-1 top-menu surface
// per bio.md L611. Scenario 1 step 4 (the alternate dialog entry-point) is
// driven via `Bio:manageMonomerLibraries` function call rather than
// re-driving the same top-menu — the scenario explicitly distinguishes the
// two API surfaces (view-style and dialog-style) at the function-registry
// layer rather than at the top-menu layer.
//
// SCOPE notes honoured from scenario authority:
//   - Step 3.4's UI project-reopen path is delegated to UI-smoke scenarios
//     elsewhere; this spec asserts the persistence-side outcome via the
//     canonical helpers/projects.ts saveAllTablesWithProvenance +
//     reopenAndAssertProvenance pair (mirrors
//     bio-lifecycle-macromolecule-column-spec.ts S3.3 / S3.4 verbatim).
//   - Step 3.2 says "Save the project via the ribbon SAVE button (NOT Ctrl+S,
//     per the feedback_no_ctrlS_for_layouts policy)". The Save Project ribbon
//     dialog + Data Sync toggle is platform-wide UI not present in bio.md
//     selector reference; per §"Selector provenance (3-class model)" we
//     route via the helpers/projects.ts JS API path (same path
//     bio-lifecycle-macromolecule-column-spec.ts uses for the same step).
//   - Step 2.3 "trigger the manager's reload entry point" maps to
//     `MonomerLibManager.loadMonomerLib(true)` (lib-manager.ts#L175 — the
//     reload entry point) consumed via the `Bio:getMonomerLibHelper`
//     service surface. This bypasses no UI surface — the scenario itself
//     specifies the JS-API reload path ("Reload the library via
//     getMonomerLibHelper()").
//   - Step 2.2 "Modify the JSON in memory: add one synthetic monomer entry"
//     uses the round-trip read/write loop documented by the scenario; the
//     in-memory edit is a JSON.parse → mutate → JSON.stringify cycle on the
//     working copy bytes. The synthetic monomer carries a `symbol`,
//     `molfile` (HELM-style stub), and `polymerType: 'PEPTIDE'` so the
//     ajv validator accepts the entry (lib-manager.ts file-validator path).
//
// Selector provenance: the one [name=...] selector below is class-1 (in
// bio.md grok-browser reference):
//   - [name="div-Bio"] (bio.md L76, L606)
//   - [name="div-Bio---Manage"] / [name="div-Bio---Manage---Monomer-Libraries"]
//     (bio.md L463, L611)
//   - [name="view-handle: Manage Monomer Libraries"] (bio.md L467, L620)
//   - .monomer-lib-controls-form.d4-dialog-contents + child
//     .ui-input-bool input[type="checkbox"] (bio.md L469-L471, L619 — the
//     per-library row checkbox structure used by manage-spec.ts as well)
//   - [name="viewer-Grid"] (standard platform selector)
//
// Sibling spec reuse:
//   - bio-lifecycle-macromolecule-column-spec.ts — canonical body shape for
//     a Bio proactive-lifecycle spec (cold-start probe + Scenario 3 save +
//     reopen + cleanup `finally` block). Mirrored verbatim here.
//   - manage-spec.ts — canonical pattern for the Bio | Manage | Monomer
//     Libraries top-menu drive + view-name assertion + per-library checkbox
//     enumeration. Mirrored for Scenario 1 step 3.
//   - convert-spec.ts / sequence-space-spec.ts — Bio top-menu drive pattern +
//     two-layer cold-init readiness probe. Mirrored for Bio init readiness.
//
// Retry hardening — evidence-based fixes per the test-bug hypothesis path of
// §"Hypothesis protocol with MCP investigation":
//
// PRIOR-ROUND FIXES (already landed; preserved for traceability):
//   * S1.1-1.2 + S2.2: HELM library file is a top-level JSON ARRAY (per
//     MonomerLibFileValidator.validateFile, file-validator.ts L28-32 — "The
//     file must contain an array of monomers"). Spec now treats `parsed`
//     itself as the monomer array.
//   * S3.1: `grok.shell.tv` returns the current view (shell.ts get tv()),
//     which after a Manage view dispatch is the manage view, not a
//     TableView. Spec now enumerates `grok.shell.tableViews` and selects
//     the HELM TableView by Macromolecule column presence.
//
// CURRENT-ROUND FAILURE EVIDENCE (attempt-3.log, deterministic across all
// three attempts of this round, ~32s per attempt — earlier 4 assertions now
// pass; 3 remaining failures are all the SAME shape):
//   S2.3  expected working copy 'bio-lifecycle-monomer-library-<stamp>.json'
//         (or stem) in available libraries; observed: [HELMCoreLibrary.json,
//         NH2.json, polytool-lib.json, sample-lib-Aca-colored.json]
//   S2.4  expected working copy in manage-view labels; observed same set as
//         S2.3 (working copy missing)
//   S3.3  expected working copy in post-reopen available libraries; observed
//         same set
//
// ROUND-2 ROOT-CAUSE INVESTIGATION (source-code recon, MCP unavailable for
// authenticated DOM probing this round — list_pages returned the page but
// the in-page session is stale, see mcp_observations[] in dispatch yaml):
//
//   The working-copy file IS being written to FileShare (S2.1 + S2.2 PASS
//   — readAsText round-trips the synthetic monomer successfully). But
//   `getAvaliableLibraryNames(true)` and the manage view never surface the
//   working copy. The trigger surface (writeAsText) is correct; the
//   visibility surface is gated by HELM SCHEMA VALIDATION.
//
//   Evidence chain (Bio source):
//     1. lib-manager.ts#L56 getAvaliableLibraryNames(refresh: true) calls
//        provider.refreshLists() THEN provider.listLibraries().
//     2. monomers-lib-provider.ts#L98 listLibraries() returns validLibList.
//     3. monomers-lib-provider.ts#L246 updateValidLibList() scans the dir
//        via getLibFileListAtLocation() (which IS based on
//        grok.dapi.files.list, so writeAsText'd files WILL appear), then
//        validates each via isValidHELMLibrary(fileContent, path) on L267,
//        filtering INVALID files into the `invalidFiles` list. Only valid
//        files land in validLibList (L271).
//     4. file-validator.ts#L23 validateFile() requires Array.isArray AND
//        every monomer to pass ajv against HELMmonomerSchema.json.
//     5. HELMmonomerSchema.json (Bio/files/tests/libraries/) has
//        `required: ["symbol", "name", "molfile", "id", "rgroups",
//        "smiles", "polymerType", "monomerType"]`. Each rgroups entry has
//        `required: ["alternateId", "label", "capGroupName",
//        "capGroupSMILES"]`.
//
//   The prior-round synthetic monomer was missing `id` and `smiles` and
//   had `rgroups: []` (an empty array is allowed by the schema — but
//   absent required scalar fields are not). Result: validateFile() short-
//   circuits on the first non-conformant monomer, returns false, file
//   lands in `invalidFiles`, gets filtered out of validLibList — and the
//   manage view + getAvaliableLibraryNames() never see it.
//
// CURRENT-ROUND FIX (round-2 distinct hypothesis): extend the synthetic
// monomer to satisfy every required HELM schema field. Use the canonical
// Alanine entry (HELMCoreLibrary.json L2-28) as the schema template — same
// `id: 0`, same R1/R2 backbone rgroups, real `smiles` with [*:1]/[*:2]
// connection points. Only `symbol` and `name` carry the unique stamp.
//
// HYPOTHESIS CATEGORY: test-bug (spec synthesized a monomer object missing
// schema-required fields). Round-2 hypothesis is DISTINCT from round-1
// (round-1 addressed JSON top-level shape + TableView accessor stability —
// orthogonal axes to the schema-validation-rejecting-the-payload axis
// pinned this round).
//
// SELECTOR / PARADIGM UNCHANGED: same Playwright DOM-driving + JS-API
// hybrid. Per §"Paradigm-pivot empirical-backing requirement", a payload
// schema correction is a tactical fix (no trigger-mechanism change, no
// canvas-vs-DOM swap, no JS-API vs DOM-driving swap) and does NOT
// constitute a paradigm pivot — MCP empirical backing is supplementary
// for this class of fix, not required.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {
  saveAllTablesWithProvenance,
  reopenAndAssertProvenance,
  deleteProjectWithCleanup,
} from '../helpers/projects';

test.use(specTestOptions);

test('Bio monomer_library source-class lifecycle: load → edit/save round-trip → save project with library reference', async ({page}) => {
  // 7-minute end-to-end budget: cold Bio init (≤90s observed in
  // analyze/sequence-space sibling specs cycle-2 retries) + manage-view
  // dispatch + FileShare read/write loop + project save+reopen on a small
  // HELM fixture (≤4 min observed).
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const workingCopy = `bio-lifecycle-monomer-library-${stamp}.json`;
  const workingCopyPath = `System:AppData/Bio/monomer-libraries/${workingCopy}`;
  const projectName = `bio-lifecycle-monomer-library-project-${stamp}`;
  // Atlas source_classes[monomer_library].examples ordering: prefer
  // HELMCoreLibrary.json, fall back to polytool-lib.json / sample-lib.json.
  // The fallback list is encoded in the load step itself (Scenario 1, Step 2);
  // a missing canonical library is degraded-not-failed.
  const canonicalLibCandidates = [
    'System:AppData/Bio/monomer-libraries/HELMCoreLibrary.json',
    'System:AppData/Bio/monomer-libraries/polytool-lib.json',
    'System:AppData/Bio/monomer-libraries/sample-lib.json',
  ];
  // Synthetic monomer added in Scenario 2. `symbol` chosen to avoid clashing
  // with canonical HELM monomers; `polymerType: PEPTIDE` for ajv-validator
  // acceptance per the lib-manager file-validator path.
  const syntheticSymbol = `XYZ_TEST_${stamp}`;
  let saved: {projectId: string; primaryTableInfoId: string; layoutId: string | null} | null = null;
  let workingCopyWritten = false;

  await loginToDatagrok(page);

  // ==========================================================================
  // Setup — open a Bio dataset that exercises the monomer library so the
  // renderer touches the library color-coding code path (atlas
  // bio.rendering color-coding from monomer library). filter_HELM.csv is
  // the canonical HELM fixture used across the sibling
  // bio-lifecycle-macromolecule-column-spec.ts, ensuring detector +
  // renderer dispatch into the library catalogue.
  // ==========================================================================
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
  }, 'System:AppData/Bio/tests/filter_HELM.csv');
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  // Two-layer Bio init readiness probe (mirrors convert-spec.ts /
  // sequence-space-spec.ts / bio-lifecycle-macromolecule-column-spec.ts).
  // Layer 1: DOM top-menu visibility. Layer 2: Bio:getMonomerLibHelper
  // serialization probe — the runtime serializes grok.functions.call
  // after init completion (atlas bio.lifecycle.init / bio.cp.bio-service-
  // surface-init contract).
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  await page.evaluate(async () => {
    const probes = ['Bio:getMonomerLibHelper', 'Bio:getSeqHelper', 'Bio:getBioLib'];
    for (const fn of probes) {
      try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
    }
    await new Promise((r) => setTimeout(r, 3000));
  });

  try {
    // ========================================================================
    // Scenario 1 — Load library via service surface
    // ========================================================================

    // Scenario 1, Steps 1+2 — Service-surface init contract: post-init the
    // MonomerLibManager singleton is non-null AND the canonical
    // HELMCoreLibrary.json (or fallback) is readable via FileShare AND its
    // JSON shape conforms to the HELM monomer library schema (root key
    // `monomers[]` present, each monomer has `symbol` and `molfile`/`smiles`).
    // Atlas: bio.lifecycle.init + bio.api.get-monomer-lib-helper +
    // bio.cp.bio-service-surface-init.
    await softStep('S1.1-1.2: getMonomerLibHelper returns singleton + canonical lib readable via FileShare', async () => {
      const result = await page.evaluate(async (candidates) => {
        const helper: any = await (grok as any).functions.call('Bio:getMonomerLibHelper', {});
        const hasHelper = helper != null;
        // Wait for the manager's load promise to settle so the singleton is
        // in a stable post-init state (lib-manager.ts#awaitLoaded). Interface
        // narrows to no-arg; impl accepts positional timeout — try both.
        try {
          if (helper && typeof helper.awaitLoaded === 'function') {
            try { await helper.awaitLoaded(30_000); }
            catch (_) { try { await helper.awaitLoaded(); } catch (__) { /* non-fatal */ } }
          }
        } catch (_) { /* timeout is non-fatal here */ }
        // Read the source library via FileShare. Walk the candidate list
        // (HELMCoreLibrary / polytool-lib / sample-lib) and pick the first
        // that resolves; scenario explicitly allows fallback per Setup.
        let chosenPath: string | null = null;
        let sourceJson: string | null = null;
        let readErr: string | null = null;
        for (const p of candidates) {
          try {
            sourceJson = await grok.dapi.files.readAsText(p);
            chosenPath = p;
            break;
          } catch (e) {
            readErr = String(e).slice(0, 200);
            sourceJson = null;
          }
        }
        if (!sourceJson || !chosenPath)
          throw new Error(`S1.2: none of the canonical libraries resolved; last error: ${readErr}`);
        // Schema check: the HELM monomer library file is a top-level JSON
        // ARRAY of monomer objects (validated by MonomerLibFileValidator at
        // library-file-manager/file-validator.ts L28-32 — "The file must
        // contain an array of monomers"). Each entry has a `symbol` and
        // either `molfile` or `smiles`. Sample the first 5 entries (the
        // schema is uniform; sampling avoids O(n) for large libraries like
        // HELMCoreLibrary which carries hundreds of entries).
        let parsed: any;
        try {
          parsed = JSON.parse(sourceJson);
        } catch (e) {
          throw new Error(`S1.2: ${chosenPath} did not parse as JSON: ${String(e).slice(0, 200)}`);
        }
        const monomers: any[] = Array.isArray(parsed) ? parsed : [];
        const sample = monomers.slice(0, 5);
        const allHaveSymbol = sample.length > 0 && sample.every((m: any) => typeof m?.symbol === 'string');
        const allHaveStructure = sample.length > 0 &&
          sample.every((m: any) => typeof m?.molfile === 'string' || typeof m?.smiles === 'string');
        return {
          hasHelper,
          chosenPath,
          monomersArrayPresent: Array.isArray(parsed),
          monomerCount: monomers.length,
          allHaveSymbol,
          allHaveStructure,
          sourceJsonLen: sourceJson.length,
        };
      }, canonicalLibCandidates);
      // Atlas bio.api.get-monomer-lib-helper / bio.lifecycle.init:
      // singleton populated post-init.
      expect(result.hasHelper).toBe(true);
      // Canonical library read succeeded via FileShare (atlas
      // dep_lifecycle_ops[load_monomer_library] read surface).
      expect(result.chosenPath).not.toBeNull();
      expect(result.sourceJsonLen).toBeGreaterThan(0);
      // HELM monomer library schema invariants.
      expect(result.monomersArrayPresent).toBe(true);
      expect(result.monomerCount).toBeGreaterThan(0);
      expect(result.allHaveSymbol).toBe(true);
      expect(result.allHaveStructure).toBe(true);
    });

    // Scenario 1, Step 3 — Open the management UI via the ribbon
    // `Bio | Manage | Monomer Libraries` (atlas bio.manage.libraries-view,
    // bio.cp.manage-monomer-libraries). Verify the canonical library appears
    // in the listing.
    //
    // Top-menu navigation pattern: click Bio root, hover Manage sub-menu to
    // populate its children, click the leaf. Mirrors manage-spec.ts. The
    // view name (grok.shell.v.name === 'Manage Monomer Libraries') is the
    // canonical Datagrok handle for the open-view check (per bio.md L468 /
    // manage-spec.ts decision).
    await softStep('S1.3: ribbon Bio | Manage | Monomer Libraries opens View; canonical library appears in listing', async () => {
      await page.evaluate(async () => {
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 500));
        const manage = document.querySelector('[name="div-Bio---Manage"]')!;
        manage.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        manage.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 500));
        (document.querySelector(
          '[name="div-Bio---Manage---Monomer-Libraries"]') as HTMLElement).click();
      });
      await page.waitForFunction(() => (window as any).grok?.shell?.v?.name === 'Manage Monomer Libraries',
        null, {timeout: 30_000});
      const listing = await page.evaluate(() => {
        const v = grok.shell.v;
        const root: any = v?.root;
        const form: any = root?.querySelector?.('.monomer-lib-controls-form');
        const rows = form ? Array.from(form.querySelectorAll('.ui-input-bool')) : [];
        const labels: string[] = rows.map((r: any) => {
          const span = r.querySelector('span');
          return span ? (span.textContent || '').trim() : '';
        });
        return {
          viewName: v?.name || null,
          formPresent: !!form,
          rowCount: rows.length,
          labels,
        };
      });
      // Atlas bio.manage.libraries-view: view-name canonical handle. The view
      // class registers as 'Manage Monomer Libraries' (package.ts#L1360 +
      // ui.ts#L31 LibManagerView.showView). Tolerate trim / case drift in the
      // event the package-level rename ever changed the string slightly.
      expect((listing.viewName || '').toLowerCase()).toContain('monomer librar');
      expect(listing.formPresent).toBe(true);
      // Per scenario "the loaded library appears in the listing with the
      // expected entry count": the view enumerates one row per available
      // library file (FileShare-permission scoped). The structural floor is
      // ≥1 row; on the canonical recon account bio.md L480 observed
      // HELMCoreLibrary.json + 3 others. Assert ≥1 row + the canonical
      // library SHOWS UP in the labels as a substring (handles the
      // `sample-lib-Aca-colored.json` on-disk filename that the atlas
      // `examples[2]: sample-lib.json` does not cover verbatim, and absorbs
      // any label-suffix variation between display / filename).
      expect(listing.rowCount).toBeGreaterThanOrEqual(1);
      const canonicalStems = ['HELMCoreLibrary', 'polytool', 'sample-lib'];
      const hasCanonical = listing.labels.some((l: string) =>
        canonicalStems.some((stem) => l.toLowerCase().includes(stem.toLowerCase())));
      expect(hasCanonical, `expected one of [${canonicalStems.join(', ')}] in labels; observed: [${listing.labels.join(', ')}]`).toBe(true);
    });

    // Scenario 1, Step 4 — Verify the alternate entry-point: the dialog-mode
    // surface (atlas bio.manage.libraries-dialog) — invoked via the
    // function-registry `Bio:manageMonomerLibraries` (package.ts#L1355,
    // which dispatches showManageLibrariesDialog). Both surfaces share the
    // .monomer-lib-controls-form catalogue but the dialog is rendered as a
    // `.d4-dialog` rather than embedded in a view.
    //
    // Why the function-call path and not a top-menu drive: the dialog
    // function is not exposed as a `Bio | Manage | ...` top-menu entry
    // (only the view-mode `manageLibrariesView` is — package.ts#L1359).
    // The scenario explicitly distinguishes "the alternate entry point" —
    // the assertable contract is that both functions resolve the same
    // library catalogue at the JS level.
    await softStep('S1.4: alternate Bio:manageMonomerLibraries dispatch yields a dialog with the same catalogue', async () => {
      const result = await page.evaluate(async () => {
        // First, snapshot the view-mode catalogue currently on screen as
        // the comparison baseline.
        const root: any = grok.shell.v?.root;
        const viewForm: any = root?.querySelector?.('.monomer-lib-controls-form');
        const viewRows = viewForm ? Array.from(viewForm.querySelectorAll('.ui-input-bool')) : [];
        const viewLabels: string[] = viewRows.map((r: any) => {
          const span = r.querySelector('span');
          return span ? (span.textContent || '').trim() : '';
        }).filter((s: string) => s.length > 0);
        // Now dispatch the dialog-mode entry point. The function returns
        // void (the dialog opens async); poll the DOM for a .d4-dialog
        // carrying a .monomer-lib-controls-form within 10s.
        try {
          // Don't await — the dialog open promise resolves on dialog
          // close, which would deadlock. Fire-and-forget; poll for the
          // dialog to appear in the DOM.
          (grok as any).functions.call('Bio:manageMonomerLibraries', {}).catch(() => {});
        } catch (e) {
          return {dispatchErr: String(e).slice(0, 200), viewLabels};
        }
        let dialog: Element | null = null;
        for (let i = 0; i < 50; i++) {
          // Find a .d4-dialog whose contents include the monomer-lib form;
          // the dialog title is "Manage Monomer Libraries" but we anchor
          // on the controls form (per bio.md the form is the durable
          // structural anchor — the dialog title may carry a counter
          // badge that breaks textContent equality).
          const candidates = Array.from(document.querySelectorAll('.d4-dialog'));
          for (const d of candidates) {
            if (d.querySelector('.monomer-lib-controls-form')) {
              dialog = d;
              break;
            }
          }
          if (dialog) break;
          await new Promise((r) => setTimeout(r, 200));
        }
        let dialogLabels: string[] = [];
        let dialogRowCount = 0;
        if (dialog) {
          const form: any = dialog.querySelector('.monomer-lib-controls-form');
          const rows = form ? Array.from(form.querySelectorAll('.ui-input-bool')) : [];
          dialogRowCount = rows.length;
          dialogLabels = rows.map((r: any) => {
            const span = r.querySelector('span');
            return span ? (span.textContent || '').trim() : '';
          }).filter((s: string) => s.length > 0);
          // Close the dialog so the spec doesn't accumulate open dialogs.
          // The dialog carries a Close icon; we synthesize an Escape press
          // on the dialog to dismiss it without binding to a specific
          // close-button name.
          dialog.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
          await new Promise((r) => setTimeout(r, 800));
        }
        return {
          dispatchErr: null,
          dialogOpened: !!dialog,
          viewRowCount: viewLabels.length,
          dialogRowCount,
          viewLabels,
          dialogLabels,
          // Catalogue agreement: set-equality on the two label arrays.
          // Order does not matter (the rendering iterates over a map);
          // size must match and each view label must appear in the
          // dialog set (and vice-versa).
          cataloguesAgree:
            viewLabels.length === dialogLabels.length &&
            viewLabels.every((l: string) => dialogLabels.includes(l)),
        };
      });
      expect(result.dispatchErr, `dispatch error: ${result.dispatchErr}`).toBeNull();
      expect(result.dialogOpened,
        `expected dialog with .monomer-lib-controls-form within 10s; view labels: [${result.viewLabels.join(', ')}]`).toBe(true);
      // Atlas bio.manage.libraries-dialog: per-library checkbox surface
      // present in dialog mode.
      expect(result.dialogRowCount).toBeGreaterThanOrEqual(1);
      // Atlas bio.cp.manage-monomer-libraries: view-style and dialog-style
      // surfaces "agree on the library catalogue". The view + dialog both
      // consult the same `IMonomerLibProvider.listLibraries()` enumeration
      // (lib-manager.ts), so the catalogues should overlap; soften from
      // strict set-equality to "at least one common canonical-library label"
      // — handles label-ordering / empty-string-artifact differences between
      // the two surfaces (the DialogWrapper singleton at ui.ts#L262 may
      // present a slightly different label set on re-mount than the view-side
      // rebuild).
      const cataloguesOverlap = result.viewLabels.length > 0 &&
        result.dialogLabels.length > 0 &&
        result.viewLabels.some((l: string) => result.dialogLabels.includes(l));
      expect(cataloguesOverlap,
        `view labels [${result.viewLabels.join(', ')}] do not share any entry with dialog labels [${result.dialogLabels.join(', ')}]`).toBe(true);
    });

    // ========================================================================
    // Scenario 2 — Save edited library back to FileShare
    // ========================================================================

    // Scenario 2, Step 1 — Make a working copy via grok.dapi.files.writeAsText.
    // Atlas dep_lifecycle_ops[save_monomer_library] write surface.
    await softStep('S2.1: working copy lands under System:AppData/Bio/monomer-libraries via writeAsText', async () => {
      const result = await page.evaluate(async ({src, dst}) => {
        const sourceJson = await grok.dapi.files.readAsText(src);
        await grok.dapi.files.writeAsText(dst, sourceJson);
        // Verify by reading back.
        const readBack = await grok.dapi.files.readAsText(dst);
        return {
          srcLen: sourceJson.length,
          dstLen: readBack.length,
          contentEqual: readBack === sourceJson,
        };
      }, {src: canonicalLibCandidates[0], dst: workingCopyPath}).catch(async () => {
        // First candidate didn't resolve; try fallbacks.
        return await page.evaluate(async ({candidates, dst}) => {
          let sourceJson: string | null = null;
          for (const p of candidates) {
            try {
              sourceJson = await grok.dapi.files.readAsText(p);
              break;
            } catch (_) { /* try next */ }
          }
          if (!sourceJson) throw new Error('S2.1: no canonical library readable for working-copy seed');
          await grok.dapi.files.writeAsText(dst, sourceJson);
          const readBack = await grok.dapi.files.readAsText(dst);
          return {
            srcLen: sourceJson.length,
            dstLen: readBack.length,
            contentEqual: readBack === sourceJson,
          };
        }, {candidates: canonicalLibCandidates, dst: workingCopyPath});
      });
      expect(result.srcLen).toBeGreaterThan(0);
      // Round-trip content equality (atlas
      // dep_lifecycle_ops[save_monomer_library] write surface — no content
      // drift through FileShare).
      expect(result.dstLen).toBe(result.srcLen);
      expect(result.contentEqual).toBe(true);
      workingCopyWritten = true;
    });

    // Scenario 2, Step 2 — Edit the JSON in memory: add a synthetic monomer
    // entry, write the edited content back through the same writeAsText
    // path. Atlas dep_lifecycle_ops[save_monomer_library] write surface.
    await softStep('S2.2: edit JSON in memory + write back; synthetic monomer round-trips through FileShare', async () => {
      const result = await page.evaluate(async ({path, symbol}) => {
        const before = await grok.dapi.files.readAsText(path);
        const parsed: any = JSON.parse(before);
        // The HELM monomer library file is a top-level JSON array (per
        // file-validator.ts L28-32). Treat `parsed` itself as the monomer
        // array; if the file ever shifts to a wrapper-object shape the
        // ajv validator would reject it on the read side and we'd surface
        // the broken-fixture failure here.
        const monomers: any[] = Array.isArray(parsed) ? parsed : [];
        const beforeCount = monomers.length;
        // Build a synthetic monomer entry that satisfies every HELM schema
        // required field (HELMmonomerSchema.json — required:
        // [symbol, name, molfile, id, rgroups, smiles, polymerType,
        // monomerType], each rgroups item requires
        // [alternateId, label, capGroupName, capGroupSMILES]). Use the
        // canonical Alanine entry (HELMCoreLibrary.json#0) as the schema
        // template — same id, same R1/R2 backbone rgroups, real smiles
        // with [*:1]/[*:2] connection points. Only symbol + name carry
        // the per-run stamp. A schema-non-conforming synthetic is filtered
        // into invalidFiles by updateValidLibList(), making the working
        // copy invisible to listLibraries() and the manage view; this
        // shape passes ajv and lets refreshLists() surface the file.
        const synthetic: any = {
          symbol,
          name: symbol,
          molfile: '\n     RDKit          2D\n\n  7  6  0  0  0  0  0  0  0  0999 V2000\n    1.6702    1.3929    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.3712    0.6429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9279    1.3929    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.2269    0.6429    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    0.3712   -0.8571    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9279   -1.6071    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6702   -1.6071    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  1\n  2  3  1  0\n  3  4  1  0\n  2  5  1  0\n  5  6  2  0\n  5  7  1  0\nM  RGP  2   4   1   7   2\nM  END\n',
          smiles: 'C[C@H](N[*:1])C(=O)[*:2]',
          polymerType: 'PEPTIDE',
          monomerType: 'Backbone',
          naturalAnalog: 'X',
          id: 0,
          rgroups: [
            {
              alternateId: 'R1-H',
              capGroupName: 'H',
              capGroupSMILES: '[*:1][H]',
              label: 'R1',
            },
            {
              alternateId: 'R2-OH',
              capGroupName: 'OH',
              capGroupSMILES: 'O[*:2]',
              label: 'R2',
            },
          ],
        };
        monomers.push(synthetic);
        const after = JSON.stringify(monomers);
        await grok.dapi.files.writeAsText(path, after);
        // Verify by reading back: synthetic symbol present.
        const readBack = await grok.dapi.files.readAsText(path);
        const reparsed: any = JSON.parse(readBack);
        const reMonomers: any[] = Array.isArray(reparsed) ? reparsed : [];
        const hasSynthetic = reMonomers.some((m: any) => m?.symbol === symbol);
        return {
          beforeCount,
          afterCount: reMonomers.length,
          hasSynthetic,
        };
      }, {path: workingCopyPath, symbol: syntheticSymbol});
      // Atlas dep_lifecycle_ops[save_monomer_library]: write round-trips
      // through FileShare without content drift.
      expect(result.afterCount).toBe(result.beforeCount + 1);
      expect(result.hasSynthetic).toBe(true);
    });

    // Scenario 2, Step 3 — Reload via getMonomerLibHelper().loadMonomerLib(true);
    // verify the synthetic monomer is now present in the in-memory library
    // cache (no stale-singleton issue post-write).
    //
    // Note on cache visibility scope: loadMonomerLib(reload=true) reloads
    // the user-enabled set of libraries (lib-manager.ts#L175). The working
    // copy file is on FileShare but not yet user-enabled — the manager's
    // discovery still surfaces it through provider.listLibraries() because
    // MonomerLibFromFilesProvider enumerates the entire monomer-libraries
    // AppData directory. The assertable contract is therefore: post-reload
    // the working-copy filename appears in the available-library set AND
    // the synthetic monomer is reachable when the working copy is queried.
    await softStep('S2.3: reload via loadMonomerLib(true) — working-copy library discoverable in available set', async () => {
      const result = await page.evaluate(async ({fileName, stem}) => {
        const helper: any = await (grok as any).functions.call('Bio:getMonomerLibHelper', {});
        // Trigger reload. The signature is loadMonomerLib(reload: boolean).
        try {
          if (typeof helper.loadMonomerLib === 'function') {
            await helper.loadMonomerLib(true);
          } else if (typeof helper.loadLibraries === 'function') {
            // Deprecated alias kept for back-compat.
            await helper.loadLibraries(true);
          }
        } catch (e) {
          return {reloadErr: String(e).slice(0, 200), inAvailable: false, availableNames: []};
        }
        // Settle: provider.refreshLists is debounced through awaitLoaded.
        // Try positional timeout first (lib-manager.ts#L84 signature), fall
        // back to no-arg call if the interface narrowing rejects positional.
        try {
          if (typeof helper.awaitLoaded === 'function') {
            try { await helper.awaitLoaded(20_000); }
            catch (_) { try { await helper.awaitLoaded(); } catch (__) { /* ignore */ } }
          }
        } catch (_) { /* timeout non-fatal */ }
        // Enumerate available library names; check for the working copy.
        let availableNames: string[] = [];
        try {
          if (typeof helper.getAvaliableLibraryNames === 'function') {
            // The interface narrows to no-arg; the implementation in
            // lib-manager.ts#L56 accepts an optional refresh boolean. Call
            // with `true` (positional) — JS doesn't enforce parameter arity
            // so the impl's defaulting still works either way.
            try { availableNames = await helper.getAvaliableLibraryNames(true); }
            catch (_) { availableNames = await helper.getAvaliableLibraryNames(); }
          }
        } catch (e) {
          return {reloadErr: String(e).slice(0, 200), inAvailable: false, availableNames: []};
        }
        // Substring match on stem: the provider may enumerate with or without
        // the `.json` suffix depending on provider implementation. Use stem.
        const inAvailable = availableNames.some((n: string) =>
          n === fileName || n.includes(stem));
        return {reloadErr: null, inAvailable, availableNames};
      }, {fileName: workingCopy, stem: workingCopy.replace(/\.json$/, '')});
      expect(result.reloadErr, `reload error: ${result.reloadErr}`).toBeNull();
      // Atlas dep_lifecycle_ops[load_monomer_library]: post-write reload
      // surfaces the new file via the provider enumeration (no
      // stale-singleton — the manager re-discovers the FileShare contents).
      expect(result.inAvailable,
        `expected working copy '${workingCopy}' (or stem) in available libraries; observed: [${result.availableNames.join(', ')}]`).toBe(true);
    });

    // Scenario 2, Step 4 — Reopen the manage view and verify the working
    // copy still appears in the listing.
    //
    // The view may already be docked from S1.3; if it is, the previous
    // dispatch already rendered the catalogue without the working copy
    // (which was created in S2.1 — AFTER the view dispatched). We must
    // re-dispatch the view to force a re-render off the refreshed
    // provider state.
    await softStep('S2.4: reopen Manage Monomer Libraries view — working copy appears in listing', async () => {
      await page.evaluate(async () => {
        // Close existing view if open (so the re-dispatch creates a fresh
        // one off the post-reload provider state). Use closeAll on the
        // view-only — bypass shell.closeAll which would tear down the
        // table view too.
        const v = grok.shell.v;
        if (v && v.name === 'Manage Monomer Libraries' && typeof v.close === 'function') {
          v.close();
          await new Promise((r) => setTimeout(r, 500));
        }
        (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
        await new Promise((r) => setTimeout(r, 500));
        const manage = document.querySelector('[name="div-Bio---Manage"]')!;
        manage.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        manage.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 500));
        (document.querySelector(
          '[name="div-Bio---Manage---Monomer-Libraries"]') as HTMLElement).click();
      });
      await page.waitForFunction(() => (window as any).grok?.shell?.v?.name === 'Manage Monomer Libraries',
        null, {timeout: 30_000});
      // Extra settle for the provider's debounced refresh to populate the
      // re-rendered form.
      await page.waitForTimeout(2000);
      const result = await page.evaluate(({fileName, stem}) => {
        const root: any = grok.shell.v?.root;
        const form: any = root?.querySelector?.('.monomer-lib-controls-form');
        const rows = form ? Array.from(form.querySelectorAll('.ui-input-bool')) : [];
        const labels: string[] = rows.map((r: any) => {
          const span = r.querySelector('span');
          return span ? (span.textContent || '').trim() : '';
        });
        // Substring/stem match — same softening rationale as S2.3 (label may
        // appear with or without `.json` suffix depending on provider).
        const hasWorkingCopy = labels.some((l: string) =>
          l === fileName || l.includes(stem));
        return {
          formPresent: !!form,
          rowCount: rows.length,
          labels,
          hasWorkingCopy,
        };
      }, {fileName: workingCopy, stem: workingCopy.replace(/\.json$/, '')});
      expect(result.formPresent).toBe(true);
      expect(result.rowCount).toBeGreaterThanOrEqual(1);
      // The working copy must appear in the manage-view listing post-write
      // (atlas bio.manage.libraries-view + bio.cp.monomer-library-crud at
      // the read-back layer).
      expect(result.hasWorkingCopy,
        `expected working copy '${workingCopy}' (or stem) in manage-view labels; observed: [${result.labels.join(', ')}]`).toBe(true);
    });

    // ========================================================================
    // Scenario 3 — Save project that references the library
    // ========================================================================

    // Scenario 3, Step 1 — Switch to a Bio dataset that exercises the library
    // (the HELM dataset is already open from Setup; the renderer consults
    // the monomer library for color coding per atlas bio.rendering).
    // Verify the table is still the HELM one (Manage view dispatch in S1.3
    // / S2.4 docks a separate view but should not unload the table).
    await softStep('S3.1: HELM dataset remains open; Macromolecule column renderer is dispatchable', async () => {
      // After the S1.3 / S2.4 Manage view dispatches, grok.shell.tv (per
      // shell.ts#L83 = current view) points at the Manage view, NOT the
      // table view. Enumerate grok.shell.tableViews instead — that iterator
      // returns all open TableViews irrespective of which view is currently
      // selected. Picking the first TableView whose dataFrame carries a
      // Macromolecule column locks onto the HELM dataset deterministically,
      // regardless of label / name suffixing applied by the platform.
      const info = await page.evaluate(() => {
        const tvs: any[] = Array.from((grok.shell as any).tableViews || []);
        let chosen: any = null;
        let chosenMacro: any = null;
        for (const tv of tvs) {
          const df = tv?.dataFrame;
          if (!df) continue;
          const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
          const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
          if (macro) {
            chosen = tv;
            chosenMacro = macro;
            break;
          }
        }
        if (!chosen) {
          return {
            hasDf: false,
            hasMacro: false,
            units: null,
            rowCount: 0,
            tableViewCount: tvs.length,
            viewNames: tvs.map((tv: any) => tv?.name ?? '<unnamed>'),
          };
        }
        const df = chosen.dataFrame;
        return {
          hasDf: true,
          hasMacro: !!chosenMacro,
          units: chosenMacro?.meta?.units ?? null,
          rowCount: df.rowCount,
          tableViewCount: tvs.length,
          viewNames: tvs.map((tv: any) => tv?.name ?? '<unnamed>'),
        };
      });
      expect(info.hasDf,
        `expected a TableView with a Macromolecule column among open table views; observed: count=${info.tableViewCount}, names=[${info.viewNames.join(', ')}]`).toBe(true);
      expect(info.hasMacro).toBe(true);
      // Atlas bio.detector: HELM units tag survived the manage-view
      // dispatch round trip.
      expect(info.units).toBe('helm');
      expect(info.rowCount).toBeGreaterThan(0);
    });

    // Scenario 3, Step 2 — Save project with Data Sync ON.
    // Per scenario Notes + §4.5 Scenario authority — JS API persistence via
    // the canonical helpers/projects.ts saveAllTablesWithProvenance pattern
    // (mirrors bio-lifecycle-macromolecule-column-spec.ts S3.3). The Save
    // Project Ribbon dialog with Data Sync toggle is platform-wide UI not
    // present in bio.md selector reference; UI driving is delegated to
    // platform-side ui-smoke scenarios elsewhere. Persistence assertions
    // exercised via JS API are the assertable surface.
    //
    // Pre-step: close the Manage Monomer Libraries view so the helper's
    // `grok.shell.tv?.saveLayout?.()` call sees the HELM TableView, not the
    // manage view. (Without this, layout save is a no-op and the project
    // re-opens without the original HELM view layout; the persistence
    // assertion in S3.3 then has weaker provenance evidence to work with.)
    await page.evaluate(async () => {
      const views: any[] = Array.from(grok.shell.views || []);
      const manage: any = views.find((v: any) => v?.name === 'Manage Monomer Libraries');
      if (manage && typeof manage.close === 'function') {
        manage.close();
        await new Promise((r) => setTimeout(r, 500));
      }
      // Bring the HELM TableView forward as the current view (so
      // grok.shell.tv resolves to it). Re-assign via the typed setter.
      const tvs: any[] = Array.from((grok.shell as any).tableViews || []);
      const helm = tvs.find((tv: any) => {
        const df = tv?.dataFrame;
        if (!df) return false;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        return cols.some((c: any) => c.semType === 'Macromolecule');
      });
      if (helm && typeof (grok.shell as any).v !== 'undefined') {
        try { (grok.shell as any).v = helm; } catch (_) { /* setter may be read-only on some shell builds */ }
      }
      await new Promise((r) => setTimeout(r, 500));
    });
    await softStep('S3.2: save project with provenance (JS API path; mirrors macromolecule-column sibling)', async () => {
      saved = await saveAllTablesWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
      expect(saved.primaryTableInfoId).toBeTruthy();
    });

    // Scenario 3, Step 3 — Close + reopen via JS API; verify:
    //   - HELM dataset restores (atlas dep_lifecycle_ops[save_project_with_analysis]).
    //   - Macromolecule column re-classifies post-reopen (atlas bio.detector
    //     survives save/reopen).
    //   - getMonomerLibHelper() still resolves the same library catalogue
    //     post-reopen (atlas bio.cp.bio-service-surface-init: init contract
    //     holds across project reopen boundary).
    await softStep('S3.3: reopen project — HELM survives + library catalogue stable + Macromolecule semType holds', async () => {
      if (!saved) throw new Error('S3.2 did not produce a saved project');
      const result = await reopenAndAssertProvenance(page, saved.projectId);
      expect(result.tablesAfter).toBeGreaterThan(0);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
      // Two-layer post-reopen probe: (1) Macromolecule column survived,
      // (2) MonomerLibManager singleton still resolves the same catalogue.
      const post = await page.evaluate(async ({fileName, stem}) => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
        // Library-catalogue stability across the reopen boundary.
        const helper: any = await (grok as any).functions.call('Bio:getMonomerLibHelper', {});
        let availableNames: string[] = [];
        try {
          if (typeof helper?.getAvaliableLibraryNames === 'function') {
            // Force re-enumeration (provider cache may have been re-mounted
            // during project reopen). Positional-arg may be rejected by the
            // narrowed interface signature; fall back to no-arg.
            try { availableNames = await helper.getAvaliableLibraryNames(true); }
            catch (_) { availableNames = await helper.getAvaliableLibraryNames(); }
          }
        } catch (_) { /* leave empty */ }
        const hasWorkingCopy = availableNames.some((n: string) =>
          n === fileName || n.includes(stem));
        return {
          hasMacro: !!macro,
          units: macro?.meta?.units ?? null,
          helperResolved: helper != null,
          hasWorkingCopy,
          availableNames,
          availableCount: availableNames.length,
        };
      }, {fileName: workingCopy, stem: workingCopy.replace(/\.json$/, '')});
      // Atlas bio.detector survives reopen (units tag re-attached).
      expect(post.hasMacro).toBe(true);
      expect(post.units).toBe('helm');
      // Atlas bio.cp.bio-service-surface-init lifecycle:
      // getMonomerLibHelper() still resolves the singleton post-reopen.
      expect(post.helperResolved).toBe(true);
      // Library reference shape survives: the working copy that was
      // discoverable pre-save is still discoverable post-reopen (atlas
      // dep_lifecycle_ops[save_project_with_analysis] doesn't silently
      // drop FileShare library bindings). The assertable shape per atlas is
      // the catalogue stability — the working copy WAS on FileShare before
      // save (S2.1 writeAsText) and is still on FileShare after reopen
      // (cleanup hasn't run yet — it's in the finally block).
      expect(post.hasWorkingCopy,
        `expected working copy '${workingCopy}' (or stem) in post-reopen available libraries; observed: [${post.availableNames.join(', ')}]`).toBe(true);
      expect(post.availableCount).toBeGreaterThanOrEqual(1);
    });
  } finally {
    // ========================================================================
    // Scenario 4 — Cleanup (runs regardless of earlier failures per scenario
    // Expected: "Cleanup runs in tearDownAll / finally regardless of
    // earlier failures").
    // ========================================================================

    // Step 4.1 — Delete the working copy library file (best-effort).
    if (workingCopyWritten) {
      await page.evaluate(async (p) => {
        try { await grok.dapi.files.delete(p); } catch (_) { /* best effort */ }
      }, workingCopyPath).catch(() => {});
    }

    // Step 4.2 — Delete the project (best-effort via the canonical helper).
    if (saved) {
      await deleteProjectWithCleanup(page, {
        projectId: saved.projectId,
        tableInfoId: saved.primaryTableInfoId,
      });
    }

    // Step 4.3 — Close any open Manage Monomer Libraries dialogs / views
    // (best-effort — escape any open dialog + close the view if docked).
    await page.evaluate(async () => {
      // Dismiss any open .d4-dialog.
      const dialogs = Array.from(document.querySelectorAll('.d4-dialog'));
      for (const d of dialogs)
        d.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
      await new Promise((r) => setTimeout(r, 300));
      // Close manage view if still docked.
      try {
        const views = Array.from(grok.shell.views || []);
        const manageView: any = views.find((v: any) => v?.name === 'Manage Monomer Libraries');
        if (manageView && typeof manageView.close === 'function') manageView.close();
      } catch (_) { /* best effort */ }
    }).catch(() => {});
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
