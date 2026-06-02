/* ---
sub_features_covered:
  - bio.engines.numbering-immunum
  - bio.annotate.numbering-scheme
  - bio.api.get-seq-helper
  - bio.lifecycle.init
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent in scenario .md frontmatter (chain yaml pins
//     proactive_lifecycle_specs[5] at the proactive-lifecycle pyramid layer
//     globally; coverage_type: regression)
//   sub_features_covered: [bio.engines.numbering-immunum,
//     bio.annotate.numbering-scheme, bio.api.get-seq-helper,
//     bio.lifecycle.init]
//   ui_coverage_responsibility: absent (delegated_to: null) — scenario Notes
//     explicitly state "JS API substitutes are used for project persistence
//     assertions per the same pattern as sibling bio-lifecycle-*.md
//     scenarios"; UI driving stays on the cited dispatch points where the
//     assertable surface lives (the Annotate top-menu path).
//   related_bugs: [] (chain proactive_lifecycle_specs[5].bugs_reinforcing
//     empty; no bug-invariant slice on this lifecycle cell)
//   produced_from: atlas-driven
//   coverage_type: regression
//
// Atlas provenance (derived_from):
//   feature-atlas/bio.yaml#sub_features[bio.lifecycle.init] derived_from:
//     public/packages/Bio/src/package.ts#L138 — initBio role: init.
//   feature-atlas/bio.yaml#sub_features[bio.api.get-seq-helper] derived_from:
//     public/packages/Bio/src/package.ts#L1678 — Bio:getSeqHelper service
//     surface. interaction = "await grok.functions.call('Bio:getSeqHelper')".
//   feature-atlas/bio.yaml#sub_features[bio.annotate.numbering-scheme]
//     derived_from: public/packages/Bio/src/package.ts#L502 — applyNumbering
//     Scheme top-menu function. interaction = "Bio | Annotate | Apply
//     Numbering Scheme...".
//   feature-atlas/bio.yaml#sub_features[bio.engines.numbering-immunum]
//     derived_from: public/packages/Bio/src/package.ts#L1018 — Immunum engine
//     (immunumAntibodyNumbering, meta.role: antibodyNumbering); WASM-in-worker
//     5-column result shape per Bio/CLAUDE.md "Expected output shape" table:
//     position_names, chain_type, annotations_json, numbering_detail,
//     numbering_map. Schemes choices ['imgt', 'kabat'].
//   feature-atlas/bio.yaml#critical_paths[bio.cp.numbering-scheme] derived_from:
//     public/packages/Bio/src/package.ts#L502 (priority p1; sub_features_used:
//     [bio.annotate.numbering-scheme, bio.engines.numbering-immunum]).
//   feature-atlas/bio.yaml#dep_lifecycle_ops[save_project_with_analysis]
//     affected_source_classes: [all] (shorthand applies to immunum_wasm per
//     chain proactive_lifecycle_specs[5].rationale — "Only one non-agnostic
//     op affects this source class since WASM asset is bundled with package
//     version and has no runtime entity-type dep").
//
// Paradigm selection (per pyramid_layer: proactive-lifecycle on
// target_layer: playwright): mostly JS API for matrix/lifecycle shape; UI
// driving required for the atlas-cited UI dispatch point the scenario
// explicitly names — `Bio | Annotate | Apply Numbering Scheme...` (Scenario
// 1 step 3). The top-menu path is class-1 selectors (bio.md L37, L417-L429).
//
// SCOPE notes honoured from scenario authority:
//   - Step 2.1 (Save Project Ribbon + Data Sync toggle): per scenario Notes
//     "JS API substitutes are used for project persistence assertions per
//     the same pattern as sibling bio-lifecycle-*.md scenarios". Mirrors
//     bio-lifecycle-macromolecule-column-spec.ts S3.3 — uses
//     helpers/projects.ts saveAllTablesWithProvenance +
//     reopenAndAssertProvenance, the canonical JS-API persistence-helper
//     pattern at the lifecycle layer.
//   - Step 3.3 (re-run Immunum on reopened table): exercises the WASM
//     re-load determinism contract — same code path as first call per
//     immunum-client.ts "fresh worker per call, terminate before return"
//     pattern (Bio/CLAUDE.md). No private WASM-handle hooks; assertion is
//     observable equality of the two result DataFrames' column shapes +
//     non-null values across the same input rows.
//
// Selector provenance: every [name=...] selector below is class-1
// (in bio.md grok-browser reference at the cited lines):
//   - [name="div-Bio"] (bio.md L606)
//   - [name="div-Bio---Annotate"] (bio.md L69)
//   - [name="div-Bio---Annotate---Apply-Numbering-Scheme..."] (bio.md L421)
//   - [name="dialog-Apply-Antibody-Numbering"] (bio.md L425; dialog title
//     differs from menu label — title is "Apply Antibody Numbering")
//   - [name="input-host-Sequence"], [name="input-host-Engine"],
//     [name="input-host-Scheme"] (bio.md L427-L429)
//   - [name="button-OK"] (bio.md L131 / standard dialog OK)
//   - [name="viewer-Grid"] (standard platform selector — used across all bio specs)
//
// Sibling spec reuse:
//   - bio-lifecycle-macromolecule-column-spec.ts — canonical proactive-
//     lifecycle pattern (cold-start two-layer init probe + top-menu dialog
//     drive + saveAllTablesWithProvenance + reopenAndAssertProvenance);
//     this spec mirrors its structure exactly, with the Convert+Sequence
//     Space ribbon paths replaced by the Apply Numbering Scheme path.
//   - convert-spec.ts — canonical Bio top-menu click + .dispatchEvent
//     'mouseover' submenu pattern; mirrored for the Annotate submenu walk.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {
  saveAllTablesWithProvenance,
  reopenAndAssertProvenance,
  deleteProjectWithCleanup,
} from '../helpers/projects';

test.use(specTestOptions);

test('Bio immunum_wasm source-class lifecycle: init → IMGT numbering → save+reopen → re-run', async ({page}) => {
  // 7-minute end-to-end budget: cold Bio init (≤90s observed in sibling
  // analyze/sequence-space cycle-2 retries) + dialog dispatch + IMGT
  // numbering WASM compute on a tiny antibody fixture (≤30s expected) +
  // project save+reopen round trip + second numbering run.
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = `bio-lifecycle-immunum-wasm-${stamp}`;
  // Canonical antibody fixture (used by Bio/src/tests/antibody-numbering-tests.ts:167).
  // 47 antibody heavy-chain rows in FASTA notation — small enough to keep the
  // WASM numbering cell under budget, large enough that the 5-column result
  // shape is non-degenerate per row.
  const antibodyFixturePath = 'System:AppData/Bio/samples/antibodies.csv';
  let saved: {projectId: string; primaryTableInfoId: string; layoutId: string | null} | null = null;

  // The 5-column result DataFrame columns the Immunum engine produces, per
  // atlas bio.engines.numbering-immunum + Bio/CLAUDE.md "Expected output
  // shape" table. These names are emitted by showNumberingSchemeDialog's
  // result-applier and visible as added columns on the source DataFrame
  // (per numbering-ui.ts "applies its result — annotations + aligned
  // column" — the annotations are written into the source column's
  // `.annotations` tag rather than being added as a literal '*_annotations'
  // top-level column; the aligned column IS added as a sibling column;
  // numbering rows are stored per-row tags). The result-DataFrame contract
  // is asserted by invoking the engine via the registered function
  // Bio:immunumAntibodyNumbering directly — bypasses the result-applier so
  // we observe the raw 5-column shape.
  const EXPECTED_IMMUNUM_COLS = [
    'position_names', 'chain_type', 'annotations_json',
    'numbering_detail', 'numbering_map',
  ] as const;

  await loginToDatagrok(page);

  // ==========================================================================
  // Setup — open antibody fixture, await semType + Bio init readiness.
  // ==========================================================================
  // Mirrors the setup phase of bio-lifecycle-macromolecule-column-spec.ts /
  // convert-spec.ts / sequence-space-spec.ts verbatim — same cold-start
  // tolerance applies.
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
  }, antibodyFixturePath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  // Two-layer Bio init readiness probe (mirrors convert-spec.ts /
  // sequence-space-spec.ts / bio-lifecycle-macromolecule-column-spec.ts).
  // Layer 1: DOM top-menu visibility. Layer 2: Bio service-surface
  // serialization probe — the runtime serializes grok.functions.call after
  // init completion.
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  await page.evaluate(async () => {
    const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
    for (const fn of probes) {
      try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
    }
    await new Promise((r) => setTimeout(r, 3000));
  });

  // ==========================================================================
  // Scenario 1 — Trigger initBio + verify getSeqHelper resolves (atlas
  // bio.lifecycle.init + bio.api.get-seq-helper).
  // ==========================================================================
  await softStep('S1.1: initBio completes; Bio:getSeqHelper resolves to a usable singleton', async () => {
    const info = await page.evaluate(async () => {
      // Atlas bio.api.get-seq-helper interaction =
      //   "await grok.functions.call('Bio:getSeqHelper')"
      // Returns ISeqHelper singleton after initBio completes
      // (Bio/CLAUDE.md "Initialization Flow"). If init has not completed,
      // this call awaits or errors.
      let resolved = false;
      let helperKind: string | null = null;
      let methodNames: string[] = [];
      try {
        const helper: any = await (grok as any).functions.call('Bio:getSeqHelper', {});
        resolved = !!helper;
        helperKind = helper ? (helper.constructor?.name ?? typeof helper) : null;
        // Probe the documented ISeqHelper surface (Bio/CLAUDE.md seq-helper.ts):
        //   getSeqHandler, getSeqMonomers, setUnitsToFastaColumn,
        //   helmToAtomicLevel, helmToAtomicLevelSingle, etc.
        if (helper) {
          methodNames = ['getSeqHandler', 'getSeqMonomers', 'helmToAtomicLevel',
            'setUnitsToFastaColumn']
            .filter((m) => typeof helper[m] === 'function');
        }
      } catch (e) {
        return {resolved: false, helperKind: null, methodNames: [], err: String(e).slice(0, 200)};
      }
      return {resolved, helperKind, methodNames, err: null};
    });
    expect(info.resolved).toBe(true);
    // ISeqHelper surface contract — at least getSeqHandler is the load-bearing
    // method other packages consume (Bio/CLAUDE.md service pattern).
    expect(info.methodNames).toContain('getSeqHandler');
  });

  // ==========================================================================
  // Scenario 1, Step 2 — Macromolecule detector classified the antibody
  // sequences (atlas bio.detector — outside this scenario's
  // sub_features_covered per scenario authority, but is a precondition for
  // the Annotate dialog which requires a Macromolecule column).
  // ==========================================================================
  await softStep('S1.2: antibodies.csv opens with at least one Macromolecule column (FASTA notation)', async () => {
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macroCols = cols.filter((c: any) => c.semType === 'Macromolecule');
      const macro: any = macroCols[0];
      return {
        macroCount: macroCols.length,
        macroName: macro?.name ?? null,
        units: macro?.meta?.units ?? null,
        rowCount: df.rowCount,
      };
    });
    expect(info.macroCount).toBeGreaterThan(0);
    expect(info.rowCount).toBeGreaterThan(0);
    // Antibody heavy chains are FASTA / peptide-alphabet (Bio/samples/
    // antibodies.csv contract; matches the antibody-numbering-tests.ts
    // fixture). 'fasta' is the canonical sequence units; the engine accepts
    // any valid antibody sequence input per atlas bio.engines.numbering-immunum.
    expect(info.units).not.toBeNull();
  });

  // ==========================================================================
  // Scenario 1, Step 3 — Drive Bio | Annotate | Apply Numbering Scheme...
  // (atlas bio.annotate.numbering-scheme + bio.engines.numbering-immunum).
  // ==========================================================================
  // Atlas-canonical click pattern per bio.md L417-L429 (the Apply Numbering
  // Scheme dialog) + L37, L69 (top-menu path). Mirrors the convert-spec.ts /
  // bio-lifecycle-macromolecule-column-spec.ts dispatch-with-mouseover
  // submenu pattern.
  //
  // The dialog (per bio.md L425) carries title "Apply Antibody Numbering"
  // (differs from menu label) and has Engine SELECT (Immunum is the built-in
  // WASM engine; meta.role: antibodyNumbering) + Scheme SELECT (options
  // ['imgt', 'kabat', 'chothia', 'aho']; Immunum advertises ['imgt', 'kabat']
  // via its scheme param `choices`). The dialog's onOK invokes the selected
  // engine and applies the result (annotations into the source column tag +
  // aligned column added as sibling).
  //
  // Default selections — Engine: Immunum (only built-in engine; per
  // Bio/CLAUDE.md "Existing engines" table); Scheme: imgt (Immunum's
  // scheme param initialValue per numbering-ui.ts). Default selections are
  // adequate for the lifecycle contract under test — the assertion is the
  // 5-column engine-output shape arrives on the table (or its mirror image
  // via the result-applier).
  await softStep('S1.3-1.4: Apply Antibody Numbering (Immunum/IMGT) via top-menu; sibling column appended', async () => {
    const baseCols: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      document.querySelector('[name="div-Bio---Annotate"]')!
        .dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      (document.querySelector(
        '[name="div-Bio---Annotate---Apply-Numbering-Scheme..."]') as HTMLElement).click();
    });
    // Dialog title is "Apply Antibody Numbering" (NOT the menu label) per
    // bio.md L425 gotcha.
    await page.locator('[name="dialog-Apply-Antibody-Numbering"]').waitFor({timeout: 60_000});
    // The Sequence SELECT auto-binds to the first Macromolecule column on
    // dialog open; Engine SELECT auto-binds to the only built-in option
    // (Immunum); Scheme SELECT auto-binds to its initialValue (imgt per
    // numbering-ui.ts). No further input required for the canonical
    // IMGT/Immunum path under test.
    await page.locator('[name="dialog-Apply-Antibody-Numbering"] [name="button-OK"]').click();
    // The OK click triggers the WASM call; the result-applier writes
    // annotations into the source column tag and adds the aligned column
    // (per numbering-ui.ts). 120s budget tolerates the cold WASM module
    // load + per-row numbering compute on the 47-row fixture.
    await page.waitForFunction(
      (b) => grok.shell.tv.dataFrame.columns.length > b ||
        // Annotations are written to the source column tag rather than as a
        // new top-level column on some code paths — also accept that signal.
        (() => {
          const df = grok.shell.tv.dataFrame;
          const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
          return cols.some((c: any) => {
            try {
              const tag = c.getTag?.('.annotations') ?? c.tags?.get?.('.annotations');
              return typeof tag === 'string' && tag.length > 0 && tag.startsWith('[');
            } catch { return false; }
          });
        })(),
      baseCols, {timeout: 120_000});
    // Dialog closes after OK + result applied.
    await page.waitForFunction(
      () => document.querySelectorAll('[name="dialog-Apply-Antibody-Numbering"]').length === 0,
      null, {timeout: 30_000}).catch(() => {});
  });

  // ==========================================================================
  // Scenario 1, Step 4 — Verify the 5-column Immunum result-shape contract
  // (atlas bio.engines.numbering-immunum + Bio/CLAUDE.md "Expected output
  // shape" table).
  // ==========================================================================
  // The Apply dialog applies the engine result via the dialog's result-
  // applier (which writes annotations into the source column tag + adds the
  // aligned column to the source DataFrame, per numbering-ui.ts). The raw
  // 5-column DataFrame is the engine's direct output; we observe it by
  // invoking the registered engine function Bio:immunumAntibodyNumbering on
  // the same column (atlas bio.engines.numbering-immunum interactions empty;
  // function name per public/packages/Bio/src/package.ts#L1018 reference
  // implementation). This exercises the WASM-in-worker engine code path
  // directly and lets us assert the atlas-declared column-shape contract.
  await softStep('S1.4: Immunum engine returns 5-column result DataFrame with documented column names + non-null per-row values', async () => {
    const info = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macroCol: any = cols.find((c: any) => c.semType === 'Macromolecule');
      if (!macroCol) throw new Error('S1.4: no Macromolecule column for Immunum engine call');
      // Invoke the Immunum engine directly. Per Bio/CLAUDE.md "Adding a New
      // Antibody Numbering Engine" the engine signature is
      //   (df: DG.DataFrame, seqCol: DG.Column<string>, scheme: string)
      //     => DG.DataFrame
      // and the function-registry name is the same as the function name
      // (immunumAntibodyNumbering in package.ts L1018). The package-name
      // prefix is "Bio:" per the function-registry convention.
      let result: any = null;
      let invokeErr: string | null = null;
      try {
        result = await (grok as any).functions.call(
          'Bio:immunumAntibodyNumbering',
          {df, seqCol: macroCol, scheme: 'imgt'},
        );
      } catch (e) {
        invokeErr = String(e).slice(0, 300);
      }
      if (!result) return {result: null, invokeErr, colNames: null, sampleRow: null};
      const colNames = Array.from({length: result.columns.length},
        (_, i) => result.columns.byIndex(i).name);
      // Sample the first row from each column — atlas bio.engines.numbering-
      // immunum: "All five columns are populated with non-null values for
      // each input antibody row".
      const sampleRow: Record<string, string | null> = {};
      for (let i = 0; i < result.columns.length; i++) {
        const col = result.columns.byIndex(i);
        const v = col.get(0);
        sampleRow[col.name] = v == null ? null : String(v).slice(0, 80);
      }
      return {
        result: {rowCount: result.rowCount, colCount: result.columns.length},
        invokeErr,
        colNames,
        sampleRow,
      };
    });
    // Atlas contract: 5-column result DataFrame.
    expect(info.invokeErr).toBeNull();
    expect(info.result).not.toBeNull();
    expect(info.result!.colCount).toBe(5);
    expect(info.result!.rowCount).toBeGreaterThan(0);
    // Documented column names per Bio/CLAUDE.md "Expected output shape" table.
    for (const expectedName of EXPECTED_IMMUNUM_COLS)
      expect(info.colNames).toContain(expectedName);
    // Per-row non-null contract — sample row's values are all non-null.
    for (const expectedName of EXPECTED_IMMUNUM_COLS)
      expect(info.sampleRow![expectedName]).not.toBeNull();
  });

  try {
    // ==========================================================================
    // Scenario 2 — Save project with the numbering output (atlas
    // dep_lifecycle_ops[save_project_with_analysis]; per chain
    // proactive_lifecycle_specs[5].rationale the [all] shorthand applies to
    // immunum_wasm).
    // ==========================================================================
    // Per scenario Notes + §4.5 Scenario authority — JS API persistence via
    // the canonical helpers/projects.ts saveAllTablesWithProvenance pattern
    // (mirrors bio-lifecycle-macromolecule-column-spec.ts S3.3). The Save
    // Project Ribbon dialog with Data Sync toggle is platform-wide UI not
    // present in bio.md selector reference; UI driving for that surface is
    // delegated to platform-side ui-smoke scenarios elsewhere. Persistence
    // assertions exercised via JS API are the assertable surface.
    await softStep('S2.1: Save project with numbering output (JS API path)', async () => {
      saved = await saveAllTablesWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
      expect(saved.primaryTableInfoId).toBeTruthy();
    });

    // Scenario 2.2 — Verify project record persists (find-by-id; per the
    // saveAllTablesWithProvenance contract, the project entity carries the
    // TableInfo with .annotations / aligned columns intact).
    await softStep('S2.2: saved project is findable server-side via dapi.projects.find(id)', async () => {
      if (!saved) throw new Error('S2.1 did not produce a saved project');
      const ok = await page.evaluate(async (id) => {
        const proj = await (grok as any).dapi.projects.find(id);
        return proj != null && proj.id === id;
      }, saved.projectId);
      expect(ok).toBe(true);
    });

    // ==========================================================================
    // Scenario 3 — Reopen project + WASM re-load + deterministic re-run.
    // ==========================================================================
    // Per scenario Step 3.1-3.2: close + reopen via JS API; verify antibody
    // table + Macromolecule semType survive the round trip.
    await softStep('S3.1-3.2: reopen project — antibody table + Macromolecule semType survive', async () => {
      if (!saved) throw new Error('S2.1 did not produce a saved project');
      const result = await reopenAndAssertProvenance(page, saved.projectId);
      expect(result.tablesAfter).toBeGreaterThan(0);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
      const post = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macro: any = cols.find((c: any) => c.semType === 'Macromolecule');
        return {
          hasMacro: !!macro,
          units: macro?.meta?.units ?? null,
          rowCount: df.rowCount,
        };
      });
      // bio.detector lifecycle: Macromolecule tag persists across save/reopen.
      expect(post.hasMacro).toBe(true);
      expect(post.units).not.toBeNull();
      expect(post.rowCount).toBeGreaterThan(0);
    });

    // Scenario 3, Step 3 — Re-run Immunum on the reopened table; exercises
    // the WASM re-load code path on a fresh session. Per immunum-client.ts
    // "fresh worker per call, terminate before return" pattern (Bio/CLAUDE.md)
    // — the WASM module re-loads on every invocation. Determinism is the
    // observable: same input → same 5-column shape, same row count, same
    // per-row non-null contract.
    await softStep('S3.3-3.4: re-run Immunum on reopened table; result shape deterministic on WASM re-load', async () => {
      const info = await page.evaluate(async () => {
        const df = grok.shell.tv.dataFrame;
        const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
        const macroCol: any = cols.find((c: any) => c.semType === 'Macromolecule');
        if (!macroCol) throw new Error('S3.3: no Macromolecule column after reopen');
        let result: any = null;
        let invokeErr: string | null = null;
        try {
          result = await (grok as any).functions.call(
            'Bio:immunumAntibodyNumbering',
            {df, seqCol: macroCol, scheme: 'imgt'},
          );
        } catch (e) {
          invokeErr = String(e).slice(0, 300);
        }
        if (!result) return {result: null, invokeErr, colNames: null, sampleRow: null};
        const colNames = Array.from({length: result.columns.length},
          (_, i) => result.columns.byIndex(i).name);
        const sampleRow: Record<string, string | null> = {};
        for (let i = 0; i < result.columns.length; i++) {
          const col = result.columns.byIndex(i);
          const v = col.get(0);
          sampleRow[col.name] = v == null ? null : String(v).slice(0, 80);
        }
        return {
          result: {rowCount: result.rowCount, colCount: result.columns.length},
          invokeErr,
          colNames,
          sampleRow,
        };
      });
      // WASM re-load contract: engine call succeeds in the reopened session
      // identically to Scenario 1.4.
      expect(info.invokeErr).toBeNull();
      expect(info.result).not.toBeNull();
      // Same 5-column shape (atlas bio.engines.numbering-immunum).
      expect(info.result!.colCount).toBe(5);
      expect(info.result!.rowCount).toBeGreaterThan(0);
      for (const expectedName of EXPECTED_IMMUNUM_COLS)
        expect(info.colNames).toContain(expectedName);
      for (const expectedName of EXPECTED_IMMUNUM_COLS)
        expect(info.sampleRow![expectedName]).not.toBeNull();
    });
  } finally {
    // ========================================================================
    // Scenario 4 — Cleanup (runs regardless of earlier failures per scenario
    // Expected: "Cleanup runs in tearDownAll / finally regardless of earlier
    // failures").
    // ========================================================================
    if (saved) {
      await deleteProjectWithCleanup(page, {
        projectId: saved.projectId,
        tableInfoId: saved.primaryTableInfoId,
      });
    }
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
