/* ---
sub_features_covered: [powerpack.dialogs.add-new-column, powerpack.dialogs.add-new-column-func, powerpack.dialogs.prepare-add-column-call, powerpack.formula.is-formula-column, powerpack.formula.widget, powerpack.dialogs]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: bug-focused
//   sub_features_covered: [powerpack.dialogs.add-new-column,
//     powerpack.dialogs.add-new-column-func,
//     powerpack.dialogs.prepare-add-column-call,
//     powerpack.formula.is-formula-column,
//     powerpack.formula.widget,
//     powerpack.dialogs]
//   ui_coverage_responsibility: [add-new-column-dialog,
//     context-panel-formula-edit,
//     add-new-column-formula-recalc-dependency,
//     save-project-with-formula-columns-persistence]
//     (delegated_to: add-new-column.md — the BASIC add-column-dialog
//      open + Recent-Activities autofill is owned by the chain's smoke
//      witness; this spec specializes on the chained-dependency
//      build, the context-panel Formula widget edit + recalc, and the
//      GROK-17109 save+reopen persistence invariant)
//   related_bugs: [GROK-17109]
//   produced_from: migrated (original at
//     PowerPack/AddNewColumn/formula-refreshing.md — split out by
//     Migrator 2026-05-20 with the Additional-Notes persistence check
//     promoted to an explicit scenario block).
//
// Bug-library cross-reference:
//   GROK-17109 (Calculated columns: columns are not saved to project
//     with data sync). The third scenario block (Save the project,
//     reopen, verify chained columns persist) is the canonical
//     GROK-17109 regression surface — exercised here for the chained
//     calc-column case (Weight2 → Weight3 → Weight4). Atlas affects[]:
//     add-new-column, add-new-column-func, prepare-add-column-call,
//     is-formula-column, formula.widget, dialogs — covered verbatim by
//     frontmatter sub_features_covered.
//   Chain witness role: per scenario Notes "Chain context" — one of
//     three GROK-17109 spans (alongside the top-level add-new-column.md
//     smoke and the multi-source AddNewColumn/add-new-column.md
//     datasync flow); this scenario contributes the dependency-chain
//     recalc invariant (Weight2 → Weight3 → Weight4) PLUS the
//     save+reopen persistence check for chained calc columns.
//
// Atlas provenance (derived_from): per feature-atlas/powerpack.yaml
//   the six consumed sub_features are at:
//     powerpack.formula                — L287
//     powerpack.formula.is-formula-column — L294
//     powerpack.formula.widget         — L301
//     powerpack.dialogs                — L587
//     powerpack.dialogs.add-new-column — L615
//       derived_from: public/packages/PowerPack/src/dialogs/add-new-column.ts#L98
//     powerpack.dialogs.add-new-column-func — L622
//     powerpack.dialogs.prepare-add-column-call — L629
//       derived_from: public/packages/PowerPack/src/dialogs/add-new-column.ts#L1638
//   Critical path consulted:
//     powerpack.cp.add-new-column-persists — L1053
//       derived_from: public/packages/PowerPack/src/dialogs/add-new-column.ts#L98
//       sub_features_used matches frontmatter exactly (5 of 6 — minus
//       the parent powerpack.dialogs container).
//
// House-style anchors:
//   - public/packages/PowerPack/src/tests/add-new-column.ts
//     (apitests-layer sibling — AddNewColumnDialog instantiated directly
//      + `dlg.codeMirror!.dispatch({changes: {...}})` pattern for
//      formula composition; canonical reference for the formula-tag
//      shape per add-new-column.ts:96-129).
//   - public/packages/UsageAnalysis/files/TestTrack/PowerPack/add-new-column-spec.ts
//     (ui-smoke sibling — basic dialog open + name input + CodeMirror
//      dispatch + [name="button-Add-New-Column---OK"] selectors).
//   - public/packages/UsageAnalysis/files/TestTrack/PowerPack/add-new-column-advanced-spec.ts
//     (bug-focused sibling — saveProjectWithProvenance + reopen +
//      formula-tag-preserved invariant for Weight2/Weight3 chain;
//      THIS spec extends the chain to Weight4 + adds the Formula info
//      panel UI surface).
//
// Reference templates:
//   - bug-focused: Projects/complex-derived-tables-spec.ts (GROK-19103
//     slice pattern — open → mutate → save-and-reopen via JS API +
//     find-by-id verification; deleteProjectWithCleanup in finally).
//   - sibling pattern: PowerPack/add-new-column-advanced-spec.ts
//     (most-direct template: same GROK-17109 invariant, same
//     openTableFromFile + saveProjectWithProvenance + reopen +
//     deleteProjectWithCleanup helper trio; THIS spec adds the
//     Formula info panel widget UI surface as the formula-edit path).
//
// Selector / API citations (all in current grok-browser/references or
// existing PowerPack/JS-API sources; no reference-file proposal needed):
//   - Toolbar icon: [name="icon-add-new-column"] — used unmodified
//     across both sibling specs; documented in
//     grok-browser/references/dialogs-menus.md "Toolbar Ribbon Icons".
//   - Dialog selectors: prepareForSeleniumTests in
//     public/packages/PowerPack/src/dialogs/add-new-column.ts:344-349
//     (name="input-Add-New-Column---Name",
//      name="button-Add-New-Column---OK",
//      name="button-Add-New-Column---CANCEL"); set unconditionally from
//     init() at L250 so reliable for both initial-dialog and
//     Formula-widget-bound-dialog instances (the widget at L194 of
//     public/packages/PowerPack/src/package.ts constructs the same
//     AddNewColumnDialog pre-bound to the column).
//   - CodeMirror dispatch: CM6 EditorView is exposed via
//     `cmContent.cmTile.view` (NOT `cmDiv.cmView.view`) — MCP recon
//     2026-05-26 cycle 03 retry round 1 confirmed: on the live page
//     `document.querySelector('.add-new-column-dialog-cm-div .cm-content').cmTile.view`
//     resolves to the EditorView with .state and .dispatch; the
//     `cmView` property does NOT exist on any ancestor of the
//     contenteditable host (chain probe across 5 ancestor levels:
//     all `cmView=undefined`). Earlier spec versions read
//     `(cmDiv as any).cmView?.view ?? (cmDiv.parentElement as any)?.cmView?.view`
//     — both always undefined → view-dispatch retry loop always
//     exhausted → fall through to keyboard fallback, which is brittle
//     on contenteditable. Root cause of cycle -02 / cycle -03 retry
//     [B-RUN-PASS, B-STAB-01] is THIS incorrect property path.
//     Empirical fix per MCP recon: read `cmContent.cmTile.view`. This
//     same property surface works for BOTH the toolbar-icon dialog
//     and the Formula info panel widget — verified live by composing
//     `${WEIGHT} + 100` into the toolbar dialog (Weight2 = 173.20 vs
//     WEIGHT=73.20, +100 exact) and `${WEIGHT} + 200` into the panel
//     widget (Weight2 = 273.20 with tag carrying `${WEIGHT} + 200`).
//     Same-paradigm tactical fix — no paradigm pivot.
//   - WIDGET-mode CM div class: empirical recon 2026-05-26 cycle 03
//     also revealed that the Formula info panel widget renders
//     `.add-new-column-dialog-cm-div` (NOT `.add-new-column-widget-cm-div`
//     as PowerPack source line 175 suggests). Reading the PowerPack
//     source closely confirms why: the constructor at
//     public/packages/PowerPack/src/dialogs/add-new-column.ts:174-190
//     evaluates `this.codeMirrorDiv.classList.add(this.widget ? ... :
//     ...)` BEFORE the subsequent `if (widget) this.widget = widget;`
//     assignment — so at the moment line 175 runs, `this.widget` is
//     `undefined` regardless of whether a widget was passed in, and the
//     `-dialog-cm-div` class wins. This is a quiescent PowerPack source
//     quirk (NOT a Datagrok bug we'd file; classList side effect
//     doesn't observably break the user flow); the spec accommodates
//     it by scoping `.add-new-column-dialog-cm-div .cm-content` to
//     `.grok-prop-panel` for the Formula info panel path. Widget mode
//     still does NOT call prepareForSeleniumTests() (guarded at
//     add-new-column.ts:217-250 by `if (!this.widget)`), so
//     `[name="button-Add-New-Column---OK"]` is unset on the widget;
//     the widget renders `ui.button('Apply', addNewColumnAction)` at
//     L262-266 and `editFormulaViaInfoPanel` clicks the Apply button
//     by text. Both adjustments are MCP-empirically-backed
//     same-paradigm tactical fixes; no paradigm pivot.
//   - Formula info panel accordion-pane toggle idempotency: MCP recon
//     2026-05-26 cycle 03 retry round 2 confirmed the Formula pane
//     header click TOGGLES expansion state (not idempotent). With the
//     pane already expanded (Context Panel auto-expands single-matching
//     widget panels; state persists across grok.shell.o column switches),
//     a blind click COLLAPSES the pane, leaving the .cm-content element
//     in the DOM but display:none (offsetWidth/Height === 0). The earlier
//     spec attempted `.cm-content { state: 'visible' }` Playwright wait
//     which then timed out at 15_000ms (the test trace shows
//     `34 × locator resolved to hidden`). Root cause of cycle -03 retry
//     round 1 [B2.3 timeout] + cascading [B3.3+4 expected "+ 50" got
//     "${Weight2} + 100"] — Weight3 edit never reached the dispatch step
//     because the pane was collapsed by step-2's blind click. Fix:
//     read `classList.contains('expanded')` / `aria-expanded` BEFORE
//     clicking; only click when NOT expanded. Empirical verification
//     via MCP: collapse-then-expand cycle reproduces "hidden in DOM"
//     state and recovers it via second click. Same-paradigm tactical
//     fix (UI driving path unchanged; idempotency check added).
//   - Apply button name= attribute: MCP recon 2026-05-26 cycle 03 retry
//     round 2 also revealed that the in-panel widget's Apply button
//     DOES carry `[name="button-Apply"]` (xamgle annotate() auto-adds
//     name= from button text — same convention as the dialog's OK
//     button). Spec now prefers `[name="button-Apply"]` over text-content
//     lookup, with text-content lookup retained as fallback.
//   - Formula info panel widget path: per
//     public/packages/PowerPack/src/package.ts:184-209, the panel is
//     registered as `@panel({name: 'Formula',
//     condition: 'PowerPack:isFormulaColumn(col)'})`. Selecting a
//     calc-column-flagged column populates the Context Panel with the
//     "Formula" accordion pane; opening that pane spawns the same
//     AddNewColumnDialog pre-bound to the column's existing formula
//     (the dialog's `addColumn: false` aux flag at L205 means OK
//     re-evaluates the existing column in-place rather than creating
//     a new one). The widget is invoked from the Context Panel via the
//     standard d4-accordion expansion of the "Formula" pane —
//     selectable by accordion-pane title text since `name=` attributes
//     on accordion panes are not stable. The direct JS API path
//     `PowerPack:formulaWidget(col)` is available as a deterministic
//     fallback when the accordion-pane click is brittle under headless
//     conditions (this still exercises the SAME widget code path —
//     the panel registration's render function — and therefore the
//     SAME powerpack.formula.widget sub_feature surface).
//
// Scope notes:
//   * Scenario steps are walked end-to-end EXCEPT for two scope
//     compressions in the Formula-info-panel block:
//     - "Modify the Weight4 formula" sub-step (block 2 step 4) is
//       executed but the assertion that Weight2/Weight3 are unaffected
//       is checked by reading their values pre-edit and post-edit
//       (we keep the pre-Weight4-edit values and confirm they did not
//       change). The scenario's "no errors are surfaced" check is
//       satisfied by absence of grok.shell.warnings additions.
//     - Formula editing for each of Weight2/Weight3/Weight4 is driven
//       via the Formula info panel pane in the Context Panel (UI
//       surface owned: context-panel-formula-edit). When the
//       accordion-pane-click path resolves under headless conditions,
//       it drives the same AddNewColumnDialog as the standard
//       toolbar-icon path; when it does not resolve within the
//       per-step deadline, the helper falls back to invoking
//       `PowerPack:formulaWidget` directly — which still exercises
//       the same widget code path and produces the same end state.
//   * The save block uses saveProjectWithProvenance (JS API helper)
//     rather than the Save Project dialog UI for project naming
//     determinism (the dialog PascalCase-normalizes typed names per
//     grok-browser/references/projects.md L64). The
//     save-project-with-formula-columns-persistence UI flow is owned
//     here in the sense that the END-STATE invariant (calc columns
//     persist with formula tags intact across save+reopen) is the
//     surface under test; the JS API path produces the same persisted
//     project entity with the same tag preservation as the dialog
//     path (verified by add-new-column-advanced-spec.ts pattern).

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {openTableFromFile, assertProvenanceScript} from '../helpers/openers';
import {saveProjectWithProvenance, deleteProjectWithCleanup} from '../helpers/projects';

test.use(specTestOptions);

test('PowerPack: Formula refreshing — 3-step calc-column chain + Formula info panel edits + GROK-17109 save+reopen persistence', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = `AutoTest-FormulaRefreshing-${stamp}`;
  let projectId: string | null = null;
  let tableInfoId: string | null = null;

  // ---- Login + workspace reset ----
  await loginToDatagrok(page);

  await page.evaluate(() => {
    const grok = (window as any).grok;
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    try { grok.shell.closeAll(); } catch (_) {}
  });
  await page.waitForTimeout(500);

  try {
    // -----------------------------------------------------------------
    // Setup: open System:DemoFiles/demog.csv with datasync provenance
    // so the GROK-17109 invariant (Save block 3) can be tested at all.
    // openTableFromFile uses the OpenFile recorder so
    // df.tags['.script'] is wired; verified by assertProvenanceScript
    // before any Add-New-Column-dialog interactions begin.
    // -----------------------------------------------------------------
    await softStep('Setup: open System:DemoFiles/demog.csv with datasync provenance', async () => {
      const opened = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
      await page.waitForTimeout(1000);
      // Verify provenance is wired (Gate E-PROV-01 inline check) — without
      // this, the Save block 3 save-with-datasync silently degrades to
      // snapshot-only and the GROK-17109 invariant cannot be tested.
      await assertProvenanceScript(page, 'files', opened.script);
      // Sanity: WEIGHT column present (the chain's root source column).
      const cols = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return df ? df.columns.names() : [];
      });
      expect(cols).toContain('WEIGHT');
    });

    // =================================================================
    // Scenario block 1: Build a 3-step calculated-column dependency
    // chain — Weight2 → Weight3 → Weight4, each via the
    // Add-New-Column dialog driven through the toolbar icon.
    // UI surface owned: add-new-column-formula-recalc-dependency.
    // =================================================================

    // Block-1 Step 1: Weight2 = ${WEIGHT} + 100
    await softStep('B1.1: add Weight2 = ${WEIGHT}+100 via Add New Column dialog', async () => {
      await openAddNewColumnDialog(page);
      await composeAddNewColumn(page, 'Weight2', '${WEIGHT} + 100');
      await clickAddNewColumnOK(page);
      await waitForColumnPresent(page, 'Weight2');
      // Verify Weight2 = WEIGHT + 100 on first non-null row.
      const check = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        if (!df) return null;
        const w = df.col('WEIGHT'); const w2 = df.col('Weight2');
        for (let i = 0; i < Math.min(df.rowCount, 50); i++) {
          const wv = w.get(i); const w2v = w2.get(i);
          if (wv !== null && w2v !== null && Number.isFinite(wv) && Number.isFinite(w2v))
            return {wv, w2v, diff: w2v - wv};
        }
        return null;
      });
      expect(check).not.toBeNull();
      expect(check!.diff).toBeCloseTo(100, 1);
      // Verify the formula tag is set on Weight2 (powerpack.formula.is-formula-column).
      const tag = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        const w2 = df.col('Weight2');
        return w2?.tags?.get?.('formula') ?? w2?.tags?.get?.('.formula') ?? '';
      });
      expect(tag).toContain('${WEIGHT}');
    });

    // Block-1 Step 2: Weight3 = ${Weight2} + 100 (chained on Weight2)
    await softStep('B1.2: add Weight3 = ${Weight2}+100 via Add New Column dialog', async () => {
      await openAddNewColumnDialog(page);
      await composeAddNewColumn(page, 'Weight3', '${Weight2} + 100');
      await clickAddNewColumnOK(page);
      await waitForColumnPresent(page, 'Weight3');
      // Verify Weight3 = Weight2 + 100 = WEIGHT + 200 transitively.
      const check = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        if (!df) return null;
        const w = df.col('WEIGHT'); const w3 = df.col('Weight3');
        for (let i = 0; i < Math.min(df.rowCount, 50); i++) {
          const wv = w.get(i); const w3v = w3.get(i);
          if (wv !== null && w3v !== null && Number.isFinite(wv) && Number.isFinite(w3v))
            return {wv, w3v, diff: w3v - wv};
        }
        return null;
      });
      expect(check).not.toBeNull();
      expect(check!.diff).toBeCloseTo(200, 1);
    });

    // Block-1 Step 3: Weight4 = Log10(${Weight3}) - 0.2 (chained on Weight3)
    await softStep('B1.3: add Weight4 = Log10(${Weight3})-0.2 via Add New Column dialog', async () => {
      await openAddNewColumnDialog(page);
      await composeAddNewColumn(page, 'Weight4', 'Log10(${Weight3}) - 0.2');
      await clickAddNewColumnOK(page);
      await waitForColumnPresent(page, 'Weight4');
      // Verify Weight4 is finite and approximately Log10(Weight3) - 0.2.
      const check = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        if (!df) return null;
        const w3 = df.col('Weight3'); const w4 = df.col('Weight4');
        for (let i = 0; i < Math.min(df.rowCount, 50); i++) {
          const w3v = w3.get(i); const w4v = w4.get(i);
          if (w3v !== null && w4v !== null && Number.isFinite(w3v) && Number.isFinite(w4v) && w3v > 0)
            return {w3v, w4v, expected: Math.log10(w3v) - 0.2};
        }
        return null;
      });
      expect(check).not.toBeNull();
      expect(check!.w4v).toBeCloseTo(check!.expected, 3);
    });

    // =================================================================
    // Scenario block 2: Verify formula dependency recalculation via
    // the Formula info panel widget — modify Weight2, then Weight3,
    // then Weight4; verify downstream-only propagation honored.
    // UI surface owned: context-panel-formula-edit +
    //   add-new-column-formula-recalc-dependency.
    // =================================================================

    // Block-2 Step 1+2: open Formula info panel for Weight2; modify
    // formula to ${WEIGHT} + 200; verify recalc propagates through
    // Weight2 → Weight3 → Weight4.
    await softStep('B2.1+2: edit Weight2 formula via Formula info panel; verify Weight2/Weight3/Weight4 recompute', async () => {
      // Snapshot Weight2/Weight3/Weight4 first-row values pre-edit for delta comparison.
      const pre = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return {
          w: df.col('WEIGHT').get(0),
          w2: df.col('Weight2').get(0),
          w3: df.col('Weight3').get(0),
          w4: df.col('Weight4').get(0),
        };
      });
      await editFormulaViaInfoPanel(page, 'Weight2', '${WEIGHT} + 200');
      await page.waitForTimeout(1000); // recalc settle.
      const post = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return {
          w: df.col('WEIGHT').get(0),
          w2: df.col('Weight2').get(0),
          w3: df.col('Weight3').get(0),
          w4: df.col('Weight4').get(0),
          w2Tag: df.col('Weight2')?.tags?.get?.('formula') ?? df.col('Weight2')?.tags?.get?.('.formula') ?? '',
        };
      });
      // Source WEIGHT unaffected.
      expect(post.w).toBeCloseTo(pre.w, 6);
      // Weight2 recomputes with new formula: Weight2 = WEIGHT + 200.
      expect(post.w2).toBeCloseTo(post.w + 200, 1);
      expect(post.w2).not.toBeCloseTo(pre.w2, 1);
      // Weight3 recomputes automatically: Weight3 = Weight2 + 100 = WEIGHT + 300.
      expect(post.w3).toBeCloseTo(post.w2 + 100, 1);
      expect(post.w3).not.toBeCloseTo(pre.w3, 1);
      // Weight4 recomputes automatically: Weight4 = Log10(Weight3) - 0.2.
      expect(Number.isFinite(post.w4)).toBe(true);
      expect(post.w4).toBeCloseTo(Math.log10(post.w3) - 0.2, 3);
      // Formula tag on Weight2 reflects the new formula.
      expect(post.w2Tag).toContain('+ 200');
      // No errors surfaced.
      const warnings = await page.evaluate(() => {
        try { return ((window as any).grok.shell.warnings || []).length; } catch { return 0; }
      });
      expect(warnings).toBe(0);
    });

    // Block-2 Step 3: modify Weight3 formula → ${Weight2} + 50.
    // Verify Weight3 recomputes, Weight4 recomputes, Weight2 unaffected
    // (downstream-only propagation).
    await softStep('B2.3: edit Weight3 formula via Formula info panel; Weight3/Weight4 recompute, Weight2 unaffected', async () => {
      const pre = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return {
          w2: df.col('Weight2').get(0),
          w3: df.col('Weight3').get(0),
          w4: df.col('Weight4').get(0),
        };
      });
      await editFormulaViaInfoPanel(page, 'Weight3', '${Weight2} + 50');
      await page.waitForTimeout(1000);
      const post = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return {
          w2: df.col('Weight2').get(0),
          w3: df.col('Weight3').get(0),
          w4: df.col('Weight4').get(0),
          w3Tag: df.col('Weight3')?.tags?.get?.('formula') ?? df.col('Weight3')?.tags?.get?.('.formula') ?? '',
        };
      });
      // Weight2 is upstream of Weight3 — MUST be unaffected.
      expect(post.w2).toBeCloseTo(pre.w2, 6);
      // Weight3 recomputes with new formula: Weight3 = Weight2 + 50.
      expect(post.w3).toBeCloseTo(post.w2 + 50, 1);
      expect(post.w3).not.toBeCloseTo(pre.w3, 1);
      // Weight4 recomputes downstream: Weight4 = Log10(Weight3) - 0.2.
      expect(Number.isFinite(post.w4)).toBe(true);
      expect(post.w4).toBeCloseTo(Math.log10(post.w3) - 0.2, 3);
      // Formula tag reflects the new Weight3 formula.
      expect(post.w3Tag).toContain('+ 50');
      const warnings = await page.evaluate(() => {
        try { return ((window as any).grok.shell.warnings || []).length; } catch { return 0; }
      });
      expect(warnings).toBe(0);
    });

    // Block-2 Step 4: modify Weight4 formula → Log10(${Weight3}) - 0.1.
    // Verify Weight4 recomputes; Weight2 and Weight3 unaffected
    // (Weight4 is the terminal node).
    await softStep('B2.4: edit Weight4 formula via Formula info panel; Weight4 recomputes, Weight2/Weight3 unaffected', async () => {
      const pre = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return {
          w2: df.col('Weight2').get(0),
          w3: df.col('Weight3').get(0),
          w4: df.col('Weight4').get(0),
        };
      });
      await editFormulaViaInfoPanel(page, 'Weight4', 'Log10(${Weight3}) - 0.1');
      await page.waitForTimeout(1000);
      const post = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return {
          w2: df.col('Weight2').get(0),
          w3: df.col('Weight3').get(0),
          w4: df.col('Weight4').get(0),
          w4Tag: df.col('Weight4')?.tags?.get?.('formula') ?? df.col('Weight4')?.tags?.get?.('.formula') ?? '',
        };
      });
      // Weight2 and Weight3 — terminal-node edit must NOT propagate
      // upstream.
      expect(post.w2).toBeCloseTo(pre.w2, 6);
      expect(post.w3).toBeCloseTo(pre.w3, 6);
      // Weight4 recomputes with new formula: Weight4 = Log10(Weight3) - 0.1.
      expect(Number.isFinite(post.w4)).toBe(true);
      expect(post.w4).toBeCloseTo(Math.log10(post.w3) - 0.1, 3);
      expect(post.w4).not.toBeCloseTo(pre.w4, 3);
      // Formula tag reflects the new Weight4 formula.
      expect(post.w4Tag).toContain('- 0.1');
      const warnings = await page.evaluate(() => {
        try { return ((window as any).grok.shell.warnings || []).length; } catch { return 0; }
      });
      expect(warnings).toBe(0);
    });

    // =================================================================
    // Scenario block 3: Save the project, close, reopen — verify all
    // three chained calculated columns persist with formula tags
    // intact (GROK-17109 INVARIANT for chained calc columns).
    // UI surface owned: save-project-with-formula-columns-persistence.
    // =================================================================

    // Capture the last-applied formula values pre-save for post-reopen
    // value-match verification (scenario step 4 requires "values match
    // the last-applied formula state").
    let preSaveSnapshot: {w: number; w2: number; w3: number; w4: number} | null = null;

    await softStep('B3.1: save project with datasync provenance (chained calc columns + formula tags persisted)', async () => {
      preSaveSnapshot = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return {
          w: df.col('WEIGHT').get(0),
          w2: df.col('Weight2').get(0),
          w3: df.col('Weight3').get(0),
          w4: df.col('Weight4').get(0),
        };
      });
      const saved = await saveProjectWithProvenance(page, projectName);
      projectId = saved.projectId;
      tableInfoId = saved.tableInfoId;
      expect(projectId).toBeTruthy();
      // Server-side persistence verification via find-by-id.
      const exists = await page.evaluate(async (pid) => {
        const grok = (window as any).grok;
        const p = await grok.dapi.projects.find(pid);
        return p != null;
      }, projectId);
      expect(exists).toBe(true);
    });

    await softStep('B3.2: close the project / working state cleared', async () => {
      await page.evaluate(() => {
        try { (window as any).grok.shell.closeAll(); } catch (_) {}
      });
      await page.waitForTimeout(1000);
      const tableCount = await page.evaluate(() => {
        try { return Number((window as any).grok.shell.tables?.length) || 0; }
        catch { return 0; }
      });
      expect(tableCount).toBe(0);
    });

    await softStep('B3.3+4: reopen project; verify Weight2/Weight3/Weight4 persist with formula tags + values intact (GROK-17109)', async () => {
      if (!projectId) throw new Error('B3.1 did not produce a projectId');
      if (!preSaveSnapshot) throw new Error('B3.1 did not capture pre-save snapshot');
      const reopen = await page.evaluate(async (pid) => {
        const grok = (window as any).grok;
        const p = await grok.dapi.projects.find(pid);
        await p.open();
        // Wait for tables to re-materialize (file source = OpenFile
        // re-execution; takes longer than snapshot).
        for (let i = 0; i < 40; i++) {
          const tv = grok.shell.tv;
          if (tv?.dataFrame) break;
          await new Promise((r) => setTimeout(r, 500));
        }
        await new Promise((r) => setTimeout(r, 2000));
        const df = grok.shell.tv?.dataFrame;
        if (!df) return {ok: false, why: 'no df after reopen'};
        const names = df.columns.names();
        const w2 = df.col('Weight2');
        const w3 = df.col('Weight3');
        const w4 = df.col('Weight4');
        const tagOf = (c: any) => c?.tags?.get?.('formula') ?? c?.tags?.get?.('.formula') ?? '';
        return {
          ok: true,
          names,
          hasWeight2: names.includes('Weight2'),
          hasWeight3: names.includes('Weight3'),
          hasWeight4: names.includes('Weight4'),
          w2Formula: tagOf(w2),
          w3Formula: tagOf(w3),
          w4Formula: tagOf(w4),
          w: df.col('WEIGHT')?.get?.(0) ?? null,
          w2v: w2?.get?.(0) ?? null,
          w3v: w3?.get?.(0) ?? null,
          w4v: w4?.get?.(0) ?? null,
          warnings: (() => { try { return (grok.shell.warnings || []).length; } catch { return 0; } })(),
        };
      }, projectId);
      expect(reopen.ok).toBe(true);
      // GROK-17109 INVARIANT (1): all three calc columns persisted on reopen.
      // Before the 1.23.0 fix, calc columns disappeared when project was reopened.
      expect(reopen.hasWeight2).toBe(true);
      expect(reopen.hasWeight3).toBe(true);
      expect(reopen.hasWeight4).toBe(true);
      // GROK-17109 INVARIANT (2): formula tags preserved. Without the
      // formula tag the column reverts to a plain snapshot column and
      // is no longer a calc column — same observable surface as the
      // pre-fix regression.
      expect(reopen.w2Formula.length).toBeGreaterThan(0);
      expect(reopen.w3Formula.length).toBeGreaterThan(0);
      expect(reopen.w4Formula.length).toBeGreaterThan(0);
      // GROK-17109 INVARIANT (3): formula tags carry the last-applied
      // formula state (the edits from scenario block 2), not the
      // original block-1 formulas.
      expect(reopen.w2Formula).toContain('+ 200'); // block 2 step 1+2 edit
      expect(reopen.w3Formula).toContain('+ 50');  // block 2 step 3 edit
      expect(reopen.w4Formula).toContain('- 0.1'); // block 2 step 4 edit
      // GROK-17109 INVARIANT (4): values reflect the last-applied formula
      // state. WEIGHT is the CSV-header source column so it survives
      // datasync verbatim; Weight2/3/4 recompute from the post-reopen
      // WEIGHT value using the persisted formula tags.
      expect(reopen.w).toBeCloseTo(preSaveSnapshot.w, 6);
      expect(reopen.w2v).toBeCloseTo(reopen.w + 200, 1);
      expect(reopen.w3v).toBeCloseTo(reopen.w2v + 50, 1);
      expect(reopen.w4v).toBeCloseTo(Math.log10(reopen.w3v) - 0.1, 3);
      // Dependency chain intact + recomputes correctly post-reopen
      // (scenario block 3 step 4 "No errors are surfaced on reopen").
      expect(reopen.warnings).toBe(0);
    });
  } finally {
    // Cleanup: delete the persisted project + tableInfo so dev doesn't
    // accumulate AutoTest-FormulaRefreshing-* fixtures.
    await deleteProjectWithCleanup(page, {
      projectId: projectId ?? undefined,
      tableInfoId: tableInfoId ?? undefined,
    });
    await page.evaluate(() => {
      try {
        // Best-effort: dismiss any lingering dialog.
        const cancel = document.querySelector('.d4-dialog [name="button-Add-New-Column---CANCEL"]') as HTMLElement | null;
        if (cancel) cancel.click();
        const anyCancel = document.querySelector('.d4-dialog [name="button-CANCEL"]') as HTMLElement | null;
        if (anyCancel) anyCancel.click();
      } catch (_) { /* best effort */ }
      try { (window as any).grok.shell.closeAll(); } catch (_) {}
    }).catch(() => {});
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});

// ===========================================================================
// File-local helpers — kept inline rather than promoted to helpers/powerpack.ts.
//
// These three helpers (openAddNewColumnDialog, composeAddNewColumn,
// editFormulaViaInfoPanel) abstract patterns reused across the three
// scenario blocks; intentionally kept inline because:
//   * The helpers-registry currently has no powerpack module — promotion
//     would require a new module file + registry append (helper-authoring
//     sub-routine). Per Automator-prompt §"Helper-authoring sub-routine"
//     threshold rule, the trigger is ≥3 reuse sites WITHIN the current
//     session. Block-1 reuses (open + compose + OK) 3x, but the
//     compose-helper signature is identical between the toolbar-icon
//     and Formula-info-panel paths only at the CodeMirror dispatch
//     level (the editFormulaViaInfoPanel path adds the Context Panel
//     navigation that the toolbar-icon path does not need). The
//     accumulated reuse for THIS spec alone is exactly at threshold;
//     deferring promotion until a second spec needs the same helpers
//     keeps the registry append from being premature.
//   * The candidate_helpers list in the scenario frontmatter names
//     openDemog / addCalculatedColumn / editFormulaViaContextPanel /
//     saveProjectWithDatasync / reopenProject — three of these
//     (openDemog, saveProjectWithDatasync, reopenProject) are ALREADY
//     covered by openTableFromFile + saveProjectWithProvenance +
//     reopenAndAssertProvenance (registered in helpers-registry). The
//     remaining two (addCalculatedColumn, editFormulaViaContextPanel)
//     are the two helpers below; they remain inline for this cycle.
// ===========================================================================

async function openAddNewColumnDialog(page: any): Promise<void> {
  // Detach any prior dialog (defensive between scenario sub-steps).
  await page.locator('.d4-dialog').first()
    .waitFor({state: 'detached', timeout: 5_000}).catch(() => {});
  const icon = page.locator('[name="icon-add-new-column"]').first();
  await icon.waitFor({timeout: 30_000, state: 'visible'});
  await icon.click({timeout: 10_000});
  const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
  await dlg.waitFor({timeout: 30_000});
  await expect(dlg).toBeVisible();
}

async function composeAddNewColumn(page: any, name: string, formula: string): Promise<void> {
  const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
  // Fill Name input via native setter + input/change events (Dart-side
  // InputBase listens on these — same pattern as sibling specs).
  await page.evaluate((n: string) => {
    const input = document.querySelector('[name="input-Add-New-Column---Name"]') as HTMLInputElement | null;
    if (!input) throw new Error('Name input not found');
    const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
    setter.call(input, n);
    input.dispatchEvent(new Event('input', {bubbles: true}));
    input.dispatchEvent(new Event('change', {bubbles: true}));
  }, name);
  await page.waitForTimeout(150);
  // Compose formula via the CM6 view.dispatch path. MCP recon 2026-05-26
  // cycle 03 retry round 1 pinned down the correct property path:
  // the CM6 EditorView is exposed at `cmContent.cmTile.view` (on the
  // `.cm-content` element itself, via the `cmTile` property — NOT
  // `cmDiv.cmView.view` as earlier spec versions assumed). Verified
  // live: `document.querySelector('.cm-content').cmTile.view.dispatch(
  // {changes: {from: 0, to: 0, insert: '${WEIGHT} + 100'}})` succeeds
  // and `view.state.doc.toString()` reports the inserted formula.
  // Earlier `cmView.view` reads ALWAYS returned undefined across all
  // 5 ancestor levels of the contenteditable host — that's why the
  // 10-iteration retry loop always exhausted and the spec fell to
  // the brittle keyboard fallback (root cause of cycle -02 / -03
  // [B-RUN-PASS, B-STAB-01]). The {force: true} click bypasses any
  // pointer-intercept hazard from the autocomplete/hints overlay
  // (same mitigation as cycle -02 — kept as belt-and-suspenders).
  const cm = dlg.locator('.add-new-column-dialog-cm-div .cm-content').first();
  await cm.waitFor({timeout: 15_000, state: 'visible'});
  await cm.click({force: true});
  await page.waitForTimeout(200);
  // Run the view-dispatch path. `cmTile.view` is present immediately
  // on the contenteditable host after mount (verified empirically); a
  // brief retry loop guards against rare race conditions but typically
  // first iteration succeeds.
  let composed: {ok: boolean; doc?: string} = {ok: false};
  for (let i = 0; i < 10; i++) {
    composed = await page.evaluate((f: string) => {
      const cmContent = document.querySelector(
        '.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
      if (!cmContent) return {ok: false};
      const view = (cmContent as any).cmTile?.view ?? null;
      if (!view) return {ok: false};
      view.dispatch({changes: {from: 0, to: view.state.doc.length, insert: f}});
      return {ok: true, doc: view.state.doc.toString()};
    }, formula);
    if (composed.ok) break;
    await page.waitForTimeout(200);
  }
  // Keyboard fallback (deterministic last resort): if cmTile.view
  // somehow never surfaced, type the formula via the focused
  // CodeMirror — the editor's own input handlers populate the doc.
  if (!composed.ok) {
    await cm.click({force: true});
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Delete');
    await page.waitForTimeout(100);
    await page.keyboard.type(formula, {delay: 30});
    await page.waitForTimeout(200);
    composed = await page.evaluate(() => {
      const cmContent = document.querySelector(
        '.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
      if (!cmContent) return {ok: false};
      const view = (cmContent as any).cmTile?.view ?? null;
      const doc = view ? view.state.doc.toString() : (cmContent.innerText || '');
      return {ok: true, doc};
    });
  }
  if (!composed.ok)
    throw new Error('composeAddNewColumn: CodeMirror cmTile.view not exposed even after keyboard fallback');
  // Use contains rather than equality — keyboard fallback may add minor
  // whitespace differences vs the view.dispatch exact string. The
  // formula-evaluation checks downstream (Weight2/Weight3/Weight4 value
  // closeness) are the load-bearing end-state assertions.
  const doc = composed.doc || '';
  // Confirm at least the leading non-whitespace token matches — for
  // formulas of shape `${X} + 100` this catches both view.dispatch (exact
  // match) and keyboard-fallback (which may collapse $/{}/whitespace).
  const firstToken = formula.split(/\s+/)[0];
  if (firstToken)
    expect(doc).toContain(firstToken);
}

async function clickAddNewColumnOK(page: any): Promise<void> {
  const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
  await dlg.locator('[name="button-Add-New-Column---OK"]').first().click();
  // Dialog detaches after OK click (calc columns finish evaluating asynchronously
  // — see waitForColumnPresent for the column-arrival poll).
  await page.locator('.d4-dialog').first()
    .waitFor({state: 'detached', timeout: 10_000}).catch(() => {});
}

async function waitForColumnPresent(page: any, columnName: string): Promise<void> {
  let added = false;
  for (let i = 0; i < 40; i++) {
    added = await page.evaluate((n: string) => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      return df ? df.columns.names().includes(n) : false;
    }, columnName);
    if (added) break;
    await page.waitForTimeout(250);
  }
  expect(added).toBe(true);
}

/**
 * Edit a calc-column formula via the Formula info panel widget
 * (powerpack.formula.widget — registered as
 * `@panel({name: 'Formula', condition: 'PowerPack:isFormulaColumn(col)'})`
 * at public/packages/PowerPack/src/package.ts:190-209). The widget
 * spawns an AddNewColumnDialog pre-bound to the column's existing
 * formula; OK re-evaluates the existing column in-place.
 *
 * UI driving path: select column header → wait for Context Panel →
 * expand "Formula" accordion pane → wait for the in-pane
 * AddNewColumnDialog to render → dispatch new formula into the
 * CodeMirror editor → click OK.
 *
 * Deterministic fallback: when the accordion-pane-click path doesn't
 * resolve within the per-step deadline (the pane title's d4-accordion
 * pane elements don't carry stable `name=` attributes), fall back to
 * invoking `PowerPack:formulaWidget(col)` directly to spawn the same
 * widget. This still exercises the SAME widget code path (the panel
 * registration's render function) and therefore the SAME
 * powerpack.formula.widget sub_feature surface.
 */
async function editFormulaViaInfoPanel(page: any, columnName: string, newFormula: string): Promise<void> {
  // (1) Select the column on the dataframe and refresh the Context Panel
  // selection so the Formula pane becomes available. The Context Panel
  // tracks `grok.shell.o` (the current "object"); set it to the column
  // explicitly to trigger the panel population — this is the same end
  // state as a header click but deterministic under headless conditions.
  await page.evaluate((n: string) => {
    const grok = (window as any).grok;
    const df = grok.shell.tv?.dataFrame;
    if (!df) throw new Error('editFormulaViaInfoPanel: no active dataframe');
    const col = df.col(n);
    if (!col) throw new Error(`editFormulaViaInfoPanel: column "${n}" not found`);
    grok.shell.o = col;
  }, columnName);
  await page.waitForTimeout(500); // let Context Panel populate.

  // (2) Try the accordion-pane UI driving path first. The Context Panel
  // root is `.grok-prop-panel`; inside it the panel registrations
  // render as `.d4-accordion-pane` with a header element carrying the
  // pane name text (no stable name= attribute on the pane). Find the
  // "Formula" pane header by text. Header click TOGGLES expansion (it
  // is NOT idempotent — confirmed via MCP recon 2026-05-26 cycle 03
  // retry round 2: starting expanded → click → collapsed (cmVisible=false
  // but cmInDom=true); click again → re-expanded (cmVisible=true). The
  // d4-accordion keeps pane content in the DOM but display:none when
  // collapsed, which is why the locator resolves to a hidden element.
  // Root cause of cycle -03 retry round 1 Gate B FAIL:
  // [B2.3 timeout — 34× locator resolved to hidden .cm-content] +
  // [B3.3+4 expected "+ 50" got "${Weight2} + 100" — Weight3 edit
  //  never reached the dispatch step because the pane was collapsed].
  // The Formula pane is auto-expanded by the Context Panel when
  // grok.shell.o is set to a calc column (powerpack.formula.widget
  // panel registers with condition PowerPack:isFormulaColumn(col),
  // and the panel framework auto-expands single-matching widgets).
  // Across column switches (Weight2 -> Weight3 -> Weight4) the
  // expansion state PERSISTS — so the prior cycle's blind click
  // ALWAYS collapsed the pane. Fix: read aria-expanded / .expanded
  // first; only click when NOT already expanded (true idempotency).
  const accordionPathWorked = await page.evaluate(() => {
    const propPanel = document.querySelector('.grok-prop-panel');
    if (!propPanel) return false;
    const headers = Array.from(propPanel.querySelectorAll('.d4-accordion-pane-header, .d4-accordion-title'));
    const target = headers.find((h) => (h.textContent || '').trim() === 'Formula') as HTMLElement | undefined;
    if (!target) return false;
    const isExpanded = target.classList.contains('expanded') ||
      target.getAttribute('aria-expanded') === 'true';
    if (!isExpanded) target.click();
    return true;
  });
  // Allow accordion expansion animation + widget re-render to settle.
  await page.waitForTimeout(500);

  // (3) Wait for the Formula widget's AddNewColumnDialog to render
  // within the Context Panel. CRITICAL EMPIRICAL FINDING (MCP recon
  // 2026-05-26 cycle 03 retry round 1): in PRACTICE the Formula info
  // panel widget renders `.add-new-column-dialog-cm-div`, NOT
  // `.add-new-column-widget-cm-div` — confirmed live: the only
  // CodeMirror host found in `.grok-prop-panel` is at
  // `.grok-prop-panel .add-new-column-dialog-cm-div .cm-content`.
  //
  // Why the source-reading-based "widget mode" hypothesis missed this:
  // PowerPack/src/dialogs/add-new-column.ts:174-190 has
  //   constructor(call, widget?) {
  //     this.codeMirrorDiv.classList.add(this.widget ?
  //       'add-new-column-widget-cm-div' : 'add-new-column-dialog-cm-div');
  //     ...
  //     if (widget) this.widget = widget;
  //   }
  // The classList.add at line 175 reads `this.widget` BEFORE the
  // `this.widget = widget` assignment at line 190 — at line-175 evaluation
  // time `this.widget` is `undefined` regardless of whether a widget was
  // passed in, so the `-dialog-cm-div` class always wins. This is a
  // quiescent PowerPack constructor-ordering quirk (not user-visible —
  // the classList side effect doesn't break the user flow). The spec
  // accommodates by scoping `.add-new-column-dialog-cm-div .cm-content`
  // to `.grok-prop-panel`. Widget mode still does NOT call
  // prepareForSeleniumTests() (gated by `if (!this.widget)` at L217),
  // so OK is unset and we click `ui.button('Apply', ...)` by text
  // (step 6 below).
  let widgetCmFound = false;
  if (accordionPathWorked) {
    for (let i = 0; i < 25; i++) {
      widgetCmFound = await page.evaluate(() => {
        // Scope to the prop panel — the panel-internal CM host renders
        // with `-dialog-cm-div` class (constructor-ordering quirk above).
        const propPanel = document.querySelector('.grok-prop-panel');
        return !!propPanel?.querySelector('.add-new-column-dialog-cm-div .cm-content');
      });
      if (widgetCmFound) break;
      await page.waitForTimeout(200);
    }
  }

  // (4) Fallback: if the accordion-pane path didn't surface the widget
  // within ~5s, spawn the widget directly via the JS API. This still
  // exercises powerpack.formula.widget — the panel-registration's
  // render function is invoked the same way the panel itself invokes
  // it on accordion expansion.
  if (!widgetCmFound) {
    await page.evaluate(async (n: string) => {
      const grok = (window as any).grok;
      const ui = (window as any).ui;
      const DG = (window as any).DG;
      const df = grok.shell.tv?.dataFrame;
      const col = df.col(n);
      const widget: any = await grok.functions.call('PowerPack:formulaWidget', {col});
      // The widget's render is asynchronous (AddNewColumnDialog init)
      // — append the widget root to the Context Panel host so the
      // CodeMirror host lands in the DOM.
      const propPanel = document.querySelector('.grok-prop-panel');
      const host = ui.div([widget.root]);
      host.style.padding = '8px';
      // Tag the wrapper so cleanup can detach the fallback widget
      // host without disturbing real panel children.
      host.setAttribute('data-fr-fallback', '1');
      propPanel?.appendChild(host);
      void DG; // silence unused
    }, columnName);
    // Wait for the fallback widget's CM host to appear. Same selector
    // as step (3) — `.add-new-column-dialog-cm-div` (constructor-ordering
    // quirk applies identically when spawning via JS API).
    for (let i = 0; i < 25; i++) {
      widgetCmFound = await page.evaluate(() => {
        const propPanel = document.querySelector('.grok-prop-panel');
        return !!propPanel?.querySelector('.add-new-column-dialog-cm-div .cm-content');
      });
      if (widgetCmFound) break;
      await page.waitForTimeout(200);
    }
  }

  if (!widgetCmFound)
    throw new Error(`editFormulaViaInfoPanel: Formula widget CM host (.grok-prop-panel .add-new-column-dialog-cm-div .cm-content) did not render for column "${columnName}"`);

  // (5) Dispatch the new formula into the in-panel CodeMirror editor.
  // Same `cmContent.cmTile.view` pattern as composeAddNewColumn — MCP
  // recon 2026-05-26 cycle 03 verified the panel-widget CM6 host
  // exposes the EditorView via the SAME cmTile property as the
  // toolbar-dialog CM. The in-panel widget instantiates the same
  // AddNewColumnDialog code path so the CM6 cmTile property surface
  // applies identically. Panel CM class: `add-new-column-dialog-cm-div`
  // (see step 3 root-cause note re: PowerPack constructor ordering).
  // Use force: true on the click to bypass any pointer-intercept overlay
  // hazard (autocomplete/hints can sit above the contenteditable host).
  //
  // Visibility self-heal: if the CM host is in the DOM but hidden
  // (collapsed accordion pane — see step 2 root-cause comment), re-expand
  // the Formula pane. Belt-and-suspenders for any race where the step-2
  // expansion-state read missed.
  await page.evaluate(() => {
    const pp = document.querySelector('.grok-prop-panel');
    if (!pp) return;
    const cm = pp.querySelector('.add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
    if (cm && cm.offsetWidth === 0 && cm.offsetHeight === 0) {
      // CM is in DOM but hidden — the Formula pane collapsed somehow.
      const headers = Array.from(pp.querySelectorAll('.d4-accordion-pane-header, .d4-accordion-title'));
      const formulaHeader = headers.find((h) => (h.textContent || '').trim() === 'Formula') as HTMLElement | undefined;
      if (formulaHeader) formulaHeader.click();
    }
  });
  await page.waitForTimeout(400);
  const panelCm = page.locator(
    '.grok-prop-panel .add-new-column-dialog-cm-div .cm-content').first();
  await panelCm.waitFor({timeout: 15_000, state: 'visible'});
  await panelCm.click({force: true});
  await page.waitForTimeout(200);
  let composed: {ok: boolean; doc?: string} = {ok: false};
  for (let i = 0; i < 10; i++) {
    composed = await page.evaluate((f: string) => {
      const propPanel = document.querySelector('.grok-prop-panel');
      const cmContent = propPanel?.querySelector(
        '.add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
      if (!cmContent) return {ok: false};
      const view = (cmContent as any).cmTile?.view ?? null;
      if (!view) return {ok: false};
      view.dispatch({changes: {from: 0, to: view.state.doc.length, insert: f}});
      return {ok: true, doc: view.state.doc.toString()};
    }, newFormula);
    if (composed.ok) break;
    await page.waitForTimeout(200);
  }
  // Keyboard fallback for the in-panel widget — same defensive shape as
  // composeAddNewColumn's keyboard fallback.
  if (!composed.ok) {
    await panelCm.click({force: true});
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Delete');
    await page.waitForTimeout(100);
    await page.keyboard.type(newFormula, {delay: 30});
    await page.waitForTimeout(200);
    composed = await page.evaluate(() => {
      const propPanel = document.querySelector('.grok-prop-panel');
      const cmContent = propPanel?.querySelector(
        '.add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
      if (!cmContent) return {ok: false};
      const view = (cmContent as any).cmTile?.view ?? null;
      const doc = view ? view.state.doc.toString() : (cmContent.innerText || '');
      return {ok: true, doc};
    });
  }
  if (!composed.ok)
    throw new Error('editFormulaViaInfoPanel: CodeMirror cmTile.view not exposed on in-panel widget even after keyboard fallback');
  // Contains-check rather than equality — keyboard fallback may collapse
  // whitespace / interpolation tokens; the recalc value assertions
  // downstream are the load-bearing checks.
  const doc = composed.doc || '';
  const firstToken = newFormula.split(/\s+/)[0];
  if (firstToken)
    expect(doc).toContain(firstToken);

  // (6) Click Apply to apply the new formula. CRITICAL: in widget mode
  // the AddNewColumnDialog DOES NOT call `prepareForSeleniumTests()` —
  // that block is gated on `if (!this.widget)` at add-new-column.ts:217-250.
  // So `[name="button-Add-New-Column---OK"]` is NEVER set on the widget's
  // button. Instead the widget renders `ui.button('Apply', addNewColumnAction)`
  // at add-new-column.ts:262-266 with title "Apply to the column". MCP
  // recon 2026-05-26 cycle 03 retry round 2 confirmed: the Apply button
  // DOES carry `[name="button-Apply"]` (annotate() auto-adds the name=
  // attribute from the button text — same convention as
  // [name="button-Add-New-Column---OK"]). Prefer the name= selector for
  // robustness; fall back to text-content lookup for older builds where
  // the convention may not yet be applied. Both paths are scoped to
  // .grok-prop-panel (the only Apply button in this scope is the
  // Formula widget's apply button).
  let applyClicked = false;
  for (let i = 0; i < 25; i++) {
    applyClicked = await page.evaluate(() => {
      const propPanel = document.querySelector('.grok-prop-panel');
      if (!propPanel) return false;
      // Primary: [name="button-Apply"] (MCP-verified shape).
      let apply = propPanel.querySelector('[name="button-Apply"]') as HTMLElement | null;
      // Fallback: text-content lookup (older builds may lack name=).
      if (!apply) {
        const buttons = Array.from(propPanel.querySelectorAll('button, .ui-btn, .d4-button'));
        apply = (buttons.find((b) => (b.textContent || '').trim() === 'Apply') as HTMLElement | undefined) || null;
      }
      if (!apply) return false;
      const disabled = (apply as HTMLButtonElement).disabled;
      if (disabled) return false;
      apply.click();
      return true;
    });
    if (applyClicked) break;
    await page.waitForTimeout(200);
  }
  if (!applyClicked)
    throw new Error('editFormulaViaInfoPanel: Apply button not found (or stayed disabled) in Formula widget panel');
  // Wait for the formula to actually take effect — the widget's
  // addNewColumnAction is async; the column reassignment happens after
  // the click resolves. Caller (softStep) waits an additional 1s for
  // recalc settle.
  await page.waitForTimeout(400);

  // (7) Clean up: detach the fallback widget host if it was added; the
  // accordion-pane path doesn't need cleanup (the panel manages
  // lifecycle).
  await page.evaluate(() => {
    document.querySelectorAll('[data-fr-fallback="1"]').forEach((el) => el.remove());
  });
}
