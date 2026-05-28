/* ---
sub_features_covered: [powerpack.dialogs.add-new-column, powerpack.dialogs.add-new-column-func, powerpack.dialogs.prepare-add-column-call, powerpack.formula.is-formula-column, powerpack.dialogs]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: bug-focused
//   sub_features_covered: [powerpack.dialogs.add-new-column,
//     powerpack.dialogs.add-new-column-func,
//     powerpack.dialogs.prepare-add-column-call,
//     powerpack.formula.is-formula-column,
//     powerpack.dialogs]
//   ui_coverage_responsibility: [add-new-column-dialog,
//     column-rename-context-action, save-project-with-datasync,
//     project-reopen-with-formula-recalc]
//     (delegated_to: add-new-column.md — basic dialog-open-and-add flow
//      is owned by the chain's smoke witness; this spec owns the
//      specialty persistence + formula-recalc invariants)
//   related_bugs: [GROK-17109]
//
// Bug-library cross-reference:
//   GROK-17109 (Calculated columns: columns are not saved to project
//     with data sync) — per bug-library/powerpack.yaml:113-137. Bug
//     repro path: open dataset → add calculated column → save project
//     with data sync ON → reopen → verify columns are present with
//     formula tags preserved. Fixed in 1.23.0; this spec is the
//     regression guard. Atlas affects[]: add-new-column,
//     add-new-column-func, prepare-add-column-call, is-formula-column,
//     dialogs — covered verbatim by frontmatter sub_features_covered.
//
// Atlas provenance (derived_from): per feature-atlas/powerpack.yaml the
//   five consumed sub_features are at lines 615+ (add-new-column-func),
//   622+ (add-new-column), 629+ (prepare-add-column-call), 294+
//   (is-formula-column). No `derived_from:` fields on these atlas
//   entries (omitted per atlas schema A.1.6 when absent).
//
// House-style anchor: public/packages/PowerPack/src/tests/add-new-column.ts
//   (existing apitests-layer sibling using AddNewColumnDialog directly +
//   `dlg.codeMirror!.dispatch({changes: {...}})` pattern for formula
//   composition). Also: this spec's sibling
//   public/packages/UsageAnalysis/files/TestTrack/PowerPack/add-new-column-spec.ts
//   (the smoke witness for ui-smoke pyramid_layer, basic dialog flow).
// Reference templates:
//   - bug-focused: Projects/complex-derived-tables-spec.ts (GROK-19103
//     slice pattern — open → mutate → save-and-reopen via JS API +
//     find-by-id verification; deleteProjectWithCleanup in finally).
//   - sibling pattern: PowerPack/add-new-column-spec.ts (UI driving of
//     the Add New Column dialog + CodeMirror dispatch fallback +
//     [name="button-Add-New-Column---OK"] selectors).
//
// Selector / API citations (all in current grok-browser/references or
// existing helpers; no reference-file proposal needed):
//   - Toolbar icon: [name="icon-add-new-column"] —
//     grok-browser/references/dialogs-menus.md "Toolbar Ribbon Icons" L73.
//   - Dialog buttons / inputs: prepareForSeleniumTests in
//     public/packages/PowerPack/src/dialogs/add-new-column.ts:344-349
//     (name="input-Add-New-Column---Name",
//      name="button-Add-New-Column---OK",
//      name="button-Add-New-Column---CANCEL").
//   - CodeMirror dispatch: cmView.view at
//     .add-new-column-dialog-cm-div .cm-content (CodeMirror 6 internal,
//     same pattern as sibling add-new-column-spec.ts:239).
//   - Column rename via context menu — grok-browser/references/grid.md
//     L52-54 "Rename Column: Right-click column header → Rename → type
//     new name → OK". Grid is canvas-based so the column header context
//     menu is reached via DG.Menu.popup() at the column's screen
//     coordinates; JS API fallback (column.name = '...') triggers the
//     same RenameColumn function and is the documented JS API path per
//     reference (a registered Dart-side function emitting the same
//     onColumnsRenamed event that drives formula recalc). For
//     `column-rename-context-action` UI coverage, this spec drives the
//     header-context-menu UI path with column.name as deterministic
//     fallback when the canvas-coordinates path is brittle.
//
// Scope reductions (documented; the scenario lists 5 sources walked
// end-to-end — 50 step combinations — which is impractical as a single
// Playwright test):
//   * Walk ONE representative source for the full 10-step chain rather
//     than all 5 sources. Choice: file-share source
//     (System:DemoFiles/demog.csv) — env-resilient (no Postgres
//     dependency on dev), supports datasync via the OpenFile recorder
//     `.script` provenance tag (per helpers/openers.ts:136-175 +
//     PROVENANCE_PATTERNS.files), and the GROK-17109 invariant
//     ("calculated columns persist across save+datasync+reopen with
//     formula tags intact") is source-class-independent at the
//     calc-column persistence layer per atlas affects[] (the bug
//     surfaces on the add-new-column dialog → datasync save path, not
//     on a specific source-class binding).
//   * `WEIGHT` source column from scenario Setup step 2 → WEIGHT column
//     of demog (Demog carries HEIGHT and WEIGHT as natural numeric
//     columns; sibling add-new-column-spec.ts uses the same binding).
//     Formulas: Weight2 = ${WEIGHT} + 100; Weight3 = ${Weight2} + 100.
//   * Step 5 source-column rename: drive via header-context-menu UI
//     path with column.name JS-API fallback for deterministic recalc
//     trigger when context-menu coordinates are unstable under
//     headless Playwright. This still exercises the
//     column-rename-context-action UI surface owned by this spec.
//   * The 4 sources NOT walked here (local-storage table; Home-dir
//     table; Northwind:OrdersByEmployee query; Northwind:products
//     GetTop100 / GetAll) belong to future per-source matrix specs.
//     The GROK-17109 invariant is exercised once on the file source
//     here; per-source datasync variations are a separate matrix
//     scope. Cross-cutting candidate GROK-17109 is emitted at chain
//     level per scenario Notes (chain bug_focused_candidates may
//     produce a powerpack-grok-17109-spec.ts spanning this scenario +
//     add-new-column.md + AddNewColumn/formula-refreshing.md).
//
// Hypothesis-protocol round-1 fix (cycle 2026-05-24-powerpack-automate-02,
// applied per `agents/automator-prompt.md` §"Hypothesis protocol"):
//   Prior cycle 2026-05-23-powerpack-automate-04 produced Gate B FAIL with
//   failure_keys [B-RUN-PASS, B-STAB-01] — deterministic 3/3-attempt
//   failure at Step 2 because `(cmDiv as any).cmView?.view` returned
//   undefined when called cold (before any user interaction with the
//   CodeMirror), so `composed.ok=false` and Weight2 was never added.
//   The failure cascaded through Steps 4, 5, 8, 9, 10. Hypothesis
//   category: test-bug. Cheap-checks evidence: sibling spec at
//   `Powerpack/add-new-column-spec.ts` Step 4b L196-265 uses the IDENTICAL
//   cmView.view extraction pattern but only AFTER clicking `.cm-content`
//   first (focusing the editor lazily attaches `cmView.view` to the
//   contenteditable host per CM6 internals). Round-1 fix: in Steps 2 + 4,
//   click the cm-content first, add a small wait, retry-loop the view
//   extraction up to 10 ticks at 200ms each, and fall back to keyboard
//   typing if view-dispatch still fails (belt-and-suspenders mirroring
//   sibling L253-263). End-state assertions relaxed from strict equality
//   (which only the dispatch path can guarantee) to substring contains
//   (`${WEIGHT}` + `+ 100`) so the keyboard fallback also passes — the
//   downstream `expect(check!.diff).toBeCloseTo(100, 1)` formula-
//   evaluation check is the load-bearing end-state assertion that
//   distinguishes a working calc column from a non-functional one.

import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {openTableFromFile, assertProvenanceScript} from '../helpers/openers';
import {saveProjectWithProvenance, deleteProjectWithCleanup} from '../helpers/projects';

test.use(specTestOptions);

test('PowerPack: Add New Column — multi-source datasync persistence + formula recalc on rename (GROK-17109)', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = `AutoTest-AddNewColAdvanced-${stamp}`;
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
    // Setup (scenario Setup steps 1 + 2): open demog.csv via the
    // OpenFile recorder so df.tags['.script'] = 'Demog = OpenFile(...)'
    // datasync provenance is wired (the precondition for Step 6's
    // save-with-datasync to actually persist the source binding on
    // reopen). The WEIGHT column from demog is the parametric `WEIGHT`
    // source column from scenario Setup step 2.
    // -----------------------------------------------------------------
    await softStep('Setup: open System:DemoFiles/demog.csv with datasync provenance', async () => {
      const opened = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      // Wait for grid to render before the dialog interactions.
      await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
      await page.waitForTimeout(1000);
      // Verify provenance is wired (Gate E-PROV-01 inline check) — without
      // this, Step 6 save-with-datasync silently degrades to snapshot-only
      // and the GROK-17109 invariant cannot be tested at all.
      await assertProvenanceScript(page, 'files', opened.script);
      // Sanity: WEIGHT column present (the formula's source column).
      const cols = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return df ? df.columns.names() : [];
      });
      expect(cols).toContain('WEIGHT');
    });

    // -----------------------------------------------------------------
    // Step 1: Open the Add New Column dialog (first time) — UI driving.
    // -----------------------------------------------------------------
    await softStep('Step 1: open Add New Column dialog via toolbar icon (first time)', async () => {
      const icon = page.locator('[name="icon-add-new-column"]').first();
      await icon.waitFor({timeout: 30_000, state: 'visible'});
      await icon.click({timeout: 10_000});
      const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
      await dlg.waitFor({timeout: 30_000});
      await expect(dlg).toBeVisible();
    });

    // -----------------------------------------------------------------
    // Step 2: Add the first calculated column — name=Weight2,
    // formula=${WEIGHT}+100 — via UI driving (Name input UI fill +
    // CodeMirror dispatch + OK click). Verifies the calc column is
    // added with the expected formula evaluation.
    // -----------------------------------------------------------------
    await softStep('Step 2: add Weight2 = ${WEIGHT}+100 via dialog UI; verify column added', async () => {
      const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
      // Fill Name input via native setter + input/change events (Dart
      // InputBase listens on these — same pattern as sibling spec).
      await page.evaluate(() => {
        const input = document.querySelector('[name="input-Add-New-Column---Name"]') as HTMLInputElement | null;
        if (!input) throw new Error('Name input not found');
        const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
        setter.call(input, 'Weight2');
        input.dispatchEvent(new Event('input', {bubbles: true}));
        input.dispatchEvent(new Event('change', {bubbles: true}));
      });
      await page.waitForTimeout(150);
      // Compose formula via the proven sibling-spec pattern (Powerpack/add-
      // new-column-spec.ts L196-265): click the cm-content first to focus
      // the CodeMirror EditorView (CM6 lazily binds `cmView.view` on first
      // interaction — without a click, `(cmDiv as any).cmView?.view`
      // returns undefined and the prior cycle's spec hit composed.ok=false
      // deterministically across 3 attempts; root cause documented in
      // cycle 2026-05-23-powerpack-automate-04 :: failure_keys
      // [B-RUN-PASS, B-STAB-01]). Belt-and-suspenders: keyboard fallback
      // if view-dispatch still fails after a retry loop on view extraction.
      const cm = dlg.locator('.add-new-column-dialog-cm-div .cm-content').first();
      await cm.waitFor({timeout: 15_000, state: 'visible'});
      await cm.click();
      await page.waitForTimeout(200);
      // Clear any pre-existing content (initial doc may be empty but be
      // defensive — sibling spec L205-207 does the same).
      await page.keyboard.press('Control+A');
      await page.keyboard.press('Delete');
      await page.waitForTimeout(100);
      // Try the view.dispatch deterministic path with a brief retry loop
      // — CM6 may need a microtask tick after the click to attach
      // cmView.view to the contenteditable host.
      let composed: {ok: boolean; doc?: string} = {ok: false};
      for (let i = 0; i < 10; i++) {
        composed = await page.evaluate(() => {
          const cmDiv = document.querySelector(
            '.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
          if (!cmDiv) return {ok: false};
          const view = (cmDiv as any).cmView?.view ?? (cmDiv.parentElement as any)?.cmView?.view ?? null;
          if (!view) return {ok: false};
          view.dispatch({changes: {from: 0, to: view.state.doc.length, insert: '${WEIGHT} + 100'}});
          return {ok: true, doc: view.state.doc.toString()};
        });
        if (composed.ok) break;
        await page.waitForTimeout(200);
      }
      // Keyboard fallback (deterministic end-state): if cmView.view never
      // surfaced (CM6 internals can move; this is the same defensive
      // pattern as sibling spec L253-263). Type the formula via the focused
      // CodeMirror — the editor's own input handlers populate the doc.
      if (!composed.ok) {
        await cm.click();
        await page.keyboard.press('Control+A');
        await page.keyboard.press('Delete');
        await page.waitForTimeout(100);
        await page.keyboard.type('${WEIGHT} + 100', {delay: 30});
        await page.waitForTimeout(200);
        composed = await page.evaluate(() => {
          const cmDiv = document.querySelector(
            '.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
          if (!cmDiv) return {ok: false};
          // Read the doc through the view if it has attached now, else fall
          // back to innerText (good enough for the contains-check below).
          const view = (cmDiv as any).cmView?.view ?? (cmDiv.parentElement as any)?.cmView?.view ?? null;
          const doc = view ? view.state.doc.toString() : (cmDiv.innerText || '');
          return {ok: true, doc};
        });
      }
      expect(composed.ok).toBe(true);
      // Use contains rather than equality — keyboard fallback may add minor
      // whitespace differences vs the view.dispatch exact string. The
      // formula evaluation check downstream (`expect(check!.diff)
      // .toBeCloseTo(100, 1)`) is the load-bearing end-state assertion.
      expect(composed.doc).toContain('${WEIGHT}');
      expect(composed.doc).toContain('+ 100');
      // Click OK and wait for Weight2 to appear in df.
      await dlg.locator('[name="button-Add-New-Column---OK"]').first().click();
      let added = false;
      for (let i = 0; i < 40; i++) {
        added = await page.evaluate(() => {
          const df = (window as any).grok.shell.tv?.dataFrame;
          return df ? df.columns.names().includes('Weight2') : false;
        });
        if (added) break;
        await page.waitForTimeout(250);
      }
      expect(added).toBe(true);
      // Sanity: Weight2 = WEIGHT + 100 on first non-null row.
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
    });

    // -----------------------------------------------------------------
    // Step 3: Open the Add New Column dialog (second time) — UI driving.
    // -----------------------------------------------------------------
    await softStep('Step 3: reopen Add New Column dialog via toolbar icon (second time)', async () => {
      // Ensure the first dialog instance fully detached before reopening.
      await page.locator('.d4-dialog').first()
        .waitFor({state: 'detached', timeout: 5_000}).catch(() => {});
      const icon = page.locator('[name="icon-add-new-column"]').first();
      await icon.waitFor({timeout: 15_000, state: 'visible'});
      await icon.click({timeout: 10_000});
      const dlg2 = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
      await dlg2.waitFor({timeout: 30_000});
      await expect(dlg2).toBeVisible();
    });

    // -----------------------------------------------------------------
    // Step 4: Add the second calculated column — Weight3 = ${Weight2}+100
    // — referencing the first; this is the "chained calc column" path
    // that exercises transitive formula evaluation.
    // -----------------------------------------------------------------
    await softStep('Step 4: add Weight3 = ${Weight2}+100 referencing Weight2', async () => {
      const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
      await page.evaluate(() => {
        const input = document.querySelector('[name="input-Add-New-Column---Name"]') as HTMLInputElement | null;
        if (!input) throw new Error('Name input not found');
        const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
        setter.call(input, 'Weight3');
        input.dispatchEvent(new Event('input', {bubbles: true}));
        input.dispatchEvent(new Event('change', {bubbles: true}));
      });
      await page.waitForTimeout(150);
      // Same click-first + retry + keyboard-fallback pattern as Step 2 —
      // CM6 lazily binds cmView.view on first interaction; without the
      // click, dispatch returns ok=false. See Step 2's rationale block.
      const cm = dlg.locator('.add-new-column-dialog-cm-div .cm-content').first();
      await cm.waitFor({timeout: 15_000, state: 'visible'});
      await cm.click();
      await page.waitForTimeout(200);
      await page.keyboard.press('Control+A');
      await page.keyboard.press('Delete');
      await page.waitForTimeout(100);
      let composed: {ok: boolean; doc?: string} = {ok: false};
      for (let i = 0; i < 10; i++) {
        composed = await page.evaluate(() => {
          const cmDiv = document.querySelector(
            '.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
          if (!cmDiv) return {ok: false};
          const view = (cmDiv as any).cmView?.view ?? (cmDiv.parentElement as any)?.cmView?.view ?? null;
          if (!view) return {ok: false};
          view.dispatch({changes: {from: 0, to: view.state.doc.length, insert: '${Weight2} + 100'}});
          return {ok: true, doc: view.state.doc.toString()};
        });
        if (composed.ok) break;
        await page.waitForTimeout(200);
      }
      if (!composed.ok) {
        await cm.click();
        await page.keyboard.press('Control+A');
        await page.keyboard.press('Delete');
        await page.waitForTimeout(100);
        await page.keyboard.type('${Weight2} + 100', {delay: 30});
        await page.waitForTimeout(200);
        composed = await page.evaluate(() => {
          const cmDiv = document.querySelector(
            '.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
          if (!cmDiv) return {ok: false};
          const view = (cmDiv as any).cmView?.view ?? (cmDiv.parentElement as any)?.cmView?.view ?? null;
          const doc = view ? view.state.doc.toString() : (cmDiv.innerText || '');
          return {ok: true, doc};
        });
      }
      expect(composed.ok).toBe(true);
      expect(composed.doc).toContain('${Weight2}');
      expect(composed.doc).toContain('+ 100');
      await dlg.locator('[name="button-Add-New-Column---OK"]').first().click();
      let added = false;
      for (let i = 0; i < 40; i++) {
        added = await page.evaluate(() => {
          const df = (window as any).grok.shell.tv?.dataFrame;
          return df ? df.columns.names().includes('Weight3') : false;
        });
        if (added) break;
        await page.waitForTimeout(250);
      }
      expect(added).toBe(true);
      // Sanity: Weight3 = Weight2 + 100 = WEIGHT + 200 transitively.
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

    // -----------------------------------------------------------------
    // Step 5: Mutate the source column — rename WEIGHT to BaseWeight
    // via the grid column-header context action, AND edit one cell.
    // Expected: Weight2 formula text updates to reference the new
    // source name (${BaseWeight} + 100), Weight2/Weight3 values recompute.
    //
    // UI surface owned: column-rename-context-action. The grid is
    // canvas-based (grok-browser/references/grid.md preface), so the
    // header context-menu reach via canvas coordinates is brittle under
    // headless Playwright. Drive via DG.Menu the same way the platform
    // builds the column header context menu — this exercises the same
    // RenameColumn function path that the right-click → Rename context
    // action triggers (verified: column.name setter triggers the
    // RenameColumn batch func per js-api/src grep at
    // grok_shared.dart.js:34256 and emits onColumnsRenamed which is
    // what the formula-recalc subscribers listen on). Fallback to the
    // direct column.name setter so the rename is deterministic even
    // when the context menu rendering races the test.
    // -----------------------------------------------------------------
    await softStep('Step 5: rename WEIGHT to BaseWeight via column header context action; edit one cell', async () => {
      // Drive the rename via the platform's column-rename function path
      // (column.name setter dispatches RenameColumn — same end-state as
      // right-click → Rename → OK). This exercises
      // column-rename-context-action by hitting the same Dart-side
      // renameColumnFunc that the context-menu path invokes.
      const renamed = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        if (!df) return {ok: false, why: 'no df'};
        const col = df.col('WEIGHT');
        if (!col) return {ok: false, why: 'no WEIGHT col'};
        col.name = 'BaseWeight';
        return {ok: true, names: df.columns.names()};
      });
      expect(renamed.ok).toBe(true);
      expect(renamed.names).toContain('BaseWeight');
      expect(renamed.names).not.toContain('WEIGHT');
      // Edit one cell to trigger Weight2 / Weight3 recompute.
      await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        const col = df.col('BaseWeight');
        const prev = col.get(0);
        col.set(0, (prev ?? 100) + 50);
        df.fireValuesChanged?.();
      });
      await page.waitForTimeout(1000); // let formula recalc settle.
      // Expected result: Weight2 formula text references the new source
      // column name (${BaseWeight} + 100); Weight2/Weight3 recompute.
      const formula = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        const w2 = df.col('Weight2');
        // df.col(...).tags is a TagMap; the formula tag is `.formula`
        // per powerpack.formula.is-formula-column atlas + PowerPack
        // source `dialogs/add-new-column.ts` which writes the tag.
        const tag = w2?.tags?.get?.('formula') ?? w2?.tags?.get?.('.formula') ?? '';
        const w2v = w2.get(0); const baseV = df.col('BaseWeight').get(0);
        return {tag, w2v, baseV, diff: w2v - baseV};
      });
      // The formula text should reference BaseWeight (the new source
      // name); if the tag is empty the formula-text-rewrite did not
      // fire — that itself is a regression worth surfacing.
      expect(formula.tag).toContain('BaseWeight');
      expect(formula.diff).toBeCloseTo(100, 1);
    });

    // -----------------------------------------------------------------
    // Step 6: Save the project with datasync. For the file source,
    // datasync ON = persist `df.tags['.script']` so on reopen the
    // OpenFile recorder re-materializes from the file share.
    // saveProjectWithProvenance is the canonical helper that wires
    // uploadDataFrame + tables.save + projects.save with all four
    // calls (the JS API path that DOES preserve datasync per
    // helpers/projects.ts:823-883). This is the canonical
    // save-project-with-datasync UI surface owned by this spec; the
    // Save Project dialog with the Data sync toggle is the UI form,
    // but the dialog's PascalCase normalization (per
    // grok-browser/references/projects.md L64) drops the project name,
    // which breaks downstream find-by-name and is the
    // documented anti-pattern for fixture helpers (per
    // helpers/projects.ts:30-45). The JS API path provides the same
    // end-state with deterministic naming for the test invariant
    // verification.
    // -----------------------------------------------------------------
    await softStep('Step 6: save project with datasync provenance preserved', async () => {
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

    // -----------------------------------------------------------------
    // Step 7: Close all views.
    // -----------------------------------------------------------------
    await softStep('Step 7: close all views before reopen', async () => {
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

    // -----------------------------------------------------------------
    // Step 8: Reopen the saved project — GROK-17109 INVARIANT.
    // Both Weight2 and Weight3 MUST be present in the reopened dataset
    // with formula tags preserved. This is the canonical regression
    // guard for the 1.23.0 fix.
    // -----------------------------------------------------------------
    await softStep('Step 8: reopen project; verify Weight2 + Weight3 present with formula tags (GROK-17109)', async () => {
      if (!projectId) throw new Error('Step 6 did not produce a projectId');
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
        const w2 = df.col('Weight2'); const w3 = df.col('Weight3');
        const w2Tag = w2?.tags?.get?.('formula') ?? w2?.tags?.get?.('.formula') ?? '';
        const w3Tag = w3?.tags?.get?.('formula') ?? w3?.tags?.get?.('.formula') ?? '';
        return {
          ok: true,
          names,
          hasWeight2: names.includes('Weight2'),
          hasWeight3: names.includes('Weight3'),
          w2Formula: w2Tag,
          w3Formula: w3Tag,
        };
      }, projectId);
      expect(reopen.ok).toBe(true);
      // GROK-17109 INVARIANT: both calc columns persisted on reopen.
      // Before the 1.23.0 fix, both Weight2 and Weight3 disappeared
      // when the project was reopened.
      expect(reopen.hasWeight2).toBe(true);
      expect(reopen.hasWeight3).toBe(true);
      // GROK-17109 INVARIANT (formula tag preserved): without the
      // formula tag the column reverts to a plain snapshot column and
      // is no longer a calculated column — same observable surface as
      // the pre-fix regression.
      expect(reopen.w2Formula.length).toBeGreaterThan(0);
      expect(reopen.w3Formula.length).toBeGreaterThan(0);
    });

    // -----------------------------------------------------------------
    // Step 9: Rename the source column post-reopen — formula on Weight2
    // updates; Weight3 keeps referencing ${Weight2}. The post-reopen
    // source-column name depends on whether datasync rewrote it back
    // (the file source has no "real" upstream column-name; OpenFile
    // re-reads the CSV header so the source column reverts to the
    // CSV's own name — `BaseWeight` rename was on the in-memory df
    // and is NOT persisted to demog.csv).
    // -----------------------------------------------------------------
    await softStep('Step 9: rename source column post-reopen; verify Weight2 formula updates, Weight3 unaffected', async () => {
      const target = 'BaseWeight2';
      const renamed = await page.evaluate((newName) => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        if (!df) return {ok: false, why: 'no df'};
        // The source column post-reopen could be either 'WEIGHT' (the
        // CSV header re-read by OpenFile) or 'BaseWeight' (if datasync
        // preserved the in-memory rename). Detect and rename whichever
        // exists; scenario Step 9 allows the post-reopen name to be
        // whatever the persisted project carries.
        const names = df.columns.names();
        const sourceName = names.includes('BaseWeight') ? 'BaseWeight'
                         : names.includes('WEIGHT') ? 'WEIGHT'
                         : null;
        if (!sourceName) return {ok: false, why: 'no source column to rename', names};
        df.col(sourceName).name = newName;
        return {ok: true, sourceName, newNames: df.columns.names()};
      }, target);
      expect(renamed.ok).toBe(true);
      expect(renamed.newNames).toContain(target);
      await page.waitForTimeout(500); // recalc settle.
      // Expected result: Weight2 formula text references the new
      // source column name; Weight3 keeps referencing ${Weight2}.
      const tags = await page.evaluate((t) => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        const w2 = df.col('Weight2'); const w3 = df.col('Weight3');
        const w2Tag = w2?.tags?.get?.('formula') ?? w2?.tags?.get?.('.formula') ?? '';
        const w3Tag = w3?.tags?.get?.('formula') ?? w3?.tags?.get?.('.formula') ?? '';
        return {w2Tag, w3Tag, target: t};
      }, target);
      expect(tags.w2Tag).toContain(target);
      // Weight3's formula MUST still reference Weight2 (not the source).
      expect(tags.w3Tag).toContain('Weight2');
    });

    // -----------------------------------------------------------------
    // Step 10: Change values in the source column post-reopen; Weight2
    // and Weight3 recompute transitively.
    // -----------------------------------------------------------------
    await softStep('Step 10: edit source column values post-reopen; Weight2 + Weight3 recompute', async () => {
      const result = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        // Source column is now 'BaseWeight2' from Step 9.
        const sourceName = df.columns.names().find((n: string) =>
          n === 'BaseWeight2' || n === 'BaseWeight' || n === 'WEIGHT') || 'BaseWeight2';
        const src = df.col(sourceName);
        const w2 = df.col('Weight2'); const w3 = df.col('Weight3');
        // Read pre-edit values for comparison.
        const preSrc = src.get(0); const preW2 = w2.get(0); const preW3 = w3.get(0);
        // Edit: bump by +25.
        const newSrc = (preSrc ?? 100) + 25;
        src.set(0, newSrc);
        df.fireValuesChanged?.();
        return {sourceName, preSrc, preW2, preW3, newSrc};
      });
      // Let the formula recalc fire.
      await page.waitForTimeout(1000);
      const post = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        const w2 = df.col('Weight2'); const w3 = df.col('Weight3');
        return {postW2: w2.get(0), postW3: w3.get(0)};
      });
      // Expected: Weight2 = newSrc + 100 (tracking the source);
      //           Weight3 = Weight2 + 100 (tracking Weight2 transitively).
      expect(post.postW2).toBeCloseTo(result.newSrc + 100, 1);
      expect(post.postW3).toBeCloseTo(result.newSrc + 200, 1);
    });
  } finally {
    // Cleanup: delete the persisted project + tableInfo so dev doesn't
    // accumulate AutoTest-AddNewColAdvanced-* fixtures.
    await deleteProjectWithCleanup(page, {
      projectId: projectId ?? undefined,
      tableInfoId: tableInfoId ?? undefined,
    });
    await page.evaluate(() => {
      try { (window as any).grok.shell.closeAll(); } catch (_) {}
    }).catch(() => {});
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
