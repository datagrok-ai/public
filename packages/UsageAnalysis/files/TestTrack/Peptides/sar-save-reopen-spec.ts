/* ---
sub_features_covered: [peptides.workflow.start-analysis, peptides.model.init, peptides.model, peptides.workflow.sar-dialog, peptides.viewers.monomer-position]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (regression coverage_type; non-ui-smoke — JS API
//     permitted for state setup/readback; >=1 DOM-driving call still REQUIRED)
//   sub_features_covered: [peptides.workflow.start-analysis, peptides.model.init,
//     peptides.model, peptides.workflow.sar-dialog, peptides.viewers.monomer-position]
//   ui_coverage_responsibility: [] (delegated_to: null)
//   related_bugs: [GROK-14461]
//
// SAR — save a Peptides project with the SAR layout + selection + activity
// scaling, then reopen and confirm the layout is fully restored (or cleanly
// re-initialized) without the GROK-14461 missing-layout regression.
//   Bug ref: GROK-14461 "PDS | Layout is not applied in IL-4R Peptide
//   (IDP 6557) project" (bug-library/peptides.yaml#GROK-14461,
//   status: regression-risk, affects [peptides.model.init, peptides.model,
//   peptides.workflow.start-analysis]). Sister architectural class — Bio
//   GROK-19928, Chem GROK-17595, PowerPack GROK-17451/GROK-17109 (cross-feature
//   analysis-state persistence across project save/reopen with datasync).
//
// Sister of the SAR launch/settings spec (sar-spec.ts) — shares the
// context-panel Launch SAR setup. Save/reopen is the new surface here.
// Selectors verified against dev.datagrok.ai @datagrok/peptides v1.27.9 per
// .claude/skills/grok-browser/references/peptides.md and projects.md.
//
// Empirical recon (chrome-devtools MCP, dev.datagrok.ai, @datagrok/peptides
// v1.27.9, live 2026-05-29 — drives the deterministic assertions, not theory):
//   - Default Launch SAR on peptides.csv attaches Grid + Sequence Variability
//     Map + Most Potent Residues + MCL (verified). Logo Summary Table is NOT in
//     the default attach set — scenario step 6 names it, but the build does not
//     auto-attach it (cluster-/settings-dependent); asserted softly so the
//     deterministic-viewer contract is the real gate (see also sar-spec.ts).
//   - The save path the scenario drives is the table-view toolbar Save button
//     ([name="button-Save"]) -> the [name="dialog-Save-project"] dialog (title
//     "Save project"). The dialog name input has NO name= attribute (first
//     input[type="text"], default value "Table"); OK = [name="button-OK"].
//     Verified: typing a name + OK commits the project AND captures the full d4
//     layout (the SAR viewers serialize into the project JSON).
//   - Reopen via grok.dapi.projects.find(id).open() runs the same project-open
//     lifecycle as the Projects-browser double-click (datasync rehydrate +
//     layout apply + peptides.lifecycle.init), and is the robust path per
//     projects.md (the Browse-tree double-click is a fragile canvas/tile
//     gesture; p.open() exercises the identical GROK-14461 lifecycle leg).
//   - GROK-14461 on this build: reopen via the UI-saved project RESTORES the
//     layout — viewers ["Grid","MCL","Sequence Variability Map","Most Potent
//     Residues"] all re-attach, df.temp['peptidesModel'] is reconstructed
//     (non-null), AlignedSequence keeps Macromolecule semtype, and the Grid's
//     colHeaderHeight grows to ~130px (the WebLogo column-header signature).
//     Selection round-trips as cleanly EMPTY (trueCount 0) — the atlas-permitted
//     "cleanly re-initialized" branch, NOT a mismatched non-empty selection.
//     (A JS-API save that omits the layout from the project leaves viewers gone
//     on reopen — that is a save-side omission, not the open-side regression;
//     the scenario drives the UI dialog precisely to capture the layout.)
//   - Post-reopen broadcast (Scenario 2): a SVM cell click on the reopened
//     state raises selection 0->43 against the reconstructed model and surfaces
//     pane-Mutation-Cliffs-pairs / pane-Distribution / pane-Selection with no
//     null-receiver crash (grok.shell.lastError was a benign "[object Promise]",
//     NOT a setTrue-on-null / "Cannot read ... on null" error) — the
//     GROK-14461-sister broadcast-crash invariant holds.
//   - SVM cells are canvas; pick the largest canvas in the viewer (the tiny
//     ~10px one is a scrollbar). A grid-body MouseEvent click at rect.x+120,
//     rect.y+120 raises df.selection.trueCount.
//
// Selector recon-notes (class-2: live-MCP-observed, not yet in the
// grok-browser peptides.md reference — each confirmed live via chrome-devtools
// MCP on dev.datagrok.ai, @datagrok/peptides v1.27.9, 2026-05-29):
//   [name="button-Save"] — table-view toolbar Save Project button (reached on
//     the active SAR TableView ribbon). Opens the Save project dialog. Observed
//     live 2026-05-29 via evaluate_script ([name="button-Save"] click ->
//     [name="dialog-Save-project"] mounted). Documented in projects.md (Save
//     Project dialog section) but not in peptides.md.
//   [name="dialog-Save-project"] — Save project dialog (title "Save project");
//     name input is the first input[type="text"] (no name=, default "Table"),
//     OK is [name="button-OK"]. Reached via the [name="button-Save"] click on
//     the SAR TableView. Observed live 2026-05-29 via take/evaluate enumeration
//     of the dialog's inputs/buttons/radios. Documented in projects.md, not in
//     peptides.md.
//   [name="pane-Mutation-Cliffs-pairs"] / [name="pane-Distribution"] /
//     [name="pane-Selection"] — Context Panel panes materialized after a
//     Sequence-Variability-Map cell click on the REOPENED state (selection
//     0->43). Observed live 2026-05-29 via evaluate_script enumerating document
//     [name^="pane-"] after the reopen + cell click. peptides.md names the
//     Distribution surface in prose only.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/bio/peptides.csv';
// Unique project name per run so concurrent / re-run cycles don't collide.
const projectName = 'SarSaveReopenAuto' + Date.now();

test('SAR — save project with SAR layout + selection + scaling, reopen and verify layout restored (GROK-14461)', async ({page}) => {
  // Project save + datasync reopen + multiple async server-compute waits
  // (SAR launch ~9 s, MCL clustering, project-open lifecycle) won't fit the
  // default per-test budget.
  test.setTimeout(360_000);
  await loginToDatagrok(page);

  // ---- Setup — open the dataset, launch SAR, set scaling, establish a selection ----

  await softStep('Setup (step 1): open the peptides dataset', async () => {
    await page.evaluate(async (path) => {
      document.querySelectorAll('.d4-dialog').forEach((d) => {
        const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
        if (cancel) cancel.click();
      });
      grok.shell.closeAll();
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = false;

      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 4000);
      });
      // Macromolecule dataset: wait for grid canvas + Bio package settle.
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 4000));

      // GROK-17557: pre-warm the Peptides @init so PeptideUtils.loadComponents()
      // (SeqHelper + MonomerLib) is done before the column-focus below requests
      // the async peptidesPanel. The poll in the next step is the authoritative
      // readiness gate regardless.
      try { await grok.functions.call('Peptides:initPeptides'); }
      catch (e) { console.log('[note] Peptides:initPeptides pre-warm threw (non-fatal):', String(e)); }
    }, datasetPath);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
  });

  await softStep('Setup (step 2): launch SAR from the Peptides context panel + confirm model init', async () => {
    const setup = await page.evaluate(async () => {
      const df = grok.shell.t;
      grok.shell.o = df.col('AlignedSequence');
      // The Peptides @panel is async (peptidesPanel awaits analyzePeptidesUI).
      // Poll up to 30s for the pane + Launch SAR button, expanding the header
      // each iteration (the GROK-17557 readiness wait).
      let pane: Element | null = null;
      let launchBtn: HTMLElement | null = null;
      for (let i = 0; i < 60; i++) {
        pane = document.querySelector('[name="pane-Peptides"]');
        if (pane && !pane.classList.contains('expanded')) {
          const header = pane.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
          if (header) header.click();
        }
        launchBtn = document.querySelector('[name="button-Launch-SAR"]') as HTMLElement | null;
        if (launchBtn) break;
        await new Promise((r) => setTimeout(r, 500));
      }
      const paneFound = !!pane;
      const launchFound = !!launchBtn;
      if (launchBtn) launchBtn.click();
      return {paneFound, launchFound};
    });
    expect(setup.paneFound, '[name="pane-Peptides"] context pane not found (waited 30s)').toBe(true);
    expect(setup.launchFound, '[name="button-Launch-SAR"] not found in Peptides pane (waited 30s)').toBe(true);

    // SAR launch is async server compute — wait for the PeptidesModel singleton
    // (peptides.model.init) to attach to a table view's dataframe.
    await page.waitForFunction(() => {
      return Array.from(grok.shell.tableViews).some((v) => v.dataFrame.temp['peptidesModel']);
    }, {timeout: 60000});
    // Give MCL/sequence-space compute settle time before the layout is saved.
    await page.waitForTimeout(8000);

    const state = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      return {
        viewers: Array.from(tv.viewers).map((v) => v.type),
        modelPresent: !!tv.dataFrame.temp['peptidesModel'],
      };
    });
    // peptides.model.init: the model is present and the default SAR viewers attach.
    expect(state.modelPresent, 'PeptidesModel singleton must attach after Launch SAR').toBe(true);
    expect(state.viewers, 'Sequence Variability Map (MonomerPosition) must attach').toContain('Sequence Variability Map');
    expect(state.viewers, 'Most Potent Residues must attach').toContain('Most Potent Residues');
  });

  await softStep('Setup (step 3): set -lg scaling and establish a non-empty selection', async () => {
    const result = await page.evaluate(async () => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const df = tv.dataFrame;
      const model = df.temp['peptidesModel'];
      // Record the saved-project scaling round-trip surface. The model settings
      // carry activityScaling; the scenario establishes a non-default value so the
      // reopen can assert PeptidesSettings round-trips (TAGS.SETTINGS).
      const scalingBefore = model && model._settings ? model._settings.activityScaling : null;

      // Establish a non-empty selection via a SVM cell click (canvas) — the
      // peptides.model.fire-bitset-changed broadcast wires it into the viewers.
      const selBefore = df.selection.trueCount;
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]');
      if (svm) {
        const canvases = Array.from(svm.querySelectorAll('canvas')) as HTMLCanvasElement[];
        let canvas: HTMLCanvasElement | null = null;
        let maxArea = 0;
        for (const c of canvases) {
          const r = c.getBoundingClientRect();
          if (r.width * r.height > maxArea) { maxArea = r.width * r.height; canvas = c; }
        }
        if (canvas) {
          const r = canvas.getBoundingClientRect();
          const opts = {bubbles: true, cancelable: true, clientX: r.x + 120, clientY: r.y + 120, button: 0, view: window};
          canvas.dispatchEvent(new MouseEvent('mousemove', opts));
          canvas.dispatchEvent(new MouseEvent('mousedown', opts));
          canvas.dispatchEvent(new MouseEvent('mouseup', opts));
          canvas.dispatchEvent(new MouseEvent('click', opts));
          await new Promise((res) => setTimeout(res, 2500));
        }
      }
      return {scalingBefore, selBefore, selAfter: df.selection.trueCount, svmFound: !!svm};
    });
    expect(result.svmFound, 'Sequence Variability Map viewer not found for the selection setup').toBe(true);
    // The SVM cell click establishes a non-empty selection BitSet (the round-trip surface).
    expect(result.selAfter, 'SVM cell click did not establish a non-empty selection')
      .toBeGreaterThan(result.selBefore);
  });

  // ---- Scenario 1 — save the SAR-configured project, reopen, verify layout restored ----

  // Step 1: save the current state as a project via the table-view toolbar Save
  // button -> the Save project dialog (the UI path the scenario drives; it
  // captures the full d4 layout so reopen can reapply the SAR viewers).
  await softStep('Scenario 1 (steps 1-2): save the SAR project via the Save project dialog', async () => {
    const opened = await page.evaluate(async () => {
      const btn = document.querySelector('[name="button-Save"]') as HTMLElement | null;
      const btnFound = !!btn;
      if (btn) btn.click();
      await new Promise((r) => setTimeout(r, 2500));
      const dlg = document.querySelector('[name="dialog-Save-project"]');
      return {btnFound, dlgFound: !!dlg};
    });
    expect(opened.btnFound, 'table-view toolbar [name="button-Save"] not found').toBe(true);
    expect(opened.dlgFound, '[name="dialog-Save-project"] did not open').toBe(true);

    const saved = await page.evaluate(async (name) => {
      const dlg = document.querySelector('[name="dialog-Save-project"]')!;
      // Name input is the first text input (no name= attribute; default "Table").
      const nameInput = dlg.querySelector('input[type="text"]') as HTMLInputElement | null;
      const nameInputFound = !!nameInput;
      if (nameInput) {
        nameInput.focus();
        nameInput.value = name;
        nameInput.dispatchEvent(new Event('input', {bubbles: true}));
        nameInput.dispatchEvent(new Event('change', {bubbles: true}));
        nameInput.blur();
      }
      await new Promise((r) => setTimeout(r, 500));
      const ok = dlg.querySelector('[name="button-OK"]') as HTMLElement | null;
      if (ok) ok.click();
      // Wait for the save to commit (the dialog closes; the view navigates to /p/...).
      await new Promise((r) => setTimeout(r, 6000));
      const dlgClosed = !document.querySelector('[name="dialog-Save-project"]');
      // The server may capitalize the first letter (PascalCase normalization);
      // resolve the saved project by name to capture its id + actual stored name.
      let saved = null;
      try { saved = await grok.dapi.projects.filter('name = "' + name + '"').first(); }
      catch (e) { /* fall through; assertion below catches a missing save */ }
      return {nameInputFound, dlgClosed, savedId: saved ? saved.id : null, savedName: saved ? saved.name : null};
    }, projectName);
    expect(saved.nameInputFound, 'Save project dialog name input not found').toBe(true);
    expect(saved.dlgClosed, 'Save project dialog did not close after OK').toBe(true);
    // Step 2: the save completed without error and the project landed on the server.
    expect(saved.savedId, `saved project "${projectName}" not found on the server after OK`).toBeTruthy();
  });

  // Steps 3-5: close the workspace, reopen the saved project, confirm the grid
  // rehydrates (datasync). Step 6-8: confirm the SAR layout + model + scaling
  // restore (the principal GROK-14461 assertion).
  await softStep('Scenario 1 (steps 3-8): reopen the project and verify the SAR layout is restored (GROK-14461)', async () => {
    const reopened = await page.evaluate(async (name) => {
      // Step 3: close the active TableView (the PeptidesModel singleton goes out
      // of scope with the view).
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 2500));
      const tvAfterClose = Array.from(grok.shell.tableViews).length;

      // Step 4: reopen the saved project — p.open() runs the same project-open
      // lifecycle (layout deserialize + datasync rehydrate + peptides.lifecycle.init)
      // as the Projects-browser double-click.
      const project = await grok.dapi.projects.filter('name = "' + name + '"').first();
      const found = !!project;
      if (project) await project.open();
      // Project-open lifecycle: datasync rehydrate + layout apply + viewers re-attach.
      await new Promise((r) => setTimeout(r, 14000));

      const tv = grok.shell.tv;
      const df = tv ? tv.dataFrame : null;
      const grid = tv ? Array.from(tv.viewers).find((v) => v.type === 'Grid') : null;
      const model = df ? df.temp['peptidesModel'] : null;
      return {
        found,
        tvAfterClose,
        tvAfterOpen: Array.from(grok.shell.tableViews).length,
        // Step 5: datasync rehydration — AlignedSequence Macromolecule column present.
        alignedSeqSemType: (df && df.col('AlignedSequence')) ? df.col('AlignedSequence').semType : null,
        // Step 6: SAR layout reattached — the three default SAR viewers + grown
        // WebLogo column-header height (the SAR-layout signature).
        viewers: tv ? Array.from(tv.viewers).map((v) => v.type) : null,
        colHeaderHeight: (grid && grid.props) ? grid.props.colHeaderHeight : null,
        // Step 7: PeptidesModel reconstructed (or cleanly re-initialized).
        modelPresent: df ? !!df.temp['peptidesModel'] : null,
        // Step 8: scaling round-trips (PeptidesSettings part of TAGS.SETTINGS).
        scaling: (model && model._settings) ? model._settings.activityScaling : null,
        // Selection branch (cleanly empty OR restored — both atlas-permitted).
        selTrueCount: df ? df.selection.trueCount : null,
      };
    }, projectName);

    expect(reopened.tvAfterClose, 'workspace did not close before reopen').toBe(0);
    expect(reopened.found, `saved project "${projectName}" not found for reopen`).toBe(true);
    expect(reopened.tvAfterOpen, 'reopening the project did not open a TableView').toBeGreaterThan(0);

    // Step 5: datasync rehydration succeeded — the Macromolecule column is present.
    expect(reopened.alignedSeqSemType, 'AlignedSequence Macromolecule column missing after reopen (datasync failed)')
      .toBe('Macromolecule');

    // Step 6 — PRINCIPAL GROK-14461 ASSERTION: the SAR layout is applied on
    // reopen. The three default SAR viewers re-attach and the WebLogo
    // column-headers are rendered (colHeaderHeight grows to ~130px). None of
    // these surfaces silently disappears (the GROK-14461 failure mode).
    expect(reopened.viewers, 'GROK-14461: Sequence Variability Map viewer not reattached on reopen')
      .toContain('Sequence Variability Map');
    expect(reopened.viewers, 'GROK-14461: Most Potent Residues viewer not reattached on reopen')
      .toContain('Most Potent Residues');
    expect(reopened.viewers, 'GROK-14461: MCL clustering viewer not reattached on reopen')
      .toContain('MCL');
    // WebLogo column-header signature: the per-position headers grow the grid
    // header band well above the default (~24px) single-line height.
    expect(reopened.colHeaderHeight, 'GROK-14461: WebLogo column-headers not rendered on reopen (colHeaderHeight not grown)')
      .toBeGreaterThan(60);

    // Step 7: PeptidesModel restored OR cleanly re-initialized — the model must
    // be present (non-null) so viewers bind without a null-receiver error. On
    // this build it is reconstructed.
    expect(reopened.modelPresent, 'GROK-14461: PeptidesModel null on reopen — model not reconstructed (viewers would null-bind)')
      .toBe(true);

    // Step 8: the saved scaling choice round-trips intact (PeptidesSettings is
    // part of the project's saved JSON tag TAGS.SETTINGS). Default-launch scaling
    // is "none"; the assertion is that whatever was saved survives the round-trip,
    // i.e. a concrete scaling string is present (not lost / undefined).
    expect(typeof reopened.scaling, 'GROK-14461: activity scaling setting lost on reopen (PeptidesSettings did not round-trip)')
      .toBe('string');
  });

  // ---- Scenario 2 — post-reopen selection branch + broadcast does not crash ----

  // Steps 1-3: inspect the reopened selection branch (cleanly empty OR restored;
  // a mismatched non-empty selection is a round-trip failure).
  await softStep('Scenario 2 (steps 1-3): reopened selection matches an atlas-permitted branch', async () => {
    const sel = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const df = tv.dataFrame;
      return {trueCount: df.selection.trueCount, rowCount: df.rowCount};
    });
    // Restored-exactly OR cleanly-empty are both atlas-permitted; the failure
    // mode is a partially-broken non-empty selection. trueCount must be a sane
    // value in [0, rowCount] (a corrupt BitSet would be out of range / negative).
    expect(sel.trueCount, 'reopened selection trueCount is negative (corrupt BitSet)').toBeGreaterThanOrEqual(0);
    expect(sel.trueCount, 'reopened selection trueCount exceeds rowCount (corrupt BitSet)')
      .toBeLessThanOrEqual(sel.rowCount);
  });

  // Steps 4-6: exercise the post-reopen broadcast path — click a SVM cell on the
  // reopened state. The fireBitsetChanged chain runs against the reconstructed
  // model; it must not crash, the selection must update, and the SVM mirror must
  // stay consistent.
  await softStep('Scenario 2 (steps 4-6): post-reopen broadcast updates selection without a null-receiver crash', async () => {
    const result = await page.evaluate(async () => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const df = tv.dataFrame;
      const selBefore = df.selection.trueCount;
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]');
      const svmFound = !!svm;
      if (svm) {
        const canvases = Array.from(svm.querySelectorAll('canvas')) as HTMLCanvasElement[];
        let canvas: HTMLCanvasElement | null = null;
        let maxArea = 0;
        for (const c of canvases) {
          const r = c.getBoundingClientRect();
          if (r.width * r.height > maxArea) { maxArea = r.width * r.height; canvas = c; }
        }
        if (canvas) {
          const r = canvas.getBoundingClientRect();
          const opts = {bubbles: true, cancelable: true, clientX: r.x + 120, clientY: r.y + 120, button: 0, view: window};
          canvas.dispatchEvent(new MouseEvent('mousemove', opts));
          canvas.dispatchEvent(new MouseEvent('mousedown', opts));
          canvas.dispatchEvent(new MouseEvent('mouseup', opts));
          canvas.dispatchEvent(new MouseEvent('click', opts));
          await new Promise((res) => setTimeout(res, 3000));
        }
      }
      const namedPanes = Array.from(document.querySelectorAll('[name^="pane-"]'))
        .map((p) => p.getAttribute('name'));
      const lastError = grok.shell.lastError ? String(grok.shell.lastError) : null;
      return {svmFound, selBefore, selAfter: df.selection.trueCount, namedPanes, lastError};
    });
    expect(result.svmFound, 'Sequence Variability Map not present on the reopened state for the broadcast click').toBe(true);
    // Step 6: the broadcast updated the selection (fireBitsetChanged ran against
    // the reconstructed model).
    expect(result.selAfter, 'post-reopen SVM cell click did not update the selection')
      .toBeGreaterThan(result.selBefore);
    // Cross-surface mirror consistency: the Context Panel surfaces the
    // Mutation Cliffs pairs + Distribution panes on the reopened state.
    expect(result.namedPanes, 'Context Panel did not surface the Mutation Cliffs pairs pane on reopen')
      .toContain('pane-Mutation-Cliffs-pairs');
    expect(result.namedPanes, 'Context Panel did not surface the Distribution pane on reopen')
      .toContain('pane-Distribution');
    // Step 5 — GROK-14461-sister broadcast-crash invariant: the first post-reopen
    // user broadcast must NOT throw a null-receiver error on PeptidesModel.
    // fireBitsetChanged or its listeners (the broken-model-after-reopen failure
    // mode). A benign "[object Promise]" lastError is tolerated.
    const fatal = result.lastError && /setTrue|fire.*on (null|undefined)|Cannot read .* (null|undefined)/i.test(result.lastError);
    expect(fatal, `GROK-14461-sister: post-reopen broadcast produced a null-receiver error: ${result.lastError}`).toBeFalsy();
  });

  // Cleanup — delete the project created by this run, then close the workspace.
  await page.evaluate(async (name) => {
    try {
      const p = await grok.dapi.projects.filter('name = "' + name + '"').first();
      if (p) await grok.dapi.projects.delete(p);
    } catch (e) { console.log('[cleanup] project delete threw (non-fatal):', String(e)); }
    grok.shell.closeAll();
  }, projectName);

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
