/* ---
sub_features_covered: [peptides.viewers.sar-base.export-mutation-cliffs, peptides.viewers.sar-base, peptides.viewers.monomer-position, peptides.workflow.start-analysis]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (coverage_type: smoke; non-ui-smoke — JS API substitution permitted)
//   sub_features_covered: [peptides.viewers.sar-base.export-mutation-cliffs,
//     peptides.viewers.sar-base, peptides.viewers.monomer-position,
//     peptides.workflow.start-analysis]
//   ui_coverage_responsibility: [] (delegated_to: null)
//   related_bugs: []
//
// Atlas provenance (derived_from):
//   feature-atlas/peptides.yaml#critical_paths[export-mutation-cliffs-from-sar-viewer]
//     derived_from: public/packages/Peptides/src/viewers/sar-viewer.ts#L560
//
// Export Mutation Cliffs from the Sequence Variability Map (SARViewer) context menu.
// Sister of the Export Invariant Map surface (export-invariant-map.md, open authoring item).
// Selectors verified against dev.datagrok.ai @datagrok/peptides v1.27.9 per
// .claude/skills/grok-browser/references/peptides.md.
//
// Empirical recon notes (live MCP, dev.datagrok.ai, 2026-05-29 — drive deterministic
// assertions, not theory):
//   - SAR launch via the context-panel Launch SAR button is the deterministic path
//     (top-menu Bio | Analyze | SAR... falls into menu overflow at the recon viewport).
//     Default attach set on peptides.csv: Grid, Sequence Variability Map,
//     Most Potent Residues, MCL.
//   - The Export actions live on the *Sequence Variability Map viewer* context menu
//     (a contextmenu dispatch on the inner SVM canvas surfaces the SARViewer-
//     contributed Export group merged with the generic grid menu). The
//     `Export Mutation Cliffs...` `.d4-menu-item-label` is present at top level; clicking
//     its `.d4-menu-item` opens `[name="dialog-Export-Mutation-Cliffs"]` (title
//     `Export Mutation Cliffs`, OK/Cancel).
//   - OK with the default (empty) Extra-columns selection opens a NEW TableView named
//     `Mutation Cliffs` (647-row demo → 6253 cliff-pair rows) with exactly:
//       Seq 1, Seq 2, Mutation (semType MacromoleculeDifference, value `seq1#seq2`),
//       Seq 1 IC50, Seq 2 IC50, Delta
//     (activity columns named `Seq 1 <activity>` / `Seq 2 <activity>`; activity col = IC50).
//   - Scenario 2 (extra-column selection): the dialog's `[name="input-Extra-columns"]`
//     is a `ui.input.columns` multiselect whose picker popup renders lazily/canvas-style
//     (`.d4-combo-popup` with a `visibility:hidden` drop-down, no enumerable per-column
//     rows; no `DG.InputBase.forInputRoot` on this build; the column value lives in an
//     unreachable closure inside SARViewer.exportMutationCliffs, NOT in dialog.inputs).
//     The column-multiselect gesture is therefore not DOM-addressable. Per the grok-browser
//     Core Principle step 3 (UI attempt failed -> JS API fallback, record reason) and the
//     non-ui-smoke paradigm (JS API substitution permitted), Scenario 2 drives the same
//     viewer context menu + Export dialog (real UI) for the trigger, then exercises the
//     production paired-column emission via the viewer's own
//     `_doExportMutationCliffs([extraCol])` (the exact code path the dialog OK invokes —
//     verified live: passing `[ID]` adds `Seq 1 ID` / `Seq 2 ID`; passing `[ID, Activity]`
//     adds both paired pairs).

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/bio/peptides.csv';

// Base columns the export always emits (activity col on the demo dataset is IC50).
const BASE_EXPORT_COLS = ['Seq 1', 'Seq 2', 'Mutation', 'Seq 1 IC50', 'Seq 2 IC50', 'Delta'];

test('Peptides — Export Mutation Cliffs to a new TableView', async ({page}) => {
  // SAR launch + MCL/sequence-space compute and the export roundtrips need more than the
  // default per-test budget.
  test.setTimeout(300_000);
  await loginToDatagrok(page);

  // ---- Setup — open the peptides dataset and launch SAR (context-panel button) ----

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
      // (SeqHelper + MonomerLib) is done before the column-focus below requests the
      // async peptidesPanel. The poll in the next step is the authoritative readiness gate.
      try { await grok.functions.call('Peptides:initPeptides'); }
      catch (e) { console.log('[note] Peptides:initPeptides pre-warm threw (non-fatal):', String(e)); }
    }, datasetPath);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
  });

  await softStep('Setup (step 2): launch SAR from the Peptides context panel', async () => {
    const setup = await page.evaluate(async () => {
      const df = grok.shell.t;
      grok.shell.o = df.col('AlignedSequence');
      // The Peptides @panel is async; poll for the pane + Launch SAR button rather than a
      // fixed wait + single query (the GROK-17557 cold-package race).
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
      const result = {paneFound: !!pane, launchFound: !!launchBtn};
      if (launchBtn) launchBtn.click();
      return result;
    });
    expect(setup.paneFound, '[name="pane-Peptides"] context pane not found (waited 30s)').toBe(true);
    expect(setup.launchFound, '[name="button-Launch-SAR"] not found in Peptides pane (waited 30s)').toBe(true);

    // SAR launch is async server compute — wait for the PeptidesModel singleton.
    await page.waitForFunction(() => {
      return Array.from(grok.shell.tableViews).some((v) => v.dataFrame.temp['peptidesModel']);
    }, {timeout: 60000});
    // Give MCL/sequence-space compute + viewer attach extra settle time so the SVM and its
    // mutation-cliffs computation are ready before the export.
    await page.waitForTimeout(9000);

    const svmReady = await page.evaluate(() => {
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]');
      return {svmFound: !!svm, svmHasCanvas: svm ? !!svm.querySelector('canvas') : false};
    });
    expect(svmReady.svmFound, 'Sequence Variability Map viewer did not attach after SAR launch').toBe(true);
    expect(svmReady.svmHasCanvas, 'Sequence Variability Map did not render its canvas').toBe(true);
  });

  // ---- Scenario 1 — Export Mutation Cliffs (default selection) opens a documented TableView ----

  // Steps 1-3: right-click the SVM viewer body, open the Export submenu, click
  // Export Mutation Cliffs, accept defaults in the column-picker dialog, click OK.
  await softStep('Scenario 1 (steps 1-3): open Export Mutation Cliffs dialog and accept defaults', async () => {
    const opened = await page.evaluate(async () => {
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]')!;
      const canvases = Array.from(svm.querySelectorAll('canvas')) as HTMLCanvasElement[];
      // The largest canvas is the SVM render surface (the small ones are scrollbars).
      let canvas: HTMLCanvasElement | null = null;
      let maxArea = 0;
      for (const c of canvases) {
        const r = c.getBoundingClientRect();
        if (r.width * r.height > maxArea) { maxArea = r.width * r.height; canvas = c; }
      }
      if (!canvas) return {canvasFound: false};
      const r = canvas.getBoundingClientRect();
      const cx = r.x + 120, cy = r.y + 120;
      // Step 1: open the viewer context menu on the SVM body.
      canvas.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: cx, clientY: cy, view: window}));
      canvas.dispatchEvent(new MouseEvent('contextmenu',
        {bubbles: true, cancelable: true, button: 2, clientX: cx, clientY: cy, view: window}));
      await new Promise((rs) => setTimeout(rs, 1200));

      // Step 2: hover the Export group, then click Export Mutation Cliffs...
      const allLabels = Array.from(document.querySelectorAll('.d4-menu-item-label'));
      const exportGroup = allLabels.find((e) => e.textContent?.trim() === 'Export');
      if (exportGroup) {
        const g = exportGroup.closest('.d4-menu-item') as HTMLElement | null;
        if (g) {
          g.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, view: window}));
          g.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true, view: window}));
        }
      }
      await new Promise((rs) => setTimeout(rs, 500));
      const exportMc = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find((e) => e.textContent?.trim() === 'Export Mutation Cliffs...');
      const menuItemFound = !!exportMc;
      if (exportMc) {
        const item = exportMc.closest('.d4-menu-item') as HTMLElement | null;
        if (item) item.click();
      }
      await new Promise((rs) => setTimeout(rs, 1500));

      // Step 3: the column-picker dialog opens; accept defaults (no extra columns) -> OK.
      const dlg = document.querySelector('[name="dialog-Export-Mutation-Cliffs"]');
      const dlgFound = !!dlg;
      const tvBefore = Array.from(grok.shell.tableViews).length;
      const ok = dlg ? dlg.querySelector('[name="button-OK"]') as HTMLElement | null : null;
      if (ok) ok.click();
      return {canvasFound: true, menuItemFound, dlgFound, tvBefore};
    });
    expect(opened.canvasFound, 'SVM render canvas not found for the context-menu trigger').toBe(true);
    expect(opened.menuItemFound,
      '`Export Mutation Cliffs...` not found on the SVM viewer context menu').toBe(true);
    expect(opened.dlgFound, '[name="dialog-Export-Mutation-Cliffs"] did not open').toBe(true);
  });

  // Steps 4-6: a new `Mutation Cliffs` TableView opens, becomes active, and carries the
  // documented column shape + a non-empty grid.
  await softStep('Scenario 1 (steps 4-6): exported TableView has the documented column shape', async () => {
    // The export builds a DataFrame + addTableView; allow a brief settle.
    await page.waitForTimeout(3000);
    const result = await page.evaluate((baseCols) => {
      const tvs = Array.from(grok.shell.tableViews);
      const mcView = tvs.find((v) => v.dataFrame && v.dataFrame.name === 'Mutation Cliffs');
      if (!mcView) return {found: false};
      const df = mcView.dataFrame;
      const mut = df.col('Mutation');
      return {
        found: true,
        cols: df.columns.names(),
        rows: df.rowCount,
        mutationSemType: mut ? mut.semType : null,
        mutationSample: mut && df.rowCount > 0 ? String(mut.get(0)) : null,
        activeViewName: grok.shell.v ? grok.shell.v.name : null,
        // Sanity that the base columns are exactly present.
        hasAllBase: baseCols.every((c: string) => df.columns.names().includes(c)),
      };
    }, BASE_EXPORT_COLS);

    // Step 4: a new TableView opened and became the active view.
    expect(result.found, 'new `Mutation Cliffs` TableView was not opened by the export').toBe(true);
    expect(result.activeViewName, 'the exported `Mutation Cliffs` view did not become active')
      .toBe('Mutation Cliffs');
    // Step 5: documented column shape — Seq 1 / Seq 2 / Mutation / activities / Delta.
    expect(result.cols, 'exported grid is missing one or more documented columns').toEqual(BASE_EXPORT_COLS);
    expect(result.hasAllBase, 'exported grid does not carry the full documented base column set').toBe(true);
    // Mutation column carries the MacromoleculeDifference semtype, value `seq1#seq2`.
    expect(result.mutationSemType, 'Mutation column must have semType MacromoleculeDifference')
      .toBe('MacromoleculeDifference');
    expect(result.mutationSample, 'Mutation value must be formatted `seq1#seq2`')
      .toContain('#');
    // Step 6: the grid is non-empty (at least one mutation-cliff pair).
    expect(result.rows, 'exported mutation-cliffs grid is empty').toBeGreaterThan(0);
  });

  // ---- Scenario 2 — Exported table includes paired columns when an extra column is selected ----

  // Step 1: return to the original SAR TableView.
  // Steps 2-3: re-open the Export Mutation Cliffs dialog from the SVM context menu (real UI),
  //   then select one extra non-default column. The column-multiselect picker is not
  //   DOM-addressable on this build (see recon notes header) — per the grok-browser
  //   Core Principle step 3, the extra-column selection falls back to the viewer's own
  //   production export path with the chosen column, exercising the paired-column emission.
  // Step 4: the new TableView carries the base columns plus a Seq-1/Seq-2 paired pair for
  //   the extra column.
  await softStep('Scenario 2 (steps 1-4): extra column selection emits paired columns', async () => {
    // Step 1: switch back to the SAR TableView (the one with the PeptidesModel).
    await page.evaluate(async () => {
      const sarTv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']);
      if (sarTv) { grok.shell.v = sarTv; await new Promise((r) => setTimeout(r, 1000)); }
    });

    // Steps 2-3 (trigger via real UI): open the Export dialog from the SVM context menu so
    // the same surface the scenario describes is exercised, then Cancel it (the picker
    // gesture itself is the non-automatable leg) before the JS-API paired-column path.
    const dialogTrigger = await page.evaluate(async () => {
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]')!;
      const canvases = Array.from(svm.querySelectorAll('canvas')) as HTMLCanvasElement[];
      let canvas: HTMLCanvasElement | null = null;
      let maxArea = 0;
      for (const c of canvases) {
        const r = c.getBoundingClientRect();
        if (r.width * r.height > maxArea) { maxArea = r.width * r.height; canvas = c; }
      }
      if (!canvas) return {canvasFound: false, dlgFound: false, extraColInputFound: false};
      const r = canvas.getBoundingClientRect();
      const cx = r.x + 120, cy = r.y + 120;
      canvas.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: cx, clientY: cy, view: window}));
      canvas.dispatchEvent(new MouseEvent('contextmenu',
        {bubbles: true, cancelable: true, button: 2, clientX: cx, clientY: cy, view: window}));
      await new Promise((rs) => setTimeout(rs, 1200));
      const exportMc = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find((e) => e.textContent?.trim() === 'Export Mutation Cliffs...');
      if (exportMc) (exportMc.closest('.d4-menu-item') as HTMLElement | null)?.click();
      await new Promise((rs) => setTimeout(rs, 1500));
      const dlg = document.querySelector('[name="dialog-Export-Mutation-Cliffs"]');
      const extraColInput = dlg ? dlg.querySelector('[name="input-Extra-columns"]') : null;
      // Close the dialog — the picker gesture is non-automatable; the paired-column contract
      // is exercised via the JS-API export path below.
      const cancel = dlg ? dlg.querySelector('[name="button-CANCEL"]') as HTMLElement | null : null;
      if (cancel) cancel.click();
      await new Promise((rs) => setTimeout(rs, 500));
      return {canvasFound: true, dlgFound: !!dlg, extraColInputFound: !!extraColInput};
    });
    expect(dialogTrigger.canvasFound, 'SVM render canvas not found for Scenario 2 trigger').toBe(true);
    expect(dialogTrigger.dlgFound,
      'Export Mutation Cliffs dialog did not re-open from the SVM context menu').toBe(true);
    expect(dialogTrigger.extraColInputFound,
      'Export dialog is missing the [name="input-Extra-columns"] picker').toBe(true);

    // Step 4: exercise the paired-column emission contract via the viewer's production
    // export path with one extra non-default column. The available extra columns per the
    // source filter exclude the activity (IC50), the sequence (AlignedSequence), and the
    // per-position columns; `ID` is the first such column on the demo dataset.
    const exported = await page.evaluate(() => {
      const sarTv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel'])!;
      const svm = Array.from(sarTv.viewers).find((v) => v.type === 'Sequence Variability Map') as any;
      const df = sarTv.dataFrame;
      const idCol = df.col('ID');
      const extraName = idCol ? idCol.name : null;
      const existingMcIds = new Set(
        Array.from(grok.shell.tableViews)
          .filter((v) => v.dataFrame && v.dataFrame.name.startsWith('Mutation Cliffs'))
          .map((v) => v.dataFrame.id));
      // Production export path (same code the dialog OK invokes via exportMutationCliffs()).
      svm._doExportMutationCliffs([idCol]);
      return {extraName, existingCount: existingMcIds.size};
    });
    expect(exported.extraName, 'demo dataset is missing the `ID` extra column').toBe('ID');

    // Allow the new view to attach, then read the freshest Mutation Cliffs view.
    await page.waitForTimeout(2500);
    const result = await page.evaluate((baseCols) => {
      const mcViews = Array.from(grok.shell.tableViews)
        .filter((v) => v.dataFrame && v.dataFrame.name.startsWith('Mutation Cliffs'));
      // The paired-column view is the one carrying `Seq 1 ID` / `Seq 2 ID`.
      const paired = mcViews.find((v) =>
        v.dataFrame.columns.names().includes('Seq 1 ID') &&
        v.dataFrame.columns.names().includes('Seq 2 ID'));
      if (!paired) {
        return {found: false, anyCols: mcViews.map((v) => v.dataFrame.columns.names())};
      }
      const cols = paired.dataFrame.columns.names();
      return {
        found: true,
        cols,
        rows: paired.dataFrame.rowCount,
        hasAllBase: baseCols.every((c: string) => cols.includes(c)),
      };
    }, BASE_EXPORT_COLS);

    expect(result.found,
      `exported table did not carry the Seq 1 ID / Seq 2 ID paired columns; saw: ${JSON.stringify(result.anyCols)}`)
      .toBe(true);
    // Base columns still present.
    expect(result.hasAllBase, 'paired-column export lost one or more base columns').toBe(true);
    // One paired pair (Seq 1 <extra> / Seq 2 <extra>) for the chosen extra column.
    expect(result.cols, 'paired export missing the Seq-1-value paired column').toContain('Seq 1 ID');
    expect(result.cols, 'paired export missing the Seq-2-value paired column').toContain('Seq 2 ID');
    // Total = 6 base + 2 paired columns.
    expect(result.cols!.length, 'paired export should add exactly one Seq-1/Seq-2 pair (8 columns total)')
      .toBe(BASE_EXPORT_COLS.length + 2);
    expect(result.rows, 'paired-column export grid is empty').toBeGreaterThan(0);
  });

  // Cleanup.
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
