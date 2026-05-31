/* ---
sub_features_covered: [peptides.workflow.sar-dialog, peptides.workflow.analyze-ui, peptides.workflow.start-analysis, peptides.model.add-sequence-space, peptides.viewers.cluster-max-activity, peptides.viewers.logo-summary-table, peptides.widgets.settings-dialog, peptides.compute.calculate-cluster-statistics]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: integration
//   sub_features_covered: [peptides.workflow.sar-dialog, peptides.workflow.analyze-ui,
//     peptides.workflow.start-analysis, peptides.model.add-sequence-space,
//     peptides.viewers.cluster-max-activity, peptides.viewers.logo-summary-table,
//     peptides.widgets.settings-dialog, peptides.compute.calculate-cluster-statistics]
//   ui_coverage_responsibility: [top-menu-bio-analyze-sar, sar-settings-wrench-button,
//     mcl-viewer-parameter-driven-rerender] (delegated_to: info-panels.md)
//   related_bugs: [GROK-19145, GROK-14357, GROK-18058, github-1549, GROK-14461]
//
// Atlas provenance (derived_from):
//   feature-atlas/peptides.yaml#critical_paths[peptide-space-sar-with-mcl] derived_from:
//     public/packages/UsageAnalysis/files/TestTrack/Peptides/peptide-space-spec.ts
//
// Peptide Space — SAR launch via the TOP-MENU (Bio | Analyze | SAR...) entry path,
// with sequence-space dim-reduction + MCL clustering. Sister of the context-panel
// Launch SAR path (sar.md / sar-spec.ts). Maps to atlas critical_path
// peptide-space-sar-with-mcl (p1) + interaction sar-with-sequence-space-and-mcl.
// Selectors verified against dev.datagrok.ai @datagrok/peptides v1.27.9 per
// .claude/skills/grok-browser/references/peptides.md.
//
// Empirical recon notes (live MCP recon 2026-05-29 — drive deterministic
// assertions, not theory):
//   - The top-menu path is deterministic at the spec viewport (1920x1080,
//     specTestOptions): [name="div-Bio"] IS rendered + visible (offsetParent
//     non-null). The grok-browser overflow caveat (Bio menu in overflow) was
//     observed at 1536px CSS only — at 1920px the menubar shows
//     Edit | View | Select | Data | ML | Bio and the path drives cleanly:
//       [name="div-Bio"] .click()
//       -> [name="div-Bio---Analyze"] mouseenter/mousemove
//       -> [name="div-Bio---Analyze---SAR..."] .click()
//       -> [name="dialog-Analyze-Peptides"] (title "Analyze Peptides") opens.
//   - Default OK on the Analyze Peptides dialog attaches (verified live):
//     Grid, Sequence Variability Map, Most Potent Residues, MCL. Logo Summary
//     Table is NOT in the default attach set (cluster-/settings-dependent) — the
//     scenario names it (sub_features_covered) but the build does not auto-attach
//     it on a default top-menu launch; recorded softly, not failed.
//   - There is NO "Cluster (MCL)" column on this build (df columns after launch:
//     Activity, AlignedSequence, ID, IC50, position cols 1..17 — no cluster
//     column). The MCL viewer (type "MCL") is a dim-reduction scatter that shares
//     the main DataFrame; cluster state lives on the model, not a df column. An
//     earlier spec asserting df.col('Cluster (MCL)') would time out — that arm is
//     not assertable here.
//   - Scenario 2 settings flow is the wrench icon
//     (i.grok-icon.fa-wrench[aria-label="Peptides analysis settings"]) ->
//     [name="dialog-Peptides-settings"] accordion (General / Viewers / Columns /
//     MCL / Sequence space) — NOT a SAR re-launch.
//   - MCL re-render contract (per scope_reductions SR-01/SR-02): changing one
//     representative MCL parameter (Inflation-Factor 1.4 -> 2.5) + OK propagates
//     to model._settings.mclSettings.inflation and re-runs the clustering
//     pipeline. The deterministic re-render observable is the applied model
//     setting (no cluster-count column to diff on this build), plus the MCL
//     viewer persisting + re-rendering its canvas, plus the GROK-19145 invariant
//     (no post-OK null-receiver crash). The scenario's "cluster count delta /
//     membership change" arm is deferred to a future spec that can capture a
//     cluster-count observable (see SR-02 + the scenario's
//     unresolved_ambiguities).

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/bio/peptides.csv';

test('Peptide Space — top-menu SAR launch with sequence-space + MCL clustering', async ({page}) => {
  // SAR launch + sequence-space embedding + MCL clustering all run async in Web
  // Workers; a settings recompute follows. Multiple multi-second waits won't fit
  // the default per-test budget.
  test.setTimeout(300_000);
  await loginToDatagrok(page);

  // Setup: open the peptides dataset (implicit prerequisite promoted to an
  // explicit Setup step in the scenario — the top-menu SAR launcher gates on an
  // active Macromolecule column + a numerical activity column). Wait for semType
  // detection + Bio package settle, then open the Grid.
  await softStep('Setup: open the peptides Macromolecule table', async () => {
    const result = await page.evaluate(async (path) => {
      document.querySelectorAll('.d4-dialog').forEach((d) => {
        const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
        if (cancel) cancel.click();
      });
      grok.shell.closeAll();
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      // Windows mode (simpleMode=false): the top-menu SAR path docks several
      // viewers; the sibling sar-spec uses the same mode.
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

      // Pre-warm the Peptides @init so PeptideUtils.loadComponents() (SeqHelper +
      // MonomerLib) is done before the SAR launch requests it — removes the
      // cold-package init leg (GROK-17557 family). The dialog/launch waits below
      // are the authoritative readiness gates regardless.
      try { await grok.functions.call('Peptides:initPeptides'); }
      catch (e) { console.log('[note] Peptides:initPeptides pre-warm threw (non-fatal):', String(e)); }

      return {
        rows: df.rowCount,
        semType: df.col('AlignedSequence')?.semType ?? null,
      };
    }, datasetPath);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
    expect(result.rows, 'peptides.csv should load 647 rows').toBe(647);
    expect(result.semType, 'AlignedSequence must be a Macromolecule column').toBe('Macromolecule');
  });

  // ---- Scenario 1 — Launch SAR from the top menu (Bio | Analyze | SAR...) ----

  // Step 1: invoke Bio | Analyze | SAR... from the top menu. The Analyze Peptides
  // config dialog opens (atlas peptides.workflow.sar-dialog -> analyze-ui).
  await softStep('Scenario 1 (step 1): launch SAR via Bio | Analyze | SAR... top menu', async () => {
    const opened = await page.evaluate(async () => {
      // Recon (2026-05-29) confirms the Bio top menu is rendered + visible at the
      // spec viewport (1920px); the overflow caveat applied only at <~1600px CSS.
      const bio = document.querySelector('[name="div-Bio"]') as HTMLElement | null;
      const bioVisible = bio ? bio.offsetParent !== null : false;
      if (bio) bio.click();
      await new Promise((r) => setTimeout(r, 700));
      const analyze = document.querySelector('[name="div-Bio---Analyze"]') as HTMLElement | null;
      if (analyze) {
        analyze.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        analyze.dispatchEvent(new MouseEvent('mousemove', {bubbles: true}));
      }
      await new Promise((r) => setTimeout(r, 700));
      const sar = document.querySelector('[name="div-Bio---Analyze---SAR..."]') as HTMLElement | null;
      if (sar) sar.click();
      await new Promise((r) => setTimeout(r, 2500));
      const dlg = document.querySelector('[name="dialog-Analyze-Peptides"]');
      return {
        bioFound: !!bio,
        bioVisible,
        analyzeFound: !!analyze,
        sarFound: !!sar,
        dialogFound: !!dlg,
      };
    });
    // Surface menu-shape drift loudly rather than silently failing the OK below.
    expect(opened.bioFound, '[name="div-Bio"] top-menu entry not found').toBe(true);
    expect(opened.bioVisible,
      '[name="div-Bio"] present but not visible (offsetParent null — menu in overflow at this viewport)').toBe(true);
    expect(opened.analyzeFound, '[name="div-Bio---Analyze"] submenu entry not found').toBe(true);
    expect(opened.sarFound, '[name="div-Bio---Analyze---SAR..."] menu item not found').toBe(true);
    expect(opened.dialogFound, '[name="dialog-Analyze-Peptides"] config dialog did not open').toBe(true);
  });

  // Step 2: accept the default config, click OK, wait for the SAR analysis layout
  // (PeptidesModel attaches; sequence-space + MCL compute run in Web Workers).
  await softStep('Scenario 1 (step 2): accept default config, run SAR, verify viewers', async () => {
    await page.evaluate(async () => {
      const dlg = document.querySelector('[name="dialog-Analyze-Peptides"]');
      const ok = (dlg?.querySelector('[name="button-OK"]')
        ?? document.querySelector('[name="button-OK"]')) as HTMLElement | null;
      if (ok) ok.click();
    });

    // Calculation-completion marker (scenario unresolved_ambiguity
    // step-2-calculation-completion-marker-needed): the analysis completes when
    // the PeptidesModel singleton attaches to a table view's DataFrame. Poll for
    // it rather than a fixed wait.
    await page.waitForFunction(() => {
      return Array.from(grok.shell.tableViews).some((v) => v.dataFrame.temp['peptidesModel']);
    }, {timeout: 90000});

    // Let the sequence-space embedding + MCL clustering settle before probing.
    await page.waitForTimeout(8000);

    const viewers = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      return Array.from(tv.viewers).map((v) => v.type);
    });
    // Deterministic default top-menu attach set (verified live 2026-05-29):
    // Sequence Variability Map + Most Potent Residues + MCL (the sequence-space +
    // MCL cluster viewer). These cover sub_features add-sequence-space +
    // cluster-max-activity + calculate-cluster-statistics.
    expect(viewers, 'Sequence Variability Map must attach after SAR launch').toContain('Sequence Variability Map');
    expect(viewers, 'Most Potent Residues must attach after SAR launch').toContain('Most Potent Residues');
    expect(viewers, 'MCL clustering viewer must attach after SAR launch').toContain('MCL');
    // Logo Summary Table: scenario names it (sub_features_covered) but the
    // current build does not auto-attach it on a default top-menu launch (it is
    // cluster-/settings-dependent). Record without failing the deterministic
    // viewer contract.
    if (!viewers.includes('Logo Summary Table'))
      console.log('[note] Logo Summary Table not in default top-menu SAR attach set on this build (cluster-/settings-dependent).');
  });

  // ---- Scenario 2 — Change an MCL setting via the wrench, verify MCL re-render ----

  // Step 3: open the SAR settings dialog via the wrench-icon button on the SAR
  // analysis toolbar (atlas getSettingsDialog(model) — 5-pane accordion).
  await softStep('Scenario 2 (step 3): open the SAR settings dialog via the wrench', async () => {
    const opened = await page.evaluate(async () => {
      // The wrench carries no name= — anchor on its aria-label.
      const wrench = document.querySelector(
        'i.grok-icon.fa-wrench[aria-label="Peptides analysis settings"]') as HTMLElement | null;
      const wrenchFound = !!wrench;
      if (wrench) wrench.click();
      await new Promise((r) => setTimeout(r, 2000));
      const dlg = document.querySelector('[name="dialog-Peptides-settings"]');
      const panes = dlg
        ? Array.from(dlg.querySelectorAll('.d4-accordion-pane'))
          .map((p) => p.querySelector('.d4-accordion-pane-header')?.textContent?.trim())
        : [];
      return {wrenchFound, dlgFound: !!dlg, panes};
    });
    expect(opened.wrenchFound, 'Peptides analysis settings wrench not found on the SAR toolbar').toBe(true);
    expect(opened.dlgFound, '[name="dialog-Peptides-settings"] did not open').toBe(true);
    // The MCL pane (where the representative parameter lives) must be present.
    expect(opened.panes, 'settings dialog must expose the MCL accordion pane').toContain('MCL');
  });

  // Step 4: in the MCL section change one representative parameter (Inflation
  // Factor 1.4 -> 2.5) to a non-default value, then OK. Per scope_reductions
  // SR-01 the exhaustive per-parameter matrix is deferred to a future
  // parameterized peptides-settings-dialog spec; this exercises the OK roundtrip
  // with a single MCL-parameter change to establish the settings -> MCL re-render
  // contract.
  await softStep('Scenario 2 (step 4): change a representative MCL parameter + OK', async () => {
    const changed = await page.evaluate(async () => {
      const dlg = document.querySelector('[name="dialog-Peptides-settings"]')!;
      // Capture the applied MCL inflation on the model BEFORE the change — this
      // is the deterministic re-render observable (no cluster-count column exists
      // on this build to diff on).
      const tvBefore = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']);
      const inflationBefore = tvBefore
        ? (tvBefore.dataFrame.temp['peptidesModel'] as any)?._settings?.mclSettings?.inflation
        : null;

      // Expand the MCL pane.
      const mclPane = Array.from(dlg.querySelectorAll('.d4-accordion-pane'))
        .find((p) => p.querySelector('.d4-accordion-pane-header')?.textContent?.trim() === 'MCL');
      if (mclPane) {
        const h = mclPane.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
        if (h && !mclPane.classList.contains('expanded')) h.click();
      }
      await new Promise((r) => setTimeout(r, 800));

      // Representative MCL parameter: Inflation Factor (governs cluster
      // granularity). Default 1.4 -> 2.5.
      const host = dlg.querySelector('[name="input-host-Inflation-Factor"]');
      const inner = host ? (host.querySelector('input') as HTMLInputElement | null) : null;
      const inputFound = !!inner;
      const prevVal = inner ? inner.value : null;
      if (inner) {
        // Drive the Dart change listener: focus, set, dispatch input+change, blur.
        inner.focus();
        inner.value = '2.5';
        inner.dispatchEvent(new Event('input', {bubbles: true}));
        inner.dispatchEvent(new Event('change', {bubbles: true}));
        inner.blur();
      }
      await new Promise((r) => setTimeout(r, 600));
      const newVal = inner ? inner.value : null;

      // Click OK.
      const ok = dlg.querySelector('[name="button-OK"]') as HTMLElement | null;
      if (ok) ok.click();
      await new Promise((r) => setTimeout(r, 2000));
      const dlgClosed = !document.querySelector('[name="dialog-Peptides-settings"]');
      return {inputFound, prevVal, newVal, inflationBefore, dlgClosed};
    });
    expect(changed.inputFound, '[name="input-host-Inflation-Factor"] not found in the MCL pane').toBe(true);
    expect(changed.newVal, 'MCL Inflation Factor did not accept the new value').toBe('2.5');
    expect(changed.inflationBefore, 'pre-change MCL inflation should be the default 1.4').toBe(1.4);
    expect(changed.dlgClosed, 'settings dialog did not close after OK').toBe(true);
  });

  // Step 5: verify the MCL Viewer re-renders after the OK roundtrip. Per
  // scope_reductions SR-02 this is sharpened to a deterministic property — the
  // applied MCL setting propagated to the model (re-render cause), the MCL viewer
  // persists + re-renders its canvas, and the model survives (GROK-19145
  // invariant: no post-OK null-receiver crash). The source's looser "MCL Viewer
  // should output different results" cluster-by-cluster diff remains an
  // unresolved ambiguity for a future spec (no cluster-count observable here).
  await softStep('Scenario 2 (step 5): verify the MCL Viewer re-renders after the settings change', async () => {
    // Let the MCL recompute (Web Worker) settle.
    await page.waitForTimeout(12000);
    const state = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const viewers = Array.from(tv.viewers).map((v) => v.type);
      const model = tv.dataFrame.temp['peptidesModel'] as any;
      const inflationAfter = model?._settings?.mclSettings?.inflation ?? null;
      const mclDom = document.querySelector('[name="viewer-MCL"]');
      // The MCL viewer re-renders its scatter onto its largest canvas — assert it
      // has rendered (non-background) content.
      let mclHasRender = false;
      if (mclDom) {
        const cs = Array.from(mclDom.querySelectorAll('canvas')) as HTMLCanvasElement[];
        let best: HTMLCanvasElement | null = null, area = 0;
        for (const c of cs) { const r = c.getBoundingClientRect(); if (r.width * r.height > area) { area = r.width * r.height; best = c; } }
        if (best && best.width > 0 && best.height > 0) {
          try {
            const ctx = best.getContext('2d')!;
            const data = ctx.getImageData(0, 0, best.width, best.height).data;
            let nonBg = 0;
            for (let i = 0; i < data.length; i += 41)
              if (data[i + 3] > 0 && (data[i] < 250 || data[i + 1] < 250 || data[i + 2] < 250)) nonBg++;
            mclHasRender = nonBg > 50;
          } catch (e) { mclHasRender = false; }
        }
      }
      const lastError = grok.shell.lastError ? String(grok.shell.lastError) : '';
      return {viewers, modelPresent: !!model, inflationAfter, mclHasRender, lastError};
    });
    // The applied MCL setting changed — proves the settings -> recompute contract
    // fired (the deterministic re-render cause).
    expect(state.inflationAfter, 'MCL Inflation Factor change did not propagate to the model').toBe(2.5);
    // The model + MCL viewer survive the recompute and re-render.
    expect(state.modelPresent, 'PeptidesModel cache lost after the MCL settings change').toBe(true);
    expect(state.viewers, 'MCL viewer must persist after the settings change').toContain('MCL');
    expect(state.mclHasRender, 'MCL viewer did not re-render visible content after the settings change').toBe(true);
    // GROK-19145 invariant: no fatal "setTrue on null" / null-receiver crash in
    // the post-OK compute path. (Benign async/Promise + resource-404 noise is
    // tolerated; only null-receiver crashes fail the invariant.)
    expect(/setTrue|null/.test(state.lastError),
      `GROK-19145 invariant: post-OK MCL compute produced a null-receiver error: ${state.lastError}`).toBe(false);
  });

  // Cleanup.
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
