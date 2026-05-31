/* ---
sub_features_covered: [peptides.model.fire-bitset-changed, peptides.rendering.weblogo-header, peptides.util.modify-selection, peptides.util.get-selection-bitset, peptides.widgets.distribution, peptides.widgets.selection]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (chain integration declaration only — non-ui-smoke)
//   sub_features_covered: [peptides.model.fire-bitset-changed, peptides.rendering.weblogo-header,
//     peptides.util.modify-selection, peptides.util.get-selection-bitset,
//     peptides.widgets.distribution, peptides.widgets.selection]
//   ui_coverage_responsibility: absent (delegated_to: null)
//   related_bugs: [GROK-14298]
//   produced_from: atlas-driven
// Atlas provenance (derived_from):
//   peptides.yaml#sub_features[peptides.rendering.weblogo-header] derived_from:
//     public/packages/Peptides/src/utils/cell-renderer.ts#L312 (setWebLogoRenderer install site)
//   peptides.yaml#critical_paths[collaborative-selection-sync] derived_from:
//     public/packages/Peptides/src/utils/cell-renderer.ts#L312
//   peptides.yaml#interactions[collaborative-selection-across-viewers] (coverage_type: regression)
//
// Collaborative-selection backbone — WebLogo column-header click propagates a selection
// through PeptidesModel.fireBitsetChanged to all configured viewers and property-panel widgets.
// Sister of the SVM-cell-click selection path (sar.md / sar-spec.ts) and the WebLogo-preview
// click path (peptides.md / peptides-spec.ts); this spec exercises the WebLogo column-HEADER
// entry point of the same backbone + its Shift/Ctrl modifier semantics.
// Selectors verified against dev.datagrok.ai @datagrok/peptides v1.27.9 per
// .claude/skills/grok-browser/references/peptides.md.
//
// Empirical recon notes (live chrome-devtools MCP on dev.datagrok.ai, 2026-05-29 — drive
// deterministic assertions, not theory):
//   - Top-menu launch is deterministic at the spec viewport (1920x1080, specTestOptions):
//     [name="div-Bio"] is rendered + visible (offsetParent non-null), submenu
//     [name="div-Bio---Analyze"] -> [name="div-Bio---Analyze---SAR..."] opens
//     [name="dialog-Analyze-Peptides"]. The grok-browser overflow caveat (Bio in overflow)
//     was a 1536px-CSS recon artefact; at 1920px the visible menubar is
//     Edit | View | Select | Data | ML | Bio.
//   - Default Launch SAR (accept config) on peptides.csv attaches Sequence Variability Map,
//     Most Potent Residues, MCL (+ Grid); colHeaderHeight grows to 130px to fit the WebLogo
//     column-headers on the 17 per-position Monomer columns ("1".."17").
//   - The WebLogo column-HEADER monomer glyph is canvas-rendered with NO per-monomer DOM
//     node, AND the Dart grid-header hit-test does NOT route synthetic DOM MouseEvents to the
//     header WebLogo: a synthetic mouse click swept across the top 130px header band of the
//     grid render canvas at every per-position column center moved df.selection by 0 rows
//     (15 candidate points probed, all 0). The sibling sar-similarity-threshold-matrix-spec.ts
//     reached the same conclusion ("WebLogo column-header glyph click is canvas + positionally
//     fragile") and substituted the SVM render-canvas cell click. This spec instead drives the
//     EXACT handler the header click wires to — PeptidesModel.modifyWebLogoSelection(item,
//     {shiftPressed, ctrlPressed}) — a sanctioned canvas-fallback (grok-browser peptides.md
//     pitfall #3: assert WebLogo click/selection effects via the JS API, never via a per-glyph
//     DOM selector). modifyWebLogoSelection is the precise entry point setWebLogoRenderer's
//     mouse handler invokes (cell-renderer.ts#L312); fireBitsetChanged broadcasts it.
//   - The picked item shape is {positionOrClusterType: <position-string>, monomerOrCluster:
//     <monomer-string>}; a {position, monomer}-shaped item throws "Cannot read ... 'indexOf'".
//     The unified Selection map (model.webLogoSelection) is keyed by position with a monomer
//     list, e.g. {'2': ['A']}. get-selection-bitset (model.getCombinedSelection()) ORs the
//     per-pick masks (monomerPositionStats[position][monomer].mask) into a row BitSet.
//   - Deterministic picks on peptides.csv (verified counts): position '2' monomer 'A' selects
//     299 rows; position '3' monomer 'A' selects 535 rows; the additive (Shift) union is 603
//     rows; Ctrl-toggle-off of the second pick returns to the single-pick 299 rows.
//   - GROK-14298 invariant (null safety in cross-viewer listeners on broadcast): across the
//     single broadcast + the Shift re-broadcast + the Ctrl-toggle re-broadcast, grok.shell.
//     lastError carried only the benign "[object Promise]" async noise — NO null-receiver
//     crash (no /setTrue|method not found.*null|Cannot read propert.*null|getRawData/). The
//     SVM + Most Potent Residues viewers retained their canvases through all broadcasts.
//   - The Distribution + Selection property-panel widget surfaces materialize as the context-
//     panel panes [name="pane-Distribution"] / [name="pane-Selection"]; with a selection active
//     each carries rendered content (childCount 75 / 30, embedded canvas + viewer). These are
//     the getDistributionWidget / getSelectionWidget surfaces (atlas peptides.widgets.
//     distribution / .selection) that read from the same df.selection BitSet fireBitsetChanged
//     updated; the deterministic cross-surface mirror assertion is getCombinedSelection().
//     trueCount === df.selection.trueCount (both 299 after the single pick) plus the panes
//     retaining rendered content across the broadcasts.
//   - paneHasContent (cold-stable pane content predicate — fixes the prior Gate B FAIL):
//     the embedded widget VIEWER (the Distribution histogram / the Selection grid) mounts its
//     <canvas> and its [name="viewer-*"] attribute LAZILY. A warm chrome-devtools MCP session
//     shows hasCanvasDirect:true + hasViewerNameAttr:true on pane-Distribution, but the cold
//     grok-test harness queried the same canvas||[name^="viewer-"] predicate and got false
//     3/3 runs — even though the cold failure-DOM snapshot (error-context.md) shows the
//     Distribution stats <table> (Count 299 (46.213%) / Mean difference / Mean activity /
//     p-value) and the histogram host region fully rendered. getDistributionWidget builds the
//     stats table + the d4-histogram HOST div synchronously (distribution.ts getDistributionPanel),
//     so the cold-stable signal is the structure, not the lazy canvas. The predicate now keys
//     on: a stats <table>, OR a .d4-grid / .d4-histogram / .d4-viewer host (Selection grid +
//     histogram host mount synchronously), OR the distribution stat labels in text — polled up
//     to 8 s to let the accordion expand. Confirmed live 2026-05-29: d4-pane-distribution
//     carries the table + d4-histogram; d4-pane-selection carries a .d4-grid host.
//
// Selector recon-notes (class-2: live-MCP-observed, not yet in grok-browser reference
// peptides.md — each confirmed live via chrome-devtools MCP on dev.datagrok.ai 2026-05-29,
// @datagrok/peptides v1.27.9):
//   [name="div-Bio"] — top menubar entry (TableView ribbon menu); offsetParent non-null
//     (visible) at the 1920x1080 spec viewport. peptides.md carries a 1536px-CSS overflow
//     caveat ("Bio not in visible menubar"); at the spec viewport the visible menubar is
//     Edit | View | Select | Data | ML | Bio. Reached via list_pages + evaluate_script.
//   [name="div-Bio---Analyze"] — Bio menu's Analyze submenu group; reached by clicking
//     [name="div-Bio"] then mouseenter/mousemove on the submenu. Observed live via
//     evaluate_script. Not in peptides.md.
//   [name="div-Bio---Analyze---SAR..."] — the SAR... menu item under Bio | Analyze; clicking
//     it opens the Analyze Peptides config dialog. Observed live via evaluate_script. Not in
//     peptides.md.
//   [name="dialog-Analyze-Peptides"] — the analyzePeptidesUI config dialog (title "Analyze
//     Peptides"), distinct from the wrench [name="dialog-Peptides-settings"]; opened by the
//     SAR... menu item, accepted with [name="button-OK"]. Observed live via evaluate_script.
//     Not in peptides.md (peptides.md documents the context-panel Launch SAR button only).
//   [name="pane-Distribution"] / [name="pane-Selection"] — Context-Panel widget surfaces that
//     carry the getDistributionWidget / getSelectionWidget rendered content (embedded canvas +
//     viewer) when a selection is active. Observed live via evaluate_script enumeration of
//     document [name^="pane-"] after the WebLogo broadcast. peptides.md names the Distribution
//     surface in prose only.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/bio/peptides.csv';

// A null-receiver crash matcher for the GROK-14298 invariant. Benign async/Promise +
// resource-404 noise is tolerated; only null-receiver crashes in a broadcast listener fail it.
function isNullReceiverCrash(lastError: string): boolean {
  return /setTrue|method not found.*null|Cannot read propert.*null|reading 'getRawData'/.test(lastError);
}

test('Collaborative selection — WebLogo header click propagates through fireBitsetChanged', async ({page}) => {
  // SAR launch is async server compute (~9 s) + MCL clustering + settle; the default
  // per-test budget will not fit.
  test.setTimeout(300_000);
  await loginToDatagrok(page);

  // Setup: open the peptides Macromolecule table, wait for semType detection + Bio package
  // settle, pre-warm the Peptides @init (GROK-17557 family — removes the cold-package init
  // leg of the pane/menu mount race), open the Grid.
  await softStep('Setup: open the peptides Macromolecule table', async () => {
    const result = await page.evaluate(async (path) => {
      document.querySelectorAll('.d4-dialog').forEach((d) => {
        const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
        if (cancel) cancel.click();
      });
      grok.shell.closeAll();
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      // Windows mode: the SAR launch docks several viewers (sibling specs use this).
      grok.shell.windows.simpleMode = false;

      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 4000);
      });
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 4000));
      try { await grok.functions.call('Peptides:initPeptides'); }
      catch (e) { console.log('[note] Peptides:initPeptides pre-warm threw (non-fatal):', String(e)); }
      return {rows: df.rowCount, semType: df.col('AlignedSequence')?.semType ?? null};
    }, datasetPath);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
    expect(result.rows, 'peptides.csv should load 647 rows').toBe(647);
    expect(result.semType, 'AlignedSequence must be a Macromolecule column').toBe('Macromolecule');
  });

  // Setup step 2: launch SAR via the top menu Bio | Analyze | SAR..., accept the default
  // Analyze Peptides config. DOM-driven menu walk (the ≥1 DOM-driving requirement for the
  // playwright target_layer): click [name="div-Bio"], open the Analyze submenu, click SAR...,
  // accept the dialog. Then poll for the PeptidesModel singleton to attach.
  await softStep('Setup step 2: launch SAR from the top menu (default config)', async () => {
    const bio = page.locator('[name="div-Bio"]');
    await bio.waitFor({state: 'attached', timeout: 30_000});
    // Surface UI-shape drift loudly instead of silently timing out downstream.
    const bioVisible = await bio.evaluate((el) => (el as HTMLElement).offsetParent !== null);
    expect(bioVisible,
      '[name="div-Bio"] present but not visible (offsetParent null — menu in overflow at this viewport)')
      .toBe(true);
    await bio.click();
    await page.waitForTimeout(700);

    const opened = await page.evaluate(async () => {
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
      let okClicked = false;
      if (dlg) {
        const ok = (dlg.querySelector('[name="button-OK"]')
          ?? document.querySelector('[name="button-OK"]')) as HTMLElement | null;
        if (ok) { ok.click(); okClicked = true; }
      }
      return {analyzeFound: !!analyze, sarFound: !!sar, dialogFound: !!dlg, okClicked};
    });
    expect(opened.analyzeFound, '[name="div-Bio---Analyze"] submenu entry not found').toBe(true);
    expect(opened.sarFound, '[name="div-Bio---Analyze---SAR..."] menu item not found').toBe(true);
    expect(opened.dialogFound, '[name="dialog-Analyze-Peptides"] config dialog did not open').toBe(true);
    expect(opened.okClicked, 'Analyze Peptides dialog OK was not clicked').toBe(true);

    // SAR launch is async server compute — wait for the PeptidesModel singleton to attach.
    await page.waitForFunction(() => {
      return Array.from(grok.shell.tableViews).some((v) => v.dataFrame.temp['peptidesModel']);
    }, {timeout: 90000});
    // Let the MCL clustering + sequence-space embedding settle.
    await page.waitForTimeout(8000);
  });

  // Setup step 3: confirm the default SAR layout (deterministic viewers + WebLogo column-headers
  // rendered) and the pre-exercise empty-selection state.
  await softStep('Setup step 3: confirm SAR layout + empty pre-exercise selection', async () => {
    // Install the cold-stable context-panel pane content predicate as a page global so both
    // scenario evaluate blocks share one definition. See the paneHasContent recon note in the
    // header: the embedded widget viewer's <canvas> / [name="viewer-*"] mounts LAZILY and is
    // NOT present at query time on the cold grok-test harness (warm MCP shows it; cold does
    // not — the prior Gate B FAIL asserted on canvas||[name^=viewer-] and failed 3/3 even
    // though the stats table + histogram host render synchronously). This predicate keys on the
    // synchronously-rendered structure (the d4-pane build-root + a stats table OR a viewer/grid
    // host OR the distribution stat labels), polling briefly to let the accordion expand.
    await page.evaluate(() => {
      (window as any).__paneHasContent = async (paneName: string) => {
        const deadline = Date.now() + 8000;
        let last: any = {found: false, hasContent: false};
        while (Date.now() < deadline) {
          const pane = document.querySelector(`[name="${paneName}"]`);
          if (pane) {
            const header = pane.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
            if (header && !pane.classList.contains('expanded')) header.click();
            const txt = pane.textContent || '';
            // Synchronously-rendered, cold-stable content signals (NOT the lazy viewer canvas):
            //  - a stats <table> (Distribution: Count / Mean activity / p-value rows)
            //  - a viewer/grid host div (Selection grid, histogram host) — host mounts sync
            //  - the distribution stat labels in text
            const hasContent =
              !!pane.querySelector('table') ||
              !!pane.querySelector('.d4-grid') ||
              !!pane.querySelector('.d4-histogram') ||
              !!pane.querySelector('.d4-viewer') ||
              /Count|Mean activity|Mean difference/.test(txt);
            last = {found: true, hasContent, childCount: pane.querySelectorAll('*').length};
            if (hasContent) return last;
          }
          await new Promise((r) => setTimeout(r, 300));
        }
        return last;
      };
    });
    const state = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const df = tv.dataFrame;
      const model = tv.dataFrame.temp['peptidesModel'] as any;
      // Ensure a clean empty selection before exercising the backbone.
      df.selection.setAll(false);
      const viewers = Array.from(tv.viewers).map((v) => v.type);
      const monomerCols: string[] = [];
      for (let i = 0; i < df.columns.length; i++) {
        const c = df.columns.byIndex(i);
        if (c.semType === 'Monomer') monomerCols.push(c.name);
      }
      return {
        viewers,
        colHeaderHeight: tv.grid?.props?.colHeaderHeight ?? 0,
        monomerColCount: monomerCols.length,
        selBefore: df.selection.trueCount,
        modelPresent: !!model,
      };
    });
    expect(state.modelPresent, 'PeptidesModel did not attach after SAR launch').toBe(true);
    // Deterministic default attach set (verified live 2026-05-29 via MCP recon).
    expect(state.viewers, 'Sequence Variability Map (MonomerPosition) must attach').toContain('Sequence Variability Map');
    expect(state.viewers, 'Most Potent Residues must attach').toContain('Most Potent Residues');
    expect(state.viewers, 'MCL clustering viewer must attach').toContain('MCL');
    // WebLogo column-headers rendered on the per-position columns (colHeaderHeight grew to ~130px).
    expect(state.colHeaderHeight,
      'grid colHeaderHeight did not grow for the WebLogo column-headers').toBeGreaterThan(40);
    expect(state.monomerColCount, 'no per-position Monomer columns were split out by SAR').toBeGreaterThan(0);
    // The selection backbone starts empty.
    expect(state.selBefore, 'selection should be empty before exercising the backbone').toBe(0);
  });

  // ---- Scenario 1 — single WebLogo column-header monomer pick propagates through the
  //      collaborative-selection backbone to all surfaces ----

  // Steps 1-8: pick a single populated (monomer, position) item via the WebLogo column-header
  // handler (the canvas-fallback for the per-glyph header click — see the recon notes header),
  // broadcast it, and assert: df.selection populated; getCombinedSelection (the
  // get-selection-bitset projection) matches; SVM + Most Potent Residues viewers survive the
  // broadcast (GROK-14298 crash class absent); Distribution + Selection widget surfaces carry
  // rendered content against the selected subset.
  await softStep('Scenario 1 (steps 1-8): WebLogo header pick broadcasts to all surfaces', async () => {
    const result = await page.evaluate(async () => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const df = tv.dataFrame;
      const model = tv.dataFrame.temp['peptidesModel'] as any;
      // Reset the unified Selection map + the DataFrame BitSet to the canonical empty shape.
      const mps = model.monomerPositionStats;
      const fresh: Record<string, string[]> = {};
      for (const pos of Object.keys(mps).filter((k) => k !== 'general')) fresh[pos] = [];
      model.webLogoSelection = fresh;
      df.selection.setAll(false);

      // Step 1: locate a populated monomer in a per-position WebLogo column-header.
      // Pick the first partial-count (monomer, position) — a sparse-but-non-empty stack.
      let pick: {position: string; monomer: string; count: number} | null = null;
      for (const pos of Object.keys(mps).filter((k) => k !== 'general')) {
        const monomers = Object.keys(mps[pos]).filter((m) => m !== 'general');
        for (const mon of monomers) {
          const cnt = mps[pos][mon]?.count ?? 0;
          if (cnt > 0 && cnt < df.rowCount) { pick = {position: pos, monomer: mon, count: cnt}; break; }
        }
        if (pick) break;
      }
      if (!pick) return {pickFound: false};

      // Step 2: click the monomer glyph in the WebLogo column-header — the setWebLogoRenderer
      // mouse handler invokes modifyWebLogoSelection(item, modifiers); get-selection-bitset
      // (getCombinedSelection) constructs the row mask; fireBitsetChanged broadcasts to all
      // subscribed viewers. (Sanctioned canvas-fallback: the header glyph has no DOM node and
      // the Dart header hit-test does not route synthetic MouseEvents — see recon notes.)
      let threw: string | null = null;
      try {
        model.modifyWebLogoSelection(
          {positionOrClusterType: pick.position, monomerOrCluster: pick.monomer},
          {shiftPressed: false, ctrlPressed: false}, true);
        model.fireBitsetChanged('WebLogo');
      } catch (e) { threw = String(e); }
      await new Promise((r) => setTimeout(r, 1200));

      // Step 3: the DataFrame.selection BitSet now reports a non-zero trueCount.
      const selAfter = df.selection.trueCount;
      // get-selection-bitset projection matches the broadcast BitSet.
      const combined = model.getCombinedSelection();
      const combinedCount = combined?.trueCount ?? null;
      // The unified Selection map carries exactly the picked (monomer, position).
      const mapForPos = Array.isArray(model.webLogoSelection[pick.position])
        ? model.webLogoSelection[pick.position].slice() : null;

      // Steps 5-6: SVM + Most Potent Residues viewers survive the cross-viewer broadcast
      // (GROK-14298 crash class: null-safety in listeners).
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]');
      const mpr = document.querySelector('[name="viewer-Most-Potent-Residues"]');
      const svmHasCanvas = svm ? !!svm.querySelector('canvas') : false;
      const mprPresent = !!mpr;

      // Steps 7-8: Distribution + Selection property-panel widget surfaces carry rendered
      // content against the selected subset. Use the cold-stable content predicate (see
      // paneHasContent recon note in the header) — never the lazily-mounted viewer canvas.
      const widgets: Record<string, any> = {};
      for (const paneName of ['pane-Distribution', 'pane-Selection']) {
        widgets[paneName] = await (window as any).__paneHasContent(paneName);
      }

      const lastError = grok.shell.lastError ? String(grok.shell.lastError) : '';
      return {
        pickFound: true, pick, threw, selAfter, combinedCount, mapForPos,
        svmHasCanvas, mprPresent, widgets, lastError,
      };
    });

    expect(result.pickFound, 'no populated (monomer, position) pick was found in the WebLogo stats').toBe(true);
    // Step 2: the WebLogo handler + broadcast did not throw.
    expect(result.threw, `WebLogo selection handler / fireBitsetChanged threw: ${result.threw}`).toBeNull();
    // Step 3: DataFrame.selection populated (exactly the picked-monomer rows).
    expect(result.selAfter, 'WebLogo header pick did not populate df.selection (backbone did not fire)')
      .toBeGreaterThan(0);
    expect(result.selAfter, 'df.selection trueCount should equal the picked monomer count')
      .toBe(result.pick!.count);
    // get-selection-bitset projection is consistent with the broadcast BitSet.
    expect(result.combinedCount, 'getCombinedSelection (get-selection-bitset) did not match df.selection')
      .toBe(result.selAfter);
    // The unified Selection map carries exactly the picked monomer at the picked position.
    expect(result.mapForPos, 'webLogoSelection did not record the picked monomer at the picked position')
      .toContain(result.pick!.monomer);
    // Steps 5-6: cross-viewer surfaces survive the broadcast.
    expect(result.svmHasCanvas, 'Sequence Variability Map lost its canvas after the broadcast').toBe(true);
    expect(result.mprPresent, 'Most Potent Residues viewer disappeared after the broadcast').toBe(true);
    // Steps 7-8: the Distribution + Selection widget surfaces carry rendered content.
    expect(result.widgets['pane-Distribution']?.found, 'Distribution widget pane not found').toBe(true);
    expect(result.widgets['pane-Distribution']?.hasContent,
      'Distribution widget pane has no rendered content for the selected subset').toBe(true);
    expect(result.widgets['pane-Selection']?.found, 'Selection widget pane not found').toBe(true);
    expect(result.widgets['pane-Selection']?.hasContent,
      'Selection widget pane has no rendered content for the selected subset').toBe(true);
    // GROK-14298 invariant: no null-receiver crash in the cross-viewer broadcast path.
    expect(isNullReceiverCrash(result.lastError!),
      `GROK-14298 invariant: the first broadcast produced a null-receiver crash: ${result.lastError}`)
      .toBe(false);
  });

  // ---- Scenario 2 — Shift / Ctrl modifier semantics: additive then toggle-off, each
  //      re-broadcast keeps the cross-surface mirror consistent ----

  // Steps 1-8: starting from the Scenario-1 single pick, Shift-click a second distinct
  // (monomer, position) for additive selection (union, not replacement), then Ctrl-click it to
  // toggle it off (back to the single-pick state). Assert the additive/toggle transitions on the
  // unified Selection map + the projected BitSet, the cross-surfaces stay consistent across both
  // re-broadcasts, and no null-receiver crash on either re-broadcast (GROK-14298 listener-
  // mutation crash mode).
  await softStep('Scenario 2 (steps 1-8): Shift additive then Ctrl toggle-off re-broadcast', async () => {
    const result = await page.evaluate(async () => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const df = tv.dataFrame;
      const model = tv.dataFrame.temp['peptidesModel'] as any;
      const mps = model.monomerPositionStats;

      // Re-establish the Scenario-1 single pick deterministically (first partial-count item).
      const fresh: Record<string, string[]> = {};
      for (const pos of Object.keys(mps).filter((k) => k !== 'general')) fresh[pos] = [];
      model.webLogoSelection = fresh;
      df.selection.setAll(false);
      const partials: {position: string; monomer: string; count: number}[] = [];
      for (const pos of Object.keys(mps).filter((k) => k !== 'general')) {
        const monomers = Object.keys(mps[pos]).filter((m) => m !== 'general');
        for (const mon of monomers) {
          const cnt = mps[pos][mon]?.count ?? 0;
          if (cnt > 0 && cnt < df.rowCount) { partials.push({position: pos, monomer: mon, count: cnt}); break; }
        }
        if (partials.length >= 2) break;
      }
      if (partials.length < 2) return {twoPicksFound: false};
      const [first, second] = partials;

      // Step 1 (re-establish Scenario-1 single pick).
      model.modifyWebLogoSelection(
        {positionOrClusterType: first.position, monomerOrCluster: first.monomer},
        {shiftPressed: false, ctrlPressed: false}, true);
      model.fireBitsetChanged('WebLogo');
      await new Promise((r) => setTimeout(r, 800));
      const selSingle = df.selection.trueCount;

      // Step 2: Shift-click a second distinct (monomer, position) — additive (union).
      let shiftThrew: string | null = null;
      try {
        model.modifyWebLogoSelection(
          {positionOrClusterType: second.position, monomerOrCluster: second.monomer},
          {shiftPressed: true, ctrlPressed: false}, true);
        model.fireBitsetChanged('WebLogo');
      } catch (e) { shiftThrew = String(e); }
      await new Promise((r) => setTimeout(r, 1000));
      const selAdditive = df.selection.trueCount;
      const mapAfterShift = {
        first: (model.webLogoSelection[first.position] ?? []).slice(),
        second: (model.webLogoSelection[second.position] ?? []).slice(),
      };
      // Step 4: SVM reflects the additive selection (viewer survives the re-broadcast).
      const svmAfterShift = document.querySelector('[name="viewer-Sequence-Variability-Map"]');
      const svmAliveAfterShift = svmAfterShift ? !!svmAfterShift.querySelector('canvas') : false;
      // Steps 5-6: widget surfaces re-render against the new selection (cold-stable predicate).
      const widgetsAfterShift: Record<string, boolean> = {};
      for (const paneName of ['pane-Distribution', 'pane-Selection']) {
        widgetsAfterShift[paneName] = (await (window as any).__paneHasContent(paneName)).hasContent;
      }
      const lastErrorAfterShift = grok.shell.lastError ? String(grok.shell.lastError) : '';

      // Step 7: Ctrl-click the second pick — toggle it off; re-broadcast the reduced selection.
      let ctrlThrew: string | null = null;
      try {
        model.modifyWebLogoSelection(
          {positionOrClusterType: second.position, monomerOrCluster: second.monomer},
          {shiftPressed: false, ctrlPressed: true}, true);
        model.fireBitsetChanged('WebLogo');
      } catch (e) { ctrlThrew = String(e); }
      await new Promise((r) => setTimeout(r, 1000));
      const selReduced = df.selection.trueCount;
      const mapAfterCtrl = {
        first: (model.webLogoSelection[first.position] ?? []).slice(),
        second: (model.webLogoSelection[second.position] ?? []).slice(),
      };
      const combinedReduced = model.getCombinedSelection()?.trueCount ?? null;
      // Step 8: cross-surfaces reflect the reduced selection.
      const svmAfterCtrl = document.querySelector('[name="viewer-Sequence-Variability-Map"]');
      const svmAliveAfterCtrl = svmAfterCtrl ? !!svmAfterCtrl.querySelector('canvas') : false;
      const widgetsAfterCtrl: Record<string, boolean> = {};
      for (const paneName of ['pane-Distribution', 'pane-Selection']) {
        widgetsAfterCtrl[paneName] = (await (window as any).__paneHasContent(paneName)).hasContent;
      }
      const lastErrorAfterCtrl = grok.shell.lastError ? String(grok.shell.lastError) : '';

      return {
        twoPicksFound: true, first, second,
        selSingle, selAdditive, selReduced, combinedReduced,
        mapAfterShift, mapAfterCtrl,
        shiftThrew, ctrlThrew,
        svmAliveAfterShift, svmAliveAfterCtrl,
        widgetsAfterShift, widgetsAfterCtrl,
        lastErrorAfterShift, lastErrorAfterCtrl,
      };
    });

    expect(result.twoPicksFound, 'need two distinct partial-count WebLogo picks for the modifier flow').toBe(true);
    // Step 2: Shift-add did not throw; additive semantics — union grows beyond the single pick
    // (additive, not replacement).
    expect(result.shiftThrew, `Shift-add WebLogo handler threw: ${result.shiftThrew}`).toBeNull();
    expect(result.selAdditive, 'Shift-click did not grow the selection (additive union expected, not replacement)')
      .toBeGreaterThan(result.selSingle!);
    // The unified Selection map carries BOTH picks after the Shift-add.
    expect(result.mapAfterShift!.first, 'first pick lost from the Selection map after Shift-add')
      .toContain(result.first!.monomer);
    expect(result.mapAfterShift!.second, 'second pick not recorded in the Selection map after Shift-add')
      .toContain(result.second!.monomer);
    // Step 4: SVM survives the re-broadcast. Steps 5-6: widget surfaces re-render.
    expect(result.svmAliveAfterShift, 'Sequence Variability Map lost its canvas after the Shift re-broadcast').toBe(true);
    expect(result.widgetsAfterShift!['pane-Distribution'],
      'Distribution widget did not re-render after the Shift re-broadcast').toBe(true);
    expect(result.widgetsAfterShift!['pane-Selection'],
      'Selection widget did not re-render after the Shift re-broadcast').toBe(true);
    // No null-receiver crash on the Shift re-broadcast.
    expect(isNullReceiverCrash(result.lastErrorAfterShift!),
      `GROK-14298 invariant: the Shift re-broadcast produced a null-receiver crash: ${result.lastErrorAfterShift}`)
      .toBe(false);

    // Step 7: Ctrl-toggle-off did not throw; selection returns to the single-pick state.
    expect(result.ctrlThrew, `Ctrl-toggle WebLogo handler threw: ${result.ctrlThrew}`).toBeNull();
    expect(result.selReduced, 'Ctrl-toggle-off did not return the selection to the single-pick count')
      .toBe(result.selSingle);
    // get-selection-bitset projection matches the reduced selection.
    expect(result.combinedReduced, 'getCombinedSelection did not match the reduced df.selection')
      .toBe(result.selReduced);
    // The second pick is removed from the unified Selection map; the first remains.
    expect(result.mapAfterCtrl!.second, 'Ctrl-toggle did not remove the second pick from the Selection map')
      .not.toContain(result.second!.monomer);
    expect(result.mapAfterCtrl!.first, 'Ctrl-toggle erroneously dropped the first (untouched) pick')
      .toContain(result.first!.monomer);
    // Step 8: cross-surfaces reflect the reduced selection.
    expect(result.svmAliveAfterCtrl, 'Sequence Variability Map lost its canvas after the Ctrl re-broadcast').toBe(true);
    expect(result.widgetsAfterCtrl!['pane-Distribution'],
      'Distribution widget did not re-render after the Ctrl re-broadcast').toBe(true);
    expect(result.widgetsAfterCtrl!['pane-Selection'],
      'Selection widget did not re-render after the Ctrl re-broadcast').toBe(true);
    // No null-receiver crash on the Ctrl re-broadcast (the listener-mutation crash mode).
    expect(isNullReceiverCrash(result.lastErrorAfterCtrl!),
      `GROK-14298 invariant: the Ctrl re-broadcast produced a null-receiver crash: ${result.lastErrorAfterCtrl}`)
      .toBe(false);
  });

  // Cleanup.
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
