/* ---
sub_features_covered: [peptides.workflow.start-analysis, peptides.util.get-selection-bitset, peptides.compute.calculate-monomer-position-statistics, peptides.rendering.weblogo-header, peptides.workflow.sar-dialog]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (coverage_type: edge — atlas edge_cases[0] anchor)
//   sub_features_covered: [peptides.workflow.start-analysis,
//     peptides.util.get-selection-bitset,
//     peptides.compute.calculate-monomer-position-statistics,
//     peptides.rendering.weblogo-header, peptides.workflow.sar-dialog]
//   ui_coverage_responsibility: [] (delegated_to: null)
//   related_bugs: [GROK-19145]
//   produced_from: atlas-driven
//
// Atlas provenance (derived_from):
//   feature-atlas/peptides.yaml#edge_cases[0] derived_from:
//     bug-library/peptides.yaml#GROK-19145
//
// SAR Similarity-threshold matrix — the GROK-19145 edge anchor. Drives the
// top-menu Bio | Analyze | SAR... config dialog across representative Similarity
// values (10 / 50 / 75 / 90) and asserts the bug invariant: SAR completes
// gracefully OR surfaces an informative message, but NEVER throws a null-receiver
// runtime crash ("method not found: 'setTrue' on null"), and the WebLogo column
// headers render with populated stats rather than silently empty. Sister of the
// context-panel Launch SAR path (sar.md / sar-spec.ts) and the top-menu launch +
// MCL re-render path (peptide-space.md / peptide-space-spec.ts).
// Selectors verified against dev.datagrok.ai @datagrok/peptides v1.27.9 per
// .claude/skills/grok-browser/references/peptides.md.
//
// Empirical recon notes (live MCP recon 2026-05-29 — drive deterministic
// assertions, not theory):
//   - The Similarity threshold IS a field of the Analyze Peptides dialog
//     ([name="dialog-Analyze-Peptides"]) — [name="input-host-Similarity-Threshold"]
//     input, default "70". This is exactly the GROK-19145 repro surface
//     ("Bio > Analyze > SAR... -> set Similarity = 90"). The host's containing
//     MCL sub-form renders display:none (advanced/clustering knobs are collapsed
//     in the dialog), so the input is offsetParent-null; setting .value +
//     dispatching input/change still drives the Dart change listener and the
//     value propagates to model._settings.mclSettings.threshold on OK (verified
//     live for 10 and 90). This is the sanctioned hidden-input driving path —
//     not a JS-API substitution: the OK button is clicked through the DOM and the
//     dialog round-trips the value into the analysis.
//   - The top-menu path is deterministic at the spec viewport (1920x1080,
//     specTestOptions): [name="div-Bio"] is rendered + visible (offsetParent
//     non-null). The grok-browser overflow caveat (Bio in overflow) was a 1536px
//     CSS recon artifact; at 1920px the menubar shows
//     Edit | View | Select | Data | ML | Bio.
//   - GROK-19145 is FIXED on this build: at Similarity=90 SAR completes,
//     MonomerPositionStats is POPULATED (17 positions, 167 monomer entries — not
//     empty/sparse), WebLogo column headers render (grid.props.colHeaderHeight
//     grows to 130 and the header canvas strip has non-background content), and
//     grok.shell.lastError is a benign "[object Promise]", never a setTrue-on-null
//     TypeError. The spec asserts the graceful-behavior invariant; it would FAIL
//     loudly if the bug regressed (null-receiver error, or empty WebLogo with no
//     placeholder).
//   - Default attach set across all thresholds: Sequence Variability Map +
//     Most Potent Residues + MCL. Logo Summary Table is cluster-result-dependent
//     (attached at threshold=90, NOT at threshold=10) — recorded softly, not
//     asserted.
//   - The selection backbone (peptides.util.get-selection-bitset /
//     modify-selection — what the WebLogo header click wires to) survives sparse
//     stats: a synthetic MouseEvent on the Sequence Variability Map render canvas
//     at (rect.x+120, rect.y+120) raises df.selection.trueCount 0 -> 43,
//     deterministic, no crash. The WebLogo column-header glyph click itself is
//     positionally fragile (hitting an exact monomer pixel on the canvas header),
//     so the SVM cell click is used as the canonical, deterministic
//     selection-backbone observable (same gesture asserted in sar-spec.ts).
//
// Selector recon-notes (class-2: live-MCP-observed, not yet in grok-browser
// reference peptides.md). Each re-confirmed live 2026-05-29 via chrome-devtools
// MCP (evaluate_script lookups against dev.datagrok.ai, @datagrok/peptides
// v1.27.9, specTestOptions viewport):
//   [name="div-Bio"] — top menubar entry (TableView ribbon menu), offsetParent
//     non-null (visible) at 1920x1080. peptides.md documents this path only with
//     a 1536px-CSS overflow caveat ("Bio not in visible menubar"); at the spec
//     viewport it IS rendered + visible. Reached: open peptides.csv table → read
//     the visible menubar (Edit | View | Select | Data | ML | Bio).
//   [name="div-Bio---Analyze"] — Bio menu's Analyze submenu group; reached by
//     clicking [name="div-Bio"] then mouseenter/mousemove on the submenu.
//   [name="div-Bio---Analyze---SAR..."] — the SAR... menu item under
//     Bio | Analyze; clicking it opens the Analyze Peptides config dialog.
//   [name="dialog-Analyze-Peptides"] — the analyzePeptidesUI config dialog
//     (title "Analyze Peptides"), distinct from the wrench
//     [name="dialog-Peptides-settings"] in peptides.md. Opened by the SAR... menu
//     item above; carries [name="button-OK"].
//   [name="input-host-Similarity-Threshold"] input — the Similarity threshold
//     field in the Analyze-Peptides dialog (default value "70"). Its containing
//     MCL sub-form is collapsed (host offsetParent-null); value+input/change
//     dispatch still round-trips into model._settings.mclSettings.threshold on OK
//     (observed: set 90 → appliedThreshold 90). peptides.md lists
//     input-Similarity-Threshold only in the wrench-settings / context-panel MCL
//     pane, not in this top-menu Analyze dialog.
//   [name="viewer-Sequence-Variability-Map"] — present in peptides.md (class-1);
//     listed here only for the selection-backbone gesture context. Largest-canvas
//     grid-body MouseEvent click at (rect.x+120, rect.y+120) raised
//     df.selection.trueCount 0 -> 43, no throw, benign lastError ([object Promise]).

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/bio/peptides.csv';

// Open the peptides Macromolecule table fresh: close everything, read the CSV,
// wait for semType detection + Bio package settle, pre-warm the Peptides @init
// (GROK-17557 family — removes the cold-package init leg), open the Grid.
async function openPeptidesTable(page: import('@playwright/test').Page) {
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
  return result;
}

// Drive Bio | Analyze | SAR... -> Analyze Peptides dialog, set the Similarity
// threshold to `value`, click OK, wait for the PeptidesModel to attach + the
// compute pipeline to settle. Returns the live SAR outcome: applied threshold,
// attached viewers, per-position stats population (the WebLogo backing data),
// whether the WebLogo header region rendered, and any null-receiver error in the
// post-OK compute path (the GROK-19145 invariant).
async function launchSarWithSimilarity(page: import('@playwright/test').Page, value: number) {
  // Open the dialog via the top menu.
  const opened = await page.evaluate(async () => {
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
    return {bioFound: !!bio, bioVisible, analyzeFound: !!analyze, sarFound: !!sar, dialogFound: !!dlg};
  });
  expect(opened.bioFound, '[name="div-Bio"] top-menu entry not found').toBe(true);
  expect(opened.bioVisible,
    '[name="div-Bio"] present but not visible (offsetParent null — menu in overflow at this viewport)').toBe(true);
  expect(opened.analyzeFound, '[name="div-Bio---Analyze"] submenu entry not found').toBe(true);
  expect(opened.sarFound, '[name="div-Bio---Analyze---SAR..."] menu item not found').toBe(true);
  expect(opened.dialogFound, '[name="dialog-Analyze-Peptides"] config dialog did not open').toBe(true);

  // Set the Similarity threshold + click OK.
  const set = await page.evaluate(async (v) => {
    const dlg = document.querySelector('[name="dialog-Analyze-Peptides"]')!;
    const inner = dlg.querySelector('[name="input-host-Similarity-Threshold"] input') as HTMLInputElement | null;
    const inputFound = !!inner;
    if (inner) {
      // Drive the Dart change listener (host is offsetParent-null — collapsed
      // MCL sub-form — but value+event dispatch still propagates on OK).
      inner.focus();
      inner.value = String(v);
      inner.dispatchEvent(new Event('input', {bubbles: true}));
      inner.dispatchEvent(new Event('change', {bubbles: true}));
      inner.blur();
    }
    await new Promise((r) => setTimeout(r, 400));
    const setVal = inner ? inner.value : null;
    const ok = (dlg.querySelector('[name="button-OK"]')
      ?? document.querySelector('[name="button-OK"]')) as HTMLElement | null;
    if (ok) ok.click();
    return {inputFound, setVal};
  }, value);
  expect(set.inputFound, '[name="input-Similarity-Threshold"] not found in the Analyze Peptides dialog').toBe(true);
  expect(set.setVal, `Similarity input did not accept the value ${value}`).toBe(String(value));

  // The analysis completes when the PeptidesModel singleton attaches to a table
  // view's DataFrame. Poll for it rather than a fixed wait (GROK-19145 invariant:
  // the workflow must reach completion gracefully even at the extreme threshold).
  await page.waitForFunction(() => {
    return Array.from(grok.shell.tableViews).some((v) => v.dataFrame.temp['peptidesModel']);
  }, {timeout: 90000});
  // Let the MCL clustering + sequence-space embedding settle.
  await page.waitForTimeout(8000);

  return await page.evaluate(() => {
    const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
    const model = tv.dataFrame.temp['peptidesModel'] as any;
    const viewers = Array.from(tv.viewers).map((v) => v.type);
    const appliedThreshold = model?._settings?.mclSettings?.threshold ?? null;

    // The WebLogo backing data — per-position monomer stats. Empty/sparse stats
    // were the GROK-19145 empty-WebLogo failure mode; assert it is populated.
    const mps = model?.monomerPositionStats ?? {};
    let positionsWithStats = 0;
    let totalMonomerEntries = 0;
    for (const k of Object.keys(mps)) {
      if (k === 'general') continue;
      const posObj = mps[k];
      if (posObj && typeof posObj === 'object') {
        const monomers = Object.keys(posObj).filter((m) => m !== 'general');
        if (monomers.length > 0) positionsWithStats++;
        totalMonomerEntries += monomers.length;
      }
    }

    // WebLogo column-header rendering: after start-analysis the per-position
    // columns carry WebLogo headers and the grid's colHeaderHeight grows (~130px)
    // to fit them. Also sample the top header strip of the grid canvas for
    // non-background content (the drawn logo glyphs) — a silently-blank header
    // (no glyph, no placeholder) is the regression.
    const colHeaderHeight = tv.grid?.props?.colHeaderHeight ?? 0;
    let headerHasRender = false;
    const gridViewer = document.querySelector('[name="viewer-Grid"]');
    if (gridViewer) {
      const cs = Array.from(gridViewer.querySelectorAll('canvas')) as HTMLCanvasElement[];
      let best: HTMLCanvasElement | null = null, area = 0;
      for (const c of cs) { const r = c.getBoundingClientRect(); if (r.width * r.height > area) { area = r.width * r.height; best = c; } }
      if (best && best.width > 0 && best.height > 0) {
        try {
          const ctx = best.getContext('2d')!;
          const h = Math.min(130, best.height);
          const data = ctx.getImageData(0, 0, best.width, h).data;
          let nonBg = 0;
          for (let i = 0; i < data.length; i += 41)
            if (data[i + 3] > 0 && (data[i] < 250 || data[i + 1] < 250 || data[i + 2] < 250)) nonBg++;
          headerHasRender = nonBg > 50;
        } catch (e) { headerHasRender = false; }
      }
    }

    const lastError = grok.shell.lastError ? String(grok.shell.lastError) : '';
    return {
      appliedThreshold, viewers, positionsWithStats, totalMonomerEntries,
      colHeaderHeight, headerHasRender, lastError,
      modelPresent: !!model,
    };
  });
}

// The GROK-19145 invariant: the post-OK compute path must NOT throw a
// null-receiver runtime error (setTrue / getRawData on a null BitSet/column).
// Benign async/Promise + resource-404 noise is tolerated; only null-receiver
// crashes fail the invariant.
function assertNoNullReceiverCrash(lastError: string, threshold: number) {
  expect(/setTrue|method not found.*null|Cannot read propert.*null|reading 'getRawData'/.test(lastError),
    `GROK-19145 invariant: SAR at Similarity=${threshold} produced a null-receiver crash: ${lastError}`).toBe(false);
}

test('SAR Similarity-threshold matrix — graceful across low/medium/high/extreme (GROK-19145)', async ({page}) => {
  // Four parametrized SAR launches, each an async server-compute pipeline
  // (~9 s + MCL clustering), plus a settle phase — won't fit the default budget.
  test.setTimeout(600_000);
  await loginToDatagrok(page);

  await softStep('Setup: open the peptides Macromolecule table', async () => {
    const result = await openPeptidesTable(page);
    expect(result.rows, 'peptides.csv should load 647 rows').toBe(647);
    expect(result.semType, 'AlignedSequence must be a Macromolecule column').toBe('Macromolecule');
  });

  // ---- Scenario 1 — Similarity matrix: SAR completes without null-receiver
  //      crash across low / medium / high / extreme values ----

  // The four representative Similarity values from the scenario. Each launch
  // re-opens the dataset fresh + re-runs SAR via the top-menu dialog (the
  // scenario re-opens the dialog per value; a clean table per launch isolates the
  // per-value compute path and avoids cross-launch model-state carryover).
  for (const threshold of [10, 50, 75, 90]) {
    const label = threshold === 10 ? 'low' : threshold === 50 ? 'medium' : threshold === 75 ? 'high' : 'extreme';
    await softStep(`Scenario 1: SAR at Similarity=${threshold} (${label}) completes gracefully`, async () => {
      // Fresh table per launch (skip the very first — Setup already opened it for
      // threshold=10).
      if (threshold !== 10)
        await openPeptidesTable(page);

      const out = await launchSarWithSimilarity(page, threshold);

      // The dialog round-tripped the Similarity value into the analysis.
      expect(out.appliedThreshold,
        `Similarity=${threshold} did not propagate to the model mclSettings.threshold`).toBe(threshold);
      // The workflow reached completion (model attached) — not a frozen/blocked UI.
      expect(out.modelPresent, `PeptidesModel did not attach after SAR at Similarity=${threshold}`).toBe(true);

      // Either the SAR viewers attached AND the WebLogo headers rendered, OR an
      // informative platform message was surfaced. On this build the workflow
      // completes with rendered output (no informative-error branch needed), so
      // assert the rendered-output arm: deterministic viewers + populated stats +
      // drawn WebLogo headers.
      expect(out.viewers, `Sequence Variability Map must attach at Similarity=${threshold}`)
        .toContain('Sequence Variability Map');
      expect(out.viewers, `Most Potent Residues must attach at Similarity=${threshold}`)
        .toContain('Most Potent Residues');
      expect(out.viewers, `MCL clustering viewer must attach at Similarity=${threshold}`)
        .toContain('MCL');
      // Logo Summary Table is cluster-result-dependent (attaches at high
      // thresholds, not low) — record without failing the deterministic contract.
      if (!out.viewers.includes('Logo Summary Table'))
        console.log(`[note] Logo Summary Table not attached at Similarity=${threshold} (cluster-result-dependent).`);

      // WebLogo column-headers render with populated stats — NOT silently empty.
      // (The GROK-19145 failure mode is empty WebLogo from empty/sparse stats.)
      expect(out.positionsWithStats,
        `MonomerPositionStats is empty at Similarity=${threshold} (WebLogo would render blank)`)
        .toBeGreaterThan(0);
      expect(out.colHeaderHeight,
        `grid colHeaderHeight did not grow for the WebLogo headers at Similarity=${threshold}`)
        .toBeGreaterThan(40);
      expect(out.headerHasRender,
        `WebLogo column-header strip drew no content at Similarity=${threshold} (silently-blank header is a regression)`)
        .toBe(true);

      // GROK-19145 invariant: no null-receiver runtime crash in the post-OK path.
      assertNoNullReceiverCrash(out.lastError, threshold);
    });
  }

  // ---- Scenario 2 — Similarity=90 single launch preserves the WebLogo rendering
  //      contract + the selection backbone survives sparse stats ----

  await softStep('Scenario 2 (steps 1-2): launch SAR at Similarity=90 (extreme), settle', async () => {
    await openPeptidesTable(page);
    const out = await launchSarWithSimilarity(page, 90);
    expect(out.appliedThreshold, 'Similarity=90 did not propagate to the model').toBe(90);

    // Step 3: the Sequence Variability Map viewer is present (attached, not skipped).
    expect(out.viewers, 'Sequence Variability Map must be attached at Similarity=90').toContain('Sequence Variability Map');

    // Step 4: the per-position WebLogo headers render — populated, not silently
    // blank. (A header with no glyph, no placeholder, no error is the regression.)
    expect(out.positionsWithStats, 'WebLogo backing stats are empty at Similarity=90').toBeGreaterThan(0);
    expect(out.headerHasRender, 'WebLogo column-headers drew nothing at Similarity=90 (silently-blank regression)').toBe(true);

    // No null-receiver runtime error regardless of how sparse the stats are.
    assertNoNullReceiverCrash(out.lastError, 90);
  });

  // Step 5: clicking a populated SAR map cell invokes the selection backbone
  // (modifySelection -> get-selection-bitset -> fireBitsetChanged -> DataFrame
  // BitSet update) without crashing — the click handler remains live even when
  // the underlying stats are sparse at the extreme threshold. The WebLogo
  // column-header glyph click is canvas + positionally fragile; the Sequence
  // Variability Map render-canvas cell click is the canonical, deterministic
  // observable of the same selection backbone (verified live: selection 0 -> 43).
  await softStep('Scenario 2 (step 5): selection backbone responds without crashing', async () => {
    const result = await page.evaluate(async () => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const df = tv.dataFrame;
      const selBefore = df.selection.trueCount;
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]');
      if (!svm) return {svmFound: false, selBefore, selAfter: selBefore, threw: null, lastError: ''};
      // Largest canvas is the render surface (the tiny one is a scrollbar).
      const canvases = Array.from(svm.querySelectorAll('canvas')) as HTMLCanvasElement[];
      let canvas: HTMLCanvasElement | null = null, maxArea = 0;
      for (const c of canvases) { const r = c.getBoundingClientRect(); if (r.width * r.height > maxArea) { maxArea = r.width * r.height; canvas = c; } }
      if (!canvas) return {svmFound: true, canvasFound: false, selBefore, selAfter: selBefore, threw: null, lastError: ''};
      const r = canvas.getBoundingClientRect();
      // Deterministic data-cell position past the row-header column + header strip.
      const opts = {bubbles: true, cancelable: true, clientX: r.x + 120, clientY: r.y + 120, button: 0, view: window};
      let threw: string | null = null;
      try {
        canvas.dispatchEvent(new MouseEvent('mousemove', opts));
        canvas.dispatchEvent(new MouseEvent('mousedown', opts));
        canvas.dispatchEvent(new MouseEvent('mouseup', opts));
        canvas.dispatchEvent(new MouseEvent('click', opts));
      } catch (e) { threw = String(e); }
      await new Promise((res) => setTimeout(res, 2500));
      const selAfter = df.selection.trueCount;
      const lastError = grok.shell.lastError ? String(grok.shell.lastError) : '';
      return {svmFound: true, canvasFound: true, selBefore, selAfter, threw, lastError};
    });
    expect(result.svmFound, 'Sequence Variability Map viewer not found for the selection-backbone click').toBe(true);
    expect(result.canvasFound, 'SVM render canvas not found for the selection-backbone click').toBe(true);
    // The click handler did not throw.
    expect(result.threw, `SVM cell click threw: ${result.threw}`).toBeNull();
    // The selection backbone fired (BitSet updated) — survives sparse stats at 90.
    expect(result.selAfter, 'SVM cell click did not update the DataFrame selection (backbone did not fire)')
      .toBeGreaterThan(result.selBefore);
    // And it did so without a null-receiver crash.
    assertNoNullReceiverCrash(result.lastError, 90);
  });

  // Cleanup.
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
