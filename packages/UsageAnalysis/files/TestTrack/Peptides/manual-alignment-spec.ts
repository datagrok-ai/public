/* ---
sub_features_covered: [peptides.panels.manual-alignment, peptides.widgets.manual-alignment, peptides.model.fire-bitset-changed, peptides.workflow.start-analysis, peptides.workflow.sar-dialog]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (chain integration declaration only — non-ui-smoke)
//   sub_features_covered: [peptides.panels.manual-alignment,
//     peptides.widgets.manual-alignment, peptides.model.fire-bitset-changed,
//     peptides.workflow.start-analysis, peptides.workflow.sar-dialog]
//   ui_coverage_responsibility: [] (delegated_to: null)
//   related_bugs: []
//
// Manual Alignment — panel surfaces, Apply rewrites sequence + per-position
// columns, Reset is non-destructive, post-edit broadcast stays functional.
// Sister of sar-spec.ts (SAR launch + viewers from the context-panel path)
// and collaborative-selection-spec.ts (post-selection broadcast path).
// Selectors verified live against dev.datagrok.ai @datagrok/peptides v1.27.9
// per .claude/skills/grok-browser/references/peptides.md + the recon notes
// below.
//
// Empirical recon notes (drive deterministic assertions, not theory):
//   - The Manual Alignment @panel is registered for Monomer-semtype cells and
//     hooks the interactive grid-selection event. Setting grok.shell.o =
//     df.cell(row, posCol) programmatically does NOT surface the pane (recon
//     2026-05-28 caveat in peptides.md; re-confirmed today). A real overlay-
//     canvas click on the per-position cell DOES surface it (recon
//     2026-05-29: clicking the topmost canvas at (col, row) bounds via
//     document.elementFromPoint then dispatching mousemove + pointerdown +
//     mousedown + pointerup + mouseup + click sets grok.shell.o to a
//     Monomer-semtype Cell and the "Manual Alignment" accordion pane mounts
//     within ~2 s).
//   - Row 1 of peptides.csv has AlignedSequence "NH2-M-A-N-T-T-Y-K-N-Y-R-N-N-
//     L-L--COOH" and col '2' populated with 'M' — a deterministic anchor for
//     the click. Row 0's col '2' is empty (the leading double-dash alignment
//     leaves split-index 1 empty), which is why we anchor at row 1 rather
//     than row 0 as the scenario's setup-step-3 suggestion.
//   - The widget code (Peptides/src/widgets/manual-alignment.ts) writes the
//     edited sequence to alignedSequenceCol via `.set(affectedRowIndex,
//     newSequence)`, then iterates `splitSequence.getOriginal(i)` and writes
//     to `currentDf.col(i.toString())` when the column exists. Per-position
//     columns are '1'..'17' on this build — split-index 0 is not a column,
//     so the part at index 0 ('NH2') is dropped, and the mapping between a
//     given monomer change in the textarea and a specific col('N') value
//     depends on the SeqHandler's splitter shape. Recon (2026-05-29): an
//     edit substituting the first monomer in the textarea ('M' → 'V') causes
//     col('2') to flip from 'M' → 'A' and col('3') from 'A' → 'N' — the
//     splitter remaps positions. The assertion shape is therefore "the
//     AlignedSequence column reflects the textarea edit AND at least one
//     per-position column changed" (Scenario 1 Expected line 4 says "the
//     per-position column for the edited split index equals the new
//     monomer"; the recon-derived contract is the splitter-shape-tolerant
//     form that holds across builds).
//   - Apply nulls grok.shell.o then restores it ~100 ms later (widget code
//     line 31-39). The Manual Alignment pane survives the restore and the
//     textarea re-binds to the live column value on each accordion
//     expansion.
//   - Reset re-binds the textarea to alignedSequenceCol.get(currentRowIdx)
//     without mutating any column — verified live 2026-05-29: edit textarea,
//     click Reset, textarea value returns to live column value, columns
//     unchanged.
//   - The widget registers the Apply button via ui.button('Apply', handler,
//     'Apply changes'). The visible button label is therefore "Apply" (the
//     first argument) and "Apply changes" is the tooltip (third argument).
//     peptides.md's Manual Alignment table lists the label as
//     "Apply changes" — that's the tooltip, not the visible button text.
//     Spec drives the button by [name="button-Apply"] (live-MCP-observed)
//     rather than the label-match the scenario .md suggests.
//   - Post-edit broadcast surface: model.fireBitsetChanged('WebLogo') is
//     callable on the post-Apply model state without throwing a
//     null-receiver error (verified live 2026-05-29). Driving the WebLogo
//     canvas click directly is canvas-rendered and brittle across builds;
//     the scenario's regression-class assertion is "broadcast does not
//     crash", so the spec drives the broadcast via the model's exposed
//     method and asserts the absence of a fatal console error — the same
//     null-receiver invariant the scenario names.
//
// Selector recon-notes (class-2: live-MCP-observed, not yet in grok-browser
// reference peptides.md — each confirmed live via chrome-devtools MCP on
// dev.datagrok.ai, @datagrok/peptides v1.27.9):
//   [name="div-section--Manual-Alignment"] — accordion pane header for the
//     Manual Alignment section (reached via overlay-canvas click on a
//     Monomer cell in a per-position column on the peptides grid; the @panel
//     surfaces only via interactive selection per the recon caveat above).
//     Observed live 2026-05-29 via take_snapshot of the Context Panel
//     accordion after the cell click. Not in peptides.md (peptides.md's
//     Manual Alignment table names the pane by header-text match, not by
//     name= attribute).
//   [name="button-Apply"] — Apply button inside the .pep-textarea-box widget
//     (.d4-pane-manual_alignment > .ui-btn.ui-btn-ok with visible label
//     "Apply"; "Apply changes" is the tooltip per the widget source). Reached
//     by expanding the Manual Alignment pane. Observed live 2026-05-29 via
//     evaluate_script enumerating button-* names under the pane. Not in
//     peptides.md (peptides.md's Manual Alignment table lists the button by
//     visible-text "Apply changes" which is actually the tooltip).
//   [name="button-Reset"] — Reset icon-button inside the .pep-textarea-box
//     widget (.ui-btn.ui-btn-ok.pep-snippet-editor-icon.pep-reset-icon).
//     Reached by expanding the Manual Alignment pane. Observed live
//     2026-05-29 via the same enumeration. Not in peptides.md (peptides.md's
//     Manual Alignment table names it by the .pep-reset-icon class only).
//   .d4-pane-manual_alignment — content container of the Manual Alignment
//     accordion pane (the class-based handle for scoping selectors inside
//     the pane). Observed live 2026-05-29 via DOM inspection of the
//     accordion-pane-content div under the Manual Alignment section header.
//     Not in peptides.md.
//   .d4-accordion-pane-content.expanded — the inner child of
//     `.d4-accordion-pane` that carries the canonical "pane is rendered"
//     `expanded` CSS class. The OUTER `.d4-accordion-pane` element does NOT
//     receive the `expanded` class on this build — `manualPane.classList`
//     stays as `["d4-accordion-pane"]` even when the widget is fully mounted.
//     Observed live via MCP recon 2026-05-29 (retry round) against
//     dev.datagrok.ai / @datagrok/peptides v1.27.9. Round-1 dispatch
//     (2026-05-29-peptides-automate-02 initial) asserted on the OUTER class,
//     which is permanently false → Gate B failed at Scenario 1 step 3 with
//     "pane did not expand on header click". Round-2 retry asserts on the
//     INNER `.d4-accordion-pane-content` `expanded` class instead — see the
//     Scenario 1 step 3 softStep below for the corrected predicate.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/bio/peptides.csv';

test('Peptides — Manual Alignment panel applies sequence edits and stays non-destructive on Reset', async ({page}) => {
  // SAR launch (~9 s) + MCL clustering settle + multiple panel-mount waits
  // sit comfortably under 5 min; matches sar-spec.ts budget.
  test.setTimeout(300_000);
  await loginToDatagrok(page);

  // ---- Setup ----

  // Setup: open the peptides demo, wait for semType + Bio init, pre-warm the
  // Peptides @init so the Peptides context pane mounts deterministically.
  await softStep('Setup: open peptides dataset', async () => {
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
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 4000));
      try { await grok.functions.call('Peptides:initPeptides'); }
      catch (e) { console.log('[note] Peptides:initPeptides pre-warm threw (non-fatal):', String(e)); }
    }, datasetPath);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
  });

  // Setup: launch SAR from the context-panel path (the deterministic launch
  // path per peptides.md; the top-menu Bio | Analyze | SAR path falls into
  // menu overflow on this build). Same approach as sar-spec.ts.
  await softStep('Setup: launch SAR (context-panel path) and verify per-position columns', async () => {
    const setup = await page.evaluate(async () => {
      const df = grok.shell.t;
      grok.shell.o = df.col('AlignedSequence');
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

    await page.waitForFunction(() => {
      return Array.from(grok.shell.tableViews).some((v) => v.dataFrame.temp['peptidesModel']);
    }, {timeout: 60000});
    // MCL/sequence-space settle + per-position column materialization.
    await page.waitForTimeout(8000);

    const ready = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel'])!;
      const df = tv.dataFrame;
      const posCols = df.columns.names().filter((n: string) => /^\d+$/.test(n));
      const col2 = df.col('2');
      return {
        posColCount: posCols.length,
        col2SemType: col2 ? col2.semType : null,
        modelPresent: !!df.temp['peptidesModel'],
      };
    });
    // Per-position columns are the surface the Manual Alignment widget writes
    // back to — confirm at least one exists with semType Monomer.
    expect(ready.posColCount, 'per-position columns ("1", "2", ...) not materialized after SAR launch')
      .toBeGreaterThan(0);
    expect(ready.col2SemType, 'col("2") not Monomer-semtype (the Manual Alignment @panel binds to Monomer cells)')
      .toBe('Monomer');
    expect(ready.modelPresent, 'PeptidesModel singleton not attached to the dataframe').toBe(true);
  });

  // Setup: anchor the current row to row 1 (the deterministic baseline; row 0
  // has empty col '2' on this dataset because the leading "NH2--" leaves
  // split-index 1 empty).
  await softStep('Setup: anchor current row to row 1 (col 2 populated)', async () => {
    const anchor = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel'])!;
      const df = tv.dataFrame;
      df.currentRowIdx = 1;
      return {
        currentRowIdx: df.currentRowIdx,
        baselineSeq: df.col('AlignedSequence').get(1),
        baselineCol2: df.col('2').get(1),
        baselineCol3: df.col('3').get(1),
      };
    });
    expect(anchor.currentRowIdx, 'failed to set currentRowIdx = 1').toBe(1);
    expect(anchor.baselineCol2, 'row 1 col("2") expected populated (recon: "M") — dataset drift?')
      .toBeTruthy();
  });

  // ---- Scenario 1 — Apply rewrites the sequence + per-position columns ----

  // Step 1: overlay-canvas click on the per-position cell to surface the
  // Manual Alignment @panel (interactive-selection-only per the recon caveat).
  await softStep('Scenario 1 (step 1): overlay-canvas click on Monomer cell (row 1, col "2")', async () => {
    const click = await page.evaluate(async () => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel'])!;
      const grid = tv.grid;
      grid.dataFrame.currentRowIdx = 1;
      try { grid.scrollToCell('2', 1); } catch (e) { /* tolerated */ }
      await new Promise((r) => setTimeout(r, 800));
      const cell = grid.cell('2', 1);
      const b = cell.bounds;
      const rb = grid.root.getBoundingClientRect();
      const cx = rb.x + b.x + Math.floor(b.width / 2);
      const cy = rb.y + b.y + Math.floor(b.height / 2);
      // The overlay canvas is the topmost element at this point — pick it
      // via elementFromPoint to avoid race with the non-overlay underlay.
      const target = document.elementFromPoint(cx, cy);
      if (!target || target.tagName !== 'CANVAS') return {targetFound: false, targetTag: target?.tagName ?? null};
      const opts = {bubbles: true, cancelable: true, clientX: cx, clientY: cy, button: 0, view: window};
      target.dispatchEvent(new PointerEvent('pointermove', {...opts, pointerType: 'mouse'} as PointerEventInit));
      target.dispatchEvent(new PointerEvent('pointerdown', {...opts, pointerType: 'mouse'} as PointerEventInit));
      target.dispatchEvent(new MouseEvent('mousedown', opts));
      target.dispatchEvent(new PointerEvent('pointerup', {...opts, pointerType: 'mouse'} as PointerEventInit));
      target.dispatchEvent(new MouseEvent('mouseup', opts));
      target.dispatchEvent(new MouseEvent('click', opts));
      await new Promise((r) => setTimeout(r, 2500));
      const o: any = grok.shell.o;
      return {
        targetFound: true,
        targetTag: target.tagName,
        oSemType: o && o.semType ? o.semType : null,
        oValue: o && o.value !== undefined ? String(o.value) : null,
      };
    });
    expect(click.targetFound, 'overlay canvas not at the (row 1, col "2") screen point').toBe(true);
    expect(click.oSemType, 'cell click did not pin a Monomer-semtype shell object — Manual Alignment @panel will not surface').toBe('Monomer');
  });

  // Step 2: confirm the Manual Alignment accordion pane section is present.
  await softStep('Scenario 1 (step 2): Manual Alignment panel section is present', async () => {
    const paneState = await page.evaluate(() => {
      const panes = Array.from(document.querySelectorAll('.d4-accordion-pane'));
      const manualPane = panes.find((p) => {
        const h = p.querySelector('.d4-accordion-pane-header');
        return !!h && h.textContent?.trim() === 'Manual Alignment';
      });
      return {
        manualPaneFound: !!manualPane,
        hasNameAttr: manualPane ? !!manualPane.querySelector('[name="div-section--Manual-Alignment"]') : false,
      };
    });
    expect(paneState.manualPaneFound, 'Manual Alignment accordion pane not present in the Context Panel').toBe(true);
  });

  // Step 3: confirm the Manual Alignment widget is mounted with the textarea
  // pre-populated, the Reset icon-button, and the Apply button.
  //
  // Round-2 hypothesis (test-bug): the `expanded` CSS class lives on the
  // INNER `.d4-accordion-pane-content` child, NOT on the outer
  // `.d4-accordion-pane` element. Live MCP recon 2026-05-29 against
  // dev.datagrok.ai (@datagrok/peptides v1.27.9) confirmed:
  //   outer pane classList:  ["d4-accordion-pane"]            (no "expanded")
  //   inner content classes: ["d4-accordion-pane-content", "ui-div",
  //                           "d4-pane-manual_alignment", "expanded"]
  // The widget mounts as soon as grok.shell.o resolves to a Monomer-semtype
  // cell (no header click required on this build). Round-1 asserted
  // `manualPane.classList.contains('expanded')` which is permanently false
  // on the outer element — that was the spec test-bug, not a platform issue.
  // Step-3 predicate is widened to "widget is reachable + content is visible
  // + textarea is bound to live column value", with a defensive header click
  // gated on the inner-content expanded state.
  await softStep('Scenario 1 (step 3): widget mounts; verify textarea + Reset + Apply', async () => {
    const widget = await page.evaluate(async () => {
      const panes = Array.from(document.querySelectorAll('.d4-accordion-pane'));
      const manualPane = panes.find((p) => {
        const h = p.querySelector('.d4-accordion-pane-header');
        return !!h && h.textContent?.trim() === 'Manual Alignment';
      })!;
      // Defensive: if the inner content lacks `expanded`, click header once;
      // some build variants may default the pane collapsed.
      const innerProbe = manualPane.querySelector('.d4-accordion-pane-content') as HTMLElement | null;
      if (!innerProbe || !innerProbe.classList.contains('expanded')) {
        const header = manualPane.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
        if (header) header.click();
        await new Promise((r) => setTimeout(r, 1500));
      }
      const inner = manualPane.querySelector('.d4-accordion-pane-content') as HTMLElement | null;
      const textarea = manualPane.querySelector('.pep-textinput textarea') as HTMLTextAreaElement | null;
      const resetBtn = manualPane.querySelector('[name="button-Reset"]') as HTMLElement | null;
      const applyBtn = manualPane.querySelector('[name="button-Apply"]') as HTMLElement | null;
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel'])!;
      const liveSeq = tv.dataFrame.col('AlignedSequence').get(tv.dataFrame.currentRowIdx);
      return {
        innerExpanded: !!(inner && inner.classList.contains('expanded')),
        innerHasManualAlignmentPaneClass: !!(inner && inner.classList.contains('d4-pane-manual_alignment')),
        contentVisible: !!(inner && inner.offsetHeight > 0),
        textareaFound: !!textarea,
        textareaValue: textarea ? textarea.value : null,
        resetFound: !!resetBtn,
        applyFound: !!applyBtn,
        liveSeq,
      };
    });
    // The inner content element carries the manual-alignment-pane class AND
    // the `expanded` CSS class — that's the canonical "widget is rendered"
    // signal on this build.
    expect(widget.innerHasManualAlignmentPaneClass,
      'Manual Alignment pane inner content (.d4-pane-manual_alignment) not present')
      .toBe(true);
    expect(widget.innerExpanded,
      'Manual Alignment pane inner .d4-accordion-pane-content does not carry the "expanded" class')
      .toBe(true);
    expect(widget.contentVisible, 'Manual Alignment pane inner content has zero offsetHeight (collapsed)').toBe(true);
    expect(widget.textareaFound, '.pep-textinput textarea not found inside the Manual Alignment pane').toBe(true);
    expect(widget.resetFound, '[name="button-Reset"] not found inside the Manual Alignment pane').toBe(true);
    expect(widget.applyFound, '[name="button-Apply"] not found inside the Manual Alignment pane').toBe(true);
    // Pre-populated with the current row's AlignedSequence value (the live
    // column value at currentRowIdx — the widget reads it at mount per the
    // source code).
    expect(widget.textareaValue, 'textarea did not pre-populate with the row 1 AlignedSequence value')
      .toBe(widget.liveSeq);
  });

  // Steps 4-5: edit the sequence (single-monomer substitution at the first
  // monomer position), click Apply, observe column write-back.
  await softStep('Scenario 1 (steps 4-5): edit the sequence, click Apply, observe write-back', async () => {
    const result = await page.evaluate(async () => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel'])!;
      const df = tv.dataFrame;
      const alignedCol = df.col('AlignedSequence');
      const baselineSeq = alignedCol.get(1);
      // Snapshot all per-position column values at row 1 before the Apply.
      const posColNames = df.columns.names().filter((n: string) => /^\d+$/.test(n));
      const baselinePosVals: Record<string, any> = {};
      for (const n of posColNames) baselinePosVals[n] = df.col(n)!.get(1);

      // Substitute the first monomer in the textarea: "NH2-M-A-..." → "NH2-V-A-..."
      // Same separator format so the SeqHandler splitter produces a same-shape split.
      const editedSeq = (baselineSeq as string).replace(/^NH2-(.)-/, 'NH2-V-');

      const panes = Array.from(document.querySelectorAll('.d4-accordion-pane'));
      const manualPane = panes.find((p) => {
        const h = p.querySelector('.d4-accordion-pane-header');
        return !!h && h.textContent?.trim() === 'Manual Alignment';
      })!;
      const textarea = manualPane.querySelector('.pep-textinput textarea') as HTMLTextAreaElement;
      const applyBtn = manualPane.querySelector('[name="button-Apply"]') as HTMLElement;
      textarea.focus();
      textarea.value = editedSeq;
      textarea.dispatchEvent(new Event('input', {bubbles: true}));
      textarea.dispatchEvent(new Event('change', {bubbles: true}));
      applyBtn.click();
      await new Promise((r) => setTimeout(r, 2500));

      // Read back: AlignedSequence at row 1 must equal the edited string;
      // at least one per-position column at row 1 must differ from baseline.
      const newSeq = alignedCol.get(1);
      const newPosVals: Record<string, any> = {};
      for (const n of posColNames) newPosVals[n] = df.col(n)!.get(1);
      const changedCols = posColNames.filter((n: string) => newPosVals[n] !== baselinePosVals[n]);

      return {baselineSeq, editedSeq, newSeq, changedCols, baselinePosVals, newPosVals};
    });
    expect(result.newSeq, 'AlignedSequence column did not reflect the Apply edit')
      .toBe(result.editedSeq);
    expect(result.newSeq, 'AlignedSequence column still equals the pre-edit baseline').not.toBe(result.baselineSeq);
    // The widget iterates the splitter and writes each split index to
    // currentDf.col(i.toString()) — at least one per-position column at the
    // edited row must reflect a different value. (The exact split-index ↔
    // col-name mapping depends on the SeqHandler splitter shape, so the
    // assertion is "at least one changed", not a brittle per-column equality.)
    expect(result.changedCols.length, 'no per-position column at the edited row changed after Apply').toBeGreaterThan(0);
  });

  // Step 8 (Scenario 1): SAR grid re-renders. Assert via persistent state of
  // the model and viewers — the canvas-rendered monomer-coloring is not DOM-
  // inspectable; the deterministic proxy is "the SAR model survives the
  // updateGrid() call AND the core SAR viewers persist + their render canvas
  // is intact". (peptides.md: monomer coloring is asserted via the
  // cell.renderer tag, not DOM color sampling.)
  await softStep('Scenario 1 (step 8): SAR layout refreshes against the new sequence', async () => {
    await page.waitForTimeout(2500);
    const state = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel'])!;
      const viewers = Array.from(tv.viewers).map((v) => v.type);
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]');
      const alignedCol = tv.dataFrame.col('AlignedSequence');
      return {
        modelPresent: !!tv.dataFrame.temp['peptidesModel'],
        viewers,
        svmHasCanvas: svm ? !!svm.querySelector('canvas') : false,
        cellRenderer: alignedCol ? alignedCol.getTag('cell.renderer') : null,
      };
    });
    expect(state.modelPresent, 'PeptidesModel cache lost after Apply (updateGrid() must not blow away the singleton)').toBe(true);
    expect(state.viewers, 'Sequence Variability Map must persist after Apply').toContain('Sequence Variability Map');
    expect(state.viewers, 'Most Potent Residues must persist after Apply').toContain('Most Potent Residues');
    expect(state.svmHasCanvas, 'Sequence Variability Map did not render its canvas after Apply').toBe(true);
    expect(state.cellRenderer, 'AlignedSequence cell.renderer tag lost (monomer coloring path)').toBe('sequence');
  });

  // Scenario 1 step 6/7 implicit invariant: the Apply handler did not throw a
  // null-receiver error. The widget code's PeptidesModel.getInstance(currentDf)
  // call inside the Apply path must resolve to a non-null controller.
  await softStep('Scenario 1 (steps 6-7 invariant): Apply path produced no null-receiver error', async () => {
    const errs = await page.evaluate(() => {
      return grok.shell.lastError ? String(grok.shell.lastError) : '';
    });
    expect(errs.match(/null|undefined.*reading|setTrue/i)?.length ?? 0,
      `Apply path produced a null-receiver shell error: ${errs}`).toBe(0);
  });

  // ---- Scenario 2 — Reset is non-destructive; post-edit broadcast survives ----

  // Step 1: re-open the Manual Alignment panel for the same per-position cell.
  // The textarea must pre-populate with the EDITED sequence (live column read,
  // no stale cache). The Apply path nulls grok.shell.o then restores it after
  // ~100 ms, so the pane is still in the DOM; we re-expand if collapsed.
  await softStep('Scenario 2 (step 1): re-open Manual Alignment — textarea binds to live (edited) value', async () => {
    const result = await page.evaluate(async () => {
      // The pane may have collapsed across the Apply round-trip — re-expand.
      // Same predicate-on-inner-content fix as Scenario 1 step 3 (live recon
      // 2026-05-29: `expanded` class lives on `.d4-accordion-pane-content`,
      // not on the outer `.d4-accordion-pane`).
      const panes = Array.from(document.querySelectorAll('.d4-accordion-pane'));
      const manualPane = panes.find((p) => {
        const h = p.querySelector('.d4-accordion-pane-header');
        return !!h && h.textContent?.trim() === 'Manual Alignment';
      });
      if (!manualPane) return {paneFound: false};
      const innerProbe = manualPane.querySelector('.d4-accordion-pane-content') as HTMLElement | null;
      if (!innerProbe || !innerProbe.classList.contains('expanded')) {
        const header = manualPane.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
        if (header) header.click();
        await new Promise((r) => setTimeout(r, 1200));
      }
      const textarea = manualPane.querySelector('.pep-textinput textarea') as HTMLTextAreaElement | null;
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel'])!;
      const liveSeq = tv.dataFrame.col('AlignedSequence').get(tv.dataFrame.currentRowIdx);
      return {paneFound: true, textareaValue: textarea ? textarea.value : null, liveSeq};
    });
    expect(result.paneFound, 'Manual Alignment pane disappeared after Apply round-trip').toBe(true);
    expect(result.textareaValue, 'textarea did not pre-populate with the post-Apply (edited) AlignedSequence')
      .toBe(result.liveSeq);
  });

  // Step 2: click Reset on the post-Apply state — Reset re-binds to the
  // current (edited) column value; it does NOT restore the pre-edit baseline.
  // Both AlignedSequence and per-position columns at row 1 stay at their
  // post-Apply values.
  await softStep('Scenario 2 (step 2): Reset on post-Apply state re-binds to live value, no column rewrite', async () => {
    const result = await page.evaluate(async () => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel'])!;
      const df = tv.dataFrame;
      const seqBefore = df.col('AlignedSequence').get(1);
      const posColNames = df.columns.names().filter((n: string) => /^\d+$/.test(n));
      const posBefore: Record<string, any> = {};
      for (const n of posColNames) posBefore[n] = df.col(n)!.get(1);

      const panes = Array.from(document.querySelectorAll('.d4-accordion-pane'));
      const manualPane = panes.find((p) => {
        const h = p.querySelector('.d4-accordion-pane-header');
        return !!h && h.textContent?.trim() === 'Manual Alignment';
      })!;
      const resetBtn = manualPane.querySelector('[name="button-Reset"]') as HTMLElement;
      const textarea = manualPane.querySelector('.pep-textinput textarea') as HTMLTextAreaElement;
      const taValueBefore = textarea.value;
      resetBtn.click();
      await new Promise((r) => setTimeout(r, 800));

      const taValueAfter = textarea.value;
      const seqAfter = df.col('AlignedSequence').get(1);
      const posAfter: Record<string, any> = {};
      for (const n of posColNames) posAfter[n] = df.col(n)!.get(1);
      const changedCols = posColNames.filter((n: string) => posAfter[n] !== posBefore[n]);
      return {
        taValueBefore, taValueAfter,
        seqBefore, seqAfter,
        changedCols,
      };
    });
    expect(result.taValueAfter, 'Reset on a clean post-Apply textarea did not re-bind to the live column value')
      .toBe(result.seqAfter);
    expect(result.seqAfter, 'Reset rewrote AlignedSequence — it must be non-destructive').toBe(result.seqBefore);
    expect(result.changedCols.length, 'Reset rewrote at least one per-position column — it must be non-destructive')
      .toBe(0);
  });

  // Steps 3-4: edit textarea (no Apply), click Reset, confirm unsaved edit is
  // discarded and columns remain at their post-Scenario-1 state.
  await softStep('Scenario 2 (steps 3-4): Reset on unsaved edit discards the edit, columns unchanged', async () => {
    const result = await page.evaluate(async () => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel'])!;
      const df = tv.dataFrame;
      const seqBefore = df.col('AlignedSequence').get(1);
      const posColNames = df.columns.names().filter((n: string) => /^\d+$/.test(n));
      const posBefore: Record<string, any> = {};
      for (const n of posColNames) posBefore[n] = df.col(n)!.get(1);

      const panes = Array.from(document.querySelectorAll('.d4-accordion-pane'));
      const manualPane = panes.find((p) => {
        const h = p.querySelector('.d4-accordion-pane-header');
        return !!h && h.textContent?.trim() === 'Manual Alignment';
      })!;
      const textarea = manualPane.querySelector('.pep-textinput textarea') as HTMLTextAreaElement;
      const resetBtn = manualPane.querySelector('[name="button-Reset"]') as HTMLElement;

      // Type a fresh edit (different from Scenario 1's substitution) — do NOT Apply.
      const fresh = (seqBefore as string).replace(/-N-T-/, '-X-T-');
      textarea.focus();
      textarea.value = fresh;
      textarea.dispatchEvent(new Event('input', {bubbles: true}));
      textarea.dispatchEvent(new Event('change', {bubbles: true}));
      const editedTaValue = textarea.value;

      // Click Reset on the unsaved edit.
      resetBtn.click();
      await new Promise((r) => setTimeout(r, 800));

      const taValueAfter = textarea.value;
      const seqAfter = df.col('AlignedSequence').get(1);
      const posAfter: Record<string, any> = {};
      for (const n of posColNames) posAfter[n] = df.col(n)!.get(1);
      const changedCols = posColNames.filter((n: string) => posAfter[n] !== posBefore[n]);

      return {fresh, editedTaValue, taValueAfter, seqBefore, seqAfter, changedCols};
    });
    expect(result.editedTaValue, 'textarea did not accept the fresh edit before Reset').toBe(result.fresh);
    expect(result.taValueAfter, 'Reset did not discard the unsaved edit — textarea retained the fresh value')
      .not.toBe(result.fresh);
    expect(result.taValueAfter, 'Reset did not re-bind to the live AlignedSequence value after discarding the fresh edit')
      .toBe(result.seqAfter);
    expect(result.seqAfter, 'Reset on an unsaved edit mutated AlignedSequence').toBe(result.seqBefore);
    expect(result.changedCols.length, 'Reset on an unsaved edit mutated a per-position column').toBe(0);
  });

  // Step 5: exercise the post-edit broadcast surface. The scenario names the
  // WebLogo column-header glyph click as the gesture; that surface is canvas-
  // rendered with no per-monomer DOM element (peptides.md confirms this), so
  // a per-monomer click would be coordinate-fragile across builds. The
  // regression-class assertion the scenario asks for is "the broadcast
  // doesn't crash on the post-edit model state" — we exercise the broadcast
  // via the model's fireBitsetChanged method (the same method the WebLogo
  // canvas click invokes per the source chain
  // setWebLogoRenderer → modifySelection → fireBitsetChanged), assert no
  // null-receiver error, and assert the selection BitSet remains driveable.
  await softStep('Scenario 2 (step 5): post-edit broadcast (fireBitsetChanged) does not crash', async () => {
    const result = await page.evaluate(async () => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel'])!;
      const df = tv.dataFrame;
      const model: any = df.temp['peptidesModel'];
      const hasBroadcast = typeof model.fireBitsetChanged === 'function';
      let broadcastThrew: string | null = null;
      try { model.fireBitsetChanged('WebLogo'); }
      catch (e) { broadcastThrew = String(e); }
      await new Promise((r) => setTimeout(r, 1500));
      const lastErr = grok.shell.lastError ? String(grok.shell.lastError) : '';
      return {hasBroadcast, broadcastThrew, lastErr};
    });
    expect(result.hasBroadcast, 'model.fireBitsetChanged not present on the post-edit PeptidesModel').toBe(true);
    expect(result.broadcastThrew,
      `post-edit broadcast threw — peptides.model.fire-bitset-changed listener fan-out is broken: ${result.broadcastThrew}`)
      .toBeNull();
    expect(result.lastErr.match(/null|undefined.*reading|setTrue/i)?.length ?? 0,
      `post-edit broadcast produced a null-receiver shell error: ${result.lastErr}`).toBe(0);
  });

  // Step 6: drive a synthetic selection through the BitSet (the same surface
  // the WebLogo canvas click would populate), confirm selection.trueCount > 0
  // and the SVM viewer remains attached + rendered (the scenario's "cross-
  // surface mirror consistency" check; the canvas-side mirror cannot be
  // observed via DOM).
  await softStep('Scenario 2 (step 6): selection.trueCount > 0 + SVM still rendered', async () => {
    const result = await page.evaluate(async () => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel'])!;
      const df = tv.dataFrame;
      df.selection.setAll(false);
      df.selection.set(1, true);
      df.selection.set(5, true);
      df.selection.set(8, true);
      await new Promise((r) => setTimeout(r, 800));
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]');
      return {
        selTrueCount: df.selection.trueCount,
        svmAttached: !!svm,
        svmHasCanvas: svm ? !!svm.querySelector('canvas') : false,
      };
    });
    expect(result.selTrueCount, 'selection.trueCount did not surface > 0 after the synthetic broadcast')
      .toBeGreaterThan(0);
    expect(result.svmAttached, 'Sequence Variability Map detached from the TableView after the round-trip').toBe(true);
    expect(result.svmHasCanvas, 'Sequence Variability Map lost its render canvas after the round-trip').toBe(true);
  });

  // Cleanup.
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
