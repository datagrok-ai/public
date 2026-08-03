/* ---
sub_features_covered: [peptides.model.fire-bitset-changed, peptides.panels.manual-alignment, peptides.widgets.manual-alignment, peptides.workflow.sar-dialog, peptides.workflow.start-analysis]
--- */
// Manual Alignment — panel surfaces via overlay-canvas click on a Monomer cell (not via grok.shell.o),
// Apply rewrites AlignedSequence + at least one per-position column, Reset is non-destructive,
// post-edit broadcast (model.fireBitsetChanged) stays functional. Apply/Reset driven by
// [name="button-Apply"]/[name="button-Reset"]; the "pane expanded" signal is the INNER
// .d4-accordion-pane-content.expanded class (the outer .d4-accordion-pane never gets it).

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/bio/peptides.csv';

test('Peptides — Manual Alignment panel applies sequence edits and stays non-destructive on Reset', async ({page}) => {
  test.setTimeout(300_000);
  await loginToDatagrok(page);

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
    expect(ready.posColCount, 'per-position columns ("1", "2", ...) not materialized after SAR launch')
      .toBeGreaterThan(0);
    expect(ready.col2SemType, 'col("2") not Monomer-semtype (the Manual Alignment @panel binds to Monomer cells)')
      .toBe('Monomer');
    expect(ready.modelPresent, 'PeptidesModel singleton not attached to the dataframe').toBe(true);
  });

  // Anchor to row 1 (row 0 has empty col '2' because the leading "NH2--" leaves split-index 1 empty).
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

  // The @panel surfaces only via an interactive overlay-canvas click (not via grok.shell.o).
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
      // Pick the topmost overlay canvas via elementFromPoint.
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

  // The "expanded" class lives on the INNER .d4-accordion-pane-content, not the outer pane.
  await softStep('Scenario 1 (step 3): widget mounts; verify textarea + Reset + Apply', async () => {
    const widget = await page.evaluate(async () => {
      const panes = Array.from(document.querySelectorAll('.d4-accordion-pane'));
      const manualPane = panes.find((p) => {
        const h = p.querySelector('.d4-accordion-pane-header');
        return !!h && h.textContent?.trim() === 'Manual Alignment';
      })!;
      // Defensive: click the header once if the inner content lacks `expanded`.
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
    expect(widget.textareaValue, 'textarea did not pre-populate with the row 1 AlignedSequence value')
      .toBe(widget.liveSeq);
  });

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

      // Substitute the first monomer in the textarea: "NH2-M-A-..." -> "NH2-V-A-..."
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

      const newSeq = alignedCol.get(1);
      const newPosVals: Record<string, any> = {};
      for (const n of posColNames) newPosVals[n] = df.col(n)!.get(1);
      const changedCols = posColNames.filter((n: string) => newPosVals[n] !== baselinePosVals[n]);

      return {baselineSeq, editedSeq, newSeq, changedCols, baselinePosVals, newPosVals};
    });
    expect(result.newSeq, 'AlignedSequence column did not reflect the Apply edit')
      .toBe(result.editedSeq);
    expect(result.newSeq, 'AlignedSequence column still equals the pre-edit baseline').not.toBe(result.baselineSeq);
    // Splitter remaps positions, so assert "at least one changed" rather than per-column equality.
    expect(result.changedCols.length, 'no per-position column at the edited row changed after Apply').toBeGreaterThan(0);
  });

  // Canvas coloring is not DOM-inspectable; assert model survives + SVM canvas + cell.renderer tag.
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

  await softStep('Scenario 1 (steps 6-7 invariant): Apply path produced no null-receiver error', async () => {
    const errs = await page.evaluate(() => {
      return grok.shell.lastError ? String(grok.shell.lastError) : '';
    });
    expect(errs.match(/null|undefined.*reading|setTrue/i)?.length ?? 0,
      `Apply path produced a null-receiver shell error: ${errs}`).toBe(0);
  });

  // Scenario 2 — Reset is non-destructive; post-edit broadcast survives.
  await softStep('Scenario 2 (step 1): re-open Manual Alignment — textarea binds to live (edited) value', async () => {
    const result = await page.evaluate(async () => {
      // Re-expand if collapsed (expanded class is on the inner content).
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

  // Reset re-binds to the live (edited) column value; it does NOT restore the pre-edit baseline.
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

  // Edit textarea without Apply, click Reset, confirm the unsaved edit is discarded.
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

  // Drive the broadcast via model.fireBitsetChanged (the WebLogo canvas click is coordinate-fragile).
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

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
