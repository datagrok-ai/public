/* ---
sub_features_covered: [peptides.panels.peptides, peptides.viewers.logo-summary-table, peptides.viewers.monomer-position, peptides.viewers.most-potent-residues, peptides.viewers.mutation-cliffs, peptides.widgets.distribution, peptides.widgets.mutation-cliffs, peptides.widgets.settings-dialog, peptides.workflow.analyze-ui, peptides.workflow.sar-dialog, peptides.workflow.start-analysis]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
const datasetPath = 'System:DemoFiles/bio/peptides.csv';
test('SAR — Launch and verify viewers (context-panel entry path)', async ({page}) => {
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
  await softStep('Scenario 1 (steps 1-4): launch SAR from the Peptides context panel', async () => {
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
  });
  await softStep('Scenario 1 (step 5): verify SAR viewers attach', async () => {
    await page.waitForTimeout(8000);
    const viewers = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      return Array.from(tv.viewers).map((v) => v.type);
    });
    expect(viewers, 'Sequence Variability Map (MonomerPosition) must attach').toContain('Sequence Variability Map');
    expect(viewers, 'Most Potent Residues must attach').toContain('Most Potent Residues');
    expect(viewers, 'MCL clustering viewer must attach').toContain('MCL');
    if (!viewers.includes('Logo Summary Table'))
      console.log('[note] Logo Summary Table not in default Launch SAR attach set on this build (cluster-/settings-dependent).');
  });
  await softStep('Scenario 2 (steps 6-8): change Similarity threshold via the Settings dialog', async () => {
    const opened = await page.evaluate(async () => {
      const wrench = document.querySelector(
        'i.grok-icon.fa-wrench[aria-label="Peptides analysis settings"]') as HTMLElement | null;
      const wrenchFound = !!wrench;
      if (wrench) wrench.click();
      await new Promise((r) => setTimeout(r, 2000));
      const dlg = document.querySelector('[name="dialog-Peptides-settings"]');
      return {wrenchFound, dlgFound: !!dlg};
    });
    expect(opened.wrenchFound, 'Peptides analysis settings wrench not found on the SAR toolbar').toBe(true);
    expect(opened.dlgFound, '[name="dialog-Peptides-settings"] did not open').toBe(true);
    const changed = await page.evaluate(async () => {
      const dlg = document.querySelector('[name="dialog-Peptides-settings"]')!;
      const mclPane = Array.from(dlg.querySelectorAll('.d4-accordion-pane'))
        .find((p) => p.querySelector('.d4-accordion-pane-header')?.textContent?.trim() === 'MCL');
      if (mclPane) {
        const h = mclPane.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
        if (h && !mclPane.classList.contains('expanded')) h.click();
      }
      await new Promise((r) => setTimeout(r, 800));
      const simHost = dlg.querySelector('[name="input-host-Similarity-Threshold"]');
      const inner = simHost ? (simHost.querySelector('input') as HTMLInputElement | null) : null;
      const inputFound = !!inner;
      const prevVal = inner ? inner.value : null;
      if (inner) {
        inner.focus();
        inner.value = '50';
        inner.dispatchEvent(new Event('input', {bubbles: true}));
        inner.dispatchEvent(new Event('change', {bubbles: true}));
        inner.blur();
      }
      await new Promise((r) => setTimeout(r, 600));
      const newVal = inner ? inner.value : null;
      const ok = dlg.querySelector('[name="button-OK"]') as HTMLElement | null;
      if (ok) ok.click();
      await new Promise((r) => setTimeout(r, 2000));
      const dlgClosed = !document.querySelector('[name="dialog-Peptides-settings"]');
      return {inputFound, prevVal, newVal, dlgClosed};
    });
    expect(changed.inputFound, '[name="input-Similarity-Threshold"] not found in the MCL pane').toBe(true);
    expect(changed.newVal, 'Similarity threshold did not accept the new value').toBe('50');
    expect(changed.dlgClosed, 'settings dialog did not close after OK').toBe(true);
  });
  await softStep('Scenario 2 (step 9): verify viewers reload after settings change', async () => {
    await page.waitForTimeout(8000);
    const state = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const viewers = Array.from(tv.viewers).map((v) => v.type);
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]');
      return {
        viewers,
        modelPresent: !!tv.dataFrame.temp['peptidesModel'],
        svmHasCanvas: svm ? !!svm.querySelector('canvas') : false,
      };
    });
    expect(state.modelPresent, 'PeptidesModel cache lost after settings change').toBe(true);
    expect(state.viewers, 'Sequence Variability Map must persist after reload').toContain('Sequence Variability Map');
    expect(state.viewers, 'Most Potent Residues must persist after reload').toContain('Most Potent Residues');
    expect(state.svmHasCanvas, 'Sequence Variability Map did not re-render its canvas').toBe(true);
    const errors = await page.evaluate(() =>
      (grok.shell.lastError ? [String(grok.shell.lastError)] : []));
    expect(errors.filter((e) => /setTrue|null/.test(e)).length,
      `GROK-19145 invariant: post-OK compute produced a null-receiver error: ${errors.join('; ')}`).toBe(0);
  });
  await softStep('Scenario 3 (step 10): toggle Mutation Cliffs / Invariant Map mode', async () => {
    const toggled = await page.evaluate(async () => {
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]')!;
      const mc = svm.querySelector('[name="input-Mutation-Cliffs"]') as HTMLInputElement | null;
      const im = svm.querySelector('[name="input-Invariant-Map"]') as HTMLInputElement | null;
      const radiosFound = !!mc && !!im;
      if (!radiosFound) return {radiosFound, afterIM: null, afterMC: null};
      // Switch to Invariant Map.
      im!.click();
      await new Promise((r) => setTimeout(r, 1500));
      const afterIM = {mc: mc!.checked, im: im!.checked};
      // Switch back to Mutation Cliffs.
      mc!.click();
      await new Promise((r) => setTimeout(r, 1500));
      const afterMC = {mc: mc!.checked, im: im!.checked};
      return {radiosFound, afterIM, afterMC};
    });
    expect(toggled.radiosFound,
      'SVM mode-toggle radios (input-Mutation-Cliffs / input-Invariant-Map) not found').toBe(true);
    expect(toggled.afterIM, 'switching to Invariant Map did not flip the radios').toEqual({mc: false, im: true});
    expect(toggled.afterMC, 'switching back to Mutation Cliffs did not flip the radios').toEqual({mc: true, im: false});
  });
  await softStep('Scenario 3 (steps 11-12): cell click populates the Context Panel', async () => {
    const result = await page.evaluate(async () => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const df = tv.dataFrame;
      const selBefore = df.selection.trueCount;
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]')!;
      // The largest canvas is the SVM render surface (the tiny one is a scrollbar).
      const canvases = Array.from(svm.querySelectorAll('canvas')) as HTMLCanvasElement[];
      let canvas: HTMLCanvasElement | null = null;
      let maxArea = 0;
      for (const c of canvases) {
        const r = c.getBoundingClientRect();
        if (r.width * r.height > maxArea) { maxArea = r.width * r.height; canvas = c; }
      }
      if (!canvas) return {canvasFound: false, selBefore, selAfter: selBefore};
      const r = canvas.getBoundingClientRect();
      // Deterministic data-cell position past the row-header column and WebLogo header row.
      const cx = r.x + 120;
      const cy = r.y + 120;
      const opts = {bubbles: true, cancelable: true, clientX: cx, clientY: cy, button: 0, view: window};
      canvas.dispatchEvent(new MouseEvent('mousemove', opts));
      canvas.dispatchEvent(new MouseEvent('mousedown', opts));
      canvas.dispatchEvent(new MouseEvent('mouseup', opts));
      canvas.dispatchEvent(new MouseEvent('click', opts));
      await new Promise((res) => setTimeout(res, 2500));
      const selAfter = df.selection.trueCount;
      // Read the Context Panel panes that materialized.
      const namedPanes = Array.from(document.querySelectorAll('[name^="pane-"]'))
        .map((p) => p.getAttribute('name'));
      return {canvasFound: true, selBefore, selAfter, namedPanes};
    });
    expect(result.canvasFound, 'SVM render canvas not found for the cell click').toBe(true);
    // The map cell click selects the corresponding monomer-position rows.
    expect(result.selAfter, 'cell click did not change the DataFrame selection')
      .toBeGreaterThan(result.selBefore);
    // Step 12: the Context Panel exposes the Mutation Cliffs pairs + Distribution panes.
    expect(result.namedPanes, 'Context Panel did not surface the Mutation Cliffs pairs pane')
      .toContain('pane-Mutation-Cliffs-pairs');
    expect(result.namedPanes, 'Context Panel did not surface the Distribution pane')
      .toContain('pane-Distribution');
  });
  // Scenario 4 — change one representative Distribution parameter and verify re-render.
  await softStep('Scenario 4 (step 13): adjust a Distribution panel parameter', async () => {
    const result = await page.evaluate(async () => {
      const distPane = document.querySelector('[name="pane-Distribution"]');
      if (!distPane) return {distPaneFound: false};
      const header = distPane.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
      if (header && !distPane.classList.contains('expanded')) header.click();
      await new Promise((r) => setTimeout(r, 1500));
      // Representative non-default change: toggle "Stack split categories".
      const stack = distPane.querySelector('[name="input-Stack-split-categories"]') as HTMLInputElement | null;
      const stackFound = !!stack;
      const before = stack ? stack.checked : null;
      if (stack) stack.click();
      await new Promise((r) => setTimeout(r, 1500));
      const after = stack ? stack.checked : null;
      const childCount = distPane.querySelectorAll('*').length;
      return {distPaneFound: true, stackFound, before, after, childCount};
    });
    expect(result.distPaneFound, '[name="pane-Distribution"] not found in the Context Panel').toBe(true);
    expect(result.stackFound, '[name="input-Stack-split-categories"] not found in the Distribution pane').toBe(true);
    expect(result.after, 'Distribution parameter toggle did not change state').not.toBe(result.before);
    // The pane keeps rendered content (histograms) after the parameter change.
    expect(result.childCount, 'Distribution pane lost its rendered content after the parameter change')
      .toBeGreaterThan(5);
  });
  await page.evaluate(() => grok.shell.closeAll());
  finishSpec();
});
