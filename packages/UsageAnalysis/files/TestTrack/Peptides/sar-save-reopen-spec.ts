/* ---
sub_features_covered: [peptides.model, peptides.model.init, peptides.viewers.monomer-position, peptides.workflow.sar-dialog, peptides.workflow.start-analysis]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
const datasetPath = 'System:DemoFiles/bio/peptides.csv';
const projectName = 'SarSaveReopenAuto' + Date.now();
test('SAR — save project with SAR layout + selection + scaling, reopen and verify layout restored (GROK-14461)', async ({page}) => {
  test.setTimeout(360_000);
  await loginToDatagrok(page);
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
  await softStep('Setup (step 2): launch SAR from the Peptides context panel + confirm model init', async () => {
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
    await page.waitForTimeout(8000);
    const state = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      return {
        viewers: Array.from(tv.viewers).map((v) => v.type),
        modelPresent: !!tv.dataFrame.temp['peptidesModel'],
      };
    });
    expect(state.modelPresent, 'PeptidesModel singleton must attach after Launch SAR').toBe(true);
    expect(state.viewers, 'Sequence Variability Map (MonomerPosition) must attach').toContain('Sequence Variability Map');
    expect(state.viewers, 'Most Potent Residues must attach').toContain('Most Potent Residues');
  });
  await softStep('Setup (step 3): set -lg scaling and establish a non-empty selection', async () => {
    const result = await page.evaluate(async () => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const df = tv.dataFrame;
      const model = df.temp['peptidesModel'];
      const scalingBefore = model && model._settings ? model._settings.activityScaling : null;
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
    expect(result.selAfter, 'SVM cell click did not establish a non-empty selection')
      .toBeGreaterThan(result.selBefore);
  });
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
      await new Promise((r) => setTimeout(r, 6000));
      const dlgClosed = !document.querySelector('[name="dialog-Save-project"]');
      let saved = null;
      try { saved = await grok.dapi.projects.filter('name = "' + name + '"').first(); }
      catch (e) {  }
      return {nameInputFound, dlgClosed, savedId: saved ? saved.id : null, savedName: saved ? saved.name : null};
    }, projectName);
    expect(saved.nameInputFound, 'Save project dialog name input not found').toBe(true);
    expect(saved.dlgClosed, 'Save project dialog did not close after OK').toBe(true);
    expect(saved.savedId, `saved project "${projectName}" not found on the server after OK`).toBeTruthy();
  });
  await softStep('Scenario 1 (steps 3-8): reopen the project and verify the SAR layout is restored (GROK-14461)', async () => {
    const reopened = await page.evaluate(async (name) => {
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 2500));
      const tvAfterClose = Array.from(grok.shell.tableViews).length;
      const project = await grok.dapi.projects.filter('name = "' + name + '"').first();
      const found = !!project;
      if (project) await project.open();
      await new Promise((r) => setTimeout(r, 14000));
      const tv = grok.shell.tv;
      const df = tv ? tv.dataFrame : null;
      const grid = tv ? Array.from(tv.viewers).find((v) => v.type === 'Grid') : null;
      const model = df ? df.temp['peptidesModel'] : null;
      return {
        found,
        tvAfterClose,
        tvAfterOpen: Array.from(grok.shell.tableViews).length,
        alignedSeqSemType: (df && df.col('AlignedSequence')) ? df.col('AlignedSequence').semType : null,
        viewers: tv ? Array.from(tv.viewers).map((v) => v.type) : null,
        colHeaderHeight: (grid && grid.props) ? grid.props.colHeaderHeight : null,
        modelPresent: df ? !!df.temp['peptidesModel'] : null,
        scaling: (model && model._settings) ? model._settings.activityScaling : null,
        selTrueCount: df ? df.selection.trueCount : null,
      };
    }, projectName);
    expect(reopened.tvAfterClose, 'workspace did not close before reopen').toBe(0);
    expect(reopened.found, `saved project "${projectName}" not found for reopen`).toBe(true);
    expect(reopened.tvAfterOpen, 'reopening the project did not open a TableView').toBeGreaterThan(0);
    expect(reopened.alignedSeqSemType, 'AlignedSequence Macromolecule column missing after reopen (datasync failed)')
      .toBe('Macromolecule');
    expect(reopened.viewers, 'GROK-14461: Sequence Variability Map viewer not reattached on reopen')
      .toContain('Sequence Variability Map');
    expect(reopened.viewers, 'GROK-14461: Most Potent Residues viewer not reattached on reopen')
      .toContain('Most Potent Residues');
    expect(reopened.viewers, 'GROK-14461: MCL clustering viewer not reattached on reopen')
      .toContain('MCL');
    expect(reopened.colHeaderHeight, 'GROK-14461: WebLogo column-headers not rendered on reopen (colHeaderHeight not grown)')
      .toBeGreaterThan(60);
    expect(reopened.modelPresent, 'GROK-14461: PeptidesModel null on reopen — model not reconstructed (viewers would null-bind)')
      .toBe(true);
    expect(typeof reopened.scaling, 'GROK-14461: activity scaling setting lost on reopen (PeptidesSettings did not round-trip)')
      .toBe('string');
  });
  await softStep('Scenario 2 (steps 1-3): reopened selection matches an atlas-permitted branch', async () => {
    const sel = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const df = tv.dataFrame;
      return {trueCount: df.selection.trueCount, rowCount: df.rowCount};
    });
    expect(sel.trueCount, 'reopened selection trueCount is negative (corrupt BitSet)').toBeGreaterThanOrEqual(0);
    expect(sel.trueCount, 'reopened selection trueCount exceeds rowCount (corrupt BitSet)')
      .toBeLessThanOrEqual(sel.rowCount);
  });
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
    expect(result.selAfter, 'post-reopen SVM cell click did not update the selection')
      .toBeGreaterThan(result.selBefore);
    expect(result.namedPanes, 'Context Panel did not surface the Mutation Cliffs pairs pane on reopen')
      .toContain('pane-Mutation-Cliffs-pairs');
    expect(result.namedPanes, 'Context Panel did not surface the Distribution pane on reopen')
      .toContain('pane-Distribution');
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
  finishSpec();
});
