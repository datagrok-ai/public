import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
const datasetPath = 'System:DemoFiles/bio/peptides.csv';
const BASE_EXPORT_COLS = ['Seq 1', 'Seq 2', 'Mutation', 'Seq 1 IC50', 'Seq 2 IC50', 'Delta'];
test('Peptides — Export Mutation Cliffs to a new TableView', async ({page}) => {
  test.setTimeout(300_000);
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
  await softStep('Setup (step 2): launch SAR from the Peptides context panel', async () => {
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
      const result = {paneFound: !!pane, launchFound: !!launchBtn};
      if (launchBtn) launchBtn.click();
      return result;
    });
    expect(setup.paneFound, '[name="pane-Peptides"] context pane not found (waited 30s)').toBe(true);
    expect(setup.launchFound, '[name="button-Launch-SAR"] not found in Peptides pane (waited 30s)').toBe(true);
    await page.waitForFunction(() => {
      return Array.from(grok.shell.tableViews).some((v) => v.dataFrame.temp['peptidesModel']);
    }, {timeout: 60000});
    await page.waitForTimeout(9000);
    const svmReady = await page.evaluate(() => {
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]');
      return {svmFound: !!svm, svmHasCanvas: svm ? !!svm.querySelector('canvas') : false};
    });
    expect(svmReady.svmFound, 'Sequence Variability Map viewer did not attach after SAR launch').toBe(true);
    expect(svmReady.svmHasCanvas, 'Sequence Variability Map did not render its canvas').toBe(true);
  });
  await softStep('Scenario 1 (steps 1-3): open Export Mutation Cliffs dialog and accept defaults', async () => {
    const opened = await page.evaluate(async () => {
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]')!;
      const canvases = Array.from(svm.querySelectorAll('canvas')) as HTMLCanvasElement[];
      let canvas: HTMLCanvasElement | null = null;
      let maxArea = 0;
      for (const c of canvases) {
        const r = c.getBoundingClientRect();
        if (r.width * r.height > maxArea) { maxArea = r.width * r.height; canvas = c; }
      }
      if (!canvas) return {canvasFound: false};
      const r = canvas.getBoundingClientRect();
      const cx = r.x + 120, cy = r.y + 120;
      canvas.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: cx, clientY: cy, view: window}));
      canvas.dispatchEvent(new MouseEvent('contextmenu',
        {bubbles: true, cancelable: true, button: 2, clientX: cx, clientY: cy, view: window}));
      await new Promise((rs) => setTimeout(rs, 1200));
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
  await softStep('Scenario 1 (steps 4-6): exported TableView has the documented column shape', async () => {
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
        hasAllBase: baseCols.every((c: string) => df.columns.names().includes(c)),
      };
    }, BASE_EXPORT_COLS);
    expect(result.found, 'new `Mutation Cliffs` TableView was not opened by the export').toBe(true);
    expect(result.activeViewName, 'the exported `Mutation Cliffs` view did not become active')
      .toBe('Mutation Cliffs');
    expect(result.cols, 'exported grid is missing one or more documented columns').toEqual(BASE_EXPORT_COLS);
    expect(result.hasAllBase, 'exported grid does not carry the full documented base column set').toBe(true);
    expect(result.mutationSemType, 'Mutation column must have semType MacromoleculeDifference')
      .toBe('MacromoleculeDifference');
    expect(result.mutationSample, 'Mutation value must be formatted `seq1#seq2`')
      .toContain('#');
    expect(result.rows, 'exported mutation-cliffs grid is empty').toBeGreaterThan(0);
  });
  await softStep('Scenario 2 (steps 1-4): extra column selection emits paired columns', async () => {
    await page.evaluate(async () => {
      const sarTv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']);
      if (sarTv) { grok.shell.v = sarTv; await new Promise((r) => setTimeout(r, 1000)); }
    });
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
      svm._doExportMutationCliffs([idCol]);
      return {extraName, existingCount: existingMcIds.size};
    });
    expect(exported.extraName, 'demo dataset is missing the `ID` extra column').toBe('ID');
    await page.waitForTimeout(2500);
    const result = await page.evaluate((baseCols) => {
      const mcViews = Array.from(grok.shell.tableViews)
        .filter((v) => v.dataFrame && v.dataFrame.name.startsWith('Mutation Cliffs'));
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
  await page.evaluate(() => grok.shell.closeAll());
  finishSpec();
});
