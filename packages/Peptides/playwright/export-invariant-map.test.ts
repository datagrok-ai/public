/* ---
sub_features_covered: [peptides.viewers.monomer-position, peptides.viewers.sar-base, peptides.viewers.sar-base.export-invariant-map, peptides.workflow.start-analysis]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';
test.use(specTestOptions);
const datasetPath = 'System:DemoFiles/bio/peptides.csv';
const MONOMER_COL = 'AAR';
test('Peptides — Export Invariant Map to a new TableView (SARViewer-base contract)', async ({page}) => {
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
    const ready = await page.evaluate(() => {
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]');
      const mpr = document.querySelector('[name="viewer-Most-Potent-Residues"]');
      return {
        svmFound: !!svm,
        svmHasCanvas: svm ? !!svm.querySelector('canvas') : false,
        mprFound: !!mpr,
      };
    });
    expect(ready.svmFound, 'Sequence Variability Map viewer did not attach after SAR launch').toBe(true);
    expect(ready.svmHasCanvas, 'Sequence Variability Map did not render its canvas').toBe(true);
    expect(ready.mprFound, 'Most Potent Residues viewer did not attach after SAR launch').toBe(true);
  });
  let originalPositionColCount = 0;
  await softStep('Setup (step 3): record the SAR per-position column count on the original table', async () => {
    originalPositionColCount = await page.evaluate(() => {
      const sarTv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel'])!;
      return sarTv.dataFrame.columns.names().filter((n) => /^\d+$/.test(n)).length;
    });
    expect(originalPositionColCount,
      'SAR analysis produced no per-position columns on the original table').toBeGreaterThan(0);
  });
  await softStep('Scenario 1 (steps 1-3): Export Invariant Map opens a new active TableView', async () => {
    const opened = await page.evaluate(async () => {
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]')!;
      const canvases = Array.from(svm.querySelectorAll('canvas')) as HTMLCanvasElement[];
      let canvas: HTMLCanvasElement | null = null;
      let maxArea = 0;
      for (const c of canvases) {
        const r = c.getBoundingClientRect();
        if (r.width * r.height > maxArea) { maxArea = r.width * r.height; canvas = c; }
      }
      if (!canvas) return {canvasFound: false, menuItemFound: false};
      const r = canvas.getBoundingClientRect();
      const cx = r.x + 120, cy = r.y + 120;
      canvas.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: cx, clientY: cy, view: window}));
      canvas.dispatchEvent(new MouseEvent('contextmenu',
        {bubbles: true, cancelable: true, button: 2, clientX: cx, clientY: cy, view: window}));
      await new Promise((rs) => setTimeout(rs, 1200));
      const exportGroup = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find((e) => e.textContent?.trim() === 'Export');
      if (exportGroup) {
        const g = exportGroup.closest('.d4-menu-item') as HTMLElement | null;
        if (g) {
          g.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, view: window}));
          g.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true, view: window}));
        }
      }
      await new Promise((rs) => setTimeout(rs, 600));
      const exportInv = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find((e) => e.textContent?.trim() === 'Export Invariant Map');
      const menuItemFound = !!exportInv;
      if (exportInv) {
        const item = exportInv.closest('.d4-menu-item') as HTMLElement | null;
        if (item) item.click();
      }
      return {canvasFound: true, menuItemFound};
    });
    expect(opened.canvasFound, 'SVM render canvas not found for the context-menu trigger').toBe(true);
    expect(opened.menuItemFound,
      '`Export Invariant Map` not found on the Sequence Variability Map context menu').toBe(true);
  });
  await softStep('Scenario 1 (steps 4-7): exported TableView has the monomer-by-position count shape', async () => {
    await page.waitForTimeout(3000);
    const result = await page.evaluate((args) => {
      const {monomerCol, expectedPosCols} = args;
      const tvs = Array.from(grok.shell.tableViews);
      const invView = tvs.find((v) => v.dataFrame && v.dataFrame.name === 'Invariant Map');
      if (!invView) return {found: false};
      const df = invView.dataFrame;
      const cols = df.columns.names();
      const monomer = df.col(monomerCol);
      const positionCols = cols.slice(1);
      const firstPos = positionCols.length > 0 ? df.col(positionCols[0]) : null;
      const firstPosSamples: number[] = [];
      let allIntegers = true;
      let anyZero = false;
      if (firstPos) {
        for (let i = 0; i < df.rowCount; i++) {
          const v = firstPos.get(i);
          if (i < 8) firstPosSamples.push(v);
          if (!Number.isInteger(v)) allIntegers = false;
          if (v === 0) anyZero = true;
        }
      }
      const firstMonomer = monomer && df.rowCount > 0 ? String(monomer.get(0)) : null;
      const monomerLooksLikeId = firstMonomer != null && /^[A-Za-z0-9()\[\]_-]{1,8}$/.test(firstMonomer);
      return {
        found: true,
        activeViewName: grok.shell.v ? grok.shell.v.name : null,
        cols,
        rows: df.rowCount,
        leadingCol: cols[0],
        leadingColType: monomer ? monomer.type : null,
        positionColCount: positionCols.length,
        allPositionColsInt: positionCols.every((c: string) => df.col(c).type === 'int'),
        firstPosSamples,
        firstPosAllIntegers: allIntegers,
        firstPosAnyZero: anyZero,
        firstMonomer,
        monomerLooksLikeId,
        expectedPosCols,
      };
    }, {monomerCol: MONOMER_COL, expectedPosCols: originalPositionColCount});
    expect(result.found, 'new `Invariant Map` TableView was not opened by the export').toBe(true);
    expect(result.activeViewName, 'the exported `Invariant Map` view did not become active')
      .toBe('Invariant Map');
    expect(result.leadingCol, 'exported grid is missing its leading monomer column').toBe(MONOMER_COL);
    expect(result.leadingColType, 'the leading monomer column must be a string column').toBe('string');
    expect(result.positionColCount,
      'exported position-column count does not match the original SAR per-position columns')
      .toBe(originalPositionColCount);
    expect(result.allPositionColsInt,
      'every per-position column must be an integer (sequence-count) column').toBe(true);
    expect(result.rows, 'exported invariant-map grid has no monomer rows').toBeGreaterThan(0);
    expect(result.firstPosAllIntegers,
      'position-column cell values must be integers (sequence counts)').toBe(true);
    expect(result.firstPosAnyZero,
      'a per-position column should carry zero counts where a monomer is absent at that position')
      .toBe(true);
    expect(result.monomerLooksLikeId,
      `leading column value '${result.firstMonomer}' is not a recognizable monomer identifier`)
      .toBe(true);
  });
  await softStep('Scenario 2 (steps 1-6): Export Invariant Map from Most Potent Residues — same shape, no error', async () => {
    await page.evaluate(async () => {
      Array.from(grok.shell.tableViews)
        .filter((v) => v.dataFrame && v.dataFrame.name === 'Invariant Map')
        .forEach((v) => v.close());
      await new Promise((r) => setTimeout(r, 600));
      const sarTv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']);
      if (sarTv) { grok.shell.v = sarTv; await new Promise((r) => setTimeout(r, 1000)); }
    });
    const consoleErrors: string[] = [];
    const onConsole = (msg: import('@playwright/test').ConsoleMessage) => {
      if (msg.type() === 'error') consoleErrors.push(msg.text());
    };
    page.on('console', onConsole);
    const triggered = await page.evaluate(async () => {
      const sarTv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel'])!;
      const mprViewer = Array.from(sarTv.viewers).find((v) => v.type === 'Most Potent Residues') as any;
      if (!mprViewer || !mprViewer.root) return {viewerFound: false, menuItemFound: false, threw: 'no MPR viewer'};
      const root = mprViewer.root as HTMLElement;
      const inner = (root.querySelector('canvas') || root.firstElementChild || root) as HTMLElement;
      inner.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, view: window, clientX: 320, clientY: 150}));
      inner.dispatchEvent(new MouseEvent('contextmenu',
        {bubbles: true, cancelable: true, button: 2, view: window, clientX: 320, clientY: 150}));
      await new Promise((rs) => setTimeout(rs, 1200));
      const exportGroup = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find((e) => e.textContent?.trim() === 'Export');
      if (exportGroup) {
        const g = exportGroup.closest('.d4-menu-item') as HTMLElement | null;
        if (g) {
          g.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, view: window}));
          g.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true, view: window}));
        }
      }
      await new Promise((rs) => setTimeout(rs, 600));
      const exportInv = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find((e) => e.textContent?.trim() === 'Export Invariant Map');
      const menuItemFound = !!exportInv;
      let threw: string | null = null;
      if (exportInv) {
        const item = exportInv.closest('.d4-menu-item') as HTMLElement | null;
        try { if (item) item.click(); }
        catch (e) { threw = String(e); }
      }
      return {viewerFound: true, menuItemFound, threw};
    });
    expect(triggered.viewerFound, 'Most Potent Residues viewer instance not found on the SAR table view').toBe(true);
    expect(triggered.menuItemFound,
      '`Export Invariant Map` not found on the Most Potent Residues context menu (SARViewer-base contract)')
      .toBe(true);
    expect(triggered.threw, `Export Invariant Map click threw from the MPR entry point: ${triggered.threw}`)
      .toBeNull();
    // Step 4-5: the exported `Invariant Map` TableView has the same documented shape.
    await page.waitForTimeout(3000);
    const result = await page.evaluate((args) => {
      const {monomerCol, expectedPosCols} = args;
      const invView = Array.from(grok.shell.tableViews)
        .find((v) => v.dataFrame && v.dataFrame.name === 'Invariant Map');
      if (!invView) return {found: false};
      const df = invView.dataFrame;
      const cols = df.columns.names();
      const positionCols = cols.slice(1);
      return {
        found: true,
        activeViewName: grok.shell.v ? grok.shell.v.name : null,
        leadingCol: cols[0],
        leadingColType: df.col(cols[0]) ? df.col(cols[0])!.type : null,
        rows: df.rowCount,
        positionColCount: positionCols.length,
        allPositionColsInt: positionCols.every((c: string) => df.col(c).type === 'int'),
        expectedPosCols,
      };
    }, {monomerCol: MONOMER_COL, expectedPosCols: originalPositionColCount});
    page.off('console', onConsole);
    // Step 4: a new TableView opened (and became active).
    expect(result.found, 'Most Potent Residues export did not open a new `Invariant Map` TableView').toBe(true);
    expect(result.activeViewName, 'the MPR-exported `Invariant Map` view did not become active')
      .toBe('Invariant Map');
    // Step 5: same shape as Scenario 1 — leading `AAR` column + one INT column per position.
    expect(result.leadingCol, 'MPR-exported grid is missing its leading monomer column').toBe(MONOMER_COL);
    expect(result.leadingColType, 'the MPR-exported leading monomer column must be a string column')
      .toBe('string');
    expect(result.positionColCount,
      'MPR-exported position-column count does not match the original SAR per-position columns')
      .toBe(originalPositionColCount);
    expect(result.allPositionColsInt,
      'every MPR-exported per-position column must be an integer (sequence-count) column').toBe(true);
    expect(result.rows, 'MPR-exported invariant-map grid has no monomer rows').toBeGreaterThan(0);
    // Step 6: no `Cannot read ... on null`-class console error during the MPR Export action.
    const nullRefErrors = consoleErrors.filter((e) => /Cannot read .* (of|on) null/i.test(e));
    expect(nullRefErrors,
      `MPR Export raised a null-reference console error (SARViewer-base contract violated): ${JSON.stringify(nullRefErrors)}`)
      .toEqual([]);
  });
  await page.evaluate(() => grok.shell.closeAll());
  finishSpec();
});
