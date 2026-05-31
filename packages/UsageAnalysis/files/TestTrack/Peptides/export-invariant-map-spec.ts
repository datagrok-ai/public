/* ---
sub_features_covered: [peptides.viewers.sar-base.export-invariant-map, peptides.viewers.sar-base, peptides.viewers.monomer-position, peptides.workflow.start-analysis]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (coverage_type: smoke; non-ui-smoke — JS API substitution permitted)
//   sub_features_covered: [peptides.viewers.sar-base.export-invariant-map,
//     peptides.viewers.sar-base, peptides.viewers.monomer-position,
//     peptides.workflow.start-analysis]
//   ui_coverage_responsibility: [] (delegated_to: null)
//   related_bugs: []
//
// Atlas provenance (derived_from):
//   feature-atlas/peptides.yaml#critical_paths[export-invariant-map-from-sar-viewer]
//     derived_from: public/packages/Peptides/src/viewers/sar-viewer.ts#L674
//   feature-atlas/peptides.yaml#sub_features[peptides.viewers.sar-base]
//     derived_from: public/packages/Peptides/src/viewers/sar-viewer.ts#L104
//
// Export Invariant Map from the SARViewer context menu (sister of export-mutation-cliffs.md).
// The Export Invariant Map command is registered on the SARViewer *base* class
// (sar-viewer.ts onTableAttached: exportGroup.item('Export Invariant Map', () =>
// this.exportInvariantMap())), so it is available from BOTH subclass entry points —
// MonomerPosition (Sequence Variability Map, Scenario 1) and MostPotentResidues
// (Most Potent Residues, Scenario 2) — with the same exported-table shape.
// Selectors verified against dev.datagrok.ai @datagrok/peptides v1.27.9 per
// .claude/skills/grok-browser/references/peptides.md.
//
// Empirical recon notes (live chrome-devtools MCP, dev.datagrok.ai, 2026-05-29 — drive
// deterministic assertions, not theory):
//   - SAR launch via the context-panel Launch SAR button is the deterministic path
//     (top-menu Bio | Analyze | SAR... falls into menu overflow at the recon viewport).
//     Default attach set on peptides.csv: Grid, Sequence Variability Map,
//     Most Potent Residues, MCL. SAR split AlignedSequence into 17 per-position columns (1..17).
//   - `Export Invariant Map` lives on the *viewer* context menu (a contextmenu dispatch on
//     the viewer body surfaces the SARViewer-contributed Export group merged with the
//     generic grid menu). The `Export Invariant Map` `.d4-menu-item-label` renders at the
//     menu top level (NO `...` suffix, unlike `Export Mutation Cliffs...`). Clicking its
//     `.d4-menu-item` calls exportInvariantMap() directly — there is NO column-picker dialog
//     (the simpler sibling of the mutation-cliffs export, which opens a dialog).
//   - exportInvariantMap() builds + addTableView's a new DataFrame named exactly
//     `Invariant Map` that becomes the active view. Verified shape:
//       leading column `AAR` (type string — this is C.COLUMNS_NAMES.MONOMER === 'AAR',
//       the monomer identifier; the scenario's "Monomer column" is conceptual, the
//       ACTUAL column name is `AAR`), then one INT column per sequence position
//       (named `1`..`17` — matching the SAR per-position columns on the original df).
//       647-row demo → 22 monomer rows; AAR values are single-letter amino acids
//       (A, C, D, E, F, G, H, I, K, ...) plus HELM-shape labels (COOH).
//       Each position column's cell values are integer sequence-counts (a position
//       column sums to the total sequence count 647; many cells are 0 where no sequence
//       carries that monomer at that position).
//   - Scenario 2 (Most Potent Residues entry point): the MPR viewer is docked zero-width
//     in the default SAR split layout, so a bounding-rect canvas click cannot reach it.
//     A `contextmenu` MouseEvent dispatched on an element INSIDE `mprViewer.root` DOES
//     surface the SARViewer-base Export group (the menu registration filters on
//     `this.root.contains(target)`, which the event path satisfies regardless of the
//     visual zero-width). Verified: the `Export Invariant Map` item is present on the MPR
//     context menu and clicking it exports the identical-shape `Invariant Map` table with
//     no thrown error — confirming the base-class Export contract holds from the MPR subclass.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/bio/peptides.csv';

// The exported invariant-map table's leading column is named `AAR`
// (C.COLUMNS_NAMES.MONOMER === 'AAR' in Peptides/src/utils/constants.ts).
const MONOMER_COL = 'AAR';

test('Peptides — Export Invariant Map to a new TableView (SARViewer-base contract)', async ({page}) => {
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
    // Give MCL/sequence-space compute + viewer attach extra settle time so the SAR viewers
    // and their monomer-position statistics are ready before the export.
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

  // The number of per-position columns the SAR analysis produced on the ORIGINAL TableView.
  // The exported invariant-map table must carry exactly this many position columns.
  let originalPositionColCount = 0;
  await softStep('Setup (step 3): record the SAR per-position column count on the original table', async () => {
    originalPositionColCount = await page.evaluate(() => {
      const sarTv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel'])!;
      // Position columns are the integer-named split columns (1, 2, 3, ...).
      return sarTv.dataFrame.columns.names().filter((n) => /^\d+$/.test(n)).length;
    });
    expect(originalPositionColCount,
      'SAR analysis produced no per-position columns on the original table').toBeGreaterThan(0);
  });

  // ---- Scenario 1 — Export Invariant Map from the Sequence Variability Map (MonomerPosition) ----

  // Steps 1-3: right-click the SVM viewer body, hover the Export submenu, click
  // Export Invariant Map; a new `Invariant Map` TableView opens and becomes active.
  await softStep('Scenario 1 (steps 1-3): Export Invariant Map opens a new active TableView', async () => {
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
      if (!canvas) return {canvasFound: false, menuItemFound: false};
      const r = canvas.getBoundingClientRect();
      const cx = r.x + 120, cy = r.y + 120;
      // Step 1: open the viewer context menu on the SVM body.
      canvas.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: cx, clientY: cy, view: window}));
      canvas.dispatchEvent(new MouseEvent('contextmenu',
        {bubbles: true, cancelable: true, button: 2, clientX: cx, clientY: cy, view: window}));
      await new Promise((rs) => setTimeout(rs, 1200));

      // Step 2: hover the Export group so its items render, then click Export Invariant Map.
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
      // `Export Invariant Map` (no `...` suffix; calls exportInvariantMap() directly, no dialog).
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

  // Steps 4-7: the new `Invariant Map` TableView carries the documented shape.
  await softStep('Scenario 1 (steps 4-7): exported TableView has the monomer-by-position count shape', async () => {
    // The export builds a DataFrame + addTableView; allow a brief settle.
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
      // First position column's cell-value sanity (integers, including zeros).
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
      // Monomer-identifier sanity: first value should look like a monomer label
      // (single-letter amino acid or a short HELM-shape label).
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

    // Step 3 (tail) + Step: a new TableView opened and became the active view.
    expect(result.found, 'new `Invariant Map` TableView was not opened by the export').toBe(true);
    expect(result.activeViewName, 'the exported `Invariant Map` view did not become active')
      .toBe('Invariant Map');
    // Step 4: a leading monomer-identifier column (named `AAR` in the platform).
    expect(result.leadingCol, 'exported grid is missing its leading monomer column').toBe(MONOMER_COL);
    expect(result.leadingColType, 'the leading monomer column must be a string column').toBe('string');
    // Step 5: one column per sequence position, matching the original SAR per-position count.
    expect(result.positionColCount,
      'exported position-column count does not match the original SAR per-position columns')
      .toBe(originalPositionColCount);
    expect(result.allPositionColsInt,
      'every per-position column must be an integer (sequence-count) column').toBe(true);
    // Step 6: at least one monomer row; position-column cell values are integers incl. zero.
    expect(result.rows, 'exported invariant-map grid has no monomer rows').toBeGreaterThan(0);
    expect(result.firstPosAllIntegers,
      'position-column cell values must be integers (sequence counts)').toBe(true);
    expect(result.firstPosAnyZero,
      'a per-position column should carry zero counts where a monomer is absent at that position')
      .toBe(true);
    // Step 7: the leading column carries a recognizable monomer identifier.
    expect(result.monomerLooksLikeId,
      `leading column value '${result.firstMonomer}' is not a recognizable monomer identifier`)
      .toBe(true);
  });

  // ---- Scenario 2 — Export Invariant Map from Most Potent Residues (SARViewer-base contract) ----

  // Step 1: return to the original SAR TableView (close the exported view).
  // Steps 2-4: right-click the Most Potent Residues viewer body, click Export Invariant Map;
  //   a new `Invariant Map` TableView opens with the SAME shape as Scenario 1.
  // Step 5-6: same shape + no error — the SARViewer-base Export contract holds from the
  //   MostPotentResidues subclass entry point, not just MonomerPosition.
  await softStep('Scenario 2 (steps 1-6): Export Invariant Map from Most Potent Residues — same shape, no error', async () => {
    // Step 1: close the prior exported view and switch back to the SAR TableView.
    await page.evaluate(async () => {
      Array.from(grok.shell.tableViews)
        .filter((v) => v.dataFrame && v.dataFrame.name === 'Invariant Map')
        .forEach((v) => v.close());
      await new Promise((r) => setTimeout(r, 600));
      const sarTv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']);
      if (sarTv) { grok.shell.v = sarTv; await new Promise((r) => setTimeout(r, 1000)); }
    });

    // Capture console errors during the MPR export to assert the SARViewer-base contract
    // holds without a `Cannot read ... on null`-class failure.
    const consoleErrors: string[] = [];
    const onConsole = (msg: import('@playwright/test').ConsoleMessage) => {
      if (msg.type() === 'error') consoleErrors.push(msg.text());
    };
    page.on('console', onConsole);

    // Steps 2-3: open the MPR context menu and click Export Invariant Map.
    // The MPR viewer is docked zero-width in the default SAR split layout, so a
    // bounding-rect canvas click cannot reach it; a `contextmenu` MouseEvent dispatched on
    // an element inside `mprViewer.root` surfaces the SARViewer-base Export group (the menu
    // registration filters on `this.root.contains(target)`, satisfied by the event path).
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
      // Hover the Export group, then click Export Invariant Map.
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

  // Cleanup.
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
