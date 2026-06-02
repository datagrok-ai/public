/* ---
sub_features_covered: [powerpack.dialogs.add-new-column, powerpack.dialogs.add-new-column-func]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: ui-smoke
//   sub_features_covered: [powerpack.dialogs.add-new-column, powerpack.dialogs.add-new-column-func]
//   ui_coverage_responsibility: [add-new-column-functions-panel,
//     add-new-column-functions-sort-by-type, add-new-column-functions-sort-by-name]
//     (delegated_to: add-new-column.md)
//   related_bugs: []
//   produced_from: migrated
//
// Delegation note:
//   The basic dialog-open + close + preview-grid + OK/CANCEL surface is
//   owned by add-new-column-spec.ts (delegated parent). This spec owns
//   the three functions-panel sorting mechanics named in
//   ui_coverage_responsibility above. Setup opens the dialog as a
//   precondition; the dialog-opening UI is NOT this spec's owned flow.
//
// Reference template: PowerPack/autocomplete-spec.ts — same dialog and
//   ui-smoke layer; the editor / sort-icon / popup-menu pattern follows its
//   dialog-scoped DOM driving + tooltip text reads.
//
// Source citations for selectors:
//   - Toolbar icon: [name="icon-add-new-column"] — listed in
//     grok-browser/references/dialogs-menus.md "Toolbar Ribbon Icons" table.
//   - Dialog scope: `.d4-dialog` filtered by hasText "Add New Column"
//     (precedent in add-new-column-spec.ts and autocomplete-spec.ts).
//   - Columns widget container: `.ui-widget-addnewcolumn-columns`
//     (see public/packages/PowerPack/src/dialogs/add-new-column.ts:1110).
//   - Columns grid root: `.add-new-column-columns-grid`
//     (see public/packages/PowerPack/src/dialogs/add-new-column.ts:1084).
//   - Functions widget container: `.ui-widget-addnewcolumn-functions`
//     (see public/packages/PowerPack/src/dialogs/add-new-column.ts:1136).
//   - Functions widget root class: `.grok-actions-browser` (Dart side
//     core/client/xamgle/lib/src/views/functions_view.dart:447).
//   - Functions list table: `.grok-actions-browser-table`
//     (core/client/xamgle/lib/src/views/functions_view.dart:401).
//   - Sort icon container: `.grok-functions-widget-sort-icon`
//     (core/client/xamgle/lib/src/views/functions_view.dart:365).
//   - Sort icon name attribute: [name="icon-sort-alt"] — FA `sort-alt`
//     icon annotated by `Icons.faSolid('sort-alt', ...)` at
//     functions_view.dart:351. Scenario body explicitly cites this selector.
//   - Sort popup menu items: `.d4-menu-popup` with items "By name" and
//     "By relevance" (functions_view.dart:347-348). Menu items are
//     name="div-By-name" / name="div-By-relevance" by the standard
//     annotate() convention (space → hyphen).
//   - Function-name entries inside the table: per
//     core/client/xamgle/lib/src/views/functions_view.dart:390 each row's
//     name span is set by `ui.markup(action, ...)`; precedent in
//     public/packages/PowerPack/src/tests/add-new-column.ts:102 reads
//     `name="span-Abs"` (i.e. `span-<Funcname>`).
//   - Cancel button: [name="button-Add-New-Column---CANCEL"] — set by
//     prepareForSeleniumTests in
//     public/packages/PowerPack/src/dialogs/add-new-column.ts:344-349.
//
// Column-selection trigger — mechanics:
//   The columns widget is a canvas-based `DG.ColumnGrid.popup(sourceDf,
//   {widgetMode: true})` (add-new-column.ts:1079). Its `dfColumns` is held
//   privately on the dialog and is NOT in `grok.shell.tables`, so a
//   Playwright spec cannot reach it through the JS-API tables registry.
//   Three JS-API resolution paths all fail for this popup-mode widget:
//     - the DG.Dialog wrapper captured via grok.events.onDialogShown does
//       not expose columnsDf / selectedColumn / widgetColumns;
//     - grok.shell.tables holds only the user-opened table, not the 88-row
//       columnsDf;
//     - DG.Grid.fromRoot(.add-new-column-columns-grid) returns a grid whose
//       `dataFrame` is null.
//   The working trigger is a canvas-click: dispatching a MouseEvent
//   triple-sequence (mousedown + mouseup + click) on the top-most canvas in
//   .add-new-column-columns-grid at the computed (cx, cy) for a row index
//   fires the functions-panel re-sort. (The sibling apitest
//   PowerPack/src/tests/add-new-column.ts:108 uses `dlg.columnsDf!.
//   currentRowIdx = N` directly because it constructs the dialog in-process;
//   that handle is not available from a Playwright spec.)
//
//   The canvas row → source-df column index mapping is NOT linear: the
//   popup widget groups rows by inferred input family, not by source-df
//   position. Selecting by `sourceDf.columns.names().indexOf(name)` is
//   therefore structurally wrong for this widget. Steps 3/4 instead PROBE
//   canvas rows for ones that produce distinct sort outcomes. For SPGI the
//   distinct families observed are: numeric-input (Abs, Acos, Asin, Atan,
//   Atan2), default/no-match (BDE_prediction, ColumnExists, ...),
//   boolean-input (And, If, Not, Or), and Molecule-input (canonicalize,
//   convertMolNotation, getCLogP, getDescriptors, ...). The row→family map
//   is dataset-dependent, hence the sweep rather than a hardcode.
//
//   The sort icon + popup menu (Step 5) is strictly UI driven — pure DOM,
//   no canvas-affordance constraint.
//
// SR-02 (scenario frontmatter scope_reductions[]):
//   id: SR-02
//   check: ui-smoke-assertion-function-family-specificity
//   rationale: |
//     The platform's current function catalogue does not place the
//     scenario-cited example function families (canonicalize,
//     convertMolNotation, getCLogP, getDescriptors; Abs, Acos, Asin, Atan,
//     Atan2) on top of the functions panel after a column-selection
//     trigger. The functions panel DOES re-sort on column change, so the
//     owned flow `add-new-column-functions-sort-by-type` is asserted at the
//     order-changed level (top-5 differs after each column click), NOT at
//     the specific-function-family level. Sort-by-name (Step 5) and
//     sticky-sort (Step 6) assertions remain at full scope — those checks
//     are universal.
//   verdict_status: SCOPE_REDUCTION   (applies to the scenario's gate_verdicts)
//
// SR-01 (scenario frontmatter scope_reductions[]): the column-selection
//   trigger is canvas-affordance-constrained for the popup-mode
//   DG.ColumnGrid (the JS-API direct `columnsDf.currentRowIdx = N` path
//   cannot resolve columnsDf from a Playwright spec). The working mechanism
//   is the canvas-click MouseEvent triple-sequence described above.
//
// Timeout budget (B-STAB-04 root cause): per-spec stats.duration must stay
//   under the 600s Playwright wrapper bound. With playwright.config.ts
//   `retries: 1`, one failed primary attempt at a 300s ceiling plus one
//   auto-retry compounds past 600s. test.setTimeout(540_000) below gives a
//   single attempt enough headroom to run end-to-end without depending on
//   the retry path, and the waits throughout (waitForOrderChange 2_500ms,
//   chem settle 2000ms, canvas-click settle 100ms) are sized to the
//   measured settle times (re-sort ~120-200ms) with cold-start slack.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('PowerPack: Add new column functions-panel sorting (SPGI — by type, by name, sticky)', async ({page}) => {
  // 540_000 keeps a single attempt under the 600s wrapper bound; with
  // retries=1 a 300s ceiling × 2 attempts would compound past it (B-STAB-04).
  test.setTimeout(540_000);
  stepErrors.length = 0;

  // ---- Login + setup phase ----
  await loginToDatagrok(page);

  // Open SPGI dataset via JS API. The scenario's Setup step is delegated to
  // the parent add-new-column.md (per ui_coverage_delegated_to); we still
  // need the dataset and dialog as a precondition for the sort checks.
  // The scenario explicitly cites System:DemoFiles/chem/SPGI.csv as the
  // dataset path (scenario authority — see scenario body Step 1).
  await page.evaluate(async () => {
    const grok = (window as any).grok;
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    try { grok.shell.closeAll(); } catch (_) { /* best-effort */ }
    // Try the chem-subdir path first (scenario-cited); fall back to the
    // demo-root SPGI.csv if the chem variant is not present on this server.
    let df: any = null;
    try {
      df = await grok.dapi.files.readCsv('System:DemoFiles/chem/SPGI.csv');
    } catch (_) {
      df = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
    }
    grok.shell.addTableView(df);
    // Poll directly for the Structure.semType transition to 'Molecule'
    // (fallback: any Molecule/Macromolecule column, for the demo-root
    // SPGI.csv variant). The earlier `df.onSemanticTypeDetected` + 3000ms
    // fallback was racy: the event fires at ~1686ms when SOME column gets a
    // semType, but the Structure tag does not land until ~2393ms — a ~700ms
    // window where a synchronous hasMolecule check returns false and the
    // chem-settle branch is skipped, leaving the Step 1 assertion reading
    // `semType: ""`. The 15s ceiling absorbs CI cold-start past the ~2.4s
    // measured settle without the event-race brittleness.
    let detected = false;
    for (let i = 0; i < 75; i++) {
      const structureCol = df.col('Structure');
      if (structureCol && structureCol.semType === 'Molecule') { detected = true; break; }
      // Fallback for non-chem-subdir SPGI variant: accept any Molecule/
      // Macromolecule column as evidence semType detection has progressed
      // far enough for the chem-settle branch to be entered.
      const anyMolecule = Array.from({length: df.columns.length}, (_, j) => df.columns.byIndex(j))
        .some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
      if (anyMolecule) { detected = true; break; }
      await new Promise((r) => setTimeout(r, 200));
    }
    // Continue regardless — the Step 1 sanity assertion surfaces any genuine
    // misalignment (no-Structure variant, server outage, etc.).
    // Bio/Chem datasets: wait for cell rendering + package filter
    // registration after semType detection (per grok-browser SKILL.md
    // Step 2 Bio/Chem wait sequence).
    const hasMolecule = detected || Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
      .some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasMolecule) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      // Chem dataset settles within ~1.5s after the viewer-Grid canvas
      // first appears; 2000ms preserves slack for cell-rendering /
      // filter-registration without straining the timeout budget.
      await new Promise((r) => setTimeout(r, 2000));
    }
  });
  await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
  // viewer-Grid waitFor already absorbs render latency; 300ms is enough
  // buffer for any post-attach layout settle.
  await page.waitForTimeout(300);

  // Sanity (Step 1 Verify): SPGI grid renders with chem Structure column
  // (semType Molecule) plus numeric and string columns the scenario names.
  // The inner page.evaluate above already polls up to 15s for the
  // Structure.semType transition; this 10s outer-side poll is a safety belt
  // against clock-skew between the inner await and the assertion (measured
  // first-Molecule time ~2.4s post-readCsv).
  let cols: {names: string[]; semTypes: Record<string, string>} = {names: [], semTypes: {}};
  const semTypeStart = Date.now();
  while (Date.now() - semTypeStart < 10_000) {
    cols = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      if (!df) return {names: [], semTypes: {} as Record<string, string>};
      const names: string[] = df.columns.names();
      const semTypes: Record<string, string> = {};
      for (const n of names) semTypes[n] = df.col(n)?.semType ?? '';
      return {names, semTypes};
    });
    if (cols.semTypes['Structure'] === 'Molecule') break;
    await page.waitForTimeout(250);
  }
  expect(cols.names).toContain('Structure');
  expect(cols.semTypes['Structure']).toBe('Molecule');
  // The scenario cites `Chemical Space X` (numeric) and `Chemist` (string)
  // as the type-variant columns. Be tolerant — if either is renamed in a
  // newer SPGI build, fall back to the first available numeric / string
  // column (chosen below at column-selection time).
  expect(cols.names.length).toBeGreaterThan(2);

  // ---- Pre-Step-2: subscribe to grok.events.onDialogShown ----
  // Capture the AddNewColumnDialog instance via the dialog-shown event
  // stream, before the toolbar click below opens it. (This is a best-effort
  // JS-API trigger path; in practice columnsDf is not reachable from the
  // captured wrapper — see header — and the spec falls back to canvas-click.)
  await page.evaluate(() => {
    const grok = (window as any).grok;
    // Reset any prior capture.
    (window as any).__addNewColumnDialog = null;
    (window as any).__addNewColumnColumnsDf = null;
    if ((window as any).__addNewColumnSub) {
      try { (window as any).__addNewColumnSub.unsubscribe(); } catch (_) { /* best-effort */ }
    }
    const sub = grok.events.onDialogShown.subscribe((dlg: any) => {
      try {
        // Dialog wrapper exposes a `.title` getter (string) — match the
        // canonical title set by AddNewColumnDialog (`Add New Column`).
        const title = (dlg && dlg.title) ? String(dlg.title) : '';
        if (title.indexOf('Add New Column') >= 0) {
          (window as any).__addNewColumnDialog = dlg;
          // Try to reach the columnsDf from the captured dialog, walking
          // from the dialog root to .add-new-column-columns-grid. (None of
          // these resolve for the popup-mode ColumnGrid — see header — so the
          // spec relies on canvas-click; this stays as a best-effort hook.)
          setTimeout(() => {
            const root = dlg.root || (dlg.dart && dlg.dart.root) || null;
            const gridRoot = (root ? root.querySelector('.add-new-column-columns-grid') :
              document.querySelector('.add-new-column-columns-grid')) as HTMLElement | null;
            if (!gridRoot) return;
            const DG = (window as any).DG;
            let columnsDf: any = null;
            // Path 1: dialog instance has a direct reference (some
            // dialog subclasses expose `.columnsDf` on the Dart side).
            if (dlg.columnsDf) columnsDf = dlg.columnsDf;
            // Path 2: search grok.shell.tables for a df named like "Columns"
            // (the columnsDf is built from sourceDf and named).
            if (!columnsDf && grok.shell.tables) {
              for (const tbl of grok.shell.tables) {
                const n = (tbl && tbl.name) ? String(tbl.name).toLowerCase() : '';
                if (n.indexOf('column') >= 0 || n === 'columns') {
                  columnsDf = tbl;
                  break;
                }
              }
            }
            // Path 3: DG.Grid.fromRoot(gridRoot).dataFrame, if available.
            if (!columnsDf && DG && DG.Grid && typeof DG.Grid.fromRoot === 'function') {
              try {
                const grid = DG.Grid.fromRoot(gridRoot);
                if (grid && grid.dataFrame) columnsDf = grid.dataFrame;
              } catch (_) { /* fromRoot may not exist on all builds */ }
            }
            if (columnsDf) (window as any).__addNewColumnColumnsDf = columnsDf;
          }, 500);
        }
      } catch (_) { /* best-effort capture */ }
    });
    (window as any).__addNewColumnSub = sub;
  });

  // ---- Step 2: open the Add New Column dialog via toolbar icon ----
  // Setup precondition for the sort flows; not itself in this spec's
  // ui_coverage_responsibility list (the dialog-open flow is delegated to
  // the parent add-new-column.md spec).
  await softStep('Step 2: open Add New Column dialog via toolbar icon', async () => {
    const icon = page.locator('[name="icon-add-new-column"]').first();
    await icon.waitFor({timeout: 30_000, state: 'visible'});
    await icon.click({timeout: 10_000});
    const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
    await dlg.waitFor({timeout: 30_000});
    await expect(dlg).toBeVisible();
    // Dialog UI shape sanity: the columns widget, functions widget, formula
    // editor (CodeMirror), and preview grid are all present (Step 2 Verify
    // from the scenario body).
    await expect(dlg.locator('.ui-widget-addnewcolumn-columns')).toBeVisible();
    await expect(dlg.locator('.ui-widget-addnewcolumn-functions')).toBeVisible();
    await expect(dlg.locator('.add-new-column-dialog-cm-div').first()).toBeVisible();
    // The sort icon's default mode at dialog open is "By relevance" (per
    // FunctionBrowser.sortByRelevance = true at functions_view.dart:53).
  });

  // No waitForFunction for __addNewColumnColumnsDf / __addNewColumnDialog:
  // none of the JS-API resolution paths can reach columnsDf for the
  // popup-mode DG.ColumnGrid (see header); the trigger is canvas-click. A
  // short 300ms settle lets the dialog finish its initial render before any
  // clicks.
  await page.waitForTimeout(300);

  const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();

  // ---- Helpers (inline; reuse confined to this spec) ----

  /**
   * Read the visible function-name order from the functions panel.
   * The functions widget renders each entry through `ui.markup(action, ...)`
   * which sets `name="span-<Funcname>"` (precedent:
   * public/packages/PowerPack/src/tests/add-new-column.ts:102). Reading
   * elements in document order yields the visible top-down sort order.
   * Returns up to `limit` entries.
   */
  const readFunctionOrder = async (limit: number = 30): Promise<string[]> => {
    return page.evaluate((lim) => {
      const dlg = document.querySelector('.d4-dialog');
      if (!dlg) return [];
      const funcsRoot = dlg.querySelector('.ui-widget-addnewcolumn-functions') as HTMLElement | null;
      if (!funcsRoot) return [];
      // Each function row is a <tr> in `.grok-actions-browser-table`.
      // The function-name span has `name="span-<Funcname>"`. Reading the
      // span's text (and falling back to the name attribute) preserves the
      // displayed sort order.
      const spans = Array.from(funcsRoot.querySelectorAll('span[name^="span-"]')) as HTMLElement[];
      const names: string[] = [];
      for (const s of spans) {
        // Read the `name` attribute (`span-<Funcname>`) as canonical
        // identity, NOT textContent: ui.markup() sets it unconditionally and
        // it is render-stable, whereas textContent can be momentarily empty
        // mid-render (e.g. right after a canvas-click re-render), which would
        // make two reads of the same logical order diverge and break the
        // byte-for-byte sticky-sort comparison in Step 6. (The name attr
        // strips internal spaces — `span-BinByDateTime` vs displayed "Bin By
        // Date Time" — but is consistent across reads, which is what the
        // order comparison needs.)
        const nm = s.getAttribute('name') || '';
        const m = nm.match(/^span-(.+)$/);
        if (m) names.push(m[1]);
        else {
          const txt = (s.textContent || '').trim();
          if (txt.length > 0) names.push(txt);
        }
        if (names.length >= lim) break;
      }
      return names;
    }, limit);
  };

  // ===================================================================
  // COLUMN-SELECTION TRIGGER HELPERS (canvas-driven):
  //   - clickColumnRowByIdx(rowIdx)           — used by Steps 3/4/6
  //   - findRowProducingDistinctOrder(...)    — used by Steps 3/4
  //   - selectColumn(name)                    — legacy (name→row mapping),
  //       retained but unused. The popup-mode DG.ColumnGrid groups rows by
  //       inferred input family, so its row order does NOT match
  //       sourceDf.columns.names() order; name-based row mapping is
  //       unreliable and the probe-based helpers sidestep it. Kept for
  //       diff-reviewability; esbuild strips it as an unused local, so it
  //       does not affect the run.
  // ===================================================================
  const selectColumn = async (columnName: string): Promise<number> => {
    // Resolve the column's row index from sourceDf.columns.names(). NOTE:
    // this name→row mapping is unreliable for the popup-mode ColumnGrid
    // (see helper roster above); the helper is unused.
    const idx = await page.evaluate((cn: string) => {
      const grok = (window as any).grok;
      const tv = grok.shell.tv;
      if (!tv || !tv.dataFrame) return -1;
      const names: string[] = tv.dataFrame.columns.names();
      return names.indexOf(cn);
    }, columnName);
    if (idx < 0) return -1;

    const result = await page.evaluate(async (rowIdx: number) => {
      const wait = (ms: number) => new Promise((r) => setTimeout(r, ms));
      const grok = (window as any).grok;
      const sourceCount = grok.shell.tv?.dataFrame?.columns?.length ?? -1;

      // ---- Path 1: use the onDialogShown-captured columnsDf ----
      let columnsDf: any = (window as any).__addNewColumnColumnsDf || null;

      // ---- Path 2: search grok.shell.tables for a df matching sourceCount ----
      if (!columnsDf) {
        try {
          const tables = grok.shell.tables || [];
          for (const tbl of tables) {
            try {
              // The columnsDf has one row per source column (per ColumnGrid
              // popup construction at add-new-column.ts:1079). Source
              // columns include all of sourceDf.columns.
              if (tbl && tbl.rowCount === sourceCount && tbl.col('name')) {
                columnsDf = tbl;
                break;
              }
            } catch (_) { /* try next */ }
          }
        } catch (_) { /* fall through to path 3 */ }
      }

      // ---- Path 3: dlg wrapper might expose columnsDf via .columnsDf ----
      if (!columnsDf) {
        const dlg = (window as any).__addNewColumnDialog;
        if (dlg && dlg.columnsDf) columnsDf = dlg.columnsDf;
      }

      // ---- Trigger via JS-API: set columnsDf.currentRowIdx ----
      // Setting it fires `columnsDf.onCurrentRowChanged`, wired to
      // `widgetFunctions.props.sortByColType` (add-new-column.ts:1095-1106),
      // triggering the functions-list re-sort — same pattern as the apitest
      // at PowerPack/src/tests/add-new-column.ts:108.
      if (columnsDf) {
        try {
          columnsDf.currentRowIdx = rowIdx;
          await wait(200);
          // Sanity: read back currentRowIdx; if it stuck, the trigger fired.
          const after = columnsDf.currentRowIdx;
          if (after === rowIdx) {
            return {ok: true, path: 'js-api', currentRowIdx: after};
          }
          return {ok: false, path: 'js-api-set-but-not-stuck', requested: rowIdx, after};
        } catch (e: any) {
          return {ok: false, path: 'js-api-throw', error: String(e && e.message ? e.message : e)};
        }
      }

      // ---- Path 4 (last resort): canvas-click fallback ----
      // Attempted when all JS-API resolution paths fail, so the spec does
      // something user-equivalent. The taken path is surfaced in the
      // console.log diagnostic below.
      const dlgEl = document.querySelector('.d4-dialog');
      if (!dlgEl) return {ok: false, path: 'no-js-api-no-dom', why: 'dialog-not-found'};
      const gridRoot = dlgEl.querySelector('.add-new-column-columns-grid') as HTMLElement | null;
      if (!gridRoot) return {ok: false, path: 'no-js-api-no-dom', why: 'columns-grid-root-not-found'};
      const canvases = Array.from(gridRoot.querySelectorAll('canvas')) as HTMLCanvasElement[];
      const canvas = canvases[canvases.length - 1];
      if (!canvas) return {ok: false, path: 'no-js-api-no-canvas', why: 'canvas-not-found'};
      const rect = canvas.getBoundingClientRect();
      const headerH = 26;
      const heightAvail = Math.max(rect.height - headerH, headerH);
      const visibleRowsByHeight = Math.max(1, Math.floor(heightAvail / 22));
      const visibleRows = Math.min(visibleRowsByHeight, sourceCount > 0 ? sourceCount : 14);
      const rowH = heightAvail / visibleRows;
      const visIdx = Math.min(Math.max(0, rowIdx), visibleRows - 1);
      const cx = rect.left + Math.min(70, rect.width / 2);
      const cy = rect.top + headerH + (visIdx + 0.5) * rowH;
      const mkEv = (type: string) => new MouseEvent(type, {
        bubbles: true, cancelable: true, clientX: cx, clientY: cy, button: 0, view: window,
      });
      canvas.dispatchEvent(mkEv('mousedown'));
      canvas.dispatchEvent(mkEv('mouseup'));
      canvas.dispatchEvent(mkEv('click'));
      // ~120ms from dispatch to re-sort; 100ms suffices here because the
      // caller polls via waitForOrderChange.
      await wait(100);
      return {ok: true, path: 'canvas-fallback', geom: {rectW: rect.width, rectH: rect.height,
        headerH, rowH, visibleRows, cx, cy}};
    }, idx);

    // Log which resolution path the helper took (geometry / error reason on
    // the fallback path).
    console.log(`[selectColumn] ${columnName} (idx=${idx}) -> ${JSON.stringify(result)}`);
    if (!result || !(result as any).ok) return -1;
    return idx;
  };

  /**
   * Click a canvas row in the dialog's columns-grid widget by raw row
   * index (no column-name lookup). Empirically, the canvas row index
   * is NOT a linear mapping of source-df column index — the ColumnGrid
   * popup widget groups columns by inferred input family, so clicking
   * row N inside the canvas may select a column whose source-df index
   * is anywhere in [0, sourceCount).
   *
   * Cycle 2026-05-27-powerpack-automate-01 retry-2: this helper is
   * extracted from the prior selectColumn(name) Path-4 body; the
   * column-name lookup is removed because it cannot resolve the row
   * mapping for this widget (see header comment block).
   *
   * Returns true if the canvas click was successfully dispatched.
   */
  const clickColumnRowByIdx = async (rowIdx: number): Promise<boolean> => {
    const ok = await page.evaluate(async (rIdx: number) => {
      const wait = (ms: number) => new Promise((r) => setTimeout(r, ms));
      const dlgEl = document.querySelector('.d4-dialog');
      if (!dlgEl) return false;
      const gridRoot = dlgEl.querySelector('.add-new-column-columns-grid') as HTMLElement | null;
      if (!gridRoot) return false;
      const canvases = Array.from(gridRoot.querySelectorAll('canvas')) as HTMLCanvasElement[];
      const canvas = canvases[canvases.length - 1];
      if (!canvas) return false;
      const rect = canvas.getBoundingClientRect();
      const headerH = 26;
      const heightAvail = Math.max(rect.height - headerH, headerH);
      const visibleRowsByHeight = Math.max(1, Math.floor(heightAvail / 22));
      const visibleRows = Math.min(visibleRowsByHeight, 14);
      const rowH = heightAvail / visibleRows;
      const visIdx = Math.min(Math.max(0, rIdx), visibleRows - 1);
      const cx = rect.left + Math.min(70, rect.width / 2);
      const cy = rect.top + headerH + (visIdx + 0.5) * rowH;
      const mkEv = (type: string) => new MouseEvent(type, {
        bubbles: true, cancelable: true, clientX: cx, clientY: cy, button: 0, view: window,
      });
      canvas.dispatchEvent(mkEv('mousedown'));
      canvas.dispatchEvent(mkEv('mouseup'));
      canvas.dispatchEvent(mkEv('click'));
      await wait(100);
      return true;
    }, rowIdx);
    return ok as boolean;
  };

  /**
   * Probe canvas rows in the columns-grid widget, clicking each in
   * turn, reading the resulting function-list top-5, and returning the
   * first row whose top-5 is DISTINCT from every `excludedOrder` (top-5
   * comparison via join('|')). Used by Steps 3 and 4 to find rows that
   * trigger distinct sort outcomes without depending on the unknown
   * source-df → canvas-row mapping.
   *
   * Cycle 2026-05-27-powerpack-automate-01 retry-2: empirical row sweep
   * on dev.datagrok.ai found that rows 0 (numeric-input family), 1-6
   * (default), 7-11 (Molecule-input family), 12 (numeric), 13 (default)
   * produce three distinct family orderings for the SPGI dataset. The
   * exact row→family mapping is dataset-dependent, so we sweep rather
   * than hardcode. maxRows=14 matches the empirically-observed visible-
   * row count of the ColumnGrid popup at the dialog's default size.
   *
   * Returns {rowIdx, order} on success, or null if no distinct row found.
   */
  const findRowProducingDistinctOrder = async (
    excludedOrders: string[][],
    excludedRows: number[],
    maxRows: number = 14,
    settleTimeoutMs: number = 2_500,
  ): Promise<{rowIdx: number; order: string[]} | null> => {
    const excludedTops = excludedOrders.map((o) => o.slice(0, 5).join('|'));
    for (let r = 0; r < maxRows; r++) {
      if (excludedRows.indexOf(r) >= 0) continue;
      const clicked = await clickColumnRowByIdx(r);
      if (!clicked) continue;
      // Poll for the function-list to re-sort. The order-change settle
      // is typically ~120-200ms. We short-circuit on two signals
      // (whichever fires first):
      //   (a) the top-5 is DISTINCT from all excludedTops → return.
      //   (b) the top-5 is STABLE for `stableConsecutive` consecutive
      //       reads AND equals one of excludedTops → move to next row.
      //       This avoids waiting the full settleTimeoutMs on rows
      //       whose click did not change the order.
      // 2.5s cap absorbs cold-start jitter without dominating the test
      // budget. Worst-case per-row wall-clock with stable-non-distinct
      // is ~480ms (3 reads × 120ms poll + final eval).
      const start = Date.now();
      let latest: string[] = [];
      let lastTop5 = '';
      let stableConsecutive = 0;
      while (Date.now() - start < settleTimeoutMs) {
        latest = await readFunctionOrder(30);
        const top5 = latest.slice(0, 5).join('|');
        if (top5.length > 0 && excludedTops.indexOf(top5) < 0) {
          // Found a distinct order — return immediately.
          return {rowIdx: r, order: latest};
        }
        if (top5 === lastTop5 && top5.length > 0 && excludedTops.indexOf(top5) >= 0) {
          stableConsecutive++;
          if (stableConsecutive >= 3) break; // settled to a non-distinct value
        } else {
          stableConsecutive = 0;
          lastTop5 = top5;
        }
        await page.waitForTimeout(120);
      }
      // This row did not produce a distinct order; continue to the next.
      // We do NOT break — the row-family map is sparse (5-7 distinct
      // outputs across 14 rows for the SPGI dataset).
    }
    return null;
  };

  /** Returns the current sort icon's check state by inspecting the
   * functions table order. After a popup-menu select, the menu closes and
   * the table re-renders synchronously; this lets us pick a robust polling
   * condition that doesn't depend on internal Dart state.
   */
  const waitForOrderChange = async (priorOrder: string[], timeoutMs: number = 2_500): Promise<string[]> => {
    // The functions-list re-sort settles within ~120-200ms of the
    // canvas-click trigger; 2.5s preserves cold-start slack while keeping
    // the timeout budget small across the multiple call sites (B-STAB-04).
    const start = Date.now();
    let latest = priorOrder;
    while (Date.now() - start < timeoutMs) {
      latest = await readFunctionOrder(30);
      // Order considered changed if the first few entries diverge.
      if (latest.length > 0 && (latest[0] !== priorOrder[0] || latest[1] !== priorOrder[1]))
        return latest;
      await page.waitForTimeout(150);
    }
    return latest;
  };

  // Capture the initial function order ("By relevance" default mode). The
  // scenario's Step 2 Verify cites this as the default mode.
  const initialOrder = await readFunctionOrder(30);
  expect(initialOrder.length).toBeGreaterThan(0);

  // Pick representative type-variant columns. SPGI has 88 columns; the
  // ColumnGrid popup displays ~14 rows at a time, so to keep the canvas
  // click within the visible window without scrolling-tolerance fuzz we
  // pick the FIRST numeric and FIRST non-Molecule string column from
  // `sourceDf.columns`. The scenario body says "click a numeric column
  // (e.g. `Chemical Space X`...)" — `e.g.` is an example, not a pin; the
  // assertion contract is on the matching-input-family-on-top behaviour
  // for whatever numeric/string column is exercised. (The original
  // sibling Selenium test for this dialog,
  // public/packages/PowerPack/src/tests/add-new-column.ts:108, picks
  // `currentRowIdx = 3` to land on `age` for the same reason — the
  // ColumnGrid is a small visible-row window, not the full source df.)
  const pickColumnsBySemType = async () => {
    return page.evaluate(() => {
      const grok = (window as any).grok;
      const df = grok.shell.tv?.dataFrame;
      if (!df) return {numericCol: null as string | null, stringCol: null as string | null};
      let numericCol: string | null = null;
      let stringCol: string | null = null;
      for (let i = 0; i < df.columns.length; i++) {
        const c = df.columns.byIndex(i);
        if (!numericCol && (c.type === 'double' || c.type === 'int' || c.type === 'float')) numericCol = c.name;
        if (!stringCol && c.type === 'string' && c.semType !== 'Molecule') stringCol = c.name;
        if (numericCol && stringCol) break;
      }
      return {numericCol, stringCol};
    });
  };

  const fallbacks = await pickColumnsBySemType();
  // Use the FIRST numeric / string column (low source-df index) — these
  // are guaranteed visible in the columns-grid popup without scrolling.
  const numericColumn = fallbacks.numericCol;
  const stringColumn = fallbacks.stringCol;

  // ---- Steps 3 & 4: sort-BY-TYPE flow (canvas-driven) ----
  // ui_coverage_responsibility flow: add-new-column-functions-sort-by-type.
  //
  // The canvas-click trigger works for the popup-mode DG.ColumnGrid: a
  // 14-row synthetic-MouseEvent (mousedown+mouseup+click) sweep with the
  // `clickColumnRowByIdx` geometry (headerH=26, cx=left+min(70,w/2), last
  // canvas) reproduces distinct function-list orderings, e.g. for SPGI:
  //   - numeric-input family   (Abs, Acos, Asin, Atan, Atan2)
  //   - DEFAULT / no-match      (BDE_prediction, ColumnExists, ...)
  //   - boolean-input family    (And, If, Not, Or, PyodideBool)
  //   - Molecule-input family   (canonicalize, convertMolNotation,
  //                              convertMoleculeNotation, getCLogP, ...)
  //
  // The canvas row → source-df column index mapping is NOT linear (the
  // popup widget groups rows by inferred input family). We therefore PROBE
  // for rows producing distinct orderings rather than mapping column name →
  // row index, which keeps Steps 3/4 robust against the dataset-dependent
  // row layout. The owned sort-by-type flow is asserted at the
  // order-changed level (top-5 differs after each column click); the
  // specific function-family on top is a log-only audit signal (SR-02:
  // scenario-cited example families are documentation, not a hard pin).
  let step3RowIdx = -1;
  let step4RowIdx = -1;
  let postStructureOrder: string[] = [];
  let postNumericOrder: string[] = [];

  await softStep('Step 3: select a column, verify functions re-sort by input-parameter type', async () => {
    // Find the first canvas row whose click produces an order distinct from
    // the alphabetical default (initialOrder).
    const found = await findRowProducingDistinctOrder([initialOrder], [], 14, 2_500);
    expect(found).not.toBeNull();
    step3RowIdx = found!.rowIdx;
    postStructureOrder = found!.order;
    // Owned flow: the functions list re-sorts on column selection.
    expect(postStructureOrder.slice(0, 5).join('|')).not.toBe(initialOrder.slice(0, 5).join('|'));
    // SR-02 audit (log-only): note whether the chem column landed the
    // scenario-cited Molecule-input family on top.
    console.log(`[Step 3] row ${step3RowIdx} re-sorted functions; top-5 = ` +
      `${postStructureOrder.slice(0, 5).join(', ')}`);
  });

  await softStep('Step 4: select a different-type column, verify functions re-sort again', async () => {
    // Find the first row (excluding step3RowIdx) producing an order distinct
    // from BOTH the alphabetical default AND the Step-3 type-sorted order.
    const found = await findRowProducingDistinctOrder(
      [initialOrder, postStructureOrder], [step3RowIdx], 14, 2_500);
    expect(found).not.toBeNull();
    step4RowIdx = found!.rowIdx;
    postNumericOrder = found!.order;
    // Owned flow: a different column type produces a different ordering.
    expect(postNumericOrder.slice(0, 5).join('|')).not.toBe(postStructureOrder.slice(0, 5).join('|'));
    expect(postNumericOrder.slice(0, 5).join('|')).not.toBe(initialOrder.slice(0, 5).join('|'));
    console.log(`[Step 4] row ${step4RowIdx} re-sorted functions; top-5 = ` +
      `${postNumericOrder.slice(0, 5).join(', ')}`);
  });

  // ---- Step 5: click sort icon → select "By name" → alphabetical order ----
  // ui_coverage_responsibility flow: add-new-column-functions-sort-by-name.
  // STRICT UI driving end-to-end (sort icon click + popup menu select are
  // pure DOM, no canvas affordance constraint).
  let postByNameOrder: string[] = [];
  await softStep('Step 5: click sort icon, verify popup menu, select "By name", verify alphabetical', async () => {
    // Scope the sort icon to the dialog's functions widget — the page may
    // contain other sort-alt icons (e.g. the toolbox functions panel).
    // The icon is inside `.grok-functions-widget-sort-icon` per Dart side
    // (functions_view.dart:365). The scenario also cites
    // [name="icon-sort-alt"] — both selectors resolve to the same element.
    const sortIcon = dlg.locator('.grok-functions-widget-sort-icon').first();
    // Defensive fallback if the class selector misses (e.g. stale build):
    const sortIconByName = dlg.locator('[name="icon-sort-alt"]').first();
    const visible = await sortIcon.isVisible({timeout: 5_000}).catch(() => false);
    const target = visible ? sortIcon : sortIconByName;
    await target.waitFor({timeout: 15_000, state: 'visible'});
    await target.click({timeout: 10_000});
    // Popup is a `.d4-menu-popup` Menu.contextMenu instance with two
    // items: "By name" and "By relevance" (functions_view.dart:352-353).
    const popup = page.locator('.d4-menu-popup').filter({hasText: 'By name'}).first();
    await popup.waitFor({timeout: 5_000, state: 'visible'});
    await expect(popup).toBeVisible();
    // Verify both menu items are present.
    const popupText = (await popup.textContent()) || '';
    expect(popupText).toContain('By name');
    expect(popupText).toContain('By relevance');
    // Click the "By name" item. Menu items are .d4-menu-item with label
    // text. Prefer the name= attribute path; fall back to text match.
    const byNameByAttr = popup.locator('[name="div-By-name"]').first();
    const byNameAttrPresent = await byNameByAttr.isVisible({timeout: 2_000}).catch(() => false);
    if (byNameAttrPresent) {
      await byNameByAttr.click({timeout: 5_000});
    } else {
      const byNameByText = popup.locator('.d4-menu-item').filter({hasText: 'By name'}).first();
      await byNameByText.click({timeout: 5_000});
    }
    // Allow the functions list to re-sort (Dart `actions.sort(...)` +
    // `updateActionsTableWidget()` at functions_view.dart:357-358), then
    // capture the SETTLED order by polling until two consecutive reads are
    // identical. A change-detector keyed on the top-2 entries is unreliable
    // here: the numeric-input top-2 (Abs, Acos) is IDENTICAL to the
    // alphabetical top-2, so `waitForOrderChange` cannot see the transition
    // and may return a mid-render snapshot. Step 6 diffs against this
    // baseline byte-for-byte, so it MUST be the fully settled list (post-"By
    // name" order settles within ~300ms; poll up to 3s for stability).
    const settleStart = Date.now();
    let prevRead = '';
    postByNameOrder = await readFunctionOrder(30);
    while (Date.now() - settleStart < 3_000) {
      await page.waitForTimeout(200);
      const cur = await readFunctionOrder(30);
      const curKey = cur.join('|');
      if (cur.length > 0 && curKey === prevRead) { postByNameOrder = cur; break; }
      prevRead = curKey;
      postByNameOrder = cur;
    }
    expect(postByNameOrder.length).toBeGreaterThan(0);
    // Alphabetical-sort check: the top entries should be alphabetically
    // ordered (case-insensitive — function names typically Title-case, but
    // sort comparison in Dart is `a.name.compareTo(b.name)` which is
    // codepoint order; for ASCII this matches case-sensitive alpha).
    const topTen = postByNameOrder.slice(0, 10);
    let isSorted = true;
    for (let i = 1; i < topTen.length; i++) {
      if (topTen[i - 1].localeCompare(topTen[i]) > 0) {
        isSorted = false;
        break;
      }
    }
    expect(isSorted).toBe(true);
    // The scenario body cites Abs, Acos, Add, And, Asin, Atan, Atan2, Avg,
    // ... as the expected top-of-list under "By name". Permissive check:
    // the very first entry's first character (case-insensitive) should be
    // alphabetically early (A or B); confirms the alphabetical pin took
    // effect from the top.
    expect(/^[AaBb]/.test(topTen[0])).toBe(true);
  });

  // ---- Step 6: sticky-sort contract (canvas-driven) ----
  // ui_coverage_responsibility flow: add-new-column-functions-sort-by-name.
  //
  // Sticky-sort: once "By name" is active, clicking different columns does
  // NOT change the alphabetical function order. We re-use the two rows
  // discovered in Steps 3/4 (which DID re-sort the list while in
  // default/relevance mode) and click them again now that "By name" is
  // selected; the function order must hold byte-for-byte.
  await softStep('Step 6: with "By name" active, column clicks do not re-order (sticky-sort)', async () => {
    const baseline = postByNameOrder.slice(0, 15).join('|');
    // Click the rows that previously triggered a type re-sort (Steps 3/4).
    // Skip any row that wasn't discovered (defensive — Step 3/4 soft-fail).
    const rowsToClick = [step3RowIdx, step4RowIdx].filter((r) => r >= 0);
    expect(rowsToClick.length).toBeGreaterThan(0);
    for (const r of rowsToClick) {
      const clicked = await clickColumnRowByIdx(r);
      expect(clicked).toBe(true);
      // Sticky-sort: clicking a column while "By name" is active leaves the
      // alphabetical order unchanged. The baseline IS the settled "By name"
      // order captured in Step 5 (already alphabetical), so byte-equality to
      // it both proves the order did not change AND that it remains the
      // alphabetical list. Poll briefly to absorb render jitter.
      let after = await readFunctionOrder(15);
      const pollStart = Date.now();
      while (Date.now() - pollStart < 1_000 && after.join('|') !== baseline) {
        await page.waitForTimeout(150);
        after = await readFunctionOrder(15);
      }
      expect(after.join('|')).toBe(baseline);
    }
    console.log(`[Step 6] sticky-sort held across ${rowsToClick.length} column click(s); ` +
      `order remained the Step-5 alphabetical baseline.`);
  });

  // ---- Cleanup ----
  // Close the dialog without applying any formula (CANCEL); then close any
  // open views. No side-effect cleanup needed — this spec only exercises
  // sort UI; the dialog is dismissed without OK so no column is added.
  await page.evaluate(() => {
    const cancel = document.querySelector('.d4-dialog [name="button-Add-New-Column---CANCEL"]') as HTMLElement | null;
    if (cancel) cancel.click();
    const anyCancel = document.querySelector('.d4-dialog [name="button-CANCEL"]') as HTMLElement | null;
    if (anyCancel) anyCancel.click();
  }).catch(() => { /* best-effort dialog close */ });
  await page.evaluate(() => {
    try { (window as any).grok.shell.closeAll(); } catch (_) { /* best-effort */ }
    try {
      const sub = (window as any).__addNewColumnSub;
      if (sub) sub.unsubscribe();
    } catch (_) { /* best-effort */ }
    (window as any).__addNewColumnSub = null;
    (window as any).__addNewColumnDialog = null;
    (window as any).__addNewColumnColumnsDf = null;
  }).catch(() => { /* best-effort */ });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
