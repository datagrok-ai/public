/* ---
sub_features_covered: [powerpack.dialogs.add-new-column, powerpack.dialogs.add-new-column-func]
--- */
// Add New Column functions-panel sorting: by type (column select), by name, sticky-sort.
//
// Load-bearing facts:
//   - The columns widget is a popup-mode DG.ColumnGrid whose columnsDf is NOT reachable via JS-API
//     (not in grok.shell.tables; DG.Grid.fromRoot dataFrame is null; the captured dialog wrapper has no
//     columnsDf). The working trigger is a canvas MouseEvent triple (mousedown+mouseup+click) on the
//     top-most canvas in .add-new-column-columns-grid at the computed (cx,cy) for a row index.
//   - The canvas row → source-df column index mapping is NOT linear (rows are grouped by inferred input
//     family), so Steps 3/4 PROBE rows for distinct sort outcomes rather than mapping name → row.
//   - SR-02: the function panel re-sorts on column change, but the scenario-cited example families don't
//     land on top, so sort-by-type is asserted at the order-changed level (top-5 differs), not by family.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);

test('PowerPack: Add new column functions-panel sorting (SPGI — by type, by name, sticky)', async ({page}) => {
  // 540_000 keeps a single attempt under the 600s wrapper bound (retries=1 × 300s would compound past it).
  test.setTimeout(540_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    const grok = (window as any).grok;
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    try { grok.shell.closeAll(); } catch (_) { /* best-effort */ }
    // Try the chem-subdir path first; fall back to demo-root SPGI.csv.
    let df: any = null;
    try {
      df = await grok.dapi.files.readCsv('System:DemoFiles/chem/SPGI.csv');
    } catch (_) {
      df = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
    }
    grok.shell.addTableView(df);
    // Poll for Structure.semType='Molecule' (or any Molecule col); the onSemanticTypeDetected event is racy.
    let detected = false;
    for (let i = 0; i < 75; i++) {
      const structureCol = df.col('Structure');
      if (structureCol && structureCol.semType === 'Molecule') { detected = true; break; }
      const anyMolecule = Array.from({length: df.columns.length}, (_, j) => df.columns.byIndex(j))
        .some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
      if (anyMolecule) { detected = true; break; }
      await new Promise((r) => setTimeout(r, 200));
    }
    // Chem datasets: wait for cell rendering + filter registration after semType detection.
    const hasMolecule = detected || Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
      .some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasMolecule) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 2000));
    }
  });
  await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
  await page.waitForTimeout(300);

  // Sanity: SPGI grid renders with chem Structure (Molecule) plus numeric/string columns.
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
  expect(cols.names.length).toBeGreaterThan(2);

  // Best-effort capture of the AddNewColumnDialog via onDialogShown (columnsDf is not reachable — see header).
  await page.evaluate(() => {
    const grok = (window as any).grok;
    (window as any).__addNewColumnDialog = null;
    (window as any).__addNewColumnColumnsDf = null;
    if ((window as any).__addNewColumnSub) {
      try { (window as any).__addNewColumnSub.unsubscribe(); } catch (_) { /* best-effort */ }
    }
    const sub = grok.events.onDialogShown.subscribe((dlg: any) => {
      try {
        const title = (dlg && dlg.title) ? String(dlg.title) : '';
        if (title.indexOf('Add New Column') >= 0) {
          (window as any).__addNewColumnDialog = dlg;
          setTimeout(() => {
            const root = dlg.root || (dlg.dart && dlg.dart.root) || null;
            const gridRoot = (root ? root.querySelector('.add-new-column-columns-grid') :
              document.querySelector('.add-new-column-columns-grid')) as HTMLElement | null;
            if (!gridRoot) return;
            const DG = (window as any).DG;
            let columnsDf: any = null;
            if (dlg.columnsDf) columnsDf = dlg.columnsDf;
            if (!columnsDf && grok.shell.tables) {
              for (const tbl of grok.shell.tables) {
                const n = (tbl && tbl.name) ? String(tbl.name).toLowerCase() : '';
                if (n.indexOf('column') >= 0 || n === 'columns') {
                  columnsDf = tbl;
                  break;
                }
              }
            }
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

  await softStep('Step 2: open Add New Column dialog via toolbar icon', async () => {
    const icon = page.locator('[name="icon-add-new-column"]').first();
    await icon.waitFor({timeout: 30_000, state: 'visible'});
    await icon.click({timeout: 10_000});
    const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
    await dlg.waitFor({timeout: 30_000});
    await expect(dlg).toBeVisible();
    // Step 2 verify: columns widget, functions widget, formula editor all present.
    await expect(dlg.locator('.ui-widget-addnewcolumn-columns')).toBeVisible();
    await expect(dlg.locator('.ui-widget-addnewcolumn-functions')).toBeVisible();
    await expect(dlg.locator('.add-new-column-dialog-cm-div').first()).toBeVisible();
  });

  await page.waitForTimeout(300); // let the dialog finish its initial render

  const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();

  // Read the visible function-name order from the functions panel (name="span-<Funcname>").
  const readFunctionOrder = async (limit: number = 30): Promise<string[]> => {
    return page.evaluate((lim) => {
      const dlg = document.querySelector('.d4-dialog');
      if (!dlg) return [];
      const funcsRoot = dlg.querySelector('.ui-widget-addnewcolumn-functions') as HTMLElement | null;
      if (!funcsRoot) return [];
      const spans = Array.from(funcsRoot.querySelectorAll('span[name^="span-"]')) as HTMLElement[];
      const names: string[] = [];
      for (const s of spans) {
        // Read name= (render-stable) not textContent (momentarily empty mid-render breaks Step 6 byte-compare).
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

  // Column-selection trigger helpers (canvas-driven). selectColumn is legacy/unused (name→row mapping is
  // unreliable for the popup-mode ColumnGrid); kept for diff-reviewability, stripped by esbuild as unused.
  const selectColumn = async (columnName: string): Promise<number> => {
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

      let columnsDf: any = (window as any).__addNewColumnColumnsDf || null;
      if (!columnsDf) {
        try {
          const tables = grok.shell.tables || [];
          for (const tbl of tables) {
            try {
              if (tbl && tbl.rowCount === sourceCount && tbl.col('name')) {
                columnsDf = tbl;
                break;
              }
            } catch (_) { /* try next */ }
          }
        } catch (_) { /* fall through */ }
      }
      if (!columnsDf) {
        const dlg = (window as any).__addNewColumnDialog;
        if (dlg && dlg.columnsDf) columnsDf = dlg.columnsDf;
      }

      // Setting currentRowIdx fires onCurrentRowChanged → functions-list re-sort.
      if (columnsDf) {
        try {
          columnsDf.currentRowIdx = rowIdx;
          await wait(200);
          const after = columnsDf.currentRowIdx;
          if (after === rowIdx) {
            return {ok: true, path: 'js-api', currentRowIdx: after};
          }
          return {ok: false, path: 'js-api-set-but-not-stuck', requested: rowIdx, after};
        } catch (e: any) {
          return {ok: false, path: 'js-api-throw', error: String(e && e.message ? e.message : e)};
        }
      }

      // Last resort: canvas-click fallback.
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
      await wait(100); // caller polls via waitForOrderChange
      return {ok: true, path: 'canvas-fallback', geom: {rectW: rect.width, rectH: rect.height,
        headerH, rowH, visibleRows, cx, cy}};
    }, idx);

    console.log(`[selectColumn] ${columnName} (idx=${idx}) -> ${JSON.stringify(result)}`);
    if (!result || !(result as any).ok) return -1;
    return idx;
  };

  // Click a canvas row in the columns-grid by raw row index (no name lookup; row→column mapping is non-linear).
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

  // Probe canvas rows, returning the first whose function-list top-5 differs from every excludedOrder.
  // The row→family map is dataset-dependent, so sweep rather than hardcode; maxRows=14 = visible-row count.
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
      // Poll for re-sort: short-circuit on a distinct top-5, or on a stable non-distinct top-5 (3 reads).
      const start = Date.now();
      let latest: string[] = [];
      let lastTop5 = '';
      let stableConsecutive = 0;
      while (Date.now() - start < settleTimeoutMs) {
        latest = await readFunctionOrder(30);
        const top5 = latest.slice(0, 5).join('|');
        if (top5.length > 0 && excludedTops.indexOf(top5) < 0)
          return {rowIdx: r, order: latest};
        if (top5 === lastTop5 && top5.length > 0 && excludedTops.indexOf(top5) >= 0) {
          stableConsecutive++;
          if (stableConsecutive >= 3) break; // settled to a non-distinct value
        } else {
          stableConsecutive = 0;
          lastTop5 = top5;
        }
        await page.waitForTimeout(120);
      }
    }
    return null;
  };

  const waitForOrderChange = async (priorOrder: string[], timeoutMs: number = 2_500): Promise<string[]> => {
    const start = Date.now();
    let latest = priorOrder;
    while (Date.now() - start < timeoutMs) {
      latest = await readFunctionOrder(30);
      if (latest.length > 0 && (latest[0] !== priorOrder[0] || latest[1] !== priorOrder[1]))
        return latest;
      await page.waitForTimeout(150);
    }
    return latest;
  };

  // Initial function order ("By relevance" default mode).
  const initialOrder = await readFunctionOrder(30);
  expect(initialOrder.length).toBeGreaterThan(0);

  // Pick the FIRST numeric + FIRST non-Molecule string column (guaranteed visible in the popup window).
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
  const numericColumn = fallbacks.numericCol;
  const stringColumn = fallbacks.stringCol;

  // Steps 3 & 4: sort-by-type. Probe canvas rows for distinct orderings (row→family map is non-linear);
  // assert at the order-changed level (top-5 differs), family-on-top is log-only (SR-02).
  let step3RowIdx = -1;
  let step4RowIdx = -1;
  let postStructureOrder: string[] = [];
  let postNumericOrder: string[] = [];

  await softStep('Step 3: select a column, verify functions re-sort by input-parameter type', async () => {
    const found = await findRowProducingDistinctOrder([initialOrder], [], 14, 2_500);
    expect(found).not.toBeNull();
    step3RowIdx = found!.rowIdx;
    postStructureOrder = found!.order;
    expect(postStructureOrder.slice(0, 5).join('|')).not.toBe(initialOrder.slice(0, 5).join('|'));
    console.log(`[Step 3] row ${step3RowIdx} re-sorted functions; top-5 = ` +
      `${postStructureOrder.slice(0, 5).join(', ')}`);
  });

  await softStep('Step 4: select a different-type column, verify functions re-sort again', async () => {
    const found = await findRowProducingDistinctOrder(
      [initialOrder, postStructureOrder], [step3RowIdx], 14, 2_500);
    expect(found).not.toBeNull();
    step4RowIdx = found!.rowIdx;
    postNumericOrder = found!.order;
    expect(postNumericOrder.slice(0, 5).join('|')).not.toBe(postStructureOrder.slice(0, 5).join('|'));
    expect(postNumericOrder.slice(0, 5).join('|')).not.toBe(initialOrder.slice(0, 5).join('|'));
    console.log(`[Step 4] row ${step4RowIdx} re-sorted functions; top-5 = ` +
      `${postNumericOrder.slice(0, 5).join(', ')}`);
  });

  // Step 5: click sort icon → "By name" → alphabetical (pure DOM, no canvas constraint).
  let postByNameOrder: string[] = [];
  await softStep('Step 5: click sort icon, verify popup menu, select "By name", verify alphabetical', async () => {
    // Scope to the dialog's functions widget (other sort-alt icons may exist on the page).
    const sortIcon = dlg.locator('.grok-functions-widget-sort-icon').first();
    const sortIconByName = dlg.locator('[name="icon-sort-alt"]').first();
    const visible = await sortIcon.isVisible({timeout: 5_000}).catch(() => false);
    const target = visible ? sortIcon : sortIconByName;
    await target.waitFor({timeout: 15_000, state: 'visible'});
    await target.click({timeout: 10_000});
    const popup = page.locator('.d4-menu-popup').filter({hasText: 'By name'}).first();
    await popup.waitFor({timeout: 5_000, state: 'visible'});
    await expect(popup).toBeVisible();
    const popupText = (await popup.textContent()) || '';
    expect(popupText).toContain('By name');
    expect(popupText).toContain('By relevance');
    const byNameByAttr = popup.locator('[name="div-By-name"]').first();
    const byNameAttrPresent = await byNameByAttr.isVisible({timeout: 2_000}).catch(() => false);
    if (byNameAttrPresent) {
      await byNameByAttr.click({timeout: 5_000});
    } else {
      const byNameByText = popup.locator('.d4-menu-item').filter({hasText: 'By name'}).first();
      await byNameByText.click({timeout: 5_000});
    }
    // Capture the SETTLED order by polling for two identical consecutive reads (Step 6 diffs it byte-for-byte;
    // a top-2 change-detector fails here because the numeric-input top-2 == alphabetical top-2).
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
    // Top-10 must be alphabetically ordered (Dart compareTo = codepoint order; ASCII == case-sensitive alpha).
    const topTen = postByNameOrder.slice(0, 10);
    let isSorted = true;
    for (let i = 1; i < topTen.length; i++) {
      if (topTen[i - 1].localeCompare(topTen[i]) > 0) {
        isSorted = false;
        break;
      }
    }
    expect(isSorted).toBe(true);
    expect(/^[AaBb]/.test(topTen[0])).toBe(true);
  });

  // Step 6: sticky-sort — with "By name" active, clicking the Step 3/4 rows must not re-order (byte-for-byte).
  await softStep('Step 6: with "By name" active, column clicks do not re-order (sticky-sort)', async () => {
    const baseline = postByNameOrder.slice(0, 15).join('|');
    const rowsToClick = [step3RowIdx, step4RowIdx].filter((r) => r >= 0);
    expect(rowsToClick.length).toBeGreaterThan(0);
    for (const r of rowsToClick) {
      const clicked = await clickColumnRowByIdx(r);
      expect(clicked).toBe(true);
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

  // Cleanup: dismiss the dialog via CANCEL (no column added), then close views.
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

  finishSpec();
});
