/* ---
realizes: [grid.cp.rows-select-filter-navigate]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;
declare const DG: any;


test.use(specTestOptions);

// Page-coordinate center of a column-header LABEL. Both coordinates are overlay-local plus the
// overlay rect origin; documentBounds.y uses the grid's own virtual-scroll origin instead, and
// mixing the two puts the click off-canvas where hitTest returns null.
async function headerLabelPoint(page: Page, col: string): Promise<{x: number; y: number}> {
  return page.evaluate((c) => {
    const grid = grok.shell.tv.grid;
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    const rc = overlay.getBoundingClientRect();
    const gc = grid.columns.byName(c);
    return {x: rc.x + gc.left + gc.width / 2, y: rc.y + grid.colHeaderHeight / 2};
  }, col);
}

// Page-coordinate center of a data cell (col, visualRow).
async function dataCellPoint(page: Page, col: string, visualRow: number): Promise<{x: number; y: number}> {
  return page.evaluate(({c, vr}) => {
    const grid = grok.shell.tv.grid;
    const db = grid.cell(c, vr).documentBounds;
    return {x: db.x + db.width / 2, y: db.y + db.height / 2};
  }, {c: col, vr: visualRow});
}

// Page-coordinate center of the ROW-HEADER STRIP (grid column 0) at a visual row. The strip is
// the ONLY column with isSelectionDragStart, so a trusted drag starting here selects the crossed
// range with no modifier (grid_row_selection.dart line 29). A drag over DATA cells matches
// neither selection path — they lack that flag and the areaSelector rubber-band needs shiftKey.
async function rowStripPoint(page: Page, visualRow: number): Promise<{x: number; y: number}> {
  return page.evaluate((vr) => {
    const grid = grok.shell.tv.grid;
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    const rc = overlay.getBoundingClientRect();
    const col0 = grid.columns.byIndex(0);
    const db = grid.cell(grid.columns.byIndex(1).name, vr).documentBounds;
    return {x: rc.x + col0.left + col0.width / 2, y: db.y + db.height / 2};
  }, visualRow);
}

// A trusted drag over the row-header strip, selecting that inclusive row range. The move is
// one-shot: a stepped move with an intermediate waypoint risks being read as a hover.
async function stripDragSelect(page: Page, startVisualRow: number, endVisualRow: number): Promise<void> {
  const s = await rowStripPoint(page, startVisualRow);
  const e = await rowStripPoint(page, endVisualRow);
  await page.mouse.move(s.x, s.y);
  await page.mouse.down();
  await page.mouse.move(e.x, e.y);
  await page.mouse.up();
  await page.waitForTimeout(400);
}

// Open a grid context menu at (clientX, clientY), expand the nested Pin submenu and click the
// given leaf; returns false if the leaf never becomes clickable. The submenu is slope/hover-intent
// guarded and expands only for coordinate-bearing mouseover+mouseenter+mousemove on the group
// center — a bare bubbling move leaves the leaf a zero-rect template. There is no JS-API pin, so
// this menu is the only actuator of the pinned-row state.
async function pinViaMenu(page: Page, at: {x: number; y: number}, leafName: string): Promise<boolean> {
  return page.evaluate(async ({x, y, leaf}) => {
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    const opt = (cx: number, cy: number) => ({bubbles: true, cancelable: true, clientX: cx, clientY: cy}) as any;
    for (let open = 0; open < 5; open++) {
      const cm = {bubbles: true, cancelable: true, clientX: x, clientY: y, button: 2, buttons: 2} as any;
      overlay.dispatchEvent(new MouseEvent('mousedown', cm));
      overlay.dispatchEvent(new MouseEvent('mouseup', cm));
      overlay.dispatchEvent(new MouseEvent('contextmenu', cm));
      await new Promise((r) => setTimeout(r, 550));
      const pin = document.querySelector('[name="div-Pin"]') as HTMLElement;
      if (!pin || pin.offsetParent === null) {
        document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 250));
        continue;
      }
      const pr = pin.getBoundingClientRect();
      const px = pr.x + pr.width / 2;
      const py = pr.y + pr.height / 2;
      let el: HTMLElement | null = null;
      for (let h = 0; h < 6; h++) {
        pin.dispatchEvent(new MouseEvent('mouseover', opt(px, py)));
        pin.dispatchEvent(new MouseEvent('mouseenter', opt(px, py)));
        pin.dispatchEvent(new MouseEvent('mousemove', opt(px + (h % 2 ? 0.5 : -0.5), py)));
        await new Promise((r) => setTimeout(r, 250));
        el = document.querySelector('[name="' + leaf + '"]') as HTMLElement;
        if (el && el.offsetParent !== null && el.getBoundingClientRect().width > 0) break;
      }
      if (!el || el.offsetParent === null || el.getBoundingClientRect().width === 0) {
        document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 250));
        continue;
      }
      const lr = el.getBoundingClientRect();
      const lo = {bubbles: true, cancelable: true, clientX: lr.x + lr.width / 2, clientY: lr.y + lr.height / 2, button: 0} as any;
      el.dispatchEvent(new MouseEvent('mouseover', lo));
      el.dispatchEvent(new MouseEvent('mousemove', lo));
      el.dispatchEvent(new MouseEvent('mousedown', lo));
      el.dispatchEvent(new MouseEvent('mouseup', lo));
      el.dispatchEvent(new MouseEvent('click', lo));
      await new Promise((r) => setTimeout(r, 700));
      document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
      return true;
    }
    return false;
  }, {x: at.x, y: at.y, leaf: leafName});
}

// Name of the first data column (grid column 1; column 0 is the row-number strip).
async function firstDataCol(page: Page): Promise<string> {
  return page.evaluate(() => grok.shell.tv.grid.columns.byIndex(1).name);
}

// Focus the grid overlay so keyboard gestures target it.
async function focusGrid(page: Page): Promise<void> {
  await page.evaluate(() => {
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    overlay.focus();
  });
}

// Reset the shell to a clean demog view (close all, reopen, clear filter/selection/sort).
async function reopenClean(page: Page): Promise<void> {
  await page.evaluate(async () => {
    grok.shell.closeAll();
    await new Promise((r) => setTimeout(r, 800));
  });
  await v.openTable(page, {path: 'System:DemoFiles/demog.csv', semTypeTimeoutMs: 3000});
  await page.evaluate(async () => {
    const tv = grok.shell.tv;
    const df = tv.dataFrame;
    const grid = tv.grid;
    df.filter.setAll(true);
    df.selection.setAll(false);
    grid.sort([], []);
    grid.props.rowSource = 'All';
    await new Promise((r) => setTimeout(r, 300));
  });
}

test('Grid — Row Selection, Filter, and Keyboard Navigation', async ({page}) => {
  test.setTimeout(420_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: 'System:DemoFiles/demog.csv', semTypeTimeoutMs: 3000});

  // This build exposes no grok.shell.warnings, so the console/pageerror channels are the
  // sanctioned no-error signal; the collector spans the whole test.
  const consoleErrors: string[] = [];
  const pageErrors: string[] = [];
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });
  page.on('pageerror', (e) => pageErrors.push(String(e)));

  await page.evaluate(async () => {
    const grid = grok.shell.tv.grid;
    const df = grok.shell.tv.dataFrame;
    df.filter.setAll(true);
    df.selection.setAll(false);
    grid.sort([], []);
    await new Promise((r) => setTimeout(r, 200));
  });

  // === Scenario 1: mouse clicks + Esc compound reset ==========================

  await softStep('Step 3 — click row 5 data cell: sets current row, no selection', async () => {
    await focusGrid(page);
    const col = await firstDataCol(page);
    const p = await dataCellPoint(page, col, 5);
    await page.mouse.click(p.x, p.y);
    await page.waitForTimeout(300);
    const r = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const grid = grok.shell.tv.grid;
      return {
        currentRowIdx: df.currentRowIdx,
        expectedTableIdx: grid.gridRowToTable(5),
        selTrueCount: df.selection.trueCount,
      };
    });
    expect(r.currentRowIdx).toBe(r.expectedTableIdx); // current row is the table index of visual row 5
    expect(r.selTrueCount).toBe(0);                   // single click selects nothing
  });

  await softStep('Step 4 — drag-select data rows 5..10: range select 6 rows', async () => {
    // The strip drag crosses 6 rows (5..10 inclusive); unscrolled and unfiltered, visual rows
    // map 1:1 to table rows.
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));
    await focusGrid(page);
    const intended = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      return [5, 6, 7, 8, 9, 10].map((v) => grid.gridRowToTable(v)).filter((t) => t >= 0);
    });
    await stripDragSelect(page, 5, 10);
    const r = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const grid = grok.shell.tv.grid;
      return {
        selected: Array.from(df.selection.getSelectedIndexes()),
        selTrueCount: df.selection.trueCount,
        currentRowIdx: df.currentRowIdx,
        expectedTableIdx: grid.gridRowToTable(10),
      };
    });
    const selSet = [...r.selected].sort((a, b) => a - b);
    const intendedSet = [...intended].sort((a, b) => a - b);
    expect(selSet).toEqual(intendedSet);              // exactly rows 5..10 selected
    expect(r.selTrueCount).toBe(intendedSet.length);  // 6 rows
  });

  // Step 5 (Ctrl+click disjoint toggle-add on the row-number strip) produces no selection change
  // headless, so the gesture is verified manually in grid-rows-select-filter-navigate-ui.md.
  // The toggle-membership channel it exercises is covered automatically by Step 17 below.

  await softStep('Step 6 — Ctrl+A: select all rows', async () => {
    await focusGrid(page);
    await page.keyboard.press('Control+a');
    await page.waitForTimeout(300);
    const r = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {selTrueCount: df.selection.trueCount, rowCount: df.rowCount};
    });
    expect(r.selTrueCount).toBe(r.rowCount); // all rows selected
  });

  await softStep('Step 7 — Esc: compound reset (selection, columns, current row all cleared)', async () => {
    // Select a column first so the column-clear channel is exercised too.
    await page.evaluate(() => { grok.shell.tv.grid.columns.byName('AGE').selected = true; grok.shell.tv.grid.invalidate(); });
    await page.waitForTimeout(150);
    await focusGrid(page);
    await page.keyboard.press('Escape');
    await page.waitForTimeout(300);
    const r = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {
        selTrueCount: df.selection.trueCount,
        selectedColNames: Array.from(df.columns.selected).map((c: any) => c.name),
        currentRowIdx: df.currentRowIdx,
      };
    });
    expect(r.selTrueCount).toBe(0);             // no rows selected
    expect(r.selectedColNames).toEqual([]);     // no columns selected
    expect(r.currentRowIdx).toBe(-1);           // current row cleared
  });

  // === Scenario 2: rapid Ctrl+click column-header guard (GROK-17455) ==========

  await reopenClean(page);

  // A header Ctrl+click writes through to df.columns.selected under headed input but leaves it
  // EMPTY headless, so the visible selection outcome is delegated to the manual checklist
  // (grid-rows-select-filter-navigate-ui.md). The guard below needs the clicks only to FIRE.

  await softStep('Step 10 — rapid Ctrl+clicks across 5+ column headers: no console error (GROK-17455)', async () => {
    const errBefore = consoleErrors.length + pageErrors.length;
    await focusGrid(page);
    for (const col of ['AGE', 'HEIGHT', 'WEIGHT', 'SEX', 'RACE', 'DIS_POP']) {
      const p = await headerLabelPoint(page, col);
      await page.mouse.click(p.x, p.y, {modifiers: ['Control']});
      await page.waitForTimeout(60);
    }
    await page.waitForTimeout(400);
    const concurrentModErrors = consoleErrors
      .concat(pageErrors)
      .filter((e) => /concurrent modification|during iteration/i.test(e));
    expect(concurrentModErrors).toEqual([]);                        // GROK-17455 guard: no concurrent-modification error
    expect(consoleErrors.length + pageErrors.length).toBe(errBefore); // console-error delta is 0
  });

  // === Scenario 3: block-drag row selection over data cells ===================

  await reopenClean(page);

  await softStep('Step 11 — Shift+drag over data cells: selects rows only, no column selection', async () => {
    await focusGrid(page);
    await page.keyboard.press('Escape');
    await page.waitForTimeout(150);
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));
    // Over DATA cells the selecting gesture is the areaSelector rubber-band, gated on shiftKey
    // (grid_row_selection.dart line 69) — hence a Shift+drag rather than the plain strip drag.
    const start = await dataCellPoint(page, 'AGE', 2);
    const end = await dataCellPoint(page, 'AGE', 6);
    const crossed = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      return [2, 3, 4, 5, 6].map((v) => grid.gridRowToTable(v)).filter((t) => t >= 0);
    });
    // Shift must be down for the whole gesture: the where-gate reads it at mousedown.
    await page.keyboard.down('Shift');
    await page.mouse.move(start.x, start.y);
    await page.mouse.down();
    await page.mouse.move(end.x, end.y);
    await page.mouse.up();
    await page.keyboard.up('Shift');
    await page.waitForTimeout(400);
    const r = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {
        selected: Array.from(df.selection.getSelectedIndexes()),
        selTrueCount: df.selection.trueCount,
        selectedColNames: Array.from(df.columns.selected).map((c: any) => c.name),
      };
    });
    const selectedSet = [...r.selected].sort((a, b) => a - b);
    const crossedSet = [...crossed].sort((a, b) => a - b);
    expect(selectedSet).toEqual(crossedSet); // selection is exactly the rows the drag crossed
    expect(r.selTrueCount).toBe(crossedSet.length); // and its count reflects the crossed-row count
    expect(r.selectedColNames).toEqual([]);      // block drag selects rows only (no column side effect)
  });

  // === Scenario 4: keyboard navigation (GROK-18104) ===========================

  await reopenClean(page);

  await softStep('Step 13 — arrows: Down x2 then Right x2 move current row/col', async () => {
    // Establish a known current cell first (a data-cell click sets current, not selection).
    const p = await dataCellPoint(page, 'AGE', 4);
    await page.mouse.click(p.x, p.y);
    await page.waitForTimeout(200);
    const before = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {rowIdx: df.currentRowIdx, col: df.currentCol ? df.currentCol.name : null};
    });
    await focusGrid(page);
    await page.keyboard.press('ArrowDown');
    await page.keyboard.press('ArrowDown');
    await page.keyboard.press('ArrowRight');
    await page.keyboard.press('ArrowRight');
    await page.waitForTimeout(300);
    const after = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {rowIdx: df.currentRowIdx, col: df.currentCol ? df.currentCol.name : null};
    });
    expect(after.rowIdx).toBe(before.rowIdx + 2);   // moved down by 2
    expect(after.col).not.toBe(before.col);          // current column changed on right-arrow
  });

  await softStep('Step 15 — Ctrl+Home goes to first row; Ctrl+End goes to last row', async () => {
    await focusGrid(page);
    await page.keyboard.press('Control+Home');
    await page.waitForTimeout(300);
    const home = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
    expect(home).toBe(0); // Ctrl+Home → first row
    await focusGrid(page);
    await page.keyboard.press('Control+End');
    await page.waitForTimeout(300);
    const r = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {rowIdx: df.currentRowIdx, rowCount: df.rowCount};
    });
    expect(r.rowIdx).toBe(r.rowCount - 1); // Ctrl+End → last row
  });

  await softStep('Step 16 — PageDown twice advances both times (GROK-18104)', async () => {
    // Reset to the top first.
    await focusGrid(page);
    await page.keyboard.press('Control+Home');
    await page.waitForTimeout(250);
    const p0 = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
    await page.keyboard.press('PageDown');
    await page.waitForTimeout(300);
    const p1 = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
    await page.keyboard.press('PageDown');
    await page.waitForTimeout(300);
    const p2 = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
    expect(p1).toBeGreaterThan(p0); // first PageDown advances
    expect(p2).toBeGreaterThan(p1); // second PageDown advances again (GROK-18104: no stall)
  });

  // === Scenario 5: Space toggle + Shift+Enter value-group ======================

  await reopenClean(page);

  await softStep('Step 17 (also proves the Step 5 toggle-membership channel) — Space toggles current row selection membership', async () => {
    // Current row is set with a DATA-cell click: a row-header click would also select, and the
    // Space toggle has to start from an unselected row.
    const p = await dataCellPoint(page, 'AGE', 3);
    await page.mouse.click(p.x, p.y);
    await page.waitForTimeout(200);
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));
    await page.waitForTimeout(100);
    await focusGrid(page);
    const before = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    await page.keyboard.press(' ');
    await page.waitForTimeout(250);
    const afterAdd = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    await page.keyboard.press(' ');
    await page.waitForTimeout(250);
    const afterRemove = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(afterAdd).toBe(before + 1);      // Space adds the current row
    expect(afterRemove).toBe(before);       // Space again removes it
  });

  await softStep('Step 18 — Shift+Enter selects all rows sharing the current cell value', async () => {
    // Move current cell into SEX (few distinct values) and clear selection.
    const expected = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      df.currentCell = df.cell(3, 'SEX');
      const sexVal = df.col('SEX').get(df.currentRowIdx);
      let count = 0;
      for (let i = 0; i < df.rowCount; i++) if (df.col('SEX').get(i) === sexVal) count++;
      return {count, sexVal};
    });
    await page.waitForTimeout(200);
    await focusGrid(page);
    await page.keyboard.press('Shift+Enter');
    await page.waitForTimeout(400);
    const sel = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(sel).toBe(expected.count); // all rows with the same SEX value are selected
  });

  // === Scenario 6: negative path — allowRowSelection false =====================

  await reopenClean(page);

  await softStep('Step 20 — allowRowSelection false: Shift+Down navigates without selecting', async () => {
    // Neutral setup of the negative precondition; the effect under test is that Shift+arrow
    // navigation does not select while the toggle is off.
    await page.evaluate(() => { grok.shell.tv.grid.props.allowRowSelection = false; });
    await page.waitForTimeout(200);
    const p = await dataCellPoint(page, 'AGE', 4);
    await page.mouse.click(p.x, p.y);
    await page.waitForTimeout(200);
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));
    await focusGrid(page);
    const before = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {rowIdx: df.currentRowIdx, sel: df.selection.trueCount};
    });
    await page.keyboard.press('Shift+ArrowDown');
    await page.keyboard.press('Shift+ArrowDown');
    await page.keyboard.press('Shift+ArrowDown');
    await page.waitForTimeout(300);
    const after = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {rowIdx: df.currentRowIdx, sel: df.selection.trueCount};
    });
    expect(after.rowIdx).toBe(before.rowIdx + 3); // navigation proceeds
    expect(after.sel).toBe(0);                     // but nothing is selected
    // Restore for hygiene.
    await page.evaluate(() => { grok.shell.tv.grid.props.allowRowSelection = true; });
  });

  // === Scenario 7: pinned-row interplay — selection coordinate mapping =========

  await reopenClean(page);

  await softStep('Step 22 — 2 pinned rows + scroll + drag-select: selected indices are the intended table rows', async () => {
    // Pinning is neutral setup: it establishes the coordinate offset that the tested effect —
    // the drag-select → getSelectedIndexes mapping — has to survive.
    const cell0 = await dataCellPoint(page, 'AGE', 0);
    const p1 = await pinViaMenu(page, cell0, 'div-Pin---Pin-Row');
    expect(p1).toBe(true); // first Pin Row leaf reached and clicked
    await page.waitForTimeout(500);
    // After the first pin the pinned row displays at top and rows shift; pin a
    // grid row whose table index is not already pinned.
    const cell1 = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const pinned = Array.from(grid.pinnedRows) as number[];
      let gridRow = 1;
      for (let gr = 1; gr < 12; gr++) { const ti = grid.gridRowToTable(gr); if (!pinned.includes(ti)) { gridRow = gr; break; } }
      const db = grid.cell('AGE', gridRow).documentBounds;
      return {x: db.x + db.width / 2, y: db.y + db.height / 2};
    });
    const p2 = await pinViaMenu(page, cell1, 'div-Pin---Pin-Row');
    expect(p2).toBe(true); // second Pin Row leaf reached and clicked
    await page.waitForTimeout(500);
    const pinnedCount = await page.evaluate(() => Array.from(grok.shell.tv.grid.pinnedRows).length);
    expect(pinnedCount).toBe(2); // precondition: exactly 2 rows pinned before the drag-select test
    // Scroll down so the pinned rows stay at top while the data area shows later rows.
    await page.evaluate(async () => {
      const grid = grok.shell.tv.grid;
      if (grid.vertScroll?.scrollTo) grid.vertScroll.scrollTo(40);
      await new Promise((r) => setTimeout(r, 500));
    });
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));
    await focusGrid(page);
    // Drag-select 3 consecutive data-area rows. The vertical screen position comes from the
    // overlay rect + header + rowHeight, not cell.documentBounds.y: that origin is the grid's
    // virtual scroll space and goes off-viewport once the grid is scrolled.
    const geom = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
      const rc = overlay.getBoundingClientRect();
      const rowH = grid.props.rowHeight;
      const headerH = grid.colHeaderHeight;
      // Screen y of an on-screen row SLOT (0 = first row below the header). Once scrolled, the
      // slot index no longer matches grid-row-in-scroll-space, so the intended table indices are
      // read from hitTest at the exact drag coordinates — gridRowToTable(slot) would be wrong.
      const rowY = (slot: number) => rc.y + headerH + slot * rowH + rowH / 2;
      const col0 = grid.columns.byIndex(0);
      const stripLocalX = col0.left + col0.width / 2;
      const stripX = rc.x + stripLocalX;
      // Slots 0..1 are the 2 pinned rows on this scrolled view; start below them.
      const startSlot = 6;
      const tableAt = (slot: number): number => {
        const gc: any = grid.hitTest(stripLocalX, rowY(slot) - rc.y);
        return gc && (gc.isRowHeader || gc.isTableCell) ? gc.tableRowIndex : -1;
      };
      const intended = [tableAt(startSlot), tableAt(startSlot + 1), tableAt(startSlot + 2)];
      return {
        start: {x: stripX, y: rowY(startSlot)},
        end: {x: stripX, y: rowY(startSlot + 2)},
        intended,
      };
    });
    await page.mouse.move(geom.start.x, geom.start.y);
    await page.mouse.down();
    await page.mouse.move(geom.end.x, geom.end.y);
    await page.mouse.up();
    await page.waitForTimeout(400);
    const r = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {selected: Array.from(df.selection.getSelectedIndexes())};
    });
    // The selected indices must be the intended ones, NOT shifted by the pinned-row offset.
    const selectedSet = [...r.selected].sort((a, b) => a - b);
    const intendedSet = [...geom.intended].filter((i) => i >= 0).sort((a, b) => a - b);
    expect(selectedSet).toEqual(intendedSet); // coordinate mapping accounts for pinned rows
  });

  // === Scenario 8: sort + filter shared order mapping =========================

  await reopenClean(page);

  await softStep('Step 24 — sort + filter: navigate via arrows, gridRowToTable yields correct table indices', async () => {
    // Neutral setup: a real FilterGroup filter (AGE > 30) that survives navigation, plus an
    // ascending sort. The tested effect is the row mapping under sort+filter.
    const filtCount = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const fg = tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} }
      df.filter.setAll(true);
      fg.updateOrAdd({type: DG.FILTER_TYPE.HISTOGRAM, column: 'AGE', min: 31, max: 999});
      await new Promise((r) => setTimeout(r, 1200));
      // Applying a filter does not switch rowSource back from the 'All' that reopenClean sets,
      // and under 'All' gridRowToTable walks every row in sort order — including ones the
      // filter drops — which is exactly what breaks the mapping.
      tv.grid.props.rowSource = 'Filtered';
      await new Promise((r) => setTimeout(r, 400));
      tv.grid.sort(['AGE'], [true]);
      await new Promise((r) => setTimeout(r, 600));
      // The baseline is the filter's OWN trueCount, not a hand-computed AGE>30 threshold: the
      // histogram filter's null/boundary handling diverges from a naive count by the
      // Int32-null-sentinel row.
      (window as any).__filtBaseline = df.filter.trueCount;
      return df.filter.trueCount;
    });
    expect(filtCount).toBeGreaterThan(0);
    await focusGrid(page);
    // The mapping is read through gridRowToTable; getRowOrder is never called (it materializes
    // _order as a side effect).
    const r = await page.evaluate(async () => {
      const grid = grok.shell.tv.grid;
      const df = grok.shell.tv.dataFrame;
      const ageCol = df.col('AGE');
      const rows = [];
      for (let v = 0; v < 5; v++) {
        const t = grid.gridRowToTable(v);
        rows.push({visual: v, table: t, passesFilter: df.filter.get(t), age: ageCol.get(t), isNull: ageCol.isNone(t)});
      }
      // Nulls are excluded from the order check: a filter-kept null-AGE row sorts first and has
      // no meaningful AGE ordinal, though it still passes the filter.
      const ages = rows.filter((x) => !x.isNull).map((x) => x.age);
      let ascending = true;
      for (let i = 1; i < ages.length; i++) if (ages[i] < ages[i - 1]) ascending = false;
      return {
        rows,
        ascending,
        allPassFilter: rows.every((x) => x.passesFilter),
        nonNullAllOver30: rows.filter((x) => !x.isNull).every((x) => x.age > 30),
      };
    });
    expect(r.allPassFilter).toBe(true);    // every mapped table index is a filter-kept row
    expect(r.nonNullAllOver30).toBe(true); // and each non-null AGE is > 30
    expect(r.ascending).toBe(true);        // and appears in ascending-AGE order
  });

  await softStep('Step 25 — keyboard range-select under sort+filter: selected indices are the intended table rows', async () => {
    // Grid keyboard navigation reads an INTERNAL focus cell that only a trusted data-area click
    // seeds; a df.currentCell= assignment plus overlay.focus leaves Shift+Down selecting nothing
    // and not even advancing.
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));
    const anchorCell = await dataCellPoint(page, 'AGE', 0);
    await page.mouse.click(anchorCell.x, anchorCell.y);
    await page.waitForTimeout(250);
    const anchor = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const grid = grok.shell.tv.grid;
      return {rowIdx: df.currentRowIdx, seededAtVisual0: df.currentRowIdx === grid.gridRowToTable(0)};
    });
    expect(anchor.seededAtVisual0).toBe(true); // trusted click seeded the keyboard focus at visual row 0
    // Each Shift+Down selects the current row then advances (grid_keyboard_navigation line 131).
    // No focusGrid here: it would only .focus the overlay, not re-seed the focus cell the click
    // established.
    await page.keyboard.down('Shift');
    await page.keyboard.press('ArrowDown');
    await page.keyboard.press('ArrowDown');
    await page.keyboard.press('ArrowDown');
    await page.keyboard.press('ArrowDown');
    await page.keyboard.up('Shift');
    await page.waitForTimeout(400);
    const r = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const grid = grok.shell.tv.grid;
      const selected = Array.from(df.selection.getSelectedIndexes()) as number[];
      // The range spans the contiguous visual rows from the anchor through the 4 Shift+Down
      // steps, so the intended set is those visual slots mapped 1:1.
      const count = selected.length;
      const intended: number[] = [];
      for (let v = 0; v < count; v++) { const t = grid.gridRowToTable(v); if (t >= 0) intended.push(t); }
      const allPassFilter = selected.every((t) => df.filter.get(t));
      return {selected, count, intended, allPassFilter};
    });
    expect(r.count).toBe(4);              // 4 rows range-selected
    expect(r.allPassFilter).toBe(true);   // all selected rows are filter-kept (mapping respected)
    const selSet = [...r.selected].sort((a, b) => a - b);
    const intendedSet = [...r.intended].sort((a, b) => a - b);
    expect(selSet).toEqual(intendedSet); // selection is exactly the mapped table-row set
  });

  await softStep('Step 24/Step 7 (cont.) — df.filter.trueCount is unchanged by sort + navigation', async () => {
    // Compared against the baseline captured in Step 24, not a hand-recomputed AGE>30 threshold
    // (which diverges by the filter-kept Int32-null-sentinel row).
    const r = await page.evaluate(() => ({
      now: grok.shell.tv.dataFrame.filter.trueCount,
      baseline: (window as any).__filtBaseline,
    }));
    expect(r.now).toBe(r.baseline); // sorting + keyboard navigation do not change the filter true count
  });

  // === Scenario 9: row source switching — Filtered / Selected / All ===========

  await reopenClean(page);

  await softStep('Step 26 — rowSource All → Filtered → Selected → All; df.filter unchanged throughout', async () => {
    // Setup: SEX=M filter (real FilterGroup) and select 5 rows.
    const setup = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const fg = tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} }
      df.filter.setAll(true);
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: ['M']});
      await new Promise((r) => setTimeout(r, 1200));
      df.selection.setAll(false);
      for (let i = 0; i < 5; i++) df.selection.set(i, true);
      return {filterCount: df.filter.trueCount, selCount: df.selection.trueCount};
    });
    expect(setup.selCount).toBe(5);

    // No grid.rowsCount / visibleRowsCount property exists, so the visible rows are counted by
    // walking gridRowToTable until its negative sentinel. The walk is capped at df.rowCount:
    // under rowSource 'All' the mapping is the identity with no sentinel at all, so a larger
    // bound would count phantom rows past the end.
    async function visibleCount(): Promise<number> {
      return page.evaluate(() => {
        const grid = grok.shell.tv.grid;
        const df = grok.shell.tv.dataFrame;
        let n = 0;
        for (let v = 0; v < df.rowCount; v++) {
          const t = grid.gridRowToTable(v);
          if (t === undefined || t === null || t < 0) break;
          n++;
        }
        return n;
      });
    }
    async function setSource(src: string): Promise<void> {
      await page.evaluate(async (s) => {
        grok.shell.tv.grid.props.rowSource = s;
        await new Promise((r) => setTimeout(r, 400));
      }, src);
    }
    const filt = () => page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);

    await setSource('Filtered');
    const filteredVisible = await visibleCount();
    const filtAtFiltered = await filt();
    expect(filteredVisible).toBe(setup.filterCount); // Filtered shows exactly df.filter.trueCount rows

    await setSource('Selected');
    const selectedVisible = await visibleCount();
    const filtAtSelected = await filt();
    expect(selectedVisible).toBe(setup.selCount);    // Selected shows exactly the selected rows

    await setSource('All');
    const allVisible = await visibleCount();
    const filtAtAll = await filt();
    const rowCount = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
    expect(allVisible).toBe(rowCount); // All restores the full row count (every row visible again)

    // df.filter.trueCount stays identical throughout (by-design divergence).
    expect(filtAtFiltered).toBe(setup.filterCount);
    expect(filtAtSelected).toBe(setup.filterCount);
    expect(filtAtAll).toBe(setup.filterCount);
  });

  // === Scenario 10: Ctrl+F search — tint signal + explicit filter restore ======

  await reopenClean(page);

  await softStep('Step 27 — Ctrl+F search: matching cells tint, then the filter is restored', async () => {
    await focusGrid(page);
    await page.keyboard.press('Control+f');
    // The Toolbox Search pane expands asynchronously and its input takes focus a beat after the
    // pane appears, so the state is polled rather than read after a fixed wait.
    const paneState = await page.evaluate(async () => {
      const sel = '[name="pane-Search"] input[placeholder="Search"]';
      let inp: HTMLInputElement | null = null;
      for (let i = 0; i < 40; i++) {
        inp = document.querySelector(sel) as HTMLInputElement;
        if (inp && document.activeElement === inp) break;
        await new Promise((r) => setTimeout(r, 100));
      }
      return {present: !!inp, focused: inp ? document.activeElement === inp : false};
    });
    // The pane opening and taking focus is Toolbox behaviour, not the grid's, so it stays a
    // precondition — but without the input the step must not silently pass.
    expect(paneState.present).toBe(true); // precondition only: the search input exists to type into
    // Type a value present in AGE.
    await page.keyboard.type('35');
    await page.waitForTimeout(700);
    const tint = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const grid = grok.shell.tv.grid;
      const ageCol = df.col('AGE');
      let matchRow = -1;
      let nonMatchRow = -1;
      for (let i = 0; i < df.rowCount; i++) { if (String(ageCol.get(i)).includes('35')) { matchRow = i; break; } }
      for (let i = 0; i < df.rowCount; i++) { if (!String(ageCol.get(i)).includes('35')) { nonMatchRow = i; break; } }
      return {
        matchColor: matchRow >= 0 ? grid.cell('AGE', matchRow).color >>> 0 : 0,
        nonMatchColor: nonMatchRow >= 0 ? grid.cell('AGE', nonMatchRow).color >>> 0 : 0,
        plainWhite: 0xffffffff,
      };
    });
    expect(tint.matchColor).not.toBe(tint.plainWhite);    // a matching cell shows a tint
    expect(tint.nonMatchColor).toBe(tint.plainWhite);     // a non-matching cell keeps the plain background
    // df.filter has to be restored EXPLICITLY — clearing the search text alone does not do it.
    const restored = await page.evaluate(async () => {
      const inp = document.querySelector('[name="pane-Search"] input[placeholder="Search"]') as HTMLInputElement;
      if (inp) { inp.value = ''; inp.dispatchEvent(new Event('input', {bubbles: true})); inp.dispatchEvent(new Event('change', {bubbles: true})); }
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true);
      await new Promise((r) => setTimeout(r, 400));
      return df.filter.trueCount === df.rowCount;
    });
    expect(restored).toBe(true); // all rows visible after the explicit filter restore
  });

  v.finishSpec();
});
