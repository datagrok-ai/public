/* ---
realizes: [grid.cp.rows-select-filter-navigate]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

async function headerLabelPoint(page: Page, col: string): Promise<{x: number; y: number}> {
  return page.evaluate((c) => {
    const grid = grok.shell.tv.grid;
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    const rc = overlay.getBoundingClientRect();
    const gc = grid.columns.byName(c);
    return {x: rc.x + gc.left + gc.width / 2, y: rc.y + grid.colHeaderHeight / 2};
  }, col);
}

async function dataCellPoint(page: Page, col: string, visualRow: number): Promise<{x: number; y: number}> {
  return page.evaluate(({c, vr}) => {
    const grid = grok.shell.tv.grid;
    const db = grid.cell(c, vr).documentBounds;
    return {x: db.x + db.width / 2, y: db.y + db.height / 2};
  }, {c: col, vr: visualRow});
}

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

async function stripDragSelect(page: Page, startVisualRow: number, endVisualRow: number): Promise<void> {
  const s = await rowStripPoint(page, startVisualRow);
  const e = await rowStripPoint(page, endVisualRow);
  await page.mouse.move(s.x, s.y);
  await page.mouse.down();
  await page.mouse.move(e.x, e.y);
  await page.mouse.up();
  await waitDfEvent(page, 'onSelectionChanged', 400); 
}

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

async function firstDataCol(page: Page): Promise<string> {
  return page.evaluate(() => grok.shell.tv.grid.columns.byIndex(1).name);
}

async function focusGrid(page: Page): Promise<void> {
  await page.evaluate(() => {
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    overlay.focus();
  });
}

const DF_STREAMS = ['onSelectionChanged', 'onCurrentRowChanged', 'onCurrentCellChanged',
  'onColumnSelectionChanged', 'onRowsFiltered'] as const;
type DfStream = typeof DF_STREAMS[number];

async function installDfStamps(page: Page): Promise<void> {
  await page.evaluate((streams) => {
    const w = window as any;
    const df = grok.shell.tv.dataFrame;
    if (w.__dfStampFor === df) return;
    w.__dfStampFor = df;
    w.__dfStamp = {};
    for (const s of streams) {
      try { df[s].subscribe(() => { w.__dfStamp[s] = Date.now(); }); } catch (_) {  }
    }
  }, DF_STREAMS as unknown as string[]);
}

async function waitDfEvent(page: Page, stream: DfStream, capMs: number): Promise<void> {
  await page.evaluate(async ({s, cap}) => {
    const w = window as any;
    const stamp = () => w.__dfStamp?.[s] ?? 0;
    const t0 = Date.now();
    if (stamp() && t0 - stamp() < 150) return; 
    const before = stamp();
    for (;;) {
      if (stamp() > before || Date.now() - t0 >= cap) return;
      await new Promise((r) => setTimeout(r, 25));
    }
  }, {s: stream, cap: capMs});
}

async function reopenClean(page: Page): Promise<void> {
  await v.closeAllAndWait(page);
  await v.openTable(page, {path: 'System:DemoFiles/demog.csv', semTypeTimeoutMs: 3000});
  await page.evaluate(() => {
    const tv = grok.shell.tv;
    const df = tv.dataFrame;
    const grid = tv.grid;

    df.filter.setAll(true);
    df.selection.setAll(false);
    grid.sort([], []);
    grid.props.rowSource = 'All';
  });
  await installDfStamps(page);
}

test('Grid — Row Selection, Filter, and Keyboard Navigation', async ({page}) => {
  test.setTimeout(420_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: 'System:DemoFiles/demog.csv', semTypeTimeoutMs: 3000});

  const consoleErrors: string[] = [];
  const pageErrors: string[] = [];
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });
  page.on('pageerror', (e) => pageErrors.push(String(e)));

  await page.evaluate(() => {
    const grid = grok.shell.tv.grid;
    const df = grok.shell.tv.dataFrame;

    df.filter.setAll(true);
    df.selection.setAll(false);
    grid.sort([], []);
  });
  await installDfStamps(page);

  await softStep('Step 3 — click row 5 data cell: sets current row, no selection', async () => {
    await focusGrid(page);
    const col = await firstDataCol(page);
    const p = await dataCellPoint(page, col, 5);
    await page.mouse.click(p.x, p.y);
    await waitDfEvent(page, 'onCurrentCellChanged', 300); 
    const r = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const grid = grok.shell.tv.grid;
      return {
        currentRowIdx: df.currentRowIdx,
        expectedTableIdx: grid.gridRowToTable(5),
        selTrueCount: df.selection.trueCount,
      };
    });
    expect(r.currentRowIdx).toBe(r.expectedTableIdx); 
    expect(r.selTrueCount).toBe(0);                   
  });

  await softStep('Step 4 — drag-select data rows 5..10: range select 6 rows', async () => {

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
    expect(selSet).toEqual(intendedSet);              
    expect(r.selTrueCount).toBe(intendedSet.length);  
  });

  await softStep('Step 6 — Ctrl+A: select all rows', async () => {
    await focusGrid(page);
    await page.keyboard.press('Control+a');
    await waitDfEvent(page, 'onSelectionChanged', 300); 
    const r = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {selTrueCount: df.selection.trueCount, rowCount: df.rowCount};
    });
    expect(r.selTrueCount).toBe(r.rowCount); 
  });

  await softStep('Step 7 — Esc: compound reset (selection, columns, current row all cleared)', async () => {

    await page.evaluate(() => { grok.shell.tv.grid.columns.byName('AGE').selected = true; grok.shell.tv.grid.invalidate(); });
    await waitDfEvent(page, 'onColumnSelectionChanged', 150);
    await focusGrid(page);
    await page.keyboard.press('Escape');
    await waitDfEvent(page, 'onSelectionChanged', 300); 
    const r = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {
        selTrueCount: df.selection.trueCount,
        selectedColNames: Array.from(df.columns.selected).map((c: any) => c.name),
        currentRowIdx: df.currentRowIdx,
      };
    });
    expect(r.selTrueCount).toBe(0);             
    expect(r.selectedColNames).toEqual([]);     
    expect(r.currentRowIdx).toBe(-1);           
  });

  await reopenClean(page);

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
    expect(concurrentModErrors).toEqual([]);                        
    expect(consoleErrors.length + pageErrors.length).toBe(errBefore); 
  });

  await reopenClean(page);

  await softStep('Step 11 — Shift+drag over data cells: selects rows only, no column selection', async () => {
    await focusGrid(page);
    await page.keyboard.press('Escape');
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false)); 

    const start = await dataCellPoint(page, 'AGE', 2);
    const end = await dataCellPoint(page, 'AGE', 6);
    const crossed = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      return [2, 3, 4, 5, 6].map((v) => grid.gridRowToTable(v)).filter((t) => t >= 0);
    });

    await page.keyboard.down('Shift');
    await page.mouse.move(start.x, start.y);
    await page.mouse.down();
    await page.mouse.move(end.x, end.y);
    await page.mouse.up();
    await page.keyboard.up('Shift');
    await waitDfEvent(page, 'onSelectionChanged', 400); 
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
    expect(selectedSet).toEqual(crossedSet); 
    expect(r.selTrueCount).toBe(crossedSet.length); 
    expect(r.selectedColNames).toEqual([]);      
  });

  await reopenClean(page);

  await softStep('Step 13 — arrows: Down x2 then Right x2 move current row/col', async () => {

    const p = await dataCellPoint(page, 'AGE', 4);
    await page.mouse.click(p.x, p.y);
    await waitDfEvent(page, 'onCurrentCellChanged', 200); 
    const before = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {rowIdx: df.currentRowIdx, col: df.currentCol ? df.currentCol.name : null};
    });
    await focusGrid(page);
    await page.keyboard.press('ArrowDown');
    await page.keyboard.press('ArrowDown');
    await page.keyboard.press('ArrowRight');
    await page.keyboard.press('ArrowRight');
    await waitDfEvent(page, 'onCurrentCellChanged', 300); 
    const after = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {rowIdx: df.currentRowIdx, col: df.currentCol ? df.currentCol.name : null};
    });
    expect(after.rowIdx).toBe(before.rowIdx + 2);   
    expect(after.col).not.toBe(before.col);          
  });

  await softStep('Step 15 — Ctrl+Home goes to first row; Ctrl+End goes to last row', async () => {
    await focusGrid(page);
    await page.keyboard.press('Control+Home');
    await waitDfEvent(page, 'onCurrentRowChanged', 300); 
    const home = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
    expect(home).toBe(0); 
    await focusGrid(page);
    await page.keyboard.press('Control+End');
    await waitDfEvent(page, 'onCurrentRowChanged', 300); 
    const r = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {rowIdx: df.currentRowIdx, rowCount: df.rowCount};
    });
    expect(r.rowIdx).toBe(r.rowCount - 1); 
  });

  await softStep('Step 16 — PageDown twice advances both times (GROK-18104)', async () => {

    await focusGrid(page);
    await page.keyboard.press('Control+Home');
    await waitDfEvent(page, 'onCurrentRowChanged', 250); 
    const p0 = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
    await page.keyboard.press('PageDown');
    await waitDfEvent(page, 'onCurrentRowChanged', 300); 
    const p1 = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
    await page.keyboard.press('PageDown');
    await waitDfEvent(page, 'onCurrentRowChanged', 300); 
    const p2 = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
    expect(p1).toBeGreaterThan(p0); 
    expect(p2).toBeGreaterThan(p1); 
  });

  await reopenClean(page);

  await softStep('Step 17 (also proves the Step 5 toggle-membership channel) — Space toggles current row selection membership', async () => {

    const p = await dataCellPoint(page, 'AGE', 3);
    await page.mouse.click(p.x, p.y);
    await waitDfEvent(page, 'onCurrentCellChanged', 200); 
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false)); 
    await focusGrid(page);
    const before = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    await page.keyboard.press(' ');
    await waitDfEvent(page, 'onSelectionChanged', 250); 
    const afterAdd = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    await page.keyboard.press(' ');
    await waitDfEvent(page, 'onSelectionChanged', 250); 
    const afterRemove = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(afterAdd).toBe(before + 1);      
    expect(afterRemove).toBe(before);       
  });

  await softStep('Step 18 — Shift+Enter selects all rows sharing the current cell value', async () => {

    const expected = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      df.currentCell = df.cell(3, 'SEX');
      const sexVal = df.col('SEX').get(df.currentRowIdx);
      let count = 0;
      for (let i = 0; i < df.rowCount; i++) if (df.col('SEX').get(i) === sexVal) count++;
      return {count, sexVal};
    });
    await focusGrid(page);
    await page.keyboard.press('Shift+Enter');
    await waitDfEvent(page, 'onSelectionChanged', 400); 
    const sel = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(sel).toBe(expected.count); 
  });

  await reopenClean(page);

  await softStep('Step 20 — allowRowSelection false: Shift+Down navigates without selecting', async () => {

    await page.evaluate(() => { grok.shell.tv.grid.props.allowRowSelection = false; }); 
    const p = await dataCellPoint(page, 'AGE', 4);
    await page.mouse.click(p.x, p.y);
    await waitDfEvent(page, 'onCurrentCellChanged', 200); 
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false)); 
    await focusGrid(page);
    const before = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {rowIdx: df.currentRowIdx, sel: df.selection.trueCount};
    });
    await page.keyboard.press('Shift+ArrowDown');
    await page.keyboard.press('Shift+ArrowDown');
    await page.keyboard.press('Shift+ArrowDown');
    await waitDfEvent(page, 'onCurrentRowChanged', 300); 
    const after = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {rowIdx: df.currentRowIdx, sel: df.selection.trueCount};
    });
    expect(after.rowIdx).toBe(before.rowIdx + 3); 
    expect(after.sel).toBe(0);                     

    await page.evaluate(() => { grok.shell.tv.grid.props.allowRowSelection = true; });
  });

  await reopenClean(page);

  await softStep('Step 22 — 2 pinned rows + scroll + drag-select: selected indices are the intended table rows', async () => {

    const cell0 = await dataCellPoint(page, 'AGE', 0);
    const p1 = await pinViaMenu(page, cell0, 'div-Pin---Pin-Row');
    expect(p1).toBe(true); 

    await page.waitForFunction(() => Array.from(grok.shell.tv.grid.pinnedRows).length >= 1,
      null, {timeout: 3000}).catch(() => {});

    const cell1 = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const pinned = Array.from(grid.pinnedRows) as number[];
      let gridRow = 1;
      for (let gr = 1; gr < 12; gr++) { const ti = grid.gridRowToTable(gr); if (!pinned.includes(ti)) { gridRow = gr; break; } }
      const db = grid.cell('AGE', gridRow).documentBounds;
      return {x: db.x + db.width / 2, y: db.y + db.height / 2};
    });
    const p2 = await pinViaMenu(page, cell1, 'div-Pin---Pin-Row');
    expect(p2).toBe(true); 

    await page.waitForFunction(() => Array.from(grok.shell.tv.grid.pinnedRows).length >= 2,
      null, {timeout: 3000}).catch(() => {});
    const pinnedCount = await page.evaluate(() => Array.from(grok.shell.tv.grid.pinnedRows).length);
    expect(pinnedCount).toBe(2); 

    await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      if (grid.vertScroll?.scrollTo) grid.vertScroll.scrollTo(40);
    });
    await page.waitForTimeout(500); 
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));
    await focusGrid(page);

    const geom = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
      const rc = overlay.getBoundingClientRect();
      const rowH = grid.props.rowHeight;
      const headerH = grid.colHeaderHeight;

      const rowY = (slot: number) => rc.y + headerH + slot * rowH + rowH / 2;
      const col0 = grid.columns.byIndex(0);
      const stripLocalX = col0.left + col0.width / 2;
      const stripX = rc.x + stripLocalX;

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
    await waitDfEvent(page, 'onSelectionChanged', 400); 
    const r = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return {selected: Array.from(df.selection.getSelectedIndexes())};
    });

    const selectedSet = [...r.selected].sort((a, b) => a - b);
    const intendedSet = [...geom.intended].filter((i) => i >= 0).sort((a, b) => a - b);
    expect(selectedSet).toEqual(intendedSet); 
  });

  await reopenClean(page);

  await softStep('Step 24 — sort + filter: navigate via arrows, gridRowToTable yields correct table indices', async () => {

    const filtCount = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const fg = tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} }
      df.filter.setAll(true);

      const filtered = new Promise<void>((resolve) => {
        let sub: any = null;
        try { sub = df.onRowsFiltered.subscribe(() => { sub.unsubscribe(); resolve(); }); }
        catch (_) {  }
        setTimeout(() => { try { sub?.unsubscribe(); } catch (_) {} resolve(); }, 1200);
      });
      fg.updateOrAdd({type: DG.FILTER_TYPE.HISTOGRAM, column: 'AGE', min: 31, max: 999});
      await filtered;

      tv.grid.props.rowSource = 'Filtered';
      await new Promise((r) => setTimeout(r, 400)); 
      tv.grid.sort(['AGE'], [true]);
      await new Promise((r) => setTimeout(r, 600)); 

      (window as any).__filtBaseline = df.filter.trueCount;
      return df.filter.trueCount;
    });
    expect(filtCount).toBeGreaterThan(0);
    await focusGrid(page);

    const r = await page.evaluate(async () => {
      const grid = grok.shell.tv.grid;
      const df = grok.shell.tv.dataFrame;
      const ageCol = df.col('AGE');
      const rows = [];
      for (let v = 0; v < 5; v++) {
        const t = grid.gridRowToTable(v);
        rows.push({visual: v, table: t, passesFilter: df.filter.get(t), age: ageCol.get(t), isNull: ageCol.isNone(t)});
      }

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
    expect(r.allPassFilter).toBe(true);    
    expect(r.nonNullAllOver30).toBe(true); 
    expect(r.ascending).toBe(true);        
  });

  await softStep('Step 25 — keyboard range-select under sort+filter: selected indices are the intended table rows', async () => {

    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false)); 
    const anchorCell = await dataCellPoint(page, 'AGE', 0);
    await page.mouse.click(anchorCell.x, anchorCell.y);
    await waitDfEvent(page, 'onCurrentCellChanged', 250); 
    const anchor = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const grid = grok.shell.tv.grid;
      return {rowIdx: df.currentRowIdx, seededAtVisual0: df.currentRowIdx === grid.gridRowToTable(0)};
    });
    expect(anchor.seededAtVisual0).toBe(true); 

    await page.keyboard.down('Shift');
    await page.keyboard.press('ArrowDown');
    await page.keyboard.press('ArrowDown');
    await page.keyboard.press('ArrowDown');
    await page.keyboard.press('ArrowDown');
    await page.keyboard.up('Shift');
    await waitDfEvent(page, 'onSelectionChanged', 400); 
    const r = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const grid = grok.shell.tv.grid;
      const selected = Array.from(df.selection.getSelectedIndexes()) as number[];

      const count = selected.length;
      const intended: number[] = [];
      for (let v = 0; v < count; v++) { const t = grid.gridRowToTable(v); if (t >= 0) intended.push(t); }
      const allPassFilter = selected.every((t) => df.filter.get(t));
      return {selected, count, intended, allPassFilter};
    });
    expect(r.count).toBe(4);              
    expect(r.allPassFilter).toBe(true);   
    const selSet = [...r.selected].sort((a, b) => a - b);
    const intendedSet = [...r.intended].sort((a, b) => a - b);
    expect(selSet).toEqual(intendedSet); 
  });

  await softStep('Step 24/Step 7 (cont.) — df.filter.trueCount is unchanged by sort + navigation', async () => {

    const r = await page.evaluate(() => ({
      now: grok.shell.tv.dataFrame.filter.trueCount,
      baseline: (window as any).__filtBaseline,
    }));
    expect(r.now).toBe(r.baseline); 
  });

  await reopenClean(page);

  await softStep('Step 26 — rowSource All → Filtered → Selected → All; df.filter unchanged throughout', async () => {

    const setup = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const fg = tv.getFiltersGroup();
      for (const f of Array.from(fg.filters as any)) { try { fg.remove(f); } catch (_) {} }
      df.filter.setAll(true);

      const filtered = new Promise<void>((resolve) => {
        let sub: any = null;
        try { sub = df.onRowsFiltered.subscribe(() => { sub.unsubscribe(); resolve(); }); }
        catch (_) {  }
        setTimeout(() => { try { sub?.unsubscribe(); } catch (_) {} resolve(); }, 1200);
      });
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: ['M']});
      await filtered;
      df.selection.setAll(false);
      for (let i = 0; i < 5; i++) df.selection.set(i, true);
      return {filterCount: df.filter.trueCount, selCount: df.selection.trueCount};
    });
    expect(setup.selCount).toBe(5);

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
      await page.evaluate((s) => { grok.shell.tv.grid.props.rowSource = s; }, src);
      await page.waitForTimeout(400); 
    }
    const filt = () => page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);

    await setSource('Filtered');
    const filteredVisible = await visibleCount();
    const filtAtFiltered = await filt();
    expect(filteredVisible).toBe(setup.filterCount); 

    await setSource('Selected');
    const selectedVisible = await visibleCount();
    const filtAtSelected = await filt();
    expect(selectedVisible).toBe(setup.selCount);    

    await setSource('All');
    const allVisible = await visibleCount();
    const filtAtAll = await filt();
    const rowCount = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
    expect(allVisible).toBe(rowCount); 

    expect(filtAtFiltered).toBe(setup.filterCount);
    expect(filtAtSelected).toBe(setup.filterCount);
    expect(filtAtAll).toBe(setup.filterCount);
  });

  await reopenClean(page);

  await softStep('Step 27 — Ctrl+F search: matching cells tint, then the filter is restored', async () => {
    await focusGrid(page);
    await page.keyboard.press('Control+f');

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

    expect(paneState.present).toBe(true); 

    await page.keyboard.type('35');

    await page.waitForFunction(() => {
      const df = grok.shell.tv.dataFrame;
      const grid = grok.shell.tv.grid;
      const ageCol = df.col('AGE');
      for (let i = 0; i < df.rowCount; i++)
        if (String(ageCol.get(i)).includes('35')) return (grid.cell('AGE', i).color >>> 0) !== 0xffffffff;
      return false;
    }, null, {timeout: 3000, polling: 100}).catch(() => {});
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
    expect(tint.matchColor).not.toBe(tint.plainWhite);    
    expect(tint.nonMatchColor).toBe(tint.plainWhite);     

    const restored = await page.evaluate(async () => {
      const inp = document.querySelector('[name="pane-Search"] input[placeholder="Search"]') as HTMLInputElement;
      if (inp) { inp.value = ''; inp.dispatchEvent(new Event('input', {bubbles: true})); inp.dispatchEvent(new Event('change', {bubbles: true})); }
      const df = grok.shell.tv.dataFrame;
      df.filter.setAll(true);
      await new Promise((r) => setTimeout(r, 400)); 
      return df.filter.trueCount === df.rowCount;
    });
    expect(restored).toBe(true); 
  });

  v.finishSpec();
});
