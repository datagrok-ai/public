/* ---
realizes: [grid.cp.edit-clipboard]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';


declare const grok: any;
declare const DG: any;


test.use(specTestOptions);


// Page-coordinate center of a data cell (col, visualRow). documentBounds is already
// in page coordinates (grid refdoc, "Coordinates for canvas gestures").
async function dataCellPoint(page: Page, col: string, visualRow: number): Promise<{x: number; y: number}> {
  return page.evaluate(({c, vr}) => {
    const grid = grok.shell.tv.grid;
    const db = grid.cell(c, vr).documentBounds;
    return {x: db.x + db.width / 2, y: db.y + db.height / 2};
  }, {c: col, vr: visualRow});
}

// Focus the grid overlay so keyboard gestures target it.
async function focusGrid(page: Page): Promise<void> {
  await page.evaluate(() => {
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    overlay.focus();
  });
}

// Trusted single click on a data cell, settled on the df channel reporting that cell as current.
async function clickCellSettled(page: Page, col: string, visualRow: number): Promise<{x: number; y: number}> {
  const p = await dataCellPoint(page, col, visualRow);
  const tableRow = await page.evaluate((vr) => grok.shell.tv.grid.gridRowToTable(vr), visualRow);
  await focusGrid(page);
  await page.mouse.click(p.x, p.y);
  await page.evaluate(async ({c, tr}) => {
    const df = grok.shell.tv.dataFrame;
    for (let i = 0; i < 40; i++) {
      if (df.currentRowIdx === tr && df.currentCol != null && df.currentCol.name === c) break;
      await new Promise((r) => setTimeout(r, 50));
    }
  }, {c: col, tr: tableRow});
  return p;
}

// Double-click a data cell to open (or, when allowEdit=false, to attempt) the inline editor. On a
// cold render the double-click can outrace the current-cell assignment the open gates on, so the
// retry re-seeds with a trusted click and opens with Enter; the read-only path passes attempts=0.
async function dblClickCell(page: Page, col: string, visualRow: number, attempts = 6): Promise<void> {
  const p = await clickCellSettled(page, col, visualRow);
  await page.waitForTimeout(150);
  await page.mouse.dblclick(p.x, p.y);
  await page.waitForTimeout(200);
  for (let i = 0; i < attempts && (await editorCount(page)) === 0; i++) {
    await clickCellSettled(page, col, visualRow);
    await page.waitForTimeout(120);
    await page.keyboard.press('Enter');
    await page.waitForTimeout(200);
  }
}

// Count of the visible DOM cell editors (0 = no editor open, 1 = editing).
async function editorCount(page: Page): Promise<number> {
  return page.evaluate(() => document.querySelectorAll('input.d4-value-editor').length);
}

// Install the clipboard-copy interceptor. The single-cell copy runs through
// document.execCommand('copy'); the copy event's clipboardData carries no payload and
// navigator.clipboard.readText is denied here, so patching execCommand is the only way to read it.
async function installCopyInterceptor(page: Page): Promise<void> {
  await page.evaluate(() => {
    const w = window as any;
    w.__copyCaptured = [];
    if (!w.__copyPatched) {
      const origExec = document.execCommand.bind(document);
      (document as any).execCommand = (cmd: string, ...rest: any[]) => {
        if (cmd === 'copy') {
          const sel = w.getSelection?.().toString();
          const ae = document.activeElement as any;
          const aeVal = ae && (ae.tagName === 'TEXTAREA' || ae.tagName === 'INPUT') ? ae.value : null;
          w.__copyCaptured.push(sel && sel.length ? sel : aeVal);
        }
        return origExec(cmd, ...rest);
      };
      w.__copyPatched = true;
    }
  });
}

async function resetCopyCapture(page: Page): Promise<void> {
  await page.evaluate(() => { (window as any).__copyCaptured = []; });
}

async function readCopyCapture(page: Page): Promise<string[]> {
  return page.evaluate(() => ((window as any).__copyCaptured ?? []).filter((x: any) => x != null));
}

// Reopen a clean demog view (close all, reopen, clear filter/selection/sort) and
// re-attach the onCellValueEdited counter + copy interceptor.
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
    const w = window as any;
    w.__editCount = 0;
    if (w.__editSub) w.__editSub.unsubscribe();
    w.__editSub = grid.onCellValueEdited.subscribe(() => w.__editCount++);
    await new Promise((r) => setTimeout(r, 300));
  });
  await installCopyInterceptor(page);
}

async function editCount(page: Page): Promise<number> {
  return page.evaluate(() => (window as any).__editCount ?? 0);
}

test('Grid — Cell Editing and Clipboard', async ({page}) => {
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
    const w = window as any;
    w.__editCount = 0;
    if (w.__editSub) w.__editSub.unsubscribe();
    w.__editSub = grid.onCellValueEdited.subscribe(() => w.__editCount++);
    await new Promise((r) => setTimeout(r, 200));
  });
  await installCopyInterceptor(page);

  // === Scenario 1: basic cell edit — commit with Enter, cancel with Esc =========

  await softStep('Step 4 — double-click AGE cell, type a new value, Enter commits it', async () => {
    const tableIdx = await page.evaluate(() => grok.shell.tv.grid.gridRowToTable(3));
    await dblClickCell(page, 'AGE', 3);
    await page.waitForTimeout(300);
    expect(await editorCount(page)).toBe(1); // Step 3: editor opened
    const editBefore = await editCount(page);
    await page.keyboard.press('Control+a');
    await page.keyboard.type('99');
    await page.keyboard.press('Enter');
    await page.waitForTimeout(400);
    const r = await page.evaluate((ti) => {
      const df = grok.shell.tv.dataFrame;
      return {age: df.col('AGE').get(ti), editorOpen: document.querySelectorAll('input.d4-value-editor').length > 0};
    }, tableIdx);
    const editAfter = await editCount(page);
    expect(r.age).toBe(99);            // committed value equals the typed value
    expect(editAfter).toBe(editBefore + 1); // onCellValueEdited fired once
    expect(r.editorOpen).toBe(false);  // editor dismissed after Enter
  });

  await softStep('Step 6 — re-open the editor, Esc leaves the value unchanged', async () => {
    const tableIdx = await page.evaluate(() => grok.shell.tv.grid.gridRowToTable(3));
    const before = await page.evaluate((ti) => grok.shell.tv.dataFrame.col('AGE').get(ti), tableIdx);
    await dblClickCell(page, 'AGE', 3);
    await page.waitForTimeout(300);
    expect(await editorCount(page)).toBe(1); // editor open
    await page.keyboard.press('Escape');
    await page.waitForTimeout(300);
    const r = await page.evaluate((ti) => ({
      age: grok.shell.tv.dataFrame.col('AGE').get(ti),
      editorOpen: document.querySelectorAll('input.d4-value-editor').length > 0,
    }), tableIdx);
    expect(r.age).toBe(before);       // Esc left the value intact
    expect(r.editorOpen).toBe(false); // editor dismissed by Esc
  });

  // === Scenario 2: Delete key clears the current cell (no modifier) =============

  await reopenClean(page);

  await softStep('Step 8 — Delete (no modifier) clears the cell, opens no editor', async () => {
    // The plain-Delete handler acts on grid.currentGridCell, so seed it with a settled click.
    // The keystroke reaches the grid only while the overlay holds focus, and a re-seat can land on
    // a different row — so the row under test is read inside the retry loop, not before it.
    await clickCellSettled(page, 'AGE', 4);
    let cleared = false;
    for (let attempt = 0; attempt < 3 && !cleared; attempt++) {
      if (attempt > 0) await clickCellSettled(page, 'AGE', 4);
      const row = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
      await page.keyboard.press('Delete');
      for (let waited = 0; waited < 3000 && !cleared; waited += 250) {
        await page.waitForTimeout(250);
        cleared = await page.evaluate((ti) =>
          grok.shell.tv.dataFrame.col('AGE').isNone(ti), row);
      }
    }
    expect(cleared).toBe(true); // the cell was cleared to null
    const editorOpen = await page.evaluate(() =>
      document.querySelectorAll('input.d4-value-editor').length > 0);
    expect(editorOpen).toBe(false); // Delete did not open an editor
  });

  // === Scenario 3: read-only grid rejects edits (GROK-20010) ====================

  await reopenClean(page);

  await softStep('Step 11 — allowEdit=false: double-click opens no editor + a read-only balloon', async () => {
    // Neutral setup of the read-only entry state; the assertable effect is the rejection.
    await page.evaluate(() => { grok.shell.tv.grid.props.allowEdit = false; });
    await page.waitForTimeout(200);
    // The read-only balloon auto-dismisses quickly, so it is observed rather than read after a
    // delay. Observe document.body, NOT .d4-balloon-container: that container is created lazily
    // with the first balloon, so observing it directly would attach to nothing.
    await page.evaluate(() => {
      const w = window as any;
      w.__balloonTexts = [];
      if (w.__balloonObs) w.__balloonObs.disconnect();
      w.__balloonObs = new MutationObserver((muts: any) => {
        for (const m of muts) for (const n of m.addedNodes) {
          if (n.nodeType !== 1) continue;
          const el = n as Element;
          const cls = typeof el.className === 'string' ? el.className : '';
          if (cls.indexOf('balloon') >= 0) w.__balloonTexts.push((el.textContent || '').trim());
          const inner = el.querySelector ? el.querySelector('[class*="balloon"]') : null;
          if (inner) w.__balloonTexts.push((inner.textContent || '').trim());
        }
      });
      w.__balloonObs.observe(document.body, {childList: true, subtree: true});
    });
    await dblClickCell(page, 'AGE', 3, 0);
    await page.waitForTimeout(500);
    const r = await page.evaluate(() => {
      const w = window as any;
      w.__balloonObs?.disconnect();
      return {
        editorOpen: document.querySelectorAll('input.d4-value-editor').length > 0,
        balloonTexts: (w.__balloonTexts ?? []).filter((t: string) => t.length > 0),
      };
    });
    expect(r.editorOpen).toBe(false); // read-only: no editor opened (GROK-20010 guard)
    expect(r.balloonTexts.some((t: string) => /read-only/i.test(t))).toBe(true); // rejection is signalled by a balloon
    await page.evaluate(() => { grok.shell.tv.grid.props.allowEdit = true; });
  });

  await softStep('Step 13 — after re-enabling edit, typing a digit replaces the cell value', async () => {
    // The digit-input handler opens the editor on grid.currentGridCell, so seed it first.
    await clickCellSettled(page, 'AGE', 5);
    const tableIdx = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
    await page.keyboard.press('7');
    await page.waitForTimeout(250);
    expect(await editorCount(page)).toBe(1); // typing started the editor
    await page.keyboard.press('Enter');
    await page.waitForTimeout(350);
    const age = await page.evaluate((ti) => grok.shell.tv.dataFrame.col('AGE').get(ti), tableIdx);
    expect(age).toBe(7); // the typed digit became the cell value
  });

  // === Scenario 4: clipboard — single-cell Ctrl+Shift+C and multi-row Ctrl+C =====

  await reopenClean(page);

  await softStep('Step 16 — Ctrl+Shift+C copies the single current-cell value', async () => {
    await clickCellSettled(page, 'AGE', 3);
    const cellVal = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return String(df.col(df.currentCol.name).get(df.currentRowIdx));
    });
    await resetCopyCapture(page);
    await page.keyboard.press('Control+Shift+C');
    await page.waitForTimeout(400);
    const captured = await readCopyCapture(page);
    expect(captured.length).toBeGreaterThan(0);        // a copy fired
    expect(captured[captured.length - 1]).toBe(cellVal); // clipboard text equals the single cell value
    expect(/[\t\n]/.test(captured[captured.length - 1] ?? '')).toBe(false); // no tab or newline
  });

  await softStep('Step 18 — Ctrl+C on a 5-row selection copies without error', async () => {
    // The copied block's content is checked manually (grid-ui.md): the Dart clipboard path leaves
    // no execCommand/setData/writeText payload, so the automated guard here is the error channel.
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      for (let i = 0; i < 5; i++) df.selection.set(i, true);
    });
    await page.waitForTimeout(150);
    const selCount = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(selCount).toBe(5); // exactly 5 rows selected before the copy
    const errBefore = consoleErrors.length + pageErrors.length;
    await focusGrid(page);
    await page.keyboard.press('Control+c');
    await page.waitForTimeout(400);
    expect(consoleErrors.length + pageErrors.length).toBe(errBefore); // multi-row copy raised no error
  });

  // === Scenario 5: Ctrl+A → Ctrl+C → Ctrl+V does not error (GROK-20010) =========

  await reopenClean(page);

  await softStep('Step 21 — Ctrl+A → Ctrl+C → Ctrl+V: no console error, row count unchanged (GROK-20010)', async () => {
    const rowsBefore = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
    const errBefore = consoleErrors.length + pageErrors.length;
    await focusGrid(page);
    await page.keyboard.press('Control+a');
    await page.waitForTimeout(150);
    await page.keyboard.press('Control+c');
    await page.waitForTimeout(200);
    await page.keyboard.press('Control+v');
    await page.waitForTimeout(600);
    const rowsAfter = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
    const concurrentModErrors = consoleErrors.concat(pageErrors)
      .filter((e) => /concurrent modification|during iteration/i.test(e));
    expect(concurrentModErrors).toEqual([]);                        // GROK-20010 guard: no concurrent-modification error
    expect(consoleErrors.length + pageErrors.length).toBe(errBefore); // console-error delta is 0
    expect(rowsAfter).toBe(rowsBefore);                              // row count unchanged
  });

  // === Scenario 6: Ctrl+V pastes into a different cell ==========================

  await reopenClean(page);

  await softStep('Step 23 — Ctrl+V pastes the copied value into a different cell', async () => {
    // The paste round-trip works because the onPaste event's clipboardData is delivered under
    // trusted key input, even though navigator.clipboard.readText is denied in this context.
    await clickCellSettled(page, 'AGE', 1);
    const srcVal = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      return df.col('AGE').get(df.currentRowIdx);
    });
    await page.keyboard.press('Control+c'); // copy the single current cell
    await page.waitForTimeout(250);
    await clickCellSettled(page, 'AGE', 10);
    const dstIdx = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
    // Source and target must hold distinct values, or a passing assert would prove nothing.
    const dstBefore = await page.evaluate((i) => grok.shell.tv.dataFrame.col('AGE').get(i), dstIdx);
    expect(dstBefore).not.toBe(srcVal);
    await page.keyboard.press('Control+v'); // paste into the different cell
    await page.waitForTimeout(500);
    const dstAfter = await page.evaluate((i) => grok.shell.tv.dataFrame.col('AGE').get(i), dstIdx);
    expect(dstAfter).toBe(srcVal); // the target cell now holds the copied value
  });

  // === Scenario 7: Shift+Del deletes selected rows; Ctrl+Z restores them ========

  await reopenClean(page);

  await softStep('Step 27 — Shift+Del deletes the 3 selected rows', async () => {
    const rowsBefore = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      df.selection.set(3, true); df.selection.set(5, true); df.selection.set(7, true);
      return df.rowCount;
    });
    await page.waitForTimeout(150);
    await focusGrid(page);
    await page.keyboard.press('Shift+Delete');
    await page.waitForTimeout(500);
    const rowsAfter = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
    expect(rowsAfter).toBe(rowsBefore - 3); // exactly the 3 selected rows removed
    (page as any).__rowsBeforeDelete = rowsBefore;
  });

  await softStep('Step 28 — Ctrl+Z restores the deleted rows', async () => {
    const rowsBefore = (page as any).__rowsBeforeDelete;
    await focusGrid(page);
    await page.keyboard.press('Control+z');
    await page.waitForTimeout(600);
    const rowsAfter = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
    expect(rowsAfter).toBe(rowsBefore); // undo restored the row count (Command.runUndoable)
  });

  // Scenario 8 (addNewRowOnLastRowEdit auto-append) is manual-only and lives in grid-ui.md: the
  // append needs the editor still attached at commit, and its post-open teardown detaches it on
  // the next redraw, ahead of every commit path automation can drive.

  v.finishSpec();
});
