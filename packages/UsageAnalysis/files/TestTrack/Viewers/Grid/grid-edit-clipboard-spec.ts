/* ---
realizes: [grid.cp.edit-clipboard]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {isLocalBootNoise, openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import * as g from './grid-helpers';

declare const grok: any;

test.use(specTestOptions);

async function armEvent(page: Page, source: 'df' | 'grid', stream: string): Promise<void> {
  await page.evaluate(({source, stream}) => {
    const w = window as any;
    const tv = grok.shell.tv;
    const obj = source === 'df' ? tv.dataFrame : tv.grid;
    if (w.__armSub) { try { w.__armSub.unsubscribe(); } catch (_) {} }
    w.__armFired = false;
    w.__armSub = obj[stream].subscribe(() => { w.__armFired = true; });
  }, {source, stream});
}

async function awaitEvent(page: Page, capMs: number): Promise<void> {
  await page.evaluate(async (capMs) => {
    const w = window as any;
    const t0 = Date.now();
    while (!w.__armFired && Date.now() - t0 < capMs) await new Promise((r) => setTimeout(r, 25));
    try { w.__armSub?.unsubscribe(); } catch (_) {}
    w.__armSub = null;
  }, capMs);
}

async function waitEditorCount(page: Page, expected: number, capMs = 2000): Promise<number> {
  return page.evaluate(async ({expected, capMs}) => {
    const t0 = Date.now();
    const read = () => document.querySelectorAll('input.d4-value-editor').length;
    while (read() !== expected && Date.now() - t0 < capMs) await new Promise((r) => setTimeout(r, 25));
    return read();
  }, {expected, capMs});
}

async function editorCount(page: Page): Promise<number> {
  return page.evaluate(() => document.querySelectorAll('input.d4-value-editor').length);
}

async function clickCellSettled(page: Page, col: string, visualRow: number): Promise<{x: number; y: number}> {
  const p = await g.cellCenter(page, col, visualRow);
  const tableRow = await page.evaluate((vr) => grok.shell.tv.grid.gridRowToTable(vr), visualRow);
  await g.focusGrid(page);
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

async function dblClickCell(page: Page, col: string, visualRow: number, attempts = 6): Promise<void> {
  const p = await clickCellSettled(page, col, visualRow);
  await page.mouse.dblclick(p.x, p.y);
  if (attempts === 0) return;
  await waitEditorCount(page, 1, 800);
  for (let i = 0; i < attempts && (await editorCount(page)) === 0; i++) {
    await clickCellSettled(page, col, visualRow);
    await page.keyboard.press('Enter');
    await waitEditorCount(page, 1, 800);
  }
}

async function installCopyInterceptor(page: Page): Promise<void> {
  await page.evaluate(() => {
    const w = window as any;
    w.__copyCaptured = [];
    if (!w.__copyPatched) {
      w.__origExecCommand = document.execCommand;
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

async function removeCopyInterceptor(page: Page): Promise<void> {
  await page.evaluate(() => {
    const w = window as any;
    if (w.__copyPatched) { (document as any).execCommand = w.__origExecCommand; w.__copyPatched = false; }
    if (w.__editSub) { w.__editSub.unsubscribe(); w.__editSub = null; }
  }).catch(() => {});
}

async function resetCopyCapture(page: Page): Promise<void> {
  await page.evaluate(() => { (window as any).__copyCaptured = []; });
}

async function primeView(page: Page): Promise<void> {
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
  });
  await installCopyInterceptor(page);
}

async function reopenClean(page: Page): Promise<void> {
  await v.closeAllAndWait(page);
  await v.openTable(page, {path: g.DEMOG, semTypeTimeoutMs: 3000});
  await primeView(page);
}

async function editCount(page: Page): Promise<number> {
  return page.evaluate(() => (window as any).__editCount ?? 0);
}

test('Grid — Cell Editing and Clipboard', async ({page}) => {
  test.setTimeout(240_000);

  await openDatagrok(page);
  const flags = await g.readShellFlags(page);
  const errors = g.trackErrors(page, isLocalBootNoise);
  try {
    await v.openTable(page, {path: g.DEMOG, semTypeTimeoutMs: 3000});
    await primeView(page);

    await softStep('Step 4 — double-click AGE cell, type a new value, Enter commits it', async () => {
      const tableIdx = await page.evaluate(() => grok.shell.tv.grid.gridRowToTable(3));
      await dblClickCell(page, 'AGE', 3);
      expect(await waitEditorCount(page, 1)).toBe(1);
      const editBefore = await editCount(page);
      await page.keyboard.press('Control+a');
      await page.keyboard.type('99');
      await armEvent(page, 'grid', 'onCellValueEdited');
      await page.keyboard.press('Enter');
      await awaitEvent(page, 1500);
      // committing the value and closing the editor are two stages of one gesture: wait on both
      const settled = await v.pollValue(() => page.evaluate((ti) => {
        const df = grok.shell.tv.dataFrame;
        return `${df.col('AGE').get(ti)}|${document.querySelectorAll('input.d4-value-editor').length}`;
      }, tableIdx), (x: string) => x === '99|0', 3000, 50);
      const r = {age: Number(settled.split('|')[0]), editorOpen: settled.split('|')[1] !== '0'};
      const editAfter = await editCount(page);
      expect(r.age).toBe(99);
      expect(editAfter).toBe(editBefore + 1);
      expect(r.editorOpen).toBe(false);
    });

    await softStep('Step 6 — re-open the editor, Esc leaves the value unchanged', async () => {
      const tableIdx = await page.evaluate(() => grok.shell.tv.grid.gridRowToTable(3));
      const before = await page.evaluate((ti) => grok.shell.tv.dataFrame.col('AGE').get(ti), tableIdx);
      await dblClickCell(page, 'AGE', 3);
      expect(await waitEditorCount(page, 1)).toBe(1);
      await page.keyboard.press('Escape');
      expect(await waitEditorCount(page, 0)).toBe(0);
      const age = await page.evaluate((ti) => grok.shell.tv.dataFrame.col('AGE').get(ti), tableIdx);
      expect(age).toBe(before);
    });

    await reopenClean(page);

    await softStep('Step 8 — Delete (no modifier) clears the cell, opens no editor', async () => {
      await clickCellSettled(page, 'AGE', 4);
      let cleared = false;
      for (let attempt = 0; attempt < 3 && !cleared; attempt++) {
        if (attempt > 0) await clickCellSettled(page, 'AGE', 4);
        const row = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
        await g.focusGrid(page);
        const focused = await page.evaluate(() =>
          document.activeElement === document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]'));
        expect(focused).toBe(true);
        await armEvent(page, 'df', 'onValuesChanged');
        await page.keyboard.press('Delete');
        await awaitEvent(page, 2000);
        cleared = await page.evaluate((ti) => grok.shell.tv.dataFrame.col('AGE').isNone(ti), row);
      }
      expect(cleared).toBe(true);
      expect(await editorCount(page)).toBe(0);
    });

    await reopenClean(page);

    await softStep('Step 11 — allowEdit=false: double-click opens no editor + a read-only balloon', async () => {
      await page.evaluate(() => { grok.shell.tv.grid.props.allowEdit = false; });

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

      const balloonTexts = await page.evaluate(async () => {
        const w = window as any;
        const t0 = Date.now();
        const has = () => (w.__balloonTexts ?? []).some((t: string) => /read-only/i.test(t));
        while (!has() && Date.now() - t0 < 1500) await new Promise((r) => setTimeout(r, 25));
        w.__balloonObs?.disconnect();
        return (w.__balloonTexts ?? []).filter((t: string) => t.length > 0);
      });
      expect(await editorCount(page)).toBe(0);
      expect(balloonTexts.some((t: string) => /read-only/i.test(t))).toBe(true);
      await page.evaluate(() => { grok.shell.tv.grid.props.allowEdit = true; });
    });

    await softStep('Step 13 — after re-enabling edit, typing a digit replaces the cell value', async () => {
      await clickCellSettled(page, 'AGE', 5);
      const tableIdx = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);
      await page.keyboard.press('7');
      expect(await waitEditorCount(page, 1)).toBe(1);
      await armEvent(page, 'grid', 'onCellValueEdited');
      await page.keyboard.press('Enter');
      await awaitEvent(page, 1500);
      const age = await page.evaluate((ti) => grok.shell.tv.dataFrame.col('AGE').get(ti), tableIdx);
      expect(age).toBe(7);
    });

    await reopenClean(page);

    await softStep('Step 16 — Ctrl+Shift+C copies the single current-cell value', async () => {
      await clickCellSettled(page, 'AGE', 3);
      const cellVal = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        return String(df.col(df.currentCol.name).get(df.currentRowIdx));
      });
      await resetCopyCapture(page);
      await page.keyboard.press('Control+Shift+C');

      const captured = await page.evaluate(async () => {
        const w = window as any;
        const t0 = Date.now();
        const read = () => (w.__copyCaptured ?? []).filter((x: any) => x != null);
        while (read().length === 0 && Date.now() - t0 < 1500) await new Promise((r) => setTimeout(r, 25));
        return read();
      });
      expect(captured.length).toBeGreaterThan(0);
      expect(captured[captured.length - 1]).toBe(cellVal);
      expect(/[\t\n]/.test(captured[captured.length - 1] ?? '')).toBe(false);
    });

    await softStep('Step 18 — Ctrl+C on a 5-row selection copies without error', async () => {
      await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        df.selection.setAll(false);
        for (let i = 0; i < 5; i++) df.selection.set(i, true);
      });
      const selCount = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
      expect(selCount).toBe(5);
      const errBefore = errors.count();
      await g.focusGrid(page);
      await page.keyboard.press('Control+c');
      // the error window after the copy is the assertion
      await page.waitForTimeout(400);
      expect(errors.count()).toBe(errBefore);
    });

    await reopenClean(page);

    await softStep('Step 21 — Ctrl+A → Ctrl+C → Ctrl+V: no console error, row count unchanged (GROK-20010)', async () => {
      const rowsBefore = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
      const errBefore = errors.count();
      await g.focusGrid(page);
      await page.keyboard.press('Control+a');
      await page.keyboard.press('Control+c');
      await page.keyboard.press('Control+v');
      // the error window after the paste is the assertion
      await page.waitForTimeout(600);
      const rowsAfter = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
      const concurrentModErrors = errors.list.filter((e) => /concurrent modification|during iteration/i.test(e));
      expect(concurrentModErrors).toEqual([]);
      expect(errors.count()).toBe(errBefore);
      expect(rowsAfter).toBe(rowsBefore);
    });

    await reopenClean(page);

    await softStep('Step 23 — Ctrl+V pastes the copied value into a different cell', async () => {
      await clickCellSettled(page, 'AGE', 1);
      const srcVal = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        return df.col('AGE').get(df.currentRowIdx);
      });
      await page.keyboard.press('Control+c');
      await clickCellSettled(page, 'AGE', 10);
      const dstIdx = await page.evaluate(() => grok.shell.tv.dataFrame.currentRowIdx);

      const dstBefore = await page.evaluate((i) => grok.shell.tv.dataFrame.col('AGE').get(i), dstIdx);
      expect(dstBefore).not.toBe(srcVal);
      await armEvent(page, 'df', 'onValuesChanged');
      await page.keyboard.press('Control+v');
      await awaitEvent(page, 2000);
      const dstAfter = await page.evaluate((i) => grok.shell.tv.dataFrame.col('AGE').get(i), dstIdx);
      expect(dstAfter).toBe(srcVal);
    });

    await reopenClean(page);

    let rowsBeforeDelete = 0;
    await softStep('Step 27 — Shift+Del deletes the 3 selected rows', async () => {
      rowsBeforeDelete = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        df.selection.setAll(false);
        df.selection.set(3, true); df.selection.set(5, true); df.selection.set(7, true);
        return df.rowCount;
      });
      await g.focusGrid(page);
      await armEvent(page, 'df', 'onRowsRemoved');
      await page.keyboard.press('Shift+Delete');
      await awaitEvent(page, 2000);
      const rowsAfter = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
      expect(rowsAfter).toBe(rowsBeforeDelete - 3);
    });

    await softStep('Step 28 — Ctrl+Z restores the deleted rows', async () => {
      await g.focusGrid(page);
      await armEvent(page, 'df', 'onRowsAdded');
      await page.keyboard.press('Control+z');
      await awaitEvent(page, 2500);
      const rowsAfter = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
      expect(rowsAfter).toBe(rowsBeforeDelete);
    });
  } finally {
    errors.stop();
    await removeCopyInterceptor(page);
    await g.leaveShellClean(page, flags);
  }
  v.finishSpec();
});
