/**
 * Playwright helpers for the PredictiveModelingView column-picker dialog
 * (`[name="dialog-Select-columns..."]`) and the column-selector popup
 * (`.d4-column-selector`).
 *
 * Why this module exists:
 * The picker renders the column list into a canvas (not DOM checkboxes)
 * and the JS-API setter `view.predict = ...` / `view.features = [...]`
 * does NOT visibly update the form — see grok-browser/references/models.md
 * "Features column-picker dialog" gotchas. Driving the canvas through
 * trusted Playwright clicks + search-filter narrowing is the only stable
 * surface; the four Models specs each re-implemented the same recipe.
 *
 * Pattern, distilled from chemprop-spec / models-validators-edge-spec:
 *   1. open the picker (mousedown on `[name="div-Features"]`)
 *   2. clear via `[name="label-None"]`
 *   3. for each column: type its name into the dialog's search input,
 *      then trusted-click the canvas overlay at a derived row coordinate
 *      (right-edge x ≈ checkbox column, swept dy across plausible row
 *      strides until the "N checked" count increments)
 *   4. clear the search filter, click OK
 *   5. as last resort fall back to `[name="label-All"]` if the per-name
 *      sweep failed to toggle (training pipeline tolerates extra cols)
 */
import {Page, expect} from '@playwright/test';

export async function setPredict(page: Page, columnName: string): Promise<void> {
  const current = await page.evaluate(() => (window as any).grok.shell.v.root
    .querySelector('[name="input-host-Predict"] .d4-column-selector-column')?.textContent?.trim());
  if (current === columnName) return;
  await page.evaluate(() => {
    const root = (window as any).grok.shell.v.root;
    const sel = root.querySelector('[name="input-host-Predict"] .d4-column-selector') as HTMLElement;
    sel.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, cancelable: true, view: window, button: 0}));
  });
  await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
    null, {timeout: 5_000});
  await page.evaluate(() => (document.querySelector('.d4-column-selector-backdrop') as HTMLElement).focus());
  await page.keyboard.type(columnName);
  await page.waitForTimeout(300);
  await page.keyboard.press('Enter').catch(() => {});
  await page.waitForTimeout(500);
  await expect.poll(async () => await page.evaluate(() =>
    (window as any).grok.shell.v.root
      .querySelector('[name="input-host-Predict"] .d4-column-selector-column')?.textContent?.trim()),
    {timeout: 10_000}).toBe(columnName);
}

async function openFeaturesPicker(page: Page): Promise<void> {
  await page.evaluate(() => {
    const editor = (window as any).grok.shell.v.root
      .querySelector('[name="div-Features"]') as HTMLElement;
    const r = editor.getBoundingClientRect();
    const ev = (t: string) => new MouseEvent(t, {
      bubbles: true, cancelable: true, view: window, button: 0,
      clientX: r.left + 10, clientY: r.top + 5,
    });
    editor.dispatchEvent(ev('mousedown'));
    editor.dispatchEvent(ev('mouseup'));
    editor.dispatchEvent(ev('click'));
  });
  await page.locator('[name="dialog-Select-columns..."]').waitFor({timeout: 10_000});
}

async function readCheckedCount(page: Page): Promise<number> {
  return await page.locator('[name="dialog-Select-columns..."]').evaluate((d) => {
    const lbl = Array.from(d.querySelectorAll('label'))
      .find((l) => /\d+ checked/.test(l.textContent || ''));
    const m = lbl?.textContent?.match(/^(\d+)/);
    return m ? parseInt(m[1], 10) : 0;
  });
}

/**
 * Toggle the given feature columns ON in the Features column-picker.
 * Idempotent re-open safe: clears via label-None first.
 *
 * Recipe (verified 2026-06-09 against dev build; canvas hit-zone is
 * viewport-independent — same coords on 1920×1080 and 1366×768):
 *   - Search-filter REORDERS the canvas (it is a real filter, not a
 *     visual decorator), so after `Search=<name>` row 0 IS the target
 *     column.
 *   - Trusted click at `(box.right - 40, box.top + 36)` toggles row 0.
 *   - One click per column — no (dx, dy) sweep needed.
 * Falls back to label-All if a click does not toggle (build drift / new
 * picker layout).
 */
export async function selectFeaturesByName(page: Page, columnNames: string[]): Promise<void> {
  await openFeaturesPicker(page);
  const dlg = page.locator('[name="dialog-Select-columns..."]');
  await dlg.locator('[name="label-None"]').click();
  await page.waitForTimeout(300);

  const search = dlg.locator('input[placeholder="Search"]');
  const overlay = dlg.locator('[name="viewer-Grid"] [name="overlay"]');

  for (const name of columnNames) {
    const before = await readCheckedCount(page);
    await search.fill('');
    await page.waitForTimeout(200);
    await search.fill(name);
    await page.waitForTimeout(500);
    const box = await overlay.boundingBox();
    if (!box) throw new Error('column-picker overlay canvas not measurable');
    await page.mouse.click(box.x + box.width - 40, box.y + 36);
    await page.waitForTimeout(300);
    const after = await readCheckedCount(page);
    if (after !== before + 1) {
      // Build drift: row-0 hit-zone moved. Bail out to label-All —
      // training pipeline tolerates extra columns; the load-bearing
      // assertion is "at least the named ones are checked".
      await search.fill('');
      await dlg.locator('[name="label-All"]').click();
      await page.waitForTimeout(300);
      break;
    }
  }
  await search.fill('');
  await page.waitForTimeout(200);
  await dlg.locator('[name="button-OK"]').click();
  await page.waitForFunction(() => !document.querySelector('[name="dialog-Select-columns..."]'),
    null, {timeout: 10_000});
}

/**
 * Open the Features picker, then uncheck a single already-checked row at
 * a fixed canvas-row geometry. Use case: train-spec Block 2 uncheck-WEIGHT
 * (Block 1 left the picker with HEIGHT and WEIGHT checked; the predict
 * column WEIGHT must be removed from features).
 *
 * Row-geometry per grok-browser/references/models.md "Canvas-overlay
 * click hit-zone (Features picker)" verified 2026-06-09 against dev:
 *   - data row stride = 28 px (build-frozen)
 *   - first data row centre = top + 36
 *   - row N centre = top + 36 + 28 * N
 *   - hit x = right - 40 (checkbox column)
 *
 * Earlier `stride = canvas.height / totalRows` derivation was WRONG —
 * the canvas is 300 px tall regardless of row count, so dividing by
 * row count puts the click target well past the actual row band (for
 * totalRows=2 the formula returned y=448, while the real row-1 centre
 * is y=287 — a 160px miss into the empty space below the last row).
 * The `totalRows` arg is kept for API compatibility but is no longer
 * load-bearing.
 *
 * Re-open invariant (verified 2026-06-09): CHECKED rows float to the
 * top of the picker on re-open, in original-column-index order.
 * Re-open needs ≥1200 ms settle before the first click — at 700 ms the
 * row-0 click misses ~30% of the time.
 */
export async function uncheckFeatureRowByPosition(
  page: Page, rowIndex: number, _totalRows?: number): Promise<void> {
  await openFeaturesPicker(page);
  await page.waitForTimeout(1_200);
  const dlg = page.locator('[name="dialog-Select-columns..."]');
  const target = await page.evaluate((idx) => {
    const overlay = document.querySelector(
      '[name="dialog-Select-columns..."] canvas[name="overlay"]') as HTMLCanvasElement;
    const r = overlay.getBoundingClientRect();
    return {x: r.left + r.width - 40, y: r.top + 36 + 28 * idx};
  }, rowIndex);
  await page.mouse.click(target.x, target.y);
  await page.waitForTimeout(400);
  await dlg.locator('[name="button-OK"]').click();
  await page.waitForFunction(() => !document.querySelector('[name="dialog-Select-columns..."]'),
    null, {timeout: 10_000});
}
