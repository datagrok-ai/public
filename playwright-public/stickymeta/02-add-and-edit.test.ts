import { test, expect } from '@playwright/test';
import * as H from './helpers';

// Add and edit sticky metadata on the SPGI molecule grid: a single cell, sticky columns, and a
// multi-row batch edit. The schema is created/deleted via the API (setup/cleanup); every metadata
// action is performed through the UI.
test.describe.configure({ mode: 'serial' });

test('Sticky Meta: add & edit metadata (cell, sticky column, batch)', async ({ page }) => {
  test.setTimeout(300_000);

  const suffix = H.uniqueSuffix();
  const schemaName = `PW_SM_Schema_${suffix}`;

  try {
    await H.gotoHome(page);
    await H.setupEnv(page);
    await H.apiDeleteAllTestSchemas(page); // defensive: no leftover molecule schema from a crashed run
    await H.apiCreateSchema(page, schemaName, [
      { name: 'rating', type: 'int' },
      { name: 'notes', type: 'string' },
      { name: 'verified', type: 'bool' },
      { name: 'review_date', type: 'datetime' },
    ]);
    await H.openSpgi(page);

    // Sanity: the molecule column exists.
    const info = await page.evaluate(() => ({
      rows: (window as any).grok.shell.t.rowCount,
      semType: (window as any).grok.shell.t.col('Structure')?.semType,
    }));
    expect(info.rows).toBe(3624);
    expect(info.semType).toBe('Molecule');

    // ---- 2.1 Add metadata to a single cell via right-click > Sticky meta > Edit for current cell ----
    const openCellEditor = () => page.evaluate(async () => {
      const g = (window as any).grok;
      const grid = g.shell.tv.grid;
      grid.scrollToCell('Structure', 0);
      await new Promise((r) => setTimeout(r, 400));
      g.shell.t.currentCell = g.shell.t.cell(0, 'Structure');
      const gridEl = document.querySelector('[name="viewer-Grid"]') as HTMLElement;
      const rect = gridEl.getBoundingClientRect();
      const b = grid.cell('Structure', 0).bounds;
      const x = rect.left + b.x + b.width / 2;
      const y = rect.top + b.y + b.height / 2;
      (document.elementFromPoint(x, y) as HTMLElement).dispatchEvent(new MouseEvent('contextmenu',
        { bubbles: true, cancelable: true, clientX: x, clientY: y, button: 2, buttons: 2 }));
      await new Promise((r) => setTimeout(r, 800));
      const menu = document.querySelector('.d4-menu-popup');
      const edit = Array.from(menu?.querySelectorAll('.d4-menu-item-label') ?? [])
        .find((l) => l.textContent?.trim() === 'Edit for current cell...');
      (edit?.closest('.d4-menu-item') as HTMLElement)?.click();
      await new Promise((r) => setTimeout(r, 1500));
      return document.querySelector('.d4-dialog .d4-dialog-header, .d4-dialog .d4-dialog-title')?.textContent?.trim();
    });

    expect(await openCellEditor()).toBe('Sticky meta');

    const readCellDialog = () => page.evaluate(() => {
      const d = document.querySelector('.d4-dialog') as HTMLElement;
      return {
        rating: (d.querySelector('[name="input-host-Rating"] input') as HTMLInputElement).value,
        notes: (d.querySelector('[name="input-host-Notes"] input') as HTMLInputElement).value,
        verified: (d.querySelector('[name="input-host-Verified"] input[type="checkbox"]') as HTMLInputElement).checked,
      };
    });
    const closeDialog = async () => {
      await page.evaluate(() => (document.querySelector('.d4-dialog [name="button-CANCEL"]') as HTMLElement | null)?.click());
      await page.waitForTimeout(500);
    };
    const expected = { rating: '5', notes: 'test note', verified: true };

    // Fill rating, notes, verified.
    await page.locator('.d4-dialog [name="input-host-Rating"] input').click();
    await page.keyboard.press('Control+a');
    await page.keyboard.type('5');
    await page.keyboard.press('Tab');
    await page.locator('.d4-dialog [name="input-host-Notes"]').first().locator('input').click();
    await page.keyboard.press('Control+a');
    await page.keyboard.type('test note');
    await page.keyboard.press('Tab');
    await page.evaluate(() => {
      const cb = document.querySelector('.d4-dialog [name="input-host-Verified"] input[type="checkbox"]') as HTMLInputElement;
      if (!cb.checked) cb.click();
    });
    await page.waitForTimeout(400);
    // Confirm the dialog holds the values before saving (catches fill failures distinctly).
    expect(await readCellDialog()).toEqual(expected);

    // Click the schema section's Save and give the server round-trip time to commit.
    await page.evaluate(() => (document.querySelector('.d4-dialog [name="button-Save"]') as HTMLElement | null)?.click());
    await page.waitForTimeout(2500);
    await closeDialog();

    // Re-open the editor and verify persistence; poll because the save round-trip can lag the reopen.
    await expect.poll(async () => {
      expect(await openCellEditor()).toBe('Sticky meta');
      const v = await readCellDialog();
      await closeDialog();
      return JSON.stringify(v);
    }, { timeout: 30_000, intervals: [2000] }).toBe(JSON.stringify(expected));

    // ---- 2.2 Create a sticky column from the Context Panel Sticky meta pane ----
    await page.evaluate(() => { (window as any).grok.shell.o = (window as any).grok.shell.t.col('Structure'); });
    await page.waitForTimeout(1500);
    const stickyHeader = page.locator('.grok-prop-panel .d4-accordion-pane-header', { hasText: 'Sticky meta' }).first();
    await stickyHeader.waitFor({ timeout: 10_000 });
    if (!(await stickyHeader.evaluate((el) => el.classList.contains('expanded')))) await stickyHeader.click();
    await page.waitForTimeout(1500);

    // Click "add all properties as columns" inside our schema's section (scoped by header text).
    await page.evaluate((schemaName) => {
      const section = Array.from(document.querySelectorAll('.grok-prop-panel .d4-build-root.ui-form'))
        .find((s) => s.querySelector('.d4-flex-row')?.textContent?.trim() === schemaName);
      (section?.querySelector('[name$="-properties-as-columns"]') as HTMLElement | null)?.click();
    }, schemaName);

    // Schema matching is async — poll for the four sticky columns to appear.
    await expect.poll(async () =>
      page.evaluate(() => {
        const names = (window as any).grok.shell.t.columns.names();
        return ['rating', 'notes', 'verified', 'review_date'].every((n) => names.includes(n));
      }), { timeout: 20_000, intervals: [500] }).toBe(true);

    // Sort ascending by the rating sticky column.
    // SCOPE NOTE: the grid is canvas-rendered with no stable DOM handle for the column header sort
    // gesture, so the sort is issued through the grid API; it is incidental to the metadata checks.
    await page.evaluate(async () => {
      (window as any).grok.shell.tv.grid.sort(['rating'], [true]);
      await new Promise((r) => setTimeout(r, 800));
    });

    // Remove the rating column, then re-add it from the pane (UI) and verify the value persists.
    // SCOPE NOTE: column removal uses the API (no clean DOM gesture on the canvas grid); the
    // metadata-preserving RE-ADD is the UI action under test.
    await page.evaluate(async () => {
      (window as any).grok.shell.t.columns.remove('rating');
      await new Promise((r) => setTimeout(r, 600));
    });
    expect(await page.evaluate(() => (window as any).grok.shell.t.columns.contains('rating'))).toBe(false);

    await page.evaluate((schemaName) => {
      const section = Array.from(document.querySelectorAll('.grok-prop-panel .d4-build-root.ui-form'))
        .find((s) => s.querySelector('.d4-flex-row')?.textContent?.trim() === schemaName);
      (section?.querySelector('[name="button-Add-rating-as-a-column"]') as HTMLElement | null)?.click();
    }, schemaName);
    await expect.poll(async () =>
      page.evaluate(() => (window as any).grok.shell.t.columns.contains('rating')),
      { timeout: 20_000, intervals: [500] }).toBe(true);

    // The cell edited in 2.1 (row 0) still carries rating = 5 after remove + re-add.
    const ratingRow0 = await page.evaluate(() => (window as any).grok.shell.t.col('rating').get(0));
    expect(ratingRow0).toBe(5);

    // ---- 2.3 Batch edit metadata on multiple rows ----
    const coords = await page.evaluate(async () => {
      const g = (window as any).grok;
      const df = g.shell.t;
      try { g.shell.tv.grid.sort([], []); } catch { /* ignore */ }
      df.selection.setAll(false);
      df.selection.set(10, true);
      df.selection.set(11, true);
      df.selection.set(12, true);
      df.currentRowIdx = 10;
      const grid = g.shell.tv.grid;
      grid.scrollToCell('Structure', 10);
      await new Promise((r) => setTimeout(r, 800));
      const gridEl = document.querySelector('[name="viewer-Grid"]') as HTMLElement;
      const rect = gridEl.getBoundingClientRect();
      const b = grid.cell('Structure', 10).bounds;
      return { x: Math.round(rect.left + b.x + b.width / 2), y: Math.round(rect.top + b.y + b.height / 2) };
    });
    await page.mouse.move(coords.x, coords.y);
    await page.mouse.click(coords.x, coords.y, { button: 'right' });
    await page.locator('.d4-menu-popup').first().waitFor({ timeout: 5000 });
    await page.evaluate(() => {
      const menu = document.querySelector('.d4-menu-popup');
      const label = Array.from(menu?.querySelectorAll('.d4-menu-item-label') ?? [])
        .find((l) => l.textContent?.trim() === 'Edit for all properties');
      (label?.closest('.d4-menu-item') as HTMLElement)?.click();
    });
    await page.locator('[name="input-host-Rows"]').waitFor({ timeout: 10_000 });

    // Set verified = true, notes = "batch note" on the selected rows.
    const verified = page.locator('input[name="prop-view-verified"]').first();
    await verified.waitFor({ timeout: 5000 });
    if (!(await verified.isChecked())) await verified.click();
    await page.waitForTimeout(400);

    // Enter edit mode on the notes value cell, then type. The cell is a view-label until clicked.
    const notesLabel = page.locator('[name="prop-view-notes"].property-grid-item-view-label').first();
    await notesLabel.click();
    const notesInput = page.locator('input.property-grid-item-editor-textbox').first();
    if (!(await notesInput.isVisible().catch(() => false))) await notesLabel.dblclick();
    await notesInput.waitFor({ timeout: 5000 });
    await notesInput.focus();
    await page.keyboard.press('Control+a');
    await page.keyboard.press('Delete');
    await page.keyboard.type('batch note');
    await page.keyboard.press('Enter');
    await page.waitForTimeout(800);

    const vals = await page.evaluate(() => {
      const df = (window as any).grok.shell.tables[0];
      const notes = df.col('notes'), verified = df.col('verified');
      const out: any[] = [];
      for (const i of [10, 11, 12]) out.push({ notes: notes.get(i), verified: verified.get(i) });
      return out;
    });
    expect(vals.length).toBe(3);
    for (const v of vals) {
      expect(v.notes).toBe('batch note');
      expect(v.verified).toBe(true);
    }
  } finally {
    await H.apiDeleteSchema(page, schemaName).catch(() => {});
    await page.evaluate(() => (window as any).grok.shell.closeAll()).catch(() => {});
  }
});
