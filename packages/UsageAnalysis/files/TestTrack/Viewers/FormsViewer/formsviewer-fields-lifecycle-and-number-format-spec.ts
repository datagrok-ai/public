/* ---
realizes: [formsviewer.cp.fields-lifecycle-and-number-format, formsviewer.int.number-format-vs-grid, formsviewer.edge.empty-fields-column-names, formsviewer.edge.tilde-columns-excluded, formsviewer.edge.over-20-columns-capped-silently, formsviewer.edge.renamed-column-follows-rename, formsviewer.edge.number-format-float-only]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import {addViewerByIcon, cleanupShell, finishSpec, openTable, openViewerProperties} from '../../helpers/viewers';
import {HOST, CURRENT, drawnLabelNames, balloonCount, withConsoleErrorCount} from '../../helpers/forms';

declare const grok: any;
declare const DG: any;

// GRID column-menu selectors, absent from grok-browser/references/viewers/formsviewer.md (the refdoc
// had no right-click primitive to reach them):
//   [name="div-Column-Properties..."] — opens the column properties dialog (root
//     [name="dialog-<ColumnName>"], New name input [name="input-New-name--"]).
//   [name="div-Remove"] — drops the column from the dataframe.

const datasetPath = 'System:DemoFiles/demog.csv';

// The field-read helpers below scope every read to the current-row card (CURRENT = green indicator).
async function cardFieldColumns(page: Page): Promise<string[]> {
  return page.evaluate((sel) => {
    const card = document.querySelector(sel) as HTMLElement | null;
    return card ? Array.from(card.querySelectorAll('[column]')).map((e) => e.getAttribute('column') ?? '') : [];
  }, CURRENT);
}

async function cardFieldText(page: Page, column: string): Promise<string | null> {
  return page.evaluate(({sel, col}) => {
    const card = document.querySelector(sel) as HTMLElement | null;
    const el = card?.querySelector(`[column="${col}"]`) as HTMLInputElement | null;
    return el ? (el.value ?? el.textContent ?? '') : null;
  }, {sel: CURRENT, col: column});
}

/**
 * Grid-side rendered text for (column, row): the raw value formatted with the GRID COLUMN's declared
 * format via the platform formatter — the source of truth, computed off the raw double, independent
 * of the viewer's own col.getString(row) path.
 *
 * HONEST-SCOPE NOTE (probe-measured, dev demog row 0): the grid column's display format and
 * col.meta.format are the SAME storage, so a divergence between "as in grid" and "as in column" is
 * UNREACHABLE today. The S2 "Same as grid" step therefore does NOT prove "matches the grid" — only
 * the modest claim that the field renders through a format (not raw precision, not empty). Auto-
 * tracks the grid if a future build ever separates the two formats.
 */
async function gridCellText(page: Page, column: string, row: number): Promise<string> {
  return page.evaluate(({col, r}) => {
    const c = grok.shell.t.col(col);
    if (c.isNone(r))
      return '';
    return DG.format(c.get(r), grok.shell.tv.grid.col(col).format);
  }, {col: column, r: row});
}

/** Centre of a grid column-header cell — the target of a trusted right-click on the canvas. */
async function gridHeaderPoint(page: Page, column: string): Promise<{x: number; y: number}> {
  return page.evaluate((col) => {
    const grid = grok.shell.tv.grid;
    const gc = grid.columns.byName(col);
    const canvas = Array.from(document.querySelectorAll('[name="viewer-Grid"] canvas'))
      .find((c: any) => c.getBoundingClientRect().width > 100) as HTMLElement;
    const r = canvas.getBoundingClientRect();
    return {x: Math.round(r.x + gc.left + gc.width / 2), y: Math.round(r.y + grid.colHeaderHeight / 2)};
  }, column);
}

/**
 * Live centre of a grid context-menu item, or null while absent OR still zero-size. A menu item is
 * in the DOM before it is clickable, so readiness is a real non-zero box, not node presence. Polls up to `ms`.
 */
async function liveMenuItemCentre(page: Page, sel: string, ms: number): Promise<{x: number; y: number} | null> {
  const deadline = Date.now() + ms;
  for (;;) {
    const c = await page.evaluate((s) => {
      const e = document.querySelector(s);
      if (!e) return null;
      const r = e.getBoundingClientRect();
      if (r.width === 0 || r.height === 0) return null;
      return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    }, sel);
    if (c || Date.now() > deadline) return c;
    await page.waitForTimeout(150);
  }
}

/**
 * Rename a dataframe column through the GRID's column-header UI — real trusted right-click on the
 * header canvas, then the Column Properties dialog's "New name" + OK. MANDATORY path: assigning
 * col.name via the JS API is forbidden by the atlas (semType-null errors and a lost field on this build).
 *
 * Exit condition: the new name is in columns.names() and the old is gone. The rename commits on a
 * shared server at unknown latency, so the whole gesture is retried until that holds; a half-open
 * menu is dismissed with Escape between attempts.
 */
async function renameColumnViaGridUI(page: Page, oldName: string, newName: string): Promise<void> {
  await expect.poll(async () => {
    const names = await page.evaluate(() => grok.shell.t.columns.names() as string[]);
    if (names.includes(newName) && !names.includes(oldName))
      return true;
    if (!names.includes(oldName))
      return false; // rename applied to the name but the new one not visible yet — let the poll wait
    await page.keyboard.press('Escape');
    const pt = await gridHeaderPoint(page, oldName);
    await page.mouse.click(pt.x, pt.y, {button: 'right'});
    const item = page.locator('[name="div-Column-Properties..."]').first();
    if (!await item.isVisible().catch(() => false))
      return false;
    await item.click();
    const nameInput = page.locator(`[name="dialog-${oldName}"] [name="input-New-name--"]`);
    if (!await nameInput.isVisible({timeout: 3000}).catch(() => false))
      return false;
    await nameInput.fill(newName);
    await page.locator(`[name="dialog-${oldName}"] [name="button-OK"]`).click();
    return false; // decided next iteration on columns.names() above
  }, {timeout: 30_000, intervals: [300, 400, 600, 800, 1000]}).toBe(true);
}

/**
 * Drop a dataframe column through the GRID's column-header context menu — trusted right-click on the
 * header canvas, then the flat "Remove" item. Exit condition: the column is gone from columns.names().
 * Readiness is a live non-zero box (liveMenuItemCentre), not node presence; the drop commits at
 * unknown latency, so the gesture is retried until it disappears, a half-open menu dismissed with Escape.
 */
async function removeColumnViaGridUI(page: Page, column: string): Promise<void> {
  await expect.poll(async () => {
    if (await page.evaluate((col) => !(grok.shell.t.columns.names() as string[]).includes(col), column))
      return true;
    await page.keyboard.press('Escape');
    const pt = await gridHeaderPoint(page, column);
    await page.mouse.click(pt.x, pt.y, {button: 'right'});
    const centre = await liveMenuItemCentre(page, '[name="div-Remove"]', 3000);
    if (centre)
      await page.mouse.click(centre.x, centre.y);
    return false; // decided next iteration on columns.names() above
  }, {timeout: 20_000, intervals: [200, 300, 400, 600, 800]}).toBe(true);
}

/**
 * Change Number Format through the property-panel choice editor — the real user route. Number Format
 * is a choice property under Misc; clicking its VALUE element (not the row) creates the <select>,
 * whose native input/change events are the real commit path.
 *
 * Two facts make this fiddly, both handled by gating on the LIVE viewer's own numberFormat prop:
 *  - A stale panel from a closed table can linger and move the wrong viewer, so each attempt drives
 *    only the VISIBLE editor and, failing that, gear-clicks the current viewer to bring ITS panel forward.
 *  - The 50 ms-debounced subscriptions rebuild the grid, so the select is re-queried live each attempt.
 */
async function setNumberFormatViaPanel(page: Page, value: string): Promise<void> {
  // Bring the CURRENT viewer's panel forward: hide any stale panel first, then open THIS viewer's
  // gear scoped through the live viewer instance so a lingering [name="viewer-Forms"] cannot shadow
  // it. Element-scoped click (position-independent), so the layout shift cannot make it miss.
  await page.evaluate(() => { grok.shell.windows.showContextPanel = false; });
  await page.waitForTimeout(300);
  await page.evaluate(() => {
    const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
    (vw?.root?.closest('.panel-base')
      ?.querySelector('.panel-titlebar [name="icon-font-icon-settings"]') as HTMLElement | null)?.click();
  });
  await page.locator('[name="prop-number-format"]').first().waitFor({state: 'visible', timeout: 8000}).catch(() => {});
  // Commit the choice on the now-visible editor. The 50 ms rebuild can replace the select, so it is
  // re-queried live each attempt; exit gate is the viewer's own prop reaching the target.
  await expect.poll(async () => page.evaluate((v) => {
    const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
    if (vw.props.numberFormat === v) return true;
    const row = document.querySelector('[name="prop-number-format"]') as HTMLElement | null;
    if (!row) return false;
    // Misc can ship COLLAPSED (rows in DOM but hidden) — expand it, else the choice editor is unreachable.
    if (row.offsetParent === null) {
      (document.querySelector('[name="prop-category-misc"]') as HTMLElement | null)?.click();
      return false;
    }
    (row.querySelector('[name="prop-view-number-format"]') as HTMLElement | null)?.click();
    const sel = row.querySelector('select') as HTMLSelectElement | null;
    if (!sel) return false;
    const opt = Array.from(sel.options).find((o) => o.value === v || (o.textContent ?? '').trim() === v);
    if (!opt) return false;
    sel.value = opt.value;
    sel.dispatchEvent(new Event('input', {bubbles: true}));
    sel.dispatchEvent(new Event('change', {bubbles: true}));
    return false;
  }, value), {timeout: 12_000, intervals: [250, 300, 400, 600, 800]}).toBe(true);
  await page.waitForTimeout(300);
}

/**
 * Remove a field via its HEADER close icon, with a REAL-mouse gesture robust to renderHeader()
 * rebuilding the header node. The icon is visibility:hidden and revealed only by a per-node
 * onmouseenter; renderHeader() rebuilds every container from scratch, so a re-render between hover
 * and click hides it again. Hence: retry the whole bundle, re-query the node each iteration, and fire
 * a real move onto the label → real click on the icon centre back-to-back against the 50 ms debounce.
 * Exit condition is the field leaving fieldsColumnNames, so a transient empty header never reports a
 * false "removed". The icon keeps a real ~8×13 px box even while hidden, so its centre is a valid target.
 */
async function removeFieldViaHeaderCloseIcon(page: Page, column: string): Promise<void> {
  await expect.poll(async () => {
    if (!await page.evaluate((col) => grok.shell.tv.viewers
      .find((x: any) => x.type === 'FormsViewer').fieldsColumnNames.includes(col), column))
      return true;

    const pts = await page.evaluate(({host, col}) => {
      const container = Array.from(document.querySelectorAll(
        `${host} .d4-multi-form-header .d4-multi-form-column-name`))
        .find((c) => c.querySelector(`div[name="div-${col}"]`)) as HTMLElement | undefined;
      if (!container)
        return null;
      const icon = container.querySelector('i.grok-icon.fal.fa-times') as HTMLElement | null;
      if (!icon)
        return null;
      const cr = container.getBoundingClientRect();
      const ir = icon.getBoundingClientRect();
      if (ir.width === 0 || ir.height === 0)
        return {zero: true as const};
      return {cx: cr.x + cr.width / 2, cy: cr.y + cr.height / 2, ix: ir.x + ir.width / 2, iy: ir.y + ir.height / 2};
    }, {host: HOST, col: column});
    if (pts === null)
      return false;
    // Tripwire: a zero-size icon box a real mouse cannot hit — fail loudly, no synthetic fallback.
    if ('zero' in pts)
      throw new Error(`[removeFieldViaHeaderCloseIcon] '${column}' close icon has a zero-size box — a real mouse cannot target it; investigate before any synthetic fallback.`);

    await page.mouse.move(pts.cx, pts.cy);
    await page.mouse.click(pts.ix, pts.iy);
    return false; // decided next iteration on fieldsColumnNames above
  }, {timeout: 20_000, intervals: [150, 150, 250, 400, 600]}).toBe(true);
}

test.use(specTestOptions);

test('Forms viewer — field lifecycle and number format', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  // #### Setup
  await softStep('Setup — open demog, add COMPUTED_H (no-format FLOAT), attach the Forms viewer', async () => {
    await openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
    // COMPUTED_H = HEIGHT * 1.0 — a FLOAT column with NO format tag, the GROK-20367 entry state
    // (demog's HEIGHT carries meta.format '0.000', invalid for the no-format repro). The guard below fails loudly.
    const entry = await page.evaluate(() => {
      const df = grok.shell.t;
      const h = df.col('HEIGHT');
      const comp = df.columns.addNewFloat('COMPUTED_H');
      for (let i = 0; i < df.rowCount; i++) comp.set(i, h.get(i) * 1.0);
      df.currentRowIdx = 0;
      return {type: comp.type, format: comp.meta.format ?? null};
    });
    expect(entry.type).toBe('double');
    expect(entry.format).toBeNull();

    await addViewerByIcon(page, 'Forms', 'Forms', 30_000, 'FormsViewer');
    await page.locator(HOST).first().waitFor({timeout: 30_000});
    await page.locator('.d4-multi-form').first().waitFor({timeout: 30_000});
    await page.evaluate(() => { grok.shell.t.currentRowIdx = 0; });
  });

  // #### Scenario 1 · Step 3 (GROK-14962)
  await softStep('Scenario 1 Step 3 — Fields render in the picked order RACE, AGE, SEX (not table order)', async () => {
    // The Select Columns picker rows are canvas-drawn (no per-checkbox handle), so the ORDERED subset
    // is set through the property API (documented refdoc fallback). Picked order is deliberately NOT
    // table order so a set-equality build can't pass by coincidence: ORDER is the GROK-14962 signal.
    await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      vw.setOptions({fieldsColumnNames: ['RACE', 'AGE', 'SEX']});
    });
    await expect.poll(() => drawnLabelNames(page), {timeout: 20_000}).toEqual(['RACE', 'AGE', 'SEX']);
    await expect.poll(() => cardFieldColumns(page), {timeout: 20_000}).toEqual(['RACE', 'AGE', 'SEX']);
  });

  // #### Scenario 1 · Step 4
  await softStep('Scenario 1 Step 4 — Removing the AGE field drops its label and elements; order kept', async () => {
    await removeFieldViaHeaderCloseIcon(page, 'AGE');

    await expect.poll(() => drawnLabelNames(page), {timeout: 20_000}).toEqual(['RACE', 'SEX']);
    expect(await page.locator(`${HOST} [column="AGE"]`).count()).toBe(0);
    await expect.poll(() => cardFieldColumns(page), {timeout: 20_000}).toEqual(['RACE', 'SEX']);
  });

  // #### Scenario 1 · Step 5
  await softStep('Scenario 1 Step 5 — Dropping RACE from the dataframe prunes its field with no error', async () => {
    const errCount = await withConsoleErrorCount(page, async () => {
      await removeColumnViaGridUI(page, 'RACE');
    });
    expect(await page.evaluate(() => grok.shell.t.columns.names().includes('RACE'))).toBe(false);
    await expect.poll(() => drawnLabelNames(page), {timeout: 20_000}).toEqual(['SEX']);
    expect(await page.locator(`${HOST} [column="RACE"]`).count()).toBe(0);
    expect(errCount).toBe(0);
    expect(await page.locator('.d4-balloon.error').count()).toBe(0);
  });

  // #### Scenario 1 · Step 6
  await softStep('Scenario 1 Step 6 — Renaming SEX→GENDER via the grid UI: the field follows the new name', async () => {
    const errCount = await withConsoleErrorCount(page, async () => {
      await renameColumnViaGridUI(page, 'SEX', 'GENDER');
    });
    await expect.poll(() => drawnLabelNames(page), {timeout: 20_000}).toEqual(['GENDER']);
    await expect.poll(() => cardFieldColumns(page), {timeout: 20_000}).toEqual(['GENDER']);
    expect(await page.evaluate(() => grok.shell.t.columns.names().includes('GENDER'))).toBe(true);
    expect(errCount).toBe(0);
  });

  // #### Scenario 1 · Step 7
  await softStep('Scenario 1 Step 7 — A >20-column table yields exactly 20 fields with no message', async () => {
    // No wide demo table on dev, so build a synthetic 26-visible-column one (25 numeric + LABEL;
    // '~hidden' excluded by prefix) — the cap claim is dataset-agnostic.
    await cleanupShell(page);
    await page.evaluate(() => {
      const N = 30;
      const cols: any[] = [];
      for (let c = 0; c < 25; c++) {
        const name = 'C' + String(c).padStart(2, '0');
        cols.push(DG.Column.fromInt32Array(name, Int32Array.from(Array.from({length: N}, (_, i) => i + c * 100))));
      }
      cols.push(DG.Column.fromStrings('LABEL', Array.from({length: N}, (_, i) => 'r' + i)));
      cols.push(DG.Column.fromStrings('~hidden', Array.from({length: N}, () => 'h')));
      const df = DG.DataFrame.fromColumns(cols);
      df.name = 'wide26';
      grok.shell.addTableView(df);
      df.currentRowIdx = 0;
    });
    await page.waitForTimeout(600);
    await addViewerByIcon(page, 'Forms', 'Forms', 30_000, 'FormsViewer');
    await page.locator(HOST).first().waitFor({timeout: 30_000});

    await expect.poll(() => drawnLabelNames(page).then((l) => l.length), {timeout: 20_000}).toBe(20);
    expect(await cardFieldColumns(page)).toHaveLength(20);
    expect((await drawnLabelNames(page)).some((n) => n.startsWith('~'))).toBe(false);
    expect(await balloonCount(page)).toBe(0);
    expect(await page.evaluate((host) => (document.querySelector(host) as HTMLElement)?.textContent?.includes('Number of columns') ?? false, HOST)).toBe(false);
  });

  // #### Scenario 1 · Step 8 (GROK-20027)
  await softStep('Scenario 1 Step 8 — Renaming C05→~SERVICE drops it from the field set (no ~ label)', async () => {
    // Exercises the column-RENAME '~'-exclusion path (distinct from initial auto-selection). Capture
    // the drawn-label count BEFORE the rename so the drop is a real one-field decrease, not only the
    // negative "no ~ / no C05" checks (an empty header mid-rebuild satisfies all negatives at once).
    const beforeCount = (await drawnLabelNames(page)).length;
    expect(beforeCount).toBeGreaterThan(0);

    const errCount = await withConsoleErrorCount(page, async () => {
      await renameColumnViaGridUI(page, 'C05', '~SERVICE');
    });
    await expect.poll(() => drawnLabelNames(page).then((l) => l.length), {timeout: 20_000}).toBe(beforeCount - 1);
    const labels = await drawnLabelNames(page);
    expect(labels.some((n) => n.startsWith('~'))).toBe(false);
    expect(labels).not.toContain('C05');
    expect(await page.locator(`${HOST} [column="~SERVICE"]`).count()).toBe(0);
    expect(errCount).toBe(0);
    expect(await page.locator('.d4-balloon.error').count()).toBe(0);
  });

  // #### Scenario 1 · Step 8b — formsviewer.edge.empty-fields-column-names
  await softStep('Scenario 1 Step 8b — Clearing every field draws nothing and raises no error', async () => {
    // Real user route: gear → Fields "..." → Select columns → None → OK (refdoc mandates the dialog;
    // picker rows are canvas-drawn). Empty field set ⇒ zero labels, zero [column], no throw, no error
    // balloon; the floor spans the clear + re-render. Scenario 2 reopens demog fresh, no restore needed.
    await openViewerProperties(page, 'Forms', '[name="prop-view-fields"]');
    const dlg = '[name="dialog-Select-columns..."]';
    const errCount = await withConsoleErrorCount(page, async () => {
      // Readiness is gated on live visibility (non-zero box), not DOM presence — controls resolve before clickable.
      const fieldsButton = page.locator('[name="prop-view-fields"] button');
      await fieldsButton.waitFor({state: 'visible', timeout: 10_000});
      await fieldsButton.click();
      await page.locator(dlg).waitFor({state: 'visible', timeout: 10_000});
      await page.locator(`${dlg} [name="label-None"]`).click();
      await page.locator(`${dlg} [name="button-OK"]`).click();
      await expect.poll(() => drawnLabelNames(page).then((l) => l.length), {timeout: 20_000}).toBe(0);
    });
    expect(await drawnLabelNames(page)).toEqual([]);
    expect(await page.locator(`${HOST} [column]`).count()).toBe(0);
    expect(errCount).toBe(0);
    expect(await page.locator('.d4-balloon.error').count()).toBe(0);
  });

  // #### Scenario 2 · Step 4 (GROK-20367)
  await softStep('Scenario 2 Step 4 — COMPUTED_H field text equals the grid column format applied to the raw value (Same as grid)', async () => {
    await openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
    await page.evaluate(() => {
      const df = grok.shell.t;
      const h = df.col('HEIGHT');
      const comp = df.columns.addNewFloat('COMPUTED_H');
      for (let i = 0; i < df.rowCount; i++) comp.set(i, h.get(i) * 1.0);
      df.currentRowIdx = 0;
    });
    await addViewerByIcon(page, 'Forms', 'Forms', 30_000, 'FormsViewer');
    await page.locator(HOST).first().waitFor({timeout: 30_000});
    // Confirm the ship default BEFORE any write — a set-then-read would be a tautology.
    const defaultFmt = await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      grok.shell.t.currentRowIdx = 0;
      return vw.props.numberFormat;
    });
    expect(defaultFmt).toBe('Same as grid');
    await page.waitForTimeout(800);

    const gridText = await gridCellText(page, 'COMPUTED_H', 0);
    expect(gridText.length).toBeGreaterThan(0);
    await expect.poll(() => cardFieldText(page, 'COMPUTED_H'), {timeout: 15_000}).toBe(gridText);
  });

  // #### Scenario 2 · Step 6
  await softStep('Scenario 2 Step 6 — An explicit float mask changes COMPUTED_H and differs from the grid', async () => {
    // The grid renders this no-format column to 2 decimals (160.48), so a "2 digits" mask would
    // coincide and can't diverge. "3 digits" (160.484) takes a distinct shape, value unchanged.
    const before = await page.evaluate(() => grok.shell.t.col('COMPUTED_H').get(0));
    await setNumberFormatViaPanel(page, '3 digits after comma');
    await expect.poll(() => cardFieldText(page, 'COMPUTED_H'), {timeout: 15_000})
      .toMatch(/^-?\d+\.\d{3}$/);
    const gridText = await gridCellText(page, 'COMPUTED_H', 0);
    expect(await cardFieldText(page, 'COMPUTED_H')).not.toBe(gridText);
    expect(await page.evaluate(() => grok.shell.t.col('COMPUTED_H').get(0))).toBe(before);
  });

  // #### Scenario 2 · Step 7
  await softStep('Scenario 2 Step 7 — INT and string fields render identically before/after the mask change', async () => {
    // demog has no QNUM column, so AGE (int) and SEX (string) are the two non-FLOAT channels here
    // (named formats touch FLOAT only). Fields set programmatically (neutral setup); the format change
    // is driven through the panel. numberFormat is baselined to 'Same as grid' first (select still shows '3 digits').
    await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      vw.setOptions({fieldsColumnNames: ['COMPUTED_H', 'AGE', 'SEX']});
      grok.shell.t.currentRowIdx = 0;
    });
    await setNumberFormatViaPanel(page, 'Same as grid');
    await page.waitForTimeout(400);
    const ageBefore = await cardFieldText(page, 'AGE');
    const sexBefore = await cardFieldText(page, 'SEX');
    // Tighten past not-null to NON-EMPTY: cardFieldText returns '' for a rendered-but-empty field,
    // degrading the before/after compare to '' === ''. Before values must carry real text.
    expect((ageBefore ?? '').length).toBeGreaterThan(0);
    expect((sexBefore ?? '').length).toBeGreaterThan(0);

    await setNumberFormatViaPanel(page, '3 digits after comma');
    // Wait for the FLOAT field to visibly change before reading the non-FLOAT ones, so "unchanged" is post-render.
    await expect.poll(() => cardFieldText(page, 'COMPUTED_H'), {timeout: 15_000}).toMatch(/^-?\d+\.\d{3}$/);
    expect(await cardFieldText(page, 'AGE')).toBe(ageBefore);
    expect(await cardFieldText(page, 'SEX')).toBe(sexBefore);
  });

  await cleanupShell(page);
  finishSpec();
});
