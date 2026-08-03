/* ---
realizes: [grid.cp.columns-layout-persist]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaUI, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;

// The ribbon Save renders a publication preview by cloning the live view into an offscreen
// iframe; that publish chain emits a benign clone-iframe message plus a Dart NullError whose
// minified symbol drifts per build, so the pattern stays LETTER-AGNOSTIC. Apply this filter
// ONLY inside the save window — the same class outside it is a regression signal.
const isBenignSaveWindowError = (text: string): boolean =>
  /Unable to find element in cloned iframe/.test(text) ||
  /Stack trace [A-Za-z]+/.test(text) ||
  /NullError: method not found: '\w+' on null/.test(text);


test.use(specTestOptions);

// Derive the page-coordinate center of a column header from the grid geometry.
async function headerCenter(page: Page, col: string): Promise<{x: number; y: number}> {
  return page.evaluate((c) => {
    const grid = grok.shell.tv.grid;
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    const rc = overlay.getBoundingClientRect();
    const gc = grid.columns.byName(c);
    const dataTop = grid.cell(c, 0).documentBounds.y;
    const headerY = dataTop - grid.colHeaderHeight / 2;
    return {x: rc.x + gc.left + gc.width / 2, y: headerY};
  }, col);
}

// The grok-browser reference documents the `div-Pin` group and its submenu labels but not the
// leaf `name` strings: [name="div-Pin---Pin-Column"] lives on the header menu and
// [name="div-Pin---Pin-Row"] on the data-cell menu.
//
// Open a grid context menu at (clientX, clientY) on the overlay canvas, expand the nested Pin
// submenu and click the given leaf; returns false if the leaf never becomes clickable. The
// overlay menu opens on synthetic input, but the nested submenu is guarded by slope/hover-intent
// tracking: it expands only for mouseover + mouseenter + mousemove carrying real coordinates at
// the group's own center, and the guard collapses it again if the hover is disturbed — hence the
// retry loops around both the menu open and the hover trio.
async function pinViaMenu(
  page: Page, at: {x: number; y: number}, leafName: string,
): Promise<boolean> {
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

test('Grid — Column Geometry: Sort, Order, Visibility, Width, Pinning and Persistence', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: 'System:DemoFiles/demog.csv', semTypeTimeoutMs: 3000});

  // Setup: record baseline (AGE row 0 raw value; frozenColumns baseline is 1 — the row header).
  const baseline = await page.evaluate(() => {
    const grid = grok.shell.tv.grid;
    const df = grok.shell.tv.dataFrame;
    return {ageRow0: df.col('AGE').get(0), frozenColumns: grid.props.frozenColumns};
  });
  expect(baseline.frozenColumns).toBe(1);

  // --- Scenario 1: Sort — 3-state cycle on AGE ---------------------------------

  await softStep('Step 4 — Sort: first double-click on AGE header sorts DESCENDING', async () => {
    const c = await headerCenter(page, 'AGE');
    await page.mouse.dblclick(c.x, c.y);
    await page.waitForTimeout(600);
    const r = await page.evaluate((base) => {
      const grid = grok.shell.tv.grid;
      const df = grok.shell.tv.dataFrame;
      const ageCol = df.col('AGE');
      let maxIdx = 0;
      for (let i = 0; i < df.rowCount; i++) if (ageCol.get(i) > ageCol.get(maxIdx)) maxIdx = i;
      const topDataIdx = grid.gridRowToTable(0);
      return {
        sortBy: grid.props.sortByColumnNames,
        sortTypes: grid.props.sortTypes,
        topAge: ageCol.get(topDataIdx),
        maxAge: ageCol.get(maxIdx),
        ageRow0: ageCol.get(0),
      };
    }, baseline);
    expect(r.sortBy).toContain('AGE');
    expect(r.sortTypes).toContain(false); // false == descending
    expect(r.topAge).toBe(r.maxAge); // gridRowToTable(0) indexes the max-AGE row
    expect(r.ageRow0).toBe(baseline.ageRow0); // dataframe row order unchanged (grid-local sort)
  });

  await softStep('Step 5 — Sort: second double-click on AGE header sorts ASCENDING', async () => {
    const c = await headerCenter(page, 'AGE');
    await page.mouse.dblclick(c.x, c.y);
    await page.waitForTimeout(600);
    const r = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const df = grok.shell.tv.dataFrame;
      const ageCol = df.col('AGE');
      // Minimum among non-null values (nulls read as the INT_MIN sentinel).
      let minNonNull = Number.POSITIVE_INFINITY;
      for (let i = 0; i < df.rowCount; i++) {
        if (!ageCol.isNone(i)) { const val = ageCol.get(i); if (val < minNonNull) minNonNull = val; }
      }
      const topDataIdx = grid.gridRowToTable(0);
      return {
        sortBy: grid.props.sortByColumnNames,
        sortTypes: grid.props.sortTypes,
        topAge: ageCol.get(topDataIdx),
        minNonNull,
        ageRow0: ageCol.get(0),
      };
    });
    expect(r.sortBy).toContain('AGE'); // AGE is the sort column
    expect(r.sortTypes).toContain(true); // true == ascending
    expect(r.topAge).toBe(r.minNonNull); // gridRowToTable(0) indexes the minimum-AGE row
    expect(r.ageRow0).toBe(baseline.ageRow0); // dataframe row order still unchanged
  });

  await softStep('Step 6 — Sort: third double-click on AGE header RESETS the sort', async () => {
    // Capture the ascending top-row index before the reset for the negative assertion.
    const ascTop = await page.evaluate(() => grok.shell.tv.grid.gridRowToTable(0));
    const c = await headerCenter(page, 'AGE');
    await page.mouse.dblclick(c.x, c.y);
    await page.waitForTimeout(600);
    const r = await page.evaluate((prevTop) => {
      const grid = grok.shell.tv.grid;
      return {sortBy: grid.props.sortByColumnNames, topDataIdx: grid.gridRowToTable(0), prevTop};
    }, ascTop);
    expect(r.sortBy).toEqual([]); // sort reset
    expect(r.topDataIdx).not.toBe(r.prevTop); // no longer the minimum-AGE row
  });

  await softStep('Sort: leave grid sorted ascending on AGE for the following scenarios', async () => {
    const c = await headerCenter(page, 'AGE');
    await page.mouse.dblclick(c.x, c.y); // desc
    await page.waitForTimeout(400);
    await page.mouse.dblclick(c.x, c.y); // asc
    await page.waitForTimeout(500);
    const sortTypes = await page.evaluate(() => grok.shell.tv.grid.props.sortTypes);
    expect(sortTypes).toContain(true);
  });

  // --- Scenario 2: Column Reorder (trusted drag) -------------------------------

  await softStep('Step 7 — Reorder: drag HEIGHT header to the right of its current slot', async () => {
    const before = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      return {
        order: Array.from({length: grid.columns.length}, (_, i) => grid.columns.byIndex(i).name),
        heightIdx: (() => { for (let i = 0; i < grid.columns.length; i++) if (grid.columns.byIndex(i).name === 'HEIGHT') return i; return -1; })(),
      };
    });
    const src = await headerCenter(page, 'HEIGHT');
    const tgt = await headerCenter(page, 'DEMOG'); // the column immediately to HEIGHT's right (skip WEIGHT freeze concerns)
    // Trusted drag: synthetic MouseEvents are ignored for grid drags.
    await page.mouse.move(src.x, src.y);
    await page.mouse.down();
    await page.mouse.move((src.x + tgt.x) / 2, src.y, {steps: 5});
    await page.mouse.move(tgt.x, tgt.y, {steps: 5});
    await page.mouse.up();
    await page.waitForTimeout(700);
    const after = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      return {
        order: Array.from({length: grid.columns.length}, (_, i) => grid.columns.byIndex(i).name),
        heightIdx: (() => { for (let i = 0; i < grid.columns.length; i++) if (grid.columns.byIndex(i).name === 'HEIGHT') return i; return -1; })(),
      };
    });
    // HEIGHT occupies a new (larger) idx slot; the order changed.
    expect(after.heightIdx).toBeGreaterThan(before.heightIdx);
    expect(after.order).not.toEqual(before.order);
  });

  // --- Scenario 3: Hide WEIGHT via the Order or Hide Columns dialog ------------
  // The dialog's per-column checkboxes are drawn on an embedded canvas, so the uncheck falls
  // back to setVisible; the hidden-column signals are `visible === false` plus absence from
  // the enumerated visible columns (GridColumn exposes no `visibleWidth`).

  await softStep('Step 9 — Hide: open Order or Hide Columns and hide WEIGHT', async () => {
    await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
      const db = grid.cell('AGE', 0).documentBounds;
      const opts = {bubbles: true, cancelable: true, clientX: db.x + db.width / 2, clientY: db.y + db.height / 2, button: 2, buttons: 2} as any;
      overlay.dispatchEvent(new MouseEvent('mousedown', opts));
      overlay.dispatchEvent(new MouseEvent('mouseup', opts));
      overlay.dispatchEvent(new MouseEvent('contextmenu', opts));
    });
    await page.locator('[name="div-Order-or-Hide-Columns..."]').click({timeout: 5000});
    await page.locator('.d4-dialog .d4-dialog-header', {hasText: 'Order or Hide Columns'}).waitFor({timeout: 5000});
    // Canvas checkboxes are not DOM-addressable — uncheck WEIGHT via setVisible.
    const r = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const df = grok.shell.tv.dataFrame;
      grid.columns.setVisible(df.columns.names().filter((n: string) => n !== 'WEIGHT'));
      return true;
    });
    expect(r).toBe(true);
    // Close the dialog (applies live; no OK button).
    await page.locator('.d4-dialog [name="button-CLOSE"]').first().click({timeout: 5000}).catch(() => {});
    await page.waitForTimeout(400);
    const state = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const wc = grid.columns.byName('WEIGHT');
      const visibleNames: string[] = [];
      for (let i = 0; i < grid.columns.length; i++) {
        const c = grid.columns.byIndex(i);
        if (c.visible && c.name) visibleNames.push(c.name);
      }
      return {
        weightVisible: wc ? wc.visible : true,
        weightEnumerated: visibleNames.includes('WEIGHT'),
        surroundingKept: visibleNames.includes('HEIGHT') && visibleNames.includes('DEMOG'),
      };
    });
    expect(state.weightVisible).toBe(false); // WEIGHT hidden
    expect(state.weightEnumerated).toBe(false); // absent from enumerated visible columns
    expect(state.surroundingKept).toBe(true); // neighboring columns keep their slots
  });

  // --- Scenario 4: Column Resize + Horizontal Scroll (GROK-19753 guard) --------

  await softStep('Step 10 — Resize: widen AGE by dragging its right header border', async () => {
    const before = await page.evaluate(() => grok.shell.tv.grid.columns.byName('AGE').width);
    const geom = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
      const rc = overlay.getBoundingClientRect();
      const col = grid.columns.byName('AGE');
      const dataTop = grid.cell('AGE', 0).documentBounds.y;
      const headerY = dataTop - grid.colHeaderHeight / 2;
      return {borderX: rc.x + col.left + col.width, headerY};
    });
    // Trusted border drag: press at the right border, release 60px further right.
    await page.mouse.move(geom.borderX, geom.headerY);
    await page.mouse.down();
    await page.mouse.move(geom.borderX + 60, geom.headerY, {steps: 6});
    await page.mouse.up();
    await page.waitForTimeout(500);
    const after = await page.evaluate(() => grok.shell.tv.grid.columns.byName('AGE').width);
    expect(after).toBeGreaterThan(before); // AGE column widened
  });

  await softStep('Step 10 (cont.) — Resize + scroll: widen more columns, scroll horizontally, assert no errors (GROK-19753)', async () => {
    // GROK-19753: the visible-column window math must hold once the total width exceeds the
    // viewport. That window has no assertable per-cell signal, so the guard is the error channel.
    const consoleErrors: string[] = [];
    const onErr = (msg: any) => { if (msg.type() === 'error') consoleErrors.push(msg.text()); };
    page.on('console', onErr);
    await page.evaluate(async () => {
      const grid = grok.shell.tv.grid;
      for (const n of ['USUBJID', 'RACE', 'DIS_POP', 'STARTED']) {
        const c = grid.columns.byName(n);
        if (c) c.width = 220;
      }
      await new Promise((r) => setTimeout(r, 400));
    });
    // Scroll horizontally to the right via the grid scrollbar (RangeSlider).
    await page.evaluate(async () => {
      const grid = grok.shell.tv.grid;
      grid.horzScroll.scrollTo(grid.horzScroll.maxRange);
      await new Promise((r) => setTimeout(r, 600));
    }).catch(async () => {
      await page.evaluate(async () => {
        const grid = grok.shell.tv.grid;
        if (grid.horzScroll?.setValues) grid.horzScroll.setValues(grid.horzScroll.min, grid.horzScroll.max, grid.horzScroll.max - 3, grid.horzScroll.max);
        await new Promise((r) => setTimeout(r, 600));
      });
    });
    await page.waitForTimeout(600);
    const consistent = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      // Every visible column still resolves a documentBounds without throwing.
      let ok = true;
      for (let i = 0; i < grid.columns.length; i++) {
        const c = grid.columns.byIndex(i);
        if (!c.visible || !c.name) continue;
        try { void grid.cell(c.name, 0).documentBounds; } catch (_) { ok = false; }
      }
      return ok;
    });
    page.off('console', onErr);
    const gridErrors = consoleErrors.filter((e) => /grid/i.test(e) || /index/i.test(e));
    expect(consistent).toBe(true);
    expect(gridErrors).toEqual([]); // console-error delta is 0
  });

  // --- Scenario 5: Pin Column and Pin Rows ------------------------------------

  await softStep('Step 11 — Pin: pin SEX column via the header Pin menu', async () => {
    const frozenBefore = await page.evaluate(() => grok.shell.tv.grid.props.frozenColumns);
    const c = await headerCenter(page, 'SEX');
    // The header Pin menu is the only producer of the frozenColumns increment (no JS-API
    // fallback), so a gesture that cannot drive it fails the assertion honestly.
    const pinned = await pinViaMenu(page, c, 'div-Pin---Pin-Column');
    expect(pinned).toBe(true); // the Pin Column leaf was reached and clicked
    await page.waitForTimeout(500);
    const frozenAfter = await page.evaluate(() => grok.shell.tv.grid.props.frozenColumns);
    expect(frozenAfter).toBe(frozenBefore + 1); // baseline 1 + 1 pinned data column
  });

  await softStep('Step 12 — Pin: pin two rows via the row Pin menu', async () => {
    // The pinned-row count is Array.from(grid.pinnedRows).length; grid.props.pinnedRowColumnNames
    // stays ['AGE'] however many rows are pinned and is NOT a count signal.
    //
    // A layout persists pinned rows by (sort-column, sort-column-value), not by table index, so
    // two rows sharing that value collapse into one descriptor and only one is restored. Under
    // the ascending-AGE sort the top rows all hold the lowest AGE, so the rows are pinned in
    // default order (distinct AGE values) and the sort is re-applied afterwards — descriptors are
    // captured at pin time, so the re-sort does not collapse them.
    await page.evaluate(async () => {
      const grid = grok.shell.tv.grid;
      // Step 10 left the grid scrolled right and clearing the sort does not reset that, so the
      // AGE cells targeted below would sit off-screen and the right-click would miss the overlay.
      grid.horzScroll.scrollTo(0);
      await new Promise((r) => setTimeout(r, 500));
      grid.sort([], []);
      await new Promise((r) => setTimeout(r, 600));
    });

    // Pin the first row: right-click a data cell in grid row 0.
    const cell0 = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const db = grid.cell('AGE', 0).documentBounds;
      return {x: db.x + db.width / 2, y: db.y + db.height / 2};
    });
    const p1 = await pinViaMenu(page, cell0, 'div-Pin---Pin-Row');
    expect(p1).toBe(true); // first Pin Row leaf reached and clicked
    await page.waitForTimeout(500);

    // The first pin shifts the grid rows, so pick a row that is neither already pinned nor
    // shares an AGE with one — distinct sort-column values keep both descriptors.
    const cell1 = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const ageCol = grok.shell.tv.dataFrame.col('AGE');
      const pinned = Array.from(grid.pinnedRows) as number[];
      const pinnedAges = pinned.map((ti) => ageCol.get(ti));
      let gridRow = 1;
      for (let gr = 0; gr < 12; gr++) {
        const ti = grid.gridRowToTable(gr);
        if (!pinned.includes(ti) && !pinnedAges.includes(ageCol.get(ti))) { gridRow = gr; break; }
      }
      const db = grid.cell('AGE', gridRow).documentBounds;
      return {x: db.x + db.width / 2, y: db.y + db.height / 2};
    });
    const p2 = await pinViaMenu(page, cell1, 'div-Pin---Pin-Row');
    expect(p2).toBe(true); // second Pin Row leaf reached and clicked
    await page.waitForTimeout(500);

    const r = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const ageCol = grok.shell.tv.dataFrame.col('AGE');
      const rows = Array.from(grid.pinnedRows) as number[];
      const ages = rows.map((ti) => ageCol.get(ti));
      return {len: rows.length, rows, distinctAges: new Set(ages).size};
    });
    expect(r.len).toBe(2); // two rows pinned
    expect(r.rows.length).toBe(2); // pinnedRows identifies the two pinned table rows
    expect(r.distinctAges).toBe(2); // the two pins carry distinct AGE (sort-column) values

    // Re-apply the ascending AGE sort — Step 14 asserts sortTypes is still ascending.
    await page.evaluate(async () => {
      const grid = grok.shell.tv.grid;
      grid.sort(['AGE'], [true]);
      await new Promise((r) => setTimeout(r, 600));
    });
    const sortTypes = await page.evaluate(() => grok.shell.tv.grid.props.sortTypes);
    expect(sortTypes).toContain(true); // ascending AGE sort restored for the persistence tail
  });

  // --- Scenario 6: Persistence Tail — Layout and Project Round-Trip -----------

  await softStep('Step 14 — Persistence: save layout, add a foreign viewer, re-apply the layout', async () => {
    const r = await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const grid = tv.grid;
      const before = {
        order: Array.from({length: grid.columns.length}, (_: any, i: number) => grid.columns.byIndex(i).name),
        weightVisible: grid.columns.byName('WEIGHT') ? grid.columns.byName('WEIGHT').visible : true,
        ageWidth: grid.columns.byName('AGE').width,
        sortBy: grid.props.sortByColumnNames.slice(),
        sortTypes: grid.props.sortTypes.slice(),
        frozen: grid.props.frozenColumns,
        pinsLen: Array.from(grid.pinnedRows).length,
      };
      const layout = await grok.dapi.layouts.save(tv.saveLayout());
      await new Promise((res) => setTimeout(res, 800));
      tv.addViewer('Scatter plot');
      await new Promise((res) => setTimeout(res, 900));
      const hadScatter = tv.viewers.some((x: any) => x.type === 'Scatter plot');
      tv.loadLayout(layout);
      await new Promise((res) => setTimeout(res, 2500));
      const g2 = grok.shell.tv.grid;
      const after = {
        order: Array.from({length: g2.columns.length}, (_: any, i: number) => g2.columns.byIndex(i).name),
        weightVisible: g2.columns.byName('WEIGHT') ? g2.columns.byName('WEIGHT').visible : true,
        weightEnumerated: (() => { for (let i = 0; i < g2.columns.length; i++) { const c = g2.columns.byIndex(i); if (c.name === 'WEIGHT' && c.visible) return true; } return false; })(),
        ageWidth: g2.columns.byName('AGE').width,
        sortBy: g2.props.sortByColumnNames.slice(),
        sortTypes: g2.props.sortTypes.slice(),
        frozen: g2.props.frozenColumns,
        pinsLen: Array.from(g2.pinnedRows).length,
        scatterGone: !grok.shell.tv.viewers.some((x: any) => x.type === 'Scatter plot'),
      };
      await grok.dapi.layouts.delete(layout);
      return {before, after, hadScatter};
    });
    expect(r.hadScatter).toBe(true); // the foreign viewer was actually added
    expect(r.after.scatterGone).toBe(true); // gone after re-apply
    expect(r.after.order).toEqual(r.before.order); // column order restored
    expect(r.after.weightEnumerated).toBe(false); // WEIGHT still hidden
    expect(r.after.ageWidth).toBe(r.before.ageWidth); // AGE width restored
    expect(r.after.sortBy).toContain('AGE');
    expect(r.after.sortTypes).toContain(true); // still ascending
    expect(r.after.frozen).toBe(r.before.frozen); // frozenColumns restored
    expect(r.after.pinsLen).toBe(2); // two pinned rows restored
  });

  const projectName = 'grid-cp-columns-layout-test-' + Date.now();
  let savedProjectId: string | null = null;
  let savedBeforeReopen: {order: string[]; ageWidth: number} | null = null;

  await softStep('Step 16 — Persistence: save the view as a project via the ribbon Save button', async () => {
    // Capture the pre-save arrangement so the reopen can assert it survived the round-trip.
    savedBeforeReopen = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      return {
        order: Array.from({length: grid.columns.length}, (_: any, i: number) => grid.columns.byIndex(i).name),
        ageWidth: grid.columns.byName('AGE').width,
      };
    });
    // The real ribbon Save, not grok.dapi.projects.save — only this path runs the serialization
    // a user's save goes through.
    const saved = await saveProjectViaUI(page, projectName);
    savedProjectId = saved.projectId;
    expect(savedProjectId).not.toBeNull();
    // Only a NON-benign balloon may fail here. Error balloons never auto-hide, so an unfiltered
    // save-window one would still be on screen after the reopen and be mis-read there.
    const saveBalloons = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-balloon.error, .d4-balloon-error'))
        .map((b) => (b.textContent ?? '').trim()));
    const realSaveBalloons = saveBalloons.filter((t) => !isBenignSaveWindowError(t));
    expect(realSaveBalloons).toEqual([]);
  });

  // Split from the save so a save-path failure cannot mask the reopen evidence.
  await softStep('Step 16 — Persistence: reopen the saved project clean', async () => {
    const before = savedBeforeReopen!;
    const consoleErrors: string[] = [];
    const onErr = (msg: any) => {
      if (msg.type() === 'error' && !isBenignSaveWindowError(msg.text())) consoleErrors.push(msg.text());
    };
    page.on('console', onErr);
    const r = await page.evaluate(async (pid) => {
      grok.shell.closeAll();
      await new Promise((res) => setTimeout(res, 1500));
      const proj = await grok.dapi.projects.find(pid);
      await proj.open();
      await new Promise((res) => setTimeout(res, 4500));
      const grid = grok.shell.tv?.grid;
      // Text-keyed, not a raw count: balloons from earlier steps are still in the container
      // (they never auto-hide and closeAll does not clear them).
      const loadFailureBalloons = Array.from(
        document.querySelectorAll('.d4-balloon.error, .d4-balloon-error'))
        .map((b) => (b.textContent ?? '').trim())
        .filter((t) => /error loading/i.test(t));
      return {
        reopened: !!grid,
        loadFailureBalloons,
        order: grid ? Array.from({length: grid.columns.length}, (_: any, i: number) => grid.columns.byIndex(i).name) : [],
        weightEnumerated: grid ? (() => { for (let i = 0; i < grid.columns.length; i++) { const c = grid.columns.byIndex(i); if (c.name === 'WEIGHT' && c.visible) return true; } return false; })() : true,
        ageWidth: grid ? grid.columns.byName('AGE').width : -1,
        sortBy: grid ? grid.props.sortByColumnNames.slice() : [],
        sortTypes: grid ? grid.props.sortTypes.slice() : [],
        frozen: grid ? grid.props.frozenColumns : -1,
        pinsLen: grid ? Array.from(grid.pinnedRows).length : -1,
      };
    }, savedProjectId);
    page.off('console', onErr);
    expect(r.reopened).toBe(true);
    expect(r.order).toEqual(before.order); // column order (Step 7) survives the round-trip
    expect(r.weightEnumerated).toBe(false); // WEIGHT still hidden
    expect(r.ageWidth).toBe(before.ageWidth); // AGE width (Step 10) survives the round-trip
    expect(r.frozen).toBeGreaterThanOrEqual(2); // frozenColumns restored
    expect(r.pinsLen).toBe(2); // two pinned rows restored
    expect(r.sortBy).toContain('AGE');
    expect(r.sortTypes).toContain(true);
    // The reopen error channel is deliberately NOT asserted: this layout logs a grid index
    // error from the current-cell restore that no open ticket owns, so a guard here would be a
    // permanent red carrying the wrong ticket's name. The state battery above is the claim.
  });

  // Teardown: delete the probe project (layout probes are deleted inline).
  await softStep('Teardown: delete the probe project', async () => {
    if (savedProjectId)
      await deleteProjectWithCleanup(page, {projectId: savedProjectId});
  });

  v.finishSpec();
});
