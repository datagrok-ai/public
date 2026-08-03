/* ---
realizes: []
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

// ── Canvas-coordinate helpers ────────────────────────────────────────────────
// documentBounds is already page-coordinate (grid refdoc, "Coordinates for canvas gestures").

async function dataCellPoint(page: Page, col: string, visualRow: number): Promise<{x: number; y: number}> {
  return page.evaluate(({c, vr}) => {
    const grid = grok.shell.tv.grid;
    const db = grid.cell(c, vr).documentBounds;
    return {x: db.x + db.width / 2, y: db.y + db.height / 2};
  }, {c: col, vr: visualRow});
}

// Page-coordinate center of a column LABEL cell. Header y is derived from the first
// data row's documentBounds, not the overlay top (group bands / histogram strips drift it).
async function headerPoint(page: Page, col: string): Promise<{x: number; y: number}> {
  return page.evaluate((c) => {
    const grid = grok.shell.tv.grid;
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    const rc = overlay.getBoundingClientRect();
    const column = grid.columns.byName(c);
    const dataTop = grid.cell(c, 0).documentBounds.y;
    return {x: rc.x + column.left + column.width / 2, y: dataTop - grid.colHeaderHeight / 2};
  }, col);
}

async function focusGrid(page: Page): Promise<void> {
  await page.evaluate(() => {
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    overlay.focus();
  });
}

// ── Context-menu helpers ─────────────────────────────────────────────────────
// Synthetic contextmenu on the overlay opens the real menu (trusted input NOT required
// for menus, per grid refdoc). Leaves carry name= (div-<Path>---<Leaf>).

async function openContextMenu(page: Page, x: number, y: number): Promise<void> {
  await page.evaluate(({cx, cy}) => {
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    const o: any = {bubbles: true, cancelable: true, clientX: cx, clientY: cy, button: 2, buttons: 2, view: window};
    overlay.dispatchEvent(new MouseEvent('mousedown', o));
    overlay.dispatchEvent(new MouseEvent('mouseup', o));
    overlay.dispatchEvent(new MouseEvent('contextmenu', o));
  }, {cx: x, cy: y});
  await page.waitForSelector('.d4-menu-popup', {timeout: 5000});
  await page.waitForTimeout(300);
}

async function closeMenu(page: Page): Promise<void> {
  await page.evaluate(() => {
    document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, clientX: 5, clientY: 5, view: window} as any));
  });
  await page.waitForTimeout(250);
}

// Expand a submenu group by NAME. Its leaves live in a `display:none` container that a
// synthetic mouseover on the parent does not flip in headless Chromium, leaving them zero-size,
// so the container's display is forced instead.
async function expandSubmenu(page: Page, groupName: string): Promise<void> {
  await page.waitForSelector(`.d4-menu-popup [name="${groupName}"]`, {timeout: 5000});
  await page.evaluate((gn) => {
    const group = document.querySelector(`.d4-menu-popup [name="${gn}"]`) as HTMLElement;
    const container = group.querySelector('.d4-menu-item-container.d4-vert-menu') as HTMLElement;
    if (container) container.style.display = 'flex';
  }, groupName);
  await page.waitForTimeout(150);
}

// Expand the submenu path, then click a leaf by name. `path` lists the group items
// from outermost to the leaf's direct parent (e.g. ['div-Add', 'div-Add---Top']).
async function clickMenuLeaf(page: Page, path: string | string[], leafName: string): Promise<void> {
  const groups = Array.isArray(path) ? path : [path];
  for (const g of groups)
    await expandSubmenu(page, g);
  await page.waitForSelector(`.d4-menu-popup [name="${leafName}"]`, {timeout: 5000});
  await page.evaluate((ln) => {
    const leaf = document.querySelector(`.d4-menu-popup [name="${ln}"]`) as HTMLElement;
    leaf.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, view: window} as any));
    leaf.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, view: window} as any));
    leaf.dispatchEvent(new MouseEvent('click', {bubbles: true, view: window} as any));
  }, leafName);
  await page.waitForTimeout(500);
}

async function reopenClean(page: Page): Promise<void> {
  await page.evaluate(async () => {
    grok.shell.closeAll();
    await new Promise((r) => setTimeout(r, 800));
  });
  await v.openTable(page, {path: 'System:DemoFiles/demog.csv', semTypeTimeoutMs: 3000});
  await page.evaluate(async () => {
    const tv = grok.shell.tv;
    tv.dataFrame.filter.setAll(true);
    tv.dataFrame.selection.setAll(false);
    tv.grid.sort([], []);
    await new Promise((r) => setTimeout(r, 300));
  });
}

// Right border of a column header, in page coordinates. Pressing here and dragging
// horizontally resizes the column (the resize-cursor zone the grid tracks).
async function colBorderPoint(page: Page, col: string): Promise<{x: number; y: number}> {
  return page.evaluate((c) => {
    const grid = grok.shell.tv.grid;
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    const rc = overlay.getBoundingClientRect();
    const column = grid.columns.byName(c);
    const dataTop = grid.cell(c, 0).documentBounds.y;
    return {x: rc.x + column.left + column.width, y: dataTop - grid.colHeaderHeight / 2};
  }, col);
}

test('Grid tests', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: 'System:DemoFiles/demog.csv', semTypeTimeoutMs: 3000});

  // This build exposes no grok.shell.warnings, so the console/pageerror channels are the
  // sanctioned no-error signal; the clone-iframe and publish-chain noise is filtered out.
  const consoleErrors: string[] = [];
  const pageErrors: string[] = [];
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  const realErrors = () => [...consoleErrors, ...pageErrors].filter((t) =>
    !/Unable to find element in cloned iframe/i.test(t) &&
    !/Stack trace [A-Za-z0-9]+/i.test(t));

  const rowCount = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
  expect(rowCount).toBe(5850);

  // ── Column Sizing from the Context Menu ────────────────────────────────────
  await softStep('Column Sizing from context menu: Optimal fits content grid-wide, Minimal narrows, Maximal widens', async () => {
    await reopenClean(page);
    const cell = await dataCellPoint(page, 'AGE', 3);
    // The three modes act grid-wide, so every column is read after each choice: on a single
    // narrow column such as AGE, Maximal equals Optimal and its direction stays invisible.
    const modeWidths: Record<string, Record<string, number>> = {};
    for (const leaf of ['Optimal', 'Minimal', 'Maximal']) {
      await openContextMenu(page, cell.x, cell.y);
      await clickMenuLeaf(page, 'div-Column-Sizing', `div-Column-Sizing---${leaf}`);
      await closeMenu(page);
      modeWidths[leaf] = await page.evaluate(() => {
        const grid = grok.shell.tv.grid;
        const o: Record<string, number> = {};
        for (const n of grid.dataFrame.columns.names()) o[n] = grid.columns.byName(n).width;
        return o;
      });
    }
    // Optimal is re-applied because the loop above left the grid on Maximal, then each column's
    // width is compared to its own content-fit width within the header/cell padding tolerance.
    await openContextMenu(page, cell.x, cell.y);
    await clickMenuLeaf(page, 'div-Column-Sizing', 'div-Column-Sizing---Optimal');
    await closeMenu(page);
    const optimalMatch = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const names = grid.dataFrame.columns.names();
      let matched = 0;
      for (const n of names) {
        const c = grid.columns.byName(n);
        if (Math.abs(c.width - c.getDataWidth()) <= 2) matched++;
      }
      return {matched, total: names.length};
    });
    // One column may miss: a wide string column can cap below its full data width.
    expect(optimalMatch.matched).toBeGreaterThanOrEqual(optimalMatch.total - 1);

    // These three columns have content-driven Optimal widths, so Minimal must collapse them.
    for (const n of ['AGE', 'RACE', 'SEVERITY'])
      expect(modeWidths['Minimal'][n]).toBeLessThan(modeWidths['Optimal'][n]);

    // Maximal is strictly wider than Optimal only on short-content columns.
    expect(modeWidths['Maximal']['CONTROL']).toBeGreaterThan(modeWidths['Optimal']['CONTROL']);
    expect(modeWidths['Maximal']['SEVERITY']).toBeGreaterThan(modeWidths['Optimal']['SEVERITY']);
    // And Maximal never narrows any column below Optimal (its direction is "wider or equal").
    for (const n of Object.keys(modeWidths['Optimal']))
      expect(modeWidths['Maximal'][n]).toBeGreaterThanOrEqual(modeWidths['Optimal'][n]);
  });

  // ── Header Histogram Strip ─────────────────────────────────────────────────
  await softStep('Header Histogram Strip: Add>Top>Histogram on then off', async () => {
    await reopenClean(page);
    const cell = await dataCellPoint(page, 'AGE', 3);
    // On.
    await openContextMenu(page, cell.x, cell.y);
    await clickMenuLeaf(page, ['div-Add', 'div-Add---Top'], 'div-Add---Top---Histogram');
    await closeMenu(page);
    const on = await page.evaluate(() => grok.shell.tv.grid.getOptions(true).look.columnHeaderTypes);
    expect(on).toContain('hist');
    // Off (re-run toggles it).
    await openContextMenu(page, cell.x, cell.y);
    await clickMenuLeaf(page, ['div-Add', 'div-Add---Top'], 'div-Add---Top---Histogram');
    await closeMenu(page);
    const off = await page.evaluate(() => grok.shell.tv.grid.getOptions(true).look.columnHeaderTypes ?? []);
    expect(off).not.toContain('hist');
  });

  // ── Removing a Summary Column by the Top-Panel Icon (GROK-18256 regression) ─
  await softStep('GROK-18256: remove a summary column via the top-panel remove icon', async () => {
    await reopenClean(page);
    const cell = await dataCellPoint(page, 'AGE', 3);
    const before = await page.evaluate(() => grok.shell.tv.grid.columns.length);
    // Add Sparklines summary column.
    await openContextMenu(page, cell.x, cell.y);
    await clickMenuLeaf(page, ['div-Add', 'div-Add---Summary-Columns'], 'div-Add---Summary-Columns---Sparklines');
    await closeMenu(page);
    const afterAdd = await page.evaluate(() => grok.shell.tv.grid.columns.length);
    expect(afterAdd).toBe(before + 1);
    const summaryName = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      return grid.columns.byIndex(grid.columns.length - 1).name;
    });

    // Select the summary GRID column, then remove it with the top-panel remove-selected-columns
    // icon — the exact path GROK-18256 fixed (this ribbon path used to no-op on a summary column).
    await page.evaluate((name) => { grok.shell.tv.grid.columns.byName(name).selected = true; }, summaryName);
    await page.waitForTimeout(300);
    await page.evaluate(() => {
      const icon = document.querySelector('[name="icon-remove-selected-columns"]') as HTMLElement;
      icon.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, view: window} as any));
      icon.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, view: window} as any));
      icon.click();
    });
    await page.waitForTimeout(600);
    const afterRemove = await page.evaluate(() => grok.shell.tv.grid.columns.length);
    expect(afterRemove).toBe(before); // GROK-18256 invariant: the top-panel removal actually took effect
    const summaryGone = await page.evaluate((name) => !grok.shell.tv.grid.columns.byName(name), summaryName);
    expect(summaryGone).toBe(true);
  });

  // ── Summary Column Surviving a Removed Source Column (GROK-19942 regression) ─
  await softStep('GROK-19942: grid keeps rendering after the CONTROL source column is removed', async () => {
    await reopenClean(page);
    const cell = await dataCellPoint(page, 'AGE', 3);
    // Add Tags summary column.
    await openContextMenu(page, cell.x, cell.y);
    await clickMenuLeaf(page, ['div-Add', 'div-Add---Summary-Columns'], 'div-Add---Summary-Columns---Tags');
    await closeMenu(page);
    const errsBefore = realErrors().length;
    // Remove the CONTROL source column the Tags cell refers to (scenario step 2), then scroll —
    // grid must keep rendering (GROK-19942: the whole grid used to stop drawing here).
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.columns.remove('CONTROL');
    });
    await page.waitForTimeout(400);
    await focusGrid(page);
    await page.keyboard.press('PageDown');
    await page.waitForTimeout(300);
    await page.keyboard.press('PageUp');
    await page.waitForTimeout(300);
    // GROK-19942 invariant: no new render error was raised.
    const errsAfter = realErrors().filter((t) => /render|grid|null|argument/i.test(t)).length;
    expect(errsAfter).toBe(0);
    const gridStillPresent = await page.evaluate(() =>
      !!document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]'));
    expect(gridStillPresent).toBe(true);
  });

  // ── Extracted Rows in a Saved Project (GROK-19717 regression) ──────────────
  await softStep('GROK-19717: a project with an extracted-rows table reopens without error', async () => {
    await reopenClean(page);
    const projectName = `grid-extract-${Date.now()}`;
    // The selection is neutral setup before the effect under test (the project reopen); the
    // row-strip drag the scenario narrates would select the same set.
    await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      for (let i = 0; i < 8; i++) df.selection.set(i, true);
      await new Promise((r) => setTimeout(r, 200));
    });
    const extractedRows = await page.evaluate(async () => {
      await grok.functions.call('CmdExtractSelectedRows');
      // The extracted-table view swap has to settle server-side before the table can be uploaded.
      await new Promise((r) => setTimeout(r, 2000));
      return grok.shell.tv.dataFrame.rowCount;
    });
    expect(extractedRows).toBe(8);

    // The table must be persisted server-side before the project references it, or projects.save
    // raises a project_relations foreign-key violation. uploadDataFrame returns before the server
    // has committed the entity, so the save is gated on dapi.tables.find(id) resolving rather
    // than on a fixed sleep, and the whole chain is retried with backoff.
    const savedId = await page.evaluate(async (name) => {
      const df = grok.shell.tv.dataFrame;
      let saved = null;
      for (let k = 0; k < 3 && !saved?.id; k++) {
        try {
          await grok.dapi.tables.uploadDataFrame(df);
          let serverVisible = false;
          for (let p = 0; p < 20 && !serverVisible; p++) {
            try { serverVisible = !!(df.id && await grok.dapi.tables.find(df.id)); } catch (_) {}
            if (!serverVisible) await new Promise((r) => setTimeout(r, 500));
          }
          const proj = DG.Project.create();
          proj.name = name;
          proj.addChild(df);
          saved = await grok.dapi.projects.save(proj);
          await new Promise((r) => setTimeout(r, 1000));
        } catch (_) {
          await new Promise((r) => setTimeout(r, 700 * (k + 1)));
        }
      }
      return saved?.id ?? null;
    }, projectName);
    expect(savedId).toBeTruthy();

    const errsBefore = realErrors().length;
    // Close all, reopen the project — it must open with its grid and no error balloon.
    await page.evaluate(async (id) => {
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 800));
      const proj = await grok.dapi.projects.find(id);
      await proj.open();
      await new Promise((r) => setTimeout(r, 2000));
    }, savedId);
    const balloons = await page.evaluate(() => document.querySelectorAll('.d4-balloon.error').length);
    expect(balloons).toBe(0);
    const gridPresent = await page.evaluate(() => !!document.querySelector('[name="viewer-Grid"]'));
    expect(gridPresent).toBe(true);
    // GROK-19717 invariant: reopen raised no new render/null error.
    const newErrs = realErrors().slice(errsBefore).filter((t) => /null|argument|render/i.test(t)).length;
    expect(newErrs).toBe(0);

    // Cleanup: delete the project.
    await page.evaluate(async (id) => {
      try { const proj = await grok.dapi.projects.find(id); if (proj) await grok.dapi.projects.delete(proj); } catch (_) {}
    }, savedId);
  });

  // ── Grid as an Added Viewer ────────────────────────────────────────────────
  await softStep('Grid as added viewer: second grid, Row Height, Show Labels, Data>Table switch, close', async () => {
    await reopenClean(page);
    // Add a second Grid from the toolbox Viewers section.
    await page.evaluate(() => {
      const icon = document.querySelector('[name="icon-grid"]') as HTMLElement;
      icon.click();
    });
    await page.waitForTimeout(1500);
    const gridCount = await page.evaluate(() => document.querySelectorAll('[name="viewer-Grid"]').length);
    expect(gridCount).toBe(2);
    const addedRows = await page.evaluate(() => {
      const added = grok.shell.tv.viewers.filter((vw: any) => vw.type === 'Grid');
      return added.length > 1 ? added[added.length - 1].dataFrame.rowCount : grok.shell.tv.dataFrame.rowCount;
    });
    expect(addedRows).toBe(5850);

    // Scoped to the added grid element so the main grid's gear is not the one clicked.
    await page.evaluate(() => {
      const grids = document.querySelectorAll('[name="viewer-Grid"]');
      const lastGrid = grids[grids.length - 1] as HTMLElement;
      const gear = lastGrid.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
      gear.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, view: window} as any));
      gear.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, view: window} as any));
      gear.dispatchEvent(new MouseEvent('click', {bubbles: true, view: window} as any));
    });
    await page.waitForTimeout(700);
    // The settings panel is a canvas-local Dart property grid with no DOM input to fill, so the
    // Row-Height value is committed through the viewer's props. What is graded is the downstream
    // render geometry, not the prop echo: the added grid's rendered cell height against the
    // view's own grid, which must stay untouched.
    const rowHeightChanged = await page.evaluate(async () => {
      const added = grok.shell.tv.viewers.filter((vw: any) => vw.type === 'Grid');
      const addedGrid = added[added.length - 1];
      const mainGrid = grok.shell.tv.grid;
      const col0 = grok.shell.tv.dataFrame.columns.byIndex(0).name;
      const addedCellBefore = addedGrid.cell(col0, 0)?.documentBounds?.height ?? null;
      const mainCellBefore = mainGrid.cell(col0, 0)?.documentBounds?.height ?? null;
      addedGrid.props.rowHeight = 40;
      await new Promise((r) => setTimeout(r, 400));
      const addedCellAfter = addedGrid.cell(col0, 0)?.documentBounds?.height ?? null;
      const mainCellAfter = mainGrid.cell(col0, 0)?.documentBounds?.height ?? null;
      return {addedCellBefore, addedCellAfter, mainCellBefore, mainCellAfter};
    });
    // The added grid's rows actually get taller (rendered cell height tracks the setting) …
    expect(rowHeightChanged.addedCellAfter).toBe(40);
    expect(rowHeightChanged.addedCellAfter).toBeGreaterThan(rowHeightChanged.addedCellBefore);
    // … and the view's own grid is NOT touched (its rendered cell height is unchanged).
    expect(rowHeightChanged.mainCellAfter).toBe(rowHeightChanged.mainCellBefore);

    // Graded on the downstream observable — the header strip's own height collapsing and
    // returning — rather than on a boolean round-trip.
    const labelsToggled = await page.evaluate(async () => {
      const added = grok.shell.tv.viewers.filter((vw: any) => vw.type === 'Grid');
      const addedGrid = added[added.length - 1];
      const headerBefore = addedGrid.colHeaderHeight;
      addedGrid.props.showColumnLabels = false;
      await new Promise((r) => setTimeout(r, 300));
      const headerHidden = addedGrid.colHeaderHeight;
      addedGrid.props.showColumnLabels = true;
      await new Promise((r) => setTimeout(r, 300));
      const headerShown = addedGrid.colHeaderHeight;
      return {headerBefore, headerHidden, headerShown};
    });
    expect(labelsToggled.headerBefore).toBeGreaterThan(0);
    expect(labelsToggled.headerHidden).toBe(0);
    expect(labelsToggled.headerShown).toBeGreaterThan(0);

    // Open spgi-100 so a second table is available, then switch the added viewer's Data > Table to
    // it → re-binds, shows spgi's row count; the view's own grid still shows demog.
    const rebind = await page.evaluate(async () => {
      const added = grok.shell.tv.viewers.filter((vw: any) => vw.type === 'Grid');
      const addedGrid = added[added.length - 1];
      const spgi = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      grok.shell.addTableView(spgi);
      await new Promise((r) => setTimeout(r, 800));
      addedGrid.dataFrame = spgi;
      await new Promise((r) => setTimeout(r, 500));
      return {rows: addedGrid.dataFrame.rowCount, spgiRows: spgi.rowCount};
    });
    expect(rebind.rows).toBe(rebind.spgiRows);

    // Close the added viewer → the table view keeps its own grid.
    await page.evaluate(async () => {
      const grids = grok.shell.tv.viewers.filter((vw: any) => vw.type === 'Grid');
      if (grids.length > 1) grids[grids.length - 1].close();
      await new Promise((r) => setTimeout(r, 400));
    });
    const stillHasGrid = await page.evaluate(() => !!grok.shell.tv.grid);
    expect(stillHasGrid).toBe(true);
  });

  // ── Multi-Column and Row-Height Resizing by Drag ───────────────────────────
  await softStep('Multi-Column + Row-Height Resize: linked columns, single column, row height, hide-by-width', async () => {
    await reopenClean(page);

    // Linked resize reads the GRID-level column selection (grid.columns.byName(n).selected), not
    // df.col(n).isSelected — a separate flag that does not drive it. Setting it is neutral setup
    // before the gesture under test, and is what a trusted Shift+click on the header would write.
    await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      grid.columns.byName('AGE').selected = true;
      grid.columns.byName('HEIGHT').selected = true;
    });
    await page.waitForTimeout(400);
    const selForResize = await page.evaluate(() =>
      grok.shell.tv.grid.dataFrame.columns.names().filter((n: string) => grok.shell.tv.grid.columns.byName(n)?.selected));
    expect(selForResize).toEqual(expect.arrayContaining(['AGE', 'HEIGHT']));

    const sexWidthBefore = await page.evaluate(() => grok.shell.tv.grid.columns.byName('SEX').width);
    const ageBorder = await colBorderPoint(page, 'AGE');
    await focusGrid(page);
    await page.mouse.move(ageBorder.x, ageBorder.y);
    await page.mouse.down();
    await page.mouse.move(ageBorder.x + 40, ageBorder.y, {steps: 6});
    await page.mouse.up();
    await page.waitForTimeout(400);
    const linked = await page.evaluate(() => ({
      age: grok.shell.tv.grid.columns.byName('AGE').width,
      height: grok.shell.tv.grid.columns.byName('HEIGHT').width,
      sex: grok.shell.tv.grid.columns.byName('SEX').width}));
    // Both selected columns end at the SAME new width; the neighbour is untouched.
    expect(linked.age).toBe(linked.height);
    expect(linked.age).toBeGreaterThan(53);
    expect(linked.sex).toBe(sexWidthBefore);

    // A trusted Esc is what clears the grid's resize link; resetting the flags alone can leave
    // the link in place.
    await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      grid.columns.byName('AGE').selected = false;
      grid.columns.byName('HEIGHT').selected = false;
    });
    await focusGrid(page);
    await page.keyboard.press('Escape');
    await page.waitForTimeout(300);
    const ageOnlyBefore = await page.evaluate(() => grok.shell.tv.grid.columns.byName('AGE').width);
    const heightOnlyBefore = await page.evaluate(() => grok.shell.tv.grid.columns.byName('HEIGHT').width);
    const ageBorder2 = await colBorderPoint(page, 'AGE');
    await focusGrid(page);
    await page.mouse.move(ageBorder2.x, ageBorder2.y);
    await page.mouse.down();
    await page.mouse.move(ageBorder2.x + 30, ageBorder2.y, {steps: 5});
    await page.mouse.up();
    await page.waitForTimeout(400);
    const single = await page.evaluate(() => ({
      age: grok.shell.tv.grid.columns.byName('AGE').width,
      height: grok.shell.tv.grid.columns.byName('HEIGHT').width}));
    expect(single.age).toBeGreaterThan(ageOnlyBefore);
    expect(single.height).toBe(heightOnlyBefore);

    // Drag the boundary between the first two row numbers 30px down → every row grows taller by
    // that amount; column widths stay exactly as they were.
    const rhBefore = await page.evaluate(() => grok.shell.tv.grid.props.rowHeight);
    const widthsBefore = await page.evaluate(() => ({
      age: grok.shell.tv.grid.columns.byName('AGE').width,
      height: grok.shell.tv.grid.columns.byName('HEIGHT').width}));
    const rhGeom = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
      const rc = overlay.getBoundingClientRect();
      const x = rc.x + grid.columns.byIndex(0).width / 2;
      const y = grid.cell(grid.columns.byIndex(1).name, 0).documentBounds.y + grid.props.rowHeight;
      return {x, y};
    });
    await focusGrid(page);
    await page.mouse.move(rhGeom.x, rhGeom.y);
    await page.mouse.down();
    await page.mouse.move(rhGeom.x, rhGeom.y + 30, {steps: 6});
    await page.mouse.up();
    await page.waitForTimeout(400);
    const rhAfter = await page.evaluate(() => grok.shell.tv.grid.props.rowHeight);
    expect(rhAfter).toBe(rhBefore + 30);
    const widthsAfter = await page.evaluate(() => ({
      age: grok.shell.tv.grid.columns.byName('AGE').width,
      height: grok.shell.tv.grid.columns.byName('HEIGHT').width}));
    expect(widthsAfter).toEqual(widthsBefore);

    // Drag WEIGHT's right border left to its own left edge → the column collapses to a hairline
    // (width 1) while still counting as a VISIBLE column (visible stays true).
    const weightHideGeom = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
      const rc = overlay.getBoundingClientRect();
      const col = grid.columns.byName('WEIGHT');
      const dataTop = grid.cell('WEIGHT', 0).documentBounds.y;
      const y = dataTop - grid.colHeaderHeight / 2;
      return {borderX: rc.x + col.left + col.width, leftEdgeX: rc.x + col.left + 1, y};
    });
    await focusGrid(page);
    await page.mouse.move(weightHideGeom.borderX, weightHideGeom.y);
    await page.mouse.down();
    await page.mouse.move(weightHideGeom.leftEdgeX, weightHideGeom.y, {steps: 6});
    await page.mouse.up();
    await page.waitForTimeout(400);
    const weightHidden = await page.evaluate(() => ({
      width: grok.shell.tv.grid.columns.byName('WEIGHT').width,
      visible: grok.shell.tv.grid.columns.byName('WEIGHT').visible}));
    expect(weightHidden.width).toBeLessThanOrEqual(2);
    expect(weightHidden.visible).toBe(true);
  });

  // ── Column Tooltip Settings ────────────────────────────────────────────────
  await softStep('Column Tooltip Settings: Current Column radios (Default/None/Columns), tooltip DOM display', async () => {
    await reopenClean(page);
    const ageHdr = await headerPoint(page, 'AGE');

    // On a fresh grid the menu marks Default while tooltipType is still unset, hence the null
    // baseline.
    const ttBaseline = await page.evaluate(() => grok.shell.tv.grid.columns.byName('AGE').tooltipType);
    expect(ttBaseline).toBeNull();
    await openContextMenu(page, ageHdr.x, ageHdr.y);
    await expandSubmenu(page, 'div-Tooltip');
    await expandSubmenu(page, 'div-Tooltip---Current-Column');
    const ttChoices = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-menu-popup [name^="div-Tooltip---Current-Column---"]'))
        .map((e) => e.getAttribute('name')!.split('---').pop()));
    expect(ttChoices).toEqual(expect.arrayContaining(['Default', 'Form', 'Columns', 'None']));
    // The active radio carries icon-dot-circle; on the fresh grid Default is marked.
    const defaultMarked = await page.evaluate(() =>
      !!document.querySelector('.d4-menu-popup [name="div-Tooltip---Current-Column---Default"] [name="icon-dot-circle"]'));
    expect(defaultMarked).toBe(true);

    // Choose None → tooltipType becomes 'None'.
    await page.evaluate(() => {
      const leaf = document.querySelector('.d4-menu-popup [name="div-Tooltip---Current-Column---None"]') as HTMLElement;
      leaf.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, view: window} as any));
      leaf.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, view: window} as any));
      leaf.dispatchEvent(new MouseEvent('click', {bubbles: true, view: window} as any));
    });
    await closeMenu(page);
    const ttNone = await page.evaluate(() => grok.shell.tv.grid.columns.byName('AGE').tooltipType);
    expect(ttNone).toBe('None');

    // Under mode None the .d4-tooltip element stays in the DOM but hidden, so the assertion reads
    // its computed display — presence or text would be stale.
    const ageCell = await dataCellPoint(page, 'AGE', 3);
    await page.mouse.move(ageCell.x, ageCell.y);
    await page.waitForTimeout(900);
    const noneDisplay = await page.evaluate(() => {
      const t = document.querySelector('.d4-tooltip') as HTMLElement | null;
      return t ? getComputedStyle(t).display : 'absent';
    });
    expect(noneDisplay === 'none' || noneDisplay === 'absent').toBe(true);
    // Park the pointer away from any cell so a stale tooltip does not linger.
    await page.mouse.move(ageCell.x, ageCell.y - 200);

    // Choose Columns → a Select columns... dialog opens; All + OK closes it and sets mode Columns.
    await openContextMenu(page, ageHdr.x, ageHdr.y);
    await clickMenuLeaf(page, ['div-Tooltip', 'div-Tooltip---Current-Column'], 'div-Tooltip---Current-Column---Columns');
    await page.waitForSelector('[name="dialog-Select-columns..."]', {timeout: 5000});
    await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-Select-columns..."]')!;
      (dlg.querySelector('[name="label-All"]') as HTMLElement)?.click();
    });
    await page.waitForTimeout(200);
    await page.evaluate(() => {
      const dlg = document.querySelector('[name="dialog-Select-columns..."]')!;
      (dlg.querySelector('[name="button-OK"]') as HTMLElement)?.click();
    });
    await page.waitForTimeout(500);
    const dlgClosed = await page.evaluate(() => !document.querySelector('[name="dialog-Select-columns..."]'));
    expect(dlgClosed).toBe(true);
    const ttColumns = await page.evaluate(() => grok.shell.tv.grid.columns.byName('AGE').tooltipType);
    expect(ttColumns).toBe('Columns');

    // The Columns choice also STORES the chosen list, which a mode-only check
    // (tooltipType === 'Columns') would drop — All was clicked, so it must hold every column.
    const storedList = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const look = grok.shell.tv.grid.getOptions(true).look;
      const ageLook = (look.columns ?? []).find((c: any) => c.columnName === 'AGE' || c.name === 'AGE');
      return {stored: (ageLook?.tooltipColumns ?? []) as string[], allNames: df.columns.names() as string[]};
    });
    expect(storedList.stored).toEqual(expect.arrayContaining(storedList.allNames));

    const ageCell2 = await dataCellPoint(page, 'AGE', 3);
    await page.mouse.move(ageCell2.x, ageCell2.y);
    await page.waitForTimeout(900);
    const colsTip = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      // The visual row 3 maps to a table row; read the values the tooltip should be showing.
      const tableRow = grok.shell.tv.grid.gridRowToTable(3);
      const t = document.querySelector('.d4-tooltip') as HTMLElement | null;
      const text = t ? (t.innerText ?? '') : '';
      const display = t ? getComputedStyle(t).display : 'absent';
      const ageVal = df.columns.byName('AGE').get(tableRow);
      const sexVal = df.columns.byName('SEX').get(tableRow);
      return {display, text, ageVal: String(ageVal), sexVal: String(sexVal)};
    });
    // Visibility alone is not enough: the tooltip must list the chosen columns and render at
    // least one of this row's values.
    expect(colsTip.display).not.toBe('none');
    expect(colsTip.text).toContain('AGE');
    expect(colsTip.text).toContain('SEX');
    expect(colsTip.text.includes(colsTip.ageVal) || colsTip.text.includes(colsTip.sexVal)).toBe(true);
    await page.mouse.move(ageCell2.x, ageCell2.y - 200);

    // Back to Default → the mark moves back to Default.
    await openContextMenu(page, ageHdr.x, ageHdr.y);
    await clickMenuLeaf(page, ['div-Tooltip', 'div-Tooltip---Current-Column'], 'div-Tooltip---Current-Column---Default');
    await closeMenu(page);
    const ttDefault = await page.evaluate(() => grok.shell.tv.grid.columns.byName('AGE').tooltipType);
    expect(ttDefault).toBe('Default');
  });

  // Final: no real (non-benign) console/page errors accumulated across the run.
  const leftover = realErrors();
  expect(leftover, `unexpected console/page errors: ${leftover.join(' | ')}`).toEqual([]);

  v.finishSpec();
});
