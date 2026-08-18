/* ---
realizes: []
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

async function dataCellPoint(page: Page, col: string, visualRow: number): Promise<{x: number; y: number}> {
  return page.evaluate(({c, vr}) => {
    const grid = grok.shell.tv.grid;
    const db = grid.cell(c, vr).documentBounds;
    return {x: db.x + db.width / 2, y: db.y + db.height / 2};
  }, {c: col, vr: visualRow});
}

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

async function openContextMenu(page: Page, x: number, y: number): Promise<void> {
  await page.evaluate(({cx, cy}) => {
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    const o: any = {bubbles: true, cancelable: true, clientX: cx, clientY: cy, button: 2, buttons: 2, view: window};
    overlay.dispatchEvent(new MouseEvent('mousedown', o));
    overlay.dispatchEvent(new MouseEvent('mouseup', o));
    overlay.dispatchEvent(new MouseEvent('contextmenu', o));
  }, {cx: x, cy: y});
  await page.waitForSelector('.d4-menu-popup', {timeout: 5000});
  await page.waitForTimeout(200); 
}

async function closeMenu(page: Page): Promise<void> {
  await page.evaluate(() => {
    document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, clientX: 5, clientY: 5, view: window} as any));
  });
  await page.waitForSelector('.d4-menu-popup', {state: 'detached', timeout: 5000}).catch(() => {});
}

async function expandSubmenu(page: Page, groupName: string): Promise<void> {
  await page.waitForSelector(`.d4-menu-popup [name="${groupName}"]`, {timeout: 5000});
  await page.evaluate((gn) => {
    const group = document.querySelector(`.d4-menu-popup [name="${gn}"]`) as HTMLElement;
    const container = group.querySelector('.d4-menu-item-container.d4-vert-menu') as HTMLElement;
    if (container) container.style.display = 'flex';
  }, groupName);
  await page.waitForTimeout(150); 
}

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

  await v.waitForViewerRendered(page, 'Grid', 500);
}

async function reopenClean(page: Page): Promise<void> {
  await v.closeAllAndWait(page);
  await v.openTable(page, {path: 'System:DemoFiles/demog.csv', semTypeTimeoutMs: 3000});
  await page.evaluate(() => {
    const tv = grok.shell.tv;
    tv.dataFrame.filter.setAll(true);
    tv.dataFrame.selection.setAll(false);
    tv.grid.sort([], []);
  });
  await v.waitForViewerRendered(page, 'Grid', 800);
}

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

  const consoleErrors: string[] = [];
  const pageErrors: string[] = [];
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  const realErrors = () => [...consoleErrors, ...pageErrors].filter((t) =>
    !/Unable to find element in cloned iframe/i.test(t) &&
    !/Stack trace [A-Za-z0-9]+/i.test(t));

  const rowCount = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
  expect(rowCount).toBe(5850);

  await softStep('Column Sizing from context menu: Optimal fits content grid-wide, Minimal narrows, Maximal widens', async () => {
    await reopenClean(page);
    const cell = await dataCellPoint(page, 'AGE', 3);

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

    expect(optimalMatch.matched).toBeGreaterThanOrEqual(optimalMatch.total - 1);

    for (const n of ['AGE', 'RACE', 'SEVERITY'])
      expect(modeWidths['Minimal'][n]).toBeLessThan(modeWidths['Optimal'][n]);

    expect(modeWidths['Maximal']['CONTROL']).toBeGreaterThan(modeWidths['Optimal']['CONTROL']);
    expect(modeWidths['Maximal']['SEVERITY']).toBeGreaterThan(modeWidths['Optimal']['SEVERITY']);

    for (const n of Object.keys(modeWidths['Optimal']))
      expect(modeWidths['Maximal'][n]).toBeGreaterThanOrEqual(modeWidths['Optimal'][n]);
  });

  await softStep('Header Histogram Strip: Add>Top>Histogram on then off', async () => {
    await reopenClean(page);
    const cell = await dataCellPoint(page, 'AGE', 3);

    await openContextMenu(page, cell.x, cell.y);
    await clickMenuLeaf(page, ['div-Add', 'div-Add---Top'], 'div-Add---Top---Histogram');
    await closeMenu(page);
    const on = await page.evaluate(() => grok.shell.tv.grid.getOptions(true).look.columnHeaderTypes);
    expect(on).toContain('hist');

    await openContextMenu(page, cell.x, cell.y);
    await clickMenuLeaf(page, ['div-Add', 'div-Add---Top'], 'div-Add---Top---Histogram');
    await closeMenu(page);
    const off = await page.evaluate(() => grok.shell.tv.grid.getOptions(true).look.columnHeaderTypes ?? []);
    expect(off).not.toContain('hist');
  });

  await softStep('GROK-18256: remove a summary column via the top-panel remove icon', async () => {
    await reopenClean(page);
    const cell = await dataCellPoint(page, 'AGE', 3);
    const before = await page.evaluate(() => grok.shell.tv.grid.columns.length);

    await openContextMenu(page, cell.x, cell.y);
    await clickMenuLeaf(page, ['div-Add', 'div-Add---Summary-Columns'], 'div-Add---Summary-Columns---Sparklines');
    await closeMenu(page);
    const afterAdd = await page.evaluate(() => grok.shell.tv.grid.columns.length);
    expect(afterAdd).toBe(before + 1);
    const summaryName = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      return grid.columns.byIndex(grid.columns.length - 1).name;
    });

    await page.evaluate((name) => { grok.shell.tv.grid.columns.byName(name).selected = true; }, summaryName);

    await page.waitForFunction(() => {
      const icon = document.querySelector('[name="icon-remove-selected-columns"]');
      return !!icon && !icon.classList.contains('d4-disabled');
    }, null, {timeout: 5000});
    const targetLen = before;
    await page.evaluate(() => {
      const icon = document.querySelector('[name="icon-remove-selected-columns"]') as HTMLElement;
      icon.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, view: window} as any));
      icon.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, view: window} as any));
      icon.click();
    });
    await page.waitForFunction((n) => grok.shell.tv.grid.columns.length === n, targetLen, {timeout: 5000}).catch(() => {});
    const afterRemove = await page.evaluate(() => grok.shell.tv.grid.columns.length);
    expect(afterRemove).toBe(before); 
    const summaryGone = await page.evaluate((name) => !grok.shell.tv.grid.columns.byName(name), summaryName);
    expect(summaryGone).toBe(true);
  });

  await softStep('GROK-19942: grid keeps rendering after the CONTROL source column is removed', async () => {
    await reopenClean(page);
    const cell = await dataCellPoint(page, 'AGE', 3);

    await openContextMenu(page, cell.x, cell.y);
    await clickMenuLeaf(page, ['div-Add', 'div-Add---Summary-Columns'], 'div-Add---Summary-Columns---Tags');
    await closeMenu(page);
    const errsBefore = realErrors().length;

    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.columns.remove('CONTROL');
    });
    await v.waitForViewerRendered(page, 'Grid', 400);
    await focusGrid(page);
    await page.keyboard.press('PageDown');
    await v.waitForViewerRendered(page, 'Grid', 300);
    await page.keyboard.press('PageUp');
    await v.waitForViewerRendered(page, 'Grid', 300);

    const errsAfter = realErrors().filter((t) => /render|grid|null|argument/i.test(t)).length;
    expect(errsAfter).toBe(0);
    const gridStillPresent = await page.evaluate(() =>
      !!document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]'));
    expect(gridStillPresent).toBe(true);
  });

  await softStep('GROK-19717: a project with an extracted-rows table reopens without error', async () => {
    await reopenClean(page);
    const projectName = `grid-extract-${Date.now()}`;

    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      for (let i = 0; i < 8; i++) df.selection.set(i, true);
    });
    await page.evaluate(async () => { await grok.functions.call('CmdExtractSelectedRows'); });

    await page.waitForFunction(() => grok.shell.tv?.dataFrame?.rowCount === 8, null, {timeout: 15000});
    const extractedRows = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
    expect(extractedRows).toBe(8);

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
    await v.closeAllAndWait(page);

    await page.evaluate(() => {
      (window as any).__grokErrBalloons = [];
      const obs = new MutationObserver((muts) => {
        for (const m of muts)
          for (const n of Array.from(m.addedNodes)) {
            if ((n as Element).nodeType !== 1) continue;
            const el = n as Element;
            if (el.matches && el.matches('.d4-balloon.error'))
              (window as any).__grokErrBalloons.push((el.textContent || '').slice(0, 80));
            el.querySelectorAll && el.querySelectorAll('.d4-balloon.error')
              .forEach((b) => (window as any).__grokErrBalloons.push((b.textContent || '').slice(0, 80)));
          }
      });
      obs.observe(document.body, {childList: true, subtree: true});
      (window as any).__grokErrBalloonObs = obs;
    });
    await page.evaluate(async (id) => {
      const proj = await grok.dapi.projects.find(id);
      await proj.open();
    }, savedId);

    await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
    await v.waitForViewerRendered(page, 'Grid', 1500);
    const balloons = await page.evaluate(() => {
      const obs = (window as any).__grokErrBalloonObs;
      if (obs) obs.disconnect();
      const captured = ((window as any).__grokErrBalloons || []).length;

      return captured + document.querySelectorAll('.d4-balloon.error').length;
    });
    expect(balloons).toBe(0);
    const gridPresent = await page.evaluate(() => !!document.querySelector('[name="viewer-Grid"]'));
    expect(gridPresent).toBe(true);

    const newErrs = realErrors().slice(errsBefore).filter((t) => /null|argument|render/i.test(t)).length;
    expect(newErrs).toBe(0);

    await page.evaluate(async (id) => {
      try { const proj = await grok.dapi.projects.find(id); if (proj) await grok.dapi.projects.delete(proj); } catch (_) {}
    }, savedId);
  });

  await softStep('Grid as added viewer: second grid, Row Height, Show Labels, Data>Table switch, close', async () => {
    await reopenClean(page);

    await page.evaluate(() => {
      const icon = document.querySelector('[name="icon-grid"]') as HTMLElement;
      icon.click();
    });

    await page.waitForFunction(() => document.querySelectorAll('[name="viewer-Grid"]').length === 2,
      null, {timeout: 10000});
    const gridCount = await page.evaluate(() => document.querySelectorAll('[name="viewer-Grid"]').length);
    expect(gridCount).toBe(2);
    const addedRows = await page.evaluate(() => {
      const added = grok.shell.tv.viewers.filter((vw: any) => vw.type === 'Grid');
      return added.length > 1 ? added[added.length - 1].dataFrame.rowCount : grok.shell.tv.dataFrame.rowCount;
    });
    expect(addedRows).toBe(5850);

    await page.evaluate(() => {
      const grids = document.querySelectorAll('[name="viewer-Grid"]');
      const lastGrid = grids[grids.length - 1] as HTMLElement;
      const gear = lastGrid.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
      gear.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, view: window} as any));
      gear.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, view: window} as any));
      gear.dispatchEvent(new MouseEvent('click', {bubbles: true, view: window} as any));
    });
    await page.waitForTimeout(700); 

    const rowHeightChanged = await page.evaluate(async () => {
      const added = grok.shell.tv.viewers.filter((vw: any) => vw.type === 'Grid');
      const addedGrid = added[added.length - 1];
      const mainGrid = grok.shell.tv.grid;
      const col0 = grok.shell.tv.dataFrame.columns.byIndex(0).name;
      const addedCellBefore = addedGrid.cell(col0, 0)?.documentBounds?.height ?? null;
      const mainCellBefore = mainGrid.cell(col0, 0)?.documentBounds?.height ?? null;

      const settled = new Promise<void>((resolve) => {
        let sub: any = null;
        try { sub = addedGrid.onViewerRendered.subscribe(() => { sub.unsubscribe(); resolve(); }); }
        catch (_) {  }
        setTimeout(() => { try { sub?.unsubscribe(); } catch (_) {} resolve(); }, 400);
      });
      addedGrid.props.rowHeight = 40;
      await settled;
      const addedCellAfter = addedGrid.cell(col0, 0)?.documentBounds?.height ?? null;
      const mainCellAfter = mainGrid.cell(col0, 0)?.documentBounds?.height ?? null;
      return {addedCellBefore, addedCellAfter, mainCellBefore, mainCellAfter};
    });

    expect(rowHeightChanged.addedCellAfter).toBe(40);
    expect(rowHeightChanged.addedCellAfter).toBeGreaterThan(rowHeightChanged.addedCellBefore);

    expect(rowHeightChanged.mainCellAfter).toBe(rowHeightChanged.mainCellBefore);

    const labelsToggled = await page.evaluate(async () => {
      const added = grok.shell.tv.viewers.filter((vw: any) => vw.type === 'Grid');
      const addedGrid = added[added.length - 1];
      const headerBefore = addedGrid.colHeaderHeight;
      const settle = (grid: any, cap: number) => new Promise<void>((resolve) => {
        let sub: any = null;
        try { sub = grid.onViewerRendered.subscribe(() => { sub.unsubscribe(); resolve(); }); }
        catch (_) {  }
        setTimeout(() => { try { sub?.unsubscribe(); } catch (_) {} resolve(); }, cap);
      });
      addedGrid.props.showColumnLabels = false;
      await settle(addedGrid, 300);
      const headerHidden = addedGrid.colHeaderHeight;
      addedGrid.props.showColumnLabels = true;
      await settle(addedGrid, 300);
      const headerShown = addedGrid.colHeaderHeight;
      return {headerBefore, headerHidden, headerShown};
    });
    expect(labelsToggled.headerBefore).toBeGreaterThan(0);
    expect(labelsToggled.headerHidden).toBe(0);
    expect(labelsToggled.headerShown).toBeGreaterThan(0);

    const rebind = await page.evaluate(async () => {
      const added = grok.shell.tv.viewers.filter((vw: any) => vw.type === 'Grid');
      const addedGrid = added[added.length - 1];
      const spgi = await grok.dapi.files.readCsv('System:AppData/Chem/tests/spgi-100.csv');
      grok.shell.addTableView(spgi);
      await new Promise((r) => setTimeout(r, 800)); 
      addedGrid.dataFrame = spgi;
      const settled = new Promise<void>((resolve) => {
        let sub: any = null;
        try { sub = addedGrid.onViewerRendered.subscribe(() => { sub.unsubscribe(); resolve(); }); }
        catch (_) {  }
        setTimeout(() => { try { sub?.unsubscribe(); } catch (_) {} resolve(); }, 500);
      });
      await settled;
      return {rows: addedGrid.dataFrame.rowCount, spgiRows: spgi.rowCount};
    });
    expect(rebind.rows).toBe(rebind.spgiRows);

    await page.evaluate(async () => {
      const grids = grok.shell.tv.viewers.filter((vw: any) => vw.type === 'Grid');
      if (grids.length > 1) grids[grids.length - 1].close();
    });
    await page.waitForFunction(() =>
      grok.shell.tv.viewers.filter((vw: any) => vw.type === 'Grid').length === 1,
      null, {timeout: 5000}).catch(() => {});
    const stillHasGrid = await page.evaluate(() => !!grok.shell.tv.grid);
    expect(stillHasGrid).toBe(true);
  });

  await softStep('Multi-Column + Row-Height Resize: linked columns, single column, row height, hide-by-width', async () => {
    await reopenClean(page);

    await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      grid.columns.byName('AGE').selected = true;
      grid.columns.byName('HEIGHT').selected = true;
    });
    await v.waitForViewerRendered(page, 'Grid', 400);
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
    await v.waitForViewerRendered(page, 'Grid', 400);
    const linked = await page.evaluate(() => ({
      age: grok.shell.tv.grid.columns.byName('AGE').width,
      height: grok.shell.tv.grid.columns.byName('HEIGHT').width,
      sex: grok.shell.tv.grid.columns.byName('SEX').width}));

    expect(linked.age).toBe(linked.height);
    expect(linked.age).toBeGreaterThan(53);
    expect(linked.sex).toBe(sexWidthBefore);

    await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      grid.columns.byName('AGE').selected = false;
      grid.columns.byName('HEIGHT').selected = false;
    });
    await focusGrid(page);
    await page.keyboard.press('Escape');
    await v.waitForViewerRendered(page, 'Grid', 300);
    const ageOnlyBefore = await page.evaluate(() => grok.shell.tv.grid.columns.byName('AGE').width);
    const heightOnlyBefore = await page.evaluate(() => grok.shell.tv.grid.columns.byName('HEIGHT').width);
    const ageBorder2 = await colBorderPoint(page, 'AGE');
    await focusGrid(page);
    await page.mouse.move(ageBorder2.x, ageBorder2.y);
    await page.mouse.down();
    await page.mouse.move(ageBorder2.x + 30, ageBorder2.y, {steps: 5});
    await page.mouse.up();
    await v.waitForViewerRendered(page, 'Grid', 400);
    const single = await page.evaluate(() => ({
      age: grok.shell.tv.grid.columns.byName('AGE').width,
      height: grok.shell.tv.grid.columns.byName('HEIGHT').width}));
    expect(single.age).toBeGreaterThan(ageOnlyBefore);
    expect(single.height).toBe(heightOnlyBefore);

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
    await v.waitForViewerRendered(page, 'Grid', 400);
    const rhAfter = await page.evaluate(() => grok.shell.tv.grid.props.rowHeight);
    expect(rhAfter).toBe(rhBefore + 30);
    const widthsAfter = await page.evaluate(() => ({
      age: grok.shell.tv.grid.columns.byName('AGE').width,
      height: grok.shell.tv.grid.columns.byName('HEIGHT').width}));
    expect(widthsAfter).toEqual(widthsBefore);

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
    await v.waitForViewerRendered(page, 'Grid', 400);
    const weightHidden = await page.evaluate(() => ({
      width: grok.shell.tv.grid.columns.byName('WEIGHT').width,
      visible: grok.shell.tv.grid.columns.byName('WEIGHT').visible}));
    expect(weightHidden.width).toBeLessThanOrEqual(2);
    expect(weightHidden.visible).toBe(true);
  });

  await softStep('Column Tooltip Settings: Current Column radios (Default/None/Columns), tooltip DOM display', async () => {
    await reopenClean(page);
    const ageHdr = await headerPoint(page, 'AGE');

    const ttBaseline = await page.evaluate(() => grok.shell.tv.grid.columns.byName('AGE').tooltipType);
    expect(ttBaseline).toBeNull();
    await openContextMenu(page, ageHdr.x, ageHdr.y);
    await expandSubmenu(page, 'div-Tooltip');
    await expandSubmenu(page, 'div-Tooltip---Current-Column');
    const ttChoices = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-menu-popup [name^="div-Tooltip---Current-Column---"]'))
        .map((e) => e.getAttribute('name')!.split('---').pop()));
    expect(ttChoices).toEqual(expect.arrayContaining(['Default', 'Form', 'Columns', 'None']));

    const defaultMarked = await page.evaluate(() =>
      !!document.querySelector('.d4-menu-popup [name="div-Tooltip---Current-Column---Default"] [name="icon-dot-circle"]'));
    expect(defaultMarked).toBe(true);

    await page.evaluate(() => {
      const leaf = document.querySelector('.d4-menu-popup [name="div-Tooltip---Current-Column---None"]') as HTMLElement;
      leaf.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, view: window} as any));
      leaf.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, view: window} as any));
      leaf.dispatchEvent(new MouseEvent('click', {bubbles: true, view: window} as any));
    });
    await closeMenu(page);
    const ttNone = await page.evaluate(() => grok.shell.tv.grid.columns.byName('AGE').tooltipType);
    expect(ttNone).toBe('None');

    const ageCell = await dataCellPoint(page, 'AGE', 3);
    await page.mouse.move(ageCell.x, ageCell.y);
    await page.waitForTimeout(900); 
    const noneDisplay = await page.evaluate(() => {
      const t = document.querySelector('.d4-tooltip') as HTMLElement | null;
      return t ? getComputedStyle(t).display : 'absent';
    });
    expect(noneDisplay === 'none' || noneDisplay === 'absent').toBe(true);

    await page.mouse.move(ageCell.x, ageCell.y - 200);

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
    await page.waitForSelector('[name="dialog-Select-columns..."]', {state: 'detached', timeout: 5000}).catch(() => {});
    const dlgClosed = await page.evaluate(() => !document.querySelector('[name="dialog-Select-columns..."]'));
    expect(dlgClosed).toBe(true);
    const ttColumns = await page.evaluate(() => grok.shell.tv.grid.columns.byName('AGE').tooltipType);
    expect(ttColumns).toBe('Columns');

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

      const tableRow = grok.shell.tv.grid.gridRowToTable(3);
      const t = document.querySelector('.d4-tooltip') as HTMLElement | null;
      const text = t ? (t.innerText ?? '') : '';
      const display = t ? getComputedStyle(t).display : 'absent';
      const ageVal = df.columns.byName('AGE').get(tableRow);
      const sexVal = df.columns.byName('SEX').get(tableRow);
      return {display, text, ageVal: String(ageVal), sexVal: String(sexVal)};
    });

    expect(colsTip.display).not.toBe('none');
    expect(colsTip.text).toContain('AGE');
    expect(colsTip.text).toContain('SEX');
    expect(colsTip.text.includes(colsTip.ageVal) || colsTip.text.includes(colsTip.sexVal)).toBe(true);
    await page.mouse.move(ageCell2.x, ageCell2.y - 200);

    await openContextMenu(page, ageHdr.x, ageHdr.y);
    await clickMenuLeaf(page, ['div-Tooltip', 'div-Tooltip---Current-Column'], 'div-Tooltip---Current-Column---Default');
    await closeMenu(page);
    const ttDefault = await page.evaluate(() => grok.shell.tv.grid.columns.byName('AGE').tooltipType);
    expect(ttDefault).toBe('Default');
  });

  const leftover = realErrors();
  expect(leftover, `unexpected console/page errors: ${leftover.join(' | ')}`).toEqual([]);

  v.finishSpec();
});
