/* ---
realizes: []
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {isLocalBootNoise, openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import * as g from './grid-helpers';

declare const grok: any;

// The summary-column steps (GROK-18256, GROK-19942: PowerGrid cell renderers) and the extracted-rows
// project step (GROK-19717) run on the server lane in grid-server-spec.ts.
test.use(specTestOptions);

async function reopenClean(page: Page): Promise<void> {
  await v.closeAllAndWait(page);
  await v.openTable(page, {path: g.DEMOG, semTypeTimeoutMs: 3000});
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

const readWidths = (page: Page) => page.evaluate(() => {
  const grid = grok.shell.tv.grid;
  const o: Record<string, number> = {};
  for (const n of grid.dataFrame.columns.names()) o[n] = grid.columns.byName(n).width;
  return o;
});

const headerTypes = (page: Page) => page.evaluate(() =>
  (grok.shell.tv.grid.getOptions(true).look.columnHeaderTypes ?? []) as string[]);

const tooltipState = (page: Page) => page.evaluate(() => {
  const t = document.querySelector('.d4-tooltip') as HTMLElement | null;
  return {display: t ? getComputedStyle(t).display : 'absent', text: t ? (t.innerText ?? '') : ''};
});

test('Grid tests', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);
  const flags = await g.readShellFlags(page);
  const errors = g.trackErrors(page, (t) => isLocalBootNoise(t) ||
    /Unable to find element in cloned iframe/i.test(t) || /Stack trace [A-Za-z0-9]+/i.test(t));
  try {
    await v.openTable(page, {path: g.DEMOG, semTypeTimeoutMs: 3000});

    const rowCount = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
    expect(rowCount).toBe(5850);

    await softStep('Column Sizing from context menu: Optimal fits content grid-wide, Minimal narrows, Maximal widens', async () => {
      await reopenClean(page);
      const cell = await g.cellCenter(page, 'AGE', 3);

      const modeWidths: Record<string, Record<string, number>> = {};
      for (const leaf of ['Optimal', 'Minimal', 'Maximal']) {
        const before = JSON.stringify(await readWidths(page));
        expect(await g.clickMenuLeaf(page, cell, ['div-Column-Sizing'], `div-Column-Sizing---${leaf}`)).toBe(true);
        modeWidths[leaf] = await v.pollValue(() => readWidths(page), (w) => JSON.stringify(w) !== before, 500, 50);
      }

      const before = JSON.stringify(await readWidths(page));
      expect(await g.clickMenuLeaf(page, cell, ['div-Column-Sizing'], 'div-Column-Sizing---Optimal')).toBe(true);
      await v.pollValue(() => readWidths(page), (w) => JSON.stringify(w) !== before, 500, 50);
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
      const cell = await g.cellCenter(page, 'AGE', 3);

      expect(await g.clickMenuLeaf(page, cell, ['div-Add', 'div-Add---Top'], 'div-Add---Top---Histogram')).toBe(true);
      const on = await v.pollValue(() => headerTypes(page), (t) => t.includes('hist'), 1500, 50);
      expect(on).toContain('hist');

      expect(await g.clickMenuLeaf(page, cell, ['div-Add', 'div-Add---Top'], 'div-Add---Top---Histogram')).toBe(true);
      const off = await v.pollValue(() => headerTypes(page), (t) => !t.includes('hist'), 1500, 50);
      expect(off).not.toContain('hist');
    });

    await softStep('Grid as added viewer: second grid, Row Height, Show Labels, Data>Table switch, close', async () => {
      await reopenClean(page);

      await page.evaluate(() => (document.querySelector('[name="icon-grid"]') as HTMLElement).click());
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
      await page.locator('.property-grid').first().waitFor({state: 'attached', timeout: 700}).catch(() => {});

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
        const w = window as any;
        const added = grok.shell.tv.viewers.filter((vw: any) => vw.type === 'Grid');
        const addedGrid = added[added.length - 1];
        const spgi = await w.__readCsv('System:AppData/Chem/tests/spgi-100.csv');
        grok.shell.addTableView(spgi);
        await w.__tableReady(3000);
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

      const widthOf = (col: string) => page.evaluate((c) => grok.shell.tv.grid.columns.byName(c).width, col);
      const sexWidthBefore = await widthOf('SEX');
      const ageWidthBefore = await widthOf('AGE');
      const ageBorder = await colBorderPoint(page, 'AGE');
      await g.focusGrid(page);
      await page.mouse.move(ageBorder.x, ageBorder.y);
      await page.mouse.down();
      await page.mouse.move(ageBorder.x + 40, ageBorder.y, {steps: 6});
      await page.mouse.up();
      const linked = await v.pollValue(() => page.evaluate(() => ({
        age: grok.shell.tv.grid.columns.byName('AGE').width,
        height: grok.shell.tv.grid.columns.byName('HEIGHT').width,
        sex: grok.shell.tv.grid.columns.byName('SEX').width})), (x) => x.age !== ageWidthBefore, 400, 50);
      expect(linked.age).toBe(linked.height);
      expect(linked.age).toBeGreaterThan(53);
      expect(linked.sex).toBe(sexWidthBefore);

      await page.evaluate(() => {
        const grid = grok.shell.tv.grid;
        grid.columns.byName('AGE').selected = false;
        grid.columns.byName('HEIGHT').selected = false;
      });
      await g.focusGrid(page);
      await page.keyboard.press('Escape');
      await v.waitForViewerRendered(page, 'Grid', 300);
      const ageOnlyBefore = await widthOf('AGE');
      const heightOnlyBefore = await widthOf('HEIGHT');
      const ageBorder2 = await colBorderPoint(page, 'AGE');
      await g.focusGrid(page);
      await page.mouse.move(ageBorder2.x, ageBorder2.y);
      await page.mouse.down();
      await page.mouse.move(ageBorder2.x + 30, ageBorder2.y, {steps: 5});
      await page.mouse.up();
      const single = await v.pollValue(() => page.evaluate(() => ({
        age: grok.shell.tv.grid.columns.byName('AGE').width,
        height: grok.shell.tv.grid.columns.byName('HEIGHT').width})), (x) => x.age !== ageOnlyBefore, 400, 50);
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
      await g.focusGrid(page);
      await page.mouse.move(rhGeom.x, rhGeom.y);
      await page.mouse.down();
      await page.mouse.move(rhGeom.x, rhGeom.y + 30, {steps: 6});
      await page.mouse.up();
      const rhAfter = await v.pollValue(() => page.evaluate(() => grok.shell.tv.grid.props.rowHeight), (h) => h !== rhBefore, 400, 50);
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
      await g.focusGrid(page);
      await page.mouse.move(weightHideGeom.borderX, weightHideGeom.y);
      await page.mouse.down();
      await page.mouse.move(weightHideGeom.leftEdgeX, weightHideGeom.y, {steps: 6});
      await page.mouse.up();
      const weightHidden = await v.pollValue(() => page.evaluate(() => ({
        width: grok.shell.tv.grid.columns.byName('WEIGHT').width,
        visible: grok.shell.tv.grid.columns.byName('WEIGHT').visible})), (x) => x.width <= 2, 400, 50);
      expect(weightHidden.width).toBeLessThanOrEqual(2);
      expect(weightHidden.visible).toBe(true);
    });

    await softStep('Column Tooltip Settings: Current Column radios (Default/None/Columns), tooltip DOM display', async () => {
      await reopenClean(page);
      const ageHdr = await g.headerCenter(page, 'AGE');

      const ttBaseline = await page.evaluate(() => grok.shell.tv.grid.columns.byName('AGE').tooltipType);
      expect(ttBaseline).toBeNull();
      await g.openGridMenu(page, ageHdr);
      const menu = await page.evaluate(() => {
        const show = (name: string) => {
          const group = document.querySelector(`.d4-menu-popup [name="${name}"]`) as HTMLElement | null;
          const container = group?.querySelector('.d4-menu-item-container.d4-vert-menu') as HTMLElement | null;
          if (container) container.style.display = 'flex';
          return !!group;
        };
        const ok = show('div-Tooltip') && show('div-Tooltip---Current-Column');
        return {
          ok,
          choices: Array.from(document.querySelectorAll('.d4-menu-popup [name^="div-Tooltip---Current-Column---"]'))
            .map((e) => e.getAttribute('name')!.split('---').pop()),
          defaultMarked: !!document.querySelector('.d4-menu-popup [name="div-Tooltip---Current-Column---Default"] [name="icon-dot-circle"]'),
        };
      });
      expect(menu.ok).toBe(true);
      expect(menu.choices).toEqual(expect.arrayContaining(['Default', 'Form', 'Columns', 'None']));
      expect(menu.defaultMarked).toBe(true);

      await page.evaluate(() => {
        const leaf = document.querySelector('.d4-menu-popup [name="div-Tooltip---Current-Column---None"]') as HTMLElement;
        leaf.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, view: window} as any));
        leaf.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, view: window} as any));
        leaf.dispatchEvent(new MouseEvent('click', {bubbles: true, view: window} as any));
      });
      await g.closeGridMenu(page);
      const ttNone = await v.pollValue(() => page.evaluate(() => grok.shell.tv.grid.columns.byName('AGE').tooltipType),
        (t) => t === 'None', 1000, 50);
      expect(ttNone).toBe('None');

      const ageCell = await g.cellCenter(page, 'AGE', 3);
      await page.mouse.move(ageCell.x, ageCell.y);
      // a tooltip that must NOT appear: the hover window is the assertion, capped at the old wait
      const none = await v.pollValue(() => tooltipState(page), (t) => t.display !== 'none' && t.display !== 'absent', 900, 100);
      expect(none.display === 'none' || none.display === 'absent').toBe(true);

      await page.mouse.move(ageCell.x, ageCell.y - 200);

      expect(await g.clickMenuLeaf(page, ageHdr, ['div-Tooltip', 'div-Tooltip---Current-Column'], 'div-Tooltip---Current-Column---Columns')).toBe(true);
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

      const ageCell2 = await g.cellCenter(page, 'AGE', 3);
      await page.mouse.move(ageCell2.x, ageCell2.y);
      const shown = await v.pollValue(() => tooltipState(page),
        (t) => t.display !== 'none' && t.display !== 'absent' && t.text.includes('AGE') && t.text.includes('SEX'), 1500, 100);
      const rowVals = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const tableRow = grok.shell.tv.grid.gridRowToTable(3);
        return {ageVal: String(df.columns.byName('AGE').get(tableRow)), sexVal: String(df.columns.byName('SEX').get(tableRow))};
      });
      expect(shown.display).not.toBe('none');
      expect(shown.text).toContain('AGE');
      expect(shown.text).toContain('SEX');
      expect(shown.text.includes(rowVals.ageVal) || shown.text.includes(rowVals.sexVal)).toBe(true);
      await page.mouse.move(ageCell2.x, ageCell2.y - 200);

      expect(await g.clickMenuLeaf(page, ageHdr, ['div-Tooltip', 'div-Tooltip---Current-Column'], 'div-Tooltip---Current-Column---Default')).toBe(true);
      const ttDefault = await v.pollValue(() => page.evaluate(() => grok.shell.tv.grid.columns.byName('AGE').tooltipType),
        (t) => t === 'Default', 1000, 50);
      expect(ttDefault).toBe('Default');
    });

    expect(errors.list, `unexpected console/page errors: ${errors.list.join(' | ')}`).toEqual([]);
  } finally {
    errors.stop();
    await g.leaveShellClean(page, flags);
  }
  v.finishSpec();
});
