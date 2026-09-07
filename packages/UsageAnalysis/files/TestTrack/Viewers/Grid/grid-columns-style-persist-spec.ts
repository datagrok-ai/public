/* ---
realizes: [grid.cp.columns-layout-persist]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import * as g from './grid-helpers';

declare const grok: any;

// Steps 14 and 16 (layout and project round-trips) run on the server lane in grid-server-spec.ts.
test.use(specTestOptions);

async function armGridEvent(page: Page, eventProp: string, key: string): Promise<void> {
  await page.evaluate(({eventProp, key}) => {
    const w = window as any;
    w.__gridEventFired = w.__gridEventFired ?? {};
    w.__gridEventFired[key] = false;
    const grid = w.grok.shell.tv.grid;
    const sub = grid[eventProp].subscribe(() => { w.__gridEventFired[key] = true; sub.unsubscribe(); });
  }, {eventProp, key});
}
async function awaitGridEvent(page: Page, key: string, capMs: number): Promise<void> {
  await page.evaluate(({key, capMs}) => new Promise<void>((resolve) => {
    const w = window as any;
    const t0 = Date.now();
    const tick = () => {
      if ((w.__gridEventFired?.[key]) || Date.now() - t0 >= capMs) { resolve(); return; }
      setTimeout(tick, 25);
    };
    tick();
  }), {key, capMs});
}

test('Grid — Column Geometry: Sort, Order, Visibility, Width, Pinning and Persistence', async ({page}) => {
  test.setTimeout(240_000);

  await openDatagrok(page);
  const flags = await g.readShellFlags(page);
  const errors = g.trackErrors(page);
  try {
    await v.openTable(page, {path: g.DEMOG, semTypeTimeoutMs: 3000});

    const baseline = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const df = grok.shell.tv.dataFrame;
      return {ageRow0: df.col('AGE').get(0), frozenColumns: grid.props.frozenColumns};
    });
    expect(baseline.frozenColumns).toBe(1);

    await softStep('Step 4 — Sort: first double-click on AGE header sorts DESCENDING', async () => {
      const c = await g.headerCenter(page, 'AGE');
      await armGridEvent(page, 'onRowsSorted', 'sort4');
      await page.mouse.dblclick(c.x, c.y);
      await awaitGridEvent(page, 'sort4', 600);
      const r = await page.evaluate(() => {
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
      });
      expect(r.sortBy).toContain('AGE');
      expect(r.sortTypes).toContain(false);
      expect(r.topAge).toBe(r.maxAge);
      expect(r.ageRow0).toBe(baseline.ageRow0);
    });

    await softStep('Step 5 — Sort: second double-click on AGE header sorts ASCENDING', async () => {
      const c = await g.headerCenter(page, 'AGE');
      await armGridEvent(page, 'onRowsSorted', 'sort5');
      await page.mouse.dblclick(c.x, c.y);
      await awaitGridEvent(page, 'sort5', 600);
      const r = await page.evaluate(() => {
        const grid = grok.shell.tv.grid;
        const df = grok.shell.tv.dataFrame;
        const ageCol = df.col('AGE');
        let minNonNull = Number.POSITIVE_INFINITY;
        for (let i = 0; i < df.rowCount; i++)
          if (!ageCol.isNone(i)) { const val = ageCol.get(i); if (val < minNonNull) minNonNull = val; }
        const topDataIdx = grid.gridRowToTable(0);
        return {
          sortBy: grid.props.sortByColumnNames,
          sortTypes: grid.props.sortTypes,
          topAge: ageCol.get(topDataIdx),
          minNonNull,
          ageRow0: ageCol.get(0),
        };
      });
      expect(r.sortBy).toContain('AGE');
      expect(r.sortTypes).toContain(true);
      expect(r.topAge).toBe(r.minNonNull);
      expect(r.ageRow0).toBe(baseline.ageRow0);
    });

    await softStep('Step 6 — Sort: third double-click on AGE header RESETS the sort', async () => {
      const ascTop = await page.evaluate(() => grok.shell.tv.grid.gridRowToTable(0));
      const c = await g.headerCenter(page, 'AGE');
      await armGridEvent(page, 'onRowsSorted', 'sort6');
      await page.mouse.dblclick(c.x, c.y);
      await awaitGridEvent(page, 'sort6', 600);
      const r = await page.evaluate((prevTop) => {
        const grid = grok.shell.tv.grid;
        return {sortBy: grid.props.sortByColumnNames, topDataIdx: grid.gridRowToTable(0), prevTop};
      }, ascTop);
      expect(r.sortBy).toEqual([]);
      expect(r.topDataIdx).not.toBe(r.prevTop);
    });

    await softStep('Sort: leave grid sorted ascending on AGE for the following scenarios', async () => {
      const c = await g.headerCenter(page, 'AGE');
      await armGridEvent(page, 'onRowsSorted', 'sortLeaveDesc');
      await page.mouse.dblclick(c.x, c.y);
      await awaitGridEvent(page, 'sortLeaveDesc', 400);
      await armGridEvent(page, 'onRowsSorted', 'sortLeaveAsc');
      await page.mouse.dblclick(c.x, c.y);
      await awaitGridEvent(page, 'sortLeaveAsc', 500);
      const sortTypes = await page.evaluate(() => grok.shell.tv.grid.props.sortTypes);
      expect(sortTypes).toContain(true);
    });

    await softStep('Step 7 — Reorder: drag HEIGHT header to the right of its current slot', async () => {
      const readOrder = () => page.evaluate(() => {
        const order = Array.from({length: grok.shell.tv.grid.columns.length}, (_: any, i: number) => grok.shell.tv.grid.columns.byIndex(i).name) as string[];
        return {order, heightIdx: order.indexOf('HEIGHT')};
      });
      const before = await readOrder();
      const src = await g.headerCenter(page, 'HEIGHT');
      const tgt = await g.headerCenter(page, 'DEMOG');

      await page.mouse.move(src.x, src.y);
      await page.mouse.down();
      await page.mouse.move((src.x + tgt.x) / 2, src.y, {steps: 2});
      await page.mouse.move(tgt.x, tgt.y, {steps: 3});
      await page.mouse.up();
      const after = await v.pollValue(readOrder, (x) => JSON.stringify(x.order) !== JSON.stringify(before.order), 700, 50);
      expect(after.heightIdx).toBeGreaterThan(before.heightIdx);
      expect(after.order).not.toEqual(before.order);
    });

    await softStep('Step 9 — Hide: open Order or Hide Columns and hide WEIGHT', async () => {
      expect(await g.clickMenuLeaf(page, await g.cellCenter(page, 'AGE', 0), [], 'div-Order-or-Hide-Columns...')).toBe(true);
      await page.locator('.d4-dialog .d4-dialog-header', {hasText: 'Order or Hide Columns'}).waitFor({timeout: 5000});

      const r = await page.evaluate(() => {
        const grid = grok.shell.tv.grid;
        const df = grok.shell.tv.dataFrame;
        grid.columns.setVisible(df.columns.names().filter((n: string) => n !== 'WEIGHT'));
        return true;
      });
      expect(r).toBe(true);

      await page.locator('.d4-dialog [name="button-CLOSE"]').first().click({timeout: 5000}).catch(() => {});
      const state = await v.pollValue(() => page.evaluate(() => {
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
      }), (x) => !x.weightVisible && !x.weightEnumerated, 400, 50);
      expect(state.weightVisible).toBe(false);
      expect(state.weightEnumerated).toBe(false);
      expect(state.surroundingKept).toBe(true);
    });

    await softStep('Step 10 — Resize: widen AGE by dragging its right header border', async () => {
      const readWidth = () => page.evaluate(() => grok.shell.tv.grid.columns.byName('AGE').width);
      const before = await readWidth();
      const geom = await page.evaluate(() => {
        const grid = grok.shell.tv.grid;
        const db = grid.cell('AGE', 0).documentBounds;
        return {borderX: db.x + db.width, headerY: db.y - grid.colHeaderHeight / 2};
      });

      await page.mouse.move(geom.borderX, geom.headerY);
      await page.mouse.down();
      await page.mouse.move(geom.borderX + 60, geom.headerY, {steps: 3});
      await page.mouse.up();
      const after = await v.pollValue(readWidth, (w) => w > before, 500, 50);
      expect(after).toBeGreaterThan(before);
    });

    await softStep('Step 10 (cont.) — Resize + scroll: widen more columns, scroll horizontally, assert no errors (GROK-19753)', async () => {
      const errBefore = errors.count();
      await page.evaluate(async () => {
        const grid = grok.shell.tv.grid;
        for (const n of ['USUBJID', 'RACE', 'DIS_POP', 'STARTED']) {
          const c = grid.columns.byName(n);
          if (c) c.width = 220;
        }
        await (window as any).__quiet('viewer:Grid.onAfterDrawContent', 200, 400);
      });

      await page.evaluate(async () => {
        const grid = grok.shell.tv.grid;
        grid.horzScroll.scrollTo(grid.horzScroll.maxRange);
        await (window as any).__quiet('viewer:Grid.onAfterDrawContent', 200, 600);
      }).catch(async () => {
        await page.evaluate(async () => {
          const grid = grok.shell.tv.grid;
          if (grid.horzScroll?.setValues) grid.horzScroll.setValues(grid.horzScroll.min, grid.horzScroll.max, grid.horzScroll.max - 3, grid.horzScroll.max);
          await (window as any).__quiet('viewer:Grid.onAfterDrawContent', 200, 600);
        });
      });
      // the error window after the scroll is the assertion; the grid's own paint closes it
      await v.waitForGridPainted(page, {gapMs: 150, capMs: 600});
      const consistent = await page.evaluate(() => {
        const grid = grok.shell.tv.grid;
        let ok = true;
        for (let i = 0; i < grid.columns.length; i++) {
          const c = grid.columns.byIndex(i);
          if (!c.visible || !c.name) continue;
          try { void grid.cell(c.name, 0).documentBounds; } catch (_) { ok = false; }
        }
        return ok;
      });
      const gridErrors = errors.list.slice(errBefore).filter((e) => /grid/i.test(e) || /index/i.test(e));
      expect(consistent).toBe(true);
      expect(gridErrors).toEqual([]);
    });

    await softStep('Step 11 — Pin: pin SEX column via the header Pin menu', async () => {
      const frozenBefore = await page.evaluate(() => grok.shell.tv.grid.props.frozenColumns);

      const geom = await page.evaluate(async () => {
        const grid = grok.shell.tv.grid;
        grid.scrollToCell('SEX', 0);
        await (window as any).__quiet('viewer:Grid.onAfterDrawContent', 200, 600);
        const db = grid.cell('SEX', 0).documentBounds;
        const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
        const orect = overlay.getBoundingClientRect();
        return {clickX: db.x + db.width / 2, overlayLeft: orect.left};
      });
      expect(geom.clickX).toBeGreaterThan(geom.overlayLeft + 20);

      const c = await g.headerCenter(page, 'SEX');
      expect(await g.pinViaMenu(page, c, 'div-Pin---Pin-Column')).toBe(true);
      const frozenAfter = await v.pollValue(() => page.evaluate(() => grok.shell.tv.grid.props.frozenColumns),
        (f) => f === frozenBefore + 1, 1000, 50);
      expect(frozenAfter).toBe(frozenBefore + 1);
    });

    await softStep('Step 12 — Pin: pin two rows via the row Pin menu', async () => {
      await page.evaluate(async () => {
        const grid = grok.shell.tv.grid;
        grid.scrollToCell('AGE', 0);
        await (window as any).__quiet('viewer:Grid.onAfterDrawContent', 200, 600);
        await new Promise<void>((resolve) => {
          const sub = grid.onRowsSorted.subscribe(() => { sub.unsubscribe(); resolve(); });
          setTimeout(() => { sub.unsubscribe(); resolve(); }, 600);
          grid.sort([], []);
        });
      });

      const cell0 = await g.cellCenter(page, 'AGE', 0);
      await armGridEvent(page, 'onPinnedRowsChanged', 'pinRow1');
      expect(await g.pinViaMenu(page, cell0, 'div-Pin---Pin-Row')).toBe(true);
      await awaitGridEvent(page, 'pinRow1', 500);

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
      await armGridEvent(page, 'onPinnedRowsChanged', 'pinRow2');
      expect(await g.pinViaMenu(page, cell1, 'div-Pin---Pin-Row')).toBe(true);
      await awaitGridEvent(page, 'pinRow2', 500);

      const r = await page.evaluate(() => {
        const grid = grok.shell.tv.grid;
        const ageCol = grok.shell.tv.dataFrame.col('AGE');
        const rows = Array.from(grid.pinnedRows) as number[];
        const ages = rows.map((ti) => ageCol.get(ti));
        return {len: rows.length, rows, distinctAges: new Set(ages).size};
      });
      expect(r.len).toBe(2);
      expect(r.rows.length).toBe(2);
      expect(r.distinctAges).toBe(2);

      await page.evaluate(async () => {
        const grid = grok.shell.tv.grid;
        await new Promise<void>((resolve) => {
          const sub = grid.onRowsSorted.subscribe(() => { sub.unsubscribe(); resolve(); });
          setTimeout(() => { sub.unsubscribe(); resolve(); }, 600);
          grid.sort(['AGE'], [true]);
        });
      });
      const sortTypes = await page.evaluate(() => grok.shell.tv.grid.props.sortTypes);
      expect(sortTypes).toContain(true);
    });
  } finally {
    errors.stop();
    await g.leaveShellClean(page, flags);
  }
  v.finishSpec();
});
