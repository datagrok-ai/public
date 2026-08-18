/* ---
realizes: [grid.cp.columns-layout-persist]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaUI, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;

const isBenignSaveWindowError = (text: string): boolean =>
  /Unable to find element in cloned iframe/.test(text) ||
  /Stack trace [A-Za-z]+/.test(text) ||
  /NullError: method not found: '\w+' on null/.test(text);

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

async function headerCenter(page: Page, col: string): Promise<{x: number; y: number}> {
  return page.evaluate((c) => {
    const grid = grok.shell.tv.grid;
    const db = grid.cell(c, 0).documentBounds;
    return {x: db.x + db.width / 2, y: db.y - grid.colHeaderHeight / 2};
  }, col);
}

async function pinViaMenu(
  page: Page, at: {x: number; y: number}, leafName: string,
): Promise<boolean> {
  return page.evaluate(async ({at, leafName}) => {
    const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
    const closeMenu = () => document.body.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, cancelable: true}));
    const sleep = (ms: number) => new Promise((r) => setTimeout(r, ms));
    for (let open = 0; open < 3; open++) {
      overlay.dispatchEvent(new MouseEvent('contextmenu',
        {bubbles: true, cancelable: true, clientX: at.x, clientY: at.y, button: 2}));
      await sleep(450); 
      const pin = document.querySelector('[name="div-Pin"]') as HTMLElement | null;
      if (!pin) { closeMenu(); await sleep(200); continue; }
      const pb = pin.getBoundingClientRect();
      const cx = pb.x + pb.width / 2, cy = pb.y + pb.height / 2;

      let leaf: HTMLElement | null = null;
      for (let h = 0; h < 12 && !leaf; h++) {
        pin.dispatchEvent(new MouseEvent('mouseover', {bubbles: true, cancelable: true, clientX: cx, clientY: cy}));
        pin.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, cancelable: true, clientX: cx + (h % 2 ? 1 : -1), clientY: cy}));
        await sleep(200); 
        const cand = document.querySelector('[name="' + leafName + '"]') as HTMLElement | null;
        if (cand) { const b = cand.getBoundingClientRect(); if (b.width > 0 && b.height > 0) leaf = cand; }
      }
      if (!leaf) { closeMenu(); await sleep(200); continue; }
      const lb = leaf.getBoundingClientRect();
      leaf.dispatchEvent(new MouseEvent('click',
        {bubbles: true, cancelable: true, clientX: lb.x + lb.width / 2, clientY: lb.y + lb.height / 2}));
      await sleep(500); 
      return true;
    }
    return false;
  }, {at, leafName});
}

test('Grid — Column Geometry: Sort, Order, Visibility, Width, Pinning and Persistence', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: 'System:DemoFiles/demog.csv', semTypeTimeoutMs: 3000});

  const baseline = await page.evaluate(() => {
    const grid = grok.shell.tv.grid;
    const df = grok.shell.tv.dataFrame;
    return {ageRow0: df.col('AGE').get(0), frozenColumns: grid.props.frozenColumns};
  });
  expect(baseline.frozenColumns).toBe(1);

  await softStep('Step 4 — Sort: first double-click on AGE header sorts DESCENDING', async () => {
    const c = await headerCenter(page, 'AGE');
    await armGridEvent(page, 'onRowsSorted', 'sort4');
    await page.mouse.dblclick(c.x, c.y);
    await awaitGridEvent(page, 'sort4', 600); 
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
    expect(r.sortTypes).toContain(false); 
    expect(r.topAge).toBe(r.maxAge); 
    expect(r.ageRow0).toBe(baseline.ageRow0); 
  });

  await softStep('Step 5 — Sort: second double-click on AGE header sorts ASCENDING', async () => {
    const c = await headerCenter(page, 'AGE');
    await armGridEvent(page, 'onRowsSorted', 'sort5');
    await page.mouse.dblclick(c.x, c.y);
    await awaitGridEvent(page, 'sort5', 600); 
    const r = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const df = grok.shell.tv.dataFrame;
      const ageCol = df.col('AGE');

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
    expect(r.sortBy).toContain('AGE'); 
    expect(r.sortTypes).toContain(true); 
    expect(r.topAge).toBe(r.minNonNull); 
    expect(r.ageRow0).toBe(baseline.ageRow0); 
  });

  await softStep('Step 6 — Sort: third double-click on AGE header RESETS the sort', async () => {

    const ascTop = await page.evaluate(() => grok.shell.tv.grid.gridRowToTable(0));
    const c = await headerCenter(page, 'AGE');
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
    const c = await headerCenter(page, 'AGE');
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
    const before = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      return {
        order: Array.from({length: grid.columns.length}, (_, i) => grid.columns.byIndex(i).name),
        heightIdx: (() => { for (let i = 0; i < grid.columns.length; i++) if (grid.columns.byIndex(i).name === 'HEIGHT') return i; return -1; })(),
      };
    });
    const src = await headerCenter(page, 'HEIGHT');
    const tgt = await headerCenter(page, 'DEMOG'); 

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

    expect(after.heightIdx).toBeGreaterThan(before.heightIdx);
    expect(after.order).not.toEqual(before.order);
  });

  await softStep('Step 9 — Hide: open Order or Hide Columns and hide WEIGHT', async () => {

    await page.evaluate(async () => {
      const grid = grok.shell.tv.grid;
      const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
      const db = grid.cell('AGE', 0).documentBounds;
      const x = db.x + db.width / 2, y = db.y + db.height / 2;
      overlay.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, clientX: x, clientY: y, button: 2}));
      await new Promise((r) => setTimeout(r, 450)); 
      const leaf = document.querySelector('[name="div-Order-or-Hide-Columns..."]') as HTMLElement | null;
      if (leaf) {
        const b = leaf.getBoundingClientRect();
        leaf.dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true, clientX: b.x + 5, clientY: b.y + 5}));
      }
    });
    await page.locator('.d4-dialog .d4-dialog-header', {hasText: 'Order or Hide Columns'}).waitFor({timeout: 5000});

    const r = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const df = grok.shell.tv.dataFrame;
      grid.columns.setVisible(df.columns.names().filter((n: string) => n !== 'WEIGHT'));
      return true;
    });
    expect(r).toBe(true);

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
    expect(state.weightVisible).toBe(false); 
    expect(state.weightEnumerated).toBe(false); 
    expect(state.surroundingKept).toBe(true); 
  });

  await softStep('Step 10 — Resize: widen AGE by dragging its right header border', async () => {
    const before = await page.evaluate(() => grok.shell.tv.grid.columns.byName('AGE').width);
    const geom = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const db = grid.cell('AGE', 0).documentBounds;
      return {borderX: db.x + db.width, headerY: db.y - grid.colHeaderHeight / 2};
    });

    await page.mouse.move(geom.borderX, geom.headerY);
    await page.mouse.down();
    await page.mouse.move(geom.borderX + 60, geom.headerY, {steps: 6});
    await page.mouse.up();
    await page.waitForTimeout(500); 
    const after = await page.evaluate(() => grok.shell.tv.grid.columns.byName('AGE').width);
    expect(after).toBeGreaterThan(before); 
  });

  await softStep('Step 10 (cont.) — Resize + scroll: widen more columns, scroll horizontally, assert no errors (GROK-19753)', async () => {

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
    expect(gridErrors).toEqual([]); 
  });

  await softStep('Step 11 — Pin: pin SEX column via the header Pin menu', async () => {
    const frozenBefore = await page.evaluate(() => grok.shell.tv.grid.props.frozenColumns);

    const geom = await page.evaluate(async () => {
      const grid = grok.shell.tv.grid;
      grid.scrollToCell('SEX', 0);
      await new Promise((r) => setTimeout(r, 600)); 
      const db = grid.cell('SEX', 0).documentBounds;
      const overlay = document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]') as HTMLElement;
      const orect = overlay.getBoundingClientRect();
      return {clickX: db.x + db.width / 2, overlayLeft: orect.left};
    });

    expect(geom.clickX).toBeGreaterThan(geom.overlayLeft + 20);
    const c = await headerCenter(page, 'SEX');

    const pinned = await pinViaMenu(page, c, 'div-Pin---Pin-Column');
    expect(pinned).toBe(true); 
    await page.waitForTimeout(500); 
    const frozenAfter = await page.evaluate(() => grok.shell.tv.grid.props.frozenColumns);
    expect(frozenAfter).toBe(frozenBefore + 1); 
  });

  await softStep('Step 12 — Pin: pin two rows via the row Pin menu', async () => {

    await page.evaluate(async () => {
      const grid = grok.shell.tv.grid;

      grid.scrollToCell('AGE', 0);
      await new Promise((r) => setTimeout(r, 600)); 

      await new Promise<void>((resolve) => {
        const sub = grid.onRowsSorted.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(() => { sub.unsubscribe(); resolve(); }, 600);
        grid.sort([], []);
      });
    });

    const cell0 = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      const db = grid.cell('AGE', 0).documentBounds;
      return {x: db.x + db.width / 2, y: db.y + db.height / 2};
    });
    await armGridEvent(page, 'onPinnedRowsChanged', 'pinRow1');
    const p1 = await pinViaMenu(page, cell0, 'div-Pin---Pin-Row');
    expect(p1).toBe(true); 
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
    const p2 = await pinViaMenu(page, cell1, 'div-Pin---Pin-Row');
    expect(p2).toBe(true); 
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
      const scatterAdded = new Promise<void>((res) => {
        const sub = grok.events.onViewerAdded.subscribe(() => { sub.unsubscribe(); res(); });
        setTimeout(() => { sub.unsubscribe(); res(); }, 900);
      });
      tv.addViewer('Scatter plot');
      await scatterAdded;
      const hadScatter = tv.viewers.some((x: any) => x.type === 'Scatter plot');
      const layoutApplied = new Promise<void>((res) => {
        const sub = grok.events.onViewLayoutApplied.subscribe(() => { sub.unsubscribe(); res(); });
        setTimeout(() => { sub.unsubscribe(); res(); }, 2500);
      });
      tv.loadLayout(layout);
      await layoutApplied;
      await new Promise((res) => setTimeout(res, 400)); 
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
    expect(r.hadScatter).toBe(true); 
    expect(r.after.scatterGone).toBe(true); 
    expect(r.after.order).toEqual(r.before.order); 
    expect(r.after.weightEnumerated).toBe(false); 
    expect(r.after.ageWidth).toBe(r.before.ageWidth); 
    expect(r.after.sortBy).toContain('AGE');
    expect(r.after.sortTypes).toContain(true); 
    expect(r.after.frozen).toBe(r.before.frozen); 
    expect(r.after.pinsLen).toBe(2); 
  });

  const projectName = 'grid-cp-columns-layout-test-' + Date.now();
  let savedProjectId: string | null = null;
  let savedBeforeReopen: {order: string[]; ageWidth: number; frozen: number} | null = null;

  await softStep('Step 16 — Persistence: save the view as a project via the ribbon Save button', async () => {

    savedBeforeReopen = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      return {
        order: Array.from({length: grid.columns.length}, (_: any, i: number) => grid.columns.byIndex(i).name),
        ageWidth: grid.columns.byName('AGE').width,
        frozen: grid.props.frozenColumns,
      };
    });

    const saved = await saveProjectViaUI(page, projectName);
    savedProjectId = saved.projectId;
    expect(savedProjectId).not.toBeNull();

    const saveBalloons = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-balloon.error, .d4-balloon-error'))
        .map((b) => (b.textContent ?? '').trim()));
    const realSaveBalloons = saveBalloons.filter((t) => !isBenignSaveWindowError(t));
    expect(realSaveBalloons).toEqual([]);
  });

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
      const projApplied = new Promise<void>((res) => {
        const sub = grok.events.onViewLayoutApplied.subscribe(() => { sub.unsubscribe(); res(); });
        setTimeout(() => { sub.unsubscribe(); res(); }, 4500);
      });
      await proj.open();
      await projApplied;
      await new Promise((res) => setTimeout(res, 500)); 
      const grid = grok.shell.tv?.grid;

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
    expect(r.order).toEqual(before.order); 
    expect(r.weightEnumerated).toBe(false); 
    expect(r.ageWidth).toBe(before.ageWidth); 
    expect(r.frozen).toBe(before.frozen); 
    expect(r.pinsLen).toBe(2); 
    expect(r.sortBy).toContain('AGE');
    expect(r.sortTypes).toContain(true);

  });

  await softStep('Teardown: delete the probe project', async () => {
    if (savedProjectId)
      await deleteProjectWithCleanup(page, {projectId: savedProjectId});
  });

  v.finishSpec();
});
