/* ---
realizes: [grid.cp.appearance-summary-persist, grid.cp.columns-layout-persist, grid.cp.dialogs-groups]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaApi, deleteProjectWithCleanup} from '../../helpers/projects';
import * as g from './grid-helpers';

declare const grok: any;
declare const DG: any;

// The server lane of the Grid section: the steps whose subject is a layout or a project surviving
// a round-trip, plus the steps that need PowerGrid's summary-column renderers (the local client
// serves no packages). The state each round-trip persists is set up through the API; the gestures
// that produce it are proved by the local-lane specs.
test.use(specTestOptions);

async function reopenProject(page: Page, projectId: string): Promise<void> {
  await page.evaluate(async (id) => {
    const w = window as any;
    grok.shell.closeAll();
    await w.__poll(() => Array.from(grok.shell.tableViews).length, (n: number) => n === 0, 3000, 50);
    await (await grok.dapi.projects.find(id)).open();
  }, projectId);
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
  await page.evaluate(() => (window as any).__tableReady(5000));
}

const loadFailureBalloons = () => Array.from(document.querySelectorAll('.d4-balloon.error, .d4-balloon-error'))
  .map((b) => (b.textContent ?? '').trim()).filter((t) => /error loading/i.test(t));

const errorBalloons = (page: Page) => page.evaluate(() =>
  Array.from(document.querySelectorAll('.d4-balloon.error, .d4-balloon-error'))
    .map((b) => (b.textContent ?? '').trim()));

// One round-trip through dapi.layouts: save, corrupt the view with a foreign viewer, re-apply.
// `readState` is serialised into the page, so it must be a closure-free function of its arguments.
async function layoutRoundTrip<T>(page: Page, readState: (grid: any, df?: any) => T): Promise<{before: T; after: T; hadScatter: boolean; scatterGone: boolean}> {
  return page.evaluate(async (readSrc) => {
    const w = window as any;
    const read = new Function('grid', 'df', `return (${readSrc})(grid, df);`);
    const tv = grok.shell.tv;
    const before = read(tv.grid, tv.dataFrame);
    const layout = await grok.dapi.layouts.save(tv.saveLayout());
    await w.__settled('grok.events.onViewerAdded', () => tv.addViewer('Scatter plot'), 900);
    const hadScatter = tv.viewers.some((x: any) => x.type === 'Scatter plot');
    await w.__settled('grok.events.onViewLayoutApplied', () => tv.loadLayout(layout), 2800);
    await w.__settledFor(() => JSON.stringify(read(grok.shell.tv.grid, grok.shell.tv.dataFrame)), 200, 1500);
    const after = read(grok.shell.tv.grid, grok.shell.tv.dataFrame);
    const scatterGone = !grok.shell.tv.viewers.some((x: any) => x.type === 'Scatter plot');
    await grok.dapi.layouts.delete(layout);
    return {before, after, hadScatter, scatterGone};
  }, readState.toString());
}

interface AppearanceRows { nullHeightRow: number; minAgeRow: number; maxAgeRow: number; }

// The appearance battery: everything the persistence steps re-read after a round-trip.
const appearanceState = (grid: any, df: any) => {
  const hc = df.col('HEIGHT'); const ac = df.col('AGE'); const scol = df.col('SEX');
  let lowRow = -1; let highRow = -1; let midRow = -1; let nullHeightRow = -1;
  for (let i = 0; i < df.rowCount; i++) {
    if (hc.isNone(i)) { if (nullHeightRow < 0) nullHeightRow = i; continue; }
    const val = hc.get(i);
    if (lowRow < 0 && val < 160) lowRow = i;
    if (highRow < 0 && val > 180) highRow = i;
    if (midRow < 0 && val >= 160 && val <= 180) midRow = i;
  }
  let minAgeRow = -1; let maxAgeRow = -1; let minV = Infinity; let maxV = -Infinity;
  for (let i = 0; i < df.rowCount; i++) {
    if (ac.isNone(i)) continue;
    const val = ac.get(i);
    if (val < minV) { minV = val; minAgeRow = i; }
    if (val > maxV) { maxV = val; maxAgeRow = i; }
  }
  let mRow = -1; let fRow = -1;
  for (let i = 0; i < df.rowCount && (mRow < 0 || fRow < 0); i++) {
    const val = scol.get(i);
    if (val === 'M' && mRow < 0) mRow = i;
    if (val === 'F' && fRow < 0) fRow = i;
  }
  return {
    cols: grid.columns.length as number,
    ageCC: df.col('AGE').getTag('.color-coding-type') as string,
    heightCC: df.col('HEIGHT').getTag('.color-coding-type') as string,
    heightCond: df.col('HEIGHT').getTag('.color-coding-conditional') as string,
    sexCC: df.col('SEX').getTag('.color-coding-type') as string,
    weightCC: df.col('WEIGHT').getTag('.color-coding-type') as string,
    weightSrc: df.col('WEIGHT').getTag('.%color-coding-linked-column-name') as string,
    ageBounds: grid.cell('AGE', 0).bounds.height as number,
    nullColor: (grid.cell('HEIGHT', nullHeightRow).color >>> 0) as number,
    cellTypes: Array.from({length: grid.columns.length}, (_: any, i: number) => grid.columns.byIndex(i).cellType) as string[],
    ageMinColor: (grid.cell('AGE', minAgeRow).color >>> 0) as number,
    ageMaxColor: (grid.cell('AGE', maxAgeRow).color >>> 0) as number,
    agePlainColor: (grid.cell('DEMOG', minAgeRow).color >>> 0) as number,
    heightLowColor: (grid.cell('HEIGHT', lowRow).color >>> 0) as number,
    heightHighColor: (grid.cell('HEIGHT', highRow).color >>> 0) as number,
    heightMidColor: (midRow >= 0 ? grid.cell('HEIGHT', midRow).color >>> 0 : -1) as number,
    sexMColor: (grid.cell('SEX', mRow).color >>> 0) as number,
    sexFColor: (grid.cell('SEX', fRow).color >>> 0) as number,
    weightMEqualsSexM: (grid.cell('WEIGHT', mRow).color === grid.cell('SEX', mRow).color) as boolean,
    weightFEqualsSexF: (grid.cell('WEIGHT', fRow).color === grid.cell('SEX', fRow).color) as boolean,
  };
};
type AppearanceState = ReturnType<typeof appearanceState>;

const readAppearance = (page: Page) => page.evaluate((src) =>
  new Function('grid', 'df', `return (${src})(grid, df);`)(grok.shell.tv.grid, grok.shell.tv.dataFrame) as AppearanceState,
appearanceState.toString());

const appearanceHolds = (s: AppearanceState) => s.ageCC === 'Linear' && s.ageMinColor !== s.ageMaxColor &&
  s.heightCC === 'Conditional' && s.heightLowColor === 0xff0000ff && s.heightHighColor === 0xffff0000 &&
  s.sexCC === 'Categorical' && s.sexMColor !== s.sexFColor && s.weightCC === 'Linked' &&
  s.weightMEqualsSexM && s.weightFEqualsSexF && s.nullColor === 0xffffaaaa;

// the summary columns raise the row height, so the round-trips compare against the height read after them
function expectAppearance(s: AppearanceState, ageBounds: number): void {
  expect(s.ageCC).toBe('Linear');
  expect(s.ageMinColor).not.toBe(s.ageMaxColor);
  expect(s.ageMinColor).not.toBe(s.agePlainColor);
  expect(s.ageMaxColor).not.toBe(s.agePlainColor);
  expect(s.heightCC).toBe('Conditional');
  expect(s.heightLowColor).toBe(0xff0000ff);
  expect(s.heightHighColor).toBe(0xffff0000);
  expect(s.heightMidColor).not.toBe(0xff0000ff);
  expect(s.heightMidColor).not.toBe(0xffff0000);
  expect(s.sexCC).toBe('Categorical');
  expect(s.sexMColor).not.toBe(s.sexFColor);
  expect(s.weightCC).toBe('Linked');
  expect(s.weightSrc).toBe('SEX');
  expect(s.weightMEqualsSexM).toBe(true);
  expect(s.weightFEqualsSexF).toBe(true);
  expect(s.ageBounds).toBe(ageBounds);
  expect(s.nullColor).toBe(0xffffaaaa);
}

const summaryTypes: {leaf: string; cellType: string}[] = [
  {leaf: 'Sparklines', cellType: 'sparkline'},
  {leaf: 'Bar-Chart', cellType: 'barchart'},
  {leaf: 'Pie-Chart', cellType: 'piechart'},
  {leaf: 'Radar', cellType: 'radar'},
  {leaf: 'Smart-Form', cellType: 'smartform'},
  {leaf: 'Tags', cellType: 'tags'},
  {leaf: 'Confidence-Interval', cellType: 'confidenceinterval'},
];

test('Grid — appearance, summary columns and stats rows persist across layout and project round-trips', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);
  const flags = await g.readShellFlags(page);
  const errors = g.trackErrors(page, g.BENIGN_NOISE);
  let projectId: string | null = null;
  let colsAfterSummary = 0;
  try {
    await v.openTable(page, {path: g.DEMOG, semTypeTimeoutMs: 3000});

    await softStep('Setup — Steps 4-10 applied through the API: four colour codings, row height 48, missing-value colour', async () => {
      await page.evaluate((cond) => {
        const df = grok.shell.tv.dataFrame; const grid = grok.shell.tv.grid;
        df.col('AGE').meta.colors.setLinear();
        df.col('HEIGHT').meta.colors.setConditional(cond);
        df.col('SEX').meta.colors.setCategorical();
        df.col('WEIGHT').tags['.color-coding-type'] = 'Linked';
        df.col('WEIGHT').tags['.%color-coding-linked-column-name'] = 'SEX';
        grid.props.rowHeight = 48;
        grid.props.missingValueColor = DG.Color.fromHtml('#FFAAAA');
        grid.invalidate();
      }, {'<160': '#0000FF', '>180': '#FF0000'});
      expectAppearance(await v.pollValue(() => readAppearance(page), (s) => appearanceHolds(s) && s.ageBounds === 48, 3000, 100), 48);
    });

    await softStep('Step 12 — Summary columns: add all seven one-click types; column count and cellType track each add', async () => {
      const before = await page.evaluate(() => grok.shell.tv.grid.columns.length);
      let count = before;
      const results: {leaf: string; added: boolean; cellType: string}[] = [];
      for (const t of summaryTypes) {
        const added = await g.clickMenuLeaf(page, await g.cellCenter(page, 'AGE', 0),
          ['div-Add', 'div-Add---Summary-Columns'], `div-Add---Summary-Columns---${t.leaf}`);
        const expectedLen = count + 1;
        const state = await page.evaluate(async ({want, cellType}) => {
          const read = () => {
            const grid = grok.shell.tv.grid;
            return {len: grid.columns.length, lastType: grid.columns.byIndex(grid.columns.length - 1).cellType};
          };
          return (window as any).__poll(read, (x: any) => x.len === want && x.lastType === cellType, 2000, 25);
        }, {want: expectedLen, cellType: t.cellType});
        results.push({leaf: t.leaf, added, cellType: state.lastType});
        if (added && state.len === count + 1) count = state.len;
      }
      colsAfterSummary = count;

      for (let i = 0; i < summaryTypes.length; i++) {
        expect(results[i].added).toBe(true);
        expect(results[i].cellType).toBe(summaryTypes[i].cellType);
      }
      expect(colsAfterSummary).toBe(before + summaryTypes.length);
    });

    await softStep('Step 13 — Stats rows: add min and max; summary columns survive and no console error (GROK-19809)', async () => {
      const errBefore = errors.count();
      for (const stat of ['min', 'max']) {
        expect(await g.clickMenuLeaf(page, await g.cellCenter(page, 'AGE', 0),
          ['div-Add', 'div-Add---Column-Stats'], `div-Add---Column-Stats---${stat}`)).toBe(true);
        // the error window after each stats row is the assertion; the grid's own paint closes it
        await v.waitForGridPainted(page, {gapMs: 120, capMs: 500});
      }
      const colsAfter = await page.evaluate(() => grok.shell.tv.grid.columns.length);
      expect(colsAfter).toBe(colsAfterSummary);
      expect(errors.list.slice(errBefore)).toEqual([]);
    });

    await softStep('Step 15 — Persistence: save layout, add a foreign viewer, re-apply the layout (GROK-19769)', async () => {
      const r = await layoutRoundTrip(page, appearanceState);
      expect(r.hadScatter).toBe(true);
      expect(r.scatterGone).toBe(true);
      expectAppearance(r.after, r.before.ageBounds);
      expect(r.after.heightCond).toBe(r.before.heightCond);
      expect(r.after.nullColor).toBe(r.before.nullColor);
      expect(r.after.cols).toBe(r.before.cols);
      expect(r.after.cellTypes).toEqual(r.before.cellTypes);
    });

    let peakCellHeight = 0;
    await softStep('Step 16 — Persistence: save the view as a project', async () => {
      peakCellHeight = (await readAppearance(page)).ageBounds;
      projectId = (await saveProjectViaApi(page, 'grid-appearance-summary-persist-' + Date.now())).projectId;
      expect(projectId).not.toBeNull();
      expect((await errorBalloons(page)).filter((t) => !g.BENIGN_NOISE(t))).toEqual([]);
    });

    await softStep('Step 17 — Persistence: Close All and reopen the project; the full battery holds, error delta 0', async () => {
      const errBefore = errors.count();
      await reopenProject(page, projectId!);
      const s = await v.pollValue(() => readAppearance(page), appearanceHolds, 5000, 100);
      const balloons = await page.evaluate(loadFailureBalloons);
      expect(balloons).toEqual([]);
      expect(errors.list.slice(errBefore)).toEqual([]);
      expectAppearance(s, peakCellHeight);
      expect(s.cols).toBe(colsAfterSummary);
      expect(s.cellTypes).toEqual(expect.arrayContaining(summaryTypes.map((t) => t.cellType)));
    });
  } finally {
    if (projectId) await deleteProjectWithCleanup(page, {projectId});
    errors.stop();
    await g.leaveShellClean(page, flags);
  }
  v.finishSpec();
});

const geometryState = (grid: any) => ({
  order: Array.from({length: grid.columns.length}, (_: any, i: number) => grid.columns.byIndex(i).name) as string[],
  weightEnumerated: (() => {
    for (let i = 0; i < grid.columns.length; i++) { const c = grid.columns.byIndex(i); if (c.name === 'WEIGHT' && c.visible) return true; }
    return false;
  })() as boolean,
  ageWidth: grid.columns.byName('AGE').width as number,
  sortBy: grid.props.sortByColumnNames.slice() as string[],
  sortTypes: grid.props.sortTypes.slice() as boolean[],
  frozen: grid.props.frozenColumns as number,
  pinsLen: Array.from(grid.pinnedRows).length as number,
});
type GeometryState = ReturnType<typeof geometryState>;

const readGeometry = (page: Page) => page.evaluate((src) =>
  new Function('grid', `return (${src})(grid);`)(grok.shell.tv.grid) as GeometryState, geometryState.toString());

test('Grid — column geometry persists across layout and project round-trips', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);
  const flags = await g.readShellFlags(page);
  const errors = g.trackErrors(page, g.BENIGN_NOISE);
  let projectId: string | null = null;
  try {
    await v.openTable(page, {path: g.DEMOG, semTypeTimeoutMs: 3000});

    await softStep('Setup — Steps 4-12 applied through the API and the Pin menu: HEIGHT moved, WEIGHT hidden, AGE widened, SEX pinned, two rows pinned, sorted by AGE', async () => {
      await page.evaluate(() => {
        const grid = grok.shell.tv.grid;
        const df = grok.shell.tv.dataFrame;
        const names = df.columns.names() as string[];
        const order = names.filter((n) => n !== 'HEIGHT');
        order.splice(order.indexOf('DEMOG') + 1, 0, 'HEIGHT');
        grid.columns.setOrder(order);
        grid.columns.setVisible(names.filter((n) => n !== 'WEIGHT'));
        grid.columns.byName('AGE').width = 160;
      });
      const frozenBefore = await page.evaluate(() => grok.shell.tv.grid.props.frozenColumns);
      expect(await g.pinViaMenu(page, await g.headerCenter(page, 'SEX'), 'div-Pin---Pin-Column')).toBe(true);
      expect(await v.pollValue(() => page.evaluate(() => grok.shell.tv.grid.props.frozenColumns),
        (f) => f === frozenBefore + 1, 1000, 50)).toBe(frozenBefore + 1);

      // pinned rows are shown again at the top, so the next pin must target a grid row whose table row is not pinned yet
      for (const n of [1, 2]) {
        const cell = await page.evaluate(() => {
          const grid = grok.shell.tv.grid;
          const pinned = Array.from(grid.pinnedRows) as number[];
          let gridRow = 0;
          for (let gr = 0; gr < 12; gr++) if (!pinned.includes(grid.gridRowToTable(gr))) { gridRow = gr; break; }
          const db = grid.cell('AGE', gridRow).documentBounds;
          return {x: db.x + db.width / 2, y: db.y + db.height / 2};
        });
        expect(await g.pinViaMenu(page, cell, 'div-Pin---Pin-Row')).toBe(true);
        await page.waitForFunction((k) => Array.from(grok.shell.tv.grid.pinnedRows).length >= k, n, {timeout: 3000});
      }
      await page.evaluate(() => new Promise<void>((resolve) => {
        const grid = grok.shell.tv.grid;
        const sub = grid.onRowsSorted.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(() => { sub.unsubscribe(); resolve(); }, 600);
        grid.sort(['AGE'], [true]);
      }));
      const s = await readGeometry(page);
      expect(s.weightEnumerated).toBe(false);
      expect(s.order.indexOf('HEIGHT')).toBeGreaterThan(s.order.indexOf('DEMOG'));
      expect(s.ageWidth).toBe(160);
      expect(s.frozen).toBe(frozenBefore + 1);
      expect(s.pinsLen).toBe(2);
      expect(s.sortBy).toContain('AGE');
      expect(s.sortTypes).toContain(true);
    });

    await softStep('Step 14 — Persistence: save layout, add a foreign viewer, re-apply the layout', async () => {
      const r = await layoutRoundTrip(page, geometryState);
      expect(r.hadScatter).toBe(true);
      expect(r.scatterGone).toBe(true);
      expect(r.after.order).toEqual(r.before.order);
      expect(r.after.weightEnumerated).toBe(false);
      expect(r.after.ageWidth).toBe(r.before.ageWidth);
      expect(r.after.sortBy).toContain('AGE');
      expect(r.after.sortTypes).toContain(true);
      expect(r.after.frozen).toBe(r.before.frozen);
      expect(r.after.pinsLen).toBe(2);
    });

    await softStep('Step 16 — Persistence: save the view as a project and reopen it clean', async () => {
      const before = await readGeometry(page);
      projectId = (await saveProjectViaApi(page, 'grid-cp-columns-layout-test-' + Date.now())).projectId;
      expect(projectId).not.toBeNull();
      expect((await errorBalloons(page)).filter((t) => !g.BENIGN_NOISE(t))).toEqual([]);

      await reopenProject(page, projectId);
      const r = await v.pollValue(() => readGeometry(page),
        (s) => JSON.stringify(s.order) === JSON.stringify(before.order) && s.pinsLen === 2, 5000, 100);
      expect(await page.evaluate(loadFailureBalloons)).toEqual([]);
      expect(r.order).toEqual(before.order);
      expect(r.weightEnumerated).toBe(false);
      expect(r.ageWidth).toBe(before.ageWidth);
      expect(r.frozen).toBe(before.frozen);
      expect(r.pinsLen).toBe(2);
      expect(r.sortBy).toContain('AGE');
      expect(r.sortTypes).toContain(true);
    });
  } finally {
    if (projectId) await deleteProjectWithCleanup(page, {projectId});
    errors.stop();
    await g.leaveShellClean(page, flags);
  }
  v.finishSpec();
});

test('Grid — column groups persist across a project round-trip (GROK-17441)', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);
  const flags = await g.readShellFlags(page);
  const errors = g.trackErrors(page, g.BENIGN_NOISE);
  let projectId: string | null = null;
  try {
    await v.openTable(page, {path: g.DEMOG, semTypeTimeoutMs: 3000});

    // what the Group columns dialog writes (xamgle column_commands.dart): the per-column group tag
    // and the .columnGroups map on the dataframe
    const groupsBeforeSave = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const groups: Record<string, any> = {
        AgeHeight: {columns: ['AGE', 'HEIGHT'], description: 'Description', color: '#1f77b4'},
        WeightSex: {columns: ['WEIGHT', 'SEX'], description: 'Description', color: '#2ca02c'},
      };
      for (const name of Object.keys(groups))
        for (const c of groups[name].columns) df.col(c).setTag('group', name);
      df.setTag('.columnGroups', JSON.stringify(groups));
      return df.getTag('.columnGroups');
    });
    expect(groupsBeforeSave).toContain('#1f77b4');
    expect(groupsBeforeSave).toContain('#2ca02c');

    await softStep('Step 33 — Save the view as a project via the ribbon Save button', async () => {
      // saveProjectViaApi (uploadDataFrame + tables.save) keeps column tags but drops the dataframe-level
      // .columnGroups tag, measured 2026-09-03; the ribbon save keeps it
      projectId = await g.saveProjectViaRibbon(page, 'grid-dialogs-groups-' + Date.now());
      expect(projectId).not.toBeNull();
      expect((await errorBalloons(page)).filter((t) => !g.BENIGN_NOISE(t))).toEqual([]);
    });

    await softStep('Step 36 — Reopen the project: group colours intact and console-error delta 0 (GROK-17441)', async () => {
      const errBefore = errors.count();
      await reopenProject(page, projectId!);
      const r = await v.pollValue(() => page.evaluate(() => {
        const df = grok.shell.tv?.dataFrame;
        let tagKeys: string[] = [];
        try { tagKeys = Object.keys(df.tags); } catch (_) { tagKeys = ['<unreadable>']; }
        return {
          reopened: !!df,
          columnGroups: df ? df.getTag('.columnGroups') : null,
          ageGroup: df ? df.col('AGE').getTag('group') : null,
          weightGroup: df ? df.col('WEIGHT').getTag('group') : null,
          tagKeys,
        };
      }), (x) => x.reopened && !!x.columnGroups, 5000, 100);
      expect(await page.evaluate(loadFailureBalloons)).toEqual([]);
      expect(r.reopened).toBe(true);
      expect(r.columnGroups, `table tags after reopen: ${r.tagKeys.join(',')}; AGE group=${r.ageGroup}`).toBe(groupsBeforeSave);
      expect(r.columnGroups).toContain('#1f77b4');
      expect(r.columnGroups).toContain('#2ca02c');
      expect(r.ageGroup).toBeTruthy();
      expect(r.weightGroup).toBeTruthy();
      expect(errors.list.slice(errBefore)).toEqual([]);
    });
  } finally {
    if (projectId) await deleteProjectWithCleanup(page, {projectId});
    errors.stop();
    await g.leaveShellClean(page, flags);
  }
  v.finishSpec();
});

test('Grid — summary columns: remove via the top panel and survive a source-column removal', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);
  const flags = await g.readShellFlags(page);
  const errors = g.trackErrors(page, g.BENIGN_NOISE);
  try {
    await softStep('GROK-18256: remove a summary column via the top-panel remove icon', async () => {
      await v.openTable(page, {path: g.DEMOG, semTypeTimeoutMs: 3000});
      const before = await page.evaluate(() => grok.shell.tv.grid.columns.length);

      expect(await g.clickMenuLeaf(page, await g.cellCenter(page, 'AGE', 3),
        ['div-Add', 'div-Add---Summary-Columns'], 'div-Add---Summary-Columns---Sparklines')).toBe(true);
      const afterAdd = await v.pollValue(() => page.evaluate(() => grok.shell.tv.grid.columns.length), (n) => n === before + 1, 2000, 50);
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
      await page.evaluate(() => {
        const icon = document.querySelector('[name="icon-remove-selected-columns"]') as HTMLElement;
        icon.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, view: window} as any));
        icon.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, view: window} as any));
        icon.click();
      });
      await page.waitForFunction((n) => grok.shell.tv.grid.columns.length === n, before, {timeout: 5000}).catch(() => {});
      const afterRemove = await page.evaluate(() => grok.shell.tv.grid.columns.length);
      expect(afterRemove).toBe(before);
      const summaryGone = await page.evaluate((name) => !grok.shell.tv.grid.columns.byName(name), summaryName);
      expect(summaryGone).toBe(true);
    });

    await softStep('GROK-19942: grid keeps rendering after the CONTROL source column is removed', async () => {
      await v.closeAllAndWait(page);
      await v.openTable(page, {path: g.DEMOG, semTypeTimeoutMs: 3000});
      const before = await page.evaluate(() => grok.shell.tv.grid.columns.length);
      expect(await g.clickMenuLeaf(page, await g.cellCenter(page, 'AGE', 3),
        ['div-Add', 'div-Add---Summary-Columns'], 'div-Add---Summary-Columns---Tags')).toBe(true);
      await v.pollValue(() => page.evaluate(() => grok.shell.tv.grid.columns.length), (n) => n === before + 1, 2000, 50);
      const errBefore = errors.count();

      await page.evaluate(() => { grok.shell.tv.dataFrame.columns.remove('CONTROL'); });
      await v.waitForViewerRendered(page, 'Grid', 400);
      await g.focusGrid(page);
      await page.keyboard.press('PageDown');
      await v.waitForViewerRendered(page, 'Grid', 300);
      await page.keyboard.press('PageUp');
      await v.waitForViewerRendered(page, 'Grid', 300);

      const errsAfter = errors.list.slice(errBefore).filter((t) => /render|grid|null|argument/i.test(t)).length;
      expect(errsAfter).toBe(0);
      const gridStillPresent = await page.evaluate(() => !!document.querySelector('[name="viewer-Grid"] canvas[name="overlay"]'));
      expect(gridStillPresent).toBe(true);
    });
  } finally {
    errors.stop();
    await g.leaveShellClean(page, flags);
  }
  v.finishSpec();
});

test('Grid — a project with an extracted-rows table reopens without error (GROK-19717)', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);
  const flags = await g.readShellFlags(page);
  const errors = g.trackErrors(page, g.BENIGN_NOISE);
  let savedId: string | null = null;
  try {
    await v.openTable(page, {path: g.DEMOG, semTypeTimeoutMs: 3000});
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

    savedId = await page.evaluate(async (name) => {
      const w = window as any;
      const df = grok.shell.tv.dataFrame;
      let saved = null;
      for (let k = 0; k < 3 && !saved?.id; k++) {
        try {
          await grok.dapi.tables.uploadDataFrame(df);
          await w.__findSaved(async () => df.id && await grok.dapi.tables.find(df.id), 10000);
          const proj = DG.Project.create();
          proj.name = name;
          proj.addChild(df);
          saved = await grok.dapi.projects.save(proj);
          await w.__findSaved(() => grok.dapi.projects.find(saved.id), 2000);
        } catch (_) {
          await new Promise((r) => setTimeout(r, 700 * (k + 1)));
        }
      }
      return saved?.id ?? null;
    }, projectName);
    expect(savedId).toBeTruthy();

    const errBefore = errors.count();
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
    await page.evaluate(async (id) => { await (await grok.dapi.projects.find(id)).open(); }, savedId);

    await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
    await page.evaluate(() => (window as any).__tableReady(5000));
    await v.waitForGridPainted(page, {gapMs: 250, capMs: 1500});
    const balloons = await page.evaluate(() => {
      const obs = (window as any).__grokErrBalloonObs;
      if (obs) obs.disconnect();
      const captured = ((window as any).__grokErrBalloons || []).length;
      return captured + document.querySelectorAll('.d4-balloon.error').length;
    });
    expect(balloons).toBe(0);
    const gridPresent = await page.evaluate(() => !!document.querySelector('[name="viewer-Grid"]'));
    expect(gridPresent).toBe(true);

    const newErrs = errors.list.slice(errBefore).filter((t) => /null|argument|render/i.test(t)).length;
    expect(newErrs).toBe(0);
  } finally {
    if (savedId) await deleteProjectWithCleanup(page, {projectId: savedId});
    errors.stop();
    await g.leaveShellClean(page, flags);
  }
});
