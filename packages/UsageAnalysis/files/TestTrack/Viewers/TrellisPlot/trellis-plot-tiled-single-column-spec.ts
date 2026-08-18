/* ---
realizes: [trellisplot.cp.tiled-single-column, trellisplot.int.tiled-single-column-layout]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const splitColumn = 'DIS_POP';

const VIEWPORT_CLAMP = 5;
const DEFAULT_TILES_PER_ROW = 4; 
const untiledCells = (cats: number) => Math.min(VIEWPORT_CLAMP, cats);
const tiledPerRow = (cats: number, tilesPerRow: number) =>
  Math.min(VIEWPORT_CLAMP, Math.min(tilesPerRow, cats));
const tiledRows = (cats: number, tilesPerRow: number) =>
  Math.min(VIEWPORT_CLAMP, Math.ceil(cats / Math.min(tilesPerRow, cats)));
const tiledCells = (cats: number, tilesPerRow: number) =>
  tiledPerRow(cats, tilesPerRow) * tiledRows(cats, tilesPerRow);

const tiledCanvasCells = (cats: number, tilesPerRow: number) =>
  Math.min(cats, tiledCells(cats, tilesPerRow));

async function gridStructure(page: Page): Promise<{
  count: number; rows: number; maxPerRow: number; withCanvas: number;
}> {
  return page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    const cells = Array.from(root.querySelectorAll('.d4-trellis-plot-cell'));
    const byTop: Record<string, number> = {};
    for (const c of cells) {
      const t = String(Math.round(c.getBoundingClientRect().top));
      byTop[t] = (byTop[t] || 0) + 1;
    }
    const counts = Object.values(byTop);
    return {
      count: cells.length,
      rows: counts.length,
      maxPerRow: counts.length ? Math.max(...counts) : 0,
      withCanvas: cells.filter((c) => !!c.querySelector('canvas')).length,
    };
  });
}

async function controlPanelState(page: Page): Promise<{
  exists: boolean; display: string; showControlPanel: boolean;
}> {
  return page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement;
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    const vs = root.querySelector('[name="viewer selector"]') as HTMLElement | null;
    const cp = vs ? vs.parentElement : null;
    return {
      exists: !!cp,
      display: cp ? getComputedStyle(cp).display : '<node missing>',
      showControlPanel: tp.props.showControlPanel as boolean,
    };
  });
}

async function setTrellisProps(page: Page, props: Record<string, any>): Promise<void> {
  await page.evaluate((p) => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    for (const k of Object.keys(p)) tp.props[k] = p[k];
  }, props);
  await v.waitForViewerRendered(page, 'Trellis plot', 900);
}

async function waitForControlPanelDisplay(
  page: Page, wantNone: boolean, capMs = 2800,
): Promise<void> {
  await page.waitForFunction((none) => {
    const root = document.querySelector('[name="viewer-Trellis-plot"]') as HTMLElement | null;
    const vs = root?.querySelector('[name="viewer selector"]') as HTMLElement | null;
    const cp = vs ? vs.parentElement : null;
    if (!cp) return false;
    return (getComputedStyle(cp).display === 'none') === none;
  }, wantNone, {timeout: capMs}).catch(() => {});
}

test('Trellis plot: tiled single-column geometry and auto-layout control-panel coupling', async ({page}) => {
  test.setTimeout(600_000);
  page.setDefaultTimeout(120_000);

  await loginToDatagrok(page);
  await page.waitForTimeout(5000); 

  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch {}
    try { grok.shell.windows.simpleMode = true; } catch {}
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    await new Promise((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(resolve, 3000);
    });
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  const splitCats = await page.evaluate(
    (name) => grok.shell.tv.dataFrame.col(name).categories.length, splitColumn);

  expect(splitCats).toBeGreaterThan(DEFAULT_TILES_PER_ROW);

  await v.addViewerByIcon(page, 'trellis-plot', 'Trellis-plot', 15000);
  await setTrellisProps(page, {xColumnNames: [splitColumn], yColumnNames: []});

  const oneColumnOnly = await page.evaluate(() => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    return tp.oneColumnOnly as boolean;
  });
  expect(oneColumnOnly).toBe(true);

  const gear = page.locator(
    '.panel-base:has([name="viewer-Trellis-plot"]) .panel-titlebar [name="icon-font-icon-settings"]',
  ).first();
  await gear.click();
  await page.waitForTimeout(1200); 

  const tiledAtEntry = await page.evaluate(() => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    return tp.props.useTiledView as boolean;
  });
  expect(tiledAtEntry).toBe(true);
  await setTrellisProps(page, {useTiledView: false});

  const untiledEntry = await gridStructure(page);
  expect(untiledEntry.rows).toBe(1);
  expect(untiledEntry.count).toBe(untiledCells(splitCats));

  await softStep('Scenario 1 — Step 4', async () => {

    await setTrellisProps(page, {useTiledView: true});
    const s = await gridStructure(page);
    expect(s.rows).toBe(tiledRows(splitCats, DEFAULT_TILES_PER_ROW));
    expect(s.maxPerRow).toBe(tiledPerRow(splitCats, DEFAULT_TILES_PER_ROW));
    expect(s.count).toBe(tiledCells(splitCats, DEFAULT_TILES_PER_ROW));
    expect(s.rows).toBeGreaterThan(untiledEntry.rows);
    expect(s.count).not.toBe(untiledEntry.count);

    expect(s.withCanvas).toBe(tiledCanvasCells(splitCats, DEFAULT_TILES_PER_ROW));
  });

  await softStep('Scenario 1 — Step 5', async () => {

    await setTrellisProps(page, {tilesPerRow: DEFAULT_TILES_PER_ROW});
    const atDefault = await gridStructure(page);
    expect(atDefault.maxPerRow).toBe(tiledPerRow(splitCats, DEFAULT_TILES_PER_ROW));
    expect(atDefault.count).toBe(tiledCells(splitCats, DEFAULT_TILES_PER_ROW));

    expect(atDefault.maxPerRow).toBeGreaterThan(1);

    await setTrellisProps(page, {tilesPerRow: 1});
    const atOne = await gridStructure(page);
    expect(atOne.maxPerRow).toBe(tiledPerRow(splitCats, 1));
    expect(atOne.rows).toBe(tiledRows(splitCats, 1));
    expect(atOne.count).toBe(tiledCells(splitCats, 1));

    expect(atOne.count).toBeLessThan(atDefault.count);
  });

  await softStep('Scenario 1 — Step 6', async () => {

    await setTrellisProps(page, {tilesPerRow: 3});
    const s = await gridStructure(page);
    expect(s.maxPerRow).toBe(tiledPerRow(splitCats, 3));
    expect(s.rows).toBe(tiledRows(splitCats, 3));
    expect(s.count).toBe(tiledCells(splitCats, 3));
    expect(s.maxPerRow).toBeLessThan(tiledPerRow(splitCats, DEFAULT_TILES_PER_ROW));
  });

  await softStep('Scenario 1 — Step 7', async () => {

    await setTrellisProps(page, {tilesPerRow: DEFAULT_TILES_PER_ROW});
    const tiled = await gridStructure(page);
    expect(tiled.rows).toBe(tiledRows(splitCats, DEFAULT_TILES_PER_ROW));
    await setTrellisProps(page, {useTiledView: false});
    const off = await gridStructure(page);
    expect(off.rows).toBe(1);
    expect(off.maxPerRow).toBe(untiledCells(splitCats));
    expect(off.count).toBe(untiledCells(splitCats));

    expect(off.rows).toBeLessThan(tiled.rows);
    expect(off.count).not.toBe(tiled.count);

    await setTrellisProps(page, {useTiledView: true});
  });

  await softStep('Scenario 2 — Step 3', async () => {

    await setTrellisProps(page, {autoLayout: true, showControlPanel: true});
    await page.setViewportSize({width: 1920, height: 1080});
    await waitForControlPanelDisplay(page, false);
    const s = await controlPanelState(page);
    expect(s.exists).toBe(true);
    expect(s.display).not.toBe('none');
  });

  await softStep('Scenario 2 — Step 4', async () => {

    const before = await controlPanelState(page);
    expect(before.exists).toBe(true);
    await page.setViewportSize({width: 900, height: 280});
    await waitForControlPanelDisplay(page, true);
    const shrunk = await controlPanelState(page);

    expect(shrunk.exists).toBe(true);
    expect(shrunk.display).toBe('none');
    expect(shrunk.showControlPanel).toBe(before.showControlPanel);
  });

  await softStep('Scenario 2 — Step 5', async () => {

    const before = await controlPanelState(page);
    await page.setViewportSize({width: 1920, height: 1080});
    await waitForControlPanelDisplay(page, false);
    const restored = await controlPanelState(page);

    expect(restored.exists).toBe(true);
    expect(restored.display).not.toBe('none');
    expect(restored.showControlPanel).toBe(before.showControlPanel);
  });

  await softStep('Scenario 2 — Step 6', async () => {

    await setTrellisProps(page, {autoLayout: false, showControlPanel: true});
    await page.setViewportSize({width: 900, height: 280});

    await waitForControlPanelDisplay(page, false);
    const explicitOn = await controlPanelState(page);
    expect(explicitOn.exists).toBe(true);
    expect(explicitOn.display).not.toBe('none');

    await setTrellisProps(page, {showControlPanel: false});
    const explicitOff = await controlPanelState(page);

    expect(explicitOff.exists).toBe(true);
    expect(explicitOff.display).toBe('none');

    await setTrellisProps(page, {showControlPanel: true, autoLayout: true});
    await page.setViewportSize({width: 1920, height: 1080});
    await waitForControlPanelDisplay(page, false);
  });

  v.finishSpec();
});
