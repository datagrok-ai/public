/* ---
realizes: [trellisplot.cp.tiled-single-column, trellisplot.int.tiled-single-column-layout]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
// The split column must be WIDER than the default Tiles Per Row, or tiled and untiled render
// the SAME grid: on the two-category SEX column both modes give one row of two
// [DOM 2026-08-11]. DIS_POP carries 6 categories on the live demog.
const splitColumn = 'DIS_POP';

// The cell grid is the scroll VIEWPORT, not the category list: oneColumnOnly clamps both ranges
// to five (trellis_plot_core.dart:855-858, getViewport() :887-889). Untiled = one row of
// min(5, cats); tiled = min(tilesPerRow, cats) wide x ceil(cats/width) rows (:79-95), same clamp. [DOM 2026-08-11]
const VIEWPORT_CLAMP = 5;
const DEFAULT_TILES_PER_ROW = 4; // trellis_plot_look.dart:70
const untiledCells = (cats: number) => Math.min(VIEWPORT_CLAMP, cats);
const tiledPerRow = (cats: number, tilesPerRow: number) =>
  Math.min(VIEWPORT_CLAMP, Math.min(tilesPerRow, cats));
const tiledRows = (cats: number, tilesPerRow: number) =>
  Math.min(VIEWPORT_CLAMP, Math.ceil(cats / Math.min(tilesPerRow, cats)));
const tiledCells = (cats: number, tilesPerRow: number) =>
  tiledPerRow(cats, tilesPerRow) * tiledRows(cats, tilesPerRow);
// A tiled grid is made of full rows, so the trailing slots of the last row are structural
// PADDING and carry no canvas [DOM 2026-08-11]. Every `count` assertion below counts SLOTS;
// where the rendered category count matters, the canvas-bearing cells are counted instead.
const tiledCanvasCells = (cats: number, tilesPerRow: number) =>
  Math.min(cats, tiledCells(cats, tilesPerRow));

// The grid has no DOM row container, so a "row" can only be inferred by grouping cells on
// their rounded top coordinate; maxPerRow is the widest such group. [DOM 2026-08-06]
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

// The control panel is the unnamed htmlDiv([viewerSelector.root, dataPropertiesDiv])
// (trellis_plot_core.dart:304), hence the parentElement hop. autoLayout only toggles its
// display (:922), so the two fields stay apart: a fused flag would pass on a deleted node. [DOM 2026-08-06]
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
  await page.evaluate(async (p) => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    for (const k of Object.keys(p)) tp.props[k] = p[k];
    await new Promise((r) => setTimeout(r, 2200));
  }, props);
}

test('Trellis plot: tiled single-column geometry and auto-layout control-panel coupling', async ({page}) => {
  test.setTimeout(600_000);
  page.setDefaultTimeout(120_000);

  await loginToDatagrok(page);
  await page.waitForTimeout(5000);

  // Setup: split by exactly ONE column (X = DIS_POP, Y empty) so oneColumnOnly engages.
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

  // Read off the live data, never hard-coded: every expected cell count, row count and
  // row width below is a clamp formula over it.
  const splitCats = await page.evaluate(
    (name) => grok.shell.tv.dataFrame.col(name).categories.length, splitColumn);
  // Precondition of the whole ladder: with cats <= Tiles Per Row the tiled grid degenerates
  // into the untiled band and Steps 4, 6 and 7 stop distinguishing the two modes.
  expect(splitCats).toBeGreaterThan(DEFAULT_TILES_PER_ROW);

  await v.addViewerByIcon(page, 'trellis-plot', 'Trellis-plot', 15000);
  await setTrellisProps(page, {xColumnNames: [splitColumn], yColumnNames: []});

  const oneColumnOnly = await page.evaluate(() => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    return tp.oneColumnOnly as boolean;
  });
  expect(oneColumnOnly).toBe(true);

  // Scenario 1 Step 2. The gear is not inside the viewer root but on the enclosing .panel-base
  // title bar, hence the :has() hop [DOM 2026-08-06]. Nothing below actuates through the panel;
  // this is the section's only exercise of that hop, so a title-bar change fails HERE.
  const gear = page.locator(
    '.panel-base:has([name="viewer-Trellis-plot"]) .panel-titlebar [name="icon-font-icon-settings"]',
  ).first();
  await gear.click();
  await page.waitForTimeout(1200);

  // Tiled View ships ON (trellis_plot_look.dart:68), so the ladder enters from the
  // non-default side: switched OFF here, Step 4 becomes a real off→on transition instead
  // of a no-op write of the value already in effect.
  const tiledAtEntry = await page.evaluate(() => {
    const tp = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Trellis plot') as any;
    return tp.props.useTiledView as boolean;
  });
  expect(tiledAtEntry).toBe(true);
  await setTrellisProps(page, {useTiledView: false});
  // Scenario 1 Step 3: DOM witness that the untiled render is alive before it is re-tiled
  // and has the SHAPE the step promises — one horizontal band of min(5, cats) cells.
  // rows = 1 is the settled untiled geometry [DOM 2026-08-11].
  const untiledEntry = await gridStructure(page);
  expect(untiledEntry.rows).toBe(1);
  expect(untiledEntry.count).toBe(untiledCells(splitCats));

  // #### Scenario 1 Step 4: Tiled View on rebuilds the band into a multi-row grid
  await softStep('Scenario 1 — Step 4', async () => {
    // The split column is wider than Tiles Per Row, so the two modes are GEOMETRICALLY
    // different; a narrower column would make this step pass on an unchanged grid.
    await setTrellisProps(page, {useTiledView: true});
    const s = await gridStructure(page);
    expect(s.rows).toBe(tiledRows(splitCats, DEFAULT_TILES_PER_ROW));
    expect(s.maxPerRow).toBe(tiledPerRow(splitCats, DEFAULT_TILES_PER_ROW));
    expect(s.count).toBe(tiledCells(splitCats, DEFAULT_TILES_PER_ROW));
    expect(s.rows).toBeGreaterThan(untiledEntry.rows);
    expect(s.count).not.toBe(untiledEntry.count);
    // Slots are not categories: the canvas-bearing cells are the category-count witness.
    expect(s.withCanvas).toBe(tiledCanvasCells(splitCats, DEFAULT_TILES_PER_ROW));
  });

  // #### Scenario 1 Step 5: Tiles Per Row = 1 re-flows the cells into a single column
  await softStep('Scenario 1 — Step 5', async () => {
    // Baseline at the default Tiles Per Row (4), where a row holds several cells, so the
    // drop to a single-column strip is a real delta and not the entry state.
    await setTrellisProps(page, {tilesPerRow: DEFAULT_TILES_PER_ROW});
    const atDefault = await gridStructure(page);
    expect(atDefault.maxPerRow).toBe(tiledPerRow(splitCats, DEFAULT_TILES_PER_ROW));
    expect(atDefault.count).toBe(tiledCells(splitCats, DEFAULT_TILES_PER_ROW));
    // Guard on the delta's precondition: a single-cell-wide baseline would make the
    // Tiles Per Row 1 read below indistinguishable from it.
    expect(atDefault.maxPerRow).toBeGreaterThan(1);

    // One cell per row, and the row count runs into the five-cell window before it runs out
    // of categories — so the slot count DROPS here (8 to 5 on DIS_POP) rather than holding.
    await setTrellisProps(page, {tilesPerRow: 1});
    const atOne = await gridStructure(page);
    expect(atOne.maxPerRow).toBe(tiledPerRow(splitCats, 1));
    expect(atOne.rows).toBe(tiledRows(splitCats, 1));
    expect(atOne.count).toBe(tiledCells(splitCats, 1));
    // The drop the step is about, stated against the baseline: 8 slots down to 5.
    expect(atOne.count).toBeLessThan(atDefault.count);
  });

  // #### Scenario 1 Step 6: Tiles Per Row = 3 reflows into rows exactly 3 cells wide
  await softStep('Scenario 1 — Step 6', async () => {
    // The row width is min(5, tilesPerRow, cats); with 6 categories the requested 3 is the
    // binding term, so the widest row is exactly 3 — the value the setting asked for, not
    // a clamp artefact. Coming from the width-1 strip of Step 5, this reads it both ways.
    await setTrellisProps(page, {tilesPerRow: 3});
    const s = await gridStructure(page);
    expect(s.maxPerRow).toBe(tiledPerRow(splitCats, 3));
    expect(s.rows).toBe(tiledRows(splitCats, 3));
    expect(s.count).toBe(tiledCells(splitCats, 3));
    expect(s.maxPerRow).toBeLessThan(tiledPerRow(splitCats, DEFAULT_TILES_PER_ROW));
  });

  // #### Scenario 1 Step 7: Tiled View off restores the single horizontal band (geometry round-trip)
  await softStep('Scenario 1 — Step 7', async () => {
    // Restore Tiles Per Row first, so turning Tiled View off is the only variable in the
    // round-trip. Untiled is a HORIZONTAL band of one row: yCategoriesViewportCount
    // degenerates to 1 with an empty Y list [DOM 2026-08-11].
    await setTrellisProps(page, {tilesPerRow: DEFAULT_TILES_PER_ROW});
    const tiled = await gridStructure(page);
    expect(tiled.rows).toBe(tiledRows(splitCats, DEFAULT_TILES_PER_ROW));
    await setTrellisProps(page, {useTiledView: false});
    const off = await gridStructure(page);
    expect(off.rows).toBe(1);
    expect(off.maxPerRow).toBe(untiledCells(splitCats));
    expect(off.count).toBe(untiledCells(splitCats));
    // A real geometry CHANGE, not shape stability: the multi-row tiled grid read a moment
    // ago collapses back into the single band.
    expect(off.rows).toBeLessThan(tiled.rows);
    expect(off.count).not.toBe(tiled.count);
    // Re-enable Tiled View for Scenario 2's entry state.
    await setTrellisProps(page, {useTiledView: true});
  });

  // #### Scenario 2 Step 3: control panel displayed with Auto Layout on and a large viewer
  await softStep('Scenario 2 — Step 3', async () => {
    // Auto-hide fires once the Trellis box falls under width 100 / height 300; at
    // 1920x1080 the box clears that threshold, so the panel must be showing.
    await setTrellisProps(page, {autoLayout: true, showControlPanel: true});
    await page.setViewportSize({width: 1920, height: 1080});
    await page.waitForTimeout(2500);
    const s = await controlPanelState(page);
    expect(s.exists).toBe(true);
    expect(s.display).not.toBe('none');
  });

  // #### Scenario 2 Step 4: shrinking below the threshold auto-hides the control panel
  await softStep('Scenario 2 — Step 4', async () => {
    // Auto-hide is a visibility override, not a property write: display goes none while
    // look.showControlPanel keeps its earlier value — hence the read BEFORE the resize, the
    // only way to witness that the auto-hide left the property alone. [DOM 2026-08-06]
    const before = await controlPanelState(page);
    expect(before.exists).toBe(true);
    await page.setViewportSize({width: 900, height: 280});
    await page.waitForTimeout(2800);
    const shrunk = await controlPanelState(page);
    // Hidden, not torn down — the two fields stay asserted apart, see controlPanelState.
    expect(shrunk.exists).toBe(true);
    expect(shrunk.display).toBe('none');
    expect(shrunk.showControlPanel).toBe(before.showControlPanel);
  });

  // #### Scenario 2 Step 5: restoring the size makes the control panel reappear
  await softStep('Scenario 2 — Step 5', async () => {
    // The panel has to come back on size alone, with no explicit toggle.
    const before = await controlPanelState(page);
    await page.setViewportSize({width: 1920, height: 1080});
    await page.waitForTimeout(2800);
    const restored = await controlPanelState(page);
    // Same node all along: it comes back out of display:none, it is not rebuilt.
    expect(restored.exists).toBe(true);
    expect(restored.display).not.toBe('none');
    expect(restored.showControlPanel).toBe(before.showControlPanel);
  });

  // #### Scenario 2 Step 6: Auto Layout off — the explicit showControlPanel value wins at any size
  await softStep('Scenario 2 — Step 6', async () => {
    // Auto Layout off disables the size gate: under the threshold it is the explicit
    // showControlPanel value that decides, in both directions.
    await setTrellisProps(page, {autoLayout: false, showControlPanel: true});
    await page.setViewportSize({width: 900, height: 280});
    await page.waitForTimeout(2800);
    const explicitOn = await controlPanelState(page);
    expect(explicitOn.exists).toBe(true);
    expect(explicitOn.display).not.toBe('none');

    await setTrellisProps(page, {showControlPanel: false});
    const explicitOff = await controlPanelState(page);
    // The explicit hide is the same display toggle as the auto-hide
    // (trellis_plot_core.dart:922) — the node survives it too.
    expect(explicitOff.exists).toBe(true);
    expect(explicitOff.display).toBe('none');

    // Restore defaults for cleanup independence.
    await setTrellisProps(page, {showControlPanel: true, autoLayout: true});
    await page.setViewportSize({width: 1920, height: 1080});
    await page.waitForTimeout(1500);
  });

  v.finishSpec();
});
