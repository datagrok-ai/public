/* ---
realizes: [pivottable.cp.configure-crosstab-values, pivottable.int.agg-types-track-agg-columns, pivottable.int.empty-aggregates-clear-pivot]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import * as proj from '../../helpers/projects';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const PIVOT = '[name="viewer-Pivot-table"]';

// Read a tag row's chip captions by row title (durable across rebuilds — the chip name
// attribute is dropped after any prop-driven rebuild, per pivot_table.md).
async function rowChips(page: Page, rowTitle: string): Promise<string[]> {
  return page.evaluate((title) => {
    const root = document.querySelector('[name="viewer-Pivot-table"]');
    const panel = Array.from(root!.querySelectorAll('.grok-pivot-column-panel'))
      .find((p) => p.querySelector('.grok-pivot-column-tags-title')?.getAttribute('d4-name') === title);
    if (!panel) return [];
    return Array.from(panel.querySelectorAll('.d4-tag')).map((t) => t.querySelector('span')?.textContent?.trim() ?? '');
  }, rowTitle);
}

async function pivotProps(page: Page) {
  return page.evaluate(() => {
    const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
    return {
      groupBy: pv.props.groupByColumnNames,
      pivot: pv.props.pivotColumnNames,
      agg: pv.props.aggregateColumnNames,
      aggTypes: pv.props.aggregateAggTypes,
    };
  });
}

// Make the view that hosts the Pivot Table viewer the current view. Clicking ADD opens a
// separate "Table aggregation" view that becomes current and backgrounds the pivot viewer,
// so any subsequent click on the pivot's + / chips times out (the element sits on a hidden
// view). Re-activate the pivot's own view before driving it again.
async function activatePivotView(page: Page) {
  await page.evaluate(() => {
    const host = Array.from(grok.shell.views).find((vw: any) =>
      Array.from(vw.viewers ?? []).some((x: any) => x.type === 'Pivot table'));
    if (host) grok.shell.v = host;
  });
  await page.waitForSelector('[name="viewer-Pivot-table"] [name="div-add-Group-by"]', {state: 'visible', timeout: 8000});
  await ensurePivotPlusClickable(page, 'div-add-Group-by');
}

// A tag-row "+" is `visible` once its view is current, but a transient overlay (a stray
// column-selector backdrop / context menu) can still intercept the pointer. Dismiss any stray
// popup, then poll until the named plus is the hit-test top element before the caller clicks it.
async function ensurePivotPlusClickable(page: Page, plusName: string) {
  const sel = `[name="viewer-Pivot-table"] [name="${plusName}"]`;
  for (let i = 0; i < 25; i++) {
    const clear = await page.evaluate((s) => {
      const plus = document.querySelector(s) as HTMLElement | null;
      if (!plus) return false;
      const r = plus.getBoundingClientRect();
      if (r.width === 0 || r.height === 0) return false;
      const top = document.elementFromPoint(r.x + r.width / 2, r.y + r.height / 2);
      return !!top && plus.contains(top);
    }, sel);
    if (clear) return;
    // Stray popup/backdrop still intercepting — Escape it and let it detach, then re-poll.
    await page.keyboard.press('Escape');
    await page.waitForTimeout(160);
  }
  throw new Error(`pivot ${plusName} + icon stayed obscured (overlay never cleared)`);
}

// Drive the tag-row "+" column picker: click the +, type the column name into the canvas
// column-grid, commit with Enter. The picker is a canvas-rendered column grid
// (.d4-column-selector-backdrop) driven by real keystrokes, not a DOM list.
async function addColumnViaPlus(page: Page, plusName: string, columnName: string) {
  await ensurePivotPlusClickable(page, plusName);
  await page.locator(`${PIVOT} [name="${plusName}"]`).click();
  await page.waitForSelector('.d4-column-selector-backdrop', {timeout: 6000});
  await page.keyboard.press(columnName[0]);
  await page.waitForTimeout(150);
  if (columnName.length > 1) await page.keyboard.type(columnName.slice(1));
  await page.waitForTimeout(150);
  await page.keyboard.press('Enter');
  await page.waitForTimeout(600);
}

// Pivot-tag menu mechanics: the chip context menu ([name="pivot-tag"]) is a d4 cascading
// vert-menu — child leaves keep a zero-size box until their parent is expanded by a sustained
// trusted hover (1px-jitter mouse moves), and the flyout collapses when the pointer leaves it.

// Right-click the first chip in the Aggregate row → opens its "pivot-tag" context menu.
async function openAggregateChipMenu(page: Page) {
  await page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Pivot-table"]');
    const panel = Array.from(root!.querySelectorAll('.grok-pivot-column-panel'))
      .find((p) => p.querySelector('.grok-pivot-column-tags-title')?.getAttribute('d4-name') === 'Aggregate');
    const chip = Array.from(panel!.querySelectorAll('.d4-tag'))[0];
    chip.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2}));
  });
  await page.locator('.d4-menu-popup[name="pivot-tag"] [name="div-Aggregation"]').waitFor({timeout: 5000});
}

// Hover a submenu parent (Aggregation / Column) until the named leaf lays out with a non-zero box.
async function expandSubmenu(page: Page, parent: 'Aggregation' | 'Column', leaf: string) {
  const box = await page.locator(`.d4-menu-popup[name="pivot-tag"] [name="div-${parent}"]`).boundingBox();
  if (!box) throw new Error(`pivot-tag ${parent} parent not found`);
  const px = box.x + box.width / 2, py = box.y + box.height / 2;
  const leafReady = () => page.evaluate((name) => {
    const el = document.querySelector(`.d4-menu-popup[name="pivot-tag"] [name="${name}"]`) as HTMLElement | null;
    if (!el) return false;
    const r = el.getBoundingClientRect();
    return r.width > 0 && r.height > 0;
  }, leaf);
  for (let i = 0; i < 25; i++) {
    await page.mouse.move(px + (i % 2), py);   // 1px jitter so each move is a distinct event
    await page.waitForTimeout(150);
    if (await leafReady()) return;
  }
  throw new Error(`pivot-tag ${parent} flyout did not reveal ${leaf}`);
}

// Click a flyout leaf while it is held open; a mid-transit collapse re-expands and retries.
async function clickFlyoutLeaf(page: Page, parent: 'Aggregation' | 'Column', leafName: string) {
  for (let attempt = 0; attempt < 3; attempt++) {
    await expandSubmenu(page, parent, leafName);
    const box = await page.locator(`.d4-menu-popup[name="pivot-tag"] [name="${leafName}"]`).boundingBox();
    if (!box) continue;
    const lx = box.x + box.width / 2, ly = box.y + box.height / 2;
    await page.mouse.move(lx, ly);
    await page.waitForTimeout(120);
    const onLeaf = await page.evaluate(({x, y, name}) => {
      const top = document.elementFromPoint(x, y);
      const el = document.querySelector(`.d4-menu-popup[name="pivot-tag"] [name="${name}"]`);
      return !!(top && el && (el.contains(top) || el === top));
    }, {x: lx, y: ly, name: leafName});
    if (!onLeaf) continue;                       // flyout collapsed mid-transit — re-expand and retry
    await page.mouse.down();
    await page.mouse.up();
    await page.waitForTimeout(400);
    return;
  }
  throw new Error(`could not click pivot-tag ${parent} leaf ${leafName}`);
}

// Click an Aggregation-type child (e.g. 'sum').
async function pickAggregation(page: Page, aggType: string) {
  await clickFlyoutLeaf(page, 'Aggregation', `div-Aggregation---${aggType}`);
}

// Click a Column child (e.g. 'HEIGHT'). The Column group is a radio group — picking a column
// REPLACES the measure's column, it does not add a second aggregate.
async function pickColumn(page: Page, colLeafName: string) {
  await clickFlyoutLeaf(page, 'Column', `div-Column---${colLeafName}`);
}

// The label of the currently-checked radio item inside the Aggregation group of the open menu.
// Requires the Aggregation submenu to be expanded (its children are zero-size otherwise); the
// caller expands it first via expandSubmenu / pickAggregation.
async function checkedAggregation(page: Page): Promise<string[]> {
  return page.evaluate(() => {
    const menu = document.querySelector('.d4-menu-popup[name="pivot-tag"]');
    if (!menu) return [];
    return Array.from(menu.querySelectorAll('[name^="div-Aggregation---"]'))
      .filter((mi) => mi.querySelector('.d4-menu-item-check i')?.className.includes('fa-dot-circle'))
      .map((mi) => mi.getAttribute('d4-name') ?? '');
  });
}

// Close any open pivot-tag context menu by pressing Escape and waiting for it to detach.
async function closePivotTagMenu(page: Page) {
  await page.keyboard.press('Escape');
  await page.locator('.d4-menu-popup[name="pivot-tag"]').waitFor({state: 'detached', timeout: 5000}).catch(() => {});
}

// Reveal a nested top-menu leaf: click the header, then keep the pointer resting (with 1px
// jitter) on the subgroup row so d4's open-delay timer flips the leaf flyout from display:none
// to laid-out — synthetic mouseover/mouseenter events do NOT drive this timer, so the
// header-click → group-hover → leaf-click sequence must use trusted CDP input throughout.
// Returns once the leaf reports a non-zero box.
async function openTopMenuLeaf(page: Page, header: string, group: string, leaf: string) {
  await page.locator(`[name="div-${header}"]`).click();
  await page.locator(`[name="div-${header}---${group}"]`).waitFor({timeout: 5000});
  const box = await page.locator(`[name="div-${header}---${group}"]`).boundingBox();
  if (!box) throw new Error(`top-menu group ${group} not found`);
  const gx = box.x + box.width / 2, gy = box.y + box.height / 2;
  const leafName = `div-${header}---${group}---${leaf}`;
  for (let i = 0; i < 25; i++) {
    await page.mouse.move(gx + (i % 2), gy);
    await page.waitForTimeout(150);
    const ready = await page.evaluate((name) => {
      const el = document.querySelector(`[name="${name}"]`) as HTMLElement | null;
      if (!el) return false;
      const r = el.getBoundingClientRect();
      return r.width > 0 && r.height > 0;
    }, leafName);
    if (ready) return leafName;
  }
  throw new Error(`top-menu leaf ${leafName} never became visible`);
}

test('Pivot Table — Configure cross-tab values', async ({page}) => {
  test.setTimeout(300_000);

  // Page-level console-error / pageerror collector — the build exposes neither
  // window.__consoleErrors nor an array grok.shell.warnings, so error guards read this
  // collector. Ignore the known harmless cloned-iframe noise.
  const consoleErrors: string[] = [];
  const IGNORE = /Unable to find element in cloned iframe/i;
  page.on('console', (m) => { if (m.type() === 'error' && !IGNORE.test(m.text())) consoleErrors.push(m.text()); });
  page.on('pageerror', (e) => { if (!IGNORE.test(e.message)) consoleErrors.push(e.message); });

  await loginToDatagrok(page);

  // --- Setup: open demog, add the Pivot Table viewer ---
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    grok.shell.addTableView(df);
    await new Promise((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
  await v.addViewerByIcon(page, 'pivot-table', 'Pivot-table', 15000);
  await page.waitForTimeout(800);

  await softStep('Setup: tag-editor header shows Group by, Aggregate and Pivot rows', async () => {
    const titles = await page.evaluate(() => Array.from(
      document.querySelectorAll('[name="viewer-Pivot-table"] .grok-pivot-column-tags-title'))
      .map((t) => t.getAttribute('d4-name')));
    expect(titles).toEqual(expect.arrayContaining(['Group by', 'Aggregate', 'Pivot']));
  });

  // === Scenario 1: Default configuration read-back and column picker guard ===

  await softStep('Scenario 1 Step 4: ADD publishes an aggregated table whose values match an independent groupBy', async () => {
    // The default auto-config reads DIS_POP / avg(AGE) / SEVERITY off the chips.
    expect(await rowChips(page, 'Group by')).toContain('DIS_POP');
    expect(await rowChips(page, 'Aggregate')).toContain('avg(AGE)');
    expect(await rowChips(page, 'Pivot')).toContain('SEVERITY');

    // Click ADD to publish the aggregated table to the workspace.
    await page.locator(`${PIVOT} .grok-pivot-counts [name="button-ADD"]`).click();
    await page.waitForTimeout(1000);

    const cmp = await page.evaluate(() => {
      const demog = grok.shell.tables.find((t: any) => t.rowCount === 5850);
      const indep = demog.groupBy(['DIS_POP']).pivot('SEVERITY').avg('AGE').aggregate();
      const indepCols = Array.from({length: indep.columns.length}, (_: any, i: number) => indep.columns.byIndex(i).name);
      const aggView = Array.from(grok.shell.views).find((vw: any) => vw.name === 'Table aggregation') as any;
      const pub = aggView?.table;
      const pubCols = pub ? Array.from({length: pub.columns.length}, (_: any, i: number) => pub.columns.byIndex(i).name) : [];
      const critName = indepCols.find((c: string) => c.startsWith('Critical'));
      return {
        published: !!pub,
        indepCols, pubCols,
        indepDisPop0: indep.col('DIS_POP').get(0),
        pubDisPop0: pub?.col('DIS_POP').get(0),
        indepCrit0: critName ? indep.col(critName).get(0) : null,
        pubCrit0: (critName && pub) ? pub.col(critName).get(0) : null,
      };
    });
    expect(cmp.published).toBe(true);
    expect(cmp.pubCols).toContain('DIS_POP');
    expect(cmp.pubCols).toEqual(cmp.indepCols);            // SEVERITY-derived headers match
    expect(cmp.pubDisPop0).toBe(cmp.indepDisPop0);          // key column matches
    expect(cmp.pubCrit0).toBe(cmp.indepCrit0);              // avg(AGE) value matches independent computation

    // ADD opened the "Table aggregation" view, which is now current and backgrounds the pivot
    // viewer. Re-activate the pivot's own view so the following steps can drive its + icons.
    await activatePivotView(page);
  });

  await softStep('Scenario 1 Step 7: opening the Group by column picker writes no console error (GROK-19114)', async () => {
    const before = consoleErrors.length;
    await ensurePivotPlusClickable(page, 'div-add-Group-by');
    await page.locator(`${PIVOT} [name="div-add-Group-by"]`).click();
    await page.waitForSelector('.d4-column-selector-backdrop', {timeout: 6000});
    await page.keyboard.press('Escape');
    await page.waitForTimeout(400);
    expect(consoleErrors.length).toBe(before);
  });

  await softStep('Scenario 1 Step 10: adding SEX in Group by yields a SEX chip and a DIS_POP+SEX grouping', async () => {
    await addColumnViaPlus(page, 'div-add-Group-by', 'SEX');
    expect(await rowChips(page, 'Group by')).toEqual(expect.arrayContaining(['DIS_POP', 'SEX']));
    const props = await pivotProps(page);
    expect(props.groupBy).toEqual(expect.arrayContaining(['DIS_POP', 'SEX']));
    // Row count of the DIS_POP×SEX grouping equals the number of distinct pairs in demog.
    const pairs = await page.evaluate(() => {
      const demog = grok.shell.tables.find((t: any) => t.rowCount === 5850);
      return demog.groupBy(['DIS_POP', 'SEX']).aggregate().rowCount;
    });
    expect(pairs).toBe(12);
  });

  await softStep('Scenario 1 Step 11: close the published aggregated table; restore default group-by', async () => {
    await page.evaluate(async () => {
      const aggView = Array.from(grok.shell.views).find((vw: any) => vw.name === 'Table aggregation') as any;
      if (aggView) aggView.close();
      await new Promise((r) => setTimeout(r, 300));
    });
    // Restore the single-key default so later scenarios start clean.
    await addColumnViaPlus(page, 'div-add-Group-by', 'DIS_POP');
  });

  // === Scenario 2: Aggregate tag context menu — multi-pick in one open menu ===

  await softStep('Scenario 2 Step 3: right-click avg(AGE) opens the menu with avg checked, no console error (GROK-17841)', async () => {
    // Reset to the default single aggregate avg(AGE) before the menu is opened.
    await page.evaluate(() => {
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      pv.props.groupByColumnNames = ['DIS_POP'];
      pv.props.aggregateColumnNames = ['AGE'];
      pv.props.aggregateAggTypes = ['avg'];
      pv.props.pivotColumnNames = ['SEVERITY'];
    });
    await page.waitForTimeout(500);
    const before = consoleErrors.length;
    await openAggregateChipMenu(page);
    // The Aggregation-type children are zero-size until the Aggregation parent is hovered.
    await expandSubmenu(page, 'Aggregation', 'div-Aggregation---avg');
    expect(await checkedAggregation(page)).toEqual(['avg']);
    expect(consoleErrors.length).toBe(before);
  });

  await softStep('Scenario 2 Step 5: picking Sum moves the mark, updates the chip, keeps the menu open (GROK-16899)', async () => {
    await pickAggregation(page, 'sum');
    expect(await page.locator('.d4-menu-popup[name="pivot-tag"]').count()).toBeGreaterThan(0);
    await expandSubmenu(page, 'Aggregation', 'div-Aggregation---sum');
    expect(await checkedAggregation(page)).toEqual(['sum']);
    expect(await rowChips(page, 'Aggregate')).toContain('sum(AGE)');
    expect((await pivotProps(page)).aggTypes).toEqual(['sum']);
  });

  await softStep('Scenario 2 Step 6: a second pick (Median) in the same open menu also takes effect (GROK-16899)', async () => {
    await pickAggregation(page, 'med');
    expect(await page.locator('.d4-menu-popup[name="pivot-tag"]').count()).toBeGreaterThan(0);
    await expandSubmenu(page, 'Aggregation', 'div-Aggregation---med');
    expect(await checkedAggregation(page)).toEqual(['med']);
    expect(await rowChips(page, 'Aggregate')).toContain('med(AGE)');
  });

  await softStep('Scenario 2 Step 8: switching column to HEIGHT rebuilds the Aggregation group to its supported types (I2)', async () => {
    await pickColumn(page, 'HEIGHT');
    expect(await page.locator('.d4-menu-popup[name="pivot-tag"]').count()).toBeGreaterThan(0);
    // Chip reflects the newly picked column; med (unsupported default carry) replaced by a type default.
    expect(await rowChips(page, 'Aggregate')).toContain('avg(HEIGHT)');
    // Expand the (rebuilt) Aggregation submenu to read exactly which types it now offers.
    await expandSubmenu(page, 'Aggregation', 'div-Aggregation---avg');
    const offered = await page.evaluate(() => Array.from(
      document.querySelectorAll('.d4-menu-popup[name="pivot-tag"] [name^="div-Aggregation---"]'))
      .map((mi) => mi.getAttribute('d4-name')));
    // I2: the Aggregation group is rebuilt to EXACTLY the numeric-type aggregation set —
    // no unsupported (string-only) item leaks in. Set-equality (sorted) fails on an extra
    // unsupported item AND on a missing supported one.
    const NUMERIC_AGG = ['first', 'count', 'values', 'unique', 'nulls', 'min', 'max', 'sum',
      'med', 'avg', 'geomean', 'stdev', 'variance', 'skew', 'kurt', 'q1', 'q2', 'q3'];
    expect([...offered].sort()).toEqual([...NUMERIC_AGG].sort());
    expect((await pivotProps(page)).agg).toEqual(['HEIGHT']);
  });

  await softStep('Scenario 2 Step 10: Remove others closes the menu, one chip remains, no console error', async () => {
    // The Column group is a radio group: picking WEIGHT REPLACES the measure's column (it does
    // not add a second aggregate), so one chip remains — now avg(WEIGHT).
    await pickColumn(page, 'WEIGHT');
    expect(await rowChips(page, 'Aggregate')).toContain('avg(WEIGHT)');
    // Changing the column via this menu rebuilds the tag row, invalidating the TagElement the
    // still-open menu captured; "Remove others" removes every chip whose TagElement != that stale
    // capture, so from the same open menu it would clear the row. Reopen the menu fresh so
    // "Remove others" binds to the live TagElement.
    await closePivotTagMenu(page);
    await openAggregateChipMenu(page);
    const before = consoleErrors.length;
    // "Remove others" is a top-level pivot-tag item (no submenu expansion needed). It removes
    // every OTHER chip; with a single chip it leaves that one and closes the menu.
    await page.locator('.d4-menu-popup[name="pivot-tag"] [name="div-Remove-others"]').click();
    await page.waitForTimeout(500);
    expect(await page.locator('.d4-menu-popup[name="pivot-tag"]').count()).toBe(0);
    expect((await rowChips(page, 'Aggregate')).length).toBe(1);
    expect(consoleErrors.length).toBe(before);
  });

  // === Scenario 3: Pivot row hidden when last aggregate is removed (I4) ===

  await softStep('Scenario 3 Step 3: removing the last aggregate hides the Pivot row and clears pivot columns (I4)', async () => {
    await page.evaluate(() => {
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      pv.props.groupByColumnNames = ['DIS_POP'];
      pv.props.aggregateColumnNames = ['AGE'];
      pv.props.aggregateAggTypes = ['avg'];
      pv.props.pivotColumnNames = ['SEVERITY'];
    });
    await page.waitForTimeout(500);
    // Remove the (only) aggregate chip via its × icon.
    await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Pivot-table"]');
      const panel = Array.from(root!.querySelectorAll('.grok-pivot-column-panel'))
        .find((p) => p.querySelector('.grok-pivot-column-tags-title')?.getAttribute('d4-name') === 'Aggregate');
      const removeIcon = panel!.querySelector('.d4-tag i, .d4-tag .grok-icon') as HTMLElement | null;
      removeIcon?.click();
    });
    await page.waitForTimeout(600);
    const state = await page.evaluate(() => {
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      const root = document.querySelector('[name="viewer-Pivot-table"]');
      const pivotPanel = Array.from(root!.querySelectorAll('.grok-pivot-column-panel'))
        .find((p) => p.querySelector('.grok-pivot-column-tags-title')?.getAttribute('d4-name') === 'Pivot') as HTMLElement | undefined;
      const pivotVisible = !!pivotPanel && pivotPanel.offsetParent !== null;
      return {aggEmpty: pv.props.aggregateColumnNames.length === 0, pivotCleared: pv.props.pivotColumnNames.length === 0, pivotVisible};
    });
    expect(state.aggEmpty).toBe(true);
    expect(state.pivotCleared).toBe(true);   // configured pivot columns silently cleared (config lost, not hidden)
    expect(state.pivotVisible).toBe(false);  // Pivot row hidden from the viewer
  });

  await softStep('Scenario 3 Step 5: re-adding an aggregate brings the Pivot row back', async () => {
    await addColumnViaPlus(page, 'div-add-Aggregate', 'AGE');
    const aggChips = await rowChips(page, 'Aggregate');
    expect(aggChips.length).toBe(1);
    expect(aggChips[0]).toContain('AGE');
    const pivotVisible = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Pivot-table"]');
      const pivotPanel = Array.from(root!.querySelectorAll('.grok-pivot-column-panel'))
        .find((p) => p.querySelector('.grok-pivot-column-tags-title')?.getAttribute('d4-name') === 'Pivot') as HTMLElement | undefined;
      return !!pivotPanel && pivotPanel.offsetParent !== null;
    });
    expect(pivotVisible).toBe(true);
  });

  // === Scenario 4: Property panel configuration path (GROK-16305 guard) ===
  //
  // The pivot column-list editors live in the VIEWER's property panel, which surfaces only when the
  // pivot viewer is the current object (the title-bar gear opens the inner GRID's props instead).
  // Each editor [name="prop-view-<row>"] is a "<n> / <total>" count label + a "..." button that
  // opens the DOM dialog [name="dialog-Select-columns..."].

  // Make the pivot viewer the current object so its property panel shows the column-list editors.
  // A prop-driven refresh steals the current-object focus asynchronously, so this is done AFTER any
  // prop write settles, and the prop rows are polled until they appear.
  async function selectViewerObject() {
    await page.evaluate(() => {
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      grok.shell.o = pv;
    });
    await page.waitForSelector('[name="prop-aggregate"]', {timeout: 8000});
  }

  // Surface the panel, then open the "..." dialog for the given row. Leaves the dialog open.
  async function openPanelColumnDialog(row: string) {
    await selectViewerObject();
    await page.waitForSelector(`[name="prop-${row}"]`, {timeout: 6000});
    await page.locator(`[name="prop-view-${row}"] button`).click();      // the "..." button (DOM)
    await page.waitForSelector('[name="dialog-Select-columns..."]', {timeout: 6000});
  }

  const panelCount = (row: string) => page.locator(`[name="prop-view-${row}"] label`).innerText();

  // Toggle the FIRST data row's checkbox in the open Select-columns dialog after typing `col` into the
  // search field. Search RE-SORTS (not filters), so `col` becomes row 1; the checkbox is a canvas-grid
  // cell (no DOM handle) at a stable offset from the live-read grid rect — driven by a trusted
  // page.mouse.click. Replace recipe: the caller searches the drop column + toggles, then the add
  // column + toggles.
  async function toggleFirstRow(col: string) {
    const dlg = page.locator('[name="dialog-Select-columns..."]');
    const search = dlg.locator('input.d4-search-input');
    await search.click();
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Delete');
    await page.keyboard.type(col);                 // keyboard type fires the Dart change listener
    await page.waitForTimeout(400);                // let the list re-sort so `col` is row 1
    const rect = await dlg.locator('.d4-grid').boundingBox();
    if (!rect) throw new Error('Select-columns grid not visible');
    // Row-1 checkbox: right edge of the grid, header band + first row. Coordinate is computed from the
    // live rect each call, so it is stable regardless of where the dialog opens.
    await page.mouse.click(rect.x + rect.width - 39, rect.y + 34);
    await page.waitForTimeout(400);
  }

  await softStep('Scenario 4 Step 3: the viewer property panel exposes the group-by / aggregate / pivot column-list editors', async () => {
    await page.evaluate(() => {
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      pv.props.groupByColumnNames = ['DIS_POP'];
      pv.props.aggregateColumnNames = ['AGE'];
      pv.props.aggregateAggTypes = ['avg'];
      pv.props.pivotColumnNames = ['SEVERITY'];
    });
    await page.waitForTimeout(700);             // let the prop-driven refresh settle (it steals focus)
    await selectViewerObject();                 // then select the viewer → its property panel appears
    const propsPresent = await page.evaluate(() => ({
      groupBy: !!document.querySelector('[name="prop-group-by"]'),
      aggregate: !!document.querySelector('[name="prop-aggregate"]'),
      pivot: !!document.querySelector('[name="prop-pivot"]'),
    }));
    expect(propsPresent.groupBy).toBe(true);
    expect(propsPresent.aggregate).toBe(true);
    expect(propsPresent.pivot).toBe(true);
  });

  await softStep('Scenario 4 Step 5: replacing AGE with WEIGHT via the aggregate editor rewrites aggregateColumnNames with no console error (GROK-16305)', async () => {
    const before = consoleErrors.length;
    await openPanelColumnDialog('aggregate');
    const dlg = page.locator('[name="dialog-Select-columns..."]');
    // Search AGE → uncheck row 1; search WEIGHT → check row 1. Each edit fires the panel change handler.
    await toggleFirstRow('AGE');                                   // drop AGE
    expect(await panelCount('aggregate')).toBe('0 / 11');          // handler fired: AGE unchecked, count → 0
    await toggleFirstRow('WEIGHT');                                // add WEIGHT
    expect(await panelCount('aggregate')).toBe('1 / 11');          // handler fired: WEIGHT checked, count → 1
    await dlg.locator('[name="button-OK"]').click();
    await page.waitForTimeout(600);
    const agg = (await pivotProps(page)).agg;
    expect(agg).toContain('WEIGHT');
    expect(agg).not.toContain('AGE');
    expect(consoleErrors.length).toBe(before);                     // GROK-16305: no console error on the panel path
    // Restore the single-aggregate default (reset only, not the tested action).
    await page.evaluate(() => {
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      pv.props.aggregateColumnNames = ['AGE']; pv.props.aggregateAggTypes = ['avg'];
    });
    await page.waitForTimeout(400);
  });

  await softStep('Scenario 4 Step 7: replacing DIS_POP with SEX via the group-by editor rewrites groupByColumnNames with no console error (GROK-16305)', async () => {
    const before = consoleErrors.length;
    await openPanelColumnDialog('group-by');
    const dlg = page.locator('[name="dialog-Select-columns..."]');
    // Search DIS_POP → uncheck row 1; search SEX → check row 1.
    await toggleFirstRow('DIS_POP');                               // drop DIS_POP
    expect(await panelCount('group-by')).toBe('0 / 11');           // handler fired: DIS_POP unchecked, count → 0
    await toggleFirstRow('SEX');                                   // add SEX
    expect(await panelCount('group-by')).toBe('1 / 11');           // handler fired: SEX checked, count → 1
    await dlg.locator('[name="button-OK"]').click();
    await page.waitForTimeout(600);
    const gb = (await pivotProps(page)).groupBy;
    expect(gb).toContain('SEX');
    expect(gb).not.toContain('DIS_POP');
    expect(consoleErrors.length).toBe(before);                     // GROK-16305: no console error on the panel path
    // Restore the single group-by default (reset only).
    await page.evaluate(() => {
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      pv.props.groupByColumnNames = ['DIS_POP'];
    });
    await page.waitForTimeout(400);
  });

  // === Scenario 5: Refresh resets configuration; layout and project persistence (github-2535) ===

  await softStep('Scenario 5 Step 4: Refresh resets the visible configuration to type-driven defaults', async () => {
    // Establish a NON-default visible configuration through the tag rows so the chips are named
    // and the reset is observable as a caption change (not a silent no-op).
    await page.evaluate(() => {
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      pv.props.groupByColumnNames = ['DIS_POP'];
      pv.props.aggregateColumnNames = ['AGE'];
      pv.props.aggregateAggTypes = ['med'];
      pv.props.pivotColumnNames = ['SEVERITY'];
    });
    await page.waitForTimeout(600);
    // Pre-Refresh the visible aggregate is the non-default med(AGE); group-by/pivot carry a chip.
    expect(await rowChips(page, 'Aggregate')).toContain('med(AGE)');
    expect((await rowChips(page, 'Group by')).length).toBeGreaterThan(0);

    await page.locator(`${PIVOT} .d4-command-bar [name="icon-redo"]`).click();
    await page.waitForTimeout(900);

    // Refresh clears group-by and pivot outright and sets the aggregations to getDefaultAggregations()
    // = avg over the first two numerical columns. The look props are NOT written back (setTags(notify:
    // false) raises no TAGS_CHANGED), so the reset is read off the visible tag chips — the honest
    // observable — never the props, which stay at their pre-Refresh values regardless of the reset.
    const gbAfter = await rowChips(page, 'Group by');
    const pivotAfter = await rowChips(page, 'Pivot');
    const aggAfter = await rowChips(page, 'Aggregate');
    expect(gbAfter).toEqual([]);
    expect(pivotAfter).toEqual([]);
    expect(aggAfter).toEqual(['avg(AGE)', 'avg(HEIGHT)']);
  });

  await softStep('Scenario 5 Step 9: a non-default aggregation survives a layout saved through the View menu and re-applied via the layout API (github-2535)', async () => {
    // Set a NON-DEFAULT aggregation, then save the layout through the real UI:
    // View -> Layout -> Save to Gallery (a documented top-menu leaf). The github-2535
    // regression lived in the UI persistence path, so the save leg is driven, not API'd.
    await page.evaluate(async () => {
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      pv.props.groupByColumnNames = ['DIS_POP'];
      pv.props.aggregateColumnNames = ['AGE'];
      pv.props.aggregateAggTypes = ['med'];
      pv.props.pivotColumnNames = ['SEVERITY'];
      await new Promise((r) => setTimeout(r, 500));
    });
    // Record the layout set present before the save so we can find the newly-saved one.
    const beforeIds = await page.evaluate(async () =>
      (await grok.dapi.layouts.getApplicable(grok.shell.tv.dataFrame)).map((l: any) => l.id));

    // Drive the save via the View -> Layout -> Save to Gallery top-menu leaf (trusted flow —
    // see openTopMenuLeaf).
    const saveLeaf = await openTopMenuLeaf(page, 'View', 'Layout', 'Save-to-Gallery');
    const lb = await page.locator(`[name="${saveLeaf}"]`).boundingBox();
    await page.mouse.move(lb!.x + lb!.width / 2, lb!.y + lb!.height / 2);
    await page.waitForTimeout(120);
    await page.mouse.down();
    await page.mouse.up();
    await page.waitForTimeout(1500);

    // Identify the layout the UI save just created (the one not present before).
    const savedId = await page.evaluate(async (prev) => {
      const applicable = await grok.dapi.layouts.getApplicable(grok.shell.tv.dataFrame);
      const fresh = applicable.filter((l: any) => !prev.includes(l.id));
      return (fresh[0] ?? applicable[applicable.length - 1])?.id ?? null;
    }, beforeIds);
    expect(savedId).not.toBeNull();                      // View -> Layout -> Save to Gallery persisted a layout

    // Perturb the visible configuration so the re-apply is observable as a real change.
    await page.evaluate(async () => {
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      pv.props.aggregateAggTypes = ['avg'];
      pv.props.groupByColumnNames = ['SEX'];
      await new Promise((r) => setTimeout(r, 400));
    });

    // Re-apply the saved layout. The individual saved-layout entry in the Toolbox Layouts
    // gallery has no captured selector, so the re-apply leg uses the JS-API apply path;
    // the SAVE leg above is the UI persistence channel the github-2535 guard exercises.
    const result = await page.evaluate(async (id) => {
      const tv = grok.shell.tv;
      const saved = await grok.dapi.layouts.find(id);
      tv.loadLayout(saved);
      await new Promise((r) => setTimeout(r, 3000));
      const pv2 = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      const restored = {
        member: !!pv2,
        groupBy: pv2?.props.groupByColumnNames,
        pivot: pv2?.props.pivotColumnNames,
        agg: pv2?.props.aggregateColumnNames,
        aggTypes: pv2?.props.aggregateAggTypes,
      };
      await grok.dapi.layouts.delete(saved);
      return restored;
    }, savedId);
    expect(result.member).toBe(true);
    expect(result.groupBy).toContain('DIS_POP');
    expect(result.pivot).toContain('SEVERITY');
    expect(result.agg).toContain('AGE');
    expect(result.aggTypes).toContain('med');            // NON-DEFAULT aggregation restored (github-2535)
  });

  await softStep('Scenario 5 Step 13: the configuration survives a project save / close / reopen round-trip (github-2535)', async () => {
    const projectName = `PivotCrosstabProj${Date.now()}`;

    // Set the NON-DEFAULT configuration to persist.
    await page.evaluate(async () => {
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      pv.props.groupByColumnNames = ['DIS_POP'];
      pv.props.aggregateColumnNames = ['AGE'];
      pv.props.aggregateAggTypes = ['med'];
      pv.props.pivotColumnNames = ['SEVERITY'];
      await new Promise((r) => setTimeout(r, 500));
    });

    // Save the workspace as a project through the real ribbon SAVE button (saveProjectViaUI clicks
    // [name="button-Save"], fills the name, dismisses the follow-up Share dialog, then polls until
    // the project is queryable). This is the ONE save path that carries the VIEW LAYOUT — a JS-API
    // project build attaches a ViewLayout that is never re-applied on reopen, so the pivot viewer +
    // its props would not come back. It is also the UI persistence channel of github-2535.
    // (Ctrl+S does not open the save dialog; the ribbon [name="button-Save"] does.)
    const saved = await proj.saveProjectViaUI(page, projectName);
    expect(saved.projectId).toBeTruthy();                 // project persisted server-side

    // Close every view, then reopen the saved project. p.open() is the UI-equivalent of
    // double-clicking the Dashboards tile (both restore tables, viewers, and layouts).
    await page.evaluate(() => grok.shell.closeAll());
    await page.waitForTimeout(800);
    const opened = await page.evaluate(async (id) => {
      const dash = await grok.dapi.projects.find(id);
      if (dash) await dash.open();
      await new Promise((r) => setTimeout(r, 4000));
      return !!grok.shell.tv;
    }, saved.projectId);
    expect(opened).toBe(true);

    const result = await page.evaluate(() => {
      const pv2 = Array.from(grok.shell.tv?.viewers ?? []).find((x: any) => x.type === 'Pivot table') as any;
      return {
        present: !!pv2,
        groupBy: pv2?.props.groupByColumnNames,
        pivot: pv2?.props.pivotColumnNames,
        agg: pv2?.props.aggregateColumnNames,
        aggTypes: pv2?.props.aggregateAggTypes,
      };
    });
    expect(result.present).toBe(true);
    expect(result.groupBy).toContain('DIS_POP');
    expect(result.pivot).toContain('SEVERITY');
    expect(result.aggTypes).toContain('med');             // NON-DEFAULT aggregation survives the project round-trip

    // Cleanup: delete the saved project.
    await proj.deleteProjectWithCleanup(page, {projectId: saved.projectId});
  });

  v.finishSpec();
});
