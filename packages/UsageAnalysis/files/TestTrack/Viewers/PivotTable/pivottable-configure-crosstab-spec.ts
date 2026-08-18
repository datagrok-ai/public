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

async function activatePivotView(page: Page) {
  await page.evaluate(() => {
    const host = Array.from(grok.shell.views).find((vw: any) =>
      Array.from(vw.viewers ?? []).some((x: any) => x.type === 'Pivot table'));
    if (host) grok.shell.v = host;
  });
  await page.waitForSelector('[name="viewer-Pivot-table"] [name="div-add-Group-by"]', {state: 'visible', timeout: 8000});
  await ensurePivotPlusClickable(page, 'div-add-Group-by');
}

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

    await page.keyboard.press('Escape');
    await page.waitForTimeout(160);   
  }
  throw new Error(`pivot ${plusName} + icon stayed obscured (overlay never cleared)`);
}

async function addColumnViaPlus(page: Page, plusName: string, columnName: string) {
  await ensurePivotPlusClickable(page, plusName);
  await page.locator(`${PIVOT} [name="${plusName}"]`).click();
  await page.waitForSelector('.d4-column-selector-backdrop', {timeout: 6000});
  await page.keyboard.press(columnName[0]);
  await page.waitForTimeout(150);   
  if (columnName.length > 1) await page.keyboard.type(columnName.slice(1));
  await page.waitForTimeout(150);   
  await page.keyboard.press('Enter');
  await v.waitForViewerRendered(page, 'Pivot table', 600);
}

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
    await page.mouse.move(px + (i % 2), py);   
    await page.waitForTimeout(150);            
    if (await leafReady()) return;
  }
  throw new Error(`pivot-tag ${parent} flyout did not reveal ${leaf}`);
}

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
    if (!onLeaf) continue;                       
    await page.mouse.down();
    await page.mouse.up();
    await v.waitForViewerRendered(page, 'Pivot table', 400);
    return;
  }
  throw new Error(`could not click pivot-tag ${parent} leaf ${leafName}`);
}

async function pickAggregation(page: Page, aggType: string) {
  await clickFlyoutLeaf(page, 'Aggregation', `div-Aggregation---${aggType}`);
}

async function pickColumn(page: Page, colLeafName: string) {
  await clickFlyoutLeaf(page, 'Column', `div-Column---${colLeafName}`);
}

async function checkedAggregation(page: Page): Promise<string[]> {
  return page.evaluate(() => {
    const menu = document.querySelector('.d4-menu-popup[name="pivot-tag"]');
    if (!menu) return [];
    return Array.from(menu.querySelectorAll('[name^="div-Aggregation---"]'))
      .filter((mi) => mi.querySelector('.d4-menu-item-check i')?.className.includes('fa-dot-circle'))
      .map((mi) => mi.getAttribute('d4-name') ?? '');
  });
}

async function closePivotTagMenu(page: Page) {
  await page.keyboard.press('Escape');
  await page.locator('.d4-menu-popup[name="pivot-tag"]').waitFor({state: 'detached', timeout: 5000}).catch(() => {});
}

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

  const consoleErrors: string[] = [];
  const IGNORE = /Unable to find element in cloned iframe/i;
  page.on('console', (m) => { if (m.type() === 'error' && !IGNORE.test(m.text())) consoleErrors.push(m.text()); });
  page.on('pageerror', (e) => { if (!IGNORE.test(e.message)) consoleErrors.push(e.message); });

  await loginToDatagrok(page);

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
  await v.waitForViewerRendered(page, 'Pivot table', 800);

  await softStep('Setup: tag-editor header shows Group by, Aggregate and Pivot rows', async () => {
    const titles = await page.evaluate(() => Array.from(
      document.querySelectorAll('[name="viewer-Pivot-table"] .grok-pivot-column-tags-title'))
      .map((t) => t.getAttribute('d4-name')));
    expect(titles).toEqual(expect.arrayContaining(['Group by', 'Aggregate', 'Pivot']));
  });

  await softStep('Scenario 1 Step 4: ADD publishes an aggregated table whose values match an independent groupBy', async () => {

    expect(await rowChips(page, 'Group by')).toContain('DIS_POP');
    expect(await rowChips(page, 'Aggregate')).toContain('avg(AGE)');
    expect(await rowChips(page, 'Pivot')).toContain('SEVERITY');

    await page.locator(`${PIVOT} .grok-pivot-counts [name="button-ADD"]`).click();
    await v.pollValue(() => page.evaluate(() =>
      Array.from(grok.shell.views).some((vw: any) => vw.name === 'Table aggregation')),
    (present) => present, 1000, 100);

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
    expect(cmp.pubCols).toEqual(cmp.indepCols);            
    expect(cmp.pubDisPop0).toBe(cmp.indepDisPop0);          
    expect(cmp.pubCrit0).toBe(cmp.indepCrit0);              

    await activatePivotView(page);
  });

  await softStep('Scenario 1 Step 7: opening the Group by column picker writes no console error (GROK-19114)', async () => {
    const before = consoleErrors.length;
    await ensurePivotPlusClickable(page, 'div-add-Group-by');
    await page.locator(`${PIVOT} [name="div-add-Group-by"]`).click();
    await page.waitForSelector('.d4-column-selector-backdrop', {timeout: 6000});
    await page.keyboard.press('Escape');
    await v.pollValue(() => page.locator('.d4-column-selector-backdrop').count(), (n) => n === 0, 400, 100);
    expect(consoleErrors.length).toBe(before);
  });

  await softStep('Scenario 1 Step 10: adding SEX in Group by yields a SEX chip and a DIS_POP+SEX grouping', async () => {
    await addColumnViaPlus(page, 'div-add-Group-by', 'SEX');
    expect(await rowChips(page, 'Group by')).toEqual(expect.arrayContaining(['DIS_POP', 'SEX']));
    const props = await pivotProps(page);
    expect(props.groupBy).toEqual(expect.arrayContaining(['DIS_POP', 'SEX']));

    const pairs = await page.evaluate(() => {
      const demog = grok.shell.tables.find((t: any) => t.rowCount === 5850);
      return demog.groupBy(['DIS_POP', 'SEX']).aggregate().rowCount;
    });
    expect(pairs).toBe(12);
  });

  await softStep('Scenario 1 Step 11: close the published aggregated table; restore default group-by', async () => {
    await page.evaluate(() => {
      const aggView = Array.from(grok.shell.views).find((vw: any) => vw.name === 'Table aggregation') as any;
      if (aggView) aggView.close();
    });
    await v.pollValue(() => page.evaluate(() =>
      Array.from(grok.shell.views).some((vw: any) => vw.name === 'Table aggregation')),
    (present) => !present, 300, 100);

    await addColumnViaPlus(page, 'div-add-Group-by', 'DIS_POP');
  });

  await softStep('Scenario 2 Step 3: right-click avg(AGE) opens the menu with avg checked, no console error (GROK-17841)', async () => {

    await v.setViewerProps(page, 'Pivot table', [{set: {
      groupByColumnNames: ['DIS_POP'],
      aggregateColumnNames: ['AGE'],
      aggregateAggTypes: ['avg'],
      pivotColumnNames: ['SEVERITY'],
    }}], 500);
    const before = consoleErrors.length;
    await openAggregateChipMenu(page);

    await expandSubmenu(page, 'Aggregation', 'div-Aggregation---avg');
    expect(await checkedAggregation(page)).toEqual(['avg']);
    expect(consoleErrors.length).toBe(before);
  });

  await softStep('Scenario 2 Step 5: picking Sum moves the mark, updates the chip, keeps the menu open (GROK-16899)', async () => {
    await pickAggregation(page, 'sum');
    expect(await page.locator('.d4-menu-popup[name="pivot-tag"]').count()).toBeGreaterThan(0);
    await expandSubmenu(page, 'Aggregation', 'div-Aggregation---sum');
    expect(await checkedAggregation(page)).toEqual(['sum']);

    expect(await v.pollValue(() => rowChips(page, 'Aggregate'),
      (c) => c.includes('sum(AGE)'), 3000, 150)).toContain('sum(AGE)');
    expect((await pivotProps(page)).aggTypes).toEqual(['sum']);
  });

  await softStep('Scenario 2 Step 6: a second pick (Median) in the same open menu also takes effect (GROK-16899)', async () => {
    await pickAggregation(page, 'med');
    expect(await page.locator('.d4-menu-popup[name="pivot-tag"]').count()).toBeGreaterThan(0);
    await expandSubmenu(page, 'Aggregation', 'div-Aggregation---med');
    expect(await checkedAggregation(page)).toEqual(['med']);
    expect(await v.pollValue(() => rowChips(page, 'Aggregate'),
      (c) => c.includes('med(AGE)'), 3000, 150)).toContain('med(AGE)');
  });

  await softStep('Scenario 2 Step 8: switching column to HEIGHT rebuilds the Aggregation group to its supported types (I2)', async () => {
    await pickColumn(page, 'HEIGHT');
    expect(await page.locator('.d4-menu-popup[name="pivot-tag"]').count()).toBeGreaterThan(0);

    expect(await v.pollValue(() => rowChips(page, 'Aggregate'),
      (c) => c.includes('avg(HEIGHT)'), 3000, 150)).toContain('avg(HEIGHT)');

    await expandSubmenu(page, 'Aggregation', 'div-Aggregation---avg');
    const offered = await page.evaluate(() => Array.from(
      document.querySelectorAll('.d4-menu-popup[name="pivot-tag"] [name^="div-Aggregation---"]'))
      .map((mi) => mi.getAttribute('d4-name')));

    const NUMERIC_AGG = ['first', 'count', 'values', 'unique', 'nulls', 'min', 'max', 'sum',
      'med', 'avg', 'geomean', 'stdev', 'variance', 'skew', 'kurt', 'q1', 'q2', 'q3'];
    expect([...offered].sort()).toEqual([...NUMERIC_AGG].sort());
    expect((await pivotProps(page)).agg).toEqual(['HEIGHT']);
  });

  await softStep('Scenario 2 Step 10: Remove others closes the menu, one chip remains, no console error', async () => {

    await pickColumn(page, 'WEIGHT');
    expect(await v.pollValue(() => rowChips(page, 'Aggregate'),
      (c) => c.includes('avg(WEIGHT)'), 3000, 150)).toContain('avg(WEIGHT)');

    await closePivotTagMenu(page);
    await openAggregateChipMenu(page);
    const before = consoleErrors.length;

    await page.locator('.d4-menu-popup[name="pivot-tag"] [name="div-Remove-others"]').click();
    await v.pollValue(() => page.locator('.d4-menu-popup[name="pivot-tag"]').count(), (n) => n === 0, 500, 100);
    expect(await page.locator('.d4-menu-popup[name="pivot-tag"]').count()).toBe(0);
    expect((await rowChips(page, 'Aggregate')).length).toBe(1);
    expect(consoleErrors.length).toBe(before);
  });

  await softStep('Scenario 3 Step 3: removing the last aggregate hides the Pivot row and clears pivot columns (I4)', async () => {
    await v.setViewerProps(page, 'Pivot table', [{set: {
      groupByColumnNames: ['DIS_POP'],
      aggregateColumnNames: ['AGE'],
      aggregateAggTypes: ['avg'],
      pivotColumnNames: ['SEVERITY'],
    }}], 500);

    await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Pivot-table"]');
      const panel = Array.from(root!.querySelectorAll('.grok-pivot-column-panel'))
        .find((p) => p.querySelector('.grok-pivot-column-tags-title')?.getAttribute('d4-name') === 'Aggregate');
      const removeIcon = panel!.querySelector('.d4-tag i, .d4-tag .grok-icon') as HTMLElement | null;
      removeIcon?.click();
    });
    await v.waitForViewerRendered(page, 'Pivot table', 600);

    const state = await v.pollValue(() => page.evaluate(() => {
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      const root = document.querySelector('[name="viewer-Pivot-table"]');
      const pivotPanel = Array.from(root!.querySelectorAll('.grok-pivot-column-panel'))
        .find((p) => p.querySelector('.grok-pivot-column-tags-title')?.getAttribute('d4-name') === 'Pivot') as HTMLElement | undefined;
      const pivotVisible = !!pivotPanel && pivotPanel.offsetParent !== null;
      return {aggEmpty: pv.props.aggregateColumnNames.length === 0, pivotCleared: pv.props.pivotColumnNames.length === 0, pivotVisible};
    }), (s) => s.aggEmpty && s.pivotCleared && !s.pivotVisible, 3000, 150);
    expect(state.aggEmpty).toBe(true);
    expect(state.pivotCleared).toBe(true);   
    expect(state.pivotVisible).toBe(false);  
  });

  await softStep('Scenario 3 Step 5: re-adding an aggregate brings the Pivot row back', async () => {
    await addColumnViaPlus(page, 'div-add-Aggregate', 'AGE');
    const aggChips = await v.pollValue(() => rowChips(page, 'Aggregate'), (c) => c.length === 1, 3000, 150);
    expect(aggChips.length).toBe(1);
    expect(aggChips[0]).toContain('AGE');
    const pivotVisible = await v.pollValue(() => page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Pivot-table"]');
      const pivotPanel = Array.from(root!.querySelectorAll('.grok-pivot-column-panel'))
        .find((p) => p.querySelector('.grok-pivot-column-tags-title')?.getAttribute('d4-name') === 'Pivot') as HTMLElement | undefined;
      return !!pivotPanel && pivotPanel.offsetParent !== null;
    }), (visible) => visible, 3000, 150);
    expect(pivotVisible).toBe(true);
  });

  async function selectViewerObject() {
    await page.evaluate(() => {
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      grok.shell.o = pv;
    });
    await page.waitForSelector('[name="prop-aggregate"]', {timeout: 8000});
  }

  async function openPanelColumnDialog(row: string) {
    await selectViewerObject();
    await page.waitForSelector(`[name="prop-${row}"]`, {timeout: 6000});
    await page.locator(`[name="prop-view-${row}"] button`).click();      
    await page.waitForSelector('[name="dialog-Select-columns..."]', {timeout: 6000});
  }

  const panelCount = (row: string) => page.locator(`[name="prop-view-${row}"] label`).innerText();

  async function toggleFirstRow(col: string) {
    const dlg = page.locator('[name="dialog-Select-columns..."]');
    const search = dlg.locator('input.d4-search-input');
    await search.click();
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Delete');
    await page.keyboard.type(col);                 
    await page.waitForTimeout(400);                
    const rect = await dlg.locator('.d4-grid').boundingBox();
    if (!rect) throw new Error('Select-columns grid not visible');

    await page.mouse.click(rect.x + rect.width - 39, rect.y + 34);
    await page.waitForTimeout(400);                
  }

  await softStep('Scenario 4 Step 3: the viewer property panel exposes the group-by / aggregate / pivot column-list editors', async () => {

    await v.setViewerProps(page, 'Pivot table', [{set: {
      groupByColumnNames: ['DIS_POP'],
      aggregateColumnNames: ['AGE'],
      aggregateAggTypes: ['avg'],
      pivotColumnNames: ['SEVERITY'],
    }}], 700);
    await selectViewerObject();                 
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

    await toggleFirstRow('AGE');                                   
    expect(await panelCount('aggregate')).toBe('0 / 11');          
    await toggleFirstRow('WEIGHT');                                
    expect(await panelCount('aggregate')).toBe('1 / 11');          
    await dlg.locator('[name="button-OK"]').click();
    await v.waitForViewerRendered(page, 'Pivot table', 600);
    const agg = (await v.pollValue(() => pivotProps(page), (p) => p.agg.includes('WEIGHT'), 3000, 150)).agg;
    expect(agg).toContain('WEIGHT');
    expect(agg).not.toContain('AGE');
    expect(consoleErrors.length).toBe(before);                     

    await v.setViewerProps(page, 'Pivot table',
      [{set: {aggregateColumnNames: ['AGE'], aggregateAggTypes: ['avg']}}], 400);
  });

  await softStep('Scenario 4 Step 7: replacing DIS_POP with SEX via the group-by editor rewrites groupByColumnNames with no console error (GROK-16305)', async () => {
    const before = consoleErrors.length;
    await openPanelColumnDialog('group-by');
    const dlg = page.locator('[name="dialog-Select-columns..."]');

    await toggleFirstRow('DIS_POP');                               
    expect(await panelCount('group-by')).toBe('0 / 11');           
    await toggleFirstRow('SEX');                                   
    expect(await panelCount('group-by')).toBe('1 / 11');           
    await dlg.locator('[name="button-OK"]').click();
    await v.waitForViewerRendered(page, 'Pivot table', 600);
    const gb = (await v.pollValue(() => pivotProps(page), (p) => p.groupBy.includes('SEX'), 3000, 150)).groupBy;
    expect(gb).toContain('SEX');
    expect(gb).not.toContain('DIS_POP');
    expect(consoleErrors.length).toBe(before);                     

    await v.setViewerProps(page, 'Pivot table', [{set: {groupByColumnNames: ['DIS_POP']}}], 400);
  });

  await softStep('Scenario 5 Step 4: Refresh resets the visible configuration to type-driven defaults', async () => {

    await v.setViewerProps(page, 'Pivot table', [{set: {
      groupByColumnNames: ['DIS_POP'],
      aggregateColumnNames: ['AGE'],
      aggregateAggTypes: ['med'],
      pivotColumnNames: ['SEVERITY'],
    }}], 600);

    expect(await rowChips(page, 'Aggregate')).toContain('med(AGE)');
    expect((await rowChips(page, 'Group by')).length).toBeGreaterThan(0);

    await page.locator(`${PIVOT} .d4-command-bar [name="icon-redo"]`).click();
    await v.waitForViewerRendered(page, 'Pivot table', 900);

    const aggAfter = await v.pollValue(() => rowChips(page, 'Aggregate'), (c) => c.length === 2, 3000, 150);
    const gbAfter = await rowChips(page, 'Group by');
    const pivotAfter = await rowChips(page, 'Pivot');
    expect(gbAfter).toEqual([]);
    expect(pivotAfter).toEqual([]);
    expect(aggAfter).toEqual(['avg(AGE)', 'avg(HEIGHT)']);
  });

  await softStep('Scenario 5 Step 9: a non-default aggregation survives a layout saved through the View menu and re-applied via the layout API (github-2535)', async () => {

    await v.setViewerProps(page, 'Pivot table', [{set: {
      groupByColumnNames: ['DIS_POP'],
      aggregateColumnNames: ['AGE'],
      aggregateAggTypes: ['med'],
      pivotColumnNames: ['SEVERITY'],
    }}], 500);

    const beforeIds = await page.evaluate(async () =>
      (await grok.dapi.layouts.getApplicable(grok.shell.tv.dataFrame)).map((l: any) => l.id));

    const saveLeaf = await openTopMenuLeaf(page, 'View', 'Layout', 'Save-to-Gallery');
    const lb = await page.locator(`[name="${saveLeaf}"]`).boundingBox();
    await page.mouse.move(lb!.x + lb!.width / 2, lb!.y + lb!.height / 2);
    await page.waitForTimeout(120);   
    await page.mouse.down();
    await page.mouse.up();
    await v.pollValue(() => page.evaluate(async (prev) =>
      (await grok.dapi.layouts.getApplicable(grok.shell.tv.dataFrame)).some((l: any) => !prev.includes(l.id)),
    beforeIds), (fresh) => fresh, 1500, 150);

    const savedId = await page.evaluate(async (prev) => {
      const applicable = await grok.dapi.layouts.getApplicable(grok.shell.tv.dataFrame);
      const fresh = applicable.filter((l: any) => !prev.includes(l.id));
      return (fresh[0] ?? applicable[applicable.length - 1])?.id ?? null;
    }, beforeIds);
    expect(savedId).not.toBeNull();                      

    await v.setViewerProps(page, 'Pivot table',
      [{set: {aggregateAggTypes: ['avg'], groupByColumnNames: ['SEX']}}], 400);

    const result = await page.evaluate(async (id) => {
      const tv = grok.shell.tv;
      const saved = await grok.dapi.layouts.find(id);
      tv.loadLayout(saved);

      const t0 = Date.now();
      let pv2: any = null;
      for (;;) {
        pv2 = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
        if (pv2?.props.groupByColumnNames?.includes('DIS_POP') || Date.now() - t0 >= 3000) break;
        await new Promise((r) => setTimeout(r, 100));
      }
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
    expect(result.aggTypes).toContain('med');            
  });

  await softStep('Scenario 5 Step 13: the configuration survives a project save / close / reopen round-trip (github-2535)', async () => {
    const projectName = `PivotCrosstabProj${Date.now()}`;

    await v.setViewerProps(page, 'Pivot table', [{set: {
      groupByColumnNames: ['DIS_POP'],
      aggregateColumnNames: ['AGE'],
      aggregateAggTypes: ['med'],
      pivotColumnNames: ['SEVERITY'],
    }}], 500);

    const saved = await proj.saveProjectViaUI(page, projectName);
    expect(saved.projectId).toBeTruthy();                 

    await v.closeAllAndWait(page);
    const opened = await page.evaluate(async (id) => {
      const dash = await grok.dapi.projects.find(id);
      if (dash) await dash.open();

      const t0 = Date.now();
      while (Date.now() - t0 < 4000) {
        if (Array.from(grok.shell.tv?.viewers ?? []).some((x: any) => x.type === 'Pivot table')) break;
        await new Promise((r) => setTimeout(r, 100));
      }
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
    expect(result.aggTypes).toContain('med');             

    await proj.deleteProjectWithCleanup(page, {projectId: saved.projectId});
  });

  v.finishSpec();
});
