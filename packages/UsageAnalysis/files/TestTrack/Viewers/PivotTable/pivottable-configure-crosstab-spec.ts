/* ---
realizes: [pivottable.cp.configure-crosstab-values, pivottable.int.agg-types-track-agg-columns, pivottable.int.empty-aggregates-clear-pivot]
--- */
import {expect} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {
  PIVOT, openPivot, rowChips, pivotProps, ensurePivotPlusClickable, addColumnViaPlus, openChipMenu,
  expandSubmenu, pickAggregation, pickColumn, checkedAggregation, closePivotTagMenu, removeChip,
} from './pivot-helpers';

declare const grok: any;

// The layout and project round-trips (Scenario 5 Steps 9 and 13) live in pivot-table-server-spec.ts.
test.use(specTestOptions);

const DEFAULT_CROSSTAB = {
  groupByColumnNames: ['DIS_POP'], aggregateColumnNames: ['AGE'], aggregateAggTypes: ['avg'], pivotColumnNames: ['SEVERITY'],
};

test('Pivot Table — Configure cross-tab values', async ({page}) => {
  test.setTimeout(300_000);

  const consoleErrors: string[] = [];
  const ignorable = (m: string) => /Unable to find element in cloned iframe/i.test(m) || isLocalBootNoise(m);
  const onConsole = (m: any) => { if (m.type() === 'error' && !ignorable(m.text())) consoleErrors.push(m.text()); };
  const onPageError = (e: Error) => { if (!ignorable(e.message)) consoleErrors.push(e.message); };
  page.on('console', onConsole);
  page.on('pageerror', onPageError);

  await openDatagrok(page);
  await openPivot(page);

  const aggViewOpen = () => page.evaluate(() =>
    Array.from(grok.shell.views).some((vw: any) => vw.name === 'Table aggregation'));

  async function activatePivotView() {
    await page.evaluate(() => {
      const host = Array.from(grok.shell.views).find((vw: any) =>
        Array.from(vw.viewers ?? []).some((x: any) => x.type === 'Pivot table'));
      if (host) grok.shell.v = host;
    });
    await page.waitForSelector(`${PIVOT} [name="div-add-Group-by"]`, {state: 'visible', timeout: 8000});
    await ensurePivotPlusClickable(page, 'div-add-Group-by');
  }

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
    expect(await v.pollValue(aggViewOpen, (present) => present, 5000, 50)).toBe(true);

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

    await activatePivotView();
  });

  await softStep('Scenario 1 Step 7: opening the Group by column picker writes no console error (GROK-19114)', async () => {
    const before = consoleErrors.length;
    await ensurePivotPlusClickable(page, 'div-add-Group-by');
    await page.locator(`${PIVOT} [name="div-add-Group-by"]`).click();
    await page.waitForSelector('.d4-column-selector-backdrop', {timeout: 6000});
    await page.keyboard.press('Escape');
    await v.pollValue(() => page.locator('.d4-column-selector-backdrop').count(), (n) => n === 0, 400, 50);
    expect(consoleErrors.length).toBe(before);
  });

  await softStep('Scenario 1 Step 10: adding SEX in Group by yields a SEX chip and a DIS_POP+SEX grouping', async () => {
    await addColumnViaPlus(page, 'div-add-Group-by', 'SEX');
    const chips = await v.pollValue(() => rowChips(page, 'Group by'), (c) => c.includes('SEX'), 3000, 50);
    expect(chips).toEqual(expect.arrayContaining(['DIS_POP', 'SEX']));
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
    await v.pollValue(aggViewOpen, (present) => !present, 2000, 50);

    await addColumnViaPlus(page, 'div-add-Group-by', 'DIS_POP');
  });

  await softStep('Scenario 2 Step 3: right-click avg(AGE) opens the menu with avg checked, no console error (GROK-17841)', async () => {
    await v.setViewerProps(page, 'Pivot table', [{set: DEFAULT_CROSSTAB}], 500);
    const before = consoleErrors.length;
    await openChipMenu(page, 'Aggregate');

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
      (c) => c.includes('sum(AGE)'), 3000, 100)).toContain('sum(AGE)');
    expect((await pivotProps(page)).aggTypes).toEqual(['sum']);
  });

  await softStep('Scenario 2 Step 6: a second pick (Median) in the same open menu also takes effect (GROK-16899)', async () => {
    await pickAggregation(page, 'med');
    expect(await page.locator('.d4-menu-popup[name="pivot-tag"]').count()).toBeGreaterThan(0);
    await expandSubmenu(page, 'Aggregation', 'div-Aggregation---med');
    expect(await checkedAggregation(page)).toEqual(['med']);
    expect(await v.pollValue(() => rowChips(page, 'Aggregate'),
      (c) => c.includes('med(AGE)'), 3000, 100)).toContain('med(AGE)');
  });

  await softStep('Scenario 2 Step 8: switching column to HEIGHT rebuilds the Aggregation group to its supported types (I2)', async () => {
    await pickColumn(page, 'HEIGHT');
    expect(await page.locator('.d4-menu-popup[name="pivot-tag"]').count()).toBeGreaterThan(0);

    expect(await v.pollValue(() => rowChips(page, 'Aggregate'),
      (c) => c.includes('avg(HEIGHT)'), 3000, 100)).toContain('avg(HEIGHT)');

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
      (c) => c.includes('avg(WEIGHT)'), 3000, 100)).toContain('avg(WEIGHT)');

    await closePivotTagMenu(page);
    await openChipMenu(page, 'Aggregate');
    const before = consoleErrors.length;

    await page.locator('.d4-menu-popup[name="pivot-tag"] [name="div-Remove-others"]').click();
    await v.pollValue(() => page.locator('.d4-menu-popup[name="pivot-tag"]').count(), (n) => n === 0, 500, 50);
    expect(await page.locator('.d4-menu-popup[name="pivot-tag"]').count()).toBe(0);
    expect((await rowChips(page, 'Aggregate')).length).toBe(1);
    expect(consoleErrors.length).toBe(before);
  });

  const pivotRowVisible = () => page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Pivot-table"]');
    const pivotPanel = Array.from(root!.querySelectorAll('.grok-pivot-column-panel'))
      .find((p) => p.querySelector('.grok-pivot-column-tags-title')?.getAttribute('d4-name') === 'Pivot') as HTMLElement | undefined;
    return !!pivotPanel && pivotPanel.offsetParent !== null;
  });

  await softStep('Scenario 3 Step 3: removing the last aggregate hides the Pivot row and clears pivot columns (I4)', async () => {
    await v.setViewerProps(page, 'Pivot table', [{set: DEFAULT_CROSSTAB}], 500);
    await removeChip(page, 'Aggregate');

    const state = await v.pollValue(async () => {
      const p = await pivotProps(page);
      return {aggEmpty: p.agg.length === 0, pivotCleared: p.pivot.length === 0, pivotVisible: await pivotRowVisible()};
    }, (s) => s.aggEmpty && s.pivotCleared && !s.pivotVisible, 3000, 100);
    expect(state.aggEmpty).toBe(true);
    expect(state.pivotCleared).toBe(true);
    expect(state.pivotVisible).toBe(false);
  });

  await softStep('Scenario 3 Step 5: re-adding an aggregate brings the Pivot row back', async () => {
    await addColumnViaPlus(page, 'div-add-Aggregate', 'AGE');
    const aggChips = await v.pollValue(() => rowChips(page, 'Aggregate'), (c) => c.length === 1, 3000, 100);
    expect(aggChips.length).toBe(1);
    expect(aggChips[0]).toContain('AGE');
    expect(await v.pollValue(pivotRowVisible, (visible) => visible, 3000, 100)).toBe(true);
  });

  async function selectViewerObject() {
    await page.evaluate(() => {
      grok.shell.o = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table');
    });
    await page.waitForSelector('[name="prop-aggregate"]', {timeout: 8000});
  }

  const SELECT_DLG = '[name="dialog-Select-columns..."]';

  async function openPanelColumnDialog(row: string) {
    await selectViewerObject();
    await page.waitForSelector(`[name="prop-${row}"]`, {timeout: 6000});
    await page.locator(`[name="prop-view-${row}"] button`).click();
    await page.waitForSelector(SELECT_DLG, {timeout: 6000});
  }

  const panelCount = (row: string) => page.locator(`[name="prop-view-${row}"] label`).innerText();

  // a cheap stamp of the dialog grid's first rows, so the click waits for the search filter to repaint them
  const dialogGridStamp = () => page.evaluate((sel) => {
    const cv = document.querySelector(`${sel} .d4-grid canvas`) as HTMLCanvasElement | null;
    const ctx = cv?.getContext('2d');
    if (!cv || !ctx || cv.width === 0 || cv.height === 0) return -1;
    const data = ctx.getImageData(0, 0, cv.width, Math.min(cv.height, 80)).data;
    let sum = 0;
    for (let i = 0; i < data.length; i += 16) sum = (sum + data[i] * 31 + i) % 1e9;
    return sum;
  }, SELECT_DLG);

  async function toggleFirstRow(row: string, col: string) {
    const dlg = page.locator(SELECT_DLG);
    const search = dlg.locator('input.d4-search-input');
    await search.click();
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Delete');
    const stampBefore = await dialogGridStamp();
    await page.keyboard.type(col);
    await v.pollValue(dialogGridStamp, (s) => s !== stampBefore, 1000, 25);
    await v.pollStable(dialogGridStamp, (a, b) => a === b, 1000, 100);
    const rect = await dlg.locator('.d4-grid').boundingBox();
    if (!rect) throw new Error('Select-columns grid not visible');

    const countBefore = await panelCount(row);
    await page.mouse.click(rect.x + rect.width - 39, rect.y + 34);
    await v.pollValue(() => panelCount(row), (t) => t !== countBefore, 2000, 50);
  }

  await softStep('Scenario 4 Step 3: the viewer property panel exposes the group-by / aggregate / pivot column-list editors', async () => {
    await v.setViewerProps(page, 'Pivot table', [{set: DEFAULT_CROSSTAB}], 700);
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
    const dlg = page.locator(SELECT_DLG);

    await toggleFirstRow('aggregate', 'AGE');
    expect(await panelCount('aggregate')).toBe('0 / 11');
    await toggleFirstRow('aggregate', 'WEIGHT');
    expect(await panelCount('aggregate')).toBe('1 / 11');
    await dlg.locator('[name="button-OK"]').click();
    const agg = (await v.pollValue(() => pivotProps(page), (p) => p.agg.includes('WEIGHT'), 3000, 100)).agg;
    expect(agg).toContain('WEIGHT');
    expect(agg).not.toContain('AGE');
    expect(consoleErrors.length).toBe(before);

    await v.setViewerProps(page, 'Pivot table',
      [{set: {aggregateColumnNames: ['AGE'], aggregateAggTypes: ['avg']}}], 400);
  });

  await softStep('Scenario 4 Step 7: replacing DIS_POP with SEX via the group-by editor rewrites groupByColumnNames with no console error (GROK-16305)', async () => {
    const before = consoleErrors.length;
    await openPanelColumnDialog('group-by');
    const dlg = page.locator(SELECT_DLG);

    await toggleFirstRow('group-by', 'DIS_POP');
    expect(await panelCount('group-by')).toBe('0 / 11');
    await toggleFirstRow('group-by', 'SEX');
    expect(await panelCount('group-by')).toBe('1 / 11');
    await dlg.locator('[name="button-OK"]').click();
    const gb = (await v.pollValue(() => pivotProps(page), (p) => p.groupBy.includes('SEX'), 3000, 100)).groupBy;
    expect(gb).toContain('SEX');
    expect(gb).not.toContain('DIS_POP');
    expect(consoleErrors.length).toBe(before);

    await v.setViewerProps(page, 'Pivot table', [{set: {groupByColumnNames: ['DIS_POP']}}], 400);
  });

  await softStep('Scenario 5 Step 4: Refresh resets the visible configuration to type-driven defaults', async () => {
    await v.setViewerProps(page, 'Pivot table', [{set: {...DEFAULT_CROSSTAB, aggregateAggTypes: ['med']}}], 600);

    expect(await rowChips(page, 'Aggregate')).toContain('med(AGE)');
    expect((await rowChips(page, 'Group by')).length).toBeGreaterThan(0);

    await page.locator(`${PIVOT} .d4-command-bar [name="icon-redo"]`).click();

    const aggAfter = await v.pollValue(() => rowChips(page, 'Aggregate'), (c) => c.length === 2, 3000, 100);
    const gbAfter = await rowChips(page, 'Group by');
    const pivotAfter = await rowChips(page, 'Pivot');
    expect(gbAfter).toEqual([]);
    expect(pivotAfter).toEqual([]);
    expect(aggAfter).toEqual(['avg(AGE)', 'avg(HEIGHT)']);
  });

  await page.evaluate(() => { grok.shell.o = null; });
  page.off('console', onConsole);
  page.off('pageerror', onPageError);
  await v.cleanupShell(page);
  v.finishSpec();
});
