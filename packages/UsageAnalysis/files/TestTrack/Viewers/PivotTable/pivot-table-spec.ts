/* ---
realizes: [pivottable.cp.chrome-history-and-drag-config, pivottable.int.history-menu-requires-existing-columns, pivottable.int.default-aggr-type-remembered]
--- */
import {expect} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {PIVOT, openPivot, rowChips, addColumnViaPlus, openChipMenu, pickAggregation, removeChip, ensurePivotPlusClickable} from './pivot-helpers';

declare const grok: any;

test.use(specTestOptions);

const HISTORY_KEY = 'grok-aggregation-history';
const isIgnorable = (m: string) => m.includes('cloned iframe') || isLocalBootNoise(m);

const pivotChips = (page: any) => Promise.all([rowChips(page, 'Group by'), rowChips(page, 'Aggregate'), rowChips(page, 'Pivot')])
  .then(([groupBy, agg, pivot]) => ({groupBy, agg, pivot}));

test('Pivot table chrome, history and drag-driven configuration', async ({page}) => {
  test.setTimeout(300_000);
  const pageErrors: string[] = [];
  const onPageError = (e: Error) => { if (!isIgnorable(e.message)) pageErrors.push(e.message); };
  page.on('pageerror', onPageError);

  await openDatagrok(page);
  await page.evaluate((k) => window.localStorage.removeItem(k), HISTORY_KEY);
  await openPivot(page);

  const hasPivot = () => page.evaluate(() => Array.from(grok.shell.tv.viewers).some((x: any) => x.type === 'Pivot table'));

  await softStep('Scenario 1 Step 2: auto-config is DIS_POP / SEVERITY / avg(AGE), counts visible', async () => {
    const r = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Pivot-table"]')!;
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      const counts = root.querySelector('.grok-pivot-counts') as HTMLElement | null;
      return {
        groupBy: pv.props.groupByColumnNames, pivot: pv.props.pivotColumnNames,
        agg: pv.props.aggregateColumnNames, aggTypes: pv.props.aggregateAggTypes,
        countsVisible: !!counts && !!counts.offsetParent,
      };
    });
    const chips = await pivotChips(page);

    expect(chips.groupBy).toEqual(['DIS_POP']);
    expect(chips.agg).toEqual(['avg(AGE)']);
    expect(chips.pivot).toEqual(['SEVERITY']);

    expect(r.groupBy).toContain('DIS_POP');
    expect(r.pivot).toContain('SEVERITY');
    expect(r.agg).toContain('AGE');
    expect(r.aggTypes).toContain('avg');
    expect(r.countsVisible).toBe(true);
  });

  await softStep('Scenario 2 Step 3: close via cross icon → viewer gone, no console-error delta (GROK-17122)', async () => {
    const errorsBefore = pageErrors.length;
    await page.evaluate(() => {
      const panel = document.querySelector('[name="viewer-Pivot-table"]')!.closest('.panel-base')!;
      (panel.querySelector('.panel-titlebar [name="Close"]') as HTMLElement)?.click();
    });
    const gone = !await v.pollValue(hasPivot, (has) => !has, 3000, 50);
    expect(gone).toBe(true);

    await v.addViewerByIcon(page, 'pivot-table', 'Pivot-table', 15000);
    const r = await page.evaluate(() => {
      const pv2 = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      return {reAdded: !!pv2, groupBy: pv2?.props.groupByColumnNames, pivot: pv2?.props.pivotColumnNames};
    });
    expect(pageErrors.length).toBe(errorsBefore);
    expect(r.reAdded).toBe(true);
    expect(r.groupBy).toContain('DIS_POP');
    expect(r.pivot).toContain('SEVERITY');
  });

  const probeChrome = async () => page.evaluate(() => {
    const root = document.querySelector('[name="viewer-Pivot-table"]')!;
    const vis = (el: Element | null) => !!el && !!(el as HTMLElement).offsetParent && getComputedStyle(el as HTMLElement).display !== 'none';
    const rowByTitle = (title: string) => [...root.querySelectorAll('.grok-pivot-column-panel')]
      .find((p) => p.querySelector(`.grok-pivot-column-tags-title[d4-name="${title}"]`)) ?? null;
    return {
      data: vis(rowByTitle('Data')), groupBy: vis(rowByTitle('Group by')),
      agg: vis(rowByTitle('Aggregate')), pivot: vis(rowByTitle('Pivot')),
      counts: vis(root.querySelector('.grok-pivot-counts')),
      cmdBar: vis(root.querySelector('.d4-command-bar')),
      history: vis(root.querySelector('.d4-command-bar [name="icon-history"]')),
    };
  });
  const setChromeProp = async (prop: string, value: boolean) =>
    v.setViewerProps(page, 'Pivot table', [{set: {[prop]: value}, wait: 450}]);

  await softStep('Scenario 3 Step 3: Show Header=false hides the Data row, tag rows and counts; they return on true', async () => {
    await setChromeProp('showHeader', false);
    const headerOff = await probeChrome();
    await setChromeProp('showHeader', true);
    const headerOn = await probeChrome();
    expect(headerOff.data).toBe(false);
    expect(headerOff.groupBy).toBe(false);
    expect(headerOff.agg).toBe(false);
    expect(headerOff.pivot).toBe(false);
    expect(headerOff.counts).toBe(false);
    expect(headerOn.groupBy).toBe(true);
    expect(headerOn.counts).toBe(true);
  });

  await softStep('Scenario 3 Step 6: Show Command Bar=false hides the command bar with history/refresh icons; it returns on true', async () => {
    await setChromeProp('showCommandBar', false);
    const cmdOff = await probeChrome();
    await setChromeProp('showCommandBar', true);
    const cmdOn = await probeChrome();
    expect(cmdOff.cmdBar).toBe(false);
    expect(cmdOff.history).toBe(false);
    expect(cmdOn.cmdBar).toBe(true);
    expect(cmdOn.history).toBe(true);
  });

  await softStep('Scenario 4 Step 4: title in the header, description Top visible, Never hides it', async () => {
    await v.setViewerProps(page, 'Pivot table', [{set: {
      showTitle: true, title: 'My Pivot', description: 'Summary stats', descriptionPosition: 'Top',
    }, wait: 500}]);
    const top = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Pivot-table"]')!;
      const panel = root.closest('.panel-base') ?? root;
      const titleText = [...document.querySelectorAll('.panel-titlebar-tabhost .panel-titlebar-text')]
        .map((e) => e.textContent!.trim()).filter(Boolean);
      const descTop = [...panel.querySelectorAll('.d4-viewer-description')]
        .map((e) => ({txt: e.textContent!.trim(), vis: !!(e as HTMLElement).offsetParent}));
      return {titleShown: titleText.includes('My Pivot'), descTopVisible: descTop.some((d) => d.txt.includes('Summary stats') && d.vis)};
    });
    await v.setViewerProps(page, 'Pivot table', [{set: {descriptionVisibilityMode: 'Never'}, wait: 500}]);
    const descNeverCount = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Pivot-table"]')!;
      const panel = root.closest('.panel-base') ?? root;
      return [...panel.querySelectorAll('.d4-viewer-description')]
        .filter((e) => !!(e as HTMLElement).offsetParent && e.textContent!.includes('Summary stats')).length;
    });

    await page.evaluate(() => {
      const pv = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
      pv.props.showTitle = false; pv.props.title = ''; pv.props.description = ''; pv.props.descriptionVisibilityMode = 'Auto';
    });
    expect(top.titleShown).toBe(true);
    expect(top.descTopVisible).toBe(true);
    expect(descNeverCount).toBe(0);
  });

  const openHistoryMenu = async () => {
    await page.evaluate(() =>
      (document.querySelector('[name="viewer-Pivot-table"] .d4-command-bar [name="icon-history"]') as HTMLElement).click());
    await page.locator('.d4-menu-popup .d4-menu-item-label').first().waitFor({timeout: 5000});
  };

  await softStep('Scenario 5 Step 4: Save parameters writes localStorage history for RACE / avg(WEIGHT)', async () => {
    await v.setViewerProps(page, 'Pivot table', [{set: {
      groupByColumnNames: ['RACE'], aggregateColumnNames: ['WEIGHT'], aggregateAggTypes: ['avg'], pivotColumnNames: [],
    }, wait: 600}]);
    await openHistoryMenu();
    await page.evaluate(() => {
      const saveItem = [...document.querySelectorAll('.d4-menu-item')]
        .find((i) => i.querySelector('.d4-menu-item-label')?.textContent?.trim() === 'Save parameters') as HTMLElement | null;
      saveItem?.click();
    });
    const r = await v.pollValue(() => page.evaluate((k) => {
      const raw = window.localStorage.getItem(k);
      let parsed: any = null;
      try { parsed = JSON.parse(raw ?? ''); } catch (_) { parsed = null; }
      const flat = Array.isArray(parsed) ? parsed.flat(2).map((a: any) => a.colName) : [];
      return {isArray: Array.isArray(parsed), len: Array.isArray(parsed) ? parsed.length : -1, names: flat};
    }, HISTORY_KEY), (x) => x.len > 0, 3000, 50);
    expect(r.isArray).toBe(true);
    expect(r.len).toBeGreaterThan(0);
    expect(r.names).toContain('RACE');
    expect(r.names).toContain('WEIGHT');
  });

  await softStep('Scenario 5 Step 6: picking the saved entry restores Group by / Aggregate (tag captions)', async () => {
    await v.setViewerProps(page, 'Pivot table', [{set: {
      groupByColumnNames: ['SEX'], aggregateColumnNames: ['AGE'], aggregateAggTypes: ['avg'],
    }, wait: 500}]);
    await openHistoryMenu();
    await page.locator('.d4-menu-popup [name="div-key(RACE),avg(WEIGHT)"]').click();
    const tags = await v.pollValue(() => pivotChips(page), (t) => t.groupBy.includes('RACE'), 3000, 50);
    expect(tags.groupBy).toEqual(['RACE']);
    expect(tags.agg).toEqual(['avg(WEIGHT)']);
  });

  await softStep('Scenario 5 Step 8: after WEIGHT is removed the history menu drops the WEIGHT entry (I8)', async () => {
    await page.evaluate(() => {
      grok.shell.tv.dataFrame.columns.remove('WEIGHT');
      (Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any).close();
    });
    expect(await v.pollValue(hasPivot, (has) => !has, 3000, 50)).toBe(false);
    await v.addViewerByIcon(page, 'pivot-table', 'Pivot-table', 15000);
    await page.locator(`${PIVOT} .grok-pivot-column-tags-title[d4-name="Group by"]`).waitFor({timeout: 15000});
    await openHistoryMenu();
    const labels = await page.evaluate(() => {
      const menu = [...document.querySelectorAll('.d4-menu-popup')].pop()!;
      const out = [...menu.querySelectorAll('.d4-menu-item-label')].map((e) => e.textContent!.trim());
      document.body.click();
      return out;
    });
    expect(labels.some((l) => l.includes('WEIGHT'))).toBe(false);
  });

  await softStep('Scenario 5 Step 9: Refresh (icon-redo) clears Group by / Pivot and re-seeds the default aggregates', async () => {
    await v.setViewerProps(page, 'Pivot table', [{set: {
      groupByColumnNames: ['RACE'], aggregateColumnNames: ['AGE'], aggregateAggTypes: ['sum'], pivotColumnNames: ['SEX'],
    }, wait: 600}]);
    await page.locator(`${PIVOT} .d4-command-bar [name="icon-redo"]`).click();
    const tags = await v.pollValue(() => pivotChips(page), (t) => t.groupBy.length === 0 && t.agg.length === 2, 3000, 50);
    expect(tags.groupBy).toEqual([]);
    expect(tags.pivot).toEqual([]);
    expect(tags.agg).toEqual(['avg(AGE)', 'avg(HEIGHT)']);
  });

  await page.evaluate((k) => window.localStorage.removeItem(k), HISTORY_KEY);
  await openPivot(page);

  await softStep('Scenario 6 Steps 1-2: add HEIGHT, choose sum → sum(HEIGHT) tag, then remove it', async () => {
    await addColumnViaPlus(page, 'div-add-Aggregate', 'HEIGHT');
    let aggChips = await v.pollValue(() => rowChips(page, 'Aggregate'), (c) => c.some((x) => x.includes('HEIGHT')), 3000, 50);
    expect(aggChips.some((c) => c.includes('HEIGHT'))).toBe(true);

    await openChipMenu(page, 'Aggregate', 'HEIGHT');
    await pickAggregation(page, 'sum');
    await page.keyboard.press('Escape');
    aggChips = await v.pollValue(() => rowChips(page, 'Aggregate'), (c) => c.some((x) => x.includes('sum(HEIGHT)')), 3000, 50);
    expect(aggChips.some((c) => c.includes('sum(HEIGHT)'))).toBe(true);

    await removeChip(page, 'Aggregate', 'HEIGHT');
    aggChips = await v.pollValue(() => rowChips(page, 'Aggregate'), (c) => !c.some((x) => x.includes('HEIGHT')), 3000, 50);
    expect(aggChips.some((c) => c.includes('HEIGHT'))).toBe(false);
  });

  await softStep('Scenario 6 Step 3: the Aggregate + popup pre-offers the remembered aggregation type (I9)', async () => {
    await ensurePivotPlusClickable(page, 'div-add-Aggregate');
    await page.locator(`${PIVOT} [name="div-add-Aggregate"]`).click();
    const backdrop = page.locator('.d4-column-selector-backdrop');
    await backdrop.waitFor({timeout: 6000});
    expect(await backdrop.count()).toBeGreaterThan(0);

    await page.keyboard.press('Escape');
    await v.pollValue(() => backdrop.count(), (n) => n === 0, 400, 50);
    const aggChips = await rowChips(page, 'Aggregate');
    expect(aggChips.some((c) => c.includes('AGE'))).toBe(true);
    expect(aggChips.some((c) => c.includes('HEIGHT'))).toBe(false);
  });

  await softStep('Scenario 8 Step 2: grouping by USUBJID makes one row per identifier, no console error (GROK-16201)', async () => {
    const errorsBefore = pageErrors.length;
    await v.setViewerProps(page, 'Pivot table', [{set: {
      groupByColumnNames: ['USUBJID'], aggregateColumnNames: ['AGE'], aggregateAggTypes: ['avg'], pivotColumnNames: [],
    }, wait: 800}]);
    const r = await v.pollValue(() => page.evaluate(() => {
      const counts = document.querySelector('[name="viewer-Pivot-table"] .grok-pivot-counts')!.textContent!.replace(/\s+/g, ' ').trim();
      const distinct = grok.shell.tv.dataFrame.col('USUBJID').categories.length;
      return {distinct, counts, rowsMatch: counts.startsWith(`${distinct} rows`)};
    }), (x) => x.rowsMatch, 3000, 50);
    expect(r.rowsMatch).toBe(true);
    expect(pageErrors.length).toBe(errorsBefore);
  });

  await softStep('Scenario 8 Step 5: ADD opens the aggregated result; key column keeps its type (GROK-16074)', async () => {
    await v.setViewerProps(page, 'Pivot table', [{set: {
      groupByColumnNames: ['DIS_POP'], aggregateColumnNames: ['AGE'], aggregateAggTypes: ['avg'], pivotColumnNames: [],
    }, wait: 600}]);
    const aggViewOpen = () => page.evaluate(() =>
      Array.from(grok.shell.views).some((vw: any) => vw.name === 'Table aggregation'));
    await page.locator(`${PIVOT} .grok-pivot-counts [name="button-ADD"]`).click();
    expect(await v.pollValue(aggViewOpen, (open) => open, 5000, 50)).toBe(true);
    const r = await page.evaluate(() => {
      const srcCol = grok.shell.tables.find((t: any) => t.rowCount === 5850).col('DIS_POP');
      const aggView = Array.from(grok.shell.views).find((vw: any) => vw.name === 'Table aggregation') as any;
      const keyCol = aggView?.dataFrame?.col('DIS_POP');
      const out = {
        opened: !!keyCol,
        keyType: keyCol?.type, srcType: srcCol.type,
        keySemType: keyCol?.semType ?? null, srcSemType: srcCol.semType ?? null,
      };
      aggView?.close();
      return out;
    });
    await v.pollValue(aggViewOpen, (open) => !open, 2000, 50);
    expect(r.opened).toBe(true);
    expect(r.keyType).toBe(r.srcType);
    expect(r.keySemType).toBe(r.srcSemType);
  });

  await softStep('Scenario 8 Step 7: switching the Data-row Table property back and forth duplicates no Data entry / header (github-3414, GROK-14995)', async () => {
    const counts = () => page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Pivot-table"]')!;
      const dataRow = [...root.querySelectorAll('.grok-pivot-column-panel')]
        .find((p) => p.querySelector('.grok-pivot-column-tags-title[d4-name="Data"]'));
      return {
        dataEntries: dataRow ? dataRow.querySelectorAll('.d4-tag').length : 0,
        headers: root.querySelectorAll('.grok-pivot-column-tags-title[d4-name="Data"]').length,
      };
    });
    const before = await counts();

    const tagBox = await page.evaluate(() => {
      const root = document.querySelector('[name="viewer-Pivot-table"]')!;
      const dataRow = [...root.querySelectorAll('.grok-pivot-column-panel')]
        .find((p) => p.querySelector('.grok-pivot-column-tags-title[d4-name="Data"]'))!;
      const r = (dataRow.querySelector('.d4-tag') as HTMLElement).getBoundingClientRect();
      return {x: r.x + r.width / 2, y: r.y + r.height / 2};
    });
    await page.mouse.click(tagBox.x, tagBox.y);
    const dlg = page.locator('.d4-dialog').last();
    await dlg.waitFor({timeout: 8000});

    await dlg.locator('[name="button-OK"]').click();
    await dlg.waitFor({state: 'detached', timeout: 5000}).catch(() => {});
    // the duplication would land within the old 700ms hold; the poll gives up at the same cap
    const after = await v.pollValue(counts, (c) => c.dataEntries !== before.dataEntries || c.headers !== 1, 700, 50);
    expect(after.dataEntries).toBe(before.dataEntries);
    expect(after.headers).toBe(1);
  });

  await page.evaluate((k) => window.localStorage.removeItem(k), HISTORY_KEY);
  page.off('pageerror', onPageError);
  await v.closeAllAndWait(page);
  v.finishSpec();
});
