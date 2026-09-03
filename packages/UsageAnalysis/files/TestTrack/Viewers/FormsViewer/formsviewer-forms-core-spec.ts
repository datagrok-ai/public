/* ---
realizes: [formsviewer.cp.forms-core, formsviewer.int.sort-mirrors-grid, formsviewer.int.selection-intersects-filter, formsviewer.edge.pin-non-unique-value-warns, formsviewer.edge.pinned-row-absent-from-ordinary-cards]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {knownOpenBug} from '../../helpers/known-open-bug';
import {
  HOST, ORDINARY, CURRENT, PINNED, PINNED_PANE,
  cardFieldValue, cardIndexByValue, balloonCount, drawnLabelNames, waitForOrderStable,
  fieldValuesByPosition, sortIndicatorLabels, sortArrow, cardContextMenu,
} from '../../helpers/forms';

declare const grok: any;

// The viewer ladder of the forms-core scenario. The layout and project round-trips (Steps 7a-7c)
// live in formsviewer-forms-core-server-spec.ts so the dev round-trips stay out of this one.
// Both halves run on the server lane: Forms is a PowerGrid package viewer, and the local client
// serves no packages, so its toolbox has no Forms icon (measured 2026-09-03, icon-Forms null).
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

async function tailUsubjids(page: Page, offset: number): Promise<string[]> {
  return (await fieldValuesByPosition(page, 'USUBJID')).slice(offset).filter((x): x is string => x !== null);
}

async function ordinaryUsubjids(page: Page): Promise<string[]> {
  return (await fieldValuesByPosition(page, 'USUBJID')).filter((x): x is string => x !== null);
}

async function ordinaryHeights(page: Page, offset: number): Promise<number[]> {
  return page.evaluate(({sel, off}) => Array.from(document.querySelectorAll(sel)).slice(off)
    .map((c) => parseFloat(((c as HTMLElement).querySelector('[column="HEIGHT"]') as HTMLInputElement)?.value))
    .filter((x) => !Number.isNaN(x)), {sel: ORDINARY, off: offset});
}

function setForms(page: Page, options: Record<string, any>): Promise<void> {
  return page.evaluate((o) => {
    grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer').setOptions(o);
  }, options);
}

// The non-unique warning is raised inside the pin itself (Step 6c sees it on its first poll), so a
// unique pin has either shown a balloon within the hold or never will.
async function expectNoBalloon(page: Page, holdMs = 1000): Promise<void> {
  expect(await v.pollValue(() => balloonCount(page), (n) => n > 0, holdMs, 100)).toBe(0);
}

test('Forms viewer — core ladder (p0)', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await softStep('Step 1 — Add the Forms viewer; default field set is the first 20 visible columns', async () => {
    await v.addViewerByIcon(page, 'Forms', 'Forms', 30_000, 'FormsViewer');
    await page.locator('.d4-multi-form').first().waitFor({timeout: 30_000});

    const expectedFields = await page.evaluate(() =>
      grok.shell.t.columns.names().filter((n: string) => !n.startsWith('~')).slice(0, 20));

    await expect.poll(() => drawnLabelNames(page), {timeout: 20_000}).toEqual(expectedFields);
    expect(expectedFields.some((n: string) => n.startsWith('~'))).toBe(false);
    expect(await v.pollValue(() => balloonCount(page), (n) => n === 0, 4000, 250)).toBe(0);
  });

  await softStep('Step 2 — The current-row card shows the grid value and follows the current row', async () => {
    const startRow = 12;
    await page.evaluate((r) => { grok.shell.t.currentRowIdx = r; }, startRow);
    const gridStart = await page.evaluate((r) => grok.shell.tv.grid.cell('HEIGHT', r).cell.valueString, startRow);
    await expect.poll(() => cardFieldValue(page, 0, 'HEIGHT', CURRENT), {timeout: 10_000}).toBe(gridStart);

    await page.evaluate(() => { grok.shell.t.currentRowIdx = 77; });
    const grid77 = await page.evaluate(() => grok.shell.tv.grid.cell('HEIGHT', 77).cell.valueString);
    await expect.poll(() => cardFieldValue(page, 0, 'HEIGHT', CURRENT), {timeout: 10_000}).toBe(grid77);
  });

  await softStep('Step 3 — Show Selected Rows renders one card per selected row beyond the two leading cards', async () => {
    const defaultOn = await page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer').props.showSelectedRows);
    expect(defaultOn).toBe(true);

    const picked = await page.evaluate(() => {
      const df = grok.shell.t;
      df.currentRowIdx = 0;
      let f = -1; let m = -1; let m2 = -1;
      for (let i = 0; i < df.rowCount; i++) {
        const s = df.col('SEX').get(i);
        if (s === 'F' && f < 0) f = i;
        else if (s === 'M' && m < 0) m = i;
        else if (s === 'M' && m2 < 0 && i !== m) m2 = i;
        if (f >= 0 && m >= 0 && m2 >= 0) break;
      }
      df.selection.setAll(false);
      df.selection.set(f, true); df.selection.set(m, true); df.selection.set(m2, true);
      return {f, m, m2, usubjids: [f, m, m2].map((r) => df.col('USUBJID').get(r))};
    });

    await expect.poll(() => page.locator(ORDINARY).count(), {timeout: 15_000}).toBe(5);
    const extras = await tailUsubjids(page, 2);
    expect(extras.length).toBe(3);
    for (const usub of extras) expect(picked.usubjids).toContain(usub);
  });

  await softStep('Step 4 — The selected-row cards equal selection ∩ filter after a filter change', async () => {
    await v.applyCategoricalFilter(page, 'SEX', ['M']);

    await expect.poll(() => page.evaluate((sel) => {
      const df = grok.shell.t;
      const inter: string[] = [];
      let selCount = 0;
      for (let i = 0; i < df.rowCount; i++) {
        if (df.selection.get(i)) selCount++;
        if (df.selection.get(i) && df.filter.get(i)) inter.push(df.col('USUBJID').get(i));
      }
      const tail = (Array.from(document.querySelectorAll(sel)).slice(2)
        .map((c) => ((c as HTMLElement).querySelector('[column="USUBJID"]') as HTMLInputElement)?.value ?? null)
        .filter((x) => x !== null)) as string[];
      return JSON.stringify({
        excludedASelectedRow: inter.length < selCount,
        match: JSON.stringify(inter) === JSON.stringify(tail),
      });
    }, ORDINARY), {timeout: 20_000})
      .toBe(JSON.stringify({excludedASelectedRow: true, match: true}));

    await v.resetFilters(page);
    await expect.poll(async () => (await tailUsubjids(page, 2)).length, {timeout: 15_000}).toBe(3);
  });

  await softStep('Step 5a — Sorting the grid by HEIGHT mirrors the card order and marks the HEIGHT label', async () => {
    await page.evaluate(() => grok.shell.tv.grid.sort(['HEIGHT'], [true]));

    await expect.poll(() => page.evaluate((sel) => {
      const df = grok.shell.t;
      const selFilter = df.selection.clone().and(df.filter);
      const order = df.getSortedOrder(['HEIGHT'], [true], selFilter);
      const expected = Array.from(order).map((r: number) => df.col('USUBJID').get(r));
      const tail = (Array.from(document.querySelectorAll(sel)).slice(2)
        .map((c) => ((c as HTMLElement).querySelector('[column="USUBJID"]') as HTMLInputElement)?.value ?? null)
        .filter((x) => x !== null)) as string[];
      return JSON.stringify({nonEmpty: expected.length > 0, match: JSON.stringify(expected) === JSON.stringify(tail)});
    }, ORDINARY), {timeout: 20_000}).toBe(JSON.stringify({nonEmpty: true, match: true}));
    expect(await sortArrow(page, 'HEIGHT')).not.toBeNull();
  });

  await softStep('Step 5b — sortByColumnName overrides the grid sort; the indicator moves to WEIGHT', async () => {
    await setForms(page, {sortByColumnName: 'WEIGHT'});
    await expect.poll(() => sortIndicatorLabels(page), {timeout: 20_000}).toEqual(['div-WEIGHT']);

    const weightArrow = await sortArrow(page, 'WEIGHT');
    expect(weightArrow).not.toBeNull();
    await expect.poll(() => page.evaluate(({sel, asc}) => {
      const df = grok.shell.t;
      const selFilter = df.selection.clone().and(df.filter);
      const order = df.getSortedOrder(['WEIGHT'], [asc], selFilter);
      const expected = (Array.from(order).map((r: number) => df.col('USUBJID').get(r))) as string[];
      const tail = (Array.from(document.querySelectorAll(sel)).slice(2)
        .map((c) => ((c as HTMLElement).querySelector('[column="USUBJID"]') as HTMLInputElement)?.value ?? null)
        .filter((x) => x !== null)) as string[];
      return JSON.stringify({nonEmpty: expected.length > 0, match: JSON.stringify(expected) === JSON.stringify(tail)});
    }, {sel: ORDINARY, asc: weightArrow === '↑'}), {timeout: 20_000}).toBe(JSON.stringify({nonEmpty: true, match: true}));
    expect(await page.evaluate(() => grok.shell.tv.grid.sortByColumns.map((c: any) => c.name)))
      .toEqual(['HEIGHT']);
  });

  await softStep('Step 5c — Turning Use Grid Sort OFF stops mirroring the grid sort (GROK-20380 known-red)', async () => {
    await setForms(page, {sortByColumnName: null});

    await expect.poll(async () => {
      const heights = await ordinaryHeights(page, 2);
      return heights.length >= 2 && heights.every((h, i, a) => i === 0 || a[i - 1] <= h);
    }, {timeout: 20_000}).toBe(true);

    const tailBefore = JSON.stringify(await ordinaryUsubjids(page));
    await v.ensurePropertyCategory(page, 'Forms', 'misc', 'use-grid-sort');
    await v.setPropertyGridCheckbox(page, 'use-grid-sort', false, 'misc');
    // the cards either leave the grid order or (GROK-20380) keep it; a change is waited for, not slept for
    await v.pollValue(async () => JSON.stringify(await ordinaryUsubjids(page)), (t) => t !== tailBefore, 2000, 100);

    const mirror = JSON.parse(await page.evaluate((sel) => {
      const df = grok.shell.t;
      const grid = grok.shell.tv.grid;
      const selFilter = df.selection.clone().and(df.filter);
      const order = df.getSortedOrder(grid.sortByColumns.map((c: any) => c.name), grid.sortTypes, selFilter);
      const mirrored = (Array.from(order).map((r: number) => df.col('USUBJID').get(r))) as string[];
      const tail = (Array.from(document.querySelectorAll(sel)).slice(2)
        .map((c) => ((c as HTMLElement).querySelector('[column="USUBJID"]') as HTMLInputElement)?.value ?? null)
        .filter((x) => x !== null)) as string[];
      return JSON.stringify({expectedLen: mirrored.length, tailLen: tail.length,
        mirrors: JSON.stringify(mirrored) === JSON.stringify(tail)});
    }, ORDINARY)) as {expectedLen: number; tailLen: number; mirrors: boolean};

    expect(mirror.expectedLen).toBeGreaterThan(0);
    expect(mirror.tailLen).toBeGreaterThan(0);
    await knownOpenBug('GROK-20380', () => { expect(mirror.mirrors).toBe(false); });
  });

  await softStep('Step 5d — Double-clicking the sort label cycles the indicator; a different label does not move it', async () => {
    await v.ensurePropertyCategory(page, 'Forms', 'misc', 'use-grid-sort');
    await v.setPropertyGridCheckbox(page, 'use-grid-sort', true, 'misc');
    await setForms(page, {sortByColumnName: 'AGE'});
    await expect.poll(() => sortArrow(page, 'AGE'), {timeout: 20_000}).not.toBeNull();

    const ageLabel = page.locator(`${HOST} .d4-multi-form-header [name="div-AGE"]`).first();

    const seq: (string | null)[] = [];
    for (let i = 0; i < 3; i++) {
      const before = await sortArrow(page, 'AGE');
      await ageLabel.dblclick();
      seq.push(await v.pollValue(() => sortArrow(page, 'AGE'), (a) => a !== before, 1500, 50));
    }
    expect(new Set(seq).size).toBe(3);
    expect(seq).toContain(null);
    expect(seq.filter((s) => s !== null).length).toBe(2);

    await setForms(page, {sortByColumnName: null});
    await v.pollValue(() => sortArrow(page, 'AGE'), (a) => a === null, 1500, 100);
    await ageLabel.dblclick();
    await expect.poll(() => sortArrow(page, 'AGE'), {timeout: 15_000}).toBe('↓');
    await expect.poll(() => sortIndicatorLabels(page), {timeout: 15_000}).toEqual(['div-AGE']);

    await page.locator(`${HOST} .d4-multi-form-header [name="div-HEIGHT"]`).first().dblclick();
    const labels = await v.pollValue(() => sortIndicatorLabels(page),
      (l) => JSON.stringify(l) !== JSON.stringify(['div-AGE']), 1000, 100);
    expect(labels).toEqual(['div-AGE']);
  });

  await softStep('Step 6a — Pin Row moves the card to the pinned pane and removes it from the ordinary set', async () => {
    await v.ensurePropertyCategory(page, 'Forms', 'misc', 'show-mouse-over-row');
    await v.setPropertyGridCheckbox(page, 'show-mouse-over-row', false, 'misc');
    await setForms(page, {sortByColumnName: null});
    await page.evaluate(() => {
      const df = grok.shell.t;
      df.mouseOverRowIdx = -1;
      df.currentRowIdx = 0;
      df.selection.setAll(false);
      df.selection.set(5, true); df.selection.set(10, true); df.selection.set(20, true);
    });

    await expect.poll(() => page.locator(ORDINARY).count(), {timeout: 15_000}).toBe(4);
    await waitForOrderStable(page);

    const byPos6a = await fieldValuesByPosition(page, 'USUBJID');
    const targetPos6a = byPos6a.findIndex((val, i) => i >= 1 && val !== null);
    expect(targetPos6a).toBeGreaterThanOrEqual(0);
    const beforeCount = await page.locator(ORDINARY).count();
    const anchor = await cardFieldValue(page, targetPos6a, 'USUBJID');
    expect(anchor).not.toBeNull();

    await cardContextMenu(page, ORDINARY, targetPos6a, 'div-Pin-Row', 'USUBJID');
    await expectNoBalloon(page);

    await expect.poll(async () =>
      page.evaluate((sel) => getComputedStyle(document.querySelector(sel) as HTMLElement).display, PINNED_PANE),
    {timeout: 15_000}).not.toBe('none');
    const pinnedValues = await page.evaluate((sel) => Array.from(document.querySelectorAll(sel))
      .map((c) => ((c as HTMLElement).querySelector('[column="USUBJID"]') as HTMLInputElement)?.value), PINNED);
    expect(pinnedValues).toEqual([anchor]);

    await expect.poll(() => page.locator(ORDINARY).count(), {timeout: 15_000}).toBe(beforeCount - 1);
    expect(await ordinaryUsubjids(page)).not.toContain(anchor);
    expect(await page.evaluate(() =>
      grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer').props.pinnedRowValues)).toEqual([anchor]);

    expect(await page.evaluate((usub) => {
      const df = grok.shell.t;
      for (let i = 0; i < df.rowCount; i++)
        if (df.col('USUBJID').get(i) === usub) return df.selection.get(i);
      return false;
    }, anchor)).toBe(true);
  });

  await softStep('Step 6b — Unpin Row returns the row to the ordinary set and hides the pinned pane', async () => {
    const anchor = await page.evaluate((sel) =>
      ((document.querySelector(sel) as HTMLElement)?.querySelector('[column="USUBJID"]') as HTMLInputElement)?.value,
    PINNED);
    const beforeCount = await page.locator(ORDINARY).count();

    await cardContextMenu(page, PINNED, 0, 'div-Unpin-Row');

    await expect.poll(async () =>
      page.evaluate((sel) => getComputedStyle(document.querySelector(sel) as HTMLElement).display, PINNED_PANE),
    {timeout: 15_000}).toBe('none');
    await expect.poll(() => page.locator(ORDINARY).count(), {timeout: 15_000}).toBe(beforeCount + 1);
    expect(await ordinaryUsubjids(page)).toContain(anchor);

    await page.evaluate(() => {
      const df = grok.shell.t;
      df.selection.setAll(false);
      df.selection.set(5, true); df.selection.set(10, true); df.selection.set(20, true);
    });
    await expect.poll(() => page.locator(ORDINARY).count(), {timeout: 15_000}).toBe(4);

    const byPos6b = await fieldValuesByPosition(page, 'USUBJID');
    const rePinPos = byPos6b.findIndex((val, i) => i >= 1 && val !== null);
    expect(rePinPos).toBeGreaterThanOrEqual(0);
    await cardContextMenu(page, ORDINARY, rePinPos, 'div-Pin-Row', 'USUBJID');
    await expectNoBalloon(page);
    await expect.poll(() => page.evaluate((sel) => Array.from(document.querySelectorAll(sel)).length, PINNED),
      {timeout: 15_000}).toBe(1);
  });

  await softStep('Step 6c — Pinning through a NON-UNIQUE field raises the exact warning; the single pin is preserved', async () => {
    const pinnedBefore = await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer').props.pinnedRowValues as string[]));
    expect(pinnedBefore.length).toBe(1);

    const target = await page.evaluate(() => {
      const df = grok.shell.t;
      const sex = df.col('SEX');
      const firstOf: Record<string, number> = {};
      for (let i = 0; i < df.rowCount; i++) {
        const val = String(sex.get(i));
        if (!(val in firstOf)) firstOf[val] = i;
      }
      let row = -1;
      for (let i = 0; i < df.rowCount; i++)
        if (firstOf[String(sex.get(i))] !== i) { row = i; break; }
      df.currentRowIdx = 0;
      if (row >= 0) df.selection.set(row, true);
      return {row, usub: df.col('USUBJID').get(row), sex: String(sex.get(row))};
    });
    expect(target.row).toBeGreaterThanOrEqual(0);

    await waitForOrderStable(page);
    const cardIdx = await cardIndexByValue(page, 'USUBJID', target.usub);
    expect(cardIdx).toBeGreaterThanOrEqual(0);

    await cardContextMenu(page, ORDINARY, cardIdx, 'div-Pin-Row', 'SEX');

    await expect.poll(() => page.evaluate(() =>
      document.querySelector('.d4-balloon.warning .d4-balloon-content')?.textContent ?? null),
    {timeout: 15_000}).toBe("You have pinned a non-unique value. It won't be applied from the layout.");

    expect(await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer').props.pinnedRowValues as string[])))
      .toContain(target.sex);

    const pinnedSexIdx = await cardIndexByValue(page, 'USUBJID', target.usub, PINNED);
    expect(pinnedSexIdx).toBeGreaterThanOrEqual(0);
    await cardContextMenu(page, PINNED, pinnedSexIdx, 'div-Unpin-Row');
    await expect.poll(() => page.evaluate(() =>
      Array.from(grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer').props.pinnedRowValues as string[])),
    {timeout: 15_000}).toEqual(pinnedBefore);
  });

  await v.cleanupShell(page);
  v.finishSpec();
});
