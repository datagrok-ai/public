/* ---
realizes: [formsviewer.cp.forms-core, formsviewer.int.sort-mirrors-grid, formsviewer.int.pinned-rows-persist-by-value, formsviewer.int.selection-intersects-filter, formsviewer.edge.pin-non-unique-value-warns, formsviewer.edge.pinned-row-absent-from-ordinary-cards]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {knownOpenBug} from '../../helpers/known-open-bug';
import * as projects from '../../helpers/projects';
import {
  HOST, ORDINARY, CURRENT, PINNED, PINNED_PANE,
  cardFieldValue, cardIndexByValue, balloonCount, drawnLabelNames, waitForOrderStable,
  fieldValuesByPosition, withConsoleErrorCount,
} from '../../helpers/forms';

declare const grok: any;

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

async function sortIndicatorLabels(page: Page): Promise<string[]> {
  return page.evaluate((host) => Array.from(document.querySelectorAll(`${host} .d4-multi-form-header .d4-multi-form-column-name`))
    .filter((l) => l.querySelector('.d4-multi-form-column-sort-indicator'))
    .map((l) => (l.querySelector('div[name^="div-"]') as HTMLElement)?.getAttribute('name') ?? '')
    .filter((n) => n.length > 0), HOST);
}

async function sortArrow(page: Page, column: string): Promise<string | null> {
  return page.evaluate(({host, col}) => {
    const label = Array.from(document.querySelectorAll(`${host} .d4-multi-form-header .d4-multi-form-column-name`))
      .find((l) => l.querySelector(`div[name="div-${col}"]`));
    const ind = label?.querySelector('.d4-multi-form-column-sort-indicator');
    return ind ? (ind.textContent ?? '').trim() : null;
  }, {host: HOST, col: column});
}

async function cardContextMenu(
  page: Page, cardSelector: string, cardIndex: number, itemName: string, columnField?: string,
): Promise<void> {
  const card = page.locator(cardSelector).nth(cardIndex);
  const target = columnField ? card.locator(`[column="${columnField}"]`).first() : card;
  await target.click({button: 'right'});
  await page.locator(`[name="${itemName}"]`).first().waitFor({timeout: 5000});
  await page.locator(`[name="${itemName}"]`).first().click();
}

async function expectNoBalloonSustained(page: Page, windowMs = 15_000): Promise<void> {
  const deadline = Date.now() + windowMs;
  let seen = 0;
  do {
    seen = await balloonCount(page);
    if (seen > 0) break;
    await page.waitForTimeout(500); 
  } while (Date.now() < deadline);
  expect(seen).toBe(0);
}

test('Forms viewer — core ladder (p0)', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await softStep('Step 1 — Add the Forms viewer; default field set is the first 20 visible columns', async () => {
    await v.addViewerByIcon(page, 'Forms', 'Forms', 30_000, 'FormsViewer');
    await page.locator('.d4-multi-form').first().waitFor({timeout: 30_000});

    const expectedFields = await page.evaluate(() =>
      grok.shell.t.columns.names().filter((n: string) => !n.startsWith('~')).slice(0, 20));

    const drawn = await drawnLabelNames(page);
    expect(drawn).toEqual(expectedFields);
    expect(drawn.some((n) => n.startsWith('~'))).toBe(false);
    expect(await balloonCount(page)).toBe(0);
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

    const defaultOn = await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      return vw.props.showSelectedRows;
    });
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
    await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      vw.setOptions({sortByColumnName: 'WEIGHT'});
    });
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

    await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      vw.setOptions({sortByColumnName: null});
    });

    await expect.poll(async () => {
      const heights = await ordinaryHeights(page, 2);
      return heights.length >= 2 && heights.every((h, i, a) => i === 0 || a[i - 1] <= h);
    }, {timeout: 20_000}).toBe(true);

    await v.ensurePropertyCategory(page, 'Forms', 'misc', 'use-grid-sort');
    await v.setPropertyGridCheckbox(page, 'use-grid-sort', false, 'misc');

    await page.waitForTimeout(3500);

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
    await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      vw.setOptions({sortByColumnName: 'AGE'});
    });
    await expect.poll(() => sortArrow(page, 'AGE'), {timeout: 20_000}).not.toBeNull();

    const ageLabel = page.locator(`${HOST} .d4-multi-form-header [name="div-AGE"]`).first();

    const seq: (string | null)[] = [];
    for (let i = 0; i < 3; i++) {
      await ageLabel.dblclick();

      await page.waitForTimeout(1500);
      seq.push(await sortArrow(page, 'AGE'));
    }
    expect(new Set(seq).size).toBe(3);
    expect(seq).toContain(null);
    expect(seq.filter((s) => s !== null).length).toBe(2);

    await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      vw.setOptions({sortByColumnName: null});
    });

    await v.pollValue(() => sortArrow(page, 'AGE'), (a) => a === null, 1500, 100);
    await ageLabel.dblclick();

    await page.waitForTimeout(1500);
    await expect.poll(() => sortArrow(page, 'AGE'), {timeout: 15_000}).toBe('↓');
    await expect.poll(() => sortIndicatorLabels(page), {timeout: 15_000}).toEqual(['div-AGE']);

    await page.locator(`${HOST} .d4-multi-form-header [name="div-HEIGHT"]`).first().dblclick();

    await page.waitForTimeout(1500);
    expect(await sortIndicatorLabels(page)).toEqual(['div-AGE']);
  });

  await softStep('Step 6a — Pin Row moves the card to the pinned pane and removes it from the ordinary set', async () => {

    await v.ensurePropertyCategory(page, 'Forms', 'misc', 'show-mouse-over-row');
    await v.setPropertyGridCheckbox(page, 'show-mouse-over-row', false, 'misc');
    await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      vw.setOptions({sortByColumnName: null});
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
    await expectNoBalloonSustained(page);

    await expect.poll(async () =>
      page.evaluate((sel) => getComputedStyle(document.querySelector(sel) as HTMLElement).display, PINNED_PANE),
    {timeout: 15_000}).not.toBe('none');
    const pinnedValues = await page.evaluate((sel) => Array.from(document.querySelectorAll(sel))
      .map((c) => ((c as HTMLElement).querySelector('[column="USUBJID"]') as HTMLInputElement)?.value), PINNED);
    expect(pinnedValues).toEqual([anchor]);

    await expect.poll(() => page.locator(ORDINARY).count(), {timeout: 15_000}).toBe(beforeCount - 1);
    expect(await ordinaryUsubjids(page)).not.toContain(anchor);
    expect(await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      return vw.props.pinnedRowValues;
    })).toEqual([anchor]);

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
    await expectNoBalloonSustained(page);
    await expect.poll(() => page.evaluate((sel) => Array.from(document.querySelectorAll(sel)).length, PINNED),
      {timeout: 15_000}).toBe(1);
  });

  await softStep('Step 6c — Pinning through a NON-UNIQUE field raises the exact warning; the single pin is preserved', async () => {

    const pinnedBefore = await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      return Array.from(vw.props.pinnedRowValues as string[]);
    });
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

    expect(await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      return Array.from(vw.props.pinnedRowValues as string[]);
    })).toContain(target.sex);

    const pinnedSexIdx = await cardIndexByValue(page, 'USUBJID', target.usub, PINNED);
    expect(pinnedSexIdx).toBeGreaterThanOrEqual(0);
    await cardContextMenu(page, PINNED, pinnedSexIdx, 'div-Unpin-Row');
    await expect.poll(() => page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      return Array.from(vw.props.pinnedRowValues as string[]);
    }), {timeout: 15_000}).toEqual(pinnedBefore);
  });

  await softStep('Step 7a / Step 7b — Re-applying the saved layout over a deliberately corrupted view restores the field set, sort-label identity and pinned row by value, and drops a foreign viewer not in the layout', async () => {

    const chosenFields = await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      const all = Array.from(vw.props.fieldsColumnNames as string[]);
      const must = ['USUBJID', 'AGE'];
      const rest = all.filter((n) => !must.includes(n));
      const picked = [...must, ...rest.slice(0, Math.max(1, rest.length - 2))].reverse();
      vw.setOptions({fieldsColumnNames: picked});
      return picked;
    });
    await expect.poll(() => drawnLabelNames(page), {timeout: 20_000}).toEqual(chosenFields);

    await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      vw.setOptions({sortByColumnName: 'AGE'});
    });
    await expect.poll(() => sortIndicatorLabels(page), {timeout: 20_000}).toEqual(['div-AGE']);

    const preLabels = await drawnLabelNames(page);
    const pre = await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      return {
        fields: Array.from(vw.props.fieldsColumnNames as string[]),
        pinnedValues: Array.from(vw.props.pinnedRowValues as string[]),
      };
    });
    const preIndicator = await sortIndicatorLabels(page);

    expect(preLabels.length).toBeGreaterThan(0);
    expect(pre.pinnedValues.length).toBeGreaterThan(0);
    expect(preIndicator.length).toBeGreaterThan(0);

    const currentUserId = await page.evaluate(() => String(grok.shell.user.id));
    const applicableLayouts = async (): Promise<{id: string; authorId: string | null}[]> =>
      page.evaluate(async () => ((await grok.dapi.layouts.getApplicable(grok.shell.t)) ?? [])
        .map((l: any) => ({id: String(l.id), authorId: l.author && l.author.id ? String(l.author.id) : null})));
    const beforeIds = (await applicableLayouts()).map((l) => l.id);

    expect(await v.driveTopMenuLeaf(page, ['View', 'Layout', 'Save to Gallery'])).toBe(true);

    let fresh: string[] = [];
    await expect.poll(async () => {

      fresh = (await applicableLayouts())
        .filter((l) => !beforeIds.includes(l.id) && l.authorId === currentUserId).map((l) => l.id);
      return fresh.length;
    }, {timeout: 20_000, intervals: [500, 1000, 2000, 3000]}).toBeGreaterThanOrEqual(1);
    if (fresh.length !== 1) {
      throw new Error(
        `Layout save produced ${fresh.length} new applicable layouts authored by the current user ` +
        `(${fresh.join(', ')}), expected exactly 1 — refusing to guess which is ours; deleting none. ` +
        `A concurrent layout creation on the same dataset is the likely cause.`);
    }
    const layoutId = fresh[0];

    try {

      await page.evaluate(() => {
        const tv = grok.shell.tv;
        const forms = tv.viewers.find((x: any) => x.type === 'FormsViewer');
        if (forms) forms.close();
        tv.addViewer('Histogram');
      });

      await expect.poll(() => page.locator('[name="viewer-Histogram"]').count(), {timeout: 20_000}).toBe(1);
      await expect.poll(() => page.locator(HOST).count(), {timeout: 20_000}).toBe(0);

      const applyErrTexts: string[] = [];
      const applyErrCount = await withConsoleErrorCount(page, async () => {
        await page.evaluate(async (id) => {
          const saved = await grok.dapi.layouts.find(id);
          grok.shell.tv.loadLayout(saved);
        }, layoutId);
        await page.locator(HOST).first().waitFor({timeout: 30_000});
      }, undefined, applyErrTexts);
      expect(applyErrCount, `layout-apply console errors: ${JSON.stringify(applyErrTexts)}`).toBe(0);

      await expect.poll(() => page.locator('[name="viewer-Histogram"]').count(), {timeout: 20_000}).toBe(0);
      expect(await page.evaluate(() => grok.shell.tv.viewers
        .filter((x: any) => x.type === 'Histogram').length)).toBe(0);
      await expect.poll(() => drawnLabelNames(page), {timeout: 20_000}).toEqual(preLabels);
      await expect.poll(() => page.evaluate(() => {
        const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
        return vw ? Array.from(vw.props.fieldsColumnNames as string[]) : null;
      }), {timeout: 20_000}).toEqual(pre.fields);

      await expect.poll(() => sortIndicatorLabels(page), {timeout: 20_000}).toEqual(preIndicator);

      await expect.poll(() => page.evaluate((sel) => Array.from(document.querySelectorAll(sel))
        .map((c) => ((c as HTMLElement).querySelector('[column="USUBJID"]') as HTMLInputElement)?.value), PINNED),
      {timeout: 20_000}).toEqual(pre.pinnedValues);
    } finally {
      await page.evaluate(async (id) => {
        const saved = await grok.dapi.layouts.find(id);
        if (saved) await grok.dapi.layouts.delete(saved);
      }, layoutId);
    }
  });

  await softStep('Step 7c — A project round-trip preserves the field set and the pinned row across a session boundary', async () => {
    const preLabels = await drawnLabelNames(page);
    const pre = await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');

      return {
        fields: Array.from(vw.props.fieldsColumnNames as string[]),
        pinnedValues: Array.from(vw.props.pinnedRowValues as string[]),
      };
    });
    const preIndicator = await sortIndicatorLabels(page);

    expect(preLabels.length).toBeGreaterThan(0);
    expect(pre.fields.length).toBeGreaterThan(0);
    expect(pre.pinnedValues.length).toBeGreaterThan(0);
    expect(preIndicator.length).toBeGreaterThan(0);

    const projName = `zz-formsviewer-core-${Date.now()}`;
    let projectId: string | null = null;
    try {

      const saved = await projects.saveProjectViaUI(page, projName);
      projectId = saved.projectId;

      await projects.reopenAndAssertProvenance(page, projectId as string);
      await page.locator(HOST).first().waitFor({timeout: 30_000});

      await expect.poll(() => drawnLabelNames(page), {timeout: 30_000}).toEqual(preLabels);
      await expect.poll(() => page.evaluate(() => {
        const vw = grok.shell.tv?.viewers?.find((x: any) => x.type === 'FormsViewer');
        return vw ? Array.from(vw.props.fieldsColumnNames as string[]) : null;
      }), {timeout: 30_000}).toEqual(pre.fields);

      await expect.poll(() => page.evaluate((sel) => Array.from(document.querySelectorAll(sel))
        .map((c) => ((c as HTMLElement).querySelector('[column="USUBJID"]') as HTMLInputElement)?.value), PINNED),
      {timeout: 30_000}).toEqual(pre.pinnedValues);

      await expect.poll(() => sortIndicatorLabels(page), {timeout: 30_000}).toEqual(preIndicator);
    } finally {
      if (projectId)
        await projects.deleteProjectWithCleanup(page, {projectId});
    }
  });

  await v.cleanupShell(page);
  v.finishSpec();
});
