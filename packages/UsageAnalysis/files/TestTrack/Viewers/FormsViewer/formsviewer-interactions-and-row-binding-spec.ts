/* ---
realizes: [formsviewer.cp.interactions-and-row-binding]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {HOST, ORDINARY, CURRENT, cardFieldValue, fieldValuesByPosition, waitForOrderStable} from '../../helpers/forms';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

async function rowOfValue(page: Page, column: string, value: string): Promise<number> {
  return page.evaluate(({col, val}) => {
    const df = grok.shell.t;
    for (let i = 0; i < df.rowCount; i++) if (df.col(col).get(i) === val) return i;
    return -1;
  }, {col: column, val: value});
}

async function establishSetup(page: Page): Promise<void> {
  await page.evaluate(() => {
    const df = grok.shell.t;
    const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
    vw.props.showSelectedRows = true;
    vw.props.showCurrentRow = true;
    vw.props.showMouseOverRow = true;
    df.mouseOverRowIdx = -1;
    df.currentRowIdx = 5;
    df.selection.setAll(false);
    [2, 5, 10, 15].forEach((r) => df.selection.set(r, true));
  });

  await expect.poll(() => page.evaluate((sel) => {
    const df = grok.shell.t;
    const usub = (r: number) => df.col('USUBJID').get(r);
    const selected: number[] = [];
    for (let i = 0; i < df.rowCount; i++) if (df.selection.get(i)) selected.push(i);
    const expectedStable = [usub(df.currentRowIdx), ...selected.map(usub)];
    const actual = Array.from(document.querySelectorAll(sel)).map((c) =>
      ((c as HTMLElement).querySelector('[column="USUBJID"]') as HTMLInputElement | null)?.value ?? null);
    const actualStable = actual.length >= 2 ? [actual[0], ...actual.slice(2)] : actual;
    return JSON.stringify(actualStable) === JSON.stringify(expectedStable);
  }, ORDINARY), {timeout: 15_000}).toBe(true);
}

async function parkPointerBelowViewer(page: Page): Promise<void> {
  const hostBox = (await page.locator(HOST).boundingBox())!;
  const viewportH = page.viewportSize()?.height ?? Math.round(hostBox.y + hostBox.height + 40);
  const x = Math.round(hostBox.x + hostBox.width / 2);
  const y = Math.round(Math.min(hostBox.y + hostBox.height + 10, viewportH - 6));
  await page.mouse.move(x, y, {steps: 10});
}

test('Forms viewer — mouse interactions and row binding (p1)', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await softStep('Setup — add the Forms viewer, confirm the three show-toggles are ON, seed the entry state', async () => {
    await v.addViewerByIcon(page, 'Forms', 'Forms', 30_000, 'FormsViewer');
    await page.locator('.d4-multi-form').first().waitFor({timeout: 30_000});

    const defaults = await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      return {
        selected: vw.props.showSelectedRows,
        current: vw.props.showCurrentRow,
        mouseOver: vw.props.showMouseOverRow,
      };
    });
    expect(defaults).toEqual({selected: true, current: true, mouseOver: true});

    await establishSetup(page);
  });

  await softStep('Scenario 1 — Clicking a non-current card sets df.currentRowIdx and the field matches df.get', async () => {
    await establishSetup(page);

    const targetUsub = await cardFieldValue(page, 2, 'USUBJID');
    expect(targetUsub).not.toBeNull();
    const targetRow = await rowOfValue(page, 'USUBJID', targetUsub as string);
    const clickedSex = await cardFieldValue(page, 2, 'SEX');

    await page.locator(ORDINARY).nth(2).click();

    await expect.poll(() => page.evaluate(() => grok.shell.t.currentRowIdx), {timeout: 10_000}).toBe(targetRow);
    const sexOfCurrent = await page.evaluate(() =>
      grok.shell.t.get('SEX', grok.shell.t.currentRowIdx));
    expect(clickedSex).toBe(sexOfCurrent);
  });

  await softStep('Scenario 2 — Clicking a card field sets df.currentCell to that column and the card row', async () => {
    await establishSetup(page);
    const cardUsub = await cardFieldValue(page, 3, 'USUBJID');
    expect(cardUsub).not.toBeNull();
    const cardRow = await rowOfValue(page, 'USUBJID', cardUsub as string);

    await page.locator(ORDINARY).nth(3).locator('[column="AGE"]').click();

    await expect.poll(() => page.evaluate(() => {
      const cc = grok.shell.t.currentCell;
      return JSON.stringify({col: cc?.column?.name ?? null, row: cc?.rowIndex ?? null});
    }), {timeout: 10_000}).toBe(JSON.stringify({col: 'AGE', row: cardRow}));
  });

  await softStep('Scenario 3 — Clicking the HEIGHT header label sets df.currentCol to HEIGHT', async () => {
    await establishSetup(page);
    await page.locator(`${HOST} .d4-multi-form-header [name="div-HEIGHT"]`).first().click();
    await expect.poll(() => page.evaluate(() => grok.shell.t.currentCol?.name ?? null),
      {timeout: 10_000}).toBe('HEIGHT');
  });

  await softStep('Scenario 4a — Ctrl+clicking a card toggles df.selection.get(row) off then back on', async () => {

    const row = 20;
    await page.evaluate((r) => {
      const df = grok.shell.t;
      df.mouseOverRowIdx = -1;
      df.currentRowIdx = r;
      df.selection.setAll(false);
      df.selection.set(r, true);
    }, row);
    await expect.poll(() => page.evaluate((r) => grok.shell.t.selection.get(r), row),
      {timeout: 10_000}).toBe(true);
    await expect.poll(() => cardFieldValue(page, 0, 'USUBJID', CURRENT),
      {timeout: 10_000}).toBe(await page.evaluate((r) => grok.shell.t.col('USUBJID').get(r), row));

    await waitForOrderStable(page);
    await page.locator(CURRENT).first().click({modifiers: ['Control']});
    await expect.poll(() => page.evaluate((r) => grok.shell.t.selection.get(r), row),
      {timeout: 10_000}).toBe(false);

    await waitForOrderStable(page);
    await page.locator(CURRENT).first().click({modifiers: ['Control']});
    await expect.poll(() => page.evaluate((r) => grok.shell.t.selection.get(r), row),
      {timeout: 10_000}).toBe(true);
  });

  await softStep('Scenario 4b — Shift+clicking the row-12 card selects exactly rows 0..12 inclusive', async () => {
    const row = 12;

    await page.evaluate((r) => {
      const df = grok.shell.t;
      df.mouseOverRowIdx = -1;
      df.currentRowIdx = r;
      df.selection.setAll(false);
      [3, 7, 12, 18, 22].forEach((i) => df.selection.set(i, true));
    }, row);
    await expect.poll(() => cardFieldValue(page, 0, 'USUBJID', CURRENT),
      {timeout: 10_000}).toBe(await page.evaluate((r) => grok.shell.t.col('USUBJID').get(r), row));
    const tailSeeded = await page.evaluate(() => [18, 22].every((i) => grok.shell.t.selection.get(i)));
    expect(tailSeeded).toBe(true);

    await waitForOrderStable(page);

    const box = (await page.locator(CURRENT).first().boundingBox())!;
    await page.keyboard.down('Shift');
    try {
      await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
    } finally {
      await page.keyboard.up('Shift');
    }

    await expect.poll(() => page.evaluate((r) => {
      const df = grok.shell.t;
      const sel: number[] = [];
      for (let i = 0; i < df.rowCount; i++) if (df.selection.get(i)) sel.push(i);
      const expected: number[] = [];
      for (let i = 0; i <= r; i++) expected.push(i);
      return JSON.stringify(sel) === JSON.stringify(expected);
    }, row), {timeout: 40_000}).toBe(true);
  });

  await softStep('Scenario 4c — Ctrl+Shift+clicking the row-8 card clears rows 0..8; rows above stay selected', async () => {

    const row = 8;
    await page.evaluate((r) => {
      const df = grok.shell.t;
      df.mouseOverRowIdx = -1;
      df.currentRowIdx = r;
    }, row);
    await expect.poll(() => cardFieldValue(page, 0, 'USUBJID', CURRENT),
      {timeout: 10_000}).toBe(await page.evaluate((r) => grok.shell.t.col('USUBJID').get(r), row));
    const nineToTwelveBefore = await page.evaluate(() =>
      [9, 10, 11, 12].every((r) => grok.shell.t.selection.get(r)));
    expect(nineToTwelveBefore).toBe(true);

    await waitForOrderStable(page);
    await page.locator(CURRENT).first().click({modifiers: ['Control', 'Shift']});

    await expect.poll(() => page.evaluate((r) => {
      const df = grok.shell.t;
      const cleared = Array.from({length: r + 1}, (_, i) => i).every((i) => !df.selection.get(i));
      const keptAbove = [9, 10, 11, 12].every((i) => df.selection.get(i));
      return JSON.stringify({cleared, keptAbove});
    }, row), {timeout: 10_000}).toBe(JSON.stringify({cleared: true, keptAbove: true}));
  });

  await softStep('Scenario 5a — Hovering an ordinary card sets df.mouseOverRowIdx to that card row', async () => {
    await establishSetup(page);

    const targetUsub = await cardFieldValue(page, 2, 'USUBJID');
    const targetRow = await rowOfValue(page, 'USUBJID', targetUsub as string);

    const beforeHover = await page.evaluate(() => grok.shell.t.mouseOverRowIdx);
    expect(beforeHover).not.toBe(targetRow);

    await page.locator(ORDINARY).nth(2).hover();
    await expect.poll(() => page.evaluate(() => grok.shell.t.mouseOverRowIdx),
      {timeout: 10_000}).toBe(targetRow);

    await parkPointerBelowViewer(page);
    await expect.poll(() => page.evaluate(() => grok.shell.t.mouseOverRowIdx),
      {timeout: 10_000}).toBe(-1);

    await expect.poll(() => cardFieldValue(page, 1, 'USUBJID'), {timeout: 10_000}).toBeNull();

    const currentUsub5a = await page.evaluate(() => grok.shell.t.col('USUBJID').get(grok.shell.t.currentRowIdx));
    await expect.poll(() => cardFieldValue(page, 0, 'USUBJID'), {timeout: 10_000}).toBe(currentUsub5a);
    expect(await page.locator(CURRENT).count()).toBe(1);
  });

  await softStep('Scenario 5b — The mouse-over card USUBJID field goes empty → hovered row value → empty', async () => {
    await establishSetup(page);
    await parkPointerBelowViewer(page);
    await page.evaluate(() => { grok.shell.t.mouseOverRowIdx = -1; });

    await expect.poll(() => cardFieldValue(page, 1, 'USUBJID'), {timeout: 10_000}).toBeNull();

    const hoverRow = 33;
    await page.evaluate((r) => { grok.shell.t.mouseOverRowIdx = r; }, hoverRow);
    const expected = await page.evaluate((r) => grok.shell.t.col('USUBJID').get(r), hoverRow);
    await expect.poll(() => cardFieldValue(page, 1, 'USUBJID'), {timeout: 10_000}).toBe(expected);

    await page.evaluate(() => { grok.shell.t.mouseOverRowIdx = -1; });
    await expect.poll(() => cardFieldValue(page, 1, 'USUBJID'), {timeout: 10_000}).toBeNull();
  });

  let baseCount = 0;
  await softStep('Scenario 6a — Turning Show Mouse Over Row OFF drops one card and stops the grid-hover fill', async () => {
    await establishSetup(page);

    baseCount = await page.locator(ORDINARY).count();

    await v.ensurePropertyCategory(page, 'Forms', 'misc', 'show-mouse-over-row');
    await v.setPropertyGridCheckbox(page, 'show-mouse-over-row', false, 'misc');
    await expect.poll(() => page.locator(ORDINARY).count(), {timeout: 15_000}).toBe(baseCount - 1);
    await expect.poll(() => page.locator(`${ORDINARY} .d4-multi-form-form-indicator-mouse-over-row`).count(),
      {timeout: 15_000}).toBe(0);

    const currentUsub = await page.evaluate(() =>
      grok.shell.t.col('USUBJID').get(grok.shell.t.currentRowIdx));
    await page.evaluate(() => { grok.shell.t.mouseOverRowIdx = 40; });

    await page.waitForTimeout(1500);
    const hovered40 = await page.evaluate(() => grok.shell.t.col('USUBJID').get(40));
    const shown = await fieldValuesByPosition(page, 'USUBJID');
    expect(shown).not.toContain(hovered40);
    expect(await cardFieldValue(page, 0, 'USUBJID', CURRENT)).toBe(currentUsub);
  });

  await softStep('Scenario 6a restore — Turning Show Mouse Over Row back ON returns the card count to the baseline', async () => {

    await page.evaluate(() => { grok.shell.t.mouseOverRowIdx = -1; });
    await v.ensurePropertyCategory(page, 'Forms', 'misc', 'show-mouse-over-row');
    await v.setPropertyGridCheckbox(page, 'show-mouse-over-row', true, 'misc');
    await expect.poll(() => page.locator(ORDINARY).count(), {timeout: 15_000}).toBe(baseCount);
  });

  await softStep('Scenario 6b — Turning Show Current Row OFF removes the current-row card and its green indicator', async () => {
    await establishSetup(page);
    await expect.poll(() => page.locator(CURRENT).count(),
      {timeout: 15_000}).toBe(1);

    await v.ensurePropertyCategory(page, 'Forms', 'misc', 'show-current-row');
    await v.setPropertyGridCheckbox(page, 'show-current-row', false, 'misc');

    await expect.poll(() => page.locator(CURRENT).count(),
      {timeout: 15_000}).toBe(0);

    await v.ensurePropertyCategory(page, 'Forms', 'misc', 'show-current-row');
    await v.setPropertyGridCheckbox(page, 'show-current-row', true, 'misc');
  });

  await v.cleanupShell(page);
  v.finishSpec();
});
