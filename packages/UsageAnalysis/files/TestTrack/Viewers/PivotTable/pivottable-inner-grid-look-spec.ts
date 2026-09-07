/* ---
realizes: [pivottable.cp.inner-grid-look-viewers, pivottable.int.viewer-columns-per-pivot-category]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {
  PIVOT, INNER_CANVAS, openPivot, pivotProps, gridLookColumn, innerRect, headerPoint, rowY,
  rightClickHeader, clickGridLeaf, applyLinearColorCoding,
} from './pivot-helpers';

// The layout round-trip and the ribbon Save entry point (Scenario 6 Steps 9-10) live in
// pivot-table-server-spec.ts.
test.use(specTestOptions);

const KEY_X = 62;
const VALUE_X = 166;
const VALUE_COL = 'Critical avg(AGE)';

async function aggregateChipCount(page: Page): Promise<number> {
  return page.evaluate(() => {
    const scope = document.querySelector('[name="viewer-Pivot-table"]');
    const row = scope && [...scope.querySelectorAll('.grok-pivot-column-panel')]
      .find((p) => p.querySelector('.grok-pivot-column-tags-title[d4-name="Aggregate"]'));
    return row ? row.querySelectorAll('.d4-tag').length : 0;
  });
}

const VIEWER_PICKER = `${PIVOT} .grok-pivot-column-tags-title[d4-name="Aggregate"] .d4-combo-popup`;
async function pickViewerColumn(page: Page) {
  const box = await page.locator(VIEWER_PICKER).first().boundingBox();
  if (!box) throw new Error('viewer-picker combo-popup not visible');
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
  await page.locator(`${VIEWER_PICKER}.d4-combo-popup-expanded .d4-list-item`).first().waitFor({timeout: 5000});
  await page.locator(`${VIEWER_PICKER}.d4-combo-popup-expanded .d4-list-item [name="icon-scatter-plot"]`).first().click();
}

test('Pivot Table — Inner grid look and viewer columns', async ({page}) => {
  test.setTimeout(360_000);

  const consoleErrors: string[] = [];
  const ignorable = (m: string) => /Unable to find element in cloned iframe/i.test(m) || isLocalBootNoise(m);
  const onConsole = (m: any) => { if (m.type() === 'error' && !ignorable(m.text())) consoleErrors.push(m.text()); };
  page.on('console', onConsole);
  const errCount = () => consoleErrors.length;

  await openDatagrok(page);
  await openPivot(page);

  await softStep('Setup: tag-editor header shows Group by, Aggregate and Pivot rows with the auto cross-tab', async () => {
    const titles = await page.evaluate(() => Array.from(
      document.querySelectorAll('[name="viewer-Pivot-table"] .grok-pivot-column-tags-title'))
      .map((t) => t.getAttribute('d4-name')));
    expect(titles).toEqual(expect.arrayContaining(['Group by', 'Aggregate', 'Pivot']));
    const props = await pivotProps(page);
    expect(props.groupBy).toEqual(['DIS_POP']);
    expect(props.agg).toEqual(['AGE']);
    expect(props.aggTypes).toEqual(['avg']);
    expect(props.pivot).toEqual(['SEVERITY']);
    await page.locator(INNER_CANVAS).first().waitFor({timeout: 15000});
  });

  await softStep('Scenario 1 Step 4: the value column is colour-coded Linear after the header-menu action, no console error', async () => {
    const errBefore = errCount();

    const baseline = await gridLookColumn(page, VALUE_COL);
    expect(baseline?.colorCodingType).toBe('Off');
    await applyLinearColorCoding(page, VALUE_X);

    const coded = await v.pollValue(() => gridLookColumn(page, VALUE_COL),
      (c) => c?.colorCodingType === 'Linear', 3000, 100);
    expect(coded?.colorCodingType).toBe('Linear');
    expect(coded?.isColorCoded).toBe(true);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 3 Step 4: Grid > Hide hides the value column in the inner grid (gridLook visible flips false)', async () => {
    const beforeCol = await gridLookColumn(page, VALUE_COL);
    expect(beforeCol?.visible).toBe(true);

    const rc = await innerRect(page);
    const vpc = page.viewportSize() ?? {width: 1920, height: 1080};
    const cell = headerPoint(rc, vpc, VALUE_X);
    await page.mouse.click(cell.px, rc.y + rowY(0));
    await rightClickHeader(page, VALUE_X);
    await clickGridLeaf(page, 'div-Grid---Hide');
    await page.keyboard.press('Escape');

    const afterCol = await v.pollValue(() => gridLookColumn(page, VALUE_COL),
      (c) => c?.visible === false, 3000, 100);
    expect(afterCol?.visible).toBe(false);
  });

  await softStep('Scenario 3 Step 6: the select-all checkbox in Order or Hide Columns restores the hidden column', async () => {
    const errBefore = errCount();
    const dialog = page.locator('.d4-dialog[name="dialog-Order-or-Hide-Columns"]');
    await rightClickHeader(page, KEY_X);
    await clickGridLeaf(page, 'div-Grid---Order-or-Hide-Columns...');
    await dialog.waitFor({timeout: 8000});
    await dialog.locator('input[type="checkbox"]').first().click();
    const restored = await v.pollValue(() => gridLookColumn(page, VALUE_COL),
      (c) => c?.visible === true, 3000, 100);
    const close = dialog.locator('[name="button-CLOSE"]');
    if (await close.count() > 0) await close.first().click();
    else await page.keyboard.press('Escape');
    await dialog.waitFor({state: 'detached', timeout: 5000}).catch(() => {});
    expect(restored?.visible).toBe(true);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 4 Step 3: with one pivot column, the viewer picker adds a viewer column, no console error', async () => {
    const errBefore = errCount();

    await v.setViewerProps(page, 'Pivot table', [{set: {pivotColumnNames: ['SEVERITY']}}], 700);
    expect((await pivotProps(page)).pivot).toEqual(['SEVERITY']);

    const chipsBefore = await aggregateChipCount(page);
    await pickViewerColumn(page);

    const chipsAfter = await v.pollValue(() => aggregateChipCount(page),
      (n) => n === chipsBefore + 1, 3000, 100);
    expect(chipsAfter).toBe(chipsBefore + 1);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 4 Step 6: with two pivot columns configured, the viewer picker runs error-free', async () => {
    const errBefore = errCount();
    await v.setViewerProps(page, 'Pivot table', [{set: {pivotColumnNames: ['SEVERITY', 'SEX']}}], 900);
    expect((await pivotProps(page)).pivot).toEqual(['SEVERITY', 'SEX']);
    const chipsBefore = await aggregateChipCount(page);
    await pickViewerColumn(page);
    await v.pollValue(() => aggregateChipCount(page), (n) => n === chipsBefore + 1, 3000, 100);
    expect(errCount()).toBe(errBefore);
  });

  page.off('console', onConsole);
  await v.closeAllAndWait(page);
  v.finishSpec();
});
