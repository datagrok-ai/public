import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const VIEWER_NAME = 'Form';
const VIEWER = `[name="viewer-${VIEWER_NAME}"]`;
const datasetPath = 'System:DemoFiles/demog.csv';

// String columns are compared verbatim; numeric and date fields are formatted
// for display, so those are only checked for "it changed with the row".
const TEXT_COLUMNS = ['USUBJID', 'SEX', 'RACE', 'DIS_POP'];

/** Value the form shows for a column — the field is an input or a plain div. */
async function fieldValue(page: Page, column: string): Promise<string> {
  return page.evaluate(({sel, col}) => {
    const root = document.querySelector(sel) as HTMLElement;
    const el = root.querySelector(`[name="input-${col.replace(/_/g, '-')}"]`) as HTMLElement | null;
    if (!el) return '';
    const input = (el as HTMLInputElement).value !== undefined ? (el as HTMLInputElement).value : null;
    const inner = el.querySelector('input') as HTMLInputElement | null;
    return (input || inner?.value || el.innerText || '').trim();
  }, {sel: VIEWER, col: column});
}

/** The same values straight from the table, for the row the form is showing. */
async function rowValues(page: Page, columns: string[]): Promise<Record<string, string>> {
  return page.evaluate((cols) => {
    const df = (window as any).grok.shell.t;
    const i = df.currentRowIdx;
    const out: Record<string, string> = {};
    for (const c of cols) out[c] = String(df.col(c).get(i));
    return out;
  }, columns);
}

const currentRow = (page: Page) => page.evaluate(() => (window as any).grok.shell.t.currentRowIdx);

test('Form viewer', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  // #### Add the viewer
  await softStep('Add Form from the Viewers toolbox', async () => {
    await page.locator('[name="icon-form"]').first().click();
    await page.locator(VIEWER).first().waitFor({timeout: 30_000});

    // Every column of the table gets a field.
    const columns = await page.evaluate(() =>
      (window as any).grok.shell.t.columns.names() as string[]);
    for (const column of columns)
      await expect(page.locator(`${VIEWER} [name="div-${column.replace(/_/g, '-')}"]`).first())
        .toBeVisible();
  });

  // #### The form shows the current row
  await softStep('The fields show the current row', async () => {
    const expected = await rowValues(page, TEXT_COLUMNS);
    for (const column of TEXT_COLUMNS)
      expect(await fieldValue(page, column)).toBe(expected[column]);
  });

  // #### Navigation arrows
  await softStep('The next and previous arrows walk the rows', async () => {
    const startRow = await currentRow(page);
    const startValues = await rowValues(page, TEXT_COLUMNS);

    await page.locator(`${VIEWER} [name="icon-chevron-right"]`).first().click();
    await expect.poll(() => currentRow(page), {timeout: 8000}).toBe(startRow + 1);

    // The form must follow the row, not keep showing the previous one.
    const nextValues = await rowValues(page, TEXT_COLUMNS);
    for (const column of TEXT_COLUMNS)
      expect(await fieldValue(page, column)).toBe(nextValues[column]);
    expect(nextValues['USUBJID']).not.toBe(startValues['USUBJID']);

    await page.locator(`${VIEWER} [name="icon-chevron-left"]`).first().click();
    await expect.poll(() => currentRow(page), {timeout: 8000}).toBe(startRow);
    expect(await fieldValue(page, 'USUBJID')).toBe(startValues['USUBJID']);
  });

  // #### The form follows the row picked in the grid
  await softStep('Picking a row in the grid updates the form', async () => {
    await page.evaluate(() => { (window as any).grok.shell.t.currentRowIdx = 42; });
    await expect.poll(() => currentRow(page), {timeout: 8000}).toBe(42);

    const expected = await rowValues(page, TEXT_COLUMNS);
    for (const column of TEXT_COLUMNS)
      expect(await fieldValue(page, column)).toBe(expected[column]);
  });

  // #### Row selector
  await softStep('The row selector selects the row on show', async () => {
    const selected = () =>
      page.evaluate(() => (window as any).grok.shell.t.selection.trueCount as number);
    await page.evaluate(() => (window as any).grok.shell.t.selection.setAll(false));
    await expect.poll(selected, {timeout: 5000}).toBe(0);

    await page.locator(`${VIEWER} [name="icon-square"]`).first().click();
    await expect.poll(selected, {timeout: 8000}).toBe(1);

    // The row that got selected must be the one the form is showing.
    const row = await currentRow(page);
    expect(await page.evaluate((i) => (window as any).grok.shell.t.selection.get(i), row)).toBe(true);

    await page.locator(`${VIEWER} [name="icon-square"]`).first().click();
    await expect.poll(selected, {timeout: 8000}).toBe(0);
  });

  // #### The navigation bar can be hidden
  await softStep('Show Navigation hides the toolbar', async () => {
    await v.openViewerProperties(page, VIEWER_NAME);
    await v.ensurePropertyCategory(page, VIEWER_NAME, 'misc', 'show-navigation');
    expect(await v.togglePropertyGridCheckbox(page, 'show-navigation')).toBe(false);
    await expect(page.locator(`${VIEWER} [name="icon-chevron-right"]`).first()).toBeHidden();

    expect(await v.togglePropertyGridCheckbox(page, 'show-navigation')).toBe(true);
    await expect(page.locator(`${VIEWER} [name="icon-chevron-right"]`).first()).toBeVisible();
  });

  // #### Individual toolbar icons
  await softStep('Show Next Row Arrow hides just that arrow', async () => {
    await v.ensurePropertyCategory(page, VIEWER_NAME, 'misc', 'show-next-row-arrow');
    expect(await v.togglePropertyGridCheckbox(page, 'show-next-row-arrow')).toBe(false);
    await expect(page.locator(`${VIEWER} [name="icon-chevron-right"]`).first()).toBeHidden();
    await expect(page.locator(`${VIEWER} [name="icon-chevron-left"]`).first()).toBeVisible();

    expect(await v.togglePropertyGridCheckbox(page, 'show-next-row-arrow')).toBe(true);
    await expect(page.locator(`${VIEWER} [name="icon-chevron-right"]`).first()).toBeVisible();
  });

  // #### Closing the viewer
  await softStep('Close the viewer from its title bar', async () => {
    await v.clickViewerTitlebarIcon(page, VIEWER_NAME, 'Close');
    await expect(page.locator(VIEWER)).toHaveCount(0);
  });

  await v.cleanupShell(page);

  v.finishSpec();
});
