import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

// The Forms viewer's root carries no `name` attribute, so it is addressed by its
// class instead of the usual `[name="viewer-…"]`.
const VIEWER = '.d4-multi-form';
const VIEWER_NAME = 'Forms';
const datasetPath = 'System:DemoFiles/demog.csv';

const TEXT_COLUMNS = ['USUBJID', 'SEX', 'RACE'];

/** How many row forms are on screen — one field name repeats once per form. */
const formCount = (page: Page) => page.locator(`${VIEWER} [name="input-SEX"]`).count();

/** Values of one form, by its index among the forms on screen. */
async function formValues(page: Page, index: number, columns: string[]): Promise<Record<string, string>> {
  return page.evaluate(({sel, i, cols}) => {
    const root = document.querySelector(sel) as HTMLElement;
    const out: Record<string, string> = {};
    for (const c of cols) {
      const all = root.querySelectorAll(`[name="input-${c.replace(/_/g, '-')}"]`);
      const el = all[i] as HTMLElement | undefined;
      const inner = el?.querySelector('input') as HTMLInputElement | null;
      out[c] = ((el as HTMLInputElement)?.value || inner?.value || el?.innerText || '').trim();
    }
    return out;
  }, {sel: VIEWER, i: index, cols: columns});
}

/** The same values straight from the table, for a given row. */
async function rowValues(page: Page, row: number, columns: string[]): Promise<Record<string, string>> {
  return page.evaluate(({r, cols}) => {
    const df = (window as any).grok.shell.t;
    const out: Record<string, string> = {};
    for (const c of cols) out[c] = String(df.col(c).get(r));
    return out;
  }, {r: row, cols: columns});
}

/** Real click on the viewer's title-bar Close button, found through its root. */
async function closeViewer(page: Page): Promise<void> {
  const point = await page.evaluate((sel) => {
    const root = document.querySelector(sel) as HTMLElement | null;
    const target = root?.closest('.panel-base')?.querySelector('[name="Close"]') as HTMLElement | null;
    if (!target) return null;
    const r = target.getBoundingClientRect();
    return {x: r.x + r.width / 2, y: r.y + r.height / 2};
  }, VIEWER);
  if (!point) throw new Error('no Close button on the Forms title bar');
  await page.mouse.click(point.x, point.y);
}

test('Forms viewer', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  // #### Add the viewer
  await softStep('Add Forms from the Viewers toolbox', async () => {
    await page.locator('[name="icon-Forms"]').first().click();
    await page.locator(VIEWER).first().waitFor({timeout: 30_000});

    const columns = await page.evaluate(() =>
      (window as any).grok.shell.t.columns.names() as string[]);
    for (const column of columns)
      await expect(page.locator(`${VIEWER} [name="div-${column.replace(/_/g, '-')}"]`).first())
        .toBeVisible();

    // Out of the box it shows the current row, and nothing else.
    await expect.poll(() => formCount(page), {timeout: 15_000}).toBe(1);
  });

  // #### The single form shows the current row
  await softStep('The form shows the current row and follows it', async () => {
    const row = await page.evaluate(() => (window as any).grok.shell.t.currentRowIdx);
    expect(await formValues(page, 0, TEXT_COLUMNS)).toEqual(await rowValues(page, row, TEXT_COLUMNS));

    await page.evaluate(() => { (window as any).grok.shell.t.currentRowIdx = 77; });
    await expect.poll(() => formValues(page, 0, TEXT_COLUMNS), {timeout: 10_000})
      .toEqual(await rowValues(page, 77, TEXT_COLUMNS));
  });

  // #### Selected rows get their own forms
  await softStep('Show Selected Rows adds a form per selected row', async () => {
    await v.openViewerProperties(page, VIEWER_NAME).catch(async () => {
      // The root has no name attribute, so the shared gear helper cannot find it.
      const point = await page.evaluate((sel) => {
        const root = document.querySelector(sel) as HTMLElement | null;
        const gear = root?.closest('.panel-base')
          ?.querySelector('.panel-titlebar [name="icon-font-icon-settings"]') as HTMLElement | null;
        if (!gear) return null;
        const r = gear.getBoundingClientRect();
        return {x: r.x + r.width / 2, y: r.y + r.height / 2};
      }, VIEWER);
      if (point) await page.mouse.click(point.x, point.y);
      await page.locator('.property-grid').first().waitFor({timeout: 10_000});
    });

    await v.ensurePropertyCategory(page, VIEWER_NAME, 'misc', 'show-selected-rows');
    // The setting ships enabled, so it is set rather than flipped.
    await v.setPropertyGridCheckbox(page, 'show-selected-rows', true);

    await page.evaluate(() => {
      const df = (window as any).grok.shell.t;
      df.selection.init((i: number) => i < 3);
    });

    // Three selected rows are shown, alongside the current-row form.
    await expect.poll(() => formCount(page), {timeout: 15_000}).toBeGreaterThanOrEqual(3);

    const shown = await formValues(page, 0, TEXT_COLUMNS);
    const candidates = await Promise.all([0, 1, 2, 77].map((r) => rowValues(page, r, TEXT_COLUMNS)));
    expect(candidates).toContainEqual(shown);

    await page.evaluate(() => (window as any).grok.shell.t.selection.setAll(false));
    await expect.poll(() => formCount(page), {timeout: 15_000}).toBe(1);
  });

  // #### The current-row form can be switched off
  await softStep('Show Current Row hides the current-row form', async () => {
    await v.ensurePropertyCategory(page, VIEWER_NAME, 'misc', 'show-current-row');
    await v.setPropertyGridCheckbox(page, 'show-current-row', false);
    await expect.poll(() => formCount(page), {timeout: 15_000}).toBe(0);

    await v.setPropertyGridCheckbox(page, 'show-current-row', true);
    await expect.poll(() => formCount(page), {timeout: 15_000}).toBe(1);
  });

  // #### Closing the viewer
  await softStep('Close the viewer from its title bar', async () => {
    await closeViewer(page);
    await expect(page.locator(VIEWER)).toHaveCount(0);
  });

  await v.cleanupShell(page);

  v.finishSpec();
});
