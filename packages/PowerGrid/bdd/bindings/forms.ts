/* The steps only the Forms viewer's features need. Everything else they use — cards, fields,
   labels and sort indicators as hit areas, the readings behind "cards" / "records shown" /
   "pinned records" / "fields" / "sort column", the properties, the selection and the filter — is
   the library's viewers tier and the platform's data steps (`grok-bdd list-steps`).

   `user renames {string} column to {string}` and `the current column should be {string}` are
   generic: promote them to the library (`bindings/platform/columns.ts`) as soon as a second
   package wants them. */
import {expect, Page} from '@playwright/test';
import {Then, When} from '@datagrok-libraries/bdd';
import {el, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/** The pinned pane the viewer reports as its `pinned` part: `display: none` until a row is pinned. */
async function pinnedPaneDisplay(page: Page): Promise<string> {
  return page.evaluate(() => {
    const found = Array.from(grok.shell.tv?.viewers ?? []).find((v: any) => v.type === 'FormsViewer' || v.type === 'Forms') as any;
    if (!found)
      throw new Error('no Forms viewer in the current table view');
    const pane = found.getWidgetStatus()?.parts?.pinned as HTMLElement | undefined;
    if (!pane)
      throw new Error('the Forms viewer reports no "pinned" part');
    return getComputedStyle(pane).display;
  });
}

export const pinnedPaneHidden = Then('the pinned pane of forms viewer should be hidden', async (page: Page) => {
  await expect.poll(() => pinnedPaneDisplay(page), {message: 'the pinned pane of forms viewer'}).toBe('none');
}, {description: 'the pane the pinned cards live in, gone once nothing is pinned'});

export const pinnedPaneShown = Then('the pinned pane of forms viewer should be shown', async (page: Page) => {
  await expect.poll(() => pinnedPaneDisplay(page), {message: 'the pinned pane of forms viewer'}).not.toBe('none');
});

export const renameColumn = When('user renames {string} column to {string}', async (page: Page, from: string, to: string) => {
  await viewers.baselineAll(page);
  await page.evaluate(([f, t]) => {
    const col = grok.shell.t.col(f);
    if (!col)
      throw new Error(`no "${f}" column in ${grok.shell.t.name}; it has: ${grok.shell.t.columns.names().join(', ')}`);
    col.name = t;
  }, [from, to]);
  await viewers.settleAll(page);
}, {tier: 'api', description: 'the UI path is Column Properties... from the grid header; a name starting with "~" hides the column'});

/** The library's `holding {key}` passes one key to Playwright, and a chord ("Control+Shift") is not
 * one — promote a two-modifier form to the library when a second viewer needs it. */
export const chordClick = When('user clicks on the {string} area of forms viewer holding {word} and {word}',
  async (page: Page, area: string, first: string, second: string) => {
    const target = el('forms viewer');
    const point = viewers.centerOf(await viewers.hitArea(page, target, area, true));
    await page.keyboard.down(first);
    await page.keyboard.down(second);
    await page.mouse.click(point.x, point.y);
    await page.keyboard.up(second);
    await page.keyboard.up(first);
  }, {tier: 'ui', description: 'a click with two modifiers held — Control+Shift on a card deselects every row up to it'});

export const currentColumnIs = Then('the current column should be {string}', async (page: Page, name: string) => {
  await expect.poll(() => page.evaluate(() => grok.shell.t.currentCol?.name ?? ''),
    {message: 'the current column of the table'}).toBe(name);
}, {description: 'dataFrame.currentCol — what a click on a field or on a header label sets'});
