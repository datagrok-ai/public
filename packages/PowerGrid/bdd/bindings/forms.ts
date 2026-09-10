/* The one step only the Forms viewer's features need: whether its pinned pane is shown. Everything
   else they use — cards, fields, labels and sort indicators as hit areas, the readings behind
   "cards" / "records shown" / "pinned records" / "fields" / "sort column", the properties, the
   selection, the filter, renaming a column, the current column, a click on an area with a chord
   held — is the library's viewers tier and the platform's data steps (`grok-bdd list-steps`). */
import {expect, Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';

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
