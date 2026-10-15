/* What only the Scripts features need: the Results table the script view puts under the code, the
   Save button's state while the core names it with a class alone, and the layouts the view leaves
   behind. Scripts on the server and their Save, the console, the alerts and the pane counts are the
   library's. */
import type {Page} from '@playwright/test';
import {Given, Then} from '@datagrok-libraries/bdd';
import {deleteLayoutsAtEnd, expect, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/** The name / value table the script view puts under the code after a run. */
export const scriptResult = Then('the script results should show {string} as {string}', async (page: Page, output: string, value: string) => {
  const read = () => page.evaluate((o) => {
    const table = Array.from(document.querySelectorAll('.d4-item-table')).find((t) => (t as HTMLElement).offsetParent !== null &&
      /name\s+value/.test((t as HTMLElement).innerText));
    if (!table)
      return null;
    const row = Array.from(table.querySelectorAll('tr')).map((r) => Array.from(r.querySelectorAll('td')).map((c) => c.textContent?.trim() ?? ''))
      .find((cells) => cells.includes(o));
    return row ? row[row.indexOf(o) + 1] ?? '' : '';
  }, output);
  await expect.poll(read, {message: `the value of "${output}" in the script results (null: no results table yet)`,
    timeout: pollMs(30000)}).toBe(value);
}, {description: 'the Results table under the editor after a run: the value column of the output\'s row'});

/* A layout saved from the script view is named after the script's dataframe output ("Df", "Df_1"). */
export const cleanLayouts = Given('the layouts saved for the script are deleted at the end', (page: Page) =>
  deleteLayoutsAtEnd(page, '^Df(_\\d+)?$'),
{tier: 'api', description: 'every "Df"-named layout of this account made since the step ran, and the project it belongs to, go when the feature ends'});

/* The Signature Editor keeps the ribbon it finds when it opens and puts that back when it is left
   (DevTools `function-signature-editor.ts`), and a view switched to a moment ago still has the
   previous render's icons in the DOM while its own panels are not on the view yet: opening the
   editor in that gap restores nothing and the ribbon stays empty. */
export const ribbonReady = Then('the ribbon of the current view should be ready', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => {
    const panels = grok.shell.v?.getRibbonPanels?.() ?? [];
    return panels.reduce((n: number, p: any[]) => n + p.length, 0);
  }), {message: "the icons the current view's own ribbon panels hold", timeout: pollMs(30000)}).toBeGreaterThan(0);
}, {tier: 'api', description: 'the view\'s own ribbon panels, not the icons the previous render left in the DOM'});
