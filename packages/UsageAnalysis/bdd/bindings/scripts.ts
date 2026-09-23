/* What only the Scripts features need: the Results table the script view puts under the code, the
   Save button's state while the core names it with a class alone, the editor's Save (DiffStudio's
   bindings own that phrase too, with state of their own), and the layouts the view leaves behind.
   Scripts on the server, the console, the alerts and the pane counts are the library's. */
import type {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, deleteChatsOf, expect, pollMs} from '@datagrok-libraries/bdd/runtime';
import {deleteLayoutsAtEnd} from '../helpers/layouts.js';

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
    timeout: pollMs(120000)}).toBe(value);
}, {description: 'the Results table under the editor after a run: the value column of the output\'s row'});

export const scriptResultListed = Then('the script results should list {string}', async (page: Page, output: string) => {
  await expect.poll(() => page.evaluate((o) => Array.from(document.querySelectorAll('.d4-item-table'))
    .filter((t) => (t as HTMLElement).offsetParent !== null && /name\s+value/.test((t as HTMLElement).innerText))
    .some((t) => Array.from(t.querySelectorAll('td')).some((c) => c.textContent?.trim() === o)), output),
  {message: `a row for "${output}" in the script results`, timeout: pollMs(180000)}).toBe(true);
}, {description: 'the run has ended with that output in the Results table under the editor, whatever its value'});

/* A layout saved from the script view is named after the script's dataframe output ("Df", "Df_1"). */
export const cleanLayouts = Given('the layouts saved for the script are deleted at the end', async (page: Page) => {
  deleteLayoutsAtEnd(page, '^Df(_\\d+)?$');
}, {tier: 'api', description: 'every "Df"-named layout of this account made since the step ran, and the project it belongs to, go when the feature ends'});

/* The editor's Save: the ribbon button, done when it reads "Saved"; the script it creates is
   deleted with its chats at feature end. Not in the library: DiffStudio's bindings own the same
   phrase with state of their own (the saved script its Model Hub claims read), and unifying the
   two is a change to that suite. */
export const saveScript = When('user saves the script', async (page: Page) => {
  const save = page.locator('[name="button-Save"]').filter({visible: true}).first();
  await save.click();
  await expect(save, 'the Save button after the save').toHaveText('Saved', {timeout: pollMs(60000)});
  let id = '';
  await expect.poll(async () => {
    id = await page.evaluate(async () => {
      const found = String(grok.shell.v?.path ?? '').match(/^\/script\/([^/?#]+)/)?.[1] ?? '';
      return found && (await grok.dapi.scripts.find(found).catch(() => null)) ? found : '';
    });
    return id !== '';
  }, {message: "the script view's script on the server", timeout: pollMs(30000)}).toBe(true);
  atFeatureEnd(page, async () => {
    await deleteChatsOf(page, id);
    await page.evaluate(async (scriptId) => {
      const script = await grok.dapi.scripts.find(scriptId).catch(() => null);
      if (script)
        await grok.dapi.scripts.delete(script);
    }, id);
  });
}, {tier: 'ui', description: 'the ribbon Save of the script view, done when it reads "Saved"; the script is deleted with its chats at feature end'});

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
