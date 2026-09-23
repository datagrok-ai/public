/* What only the Scripts features need: the Results table the script view puts under the code, the
   Save button's state while the core names it with a class alone, the editor's Save (DiffStudio's
   bindings own that phrase too, with state of their own), and the layouts a failed save leaves
   behind. Scripts on the server, the console, the alerts and the pane counts are the library's. */
import type {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, expect, pollMs} from '@datagrok-libraries/bdd/runtime';

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

/* A layout saved from the script view is named after the script's dataframe output ("Df", "Df_1"),
   and the save that fails on the share still leaves it — with the project it made. */
export const cleanLayouts = Given('the layouts saved for the script are deleted at the end', async (page: Page) => {
  const since = Date.now() - 60 * 1000;
  atFeatureEnd(page, async () => {
    const left: string[] = await page.evaluate(async (from) => {
      const me = (await grok.dapi.users.current()).id;
      const gone: string[] = [];
      for (const layout of await grok.dapi.layouts.list({pageSize: 1000})) {
        const made = layout.createdOn ? new Date(layout.createdOn.toString()).getTime() : 0;
        if (made < from || String(layout.author?.id ?? '') !== String(me) || !/^Df(_\d+)?$/.test(String(layout.name)))
          continue;
        const project = (await grok.dapi.projects.list({pageSize: 1000})).find((p: any) => p.name === layout.name);
        if (project)
          await grok.dapi.projects.delete(project);
        await grok.dapi.layouts.delete(layout);
        gone.push(String(layout.name));
      }
      return gone;
    }, since);
    if (left.length > 0)
      console.warn(`bdd: deleted the layouts the script view left behind: ${left.join(', ')}`);
  });
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
    await page.evaluate(async (scriptId) => {
      const headers = {Authorization: String(grok.dapi.token)};
      const root = grok.dapi.root;
      const chats = await (await fetch(`${root}/chats?entityId=${scriptId}`, {headers})).json();
      for (const chat of Array.isArray(chats) ? chats : [])
        await fetch(`${root}/chats/${chat.id}`, {method: 'DELETE', headers});
      const script = await grok.dapi.scripts.find(scriptId).catch(() => null);
      if (script)
        await grok.dapi.scripts.delete(script);
    }, id);
  });
}, {tier: 'ui', description: 'the ribbon Save of the script view, done when it reads "Saved"; the script is deleted with its chats at feature end'});
