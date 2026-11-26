/* Names and claims the Scripting and automation how-to clips (features/howto/13-scripting) need and the Tutorials
   vocabulary does not give. The console's command line is a bare <input class="ui-input-editor"> under
   `.d4-console-input`, outside any `.ui-input-root`, so no input kind reaches it, and the library's console steps
   press the backquote and fill it off camera: a clip clicks it and types into it on camera. The Results table the
   script view draws under the code after a run is claimed by UsageAnalysis' own step, which this project does not
   load: the same reading is here, and for the Variables panel, where a run from the context panel puts its outputs. */
import type {Page} from '@playwright/test';
import {Given, Then, element} from '@datagrok-libraries/bdd';
import {atFeatureEnd, expect, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/** The Scripts gallery's search text is an account setting: a clip that searches it leaves the account's next
 * visit filtered unless it is emptied again (the library's `user opens the Scripts view` empties it the same way). */
export const scriptsSearchPutBack = Given('the Scripts gallery search is emptied at feature end', async (page: Page) => {
  atFeatureEnd(page, async () => {
    await page.evaluate(() => {
      const open = Array.from(grok.shell.views as Iterable<any>).find((v) => v.type === 'scripts');
      if (open)
        grok.shell.v = open;
      else
        grok.shell.route('/scripts');
    });
    await page.waitForFunction(() => grok.shell.v?.type === 'scripts', null, {timeout: pollMs(30000)});
    const search = page.locator('.grok-gallery-search-bar input').filter({visible: true}).first();
    await search.waitFor({timeout: pollMs(30000)});
    if (await search.inputValue() !== '')
      await search.fill('');
    await expect.poll(() => page.locator('.grok-items-view-counts').filter({visible: true}).first().textContent()
      .then((t) => (t ?? '').trim()).catch(() => ''), {message: 'the Scripts gallery counter once its search is empty',
      timeout: pollMs(60000)}).toMatch(/^\d+ of \d+$|^\d{2,}$/);
    await page.evaluate(() => grok.shell.v?.type === 'scripts' && grok.shell.v.close());
  });
}, {tier: 'api', description: 'the Scripts view\'s search box is emptied when the feature ends, as the account had it'});

element('console command line', {selector: '.d4-console-wrapper .d4-console-input input',
  aliases: ['console input'],
  description: 'the "> Enter command" line at the bottom of the console: Enter runs what it holds'});

/* A function's Run pane in the context panel ends with a RUN button; a button phrase matches text without case, so
   "RUN button in Run pane" also finds the pane's own header ("Run", role=button) and is ambiguous. */
element('context panel run button', {selector: '.grok-prop-panel [name="button-RUN"]',
  description: 'the RUN button of the Run pane a function or script shows in the context panel: runs it with the pane\'s inputs'});

/** A script's Script pane in the context panel shows its code in a plain text area (no code editor), whose text is
 * its value, not its text content. */
export const scriptPaneShows = Then('the Script pane of the context panel should show {string}', async (page: Page, text: string) => {
  await expect.poll(() => page.evaluate(() => {
    const content = document.querySelector('.grok-prop-panel [name="div-section--Script"]')?.nextElementSibling;
    if (!content || (content as HTMLElement).offsetParent === null)
      return null;
    const area = content.querySelector('textarea');
    return area ? area.value : content.textContent;
  }), {message: 'the code the Script pane of the context panel shows (null: the pane is not open)', timeout: pollMs(15000)})
    .toContain(text);
}, {description: 'the code in the open Script pane of a script\'s context panel, by substring'});

/** The value beside a name in a visible name / value table (`.d4-item-table`), on the whole page or inside the dock
 * panel titled `panel`: null while there is no such table. */
function readNameValue(page: Page, panel: string | null, name: string): Promise<string | null> {
  return page.evaluate(([p, n]) => {
    const scopes = p === null ? [document.body] : Array.from(document.querySelectorAll('.panel-base'))
      .filter((e) => e.querySelector('.panel-titlebar-text')?.textContent?.trim() === p);
    const tables = scopes.flatMap((e) => Array.from(e.querySelectorAll('.d4-item-table')));
    const table = tables.find((t) => (t as HTMLElement).offsetParent !== null && /name\s+value/.test((t as HTMLElement).innerText));
    if (!table)
      return null;
    const row = Array.from(table.querySelectorAll('tr')).map((r) => Array.from(r.querySelectorAll('td')).map((c) => c.textContent?.trim() ?? ''))
      .find((cells) => cells.includes(n));
    return row ? row[row.indexOf(n) + 1] ?? '' : '';
  }, [panel, name] as [string | null, string]);
}

/** The name / value table the script view puts under the code after a run. */
export const scriptResultReads = Then('the script result {string} should read {string}', async (page: Page, output: string, value: string) => {
  await expect.poll(() => readNameValue(page, null, output),
    {message: `the value of "${output}" in the script results (null: no results table yet)`, timeout: pollMs(30000)}).toBe(value);
}, {description: 'the Results table under the script editor after a run: the value shown beside the output\'s name'});

/** The Variables panel (the status bar's Variables toggle) lists the session's variables, where a function run
 * from the context panel puts its outputs. */
export const variableReads = Then('the Variables panel should show {string} as {string}', async (page: Page, name: string, value: string) => {
  await expect.poll(() => readNameValue(page, 'Variables', name),
    {message: `the value of "${name}" in the Variables panel (null: no table yet)`, timeout: pollMs(30000)}).toBe(value);
}, {description: 'the value beside the variable\'s name in the Variables panel'});
