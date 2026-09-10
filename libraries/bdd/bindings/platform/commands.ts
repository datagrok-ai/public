/* The top menu and what a command of it does: a path picked by real pointer moves, the
   function call it starts awaited through the platform's own call events, and the columns it
   added to the table read against the columns it started with. A dialog a command shows in
   between is driven with the base steps (`OK button in "Sequence Space" dialog`). */
import {Page} from '@playwright/test';
import {expect} from '../../src/runtime/patience.js';
import {Then, When} from '../../src/registry.js';
import {closeTopMenu, columnsSince, menuNames, openTopMenu, pickTopMenu, visibleLabels, waitCommand} from '../../src/runtime/menus.js';

const PATH = '"Bio > Analyze > Sequence Space..." — the labels as the menu shows them, separated by ">" or "|"';
// the time a command may run before that is a platform failure (a WASM alignment, a dimensionality reduction)
const COMMAND_CAP = 120000;

export const pickFromTopMenu = When('user picks {string} from the top menu', (page: Page, path: string) => pickTopMenu(page, path),
  {tier: 'ui', description: `${PATH}; the groups open under the pointer, the leaf is clicked; the function call it starts is watched for "the top menu command should have completed"`});

export const openInTopMenu = When('user opens {string} in the top menu', (page: Page, path: string) => openTopMenu(page, path, false),
  {tier: 'ui', description: 'a group of the top menu, left open — its items are then "\\"Bio > Analyze > MSA...\\" menu item"'});

export const closeTheTopMenu = When('user closes the top menu', (page: Page) => closeTopMenu(page), {tier: 'ui'});

/** The paths grouped by their parent: a group opens once, its leaves are checked, the menu
 * closes — as many walks as groups, not as leaves. */
export const topMenuLists = Then('the top menu should list:', async (page: Page, rows: string[][]) => {
  const groups = new Map<string, string[]>();
  for (const [path] of rows) {
    const {segments} = menuNames(path);
    if (segments.length < 2)
      throw new Error(`"${path}" names a group, not a command: a command is "Group > Item"`);
    const group = segments.slice(0, -1).join(' > ');
    groups.set(group, [...(groups.get(group) ?? []), path]);
  }
  for (const [group, paths] of groups) {
    await openTopMenu(page, group, false);
    for (const path of paths) {
      const {segments, names} = menuNames(path);
      await page.locator(`[name="${names[names.length - 1]}"]`).first().waitFor({state: 'visible', timeout: 5000}).catch(async () => {
        throw new Error(`no "${segments[segments.length - 1]}" in the ${group} menu; it shows: ${await visibleLabels(page, names[names.length - 2]) || 'nothing'}`);
      });
    }
    await closeTopMenu(page);
  }
}, {tier: 'ui', description: `one path per row (${PATH}); every group opens once and each of its leaves is found visible, then the menu closes`});

export const commandCompleted = Then('the top menu command should have completed', async (page: Page) => {
  await waitCommand(page, COMMAND_CAP);
}, {description: 'the function call the last picked menu item started has ended (the platform\'s onAfterRunAction), dialog and all — up to two minutes'});

/** The columns the current table gained since the last menu command; polled, since a command
 * that returned at once (a dialog of its own) may still be at work. */
async function expectAdded(page: Page, predicate: (added: string[]) => string | null, what: string): Promise<void> {
  let shown = '';
  try {
    await expect.poll(async () => {
      const c = await columnsSince(page);
      if (c.before === null)
        throw new Error('no menu command has been picked in this scenario, so there is no "before"');
      if (!c.same)
        throw new Error(`the current table is not the one the command started on; it has: ${c.now.join(', ')}`);
      const added = c.now.filter((n) => !c.before!.includes(n));
      shown = added.length === 0 ? 'no column was added' : `added: ${added.join(', ')}`;
      return predicate(added) ?? 'ok';
    }, {timeout: 60000}).toBe('ok');
  }
  catch (e) {
    if (!shown)
      throw e;
    throw new Error(`${what} — ${shown} (within 60 s of the command)`);
  }
}

export const newColumnsCount = Then('{int} new column(s) should have been added', (page: Page, count: number) =>
  expectAdded(page, (added) => added.length === count ? null : `${added.length} added`, `${count} new column(s)`),
  {description: 'to the current table since the last menu command, by name; waits for a command still at work'});

export const newColumnNamed = Then('a new column {string} should have been added', (page: Page, name: string) =>
  expectAdded(page, (added) => added.includes(name) ? null : 'missing', `a new column "${name}"`));

export const newColumnMatching = Then('a new column matching {string} should have been added', (page: Page, pattern: string) =>
  expectAdded(page, (added) => added.some((n) => new RegExp(pattern).test(n)) ? null : 'missing', `a new column matching /${pattern}/`),
  {description: 'a regular expression over the names added — for a name the platform suffixes (getUnusedName)'});

export const noNewColumn = Then('no new column should have been added', async (page: Page) => {
  const c = await columnsSince(page);
  if (c.before === null)
    throw new Error('no menu command has been picked in this scenario, so there is no "before"');
  expect(c.now.filter((n) => !c.before!.includes(n)), 'columns added since the menu command').toEqual([]);
}, {description: 'read once the command has completed: put "the top menu command should have completed" first'});
