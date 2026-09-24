/* The top menu and what a command of it does: a path picked by real pointer moves, the
   function call it starts awaited through the platform's own call events, and the columns it
   added to the table read against the columns it started with. A dialog a command shows in
   between is driven with the base steps (`OK button in "Sequence Space" dialog`). */
import {Page} from '@playwright/test';
import {expect, pollMs} from '../../src/runtime/patience.js';
import {Then, When} from '../../src/registry.js';
import {closeTopMenu, columnsSince, menuNames, openTopMenu, pickTopMenu, visibleLabels, waitCommand} from '../../src/runtime/menus.js';

const PATH = '"Bio > Analyze > Sequence Space..." — the labels as the menu shows them, separated by ">" or "|"';
// the time a command may run before that is a platform failure (a WASM alignment, a dimensionality reduction)
const COMMAND_CAP = 120000;

export const pickFromTopMenu = When('user picks {string} from the top menu', (page: Page, path: string) => pickTopMenu(page, path),
  {tier: 'ui', description: `${PATH}; the groups open under the pointer, the leaf is clicked; the function call it starts is watched for "the top menu command should have completed"`});

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
}, {description: 'the function call the last picked menu item started has ended (the platform\'s onAfterRunAction), dialog and all — up to two minutes; fails when that call had ended before its dialog\'s OK, which only a package\'s own dialog opener does — claim what the OK produces instead'});

/** The columns the current table gained since the last menu command; polled, since a command
 * that returned at once (a dialog of its own) may still be at work. */
async function expectAdded(page: Page, predicate: (added: string[]) => string | null, what: string): Promise<void> {
  let shown = '';
  // a command that computes in the browser (a descriptor batch, a toxicity run over a thousand
  // molecules) keeps working well past the shared expect budget
  const budget = pollMs(Number(process.env.BDD_COMMAND_TIMEOUT ?? 180000));
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
    }, {timeout: budget}).toBe('ok');
  }
  catch (e) {
    if (!shown)
      throw e;
    throw new Error(`${what} — ${shown} (within ${Math.round(budget / 1000)} s of the command)`);
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

export const newColumnsMatching = Then('{int} new columns matching {string} should have been added', (page: Page, count: number, pattern: string) =>
  expectAdded(page, (added) => {
    const n = added.filter((name) => new RegExp(pattern).test(name)).length;
    return n === count ? null : `${n} matching`;
  }, `${count} new columns matching /${pattern}/`),
{description: 'counted among the columns added since the last menu command — a second run of a command that suffixes its column name'});

/** The last column of the current table whose name matches: its name, distinct values and missing count. */
async function newestMatching(page: Page, pattern: string): Promise<{name: string; distinct: number; missing: number}> {
  const facts = await page.evaluate((p) => {
    const t = (window as any).grok.shell.t;
    const names: string[] = t ? t.columns.names().filter((n: string) => new RegExp(p).test(n)) : [];
    if (names.length === 0)
      return {name: '', distinct: -1, missing: 0};
    const col = t.col(names[names.length - 1]);
    const values = new Set<string>();
    let missing = 0;
    for (let i = 0; i < t.rowCount; i++) {
      if (col.isNone(i))
        missing++;
      else
        values.add(String(col.get(i)));
    }
    return {name: col.name as string, distinct: values.size, missing};
  }, pattern);
  if (facts.distinct < 0)
    throw new Error(`the current table has no column matching /${pattern}/`);
  return facts;
}

export const newestMatchingDistinct = Then('the newest column matching {string} should have {int} distinct values', async (page: Page, pattern: string, count: number) => {
  const facts = await newestMatching(page, pattern);
  expect(facts.distinct, `distinct values of "${facts.name}" (${facts.missing} missing)`).toBe(count);
}, {description: 'the last column of the current table whose name matches the regular expression; missing values are not a value'});

export const newestMatchingFilled = Then('the newest column matching {string} should have no missing values', async (page: Page, pattern: string) => {
  let seen = '';
  await expect.poll(async () => {
    const facts = await newestMatching(page, pattern);
    seen = `"${facts.name}"`;
    return facts.missing;
  }, {message: `missing values in the newest column matching /${pattern}/`}).toBe(0).catch(() => {
    throw new Error(`missing values in ${seen}`);
  });
}, {description: 'polled: a column a command adds is filled a moment after it appears'});

export const noNewColumn = Then('no new column should have been added', async (page: Page) => {
  const c = await columnsSince(page);
  if (c.before === null)
    throw new Error('no menu command has been picked in this scenario, so there is no "before"');
  expect(c.now.filter((n) => !c.before!.includes(n)), 'columns added since the menu command').toEqual([]);
}, {description: 'read once the command has completed: put "the top menu command should have completed" first'});
