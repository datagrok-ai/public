/* The shell's workspace and what surrounds it: the left sidebar's Dashboards panel, the tables and
   table views that are open, how a project open ended, what a dialog says, the time between two
   gestures, and files on the server (a path, a folder of the user's own files). */
import {type Page} from '@playwright/test';
import {expect, pollMs} from '../../src/runtime/patience.js';
import {Given, Then, When} from '../../src/registry.js';
import {atFeatureEnd} from '../../src/runtime/harness.js';
import {exactText} from '../../src/runtime/locate.js';

declare const grok: any;

const listOf = (text: string): string[] => text.split(',').map((s) => s.trim()).filter(Boolean);

/* --- the Dashboards panel ------------------------------------------------------------------------
   The left sidebar's Dashboards tab lists "New Dashboard" (the scratchpad) and every open project,
   each with its tables; tables are added to an open project by a drop on its row. The tab toggles
   its panel, so it is clicked only when the panel is not shown. The panel stays open beside the
   Browse panel into the next feature (its New Dashboard row then brings a second "Save" button to
   every table view), so it is closed again when the feature ends. */
const dashboardsClosedAtEnd = new WeakSet<Page>();

export const dashboardsPanelOpen = Given('the dashboards panel of the left sidebar is open', async (page: Page) => {
  const tab = page.locator('[name="sidebar"] [name="Dashboards"]').first();
  const row = page.locator('.grok-view-browse [name="tree-New-Dashboard"]').filter({visible: true}).first();
  if (!(await tab.evaluate((e) => e.classList.contains('selected'))) || await row.count() === 0)
    await tab.click();
  await expect(row, 'the New Dashboard row of the Dashboards panel').toBeVisible({timeout: pollMs(15000)});
  if (dashboardsClosedAtEnd.has(page))
    return;
  dashboardsClosedAtEnd.add(page);
  // first: before the browse panel's own restore puts the shell back in simple mode
  atFeatureEnd(page, async () => {
    dashboardsClosedAtEnd.delete(page);
    await page.evaluate(() => {
      const t = document.querySelector('[name="sidebar"] [name="Dashboards"]') as HTMLElement | null;
      if (t?.classList.contains('selected'))
        t.click();
    });
    await expect(page.locator('[name="sidebar"] [name="Dashboards"]').first(), 'the Dashboards tab of the left sidebar, closed at feature end')
      .not.toHaveClass(/\bselected\b/);
  }, true);
}, {tier: 'ui', description: 'clicks the Dashboards tab of the left sidebar unless its panel is already shown; done when the New Dashboard row is visible; the panel is closed again when the feature ends'});

/* --- the tables and table views that are open ------------------------------------------------------ */

export const openTablesExactly = Then('the open tables should be exactly {string}', async (page: Page, names: string) => {
  await expect.poll(() => page.evaluate(() => (grok.shell.tables as any[]).map((t) => String(t.name)).sort().join(', ')),
    {message: 'the names of the open tables, sorted', timeout: pollMs(30000)}).toBe(listOf(names).sort().join(', '));
}, {description: 'grok.shell.tables by name (comma-separated), whatever view shows them: these and no others'});

export const noTableLeft = Then('no table should be left in the workspace', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => (grok.shell.tables as any[]).map((t) => t.name).join(' | ')),
    {message: 'the tables still open in the workspace'}).toBe('');
}, {description: 'grok.shell.tables is empty, polled — the claim after a Close All made through the UI, before anything is reopened'});

export const tableViewsOpen = Then('the table views {string} should be open', async (page: Page, names: string) => {
  const wanted = listOf(names);
  await expect.poll(() => page.evaluate((w) => {
    const open = (Array.from(grok.shell.views) as any[]).filter((v) => v.type === 'TableView').map((v) => String(v.name));
    const missing = w.filter((n) => !open.includes(n));
    return missing.length === 0 ? 'all open' : `missing ${missing.join(' | ')}; open: ${open.join(' | ') || 'none'}`;
  }, wanted), {message: 'the open table views', timeout: pollMs(60000)}).toBe('all open');
}, {description: 'every named table view is open (comma-separated), whatever else is'});

/* --- how an open ended -----------------------------------------------------------------------------
   A project open ends one of two ways: its table view in front with the rows its creation script
   produced (the data-sync mark), or an error dialog. The barrier waits for either and claims neither,
   so a scenario after it reads a settled state — a @known-failure one included, whose narrowed
   budget would otherwise be spent on the open itself. */
async function openOutcome(page: Page, view: string, dialog: string): Promise<string> {
  const shown = await page.locator('.d4-dialog').filter({visible: true})
    .filter({has: page.locator('.d4-dialog-title', {hasText: exactText(dialog)})}).count();
  if (shown > 0)
    return `the "${dialog}" dialog`;
  return page.evaluate((v) => {
    const df = grok.shell.tv?.dataFrame;
    return grok.shell.v?.name === v && df?.rowCount > 0 && df.getTag('.data-sync') === 'success' ?
      `the "${v}" view` : `neither yet (current view "${grok.shell.v?.name}", ${df?.rowCount ?? 'no'} rows, mark ${df?.getTag('.data-sync') ?? 'none'})`;
  }, view);
}

export const openEndsIn = Then('the open should end in the {string} view or the {string} dialog', async (page: Page, view: string, dialog: string) => {
  await expect.poll(() => openOutcome(page, view, dialog), {message: 'how the open ended', timeout: pollMs(90000)}).toMatch(/^the "/);
}, {description: 'a barrier: either the view is current with rows re-read by data sync, or the dialog is showing — whichever comes first'});

export const noDialogShows = Then('no {string} dialog should show {string}', async (page: Page, dialog: string, text: string) => {
  const hits = page.locator('.d4-dialog').filter({visible: true})
    .filter({has: page.locator('.d4-dialog-title', {hasText: dialog})}).filter({hasText: text});
  await expect(hits, `a visible "${dialog}" dialog showing ${text}`).toHaveCount(0);
}, {description: 'no visible dialog of that title holds the text; read once, so it follows a barrier such as "the open should end in …"'});

/* --- the time between two gestures -----------------------------------------------------------------
   A defect that shows only when two gestures come close enough together: the time of the first is
   noted, and the second claims the gap, so a slow run fails as itself instead of hiding the defect. */
const notedTimes = new Map<string, number>();

export const noteTime = When('user notes the time as {string}', async (_page: Page, label: string) => {
  notedTimes.set(label, Date.now());
}, {tier: 'api', description: 'the wall-clock time under a label, for "at most … seconds should have passed since"'});

export const atMostSince = Then('at most {int} seconds should have passed since {string}', async (_page: Page, seconds: number, label: string) => {
  const at = notedTimes.get(label);
  if (at === undefined)
    throw new Error(`no time noted as "${label}"`);
  expect(Math.round((Date.now() - at) / 1000), `seconds since "${label}"`).toBeLessThanOrEqual(seconds);
}, {description: 'a read of the gap, not a wait: it fails when the run was slower than the claim that follows needs'});

/* --- files on the server ---------------------------------------------------------------------------- */

export const fileOnServer = Then('the file {string} should be on the server', async (page: Page, path: string) => {
  await expect.poll(() => page.evaluate((p) => grok.dapi.files.exists(p), path),
    {message: `the file ${path} on the server`, timeout: pollMs(30000)}).toBe(true);
}, {tier: 'api', description: 'grok.dapi.files.exists on the path (<namespace>:<share>/<file>)'});

/* Browse > Files > My files is the account's Home storage, "<namespace>:Home" (the namespace is the
   user's own project, "Admin" for the admin login). A feature's folder there is written through the
   JS API — the files are the scene, not the subject — and removed with everything in it before the
   feature starts and when it ends, read back from the server. */
function homeFolder(page: Page, folder: string): Promise<string> {
  return page.evaluate((f) => `${grok.shell.user.project.name}:Home/${f}/`, folder);
}

/* The folder named, and any folder of its family (the same letters before a 13-digit time) whose
   time is over an hour old: what a killed run left. */
async function removeHomeFolder(page: Page, folder: string): Promise<void> {
  const path = await homeFolder(page, folder);
  const home = path.slice(0, path.indexOf('/') + 1);
  const family = /^(.*[A-Za-z])\d{13,}$/.exec(folder)?.[1] ?? null;
  const doomed = (names: string[]): string[] => names.filter((n) => n === folder ||
    (family != null && n.startsWith(family) && /^\d{13,}$/.test(n.slice(family.length)) && Date.now() - Number(n.slice(family.length)) > 60 * 60 * 1000));
  const folders = () => page.evaluate(async (h) =>
    (await grok.dapi.files.list(h, false)).filter((f: any) => f.isDirectory).map((f: any) => String(f.name ?? f.fileName)), home);
  for (const name of doomed(await folders()))
    await page.evaluate(async (p) => {
      for (const f of await grok.dapi.files.list(p, true))
        if (!f.isDirectory)
          await grok.dapi.files.delete(f.fullPath);
      await grok.dapi.files.delete(p.slice(0, -1));
    }, `${home}${name}/`);
  await expect.poll(async () => doomed(await folders()).join(' | ') || 'gone',
    {message: `the folder ${path} (and old ones of its family) in the user's files`, timeout: pollMs(15000)}).toBe('gone');
}

export const noHomeFolder = Given('no folder {string} is in the user\'s files', async (page: Page, folder: string) => {
  atFeatureEnd(page, () => removeHomeFolder(page, folder));
  await removeHomeFolder(page, folder);
}, {tier: 'api', description: 'removes that folder of Browse > Files > My files with its files, read back from the server; again when the feature ends'});

export const homeFileWritten = Given('the file {string} of the folder {string} in the user\'s files holds:', async (page: Page, file: string, folder: string, text: string) => {
  const path = await homeFolder(page, folder);
  const listed = await page.evaluate(async ([p, f, t]) => {
    await grok.dapi.files.writeAsText(p + f, t);
    return (await grok.dapi.files.list(p, false)).map((x: any) => x.fileName).sort();
  }, [path, file, text] as [string, string, string]);
  expect(listed, `the files of ${path}`).toContain(file);
}, {tier: 'api', description: 'written through the JS API into Browse > Files > My files/<folder> (made if missing); "no folder … is in the user\'s files" removes it'});
