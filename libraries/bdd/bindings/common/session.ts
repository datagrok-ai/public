/* The session is the storage state global setup writes (@datagrok-libraries/test); this step only
   lands in the shell — once per browser page, the scenarios of a feature share it — and applies the
   automation-friendly shell settings every suite applies, plus the viewer runtime: every viewer the
   page holds or adds later renders immediately (no debounce, no animation-frame wait). */
import type {BrowserContext, Page} from '@playwright/test';
import {Given, When} from '../../src/registry.js';
import {takeErrors} from '../../src/runtime/harness.js';
import {expect, pollMs} from '../../src/runtime/patience.js';
import {installViewerRuntime, takeBalloons} from '../../src/runtime/viewers.js';
import * as guide from '../../src/runtime/guide.js';

declare const grok: any;

const homeNotWaitedFor = new WeakSet<Page>();

/* What a PowerPack Home widget logs while it loads would land on whichever scenario runs by then. The
   widgets host shows up within a second of the shell; a stand without it is not waited for again. */
async function homeWidgetsSettled(page: Page): Promise<void> {
  if (homeNotWaitedFor.has(page))
    return;
  const hasHost = await page.waitForFunction(() => grok.shell.v?.root?.querySelector('.power-pack-widgets-host') != null,
    null, {timeout: 10000}).then(() => true, () => false);
  const settled = hasHost && await page.waitForFunction(() => {
    const contents = [...grok.shell.v.root.querySelectorAll('.power-pack-widgets-host .power-pack-widget-content')] as HTMLElement[];
    return contents.length > 0 && contents.every((c) => c.children.length > 0 && c.querySelector('.grok-loader') == null);
  }, null, {timeout: 30000}).then(() => true, () => false);
  if (!settled) {
    homeNotWaitedFor.add(page);
    console.warn(`bdd: ${hasHost ? 'the Home widgets did not finish loading in 30 s' : 'no PowerPack Home widgets'}; not waited for again on this page`);
  }
}

/** A freshly loaded shell made ready for a feature: nothing open, the Home view settled and the in-page
 * runtime installed. What the page logged while it booted stays on the error floor: a widget that fails
 * to load is the claim of a scenario that reloads or signs in. */
export async function resetShellAfterLoad(page: Page, simpleMode = guide.shellSimpleMode()): Promise<void> {
  await page.evaluate((simple) => {
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = simple;
    // the console and help panels a feature leaves open share the right column with the context panel and can
    // squeeze it to its title bar; a feature that needs one opens it
    grok.shell.windows.showConsole = false;
    grok.shell.windows.showHelp = false;
  }, simpleMode);
  // closeAll re-adds the Home view asynchronously; a table opened before it lands ends up behind it
  await page.waitForFunction(() => grok.shell.v?.type === 'datagrok', null, {timeout: 60000});
  await homeWidgetsSettled(page);
  await installViewerRuntime(page);
}

/** The shell mode the page runs in now: a reload or a sign-in mid-feature keeps what the feature set. */
export const currentSimpleMode = (page: Page): Promise<boolean> =>
  page.evaluate(() => grok.shell.windows.simpleMode as boolean).catch(() => guide.shellSimpleMode());

/** Another session on the same page, swapped the way a sign-out and a sign-in swap it: the auth cookie
 * and the stored token replaced, and what a sign-out clears (Auth.cleanup in user.dart) cleared — the
 * function list the client cached for the account before, which the next account would otherwise start
 * from. The cache is deleted from a document of the same origin the platform does not run in, so no open
 * connection holds the deletion up. */
export async function signInWithSession(page: Page, token: string, login: string, simple?: boolean): Promise<void> {
  simple ??= await currentSimpleMode(page);
  const origin = new URL(page.url()).origin;
  await page.goto(`${origin}/favicon.ico`, {waitUntil: 'load', timeout: 60000});
  await page.context().clearCookies({name: 'auth'});
  await page.context().addCookies([{name: 'auth', value: token, domain: new URL(origin).hostname, path: '/'}]);
  await page.evaluate((t) => new Promise<void>((resolve, reject) => {
    localStorage.setItem('auth', t);
    const request = indexedDB.deleteDatabase('CachedFuncs');
    request.onsuccess = () => resolve();
    request.onerror = () => reject(new Error(`the client's function cache was not deleted: ${request.error}`));
  }), token);
  await page.goto(`${origin}/`, {waitUntil: 'domcontentloaded', timeout: 180000});
  await page.locator('[name="Browse"]').first().waitFor({timeout: 180000});
  await expect.poll(() => page.evaluate(() => grok.shell.user?.login ?? ''),
    {message: 'the login of the account the shell runs as', timeout: pollMs(30000)}).toBe(login);
  await resetShellAfterLoad(page, simple);
}

// the account each browser context started with: a feature that signed in as another and could not
// sign back leaves its context there, and the next feature on the page must not run as that account
const startedAs = new WeakMap<BrowserContext, {token: string; login: string}>();

export const loggedIn = Given('user is logged in', async (page: Page) => {
  // a worker runs one spec after another on the same page, so what a feature leaves behind (an open
  // dialog, a docked panel, a sticky option) reaches the next one; BDD_FRESH_PAGE starts each
  // feature from a reload, at the cost of a shell load per feature
  guide.silent(page);
  const inShell = process.env.BDD_FRESH_PAGE !== '1' &&
    await page.evaluate(() => typeof (window as any).grok?.shell?.closeAll === 'function').catch(() => false);
  if (!inShell) {
    // a dev stand's pub serve can take minutes to hand out the bundle while it recompiles or is
    // starved: that is a delay once per page, not a failure of the feature
    const start = Date.now();
    await page.goto('/', {waitUntil: 'domcontentloaded', timeout: 180000});
    await page.locator('[name="Browse"]').first().waitFor({timeout: 180000});
    const seconds = Math.round((Date.now() - start) / 1000);
    if (seconds >= 30)
      console.warn(`bdd: the shell took ${seconds} s to load (a dev stand serving a bundle it is recompiling?)`);
  }
  const now = await page.evaluate(() => ({token: localStorage.getItem('auth') ?? String(grok.dapi.token ?? ''),
    login: String(grok.shell.user?.login ?? '')}));
  const first = startedAs.get(page.context());
  if (!first)
    startedAs.set(page.context(), now);
  if (first && first.login !== now.login) {
    console.warn(`bdd: the page was left signed in as ${now.login}; signing ${first.login} back in`);
    await signInWithSession(page, first.token, first.login, guide.shellSimpleMode());
  }
  else
    await resetShellAfterLoad(page);
  // what the stand logs or shows while booting (a broken package's autostart, "Debugging packages") is not the scenario's
  takeErrors(page);
  await takeBalloons(page);
}, {tier: 'ui', description: 'the error floor starts here: "no errors should have been logged" counts from this step; a page an earlier feature left signed in as another account signs the running one back in first'});

/** The browser's reload: a new session of the shell with nothing in memory, for a reopen that has
 * to come from the server alone. The panels and the gallery a feature opened are gone with it. */
export const reloadPage = When('user reloads the page', async (page: Page) => {
  const simple = await currentSimpleMode(page);
  await page.reload({waitUntil: 'domcontentloaded', timeout: 180000});
  await page.locator('[name="Browse"]').first().waitFor({timeout: 180000});
  await resetShellAfterLoad(page, simple);
}, {tier: 'ui', description: 'page.reload, then the shell set up as "user is logged in" sets it up, in the mode the feature ran in; what the page logs while it loads stays on the error floor'});
