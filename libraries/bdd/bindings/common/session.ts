/* The session is the storage state global setup writes (@datagrok-libraries/test); this step only
   lands in the shell — once per browser page, the scenarios of a feature share it — and applies the
   automation-friendly shell settings every suite applies, plus the viewer runtime: every viewer the
   page holds or adds later renders immediately (no debounce, no animation-frame wait). */
import type {Page} from '@playwright/test';
import {Given, Then} from '../../src/registry.js';
import {atFeatureEnd, takeErrors} from '../../src/runtime/harness.js';
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

/** Lands in the shell (reloading it when the page is not there yet, or `reload` asks) and applies the
 * automation settings; the error and balloon floors start here. */
async function enterShell(page: Page, reload = false): Promise<void> {
  // a worker runs one spec after another on the same page, so what a feature leaves behind (an open
  // dialog, a docked panel, a sticky option) reaches the next one; BDD_FRESH_PAGE starts each
  // feature from a reload, at the cost of a shell load per feature
  guide.silent(page);
  const inShell = !reload && process.env.BDD_FRESH_PAGE !== '1' &&
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
  await page.evaluate((simple) => {
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = simple;
    // the console and help panels a feature leaves open share the right column with the context panel and can
    // squeeze it to its title bar; a feature that needs one opens it
    grok.shell.windows.showConsole = false;
    grok.shell.windows.showHelp = false;
  }, guide.shellSimpleMode());
  // closeAll re-adds the Home view asynchronously; a table opened before it lands ends up behind it
  await page.waitForFunction(() => grok.shell.v?.type === 'datagrok', null, {timeout: 60000});
  await homeWidgetsSettled(page);
  await installViewerRuntime(page);
  // what the stand logs or shows while booting (a broken package's autostart, "Debugging
  // packages") is not the scenario's — and a boot balloon left on screen covers the top right corner
  takeErrors(page);
  await takeBalloons(page);
  await page.evaluate(() => document.querySelectorAll('.d4-balloon').forEach((b) => b.remove()));
}

export const loggedIn = Given('user is logged in', (page: Page) => enterShell(page), {tier: 'ui', description: 'the error floor starts here: "no errors should have been logged" counts from this step'});

/* --- the second account ---------------------------------------------------------------------------
   A feature that looks at the platform through another user's eyes (what is shared with them, what a
   user without a privilege is refused) signs that user in on the same page — never a second page —
   and the first account comes back at feature end. The token is the one global setup resolved into
   DATAGROK_AUTH_TOKEN_2 (a CI runner's own, or a password login of DATAGROK_SHARING_LOGIN); the
   account is the one the sharing steps share with. */

/** The first account's token, taken from the page before it signs the second one in. */
const firstToken = new WeakMap<Page, string>();
const firstLogin = new WeakMap<Page, string>();

async function signInWith(page: Page, token: string): Promise<void> {
  const origin = new URL(page.url() && page.url().startsWith('http') ? page.url() : (process.env.DATAGROK_URL ?? 'http://localhost:8888'));
  await page.context().addCookies([{name: 'auth', value: token, domain: origin.hostname, path: '/'}]);
  await page.evaluate((t) => window.localStorage.setItem('auth', t), token);
  await enterShell(page, true);
}

async function currentLogin(page: Page): Promise<string> {
  return page.evaluate(async () => (await grok.dapi.users.current()).login as string);
}

export const signInAsSecond = Given('user signs in as the second user', async (page: Page) => {
  const token = process.env.DATAGROK_AUTH_TOKEN_2;
  if (!token)
    throw new Error('no second account to sign in as: set DATAGROK_AUTH_TOKEN_2, or DATAGROK_SHARING_LOGIN with DATAGROK_SHARING_PASSWORD');
  if (!firstToken.has(page)) {
    const own = await page.evaluate(() => window.localStorage.getItem('auth'));
    if (!own)
      throw new Error('the page holds no auth token of its own to come back to');
    firstToken.set(page, own);
    firstLogin.set(page, await currentLogin(page));
    atFeatureEnd(page, async () => {
      if (firstToken.has(page))
        await signBackIn(page);
    });
  }
  await signInWith(page, token);
  const login = await currentLogin(page);
  if (process.env.DATAGROK_SHARING_LOGIN && login !== process.env.DATAGROK_SHARING_LOGIN)
    throw new Error(`signed in as "${login}", not the second account "${process.env.DATAGROK_SHARING_LOGIN}"`);
}, {tier: 'ui', description: 'the same page reloads under the second account (the one sharing steps share with); the first comes back at feature end'});

async function signBackIn(page: Page): Promise<void> {
  const token = firstToken.get(page);
  if (!token)
    throw new Error('the second user was never signed in on this page');
  await signInWith(page, token);
  firstToken.delete(page);
}

export const signBackInAsFirst = Given('user signs in again as the first user', (page: Page) => signBackIn(page),
  {tier: 'ui', description: 'the page reloads under the account the run started with'});

export const signedInAs = Then('the second user should be signed in', async (page: Page) => {
  const login = await currentLogin(page);
  if (login !== process.env.DATAGROK_SHARING_LOGIN)
    throw new Error(`the signed-in user is "${login}", not the second account "${process.env.DATAGROK_SHARING_LOGIN}"`);
}, {tier: 'api', description: 'the account the server sees behind the page'});

export const firstSignedIn = Then('the first user should be signed in', async (page: Page) => {
  const login = await currentLogin(page);
  const first = firstLogin.get(page);
  if (!first)
    throw new Error('the second user was never signed in on this page, so there is no first user to compare with');
  if (login !== first)
    throw new Error(`the signed-in user is "${login}", not the first account "${first}"`);
}, {tier: 'api', description: 'the account the run started with, as the server sees it behind the page'});
