import {Page} from '@playwright/test';
import {baseUrl, getSecondUserLogin, resolveSecondUserToken} from '../spec-login';

/**
 * The second user, on a page of its own, booted once per worker.
 *
 * `loginAsSecondUser` re-injects a token into the shared owner page — two navigations and a
 * boot — and `loginToDatagrok` cannot switch back at all: its `alreadyUp` early return
 * (spec-login.ts) sees a booted shell and returns without injecting the owner token. Every
 * "owner" block after the first recipient block therefore ran as the recipient, which is what
 * "You don't have a permission to share this object" and the 90s `waitForIdentity` waits were.
 * Holding the recipient on a second context keeps both identities live at once and costs one
 * boot for the whole section instead of one per block.
 */
let recipient: Page | null = null;

/** The dev-key exchange is a bare Node fetch with no retry, and it failed twice in one run. */
async function retry<T>(what: () => Promise<T>): Promise<T> {
  for (let attempt = 0; ; attempt++) {
    try { return await what(); }
    catch (e) {
      if (attempt === 2) throw e;
      await new Promise((r) => setTimeout(r, 2000));
    }
  }
}

export const secondUserLogin = (): Promise<string> => retry(getSecondUserLogin);

export async function recipientPage(owner: Page): Promise<Page> {
  if (recipient && !recipient.isClosed()) return recipient;
  const browser = owner.context().browser();
  if (!browser)
    throw new Error('recipientPage: the owner page has no browser to open a second context on');
  const token = await retry(resolveSecondUserToken);
  const context = await browser.newContext({viewport: {width: 1920, height: 1080}});
  context.setDefaultTimeout(15_000);
  context.setDefaultNavigationTimeout(120_000);
  const p = await context.newPage();
  await p.goto(baseUrl + '/oauth/');
  await context.addCookies([{name: 'auth', value: token, domain: new URL(baseUrl).hostname, path: '/'}]);
  await p.evaluate((t) => window.localStorage.setItem('auth', t), token);
  await p.goto(baseUrl);
  await p.waitForFunction(() => document.querySelector('.grok-preloader') == null, null, {timeout: 120_000});
  await p.locator('[name="Browse"]').waitFor({timeout: 60_000});
  recipient = p;
  return p;
}

/**
 * Opens the PermissionsView for an entity.
 *
 * GROK-20322 removed the Share dialog's "Advanced editor..." link
 * (core/client/xamgle/lib/src/commands/file/share_dataset.dart, the deleted
 * `htmlLink('Advanced editor...')`), so the matrix is reached by its route instead.
 */
export async function openPermissionsView(page: Page, entityId: string): Promise<void> {
  await page.goto(`${baseUrl}/permissions/${entityId}`);
  await page.waitForFunction(() => /\/permissions\/[0-9a-f-]+/.test(window.location.href),
    null, {timeout: 30_000});
}
