/**
 * Playwright session helpers.
 *
 * Namespace per helpers-candidates.yaml: `helpers.playwright.session.*`.
 * Imported as: `import {logoutAndLoginAs} from '../helpers/session';`.
 *
 * Canonical home for session-related helpers (re-auth, future switchUser /
 * refreshAuth). Re-authentication is token-based: the same token-injection
 * path used by initial login (`spec-login.ts`) is re-run with a different
 * token. No login form, no password, runs headless.
 *
 * UI-layer reference: `grok-browser/references/projects.md` "Re-auth pattern".
 */

import {Page, expect} from '@playwright/test';
import {loginToDatagrok, loginAsSecondUser} from '../spec-login';

const currentLogin = (page: Page): Promise<string | null> =>
  page.evaluate(() => (window as any).grok?.shell?.user?.login ?? null);

/**
 * Re-authenticate as the primary or second user via token injection.
 *
 * The session is switched by re-injecting the corresponding auth token
 * (`DATAGROK_AUTH_TOKEN` for primary, `DATAGROK_AUTH_TOKEN_2` for second),
 * which overwrites the existing cookie / localStorage entry and reloads.
 * This is fully automatable — no UI logout form and no manual credential
 * entry. Switching to the second user requires `DATAGROK_AUTH_TOKEN_2`
 * to be set (see `loginAsSecondUser` in spec-login.ts).
 *
 * @param page - Playwright Page (already logged in).
 * @param target.as - Which user to re-authenticate as: 'primary' or 'second'.
 * @param options.restoreOriginal - After switching, re-authenticate back as
 *   the other user to restore the baseline session. Default: false.
 */
export async function logoutAndLoginAs(
  page: Page,
  target: {as: 'primary' | 'second'},
  options?: {restoreOriginal?: boolean},
): Promise<void> {
  const before = await currentLogin(page);

  await (target.as === 'second' ? loginAsSecondUser : loginToDatagrok)(page);
  await expect.poll(() => currentLogin(page),
    {timeout: 30000, intervals: [1000, 2000, 3000]}).not.toBe(before);

  if (options?.restoreOriginal) {
    await (target.as === 'second' ? loginToDatagrok : loginAsSecondUser)(page);
    await expect.poll(() => currentLogin(page),
      {timeout: 30000, intervals: [1000, 2000, 3000]}).toBe(before);
  }
}
