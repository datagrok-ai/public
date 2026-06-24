/**
 * Playwright session helpers.
 *
 * Imported as: `import {logoutAndLoginAs} from '../helpers/session';`.
 */

import {Page, expect} from '@playwright/test';

// ---------------------------------------------------------------------------
// logoutAndLoginAs — UI logout + login-as-different-user (manual password).
// ---------------------------------------------------------------------------

/**
 * Drive the UI logout + login-as-different-user sequence with optional
 * restore-original.
 *
 * **Manual password entry required.** This helper accepts only a `username`
 * — passwords are NEVER passed as parameters. The helper drives the flow
 * up to the point where the password input has focus, then **pauses** for
 * the human operator to type the password into the browser. Once the
 * post-login Browse panel appears, the helper resumes.
 *
 * The signature deliberately excludes a `password` field for security. For
 * automated end-to-end credential flows, a future `refreshAuth` helper backed
 * by `auth-<env>.json` storage state would be the right path.
 *
 * Flow:
 *   1. Locate user-avatar / user-icon → trigger logout
 *   2. Wait for login form to render
 *   3. Fill `username` (helper-driven)
 *   4. **Pause for manual password entry** (waits up to 5 minutes for
 *      `[name="Browse"]` to appear, indicating successful login)
 *   5. If `restoreOriginal: true` — repeat steps 1-4 with the original
 *      user (the human must remember the original credentials and enter
 *      that password too)
 *
 * **Headed mode required.** This helper cannot run headlessly because
 * the password entry is a manual step. Set up the test to run via
 * `npx playwright test --project=dev <test> --headed` (omit `CI=1` env).
 *
 * **Headed-browser conflict warning.** Headed Playwright runs can collide
 * with another browser instance attached to the same machine — ensure no
 * conflicting browser session is active before running.
 *
 * @param page - Playwright Page (already logged in).
 * @param credentials.username - Target login (the password will be typed
 *   manually by the operator).
 * @param options.restoreOriginal - After logging in as the new user and
 *   verifying the helper's caller has finished its flow, log out again and
 *   log back in as the original user. Default: false.
 */
export async function logoutAndLoginAs(
  page: Page,
  credentials: {username: string},
  options?: {restoreOriginal?: boolean},
): Promise<void> {
  // Step 1 — Click user avatar / icon to open the user menu, then click Logout.
  // The user avatar lives on the sidebar ([name="User"]); we locate the
  // "Logout" item by text as a defensive fallback.
  const originalUsername = await page.evaluate(() => {
    const grok = (window as any).grok;
    return grok?.shell?.user?.login ?? null;
  });

  await driveLogout(page);
  await fillUsernameAndAwaitManualPassword(page, credentials.username);

  // Verify the new login took effect (read shell.user.login).
  await expect.poll(async () => page.evaluate(() => {
    const grok = (window as any).grok;
    return grok?.shell?.user?.login ?? null;
  }), {timeout: 30000, intervals: [1000, 2000, 3000]}).not.toBe(originalUsername);

  if (options?.restoreOriginal && originalUsername) {
    await driveLogout(page);
    await fillUsernameAndAwaitManualPassword(page, originalUsername);
    await expect.poll(async () => page.evaluate(() => {
      const grok = (window as any).grok;
      return grok?.shell?.user?.login ?? null;
    }), {timeout: 30000, intervals: [1000, 2000, 3000]}).toBe(originalUsername);
  }
}

/**
 * Drive the logout side of the flow. The Datagrok user avatar / icon lives
 * in the sidebar; clicking it opens a panel with a Logout option.
 */
async function driveLogout(page: Page): Promise<void> {
  // Try the sidebar User tab first.
  const userTab = page.locator('[name="User"]');
  if (await userTab.isVisible({timeout: 5000}).catch(() => false))
    await userTab.click({force: true});
  // Click any visible Logout / Log out item.
  const logoutItem = page.locator('text=/^(Log out|Logout)$/i').first();
  await logoutItem.waitFor({timeout: 10000});
  await logoutItem.click({force: true});
  // Wait for the login form (one of its inputs) to appear.
  const loginInput = page.locator('#signup-login-fields input[placeholder="Login or Email"]');
  await loginInput.waitFor({timeout: 30000});
}

/**
 * Fill the username field, focus the password field, then wait up to 5
 * minutes for the post-login Browse panel to appear (indicating the human
 * operator has typed the password and submitted).
 */
async function fillUsernameAndAwaitManualPassword(page: Page, username: string): Promise<void> {
  const loginInput = page.locator('#signup-login-fields input[placeholder="Login or Email"]');
  await loginInput.focus();
  await page.keyboard.type(username);
  await page.keyboard.press('Tab');
  // Password input now has focus; the human operator types the password and
  // presses Enter. Wait for the post-login indicator: the Browse sidebar tab.
  // 5-minute window so Olena has time to type.
  await page.locator('[name="Browse"]').waitFor({timeout: 300000});
}
