/**
 * Runnable test for `helpers/session.ts` helper `logoutAndLoginAs`.
 *
 * Token-based re-auth — runs headless. Switching to the second user requires
 * a second-user token; the test skips when `DATAGROK_AUTH_TOKEN_2` is not set.
 * Run via:
 *
 *   cd public/packages/UsageAnalysis/files/TestTrack
 *   DATAGROK_URL=https://dev.datagrok.ai \
 *     npx playwright test --project=dev helpers/session-spec.ts
 *
 * (`DATAGROK_AUTH_TOKEN` / `DATAGROK_AUTH_TOKEN_2` are provided by the runner.)
 *
 * Test flow:
 *   1. Login as the primary user (DATAGROK_AUTH_TOKEN).
 *   2. Call logoutAndLoginAs({as: 'second'}) to re-authenticate as the second
 *      user via token injection.
 *   3. Verify shell.user.login changed.
 *   4. Call logoutAndLoginAs({as: 'primary'}) to restore the baseline session.
 *   5. Verify shell.user.login is back to the original.
 *
 * No fixture state created — operates purely on the authenticated session.
 */

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {logoutAndLoginAs} from './session';

test.use({...specTestOptions, navigationTimeout: 180_000});

// Two-user test: loginAsSecondUser resolves the second token from
// DATAGROK_AUTH_TOKEN_2 or the config `key2:`, and throws if neither exists —
// so this fails (rather than skips) when no second user is configured.
const readLogin = (page: any) =>
  page.evaluate(() => (window as any).grok?.shell?.user?.login ?? null);

test('helpers.playwright.session — logoutAndLoginAs end-to-end', async ({page}) => {
  // Generous timeout — each re-auth reloads and waits for the Browse panel.
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  const originalLogin = await readLogin(page);
  expect(originalLogin).toBeTruthy();

  await softStep('logoutAndLoginAs — switch to second user (token injection)', async () => {
    await logoutAndLoginAs(page, {as: 'second'});
    const newLogin = await readLogin(page);
    expect(newLogin).not.toBe(originalLogin);
    expect(newLogin).toBeTruthy();
  });

  await softStep('logoutAndLoginAs — restore primary user (token injection)', async () => {
    await logoutAndLoginAs(page, {as: 'primary'});
    const restored = await readLogin(page);
    expect(restored).toBe(originalLogin);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} session step(s) failed:\n${summary}`);
  }
});
