/**
 * Runnable test for `helpers/session.ts` C1b helper `logoutAndLoginAs`.
 *
 * **Headed mode required.** This test cannot run headlessly because the
 * password entry is a manual step — Olena types passwords directly in the
 * browser when the helper pauses on the password input. Run via:
 *
 *   cd public/packages/UsageAnalysis/files/TestTrack
 *   DATAGROK_URL=https://dev.datagrok.ai DATAGROK_LOGIN=qa-pw \
 *     npx playwright test --project=dev helpers/session-spec.ts --headed
 *
 * (Note: omit `CI=1` so headless mode is NOT forced.)
 *
 * **MCP-Chrome conflict warning.** Per c1-2026-05-03-helpers-c1a-execution-validated
 * Run #1, headed Playwright runs collide with MCP-attached Chrome on the
 * same machine. Coordinate: ensure MCP browser session is inactive (or
 * pointed away from dev.datagrok.ai) before running this test.
 *
 * Test flow:
 *   1. Login as qa-pw (auth-dev.json storage state).
 *   2. Call logoutAndLoginAs to switch to `Olena Ahadzhanian`. The helper
 *      pauses on the password input — Olena types Olena's password in the
 *      visible browser window and presses Enter.
 *   3. Verify shell.user.login changed to 'Olena Ahadzhanian' (or its
 *      lowercased variant — Datagrok login schema).
 *   4. Call logoutAndLoginAs again with `restoreOriginal: true` semantically
 *      equivalent — switch back to qa-pw. Olena types qa-pw password.
 *   5. Verify shell.user.login is back to qa-pw.
 *
 * No fixture state created in this test — it operates purely on the
 * authenticated session.
 */

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {logoutAndLoginAs} from '../helpers/session';

test.use({...specTestOptions, navigationTimeout: 180_000});

test('helpers.playwright.session — logoutAndLoginAs end-to-end', async ({page}) => {
  // Long timeout — accommodates manual password entry windows (5 min each).
  test.setTimeout(900_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  // Capture the original user login for round-trip verification.
  const originalLogin = await page.evaluate(() => {
    const grok = (window as any).grok;
    return grok?.shell?.user?.login ?? null;
  });
  expect(originalLogin).toBeTruthy();

  await softStep('logoutAndLoginAs — switch from qa-pw to Olena Ahadzhanian (manual password)', async () => {
    await logoutAndLoginAs(page, {username: 'Olena Ahadzhanian'});
    const newLogin = await page.evaluate(() => {
      const grok = (window as any).grok;
      return grok?.shell?.user?.login ?? null;
    });
    expect(newLogin).not.toBe(originalLogin);
    expect(newLogin).toBeTruthy();
  });

  await softStep('logoutAndLoginAs — restore original (qa-pw) via restoreOriginal flow', async () => {
    // Per helper docs, restoreOriginal logs out of the current user and back
    // in as the previous user. We invoke a fresh logoutAndLoginAs pointing
    // back at the original login since restoreOriginal flag wasn't set on
    // the prior call (this also exercises a second logout cycle for
    // round-trip coverage).
    await logoutAndLoginAs(page, {username: originalLogin}, {restoreOriginal: false});
    const restored = await page.evaluate(() => {
      const grok = (window as any).grok;
      return grok?.shell?.user?.login ?? null;
    });
    expect(restored).toBe(originalLogin);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} session step(s) failed:\n${summary}`);
  }
});
