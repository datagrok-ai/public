import {test, Page} from '@playwright/test';

export const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
export const login = process.env.DATAGROK_LOGIN ?? 'admin';
export const password = process.env.DATAGROK_PASSWORD ?? 'admin';

export const specTestOptions = {
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
};

export interface StepError { step: string; error: string; }

export const stepErrors: StepError[] = [];

export async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) {
    stepErrors.push({step: name, error: e?.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e?.message ?? e}`);
  }
}

export async function loginToDatagrok(page: Page) {
  await page.goto(baseUrl);
  // Fresh-context flake: the login form sometimes appears while Dart is still
  // mounting listeners — typing + Enter are then silently dropped. Give the
  // page a moment to settle before probing the form. `networkidle` is the
  // cheapest reliable signal that Dart boot is done.
  await page.waitForLoadState('networkidle', {timeout: 30000}).catch(() => {});

  const submitLogin = async () => {
    // The active login inputs live inside `#signup-login-fields` (there are
    // several hidden duplicates for Reset / Sign-Up flows elsewhere on the
    // page). Scope directly to that container — and use `.focus()` instead of
    // `.click()` to avoid the "signup-container intercepts pointer events"
    // flake, where Playwright's actionability check picks the wrong stacking
    // context on fresh page loads.
    const loginInput = page.locator('#signup-login-fields input[placeholder="Login or Email"]');
    const passwordInput = page.locator('#signup-login-fields input[placeholder="Password"]');
    if (!(await loginInput.isVisible({timeout: 15000}).catch(() => false)))
      return;
    await loginInput.focus();
    await page.keyboard.type(login);
    // Dart re-renders #signup-login-fields on input — pre-resolved locators
    // detach. Tab to the next field so we don't need to re-locate the
    // password input through Playwright's actionability pipeline.
    await page.keyboard.press('Tab');
    await page.keyboard.type(password);
    await page.keyboard.press('Enter');
    // Dev/public login markup has no <form>, so Enter is sometimes a no-op —
    // as a backup, try to click the "Login" button (filtered to exclude "Login
    // with Google"). Race it against the post-login Browse tab so we don't
    // block if login already succeeded or the button is mid-detach.
    const loginBtn = page.locator('button').filter({hasText: /^Login$/}).and(page.locator(':visible'));
    await Promise.race([
      page.locator('[name="Browse"]').waitFor({timeout: 5000}).catch(() => {}),
      (async () => {
        if (await loginBtn.isVisible({timeout: 2000}).catch(() => false))
          await loginBtn.click({force: true, timeout: 3000}).catch(() => {});
      })(),
    ]);
  };

  await submitLogin();
  // One retry if Browse still hasn't appeared — covers cases where the first
  // submit was dropped because Dart hadn't wired up the change listener yet.
  if (!(await page.locator('[name="Browse"]').isVisible({timeout: 15000}).catch(() => false)))
    await submitLogin();

  await page.locator('[name="Browse"]').waitFor({timeout: 120000});
}
