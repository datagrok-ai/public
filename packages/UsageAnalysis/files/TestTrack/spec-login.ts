import {test, Page} from '@playwright/test';

export const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';

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

async function injectToken(page: Page, token: string) {
  // Navigate to the origin first so the cookie/localStorage entries are
  // attached to the right host. The `/oauth/` path matches what `grok test`
  // does for Puppeteer (test-utils.ts:135).
  await page.goto(baseUrl + '/oauth/');
  const u = new URL(baseUrl);
  await page.context().addCookies([{name: 'auth', value: token, domain: u.hostname, path: '/'}]);
  await page.evaluate((t) => window.localStorage.setItem('auth', t), token);
  await page.goto(baseUrl);
  await page.waitForFunction(() => document.querySelector('.grok-preloader') == null, null, {timeout: 120_000});
  await page.locator('[name="Browse"]').waitFor({timeout: 60_000});
}

export async function loginToDatagrok(page: Page) {
  const token = process.env.DATAGROK_AUTH_TOKEN;
  if (!token || token.length === 0)
    throw new Error('DATAGROK_AUTH_TOKEN is not set. Run via `grok test`, which derives the token from ~/.grok/config.yaml.');
  await injectToken(page, token);
}

export async function loginAsSecondUser(page: Page) {
  const token2 = process.env.DATAGROK_AUTH_TOKEN_2;
  if (!token2 || token2.length === 0)
    throw new Error('DATAGROK_AUTH_TOKEN_2 is not set. Provide a second-user dev key via DATAGROK_DEV_KEY_2 (the runner exchanges it for a token).');
  await injectToken(page, token2);
}
