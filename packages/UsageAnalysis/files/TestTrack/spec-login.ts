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
  const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
  if (await loginInput.isVisible({timeout: 15000}).catch(() => false)) {
    await loginInput.click();
    await page.keyboard.type(login);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).click();
    await page.keyboard.type(password);
    await page.keyboard.press('Enter');
  }
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});
}
