import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {logoutAndLoginAs} from './session';

test.use({...specTestOptions, navigationTimeout: 180_000});

const readLogin = (page: any) =>
  page.evaluate(() => (window as any).grok?.shell?.user?.login ?? null);

test('helpers.playwright.session — logoutAndLoginAs end-to-end', async ({page}) => {

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
