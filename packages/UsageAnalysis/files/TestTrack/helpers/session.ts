import {Page, expect} from '@playwright/test';
import {loginToDatagrok, loginAsSecondUser} from '../spec-login';

const currentLogin = (page: Page): Promise<string | null> =>
  page.evaluate(() => (window as any).grok?.shell?.user?.login ?? null);

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
