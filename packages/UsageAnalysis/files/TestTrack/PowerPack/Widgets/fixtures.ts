import { test as base, expect, Page, BrowserContext } from '@playwright/test';

type WorkerFixtures = { _app: { context: BrowserContext; page: Page } };
type TestFixtures = { homePage: Page };

export const test = base.extend<TestFixtures, WorkerFixtures>({
  _app: [async ({ browser }, use, workerInfo) => {
    const storageState = (workerInfo.project.use as { storageState?: string }).storageState;
    const context = await browser.newContext({ storageState });
    const page = await context.newPage();
    await use({ context, page });
    await context.close();
  }, { scope: 'worker' }],

  homePage: async ({ _app }, use) => {
    await use(_app.page);
  },
});

export { expect };
