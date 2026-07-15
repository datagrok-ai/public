import { test as base, expect, Page, BrowserContext } from '@playwright/test';

// A worker-scoped page shared by every test in the worker. The app shell is booted once
// (an expensive, occasionally-flaky operation on dev); between tests we reset in-app with
// `grok.shell.closeAll()` instead of reloading — see resetHome() in helpers.ts. This turns
// ~N full SPA boots per file into a single one.
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
