/* One browser page for all the scenarios of a feature. The first scenario opens it (inside a test,
   so Playwright applies the project's context options — login storage state, viewport — and
   records traces and screenshots per scenario as usual), every scenario ends with the shell reset
   (dialogs and popups gone, `grok.shell.closeAll()`, back on the Home view), the last one closes
   it. The generated spec keeps calling `test()` itself, so reports point at the spec line. */
import type {Browser, Page, PlaywrightTestArgs, PlaywrightTestOptions, PlaywrightWorkerArgs, PlaywrightWorkerOptions,
  TestType} from '@playwright/test';
import {leave} from './args.js';

type Test = TestType<PlaywrightTestArgs & PlaywrightTestOptions, PlaywrightWorkerArgs & PlaywrightWorkerOptions>;

export interface FeatureSession {
  /** The feature's page — opened on first use, shared by the scenarios that follow. */
  page(browser: Browser): Promise<Page>;
}

const HOME_VIEW = 'datagrok';
const LEFTOVERS = '[data-u2="dialog"], [data-u2="menu"], [data-u2="tooltip"], [data-u2="notify"] > *, ' +
  '.d4-dialog, .d4-menu-popup, .d4-balloon, .d4-tooltip';

export function feature(test: Test): FeatureSession {
  let page: Page | undefined;
  test.afterEach(async () => {
    if (page && !page.isClosed()) {
      leave(page);
      await resetShell(page);
    }
  });
  test.afterAll(async () => {
    await page?.context().close().catch(() => undefined);
    page = undefined;
  });
  return {
    async page(browser: Browser): Promise<Page> {
      if (!page || page.isClosed())
        page = await browser.newPage();
      return page;
    },
  };
}

/** Everything closed and the Home view current — the state the next scenario starts from. A page
 * that is not in the shell (about:blank, the login page) is left alone. */
export async function resetShell(page: Page): Promise<void> {
  const inShell = await page.evaluate(() => typeof (window as any).grok?.shell?.closeAll === 'function').catch(() => false);
  if (!inShell)
    return;
  await page.keyboard.press('Escape').catch(() => undefined);
  await page.evaluate((leftovers) => {
    (window as any).grok.shell.closeAll();
    for (const e of document.querySelectorAll(leftovers))
      e.remove();
  }, LEFTOVERS).catch(() => undefined);
  await page.waitForFunction((home) => (window as any).grok?.shell?.v?.type === home, HOME_VIEW, {timeout: 60000})
    .catch(() => undefined);
}
