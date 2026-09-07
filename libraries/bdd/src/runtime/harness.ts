/* One browser page for all the scenarios of a feature. The first scenario opens it (inside a test,
   so Playwright applies the project's context options — login storage state, viewport — and
   records traces and screenshots per scenario as usual), every scenario ends with the shell reset
   (dialogs and popups gone, `grok.shell.closeAll()`, back on the Home view), the last one closes
   it. The generated spec keeps calling `test()` itself, so reports point at the spec line.
   The page's console errors and uncaught exceptions are collected from the moment it opens, so a
   scenario can assert a zero-error floor (`no errors should have been logged`). */
import type {Browser, Page, PlaywrightTestArgs, PlaywrightTestOptions, PlaywrightWorkerArgs, PlaywrightWorkerOptions,
  TestType} from '@playwright/test';
import {expect} from '@playwright/test';
import {leave} from './args.js';

type Test = TestType<PlaywrightTestArgs & PlaywrightTestOptions, PlaywrightWorkerArgs & PlaywrightWorkerOptions>;

export interface FeatureSession {
  /** The feature's page — opened on first use, shared by the scenarios that follow. */
  page(browser: Browser): Promise<Page>;
}

export interface Journey {
  /** Runs a scenario as a soft step: a failure is recorded and the next scenario still runs. */
  scenario(name: string, body: () => Promise<void>): Promise<void>;
  /** Fails the test when any scenario failed, listing them. */
  finish(): void;
}

/** A `@journey` feature: one test, the Background once, the scenarios in order on the same shell
 * state — the way a hand-written spec chains `softStep`s. Each scenario leaves the state it changed
 * as it found it, so the next one starts where the Background left off. The test's budget is the
 * per-test timeout times the scenario count. */
export function journey(test: Test, scenarios: number): Journey {
  test.setTimeout(test.info().timeout * scenarios);
  const failed: string[] = [];
  return {
    async scenario(name: string, body: () => Promise<void>): Promise<void> {
      try {
        await test.step(name, body);
      }
      catch (e) {
        failed.push(`${name}: ${(e as Error).message ?? e}`);
      }
    },
    finish(): void {
      expect(failed, `${failed.length} of ${scenarios} scenarios failed`).toEqual([]);
    },
  };
}

const HOME_VIEW = 'datagrok';
const LEFTOVERS = '[data-u2="dialog"], [data-u2="menu"], [data-u2="tooltip"], [data-u2="notify"] > *, ' +
  '.d4-dialog, .d4-menu-popup, .d4-balloon, .d4-tooltip';

const errors = new WeakMap<Page, string[]>();

/** Starts collecting the page's console errors and uncaught exceptions. */
export function watchErrors(page: Page): void {
  if (errors.has(page))
    return;
  const list: string[] = [];
  errors.set(page, list);
  page.on('console', (m) => {
    // a resource the stand does not serve (a help page) is logged as a console error by the
    // browser, not raised by the platform's code — not part of the error floor
    if (m.type() === 'error' && !m.text().startsWith('Failed to load resource'))
      list.push(m.location().url ? `${m.text()} (${m.location().url})` : m.text());
  });
  page.on('pageerror', (e) => list.push(String(e)));
}

/** The errors logged since the last call (or since the page opened), and clears them. */
export function takeErrors(page: Page): string[] {
  const list = errors.get(page) ?? [];
  const out = [...list];
  list.length = 0;
  return out;
}

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
      if (!page || page.isClosed()) {
        page = await browser.newPage();
        watchErrors(page);
      }
      return page;
    },
  };
}

/** Everything closed and the Home view current — the state the next scenario starts from. A page
 * that is not in the shell (about:blank, the login page) is left alone. Errors the teardown itself
 * raises (work cancelled by `closeAll`) are dropped, so they are not charged to the next scenario. */
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
  takeErrors(page);
}
