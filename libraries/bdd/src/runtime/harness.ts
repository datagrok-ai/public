/* One browser page for all the scenarios of a feature. The first scenario opens it (inside a test,
   so Playwright applies the project's context options — login storage state, viewport — and
   records traces and screenshots per scenario as usual), every scenario ends with the shell reset
   (dialogs and popups gone, `grok.shell.closeAll()`, back on the Home view), the last one closes
   it. The generated spec keeps calling `test()` itself, so reports point at the spec line.
   The page's console errors and uncaught exceptions are collected from the moment it opens, so a
   scenario can assert a zero-error floor (`no errors should have been logged`). */
import {join, sep} from 'node:path';
import {fileURLToPath} from 'node:url';
import type {Browser, Page, PlaywrightTestArgs, PlaywrightTestOptions, PlaywrightWorkerArgs, PlaywrightWorkerOptions,
  TestType} from '@playwright/test';
import {leave} from './args.js';
import {failure, isWaitFailure, journeyFailure} from './failure.js';
import {explain} from './locate.js';
import {whileExpectedToFail} from './patience.js';
import {takeBalloons} from './viewers.js';

type Test = TestType<PlaywrightTestArgs & PlaywrightTestOptions, PlaywrightWorkerArgs & PlaywrightWorkerOptions>;

export interface FeatureSession {
  /** The feature's page — opened on first use, shared by the scenarios that follow. */
  page(browser: Browser): Promise<Page>;
  /** One Gherkin step: a Playwright step located at the feature line, whose failure names the
   * line, the step as written, the reason, and — when Playwright gave up on an element — what the
   * page shows where the phrase looked. */
  step(line: number, title: string, body: () => Promise<unknown>): Promise<void>;
}

export interface Journey {
  /** Runs a scenario as a soft step: a failure is recorded and the next scenario still runs.
   * `knownFailure` inverts that scenario: its failure is expected and does not fail the test,
   * while its passing does — the bug it describes is fixed and the tag has to go. */
  scenario(name: string, body: () => Promise<void>, options?: {knownFailure?: boolean}): Promise<void>;
  /** Fails the test when any scenario failed, listing them. */
  finish(): void;
}

/** What a journey scenario may add to the test's budget: a scenario takes a few seconds, and a
 * hung one should not hold a worker for the whole per-test timeout times the scenario count. */
const SCENARIO_BUDGET_MS = 20000;

/** A `@journey` feature: one test, the Background once, the scenarios in order on the same shell
 * state — the way a hand-written spec chains `softStep`s. Each scenario leaves the state it changed
 * as it found it, so the next one starts where the Background left off, and owns its error and
 * balloon floors: what an earlier scenario logged is not charged to it. */
export function journey(test: Test, scenarios: number, page?: Page): Journey {
  test.setTimeout(test.info().timeout + SCENARIO_BUDGET_MS * scenarios);
  const failed: {name: string; error: unknown}[] = [];
  return {
    async scenario(name: string, body: () => Promise<void>, options?: {knownFailure?: boolean}): Promise<void> {
      try {
        if (page) {
          takeErrors(page);
          await takeBalloons(page).catch(() => undefined);
        }
        await test.step(name, options?.knownFailure ? () => whileExpectedToFail(body) : body);
      }
      catch (e) {
        if (!options?.knownFailure)
          failed.push({name, error: e});
        return;
      }
      if (options?.knownFailure)
        failed.push({name, error: new Error('tagged @known-failure and passed — the bug it describes is fixed, so the tag has to go')});
    },
    finish(): void {
      if (failed.length > 0)
        throw journeyFailure(failed, scenarios);
    },
  };
}

/** `<root>/generated/x/y.test.ts` + `features/x/y.feature` → the feature file (the layout
 * `outFileFor` writes). */
function featureFile(specUrl: string, path: string): string {
  const spec = fileURLToPath(specUrl);
  const i = spec.lastIndexOf(`${sep}generated${sep}`);
  return join(i < 0 ? spec : spec.slice(0, i), ...path.split('/'));
}

const HOME_VIEW = 'datagrok';
// what Escape closes: dialogs and popup menus of both UI generations
const CLOSABLE = '[data-u2="dialog"], [data-u2="menu"], .d4-dialog, .d4-menu-popup';
// transient notifications, taken away as their close icons would
const NOTICES = '[data-u2="notify"] > *, .d4-balloon';

const errors = new WeakMap<Page, string[]>();
const cleanups = new WeakMap<Page, (() => Promise<void>)[]>();

/** Runs when the feature's page closes, whatever its scenarios did — for state a step created on
 * the server (a project, an uploaded table). */
export function atFeatureEnd(page: Page, cleanup: () => Promise<void>): void {
  const list = cleanups.get(page) ?? [];
  cleanups.set(page, list);
  list.push(cleanup);
}

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

/** One page per worker: every feature the worker runs uses the page the first one opened (the
 * shell boots once, ~4 s, and a package initializes once; the next feature starts from `user is
 * logged in` on the shell it finds, reset). A page that closed (a crash, a failed test restarting
 * the worker) is replaced in the same context, which keeps the storage state and the HTTP cache.
 * The browser fixture closes the context with the worker. */
let shared: Page | undefined;

export function feature(test: Test, path = '', specUrl = ''): FeatureSession {
  let page: Page | undefined;
  const file = path && specUrl ? featureFile(specUrl, path) : undefined;
  test.afterEach(async () => {
    if (page && !page.isClosed()) {
      leave(page);
      await resetShell(page);
    }
  });
  test.afterAll(async () => {
    if (page && !page.isClosed()) {
      for (const cleanup of cleanups.get(page) ?? [])
        await cleanup().catch((e) => console.warn(`cleanup failed: ${(e as Error).message}`));
      cleanups.delete(page);
    }
    page = undefined;
  });
  return {
    async page(browser: Browser): Promise<Page> {
      if (!page || page.isClosed()) {
        if (shared && shared.context().browser() !== browser) {
          await shared.context().close().catch(() => undefined);
          shared = undefined;
        }
        if (shared && shared.isClosed())
          shared = await shared.context().newPage();
        shared ??= await (await browser.newContext()).newPage();
        watchErrors(shared);
        page = shared;
      }
      return page;
    },
    async step(line: number, title: string, body: () => Promise<unknown>): Promise<void> {
      await test.step(title, async () => {
        try {
          await body();
        }
        catch (e) {
          const shown = isWaitFailure(e) && page && !page.isClosed() ? await explain(page).catch(() => '') : '';
          throw failure(`${path || 'feature'}:${line}`, title, e, shown, file ? `${file}:${line}:1` : '');
        }
      }, {location: file ? {file, line, column: 1} : undefined});
    },
  };
}

/** Everything closed and the Home view current — the state the next scenario starts from. A page
 * that is not in the shell (about:blank, the login page) is left alone. Dialogs and menus are
 * closed the platform's way (Escape, as many times as there are open ones), the tooltip through
 * its API, notifications as their close icons would; a dialog that survives that is reported in
 * the run's output rather than pulled out of the DOM behind the platform's back. Errors the
 * teardown itself raises (work cancelled by `closeAll`) are dropped. */
export async function resetShell(page: Page): Promise<void> {
  const inShell = await page.evaluate(() => typeof (window as any).grok?.shell?.closeAll === 'function').catch(() => false);
  if (!inShell)
    return;
  const open = (): Promise<number> => page.locator(CLOSABLE).filter({visible: true}).count().catch(() => 0);
  for (let i = 0; i < 3 && await open() > 0; i++)
    await page.keyboard.press('Escape').catch(() => undefined);
  const left: string = await page.evaluate((notices) => {
    const w = window as any;
    w.ui?.tooltip?.hide?.();
    for (const e of document.querySelectorAll(notices))
      e.remove();
    w.grok.shell.closeAll();
    return Array.from(document.querySelectorAll('.d4-dialog, [data-u2="dialog"]'))
      .filter((e) => (e as HTMLElement).offsetParent !== null).map((e) => e.getAttribute('name') ?? e.tagName).join(', ');
  }, NOTICES).catch(() => '');
  if (left)
    console.warn(`bdd: a dialog is still open after the shell reset: ${left}`);
  // a view that closeAll leaves (the Model Hub's card view survives it) would otherwise keep the
  // next scenario off the Home view: closed one by one, then reported rather than waited for
  const stayed: string = await page.waitForFunction((home) => (window as any).grok?.shell?.v?.type === home, HOME_VIEW, {timeout: 5000})
    .then(() => '')
    .catch(() => page.evaluate((home) => {
      const w = window as any;
      const others = Array.from(w.grok.shell.views as Iterable<any>).filter((v) => v.type !== home);
      for (const v of others)
        v.close();
      const still = Array.from(w.grok.shell.views as Iterable<any>).filter((v) => v.type !== home);
      const homeView = Array.from(w.grok.shell.views as Iterable<any>).find((v) => v.type === home);
      if (homeView)
        w.grok.shell.v = homeView;
      return still.map((v) => `${v.type}:${v.name}`).join(', ');
    }, HOME_VIEW).catch(() => ''));
  if (stayed)
    console.warn(`bdd: views that survived closeAll and their own close() at the shell reset: ${stayed}`);
  await page.waitForFunction((home) => (window as any).grok?.shell?.v?.type === home, HOME_VIEW, {timeout: 60000})
    .catch(() => console.warn('bdd: the Home view is not current after the shell reset'));
  takeErrors(page);
}
