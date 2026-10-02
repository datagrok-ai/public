/* One browser page for all the scenarios of a feature. The first scenario opens it (inside a test,
   so Playwright applies the project's context options — login storage state, viewport — and
   records traces and screenshots per scenario as usual), every scenario ends with the shell reset
   (dialogs and popups gone, `grok.shell.closeAll()`, back on the Home view), the last one closes
   it. The generated spec keeps calling `test()` itself, so reports point at the spec line.
   The page's console errors and uncaught exceptions are collected from the moment it opens, so a
   scenario can assert a zero-error floor (`no errors should have been logged`). */
import {join, sep} from 'node:path';
import {fileURLToPath} from 'node:url';
import {randomUUID} from 'node:crypto';
import type {Browser, Page, PlaywrightTestArgs, PlaywrightTestOptions, PlaywrightWorkerArgs, PlaywrightWorkerOptions,
  TestType} from '@playwright/test';
import {leave} from './args.js';
import {failure, isSkip, isWaitFailure, journeyFailure, reasonOf} from './failure.js';
import * as guide from './guide.js';
import {explain} from './locate.js';
import {logMemory, rendererMb} from './memory.js';
import {whileExpectedToFail} from './patience.js';
import {takeBalloons} from './viewers.js';

type Test = TestType<PlaywrightTestArgs & PlaywrightTestOptions, PlaywrightWorkerArgs & PlaywrightWorkerOptions>;

export interface FeatureSession {
  /** Resolves {run} to this feature instance's unique suffix and {time} to its start in epoch ms. */
  text(value: string): string;
  /** The feature's page — opened on first use, shared by the scenarios that follow. */
  page(browser: Browser): Promise<Page>;
  /** One Gherkin step: a Playwright step located at the feature line, whose failure names the
   * line, the step as written, the reason, and — when Playwright gave up on an element — what the
   * page shows where the phrase looked. */
  step(line: number, title: string, body: () => Promise<unknown>, table?: string[][]): Promise<void>;
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
      const overlays = page ? await openOverlays(page) : 0;
      try {
        if (page) {
          takeErrors(page);
          await takeBalloons(page).catch(() => undefined);
        }
        await test.step(name, options?.knownFailure ? () => whileExpectedToFail(body) : body);
      }
      catch (e) {
        // a capability gate's skip ends the journey as skipped, unless an earlier scenario failed: a skip
        // must not hide that failure
        if (isSkip(e))
          throw failed.length > 0 ? journeyFailure(failed, scenarios, `the rest was skipped at "${name}": ${reasonOf(e)}`) : e;
        if (!options?.knownFailure)
          failed.push({name, error: e});
        // a scenario that stopped midway leaves its dialog or menu over the viewers the next one uses;
        // those open before it (a dialog the Background opened for every scenario) stay
        if (page)
          await closeOverlays(page, overlays);
        return;
      }
      if (options?.knownFailure)
        failed.push({name, error: new Error(KNOWN_FAILURE_PASSED)});
    },
    finish(): void {
      if (failed.length > 0)
        throw journeyFailure(failed, scenarios);
    },
  };
}

const KNOWN_FAILURE_PASSED = 'tagged @known-failure and passed — the bug it describes is fixed, so the tag has to go';

/** A `@known-failure` scenario outside a journey: its steps failing is the defect it describes and
 * passes the test; its steps passing fails it — the bug is fixed and the tag has to go. */
export async function knownFailure(body: () => Promise<void>): Promise<void> {
  try {
    await whileExpectedToFail(body);
  }
  catch (e) {
    if (isSkip(e))
      throw e;
    return;
  }
  throw new Error(KNOWN_FAILURE_PASSED);
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
const NOTICES = '[data-u2="notify"] > *, .d4-balloon, .ui-hint-popup';

const errors = new WeakMap<Page, string[]>();
const cleanups = new WeakMap<Page, (() => Promise<void>)[]>();
// a page whose renderer crashed stays open, and every call on it fails
const crashed = new WeakSet<Page>();
const usable = (page: Page): boolean => !page.isClosed() && !crashed.has(page);

/** Closes a crashed page and drops its feature-end cleanups, which call into it; how many it dropped. */
async function abandon(page: Page): Promise<number> {
  const dropped = cleanups.get(page)?.length ?? 0;
  cleanups.delete(page);
  await page.close().catch(() => undefined);
  return dropped;
}

/** Runs when the feature ends, whatever its scenarios did — for state a step created on the server
 * (a project, an uploaded table); a page that crashed runs none, and the feature fails for them. */
export function atFeatureEnd(page: Page, cleanup: () => Promise<void>, first = false): void {
  const list = cleanups.get(page) ?? [];
  cleanups.set(page, list);
  if (first)
    list.unshift(cleanup);
  else
    list.push(cleanup);
}

/** The two console errors the browser raises about something that is not the platform's code.
 * Both are matched on the message AND on where it came from — a broad pattern here is how a
 * suite ends up silencing the failures it exists to catch. */
function ignoredError(text: string, url: string): boolean {
  // a resource the stand does not serve (a help page), logged by the browser rather than raised
  if (text.startsWith('Failed to load resource'))
    return true;
  // an embedded third-party player refusing a feature policy of the page it is framed in:
  // a card of the Projects gallery carries a YouTube iframe, and its complaint is not ours
  return text.startsWith('Permissions policy violation') && /^https:\/\/(www\.)?youtube\.com\//.test(url);
}

/** Starts collecting the page's console errors and uncaught exceptions. */
export function watchErrors(page: Page): void {
  if (errors.has(page))
    return;
  const list: string[] = [];
  errors.set(page, list);
  // "Stack trace X" arrives seconds after its "Look below, ID = X" error: joined while unreported, else dropped
  const announced = new Set<string>();
  page.on('console', (m) => {
    const text = m.text();
    if (m.type() !== 'error' || ignoredError(text, m.location().url))
      return;
    const continuation = /^Stack trace (\S+)/.exec(text);
    if (continuation && announced.has(continuation[1])) {
      const parent = list.findIndex((e) => new RegExp(`Look below, ID = ${continuation[1]}(\\s|$)`).test(e));
      if (parent >= 0)
        list[parent] += `\n${text}`;
      return;
    }
    const id = /Look below, ID = (\S+)/.exec(text);
    if (id)
      announced.add(id[1]);
    list.push(m.location().url ? `${text} (${m.location().url})` : text);
  });
  page.on('pageerror', (e) => list.push(String(e)));
  page.on('crash', () => crashed.add(page));
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
 * logged in` on the shell it finds, reset). A page that closed or crashed (a failed test restarting
 * the worker, a page past its feature count below) is replaced in the same context, which keeps the
 * storage state and the HTTP cache. The browser fixture closes the context with the worker. */
let shared: Page | undefined;
let lastTime = 0;

/** A page is replaced after `BDD_PAGE_MAX_FEATURES` features or once its renderer holds more than
 * `BDD_PAGE_MAX_MB`: what it keeps meanwhile and what a new one costs are in CLAUDE.md. */
const PAGE_MAX_FEATURES = Number(process.env.BDD_PAGE_MAX_FEATURES ?? 25);
const PAGE_MAX_MB = Number(process.env.BDD_PAGE_MAX_MB ?? 3000);
// the operating system is asked for the renderer's memory, ~0.2 s a reading: every third feature
const PAGE_MB_EVERY = 3;
const featuresRun = new WeakMap<Page, number>();
let mbUnreadable = false;

async function afterFeature(page: Page, path: string): Promise<void> {
  const count = (featuresRun.get(page) ?? 0) + 1;
  featuresRun.set(page, count);
  const logged = await logMemory(page, {feature: path, onPage: count});
  let why = PAGE_MAX_FEATURES > 0 && count >= PAGE_MAX_FEATURES ? `${count} features` : '';
  if (!why && PAGE_MAX_MB > 0 && (logged !== undefined || count % PAGE_MB_EVERY === 0)) {
    const mb = logged ?? await rendererMb(page).catch((error) => {
      if (!mbUnreadable)
        console.warn(`bdd: the renderer's memory cannot be read, so BDD_PAGE_MAX_MB does not apply: ${error}`);
      mbUnreadable = true;
      return 0;
    });
    why = mb > PAGE_MAX_MB ? `${mb} MB in its renderer` : '';
  }
  if (!why)
    return;
  console.warn(`bdd: a new page after ${path}: this one holds ${why}`);
  await page.close().catch(() => undefined);
}

export function feature(test: Test, path = '', specUrl = ''): FeatureSession {
  let page: Page | undefined;
  // the cleanups of a page that crashed during the feature: none of them ran
  let lost = 0;
  const runId = randomUUID();
  // a login takes only [a-z0-9._-], and a user can never be deleted, so a fixture user is named by
  // when it was made; two features of one worker never start in the same millisecond
  const time = String(lastTime = Math.max(Date.now(), lastTime + 1));
  const text = (value: string): string => value.replaceAll('{run}', runId).replaceAll('{time}', time);
  const file = path && specUrl ? featureFile(specUrl, path) : undefined;
  test.afterEach(async () => {
    if (page && usable(page)) {
      leave(page);
      await resetShell(page);
    }
  });
  test.afterAll(async () => {
    const failures: unknown[] = [];
    if (page && crashed.has(page))
      lost += await abandon(page);
    else if (page && !page.isClosed()) {
      const list = cleanups.get(page) ?? [];
      cleanups.delete(page);
      // a cleanup deletes what is still there, so one that failed on a request the stand dropped under load (nginx
      // answering 502 when its connection to Datlas fails) is run again before the feature fails for it
      for (const [i, cleanup] of list.entries()) {
        if (crashed.has(page)) {
          lost += list.length - i;
          break;
        }
        for (let attempt = 1; ; attempt++) {
          try {
            await cleanup();
            break;
          }
          catch (error) {
            if (attempt === 3 || crashed.has(page)) {
              failures.push(error);
              break;
            }
            await page.waitForTimeout(2000).catch(() => undefined);
          }
        }
      }
      if (crashed.has(page))
        await page.close().catch(() => undefined);
      else
        await afterFeature(page, path);
    }
    page = undefined;
    if (lost) {
      failures.push(new Error(`the page crashed: ${lost} feature-end cleanup(s) could not run, ` +
        'what they remove stays on the server'));
    }
    if (failures.length)
      throw new AggregateError(failures, 'Feature cleanup failed');
  });
  return {
    text,
    async page(browser: Browser): Promise<Page> {
      if (!page || !usable(page)) {
        if (page && crashed.has(page))
          lost += await abandon(page);
        // a crashed page's context may be signed in as another account (a feature's sign-in, whose
        // switch back went with the dropped cleanups): a new context starts from the configured state
        if (shared && (shared.context().browser() !== browser || crashed.has(shared))) {
          await shared.context().close().catch(() => undefined);
          shared = undefined;
        }
        if (shared?.isClosed())
          shared = await shared.context().newPage();
        shared ??= await (await browser.newContext()).newPage();
        watchErrors(shared);
        await guide.attach(shared);
        page = shared;
      }
      return page;
    },
    async step(line: number, title: string, body: () => Promise<unknown>, table?: string[][]): Promise<void> {
      title = text(title);
      await test.step(title, async () => {
        await guide.begin(page, test.info(), line, title, table);
        try {
          await body();
        }
        catch (e) {
          if (isSkip(e))
            throw e;
          const shown = isWaitFailure(e) && page && usable(page) ? await explain(page).catch(() => '') : '';
          throw failure(`${path || 'feature'}:${line}`, title, e, shown, file ? `${file}:${line}:1` : '');
        }
        finally {
          await guide.end(page);
        }
      }, {location: file ? {file, line, column: 1} : undefined});
    },
  };
}

/** Task bar entries a reset already waited out: a job that never ends is waited for once, not by
 * every later scenario of the page. */
const stuckEntries = new WeakMap<Page, Set<string>>();
const SETTLE_MS = 25000;
const COMMAND_SETTLE_MS = 60000;

/** Work a scenario started and never awaited (a known failure ends at its first failing claim) must
 * not finish on the next feature's shell: an analysis that ends reopens its table and makes it
 * current. The command the scenario armed says exactly when its call is over and gets the longer
 * wait (an embedding under another worker's load takes over 25 s); the platform's progress entries
 * cover work a command's onAfterRunAction comes before, on their own budget. */
async function settleWork(page: Page): Promise<void> {
  const settled = await page.evaluate((ms) => (window as any).__bdd?.settleCommand?.(ms) ?? true, COMMAND_SETTLE_MS).catch(() => true);
  if (!settled)
    console.warn(`bdd: a menu command the scenario started was still running ${COMMAND_SETTLE_MS / 1000} s into the shell reset`);
  const deadline = Date.now() + SETTLE_MS;
  const stuck = stuckEntries.get(page) ?? new Set<string>();
  stuckEntries.set(page, stuck);
  const running = (known: string[]): string[] => Array.from(document.querySelectorAll('.d4-task-bar-entry'))
    .filter((e) => (e as HTMLElement).offsetParent !== null)
    .map((e) => (e.textContent ?? '').trim())
    .filter((t) => !known.includes(t));
  const known = [...stuck];
  // the predicate runs in the page, where only its own source exists: `running` goes in as text
  const quiet = await page.waitForFunction(`(${running})(${JSON.stringify(known)}).length === 0`, undefined,
    {timeout: Math.max(1, deadline - Date.now()), polling: 100}).then(() => true).catch(() => false);
  const left: string[] = quiet ? [] : await page.evaluate(running, known).catch(() => []);
  if (left.length > 0) {
    for (const t of left)
      stuck.add(t);
    console.warn(`bdd: still running ${SETTLE_MS / 1000} s into the shell reset (not waited for again): ${left.join(' | ')}`);
  }
}

/** Dialogs and menus closed the platform's way: Escape, as many times as there are open ones. */
const openOverlays = (page: Page): Promise<number> =>
  page.locator(CLOSABLE).filter({visible: true}).count().catch(() => 0);

/** Escape, the topmost first, while more dialogs and menus are open than `keep`. */
async function closeOverlays(page: Page, keep = 0): Promise<void> {
  for (let i = 0; i < 3 && await openOverlays(page) > keep; i++)
    await page.keyboard.press('Escape').catch(() => undefined);
}

/** Everything closed and the Home view current — the state the next scenario starts from. A page
 * that is not in the shell (about:blank, the login page) is left alone. Work the scenario left
 * running is waited out first; then dialogs and menus are closed the platform's way (Escape, as
 * many times as there are open ones), the tooltip through its API, notifications as their close
 * icons would; a dialog that survives that is reported in the run's output rather than pulled out
 * of the DOM behind the platform's back. Errors the teardown itself raises (work cancelled by
 * `closeAll`) are dropped. */
export async function resetShell(page: Page): Promise<void> {
  const inShell = await page.evaluate(() => typeof (window as any).grok?.shell?.closeAll === 'function').catch(() => false);
  if (!inShell)
    return;
  await settleWork(page);
  await closeOverlays(page);
  const left: string = await page.evaluate((notices) => {
    const w = window as any;
    w.ui?.tooltip?.hide?.();
    for (const e of document.querySelectorAll(notices))
      e.remove();
    // a demo's script panel is docked, not a view, and its script keeps going: its own Back button closes and cancels it
    for (const script of document.querySelectorAll('.demo-app-script'))
      (script.querySelector('.tutorials-root-header > button') as HTMLElement | null)?.click();
    if (w.grok.shell.windows.presentationMode)
      w.grok.shell.windows.presentationMode = false;
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
