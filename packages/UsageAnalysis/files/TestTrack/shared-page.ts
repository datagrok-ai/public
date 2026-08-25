import {test as base, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, stepErrors} from './spec-login';

/**
 * A Datagrok that is booted once per worker instead of once per spec.
 *
 * The platform costs ~10s to come up (3.5s for the document, 5.8s registering ~3300
 * functions) and every spec was paying it, which is more wall-clock than the assertions
 * themselves. The page below is created and logged into once, then handed to every test
 * in the worker; `resetShell` puts it back to a first-run state afterwards.
 *
 * What this trades away is process-level isolation: a test that corrupts the shell in a
 * way `resetShell` does not undo will be seen by the next one. Anything a test creates on
 * the SERVER (projects, layouts) is still its own to delete — that has always been true.
 */

const OVERLAYS = ['.d4-menu-popup', '.d4-column-selector-backdrop', '.d4-tooltip', '.d4-dialog'];

export async function resetShell(page: Page): Promise<void> {
  // Dismiss the platform's way first, and only rip nodes out if that fails: an overlay
  // removed from the DOM behind the platform's back leaves it believing a popup is still
  // open, and the next test's clicks land on a modal layer that is no longer visible.
  // These were the shared-page failures — "popup did not open" and click timeouts, in
  // specs that pass in isolation.
  await page.keyboard.press('Escape').catch(() => {});
  await page.keyboard.press('Escape').catch(() => {});
  // hover state is per-page too: the on-canvas column selectors only appear under the
  // pointer, so a pointer parked on the last test's viewer changes what the next one sees
  await page.mouse.move(0, 0).catch(() => {});

  await page.evaluate(async (overlays: string[]) => {
    const w = window as any;
    const grok = w.grok;

    for (const d of Array.from(grok.shell.dialogs ?? [])) {
      try { (d as any).close(); } catch (_) {}
    }
    const stuck = () => overlays.flatMap((s) => Array.from(document.querySelectorAll(s)));
    await new Promise<void>((resolve) => {
      const t0 = Date.now();
      const tick = () => {
        if (stuck().length === 0 || Date.now() - t0 > 1000) return resolve();
        setTimeout(tick, 50);
      };
      tick();
    });
    for (const el of stuck()) el.remove();

    grok.shell.closeAll();
    await new Promise<void>((resolve) => {
      const t0 = Date.now();
      const tick = () => {
        if (Array.from(grok.shell.tableViews).length === 0 || Date.now() - t0 > 3000) return resolve();
        setTimeout(tick, 50);
      };
      tick();
    });

    // shell settings the specs flip: openTable sets both, and leaving them on changes
    // what the NEXT spec sees before it has run a line
    try { grok.shell.windows.simpleMode = false; } catch (_) {}
    try { grok.shell.settings.showFiltersIconsConstantly = false; } catch (_) {}
    try { grok.shell.o = null; } catch (_) {}
    document.body.classList.remove('selenium');

    // the event-wait layer stamps state on window; a stale __lastRender makes the next
    // spec's first waitForViewerRendered resolve against a render that already happened
    delete w.__lastRender;
    delete w.__canvasColorSnap;
    delete w.__canvasQuiet;
  }, OVERLAYS);
}

export const test = base.extend<{page: Page}, {shared: {page: Page | null}}>({
  // Worker-scoped holder rather than a worker-scoped page: the context has to be built
  // from `contextOptions`, which is test-scoped, and building it by hand instead dropped
  // every project-level `use` — the Desktop Chrome device settings among them, which
  // changed how viewers rendered and failed specs that pass on their own.
  shared: [async ({}, use) => {
    const holder: {page: Page | null} = {page: null};
    await use(holder);
    await holder.page?.context().close().catch(() => {});
  }, {scope: 'worker'}],

  page: async ({browser, contextOptions, shared}, use) => {
    if (!shared.page) {
      const context = await browser.newContext(contextOptions);
      context.setDefaultTimeout(specTestOptions.actionTimeout);
      context.setDefaultNavigationTimeout(specTestOptions.navigationTimeout);
      shared.page = await context.newPage();
      await loginToDatagrok(shared.page);
    }
    stepErrors.length = 0;
    await revive(shared.page);
    await use(shared.page);
    await revive(shared.page);
  },
});

/**
 * Puts the shared page back into a usable state, re-booting it if it is past saving.
 *
 * Without this a single network blip ends the worker: the page dies, resetShell's evaluate
 * throws inside the fixture, and every remaining test in the file fails on a corpse. That
 * is the cost sharing adds over a per-test context, so it has to be paid back here.
 */
async function revive(page: Page): Promise<void> {
  try {
    await resetShell(page);
    await drainErrors(page);
    return;
  } catch (_) { /* fall through to a re-boot */ }
  await loginToDatagrok(page).catch(() => {});
  await resetShell(page).catch(() => {});
}

/**
 * Waits until the page stops emitting errors, so the previous test's fallout is not
 * charged to the next one.
 *
 * closeAll() cancels work that is still in flight, and the resulting console errors land
 * a beat later — after the next test has registered its own `page.on('console')` guard.
 * Specs then failed on `expect(errCount()).toBe(errBefore)` with counts they never caused,
 * which is what "Expected: 0, Received: 2" was in the shared-page run.
 */
async function drainErrors(page: Page, quietMs = 500, capMs = 3000): Promise<void> {
  let last = Date.now();
  const bump = () => { last = Date.now(); };
  page.on('console', bump);
  page.on('pageerror', bump);
  try {
    const deadline = Date.now() + capMs;
    while (Date.now() < deadline && Date.now() - last < quietMs)
      await page.waitForTimeout(50);
  } finally {
    page.off('console', bump);
    page.off('pageerror', bump);
  }
}

export {expect} from '@playwright/test';
