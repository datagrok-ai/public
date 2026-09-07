import {test as base, Page} from '@playwright/test';
import {installLedger, ledgerAnnotations, openDatagrok, setLane, specTestOptions, stepErrors} from './spec-login';
import {drainPendingDeletes} from './helpers/projects';

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

// .d4-balloon earns its place: a balloon raised through Balloon.error/warning is STICKY
// (helpers/balloons.ts), so one test's error toast sits over the next test's UI and its
// container eats the clicks — "d4-balloon-container intercepts pointer events" was the
// shared-page click timeout in scatter-plot, statistics and tile-viewer.
const OVERLAYS = ['.d4-menu-popup', '.d4-tooltip', '.d4-dialog', '.d4-balloon'];


// The column-selector backdrop is deliberately absent from OVERLAYS: the platform puts that
// class on a wrapper AROUND the viewer, so removing the element removes the viewer. Left
// stuck it also makes pickColumnViaSelectorTrusted believe a popup opened when none did,
// so the class is stripped below instead.

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
    // visible ones only: the tooltip is a persistent singleton that merely hides, and treating
    // it as stuck burned the full 1s here on every test and then ripped the node out
    const stuck = () => overlays.flatMap((s) => Array.from(document.querySelectorAll(s)))
      .filter((e) => (e as HTMLElement).offsetParent !== null);
    await new Promise<void>((resolve) => {
      const t0 = Date.now();
      const tick = () => {
        if (stuck().length === 0 || Date.now() - t0 > 1000) return resolve();
        setTimeout(tick, 50);
      };
      tick();
    });
    for (const el of stuck()) el.remove();
    // the hidden tooltip singleton keeps its last text, and a later "no tooltip" read sees it
    try { w.ui.tooltip.hide(); } catch (_) {}
    for (const t of Array.from(document.querySelectorAll('.d4-tooltip')))
      if ((t as HTMLElement).offsetParent === null) (t as HTMLElement).innerHTML = '';
    for (const e of Array.from(document.querySelectorAll('.d4-column-selector-backdrop')))
      e.classList.remove('d4-column-selector-backdrop');
    // the container itself is what intercepts, and it outlives its children
    for (const c of Array.from(document.querySelectorAll('.d4-balloon-container')))
      (c as HTMLElement).innerHTML = '';

    grok.shell.closeAll();
    await new Promise<void>((resolve) => {
      const t0 = Date.now();
      const tick = () => {
        if (Array.from(grok.shell.tableViews).length === 0 || Date.now() - t0 > 3000) return resolve();
        setTimeout(tick, 50);
      };
      tick();
    });

    // a closed viewer's property grid outlives it in the context panel and `shell.o = null` is
    // ignored, so the next spec's property helpers edit a dead grid; the accordion section state
    // is persisted per column and gives later specs a live filter widget they never opened
    for (const g of Array.from(document.querySelectorAll('.property-grid'))) g.remove();
    try { for (const k of Object.keys(localStorage)) if (k.startsWith('Accordion:')) localStorage.removeItem(k); } catch (_) {}

    try { grok.shell.windows.simpleMode = false; } catch (_) {}
    // with the context panel hidden, `grok.shell.o = viewer` is ignored and no property grid is
    // built, which is how a spec that hides it starved the next one's property-panel step
    try { grok.shell.windows.showContextPanel = true; } catch (_) {}
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

// fixture-side steps are dropped by the JSON reporter, so their cost is recorded as an annotation
async function timed(title: string, fn: () => Promise<void>): Promise<void> {
  const t0 = Date.now();
  try { await fn(); }
  finally { try { base.info().annotations.push({type: 'fixture', description: `${Date.now() - t0}ms ${title}`}); } catch (_) {} }
}

function laneTest(lane: 'local' | 'server') {
  return base.extend<{page: Page}, {shared: {page: Page | null; dirty: boolean}}>({
    // Worker-scoped holder rather than a worker-scoped page: the context has to be built
    // from `contextOptions`, which is test-scoped, and building it by hand instead dropped
    // every project-level `use` — the Desktop Chrome device settings among them, which
    // changed how viewers rendered and failed specs that pass on their own.
    shared: [async ({}, use) => {
      const holder: {page: Page | null; dirty: boolean} = {page: null, dirty: true};
      await use(holder);
      if (holder.page) await drainPendingDeletes(holder.page);
      await holder.page?.context().close().catch(() => {});
    }, {scope: 'worker'}],

    page: async ({browser, contextOptions, shared}, use) => {
      if (!shared.page) {
        // the boot runs inside the first test's budget; a slow stand (60s+ at two workers) must not
        // fail that test at the 60s config timeout before its first step
        base.info().setTimeout(base.info().timeout + 120_000);
        const context = await browser.newContext(contextOptions);
        context.setDefaultTimeout(specTestOptions.actionTimeout);
        context.setDefaultNavigationTimeout(specTestOptions.navigationTimeout);
        shared.page = await context.newPage();
        setLane(shared.page, lane);
        installLedger(shared.page);
        await timed('fixture: boot ' + lane, () => openDatagrok(shared.page!));
      }
      // a spec may raise the page's default timeouts for itself (trellis: 120s); they must not
      // outlive it, or every failed locator in the next spec waits 120s instead of 15s
      shared.page.setDefaultTimeout(specTestOptions.actionTimeout);
      shared.page.setDefaultNavigationTimeout(specTestOptions.navigationTimeout);
      stepErrors.length = 0;
      // the page is clean unless the previous test's teardown never completed
      if (shared.dirty) await timed('fixture: revive before', () => revive(shared.page!, true));
      shared.dirty = true;
      // the error drain after the test only earns its 300ms quiet window when the test produced
      // errors; a clean test had nothing in flight to wait out (155 x 0.3s in the final run)
      let noisy = false;
      const mark = (m: any) => { if (typeof m?.type !== 'function' || m.type() === 'error') noisy = true; };
      shared.page.on('console', mark);
      shared.page.on('pageerror', mark);
      try { await use(shared.page); }
      finally { shared.page.off('console', mark); shared.page.off('pageerror', mark); }
      await timed('fixture: revive after', () => revive(shared.page!, noisy));
      shared.dirty = false;
      for (const a of ledgerAnnotations((shared.page as any).__ledger)) base.info().annotations.push(a);
    },
  });
}

/**
 * The server lane: an authenticated client, one boot per worker. What every spec used
 * before lanes existed, and what a spec whose subject is server state still needs.
 */
export const test = laneTest('server');

/**
 * The local lane: `?mode=local`, no session, no server (core/docs/features/ui2/LOCAL_MODE.md).
 * For a spec whose subject is client behaviour. Its holder is separate from the server lane's,
 * so a worker that runs both boots one page of each rather than re-navigating between them.
 */
export const localTest = laneTest('local');

/**
 * Puts the shared page back into a usable state, re-booting it if it is past saving.
 *
 * Without this a single network blip ends the worker: the page dies, resetShell's evaluate
 * throws inside the fixture, and every remaining test in the file fails on a corpse. That
 * is the cost sharing adds over a per-test context, so it has to be paid back here.
 */
async function revive(page: Page, drain: boolean): Promise<void> {
  try {
    await resetShell(page);
    if (drain) await drainErrors(page);
    return;
  } catch (_) { /* fall through to a re-boot */ }
  await openDatagrok(page).catch(() => {});
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
async function drainErrors(page: Page, quietMs = 300, capMs = 3000): Promise<void> {
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
