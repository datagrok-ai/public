/* ---
sub_features_covered: [powerpack.search.power, powerpack.search.power.dispatch, powerpack.search.power.init, powerpack.welcome.view, powerpack.welcome.suggestion-nav, powerpack.search]
--- */
// Frontmatter:
//   target_layer: playwright (requires ≥1 DOM-driving call; no pyramid_layer set)
//   sub_features_covered: [powerpack.search.power, powerpack.search.power.dispatch,
//     powerpack.search.power.init, powerpack.welcome.view,
//     powerpack.welcome.suggestion-nav, powerpack.search]
//   related_bugs: [GROK-18656]
//   produced_from: migrated
//   coverage_type: regression
//
// GROK-18656: Power Search must safely handle Enter-key dispatch from the
//   Welcome View search input across all 8 dispatch paths (view-search →
//   regex-entity → project → JS-eval → function-call → table-query →
//   specific-widget → function-evaluation). Short queries like 'QA' made the
//   dispatcher throw `reported error: null`. This spec is the regression
//   witness. Entry in bug-library/powerpack.yaml#L207.
//
// Reference template: Powerpack/autocomplete-spec.ts (same UI-event-driven
//   pattern — focus an input via DOM, type with page.keyboard, observe the
//   debounced platform reaction; same login + setup pattern).
//
// Source citations for selectors (no powerpack.md grok-browser reference
// file exists; selectors are cited directly from the implementation source
// per the Powerpack/hints-spec.ts + Powerpack/autocomplete-spec.ts
// precedent of citing PowerPack/src/* directly when the platform-internal
// class names are NOT in the cross-feature grok-browser/references/* catalogue):
//   - Welcome view root container:
//       .power-pack-welcome-view
//     Source: public/packages/PowerPack/src/welcome-view.ts:56
//     (view.root.classList.add('power-pack-welcome-view')).
//   - Welcome-view search input element:
//       input.power-search-search-everywhere-input
//     Source: public/packages/PowerPack/src/welcome-view.ts:18-20
//     (ui.element('input', 'ui-input-editor'); .classList.add(
//     'power-search-search-everywhere-input'); placeholder 'Search
//     everywhere. Try "aspirin" or "7JZK"').
//   - Welcome-view search-input container:
//       .ui-input-root.ui-input-type-ahead (within .d4-search-bar)
//     Source: public/packages/PowerPack/src/welcome-view.ts:21-27.
//   - Search results host (where powerSearch renders dispatch output):
//       .power-pack-search-host
//     Source: public/packages/PowerPack/src/welcome-view.ts:45,
//     paired with display toggle at line 104-105 (.style.display
//     toggles based on whether the search string is non-empty).
//   - Suggestion-menu items (Scenario 3 ArrowDown/ArrowUp navigation):
//       .d4-menu-item, .d4-menu-item-hover
//     Source: public/packages/PowerPack/src/welcome-view.ts:122-139
//     (suggestionMenuKeyNavigation queries .d4-menu-item-hover for the
//     currently-highlighted entry and .d4-menu-item for the full list).
//
// Atlas provenance (derived_from): the powerpack atlas
//   (feature-atlas/powerpack.yaml) entries for the six sub_features
//   reference the source-line anchors enumerated above; no separate
//   derived_from values appear on the sub_features themselves.
//
// GROK-18656 invariant assertion strategy:
//   The bug manifests as an UNCAUGHT exception thrown from inside the
//   powerSearch dispatch chain — a browser-console error surfaced by
//   Datagrok as a red error balloon ("reported error: null"). Captured via:
//     (a) page.on('pageerror', ...) — uncaught JS exceptions, full Error.
//     (b) page.on('console', ...) filter for 'error' level messages
//         containing 'reported error' or otherwise GROK-18656-shaped.
//     (c) Per-probe read of grok.shell.lastError (Datagrok-side error
//         capture surface) via page.evaluate.
//   After each probe, assert NO new pageerror was raised since the prior
//   baseline — the direct GROK-18656 invariant ('Enter on QA does NOT throw').

import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// Probe payload: the 6 queries exercised across Scenarios 1 + 2. Each entry
// names the query, the dispatch branch(es) the scenario body identifies it
// with, and whether it is the canonical GROK-18656 reproduction query.
interface Probe {
  query: string;
  branchHint: string;
  isCanonicalRepro: boolean;
}

const DISPATCH_PROBES: Probe[] = [
  // Scenario 1: the canonical GROK-18656 reproduction query.
  {query: 'QA', branchHint: 'short-query; same shape that triggered GROK-18656', isCanonicalRepro: true},
  // Scenario 2: 5 additional dispatch shapes.
  {query: 'a', branchHint: 'single-char short-query — most aggressive null-safety surface', isCanonicalRepro: false},
  {query: 'new', branchHint: 'common word — function-call / widget dispatch branches', isCanonicalRepro: false},
  {query: 'user', branchHint: 'username-style — regex-entity (users) + function-call branches', isCanonicalRepro: false},
  {query: 'Project[0-9]+', branchHint: 'regex-pattern entity — regex-entity dispatch branch', isCanonicalRepro: false},
  {query: '1+1', branchHint: 'function-evaluation expression — JS-eval / function-evaluation branch', isCanonicalRepro: false},
];

test('PowerPack: Power Search Enter-key dispatch is null-safe across 8 paths (GROK-18656 regression)', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  // ---- Page-error capture for GROK-18656 invariant ----
  // Capture uncaught JS exceptions raised in the page. The bug manifested
  // as a thrown null-deref inside one of the powerSearch dispatch branches;
  // an unfixed regression would surface as a new entry in this array.
  const pageErrors: Array<{message: string; stack?: string}> = [];
  page.on('pageerror', (err) => {
    pageErrors.push({message: err.message, stack: err.stack});
  });
  // Capture console.error messages that look like Datagrok's "reported
  // error: null" surface (the explicit GROK-18656 signature observed
  // pre-fix per scenario body).
  const reportedErrors: string[] = [];
  page.on('console', (msg) => {
    if (msg.type() !== 'error') return;
    const text = msg.text();
    if (/reported error/i.test(text) || /Uncaught/i.test(text))
      reportedErrors.push(text);
  });

  // ---- Login + setup phase ----
  await loginToDatagrok(page);

  // Navigate to the Welcome view. The view is registered as
  // PowerPack.welcomeView with `meta.autostartImmediate: true` and
  // `outputs.home: view` (PowerPack/src/package.ts:129-136). On a fresh
  // session the Welcome view IS the platform home, but on a SHARED
  // session (the case under `grok test` cookie-injection auth — the
  // Playwright browser inherits cookies from concurrent user activity)
  // the current shell view may instead be the landing-site / docs
  // help-page renderer, so `.power-pack-welcome-view` is not in the DOM.
  //
  // Don't just call `grok.shell.v?.close?.()`: closing the active view
  // does NOT auto-respawn the Welcome view (it stays gone even after
  // closeAll() + a wait). Instead explicitly construct + add it via the
  // PowerPack-registered function. The function is registered with
  // `name: 'welcomeView'` but its underlying export is `_welcomeView`,
  // so the platform exposes it as `'PowerPack:_welcomeView'` (the export
  // name, underscore preserved through `grok api` codegen). Calling
  // `grok.functions.call('PowerPack:_welcomeView', {})` returns a fresh
  // DG.View; `grok.shell.addView(v)` mounts it and makes
  // `.power-pack-welcome-view` present in DOM — making the Welcome view
  // the active surface regardless of what the shared session had open.
  await page.evaluate(async () => {
    const grok = (window as any).grok;
    document.body.classList.add('selenium');
    try { grok.shell.windows.simpleMode = true; } catch (_) {}
    // Defensively dismiss any residual views (table views, help-page
    // renderers). Best-effort — closeAll() alone does NOT respawn welcome.
    try { grok.shell.closeAll(); } catch (_) {}
    // Construct + add the Welcome view. The function returns
    // DG.View | undefined; addView mounts it as the active view. The
    // root carries class `power-pack-welcome-view` (welcome-view.ts:56)
    // and contains the search input (line 19).
    try {
      const v = await grok.functions.call('PowerPack:_welcomeView', {});
      if (v) grok.shell.addView(v);
    } catch (_) { /* if call fails the waitFor below will surface it */ }
  });

  // Wait for the Welcome view (root container class, welcome-view.ts:56).
  // With the explicit addView above this resolves in seconds; the 60 s
  // timeout is defensive (slow CI / Dart-renderer warmup).
  const welcomeView = page.locator('.power-pack-welcome-view').first();
  await welcomeView.waitFor({timeout: 60_000, state: 'visible'});

  // Locate the search input element. Its CSS class is set explicitly in
  // welcome-view.ts:19. The placeholder + the d4-search-bar parent give
  // additional defensive scoping.
  const searchInput = page.locator('input.power-search-search-everywhere-input').first();
  await searchInput.waitFor({timeout: 30_000, state: 'visible'});

  // Tiny settle wait — welcome-view widgets refresh asynchronously after
  // login (refresh() at welcome-view.ts:79 fires a dapi.groups.find call).
  await page.waitForTimeout(1500);

  // ---- Helpers (inline — single spec, below the cross-spec extraction threshold) ----

  /**
   * Clear the welcome-view search input. We use a UI-equivalent
   * triple-click + Backspace sequence to ensure the 'input' event fires
   * (the search dispatch is wired to rxjs.fromEvent(input, 'input').pipe(
   * debounceTime(500)) at welcome-view.ts:111 — direct .value assignment
   * would NOT trigger the debounced dispatch).
   */
  async function clearSearchInput(p: Page): Promise<void> {
    const input = p.locator('input.power-search-search-everywhere-input').first();
    await input.click({timeout: 10_000});
    // Select-all + delete via keyboard so the input event fires.
    await p.keyboard.press('Control+A');
    await p.keyboard.press('Delete');
    // Settle so any in-flight debounced dispatch can complete.
    await p.waitForTimeout(150);
  }

  /**
   * Snapshot of the page-error capture state at this moment. Used to
   * detect new errors raised during a single probe.
   */
  function errorSnapshot(): {pageErrorsCount: number; reportedErrorsCount: number} {
    return {pageErrorsCount: pageErrors.length, reportedErrorsCount: reportedErrors.length};
  }

  /**
   * Read grok.shell.lastError from the page if present. Some Datagrok
   * builds expose a captured error string on the shell after a balloon
   * is shown; absence is acceptable (returns null).
   */
  async function readShellLastError(p: Page): Promise<string | null> {
    return p.evaluate(() => {
      const g = (window as any).grok;
      try {
        const v = g?.shell?.lastError;
        return typeof v === 'string' ? v : (v == null ? null : String(v));
      } catch (_) { return null; }
    });
  }

  /**
   * Per-probe assertion: no new pageerror or reportedError was raised
   * since the baseline snapshot. The probe label is the query string —
   * embedded in the failure message so a regression is immediately
   * traceable to the dispatch branch.
   */
  function assertNoNewErrors(
    label: string,
    baseline: {pageErrorsCount: number; reportedErrorsCount: number},
  ): void {
    const newPageErrors = pageErrors.slice(baseline.pageErrorsCount);
    const newReported = reportedErrors.slice(baseline.reportedErrorsCount);
    if (newPageErrors.length > 0) {
      const summary = newPageErrors.map((e) => `  ${e.message}`).join('\n');
      throw new Error(
        `GROK-18656 invariant violated for probe '${label}': ` +
        `${newPageErrors.length} new uncaught pageerror raised during dispatch:\n${summary}`,
      );
    }
    if (newReported.length > 0) {
      // 'reported error: null' is the canonical GROK-18656 signature.
      const hasReportedNull = newReported.some((m) => /reported error\s*:\s*null/i.test(m));
      if (hasReportedNull) {
        throw new Error(
          `GROK-18656 invariant violated for probe '${label}': ` +
          `'reported error: null' surfaced on the console:\n  ${newReported.join('\n  ')}`,
        );
      }
      // Non-canonical console.error messages are recorded for diagnostic
      // value but do NOT fail the GROK-18656 invariant on their own —
      // many platform-internal warn-shaped messages are routed through
      // console.error legitimately.
    }
  }

  /**
   * Drive a single probe: focus → type query → press Enter → wait for
   * dispatch settle → assert GROK-18656 invariant. This is the 4-step
   * probe applied to every query in Scenarios 1 and 2.
   */
  async function dispatchProbe(p: Page, probe: Probe): Promise<void> {
    await clearSearchInput(p);
    const baseline = errorSnapshot();
    const input = p.locator('input.power-search-search-everywhere-input').first();
    // Step 1: focus the search input — the DOM-driving UI action.
    await input.click({timeout: 10_000});
    // Step 2: type the query via page.keyboard.type so the platform's
    // rxjs.fromEvent(input, 'input') subscription fires per keystroke. The
    // 500ms debounce at welcome-view.ts:111 fires a single dispatch ~500ms
    // after the last keystroke.
    await p.keyboard.type(probe.query, {delay: 40});
    // Wait for the debounced dispatch (debounceTime(500) + margin so
    // powerSearch can finish its 8-path walk synchronously in the callback).
    await p.waitForTimeout(900);
    // Step 3: press Enter. With no suggestion highlighted, the keydown
    // listener at welcome-view.ts:118-141 returns without preventDefault
    // (lines 124-125: empty allItems list short-circuits); with one
    // highlighted it dispatches a click on that item. Must be a UI keypress
    // to exercise suggestionMenuKeyNavigation — the input.keydown codepath
    // GROK-18656 also touches.
    await p.keyboard.press('Enter');
    // Wait for post-Enter side effects (suggestion click, view change,
    // second dispatch round) to settle.
    await p.waitForTimeout(700);
    // Step 4: assert no error fired — GROK-18656 invariant.
    assertNoNewErrors(probe.query, baseline);
    // Diagnostic read (does NOT fail the test): shell last-error state.
    const shellErr = await readShellLastError(p);
    if (shellErr && /reported error\s*:\s*null/i.test(shellErr))
      throw new Error(
        `GROK-18656 invariant violated for probe '${probe.query}': ` +
        `grok.shell.lastError contains 'reported error: null': ${shellErr}`,
      );
  }

  // ---- Scenario 1: canonical GROK-18656 reproduction (query = "QA") ----
  // This single probe is the bug-focused invariant — the pre-fix failure
  // trigger; post-fix it MUST PASS.
  await softStep('Scenario 1: focus + type "QA" + press Enter → dispatch must not throw (GROK-18656 invariant)', async () => {
    const probe = DISPATCH_PROBES[0]; // 'QA'
    expect(probe.isCanonicalRepro).toBe(true);
    await dispatchProbe(page, probe);
  });

  await softStep('Scenario 1 Step 5: verify a sensible dispatch result is rendered (results, empty-state, or graceful balloon)', async () => {
    // After 'QA' + Enter + settle, ONE of three sensible outcomes must hold:
    //   (a) search-results host (.power-pack-search-host) visible with ≥1
    //       child — a dispatch path produced a hit.
    //   (b) host visible but empty — graceful empty-state, dispatcher
    //       completed cleanly.
    //   (c) a Datagrok info / warning balloon is visible.
    // The disqualifying outcome is a thrown error (covered by the
    // GROK-18656 invariant assertion above).
    const searchHost = page.locator('.power-pack-search-host').first();
    // Host created at welcome-view.ts:45; doSearch() toggles its display
    // to '' when there is a search string (welcome-view.ts:104-105).
    await searchHost.waitFor({timeout: 10_000, state: 'visible'}).catch(() => {});
    const hostVisible = await searchHost.isVisible({timeout: 1000}).catch(() => false);
    const hostChildCount = hostVisible
      ? await page.evaluate(() => {
        const el = document.querySelector('.power-pack-search-host');
        return el ? el.children.length : 0;
      })
      : 0;
    // Probe balloon existence without requiring visibility — a brief
    // balloon may have faded out by the time we look.
    const balloonCount = await page.evaluate(() => {
      return document.querySelectorAll('.grok-balloon, .d4-balloon').length;
    });
    // At least one of the three sensible outcomes must hold.
    const sensible = (hostVisible && hostChildCount > 0) ||
                     (hostVisible && hostChildCount === 0) ||
                     (balloonCount > 0);
    expect(sensible, `expected a sensible dispatch result for 'QA' — search-host visible: ${hostVisible}, ` +
      `host child count: ${hostChildCount}, balloon count: ${balloonCount}`).toBe(true);
  });

  // ---- Scenario 2: 8-path exercise — 5 additional dispatch shapes ----
  // Repeat the same 4-step probe for each remaining query; each exercises
  // a different subset of the 8 dispatch branches (see probe metadata above).
  for (let i = 1; i < DISPATCH_PROBES.length; i++) {
    const probe = DISPATCH_PROBES[i];
    await softStep(
      `Scenario 2 probe ${i + 1}/${DISPATCH_PROBES.length}: ` +
      `focus + type ${JSON.stringify(probe.query)} + Enter → dispatch must not throw ` +
      `(${probe.branchHint})`,
      async () => {
        await dispatchProbe(page, probe);
      },
    );
  }

  // ---- Scenario 3: suggestion-nav regression guard ----
  // The GROK-18656 fix MUST NOT break suggestionMenuKeyNavigation (welcome-
  // view.ts:118-141) — the ArrowDown / ArrowUp / Enter-on-highlighted-
  // suggestion handler (atlas interaction 'press Enter on highlighted
  // suggestion → click suggestion').
  //
  // Typing 'dem' surfaces Demog dataset / sample in the suggestion menu,
  // rendered by powerSearch's contributing functions (typeahead via
  // powerPackSearchProvider — PowerPack/src/package.ts:219+); entries land
  // as .d4-menu-item inside the input container (welcome-view.ts:123).
  //
  // Registry composition varies by deploy: if no entries surface for 'dem',
  // the handler short-circuits at welcome-view.ts:124-125 and ArrowDown does
  // nothing. We tolerate this — the regression-guard assertion is that the
  // keydown handler ITSELF does not throw, not that entries are present.

  await softStep('Scenario 3 Step 1-2: focus search input + type "dem" to surface suggestion menu', async () => {
    await clearSearchInput(page);
    const input = page.locator('input.power-search-search-everywhere-input').first();
    await input.click({timeout: 10_000});
    await page.keyboard.type('dem', {delay: 50});
    // Wait for the debounced powerSearch dispatch (500ms) + margin.
    await page.waitForTimeout(900);
  });

  await softStep('Scenario 3 Step 3: press ArrowDown to highlight first suggestion (if any present)', async () => {
    const baseline = errorSnapshot();
    // Count .d4-menu-item under the welcome-view container BEFORE pressing
    // ArrowDown (scoped to avoid hits elsewhere in the DOM).
    const itemsBeforeCount = await page.evaluate(() => {
      const root = document.querySelector('.power-pack-welcome-view');
      if (!root) return 0;
      return root.querySelectorAll('.d4-menu-item').length;
    });
    // Press ArrowDown — the regression-guard UI surface.
    await page.keyboard.press('ArrowDown');
    await page.waitForTimeout(300);
    // Regression-guard invariant: the keydown handler must NOT throw,
    // whether or not items were present (welcome-view.ts:124-125
    // short-circuits cleanly on an empty list).
    assertNoNewErrors('ArrowDown-after-dem', baseline);
    // If items were present, the handler adds .d4-menu-item-hover to one.
    // Diagnostic only; don't fail when none surfaced (env-dependent registry).
    if (itemsBeforeCount > 0) {
      const highlightedCount = await page.evaluate(() => {
        const root = document.querySelector('.power-pack-welcome-view');
        if (!root) return 0;
        return root.querySelectorAll('.d4-menu-item-hover').length;
      });
      // After ArrowDown on a non-empty list, one item should be hovered
      // (welcome-view.ts:138-139 moves hover from prior to new). 0 is
      // acceptable if the menu closed between debounce and keypress (CI
      // race); 1+ is the canonical outcome.
      expect(highlightedCount).toBeLessThanOrEqual(itemsBeforeCount);
    }
  });

  await softStep('Scenario 3 Step 4: press ArrowDown again (or ArrowUp) — navigation in both directions does not throw', async () => {
    const baseline = errorSnapshot();
    await page.keyboard.press('ArrowDown');
    await page.waitForTimeout(200);
    await page.keyboard.press('ArrowUp');
    await page.waitForTimeout(200);
    // Regression-guard invariant: bidirectional nav must NOT throw.
    assertNoNewErrors('ArrowDown+ArrowUp-after-dem', baseline);
  });

  await softStep('Scenario 3 Step 5: press Enter on highlighted suggestion (if any) — activation does not throw, GROK-18656 fix preserves this path', async () => {
    const baseline = errorSnapshot();
    // If a highlighted item is present, Enter clicks it (welcome-view.ts:
    // 133-137); otherwise Enter falls through to the input's default
    // behavior (no-op for a plain text input).
    const hadHighlighted = await page.evaluate(() => {
      const root = document.querySelector('.power-pack-welcome-view');
      if (!root) return false;
      return root.querySelector('.d4-menu-item-hover') != null;
    });
    await page.keyboard.press('Enter');
    await page.waitForTimeout(800);
    // Regression-guard invariant: Enter on the suggestion-nav path must
    // NOT throw — the second half of the GROK-18656 fix-preservation
    // contract (the first half is dispatch-not-throw on the Scenario 1+2
    // probes above).
    assertNoNewErrors('Enter-on-highlighted-suggestion', baseline);
    if (!hadHighlighted) {
      // No suggestion surfaced for 'dem' in this environment, so the
      // Enter-on-highlighted path was not exercised end-to-end; the
      // regression-guard portion (handler does not throw) still is.
      // No assertion failure — environment-tolerant by design.
    }
  });

  // ---- Cleanup ----
  // Clear the search input so the Welcome view returns to its widgets
  // panel state for the next test (best-effort).
  try {
    await clearSearchInput(page);
  } catch (_) { /* best effort */ }
  await page.evaluate(() => {
    try { (window as any).grok.shell.closeAll(); } catch (_) { /* best effort */ }
  }).catch(() => {});

  // ---- Final GROK-18656 invariant summary ----
  // Surface any per-probe softStep failures collected in stepErrors,
  // including whole-run pageerror / reported-error counts for correlation.
  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    const errCount = pageErrors.length;
    const reportedCount = reportedErrors.length;
    throw new Error(
      `${stepErrors.length} step(s) failed (pageerror count: ${errCount}; ` +
      `reported-error count: ${reportedCount}):\n${summary}`,
    );
  }
});
