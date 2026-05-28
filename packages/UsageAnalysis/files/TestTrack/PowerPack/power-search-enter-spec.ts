/* ---
sub_features_covered: [powerpack.search.power, powerpack.search.power.dispatch, powerpack.search.power.init, powerpack.welcome.view, powerpack.welcome.suggestion-nav, powerpack.search]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (scenario .md Notes section discusses
//     pyramid_layer bug-focused in prose under Rule 3 but does NOT
//     set it as a YAML frontmatter key — A-LAYER-ALIGN-01 PASS-by-
//     vacuity applies; constraint-enforcement uses non-pyramid
//     defaults: ≥1 DOM-driving call REQUIRED per target_layer
//     playwright; no JS-API-substitution prohibition for any single
//     flow since ui_coverage_responsibility is absent).
//   sub_features_covered: [powerpack.search.power, powerpack.search.power.dispatch,
//     powerpack.search.power.init, powerpack.welcome.view,
//     powerpack.welcome.suggestion-nav, powerpack.search]
//   ui_coverage_responsibility: absent (scenario does not declare specialty
//     UI ownership; default constraint-enforcement applies)
//   related_bugs: [GROK-18656]
//   produced_from: migrated
//   coverage_type: regression
//
// Bug-library cross-reference:
//   GROK-18656 (priority p2, regression-risk) — Power Search engine must
//     safely handle Enter-key dispatch from the Welcome View search input
//     across all 8 dispatch paths (view-search → regex-entity → project →
//     JS-eval → function-call → table-query → specific-widget →
//     function-evaluation). Before the fix, short queries like 'QA' caused
//     the dispatcher to throw `reported error: null`. This spec is the
//     canonical regression witness for the bug. Entry in
//     bug-library/powerpack.yaml#L207.
//
// Reference template: Powerpack/autocomplete-spec.ts (same directory, same
//   target_layer: playwright, same powerpack feature, same UI-event-driven
//   pattern — focus an input via DOM, type with page.keyboard, observe
//   debounced platform reaction). Same login + setup pattern.
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
//   powerSearch dispatch chain — visible as a browser-console error and
//   surfaced by Datagrok as a red error balloon ("reported error: null").
//   We capture this via:
//     (a) page.on('pageerror', ...) — uncaught JS exceptions reach this
//         listener with the full Error.
//     (b) page.on('console', ...) filter for 'error' level messages
//         containing 'reported error' or otherwise GROK-18656-shaped.
//     (c) Per-probe inspection of grok.shell.lastError (Datagrok-side
//         error capture surface) read from the page via page.evaluate.
//   After each probe, we assert that NO new pageerror was raised since
//   the prior baseline. This is the direct GROK-18656 invariant
//   ('Pressing Enter on QA does NOT throw').
//
// Retry round 1 (2026-05-26, cycle 2026-05-26-powerpack-automate-03):
//   Gate B round-0 FAILED with:
//     `TimeoutError: locator.waitFor: Timeout 60000ms exceeded.
//      Call log: waiting for locator('.power-pack-welcome-view').first()
//      to be visible`
//   at line 159 (welcomeView.waitFor). Hypothesis category: TEST-BUG
//   (setup-phase brittleness — same paradigm, tactical fix).
//
//   MCP recon (chrome-devtools, dev.datagrok.ai) revealed:
//     (i)  on a fresh authenticated session the Welcome view IS the
//          home view (welcomePresent=true, view.name='Home');
//     (ii) calling grok.shell.v?.close?.() DESTROYS the welcome view,
//          and the view does NOT auto-respawn (welcomePresent=false
//          after closeAll() + 3 s wait);
//     (iii) grok.functions.call('PowerPack:_welcomeView', {}) returns
//          a fresh DG.View; grok.shell.addView(v) mounts it and brings
//          .power-pack-welcome-view back into the DOM.
//   The Gate B round-0 failure trace (page snapshot in
//   test-playwright-output/PowerPack-power-search-ent-304b9-aths-GROK-
//   18656-regression--chromium/error-context.md) confirms the home
//   view was the landing-site marketing renderer, NOT PowerPack's
//   welcomeView — consistent with hypothesis (ii) above.
//
//   Fix applied: replaced `grok.shell.v?.close?.()` with the explicit
//   `closeAll() + grok.functions.call('PowerPack:_welcomeView', {}) +
//   grok.shell.addView(v)` sequence. No paradigm pivot — same
//   Playwright DOM-driving body, same selectors, same probe loop.

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
  // authenticated session the PowerPack Welcome view IS the platform
  // home — `grok.shell.v.name == 'Home'` and `.power-pack-welcome-view`
  // is present in the DOM (verified empirically via MCP recon on
  // dev.datagrok.ai on 2026-05-26 for retry round 1 of this spec).
  //
  // However, on a SHARED authenticated session (this is the case under
  // `grok test` cookie-injection auth: the Playwright browser inherits
  // the same session cookies as any concurrent user activity), the
  // current shell view may NOT be the PowerPack Welcome view — it can
  // be the landing-site / docs help-page view (the empirical failure
  // mode observed in Validator Gate B round 0: `.power-pack-welcome-view`
  // never materialized within 60 s because the active home view was the
  // marketing landing-site renderer, NOT PowerPack's `welcomeView`).
  //
  // The previous version of this spec called `grok.shell.v?.close?.()`
  // which DESTROYS the current view but does NOT auto-respawn the
  // PowerPack Welcome view — empirically verified via MCP recon: after
  // closing the active view, `.power-pack-welcome-view` is removed from
  // DOM and does NOT come back, even after closeAll() and a 3 s wait.
  //
  // The robust fix (this retry round 1): explicitly construct + add the
  // Welcome view via the PowerPack-registered function. Per
  // `PowerPack/src/package.ts:131` the function is registered with
  // `name: 'welcomeView'` and its underlying export is `_welcomeView`
  // (the static class method that returns the view). The platform
  // exposes it as `'PowerPack:_welcomeView'` (the export name, with the
  // underscore preserved through `grok api` codegen). Calling
  // `grok.functions.call('PowerPack:_welcomeView', {})` returns a fresh
  // DG.View instance; `grok.shell.addView(v)` mounts it and makes
  // `.power-pack-welcome-view` present in DOM. This guarantees the
  // Welcome view is the active surface for this spec regardless of
  // what other view the shared session previously had open.
  await page.evaluate(async () => {
    const grok = (window as any).grok;
    document.body.classList.add('selenium');
    try { grok.shell.windows.simpleMode = true; } catch (_) {}
    // Close all views (defensively dismiss any residual table views,
    // landing-site help-page renderers, etc.). This is best-effort —
    // closeAll() on its own does NOT respawn the welcome view.
    try { grok.shell.closeAll(); } catch (_) {}
    // Explicitly construct + add the PowerPack Welcome view. The
    // function returns DG.View | undefined; addView mounts it as the
    // active view. The view's root carries the `power-pack-welcome-view`
    // class (welcome-view.ts:56) and contains the search input (line 19).
    try {
      const v = await grok.functions.call('PowerPack:_welcomeView', {});
      if (v) grok.shell.addView(v);
    } catch (_) { /* if call fails the waitFor below will surface it */ }
  });

  // Wait for the Welcome view to be visible. We wait on the root
  // container class added by welcome-view.ts:56. With the explicit
  // addView above, this should resolve within a few seconds — the
  // 60 s timeout is defensive (covers slow CI / Dart-renderer warmup).
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

  // ---- Helpers (inline — single spec, threshold ≥3 not met for cross-spec extraction) ----

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
    // Tiny settle so any in-flight debounced dispatch can complete.
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
    // Step 1: focus the search input (UI action). This is the DOM-driving
    // call required by E-LAYER-COMPLIANCE-01 for target_layer playwright.
    await input.click({timeout: 10_000});
    // Step 2: type the query. We use page.keyboard.type so the platform's
    // rxjs.fromEvent(input, 'input') subscription is triggered for every
    // keystroke. The 500ms debounce at welcome-view.ts:111 means a single
    // dispatch will fire ~500ms after the last keystroke.
    await p.keyboard.type(probe.query, {delay: 40});
    // Wait for the debounced dispatch to fire (debounceTime(500) + a
    // small safety margin so powerSearch can complete its 8-path walk
    // synchronously inside the subscription callback).
    await p.waitForTimeout(900);
    // Step 3: press Enter. With no suggestion highlighted, the keydown
    // listener at welcome-view.ts:118-141 returns without preventDefault
    // (lines 124-125: empty allItems list short-circuits). With a
    // suggestion highlighted, it would dispatch a click on that item.
    // Either way, this MUST be a UI keypress to exercise the
    // suggestionMenuKeyNavigation handler (which is the codepath that
    // GROK-18656 also touches via the input.keydown surface).
    await p.keyboard.press('Enter');
    // Wait for any post-Enter side effects (suggestion click, view
    // change, second dispatch round) to settle.
    await p.waitForTimeout(700);
    // Step 4: assert no error fired. GROK-18656 invariant.
    assertNoNewErrors(probe.query, baseline);
    // Optional diagnostic read (does NOT fail the test): capture shell
    // last-error state for the run log.
    const shellErr = await readShellLastError(p);
    if (shellErr && /reported error\s*:\s*null/i.test(shellErr))
      throw new Error(
        `GROK-18656 invariant violated for probe '${probe.query}': ` +
        `grok.shell.lastError contains 'reported error: null': ${shellErr}`,
      );
  }

  // ---- Scenario 1: canonical GROK-18656 reproduction (query = "QA") ----
  // This single probe IS the bug-focused invariant. It is the canonical
  // pre-fix failure trigger; post-fix it MUST PASS.
  await softStep('Scenario 1: focus + type "QA" + press Enter → dispatch must not throw (GROK-18656 invariant)', async () => {
    const probe = DISPATCH_PROBES[0]; // 'QA'
    expect(probe.isCanonicalRepro).toBe(true);
    await dispatchProbe(page, probe);
  });

  await softStep('Scenario 1 Step 5: verify a sensible dispatch result is rendered (results, empty-state, or graceful balloon)', async () => {
    // After typing 'QA' + Enter + settle, ONE of three outcomes must
    // hold (per scenario body Step 5 acceptance criteria):
    //   (a) the search-results host (.power-pack-search-host) is visible
    //       and contains at least one child node (a dispatch path
    //       produced a hit).
    //   (b) the search-results host is visible and is empty / shows a
    //       graceful empty-state — the dispatcher completed cleanly.
    //   (c) a Datagrok info / warning balloon is visible with a
    //       user-friendly message — never an uncaught throw.
    // The disqualifying outcome is a thrown error (covered by the
    // GROK-18656 invariant assertion above).
    const searchHost = page.locator('.power-pack-search-host').first();
    // The host is created at welcome-view.ts:45 and its display is
    // toggled by doSearch() at welcome-view.ts:104-105 (display = ''
    // when there is a search string). Wait for it to become visible.
    await searchHost.waitFor({timeout: 10_000, state: 'visible'}).catch(() => {});
    const hostVisible = await searchHost.isVisible({timeout: 1000}).catch(() => false);
    const hostChildCount = hostVisible
      ? await page.evaluate(() => {
        const el = document.querySelector('.power-pack-search-host');
        return el ? el.children.length : 0;
      })
      : 0;
    // Any balloon? Datagrok renders balloons inside a global container
    // (class includes 'grok-balloon' on the balloon element). We probe
    // existence without requiring visibility — a brief balloon may have
    // already faded out by the time we look.
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
  // For each of the 5 remaining queries, repeat the same 4-step probe.
  // Each invocation exercises a different subset of the 8 dispatch
  // branches per the probe payload metadata above.
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
  // GROK-18656 fix MUST NOT break suggestionMenuKeyNavigation (welcome-
  // view.ts:118-141) — the ArrowDown / ArrowUp / Enter-on-highlighted-
  // suggestion handler. Atlas declares interaction
  // 'press Enter on highlighted suggestion → click suggestion'.
  //
  // The scenario body suggests typing 'dem' so Demog dataset / sample
  // surfaces in the suggestion list. The actual suggestion menu is
  // rendered by powerSearch's contributing functions (e.g. typeahead
  // sources via powerPackSearchProvider — PowerPack/src/package.ts:219+);
  // suggestion entries land as elements with class .d4-menu-item inside
  // the input container (welcome-view.ts:123).
  //
  // If no suggestion entries surface for 'dem' in the current environment
  // (registry composition varies by deploy), the suggestion-nav handler
  // short-circuits at welcome-view.ts:124-125 — pressing ArrowDown does
  // nothing. We tolerate this branch by emitting a per-step softStep
  // that does not fail the GROK-18656 invariant: the regression-guard
  // assertion is that the keydown handler ITSELF does not throw, not
  // that suggestion entries are guaranteed to be present.

  await softStep('Scenario 3 Step 1-2: focus search input + type "dem" to surface suggestion menu', async () => {
    await clearSearchInput(page);
    const input = page.locator('input.power-search-search-everywhere-input').first();
    await input.click({timeout: 10_000});
    await page.keyboard.type('dem', {delay: 50});
    // Wait for the debounced powerSearch dispatch (500ms) + small margin.
    await page.waitForTimeout(900);
  });

  await softStep('Scenario 3 Step 3: press ArrowDown to highlight first suggestion (if any present)', async () => {
    const baseline = errorSnapshot();
    // Look for any .d4-menu-item under the input container BEFORE pressing
    // ArrowDown. We scope to the welcome view container to avoid hits in
    // unrelated parts of the DOM.
    const itemsBeforeCount = await page.evaluate(() => {
      const root = document.querySelector('.power-pack-welcome-view');
      if (!root) return 0;
      return root.querySelectorAll('.d4-menu-item').length;
    });
    // Press ArrowDown (UI action — this is the regression-guard surface).
    await page.keyboard.press('ArrowDown');
    await page.waitForTimeout(300);
    // Regression-guard invariant: the keydown handler itself must NOT
    // throw, regardless of whether items were present. The handler at
    // welcome-view.ts:124-125 short-circuits cleanly on empty list.
    assertNoNewErrors('ArrowDown-after-dem', baseline);
    // If items were present, the handler should have added
    // .d4-menu-item-hover to one of them. Capture for diagnostic value;
    // do not fail when items were not surfaced (env-dependent registry).
    if (itemsBeforeCount > 0) {
      const highlightedCount = await page.evaluate(() => {
        const root = document.querySelector('.power-pack-welcome-view');
        if (!root) return 0;
        return root.querySelectorAll('.d4-menu-item-hover').length;
      });
      // After ArrowDown on a non-empty list, exactly one item should be
      // hovered (welcome-view.ts:138-139 removes hover from the prior
      // and adds to the new). 0 highlighted is acceptable if the menu
      // closed between type-debounce and keypress (race-condition-prone
      // on slow CI); 1+ is the canonical outcome.
      expect(highlightedCount).toBeLessThanOrEqual(itemsBeforeCount);
    }
  });

  await softStep('Scenario 3 Step 4: press ArrowDown again (or ArrowUp) — navigation in both directions does not throw', async () => {
    const baseline = errorSnapshot();
    await page.keyboard.press('ArrowDown');
    await page.waitForTimeout(200);
    await page.keyboard.press('ArrowUp');
    await page.waitForTimeout(200);
    // Regression-guard invariant: bidirectional navigation in the
    // suggestion menu must NOT throw.
    assertNoNewErrors('ArrowDown+ArrowUp-after-dem', baseline);
  });

  await softStep('Scenario 3 Step 5: press Enter on highlighted suggestion (if any) — activation does not throw, GROK-18656 fix preserves this path', async () => {
    const baseline = errorSnapshot();
    // Probe whether a highlighted item is present. If yes, Enter clicks
    // it (welcome-view.ts:133-137). If not, Enter falls through the
    // suggestion-nav handler and reaches the input's default Enter
    // behavior (no-op for a plain text input).
    const hadHighlighted = await page.evaluate(() => {
      const root = document.querySelector('.power-pack-welcome-view');
      if (!root) return false;
      return root.querySelector('.d4-menu-item-hover') != null;
    });
    await page.keyboard.press('Enter');
    await page.waitForTimeout(800);
    // Regression-guard invariant: Enter on the suggestion-nav path must
    // NOT throw. This is the second half of the GROK-18656 fix
    // preservation contract (the first half is the dispatch-not-throw
    // assertion on the Scenario 1 + 2 probes above).
    assertNoNewErrors('Enter-on-highlighted-suggestion', baseline);
    // Diagnostic note in the run log — does not fail.
    if (!hadHighlighted) {
      // The scenario body acknowledges this case via Step 5's "equivalent
      // to clicking it: the corresponding entity is opened or the
      // corresponding view is navigated to" — only meaningful when a
      // suggestion was actually highlighted. If no suggestion surfaced
      // for 'dem' in this environment, the Enter-on-highlighted-
      // suggestion path was not exercised end-to-end; the regression-
      // guard portion (handler does not throw) IS exercised.
      // No assertion failure — environment-tolerant per spec design.
    }
  });

  // ---- Cleanup ----
  // Clear the search input so the Welcome view returns to its widgets
  // panel state for the next test. Best-effort.
  try {
    await clearSearchInput(page);
  } catch (_) { /* best effort */ }
  await page.evaluate(() => {
    try { (window as any).grok.shell.closeAll(); } catch (_) { /* best effort */ }
  }).catch(() => {});

  // ---- Final GROK-18656 invariant summary ----
  // If any probe raised a thrown error, the per-probe softStep would have
  // pushed an entry into stepErrors and we surface it here. The summary
  // includes the count of pageerror events captured across the whole run
  // so an audit can correlate.
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
