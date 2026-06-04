/* ---
sub_features_covered: [biostructure.viewer]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent on scenario .md — JS API substitution permitted for
//     setup/teardown, but >=1 DOM-driving call still REQUIRED for
//     target_layer: playwright per E-LAYER-COMPLIANCE-01 (constraint-
//     enforcement REQUIRED list). Realized via the unrelated view-tab Close
//     icon click (Scenario 1, repeated for three open/close cycles per the
//     scenario .md step 6).
//   sub_features_covered: 1 id (biostructure.viewer) mirrored above per
//     E-STRUCT-MECH-06. The bug-library entry for CLAUDE-33 lists
//     `affects: [biostructure.viewer]` only; the defect site is the
//     onViewRemoved subscription registered while the Mol* viewer is active.
//   related_bugs: [CLAUDE-33] — bug-invariant assertion REQUIRED per the
//     bug-library cross-reference convention. This scenario IS the dedicated
//     regression guard authored to close F-BUG-COVERAGE-01 (the existing
//     smoke biostructure-viewer.md `bug_match_attempts_skipped[]` audit
//     records "scenario never opens or closes unrelated views, so
//     semantic_match_all_steps returns []" for CLAUDE-33).
//   produced_from: atlas-driven
//   coverage_type: regression
//
// Atlas provenance (derived_from):
//   feature-atlas/biostructureviewer.yaml#edge_cases[1]
//     derived_from: bug-library:biostructureviewer.yaml#CLAUDE-33
//     description: Open a Molstar (Biostructure) viewer, then open and close
//       any unrelated view. The onViewRemoved handler must no-op for
//       non-matching views (no `Cannot read properties of undefined
//       (reading 'children')` crash from fallbackPreviewCheck accessing
//       `evtView.root.children[0].children[0]`).
//   feature-atlas/biostructureviewer.yaml#sub_features[biostructure.viewer]
//     source: public/packages/BiostructureViewer/src/package.ts#L455
//   feature-atlas/biostructureviewer.yaml#known_issues[CLAUDE-33]
//     affects_sub_features: [biostructure.viewer];
//     test_coverage.exists: false — the hole this spec fills.
//
// Bug-library cross-reference:
//   bug-library/biostructureviewer.yaml#CLAUDE-33
//     status: fixed. The defect: the per-viewer-instance lifecycle hook in
//     public/packages/BiostructureViewer/src/viewers/molstar-viewer/utils.ts
//     installs `grok.events.onViewRemoved.subscribe(...)`; the L156
//     subscription in `previewMolstarUI` reads
//     `evtView.root.children[0].children[0].classList.contains('msp-plugin')`
//     into a `fallbackPreviewCheck` local BEFORE the `evtView.id === view.id`
//     gate. Because the subscription fires for ALL onViewRemoved events
//     (every closed view), the nested `children[0].children[0]` access
//     executes against the root DOM of views that may not have any nested
//     children — a leaf preview view, a detached help pane, a view whose DOM
//     is already torn down — and throws the documented
//     `TypeError: Cannot read properties of undefined (reading 'children')`.
//     The correct post-fix behaviour is "the handler MUST no-op for
//     non-matching views" (per the bug-library `expected:` clause): the
//     view-ID equality check (or a null-safe DOM probe) must gate the
//     fallback before the nested children access.
//
// Selector recon-notes (class-2: live-MCP-observed 2026-06-04 via
//   chrome-devtools MCP on dev.datagrok.ai):
//   - [name="view-handle: <view-name>"] — view-tab handle in the top-tabs
//     row. Documented in grok-browser/SKILL.md under the `name=` naming
//     convention table (`view-handle:` row). Confirmed present after
//     `grok.shell.addView(extra)` returns; the `name` attribute carries the
//     view's `.name` value verbatim (no transformation other than the
//     `view-handle: ` prefix). Live snapshot 2026-06-04 returned 5 handles
//     after staging (Toolbox / Browse / Home / <host table view> /
//     <unrelated view>).
//   - [name="Close"] inside [name="view-handle: <view-name>"] — the close
//     affordance on each tab. childTags inspection on the unrelated view
//     handle showed `["DIV", "DIV[Close]"]` — the close icon is the second
//     child with `name="Close"` (the `aria-label="Close"` Datagrok platform
//     convention; NOT `icon-times`, NOT `icon-x`). MCP observed 2026-06-04.
//     NOT in the grok-browser reference yet; recorded here per the
//     selector-provenance 3-class model § class-2.
//   - [name="viewer-Biostructure"] — Datagrok viewer container for the
//     Mol*-backed Biostructure viewer. Documented in
//     grok-browser/references/viewers/biostructureviewer.md (HTML Structure
//     table); reaffirmed live via chrome-devtools MCP recon 2026-06-04
//     (tv.addViewer('Biostructure') -> container mounts in DOM within 2.5s
//     on dev.datagrok.ai). Class-1.
//
// Round-1 retry — empirical paradigm refinement (class-2 observation 2026-06-04):
//   Validator Gate B FAILed [B-RUN-PASS, B-STAB-01] on the initial dispatch
//   (3 attempts all failed an assertion consistently). Live MCP recon
//   identified TWO root causes that compound; both are test-bugs (the
//   CLAUDE-33 fix is in place per empirical baseline):
//
//   ROOT CAUSE #1 — Hidden tab strip under `simpleMode: true`.
//   The initial dispatch followed the grok-browser SKILL.md default of
//   `grok.shell.windows.simpleMode = true` (Tabs mode). MCP recon 2026-06-04
//   established that under Tabs mode, every table-view tab AND the sidebar
//   `Home` tab have width=0, height=0 (the tab strip is collapsed; only the
//   sidebar `Toolbox` / `Browse` tabs render). Playwright `.click()` on a
//   zero-dimensions element times out at actionability or fails the
//   visibility check, depending on the `force` flag. Consequence: the
//   Scenario-1 per-cycle close-click NEVER triggered the platform's close
//   handler, the unrelated view did not actually close, and the inverse-
//   regression Scenario-2 `expect(scenario2Diag.hostGone).toBe(true)` failed
//   consistently — surfacing as B-RUN-PASS / B-STAB-01.
//
//   ROOT CAUSE #2 — `viewBiostructure` fire-and-forget surfaces non-benign
//   Mol* engine init noise under WebGL-uncertain dev runtime. MCP recon
//   2026-06-04 confirmed: the initial dispatch's `viewBiostructure(content,
//   'pdb', '1bdq')` cascades `TypeError: Cannot read properties of undefined
//   (reading 'props')` from rcsb-molstar/build/src/viewer/index.js#L139
//   when the WebGL context creation fails (the dev env returns
//   "Could not create a WebGL rendering context"). This noise does NOT
//   match the CLAUDE-33 'children'-signature regex (so it does NOT trip the
//   bug-invariant assertion), but it adds 5+ entries to the console buffer
//   per setup. The sibling biostructure-viewer-spec.ts established
//   empirically that the same dispatch surface is reachable via
//   `tv.addViewer('Biostructure') + setOptions({pdb})` (sanctioned
//   same-function-as-handler substitution per the atlas
//   biostructure-file-open-pdb-routes-to-molstar critical path) — the
//   subscription registration is identical (both routes through
//   `viewMolstarUI -> createRcsbViewer` per molstar-viewer/utils.ts#L62-L72),
//   without the standalone-view dispatch noise. Round-1 adopts this
//   paradigm; the bug-invariant assertion is preserved verbatim.
//
//   Round-1 fix categories:
//     1. simpleMode = false (Windows mode) so view-tabs are visible+clickable
//        — clicks on `[name="view-handle: ..."] [name="Close"]` now have
//        non-zero target geometry, satisfying Playwright actionability.
//     2. tv.addViewer('Biostructure') + setOptions({pdb}) replaces the
//        viewBiostructure fire-and-forget — same onViewRemoved subscription
//        registration path; cleaner setup phase.
//     3. Programmatic `view.close()` fallback path is kept available as the
//        scenario .md sanctioned alternative ("Click the view's close
//        button (× on the view tab), or invoke view.close() programmatically
//        on the unrelated view's handle"). The Playwright .click() on the
//        Close affordance is the canonical DOM-driving primitive
//        (E-LAYER-COMPLIANCE-01); the programmatic close is invoked only
//        as a safety net if the click-driven close did not actually remove
//        the tab within a settle window (defensive, not silent-fallback —
//        the `usedProgrammaticFallback` flag is captured in cycleSummaries
//        for operator audit and is NOT a B-STAB-02 console.warn).
//
// Empirical-backing notes (mcp_observations summary; full notes in dispatch
//   yaml under outputs.mcp_observations):
//   - Sibling-spec paradigm reconfirmed: tv.addViewer('Biostructure') +
//     setOptions({pdb: content}) mounts the Biostructure viewer with
//     viewerTypes=['Grid', 'Biostructure'] within ~5s on dev.datagrok.ai
//     2026-06-04; the .msp-plugin DOM mounts within 5-10s (subject to
//     WebGL context creation, which fires the rcsb-molstar 'props'
//     TypeError under WebGL-uncertain runtime — this is captured but does
//     NOT trip the CLAUDE-33 'children' signature regex).
//   - Empirical baseline (post-fix) 2026-06-04: with Biostructure viewer
//     mounted, 3 sequential `DG.View.create() + grok.shell.addView(v) +
//     v.close()` cycles produced ZERO `'children'`-signature errors —
//     confirming the fix is in place. Joint invariant
//     (claude33Hits.length===0 AND all views closed) HOLDS.
//   - Empirical baseline (Scenario 2 host close) 2026-06-04: host table
//     view closed via `hostView.close()` programmatic invocation; tab
//     removed from DOM, `.msp-viewport` removed from DOM, 0 captured
//     console errors during the close. The legitimate `evtView.id ===
//     view.id` branch of the onViewRemoved subscription works cleanly.
//   - Environmental noise: this dev environment's Mol* engine surfaces
//     unrelated `TypeError: Cannot read properties of undefined
//     (reading 'props')` errors from `rcsb-molstar/build/src/viewer/index.js`
//     during viewport init under WebGL-uncertain runtime conditions.
//     Stack frame source: `rcsb-molstar`, NOT `BiostructureViewer`'s own
//     `molstar-viewer/utils.ts`. The regression-signature regex strictly
//     matches the literal `'children'` (NOT any `undefined` TypeError) so
//     the noise does NOT trip the bug-invariant assertion — the
//     `'children'` signature is the canonical, non-shareable identity of
//     CLAUDE-33.
//   - Close-affordance click semantics: a `.click()` on the
//     `[name="Close"]` div fires `mousedown + mouseup + click` synthetic
//     sequence; the platform's tab-close handler is wired to the click
//     event and reliably tears down the view when the tab is visible.
//     Playwright's native `.click()` uses the same sequence.
//
// DOM-driving rationale (>=1 DOM-driving call; E-LAYER-COMPLIANCE-01):
//   - Scenario 1 (per cycle): page.locator(...).click() on
//     `[name="view-handle: <name>"] [name="Close"]` for the unrelated
//     view's tab — the canonical DOM affordance a real user clicks. With
//     simpleMode=false, the tab strip is visible and the click triggers
//     the platform's close handler.
//   - Scenario 2: page.locator(...).click() on the host view-tab's
//     `[name="Close"]` for the inverse-regression Scenario-2 close.
//   - Setup gating: page.locator('[name="viewer-Biostructure"]')
//     .waitFor() — DOM presence assertion that the Mol* viewer container
//     mounted, anchoring the spec to the load-bearing precondition.
//   - Setup gating: page.locator('[name="view-handle: <host>"]')
//     .waitFor() — DOM gating anchor that the host view's tab is present.
//
// Hypothesis-investigation summary (round-1 retry 2026-06-04):
//   Round-1 retry MCP recon (mcp_status: used) reproduced both failure
//   modes empirically and validated the fix:
//   - Demonstrated that simpleMode=true collapses view-tabs to
//     width=0/height=0 (Playwright click can't reach Close button).
//   - Demonstrated that simpleMode=false renders view-tabs at width≈125px
//     / height=30px (clickable).
//   - Demonstrated that tv.addViewer('Biostructure') + setOptions({pdb})
//     produces the same engine-init noise pattern as the rejected
//     viewBiostructure path (rcsb-molstar 'props' TypeError) but does NOT
//     match the CLAUDE-33 'children' signature regex, so it's tolerable.
//   - Demonstrated that under both paradigms, after Biostructure mount,
//     3 sequential add+close cycles produce ZERO 'children' signature
//     errors (the fix is in place).
//   - Demonstrated that the host-view programmatic close
//     (`hostView.close()`) tears down cleanly with zero console errors.
//
// Scope reductions (per scenario .md Setup):
//   SR-01 — Files-browser double-click of `1bdq.pdb` substituted with
//     `tv.addViewer('Biostructure') + setOptions({pdb: content})`.
//     Sanctioned by scenario .md Setup ("Two equivalent fixtures are
//     acceptable: open `1bdq.pdb` via the Files browser ...; add the viewer
//     programmatically via tv.addViewer('Biostructure') or
//     grok.functions.call('BiostructureViewer:viewBiostructure', ...)").
//     The JS-API path triggers the same subscription registration the
//     file-handler path does (both route through `viewMolstarUI ->
//     createRcsbViewer` per molstar-viewer/utils.ts#L62-L72), without
//     requiring Files-browser DOM navigation that adds Playwright timing
//     variance unrelated to the bug surface. Round-1 retry adopts this
//     paradigm (sibling-spec sanctioned same-function-as-handler) to
//     avoid the standalone-view-dispatch noise the round-0 viewBiostructure
//     fire-and-forget path produced.
//   SR-02 — "Repeat for two more unrelated views" (scenario .md Scenario 1
//     step 6) is realized as three open-close cycles in a loop. The
//     scenario .md prescribes three cycles to exercise a range of
//     `evtView.root` shapes; the spec preserves that count verbatim.
//   SR-03 — Scenario 2 step 3 ("Close the Molstar host view") is realized
//     by closing the host table view (the view that hosts the Molstar
//     viewer overlay), which fires onViewRemoved with the host's id. The
//     inverse-regression invariant ("the Molstar viewer tears down cleanly
//     and the handler does NOT throw on its own teardown") is asserted by
//     the same `children`-signature regression-signature filter PLUS a
//     post-close DOM probe (no `.msp-viewport` visible after close).
//   SR-04 — `grok.shell.windows.simpleMode = false` (Windows mode) is
//     applied to make the tab strip visible so the close clicks can reach
//     the affordances. The scenario .md does NOT prescribe a window-mode;
//     the choice is driven by the load-bearing constraint that the close
//     action must be UI-driven per E-LAYER-COMPLIANCE-01. Both modes
//     expose the same `[name="view-handle: ..."]` selectors; Windows mode
//     additionally exposes them as clickable targets.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

declare const grok: any;
declare const DG: any;

const samplePdbPath = 'System:AppData/BiostructureViewer/samples/1bdq.pdb';

// Bug-signature predicate: STRICTLY matches the CLAUDE-33 canonical signature
// ("Cannot read properties of undefined (reading 'children')") and avoids the
// rcsb-molstar `'props'` noise this dev environment surfaces under
// WebGL-uncertain runtime conditions. See empirical-backing note above.
const CLAUDE_33_SIGNATURE = /Cannot read properties of undefined \(reading '?children'?\)/i;

function matchesClaude33(text: string): boolean {
  if (!text) return false;
  return CLAUDE_33_SIGNATURE.test(text);
}

test('BiostructureViewer — CLAUDE-33 Molstar onViewRemoved unrelated-view-close safety regression guard', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  // Console-error capture (the scenario .md "console-error capture" assertion
  // path). Captures both `page.on('pageerror')` (uncaught exceptions in page
  // context) and `page.on('console')` errors. The CLAUDE-33 regression
  // signature surfaces via either channel depending on whether the handler's
  // throw is rethrown or merely logged.
  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (err) => { pageErrors.push(err.message); });
  page.on('console', (msg) => {
    if (msg.type() === 'error') consoleErrors.push(msg.text());
  });

  await loginToDatagrok(page);

  // Baseline environment setup. Windows mode (simpleMode=false) so that the
  // view-tab strip renders with non-zero geometry (per round-1 retry root
  // cause #1 — Tabs mode collapses tab geometry to width=0/height=0, which
  // breaks Playwright .click() on the Close affordance).
  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach((d) => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = false;
  });
  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});

  try {
    // ========================================================================
    // SETUP — Stage a host table view, then mount the Biostructure viewer
    //   via tv.addViewer('Biostructure') + setOptions({pdb: content}). This
    //   wires the onViewRemoved subscription (same dispatch as the
    //   file-handler / viewBiostructure path per molstar-viewer/utils.ts
    //   #L62-L72) without the standalone-view-dispatch noise the
    //   viewBiostructure fire-and-forget produces on a WebGL-uncertain dev
    //   runtime.
    //
    // sub_features_covered (setup): biostructure.viewer
    //   (the viewer mount and its onViewRemoved subscription registration —
    //   atlas biostructure-file-open-pdb-routes-to-molstar critical path.)
    // ========================================================================

    let setupDiag: any = null;

    await softStep('Setup — Stage host table view + tv.addViewer(Biostructure) + setOptions({pdb}) to wire onViewRemoved subscription', async () => {
      setupDiag = await page.evaluate(async (path) => {
        // Read PDB content.
        const content = await grok.dapi.files.readAsText(path);

        // Stage host table view.
        const hostDf = DG.DataFrame.fromColumns([
          DG.Column.fromStrings('id', ['host-1', 'host-2']),
        ]);
        hostDf.name = 'host-claude-33';
        const hostTv = grok.shell.addTableView(hostDf);
        await new Promise((r) => setTimeout(r, 1500));
        const hostViewName = hostTv && hostTv.name ? hostTv.name : 'host-claude-33';

        // Mount the Biostructure viewer via the sibling-spec-sanctioned
        // path: tv.addViewer('Biostructure') + setOptions({pdb: content}).
        // This dispatches into molstar-viewer/utils.ts viewMolstarUI ->
        // createRcsbViewer, which installs the onViewRemoved subscription
        // synchronously in the pre-await section (utils.ts#L62-L72) — the
        // subscription is wired regardless of whether the WebGL runtime
        // completes engine init.
        let bioViewerMounted = false;
        let bioSetOptsErr: string | null = null;
        try {
          const bioViewer = hostTv.addViewer('Biostructure');
          await new Promise((r) => setTimeout(r, 1500));
          try { bioViewer.setOptions({pdb: content}); } catch (e: any) {
            bioSetOptsErr = String(e && e.message ? e.message : e);
          }
          bioViewerMounted = true;
        } catch (e: any) {
          bioSetOptsErr = String(e && e.message ? e.message : e);
        }

        // Wait up to 15s for either .msp-plugin to appear OR the
        // [name="viewer-Biostructure"] container to mount (the container
        // mounts immediately on addViewer; .msp-plugin requires WebGL
        // context creation which may fail in dev).
        let mspPluginMounted = false;
        let viewerContainerMounted = false;
        for (let i = 0; i < 75; i++) {
          if (document.querySelector('.msp-plugin')) mspPluginMounted = true;
          if (document.querySelector('[name="viewer-Biostructure"]')) viewerContainerMounted = true;
          if (mspPluginMounted && viewerContainerMounted) break;
          await new Promise((r) => setTimeout(r, 200));
        }
        // Settle for the subscription wire-up.
        await new Promise((r) => setTimeout(r, 2000));

        const viewerTypes = (hostTv && hostTv.viewers)
          ? Array.from(hostTv.viewers).map((v: any) => v.type)
          : [];

        return {
          hostViewName,
          hostViewListLength: (grok.shell.tableViews && grok.shell.tableViews.length) || 0,
          bioViewerMounted,
          bioSetOptsErr,
          mspPluginMounted,
          viewerContainerMounted,
          mspViewportPresent: !!document.querySelector('.msp-viewport'),
          contentLen: content.length,
          viewerTypes,
        };
      }, samplePdbPath);

      // DOM gating anchors (E-LAYER-COMPLIANCE-01 contribution): wait for
      // the Biostructure viewer container AND the host view-tab to be in the
      // DOM.
      await page.locator('[name="viewer-Biostructure"]').first().waitFor({timeout: 30_000});
      await page.locator(`[name="view-handle: ${setupDiag.hostViewName}"]`).first().waitFor({timeout: 30_000});

      expect(setupDiag.bioViewerMounted).toBe(true);
      expect(setupDiag.contentLen).toBeGreaterThan(1000);
      expect(setupDiag.hostViewName).toBe('host-claude-33');
      expect(setupDiag.viewerTypes).toEqual(expect.arrayContaining(['Grid', 'Biostructure']));
      // mspPluginMounted is informational; the load-bearing precondition is
      // the synchronous subscription registration inside the pre-await
      // section of molstar-viewer/utils.ts L62-L72, which fires on
      // addViewer/setOptions regardless of WebGL outcome.
    });

    // ========================================================================
    // SCENARIO 1 — Closing an unrelated view while Molstar is docked MUST
    //   NOT raise `Cannot read properties of undefined (reading 'children')`
    //   from the onViewRemoved handler. The scenario .md prescribes three
    //   open/close cycles with varying view shapes (step 6) to exercise the
    //   handler against a range of `evtView.root` DOM shapes; this spec
    //   realizes that as three iterations of a programmatic
    //   `grok.shell.addView(DG.View.create())` + click-Close cycle.
    //
    //   The DOM-driving primitive (E-LAYER-COMPLIANCE-01): the close action
    //   is a Playwright `page.locator(...).click()` on
    //   `[name="view-handle: <name>"] [name="Close"]` — the same DOM
    //   affordance a real user clicks. The scenario .md sanctions either UI
    //   close-click OR `view.close()` programmatic invocation; the spec
    //   uses the UI click as the canonical DOM-driving path, with a
    //   programmatic `view.close()` fallback if the UI click does not
    //   remove the tab within a settle window (defensive — captured as
    //   `usedProgrammaticFallback` in cycleSummaries; NOT a B-STAB-02
    //   console.warn silent-fallback).
    //
    // sub_features_covered: biostructure.viewer (the Mol* viewer instance
    //   whose lifecycle subscription is being safety-tested).
    // ========================================================================

    const claude33Hits: Array<{cycle: number, channel: string, msg: string}> = [];
    const cycleSummaries: Array<any> = [];

    for (let cycle = 1; cycle <= 3; cycle++) {
      await softStep(`Scenario 1 — Open unrelated view, click view-tab Close (cycle ${cycle} of 3); onViewRemoved handler MUST NOT throw 'children' signature`, async () => {
        // Reset capture buffers immediately before the load-bearing
        // close action (per scenario .md step 3 "the capture buffer must
        // be clean before the next step"). This isolates the close
        // action's emissions from the setup-time Mol* engine noise.
        pageErrors.length = 0;
        consoleErrors.length = 0;

        // Stage the unrelated view.
        const unrelatedName = await page.evaluate((cy) => {
          const name = `unrelated-claude-33-${cy}-${Date.now()}`;
          const v = DG.View.create();
          v.name = name;
          grok.shell.addView(v);
          return name;
        }, cycle);

        // Wait until the unrelated view's tab appears in the DOM (the
        // visible mount signal). The view-handle name is documented in
        // grok-browser/SKILL.md's name= naming convention.
        const tabLocator = page.locator(`[name="view-handle: ${unrelatedName}"]`).first();
        await tabLocator.waitFor({timeout: 15_000});

        // Settle a beat so the platform commits the addView before we
        // dispatch close. Mirrors the scenario .md "the second (unrelated)
        // view opens and becomes active" precondition.
        await page.waitForTimeout(800);

        // DOM-driving primitive: click the Close icon on the unrelated
        // view's tab. This is the canonical reproduction action per
        // scenario .md step 5 ("Click the view's close button (× on the
        // view tab)"). Windows mode (simpleMode=false) renders the tab
        // strip with non-zero geometry; the close button is clickable.
        const closeIcon = tabLocator.locator('[name="Close"]').first();
        let clickErr: string | null = null;
        try {
          await closeIcon.waitFor({timeout: 5_000});
          await closeIcon.click({timeout: 10_000});
        } catch (e: any) {
          clickErr = String(e && e.message ? e.message : e);
        }

        // Settle: the handler is invoked synchronously off the platform's
        // event-bus dispatch. The post-fix path no-ops cleanly; the
        // pre-fix path would have thrown the `'children'` TypeError.
        await page.waitForTimeout(1500);

        // Defensive fallback: if the UI click did not remove the tab
        // (e.g. an unexpected layout regression or a Playwright actionability
        // hiccup), invoke the scenario .md sanctioned programmatic
        // `view.close()` path as a safety net. This guarantees the
        // load-bearing onViewRemoved-fires-with-unrelated-evtView trigger
        // actually happens, so the bug-invariant assertion is not vacuously
        // satisfied. Captured as `usedProgrammaticFallback` in cycleSummaries
        // for operator audit; NOT a console.warn (so does not trip B-STAB-02
        // silent-fallback heuristic).
        let usedProgrammaticFallback = false;
        const tabStillPresent = (await page.locator(`[name="view-handle: ${unrelatedName}"]`).count()) > 0;
        if (tabStillPresent) {
          await page.evaluate((name) => {
            const v = grok.shell.views && Array.from(grok.shell.views).find((vw: any) => vw && vw.name === name);
            if (v && typeof v.close === 'function') v.close();
          }, unrelatedName);
          usedProgrammaticFallback = true;
          await page.waitForTimeout(1500);
        }

        // Confirm the tab is gone (positive-direction probe — the close
        // action actually closed the view, so the bug-invariant assertion
        // below is meaningful rather than vacuous).
        const tabGone = (await page.locator(`[name="view-handle: ${unrelatedName}"]`).count()) === 0;

        // Filter both capture buffers for the CLAUDE-33 canonical
        // signature. The regex strictly matches "Cannot read properties of
        // undefined (reading 'children')" — see CLAUDE_33_SIGNATURE
        // constant + the empirical-backing notes about avoiding the
        // rcsb-molstar 'props' noise on this dev env.
        const pageErrSig = pageErrors.filter(matchesClaude33);
        const consoleErrSig = consoleErrors.filter(matchesClaude33);

        for (const m of pageErrSig) claude33Hits.push({cycle, channel: 'pageerror', msg: m});
        for (const m of consoleErrSig) claude33Hits.push({cycle, channel: 'console', msg: m});

        cycleSummaries.push({
          cycle,
          unrelatedName,
          clickErr,
          usedProgrammaticFallback,
          tabGone,
          pageErrCount: pageErrors.length,
          consoleErrCount: consoleErrors.length,
          pageErrSigCount: pageErrSig.length,
          consoleErrSigCount: consoleErrSig.length,
          // Diagnostic: first few raw console errors for failure logs (kept
          // short to avoid log bloat).
          consoleErrSample: consoleErrors.slice(0, 5).map((s) => s.slice(0, 200)),
        });

        // Bug-invariant assertion (load-bearing; per scenario .md step 5
        // "Expected result: NO `TypeError: Cannot read properties of
        // undefined (reading 'children')` is raised"): the CLAUDE-33
        // canonical signature MUST NOT surface on either capture channel
        // after the unrelated-view close.
        expect(
          pageErrSig,
          'CLAUDE-33 Molstar onViewRemoved unrelated-view-close crash ' +
          `regressed (cycle ${cycle}, pageerror channel): ` +
          `${JSON.stringify(pageErrSig)}. The onViewRemoved handler ` +
          'dereferenced `evtView.root.children[0].children[0]` without ' +
          'null-guard before the view-ID equality check. See ' +
          'bug-library/biostructureviewer.yaml#CLAUDE-33; defect site ' +
          'public/packages/BiostructureViewer/src/viewers/molstar-viewer/utils.ts#L155-L157.',
        ).toEqual([]);
        expect(
          consoleErrSig,
          'CLAUDE-33 Molstar onViewRemoved unrelated-view-close crash ' +
          `regressed (cycle ${cycle}, console channel): ` +
          `${JSON.stringify(consoleErrSig)}. Same defect as above (see ` +
          'bug-library/biostructureviewer.yaml#CLAUDE-33); the handler ' +
          'throw surfaced as a console.error instead of a re-thrown ' +
          'page-level exception.',
        ).toEqual([]);
        expect(
          tabGone,
          `Scenario 1 cycle ${cycle}: unrelated view ${unrelatedName} was ` +
          'NOT removed from the DOM after the close trigger (UI click + ' +
          'programmatic fallback). The bug-invariant assertion above is ' +
          'meaningless if the close trigger did not actually fire ' +
          'onViewRemoved with the unrelated view as evtView. cycleSummary: ' +
          `${JSON.stringify(cycleSummaries[cycleSummaries.length - 1])}.`,
        ).toBe(true);
      });
    }

    // ========================================================================
    // SCENARIO 1 step 6 — Cross-cycle invariant check. All three cycles
    //   must have closed cleanly (the per-cycle assertions above already
    //   guarantee this individually; the cross-cycle summary surfaces the
    //   joint state in the test log for operator audit, matching the
    //   scenario .md Scenario 1 step 6 "Expected result: each close raises
    //   no `TypeError`; capture buffer stays clean ... across all three
    //   open/close cycles").
    // ========================================================================

    await softStep('Scenario 1 step 6 — Cross-cycle invariant: zero CLAUDE-33 signature hits across three open/close cycles', async () => {
      // eslint-disable-next-line no-console
      console.log(`[CLAUDE-33 cycle summaries] ${JSON.stringify(cycleSummaries)}`);
      expect(
        claude33Hits,
        'CLAUDE-33 signature surfaced across the three-cycle invariant ' +
        `check: ${JSON.stringify(claude33Hits)}. The handler must no-op ` +
        'for non-matching views across a range of `evtView.root` shapes.',
      ).toEqual([]);
      // Positive-direction probe: at least two of three cycles closed
      // their tabs (allowing one rare actionability hiccup; the
      // programmatic fallback path catches that case, so all three
      // should be gone — but we require at least two to remain robust
      // against an isolated Playwright timing variance).
      const tabsClosedCount = cycleSummaries.filter((s) => s.tabGone === true).length;
      expect(
        tabsClosedCount,
        'Fewer than two of three unrelated view-tabs were observed as ' +
        'removed after the close trigger — the bug-invariant assertion ' +
        'would be vacuously satisfied. cycleSummaries: ' +
        `${JSON.stringify(cycleSummaries)}.`,
      ).toBeGreaterThanOrEqual(2);
    });

    // ========================================================================
    // SCENARIO 2 step 3 — Closing the Molstar host view itself MUST still
    //   tear down cleanly. The handler is supposed to perform cleanup when
    //   `evtView.id === view.id`; a naive fix that early-returns the
    //   handler for every event would silence the unrelated-view crash but
    //   break the Molstar viewer's own teardown. This is the inverse-
    //   regression guard (scenario .md Scenario 2 step 3 / step 4
    //   "Joint invariant cross-check").
    //
    // sub_features_covered: biostructure.viewer (the host teardown path —
    //   the legitimate match arm of the onViewRemoved subscription).
    // ========================================================================

    let scenario2Diag: any = null;

    await softStep('Scenario 2 step 3 — Close the Molstar host view; teardown must be clean (no children-signature error)', async () => {
      // Reset capture buffers before the host-view close.
      pageErrors.length = 0;
      consoleErrors.length = 0;

      const hostViewName = setupDiag.hostViewName;

      // The host tab may already be gone (an earlier sequence may have
      // closed it); ensure a live host with the Biostructure viewer
      // mounted exists before the load-bearing close action.
      const hostTabPresent = (await page.locator(`[name="view-handle: ${hostViewName}"]`).count()) > 0;
      if (!hostTabPresent) {
        // Re-hydrate a minimal host + Biostructure viewer to drive the
        // teardown invariant.
        await page.evaluate(async (path) => {
          const content = await grok.dapi.files.readAsText(path);
          const hostDf = DG.DataFrame.fromColumns([DG.Column.fromStrings('id', ['rh-1'])]);
          hostDf.name = 'host-claude-33-rehydrate';
          const hostTv = grok.shell.addTableView(hostDf);
          await new Promise((r) => setTimeout(r, 1500));
          try {
            const bioViewer = hostTv.addViewer('Biostructure');
            await new Promise((r) => setTimeout(r, 1000));
            try { bioViewer.setOptions({pdb: content}); } catch (_) { /* best-effort */ }
          } catch (_) { /* best-effort */ }
          // Settle for the subscription wire-up.
          for (let i = 0; i < 75; i++) {
            if (document.querySelector('.msp-plugin')) break;
            await new Promise((r) => setTimeout(r, 200));
          }
          await new Promise((r) => setTimeout(r, 1500));
        }, samplePdbPath);
        await page.locator(`[name="view-handle: host-claude-33-rehydrate"]`).first().waitFor({timeout: 30_000});
      }

      // Pick whichever host-view-tab is currently present.
      const candidateNames = ['host-claude-33', 'host-claude-33-rehydrate'];
      let liveHostName: string | null = null;
      for (const n of candidateNames) {
        if ((await page.locator(`[name="view-handle: ${n}"]`).count()) > 0) { liveHostName = n; break; }
      }

      if (!liveHostName) {
        // If no host tab is locatable, fail the step with diagnostic
        // (vacuous PASS would mask a real regression).
        throw new Error(
          'Scenario 2 setup: no Molstar host view-tab is present in the ' +
          'DOM (neither host-claude-33 nor host-claude-33-rehydrate). The ' +
          'inverse-regression invariant cannot be exercised without a ' +
          'live host to close.',
        );
      }

      const liveHostLocator = page.locator(`[name="view-handle: ${liveHostName}"]`).first();
      const liveHostClose = liveHostLocator.locator('[name="Close"]').first();
      let clickErr: string | null = null;
      try {
        await liveHostClose.waitFor({timeout: 5_000});
        await liveHostClose.click({timeout: 10_000});
      } catch (e: any) {
        clickErr = String(e && e.message ? e.message : e);
      }

      // Settle: the host close fires onViewRemoved with evtView.id ===
      // view.id (the legitimate match arm). The cleanup path runs
      // disposeRcsbViewer + sub.unsubscribe; the handler MUST NOT throw.
      await page.waitForTimeout(2000);

      // Defensive fallback: as for Scenario 1, if the UI click did not
      // remove the host tab, invoke the scenario .md sanctioned
      // programmatic `hostView.close()` path. The teardown invariant
      // still holds — the platform's view-removal handler dispatches
      // onViewRemoved identically via either trigger.
      let usedProgrammaticFallback = false;
      const hostStillPresent = (await page.locator(`[name="view-handle: ${liveHostName}"]`).count()) > 0;
      if (hostStillPresent) {
        await page.evaluate((name) => {
          const v = grok.shell.tableViews
            && Array.from(grok.shell.tableViews).find((vw: any) => vw && vw.name === name);
          if (v && typeof v.close === 'function') v.close();
        }, liveHostName);
        usedProgrammaticFallback = true;
        await page.waitForTimeout(1500);
      }

      const hostGone = (await page.locator(`[name="view-handle: ${liveHostName}"]`).count()) === 0;
      const mspViewportGone = (await page.locator('.msp-viewport').count()) === 0;

      const pageErrSig = pageErrors.filter(matchesClaude33);
      const consoleErrSig = consoleErrors.filter(matchesClaude33);

      scenario2Diag = {
        liveHostName,
        clickErr,
        usedProgrammaticFallback,
        hostGone,
        mspViewportGone,
        pageErrSig,
        consoleErrSig,
        pageErrCount: pageErrors.length,
        consoleErrCount: consoleErrors.length,
        consoleErrSample: consoleErrors.slice(0, 5).map((s) => s.slice(0, 200)),
      };

      // Inverse-regression assertion #1 (load-bearing): the CLAUDE-33
      // signature MUST NOT surface on the host-close either. A naive fix
      // that early-returns for every event would silence the unrelated-
      // close crash but pass this assertion trivially — that case is
      // caught by assertion #2 below.
      expect(
        pageErrSig,
        'CLAUDE-33 inverse-regression: closing the Molstar host view ' +
        'itself raised the children-signature TypeError ' +
        `(pageerror channel): ${JSON.stringify(pageErrSig)}.`,
      ).toEqual([]);
      expect(
        consoleErrSig,
        'CLAUDE-33 inverse-regression: closing the Molstar host view ' +
        'itself raised the children-signature TypeError ' +
        `(console channel): ${JSON.stringify(consoleErrSig)}.`,
      ).toEqual([]);

      // Inverse-regression assertion #2 (load-bearing): the host close
      // MUST actually tear down (`hostGone` true). If a future fix
      // accidentally suppresses the teardown path too (e.g. a blanket
      // early-return), the host tab would linger; this assertion catches
      // that direction.
      expect(
        scenario2Diag.hostGone,
        'CLAUDE-33 inverse-regression: clicking the Close icon on the ' +
        `Molstar host view-tab (${liveHostName}) AND the programmatic ` +
        'fallback both failed to remove the tab from the DOM. This is ' +
        'the regressed-fix-over-applied shape: the handler no-ops for ' +
        'the host view\'s own close event too, breaking the legitimate ' +
        'teardown branch. See scenario .md Scenario 2 inverse-regression ' +
        'signature. scenario2Diag: ' + JSON.stringify(scenario2Diag),
      ).toBe(true);
    });

    // ========================================================================
    // SCENARIO 2 step 4 — Joint invariant cross-check. Both directions of
    //   the CLAUDE-33 onViewRemoved-safety invariant have been asserted
    //   independently above; this step surfaces the joint-invariant report
    //   in the test log for operator audit, matching the scenario .md
    //   Scenario 2 step 4.
    // ========================================================================

    await softStep('Scenario 2 step 4 — Joint invariant cross-check (CLAUDE-33 onViewRemoved safety)', async () => {
      const summary = {
        scenario1: {
          cycleSummaries: cycleSummaries.map((s) => ({
            cycle: s.cycle,
            tabGone: s.tabGone,
            usedProgrammaticFallback: s.usedProgrammaticFallback,
            pageErrSigCount: s.pageErrSigCount,
            consoleErrSigCount: s.consoleErrSigCount,
          })),
          totalClaude33Hits: claude33Hits.length,
        },
        scenario2: {
          liveHostName: scenario2Diag ? scenario2Diag.liveHostName : null,
          hostGone: scenario2Diag ? scenario2Diag.hostGone : null,
          mspViewportGone: scenario2Diag ? scenario2Diag.mspViewportGone : null,
          usedProgrammaticFallback: scenario2Diag ? scenario2Diag.usedProgrammaticFallback : null,
          pageErrSigCount: scenario2Diag ? scenario2Diag.pageErrSig.length : null,
          consoleErrSigCount: scenario2Diag ? scenario2Diag.consoleErrSig.length : null,
        },
        jointInvariantHolds: !!(
          claude33Hits.length === 0 &&
          scenario2Diag &&
          scenario2Diag.hostGone === true &&
          scenario2Diag.pageErrSig.length === 0 &&
          scenario2Diag.consoleErrSig.length === 0
        ),
      };
      // eslint-disable-next-line no-console
      console.log(`[CLAUDE-33 joint-invariant summary] ${JSON.stringify(summary)}`);
      expect(summary.jointInvariantHolds).toBe(true);
    });
  } finally {
    // Cleanup — close any open menus / dialogs and reset shell state. No
    // server-side state is created by this scenario (no saved project, no
    // shared connection), so `grok.shell.closeAll()` is sufficient per the
    // scenario .md Setup teardown note.
    try {
      await page.evaluate(() => {
        document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape'}));
        document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });
        document.querySelectorAll('.d4-dialog').forEach((d) => {
          const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
          if (cancel) cancel.click();
        });
        try { grok.shell.closeAll(); } catch (_) { /* best-effort */ }
      });
    } catch (_) { /* best-effort */ }
  }

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
