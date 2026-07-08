---
feature: powerpack
target_layer: playwright
coverage_type: regression
priority: p1
realizes: [direct-link-loading-window]
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Powerpack/direct-link-loading.md
migration_date: 2026-05-23
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
related_bugs:
  - GROK-18721
realized_as:
  - direct-link-loading-spec.ts
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-05-23-powerpack-migrate-02
    timestamp: 2026-05-23T00:00:00Z
    review_round: 1
    failure_keys: []
  d:
    verdict: PASS
    cycle_id: 2026-05-23-powerpack-migrate-02
    timestamp: 2026-05-23T00:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T11:05:01Z
    spec_runs:
      - spec: direct-link-loading-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 52
        failure_keys: []
---

# PowerPack — Direct-link entry path loading window rendering (GROK-18721 regression)

Regression test for GROK-18721: when a Datagrok app or dashboard is
opened via a direct link (URL paste in a fresh browser, bookmark, or
external deeplink), PowerPack's loading window must render fully and
properly — the same visual quality as when navigating to the same app
or dashboard from inside the platform. Before the fix, the direct-link
entry path produced a cropped or incomplete loading-window layout,
because PowerPack's initialization ran before the page's layout
dimensions were established.

## Setup

1. Open Datagrok with PowerPack installed (default platform load).
   Confirm the user has access to the platform.
2. Ensure a saved project or app with a known direct-link URL exists.
   If no suitable project / app is available, create one first:
   - Open `System:DemoFiles/demog.csv` so a Demog table view is the
     active view.
   - `File | Save As Project...` and save with name
     `direct-link-loading-test`.
   - Record the direct-link URL surfaced by Datagrok for this saved
     project (the form is `https://<server>/p/<owner>.<project>`,
     e.g. `https://dev.datagrok.ai/p/<user>.direct-link-loading-test`
     — the same shape as the bug report's
     `https://dev.datagrok.ai/p/Opavlenko.Demog_95`).
3. Have a fresh browser context ready (incognito / private window OR
   a new browser profile with no warm Datagrok session). Direct-link
   entry must be exercised against an uncached / non-warm session so
   the lifecycle.init race condition that produced the bug is on
   the entry path.

## Scenarios

### Scenario 1: Direct-link entry renders PowerPack loading window fully (regression repro)

1. **Open a fresh browser context.** Use an incognito / private
   window or a browser profile with no warm Datagrok session. The
   address bar must not auto-complete the project URL; this ensures
   PowerPack's `powerPackInit` runs from scratch on this navigation.
2. **Paste the direct-link URL.** Paste the URL recorded in Setup
   step 2 (e.g.
   `https://dev.datagrok.ai/p/<user>.direct-link-loading-test`) into
   the address bar and press Enter.
3. **Observe the PowerPack loading window during page load.** While
   the page is still loading (between the initial blank state and
   the fully rendered project view), capture / inspect the
   PowerPack loading window. Verify that:
   - The loading window renders with proper dimensions — width and
     height match the available viewport (or the platform's
     configured loading affordance size), not a degenerate /
     zero-dimension box.
   - The loading window's contents (PowerPack spinner / branding /
     progress affordance, depending on platform skin) are fully
     visible — no cropping along any edge, no clipping of icons or
     text.
   - No incomplete-layout artifacts are visible: no overlapping
     panels, no orphan background fragments, no partially-painted
     welcome-view skeleton beneath the loading window.
4. **Wait for the load to complete.** Allow the page to finish
   loading; the project's table view (the Demog dataset under the
   `direct-link-loading-test` project layout) should become the
   active view.
5. **Verify post-load rendering.** Once load completes:
   - The project's table view renders correctly: table grid is
     visible, columns are labeled, the breadcrumb / view title
     reflects `direct-link-loading-test`.
   - The home / welcome view that PowerPack briefly hosted during
     load has yielded cleanly to the project view (no leftover
     home-view fragments, no zombie loading affordance).
   - PowerPack dashboards / widgets, if visible at any moment of
     the load (the welcome view's home-dashboard host can flash
     during transition), did not render cropped or incomplete.

Expected:
- The PowerPack loading window renders fully during direct-link
  page load — proper dimensions, no cropping, complete layout.
- After load completes, the project's view renders correctly with
  the same visual quality as in-app navigation.
- No visual glitches (cropped loading UI, incomplete welcome
  skeleton, overlapping panels) are observed at any point of the
  direct-link entry path.

### Scenario 2: In-app navigation renders the same loading window state (control case)

This scenario is the regression-guard control. It exercises the
in-app navigation path that — per the bug report — never reproduced
the visual glitch. Comparing direct-link entry (Scenario 1) against
in-app navigation (this Scenario 2) verifies the fix did not
introduce a regression on the previously-working path AND that the
direct-link path now matches the same visual quality.

1. **Start from inside Datagrok.** Open Datagrok in any window
   (warm session is fine; this scenario does not require a fresh
   browser context). Navigate to the home / welcome view.
2. **Open the same project via in-app navigation.** Use one of
   the in-app entry paths the platform supports:
   - From the Recent Projects widget on the home view, click the
     `direct-link-loading-test` card.
   - Or from the Projects browser (`Browse | Projects`), locate
     and open the project.
   - Or from the global search bar, type the project name and
     select the result.
3. **Observe the PowerPack loading window during in-app open.**
   The loading window may briefly surface during the project open;
   capture / inspect its state. Verify the same render quality as
   the expected Scenario 1 post-fix behavior:
   - Proper dimensions, no cropping, complete layout.
4. **Wait for the open to complete.** The project's table view
   becomes the active view.
5. **Compare against Scenario 1 outcome.** Cross-check that the
   final rendered state of the project view matches Scenario 1's
   step 5 outcome — i.e. the direct-link entry path produces the
   same visual quality as in-app navigation.

Expected:
- In-app navigation renders the PowerPack loading window fully and
  properly (preserved control behavior — never broken by the bug).
- The direct-link entry path (Scenario 1) matches this visual
  quality after the fix.
- No regression on the in-app navigation path is introduced by the
  fix.

## Notes

- **Related bug.** GROK-18721 (fixed): copying a direct link to an
  app or dashboard (e.g. `https://dev.datagrok.ai/p/Opavlenko.Demog_95`)
  and opening it in a fresh browser caused PowerPack's loading window
  to render cropped or incomplete during page load. Direct-link
  navigation must render the loading window with the same visual
  quality as in-app navigation.
