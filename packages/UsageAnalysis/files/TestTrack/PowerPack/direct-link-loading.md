---
feature: powerpack
sub_features_covered:
  - powerpack.lifecycle.init
  - powerpack.welcome.view
  - powerpack.dashboards
target_layer: playwright
coverage_type: regression
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
    claims:
      - check: A-STRUCT-MECH-01
        status: PASS
        evidence: |
          Frontmatter delimited by leading and trailing `---` fences
          parses as well-formed YAML. All four required fields are
          present and well-typed. feature is the string powerpack.
          sub_features_covered is a list of three ids — each one is
          resolvable in references/feature-atlas/powerpack.yaml
          (powerpack.lifecycle.init at line 58, powerpack.welcome.view
          at line 87, powerpack.dashboards at line 102). target_layer
          is the string playwright. coverage_type is the string
          regression. The four Phase 1 sidecar-less keys
          (source_text_fixes, candidate_helpers, unresolved_ambiguities,
          scope_reductions) are present as empty inline-flow lists; that
          is not required by A-STRUCT-MECH-01 but is consistent with
          the sidecar-less emission contract.
      - check: A-STRUCT-MECH-02
        status: PASS
        evidence: |
          Body contains three top-level `## ` headings (Setup,
          Scenarios, Notes) and two H3 scenario blocks under
          `## Scenarios`: Scenario 1 direct-link regression repro and
          Scenario 2 in-app navigation control. A-STRUCT-MECH-02 only
          requires that at least one `## ` heading exist in body; that
          is satisfied. The H3 scenario nesting is the established
          TestTrack convention used by sibling Powerpack files.
      - check: A-STRUCT-MECH-03
        status: PASS
        evidence: |
          Both scenario sub-headings carry numbered steps. Scenario 1
          opens with `1. **Open a fresh browser context.**` and runs
          through `5. **Verify post-load rendering.**`. Scenario 2
          opens with `1. **Start from inside Datagrok.**` and runs
          through `5. **Compare against Scenario 1 outcome.**`. No
          scenario heading lacks numbered steps.
      - check: A-STRUCT-MECH-04
        status: PASS
        evidence: |
          Neither scenario is empty. Scenario 1 exercises fresh-browser
          direct-link entry with explicit observation of the loading
          window during page load plus post-load rendering checks.
          Scenario 2 is the in-app navigation regression-guard control
          with explicit comparison to Scenario 1's outcome. Both
          scenarios contain executable verification content far beyond
          a bare heading.
      - check: A-STRUCT-MECH-05
        status: PASS
        evidence: |
          target_layer value is playwright, which is one of the three
          canonical enum members (playwright, apitest, manual-only).
      - check: A-STRUCT-MECH-06
        status: PASS
        evidence: |
          coverage_type value is regression, one of the four canonical
          enum members (smoke, regression, edge, perf). Not a
          severity-axis value (p0..p3 belong to the
          critical_paths/priority axis only and would FAIL the unified
          test-kind enum).
      - check: A-STRUCT-03
        status: PASS
        evidence: |
          coverage_type label is declared at file frontmatter level
          (regression) and applies uniformly to both scenarios; the
          mode file explicitly permits file-frontmatter-level
          declaration. The Notes section justifies the choice
          (bug-focused regression guard for GROK-18721; smoke lives
          elsewhere per chain ui_coverage_plan; direct-link entry is a
          primary user path, not a boundary value, so not edge; not
          perf).
      - check: A-STRUCT-04
        status: PASS
        evidence: |
          Common preconditions are factored into `## Setup` (3 numbered
          steps): PowerPack installed and user access confirmed; a
          saved project with a known direct-link URL exists (with a
          fixture creation fallback via System:DemoFiles/demog.csv and
          File menu Save As Project); a fresh browser context is
          available. Neither scenario re-states this material; each
          starts from its own entry-path-specific step 1.
      - check: A-LAYER-ALIGN-01
        status: PASS
        evidence: |
          Frontmatter has no pyramid_layer key (verified across the
          full preamble block). PASS-by-vacuity applies per the mode
          file. The prose mention of "bug-focused per Rule 3" in the
          Notes section is documentation, not a frontmatter field, and
          the hard alignment rule applies only to the ui-smoke value
          paired with coverage_type smoke.
      - check: A-CONT-01
        status: PASS
        evidence: |
          Concrete names used throughout. Platform entities are cited
          by name (System:DemoFiles/demog.csv, File menu Save As
          Project, Browse Projects, Recent Projects widget, project
          name direct-link-loading-test). Source-file citations were
          verified against repo state: package.ts line 134 hosts the
          welcomeView registration with autostartImmediate true, line
          138 hosts the first `@grok.decorators.dashboard` entry
          (Spotlight), and line 315 hosts powerPackInit under
          `@grok.decorators.init`. The URL template uses `<server>`
          and `<user>` as environment-parametric markers paired with
          the bug report's concrete example
          dev.datagrok.ai/p/Opavlenko.Demog_95 — these are not generic
          stand-ins like `<column>` or `<TODO>`, and no
          atlas-hallucination patterns are present.
      - check: A-BUG-01
        status: PASS
        evidence: |
          Atlas powerpack.yaml known_issues (starts at line 1352) uses
          the modern schema where each entry carries test_coverage with
          exists false and paths empty list, rather than the legacy
          literal `test_coverage: needed`. Under strict literal
          predicate no entry qualifies, so PASS-by-vacuity. Under the
          equivalent-modern reading (exists false equates to needed,
          per the file-level comment at lines 1349 to 1351), the
          single entry scoped to this scenario is GROK-18721 at lines
          1373 to 1381, whose affects_sub_features list exactly
          matches this scenario's sub_features_covered list
          (powerpack.lifecycle.init, powerpack.welcome.view,
          powerpack.dashboards). GROK-18721 is addressed via both
          clause (a) — related_bugs lists GROK-18721 in frontmatter —
          and clause (b) — body references at the H1 title, the
          opening prose paragraph, the Related bug Notes entry, and
          the Chain context Notes entry. The other 8 atlas bugs
          target different sub-features (formula-lines, add-new-column,
          search.power, io.xlsx, widgets.recent-projects, …) and are
          chain-level concerns owned by F-BUG-COVERAGE-01, outside
          single-scenario Gate A scope.
      - check: A-MERIT-01
        status: PASS
        evidence: |
          No effort or complexity opt-outs anywhere. The Notes section
          closes with an explicit "Deferrals. None." sentence. Both
          scenarios are fully authored with concrete steps and
          expected results.
      - check: A-MERIT-02
        status: PASS
        evidence: |
          No TODO / later / next-phase markers in body or frontmatter
          (verified across the full file). The Notes section's
          Deferrals item explicitly states no deferral is required,
          citing that the two scenarios cover the bug's reproduction
          surface plus the regression-guard control.
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

Bug-focused regression scenario for `GROK-18721`: when a Datagrok app
or dashboard is opened by following a direct link (URL paste in a
fresh browser context, bookmark, or external deeplink), PowerPack's
loading window must render fully and properly — the same visual
quality as when navigating to the same app or dashboard from inside
the platform. Before the fix, the direct-link entry path produced a
cropped / incomplete loading-window layout: PowerPack's
`powerPackInit` (the package's `@init` handler at
`public/packages/PowerPack/src/package.ts#L315`) ran before the page's
layout dimensions were established, and the welcome / home dashboard
host rendered against undefined viewport sizing.

Atlas surface exercised:
- `powerpack.lifecycle.init` — `powerPackInit` runs at package init
  on a fresh load; the bug is entry-path-dependent because in-app
  navigation does not re-run init while direct-link entry does.
- `powerpack.welcome.view` — `welcomeView`
  (`public/packages/PowerPack/src/package.ts#L134`,
  `meta.autostartImmediate: true`) replaces the platform home view on
  direct-link entry; the loading window surfaces inside the welcome
  view host.
- `powerpack.dashboards` — the home-dashboard widgets
  (`public/packages/PowerPack/src/package.ts#L138`) host the loading
  affordance during direct-link page load; the cropping bug
  manifests against this surface.

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

- **target_layer rationale.** `playwright`. The bug is a visual
  rendering glitch in the PowerPack loading window during a specific
  entry path (direct-link). A headless JS-API exercise (apitest)
  cannot verify the loading window's rendering state — viewport
  dimensions, painting completeness, and visual cropping are
  browser-level observations. Pixel-precision drag and native file
  picker are not involved, so `manual-only` is not required.
- **coverage_type rationale.** `regression`. Bug-focused
  (`pyramid_layer: bug-focused` per Rule 3 — canonical GROK-18721
  reproduction surface), guards against re-regression on the
  direct-link entry path. Not `smoke` (this is not the section's
  golden path — the section's smoke is the top-level
  `add-new-column.md` per chain `ui_coverage_plan`). Not `edge` (the
  bug surfaces on a primary entry path — direct-link navigation is a
  common user action via URL paste, bookmark, deeplink — not on a
  boundary value or unusual input). Not `perf`.
- **Pyramid layer.** `bug-focused` per Rule 3 — discriminator test:
  GROK-18721 fails Scenario 1 before fix (PowerPack loading window
  renders cropped / incomplete during direct-link page load). After
  the fix, both scenarios pass and direct-link entry matches in-app
  navigation visual quality.
- **Atlas sub_features traceability.**
  - `powerpack.lifecycle.init` — `powerPackInit`
    (`public/packages/PowerPack/src/package.ts#L315`); the bug is
    entry-path-dependent on lifecycle init timing. Direct-link
    bypasses the warm-session shortcut and re-runs init against an
    unestablished viewport.
  - `powerpack.welcome.view` — `welcomeView`
    (`public/packages/PowerPack/src/package.ts#L134`,
    `meta.autostartImmediate: true`); the welcome view hosts the
    loading affordance on fresh entry.
  - `powerpack.dashboards` — the home-dashboard widgets
    (`public/packages/PowerPack/src/package.ts#L138`); cropping
    manifests against the dashboard host area during the welcome /
    loading transition.
- **Related bug.** `GROK-18721` (p2, status fixed).
  Reproduction: copy a direct link to an app or dashboard (e.g.
  `https://dev.datagrok.ai/p/Opavlenko.Demog_95`) → paste in a
  fresh browser → PowerPack loading window renders cropped /
  incomplete during page load. Expected: direct-link navigation
  must render PowerPack's loading window fully and properly with
  the same visual quality as in-app navigation.
- **Edge-case-for-atlas note (from bug-library).**
  entry-path-dependent rendering: lifecycle init timing on
  direct-link bypasses normal view-change handler; loading window
  renders before layout dimensions are established. The two
  scenarios above pair the broken direct-link path with the
  preserved in-app path to catch any future regression on either
  side.
- **Chain context.** This scenario is the section's bug-focused
  witness for GROK-18721 — Critic F's `bug_focused_candidates[]`
  entry for the bug (added in
  `cycle-2026-05-20-powerpack-coverage`) had empty `spans[]`
  pending this scenario's authoring. The chain's
  `order_from_files[]` is updated to include this scenario.
- **Deferrals.** None. The two scenarios cover the bug's
  reproduction surface (direct-link entry → cropped loading
  window) plus the regression-guard control (in-app navigation
  with the same target). No deferral required.
- **Coverage map.** Coverage map for PowerPack
  (`references/coverage-map/powerpack.yaml`) is not present at
  authoring time — gap-vs-coverage cross-check skipped per
  STEP B fallback. The Critic F coverage-gap dispatch
  (`gap: bug-uncovered :: GROK-18721`) drove this scenario's
  authoring directly.
