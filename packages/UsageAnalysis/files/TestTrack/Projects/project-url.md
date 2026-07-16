---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: [views.projects]
realized_as:
  - project-url-spec.ts
pyramid_layer: integration
ui_coverage_responsibility:
  - context-panel-links-url-copy
  - new-tab-open-url
ui_coverage_delegated_to: null
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/projects/project-url.md
migration_date: 2026-05-20
source_text_fixes: []
candidate_helpers:
  - helpers.playwright.projects.buildVariantsComposite
unresolved_ambiguities:
  - click-the-project-single-click-vs-double-click
  - context-panel-links-location
  - url-format-query-parameters
  - new-tab-vs-incognito-tab
  - order-vs-dependency-contradiction
  - source-text-correction-recurrence
scope_reductions:
  - id: SR-01
    check: E-SCENARIO-RUNTIME-ALIGNMENT
    rationale: |
      The existing `project-url-spec.ts` exercises only the `demog`
      representative source (1 of the 4 project variants) via direct URL
      navigation; the copied-with-link / copied-with-clone /
      personal-view-customizations variants are not deep-linked in the spec.
      The URL build/apply/shell-open contract is source-agnostic, so the
      single-variant walk preserves the invariant. NOTE: gate_verdicts.b is
      FAIL (spec unstable, [B-RUN-PASS, B-STAB-01]) — a separate open issue,
      not addressed by this scope reduction.
    verdict_status: SCOPE_REDUCTION
related_bugs: []
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: batch-6.6-pilot-2026-05-20-projects-project-url
    timestamp: 2026-05-20T00:00:00Z
    review_round: 1
    failure_keys: []
  d:
    verdict: PASS
    cycle_id: batch-6.6-pilot-2026-05-20-projects-project-url
    timestamp: 2026-05-20T00:00:00Z
    failure_keys: []
  e:
    verdict: SCOPE_REDUCTION
    cycle_id: batch-6.6-pilot-2026-05-20-projects-project-url
    timestamp: 2026-05-20T00:00:00Z
    review_round: 1
    failure_keys: []
  b:
    verdict: FAIL
    cycle_id: batch-6.6-pilot-2026-05-20-projects-project-url
    timestamp: 2026-05-20T12:35:00Z
    spec_runs:
      - spec: project-url-spec.ts
        result: failed
        attempts: 3
        duration_seconds: 70
        failure_keys: [B-RUN-PASS, B-STAB-01]
---

# Project URL — Deep-link reopen for saved project variants

For the `demog` project (a file-share source, used here as the
representative case since the URL deep-link flow is identical
regardless of the underlying data source), navigate to Browse >
Dashboards, locate each of its four variants — the original, the copy
made with Link mode, the copy made with Clone mode, and the copy
saved with Personal View Customizations — copy the deep-link URL from
Context Panel > Links for each, open that URL in a new browser tab,
and verify the corresponding project loads correctly.

This scenario depends on the three copy-mode variants produced by
`projects-copy-clone.md`, plus the original `demog` project (the
`Test_Case<N>` projects from `uploading.md` can substitute for the
"original" role when needed).

## Setup

1. **Fixture prerequisite:** the four project variants must exist on
   the server:
   - **original**: the saved `demog` project
     (Variant C representative source — file-share). The
     `Test_Case<N>_Sync` projects from `uploading.md` are an
     acceptable substitute when `demog` is unavailable in the chain
     run.
   - **copied-with-link**: the `<original>-link` Save-Copy-with-Link
     variant produced by `projects-copy-clone.md`.
   - **copied-with-clone**: the `<original>-clone` Save-Copy-with-Clone
     variant produced by `projects-copy-clone.md`.
   - **saved-with-personal-view-customizations**: the
     `<original>-personal-view-customizations` Save-personal-view-
     customizations variant produced by `projects-copy-clone.md`.
   In `end_to_end_fixtures` strategy, the Automator stage builds the
   composite `copy-clone-customizations-variants` fixture in a
   `beforeAll` block by either reusing fixtures from prior chain
   runs or replaying the relevant save-flow cases via `js-api-replay`
   on top of the `demog-project-with-viewers` baseline fixture.
2. The browser session is authenticated as the project owner (the
   user who produced the four variants upstream). Cross-user share
   verification is OUT of scope here.

## Scenarios

### URL deep-link reopen for each variant

The following workflow is parameterised over the 4 project variants.
Implemented per `end_to_end_fixtures` strategy at the Playwright layer
(per `output_plan.project-url.md.strategy = end_to_end_fixtures` in
`scenario-chains/projects.yaml` rev 3). Variant C source-agnostic
discipline: the URL-build / URL-apply / shell-open path is identical
across all four variants; the assertion verifies uniformity, not
source-class-specific behavior.

For each `<variant>` in the 4-variant list:

1. Go to **Browse** > **Dashboards**.
2. Click the project corresponding to `<variant>`:
   - original project
   - copied with the link
   - copied with clone
   - saved with personal view customizations
3. Go to **Context Panel** > **Links** and copy the URL shown there.
4. Open a new tab in the browser and paste/navigate to the copied URL.
5. **Verify on URL load (new tab):** the project corresponding to
   `<variant>` opens — its tables, viewers, and view layout render in
   the new tab; no errors in the browser console (F12).

## Notes

- **Why `demog` (file-share) is the representative source.** The
  Context-Panel-Links URL copy + new-tab open flow works identically
  regardless of the underlying data source, so testing it once is
  enough. `demog` was picked because it is the smallest, cheapest
  baseline project, and
  `projects-copy-clone.md` derives its three copy-mode variants from
  that same source. Other source classes (query, script, Spaces, DB
  table, derived) are covered by the separate
  `projects-lifecycle-*.md` scenarios, not by this one.
- **UI coverage owned here.** The clipboard-copy of the deep-link URL
  from Context Panel > Links, and the open-in-new-tab flow, are not
  covered by any other scenario in this section. Save / share / open
  / delete dialogs are owned elsewhere and are not exercised here.
- **Source-text correction.** The original scenario listed a
  "with layout" variant that had no actual producer in this section;
  it was corrected to "saved with personal view customizations" to
  match the mode that `projects-copy-clone.md` actually produces.
- **Cross-fixture dependency.** This is the only scenario in the
  section that needs a fixture composed from two different
  producers — the three copy-mode variants from
  `projects-copy-clone.md`, plus the original `demog` project.
  Coordinating that setup is more involved than
  the section's other scenarios.
- **No Share operation tested.** Step 3 copies the URL from Context
  Panel > Links — a clipboard copy of a deep-link, not the
  right-click Share dialog. Sharing is not exercised here.
- **Existing spec covers one variant only.** The existing
  `project-url-spec.ts` is a scope reduction that exercises only the
  `demog` representative source via direct URL navigation
  (`page.goto({BASE}/p/{nqName})`), not the full 4-variant loop
  described above.
