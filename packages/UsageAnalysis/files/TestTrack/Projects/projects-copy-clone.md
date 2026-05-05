---
feature: projects
sub_features_covered: [projects.view.browse, projects.view.single, projects.shell.open, projects.shell.share-via-context-menu, projects.api.save, projects.add-relation, projects.add-link]
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/projects-copy-clone.md
migration_date: 2026-05-04
migration_report: projects-copy-clone-migration-report.md
ui_companion: projects-copy-clone-ui.md
pyramid_layer: bug-focused
ui_coverage_responsibility:
  - save-copy-with-link-dialog
  - save-copy-with-clone-dialog
  - save-personal-view-customizations-dialog
  - pcmdShareProject
  - share-dialog-recipients
  - context-panel-sharing-tab
ui_coverage_delegated_to: null
related_bugs: [GROK-19750, GROK-19103, GROK-19403]
---

# Projects Copy Clone

For each previously-created project (consumed from `upload-project.md`'s
`demog` and `uploading.md`'s `Test_Case<N>_Sync` / `_NoSync` projects),
preview in Browse, share with a recipient, open, and run a sequence of
edit-and-save operations across three save modes (Link, Clone, Personal
View Customizations). Verify the GROK-19750 invariant: after Save-Copy-
with-Link, the original project's viewers must remain intact. Re-share
the newly-created variant projects.

This scenario is the producer of three variant fixtures consumed by
`project-url.md`: `<original>-link`, `<original>-clone`,
`<original>-personal-view-customizations`. Implemented as
`chained_tests` per `scenario-chains/projects.yaml` rev 3 — same
fixture (`multi-source-saved-projects`) reused across all chained steps;
no tear-down between steps within a chain.

This scenario is also the chain's `bug-focused` witness for
**GROK-19750** (sub-flow 4b walks the bug's reproduction path; step 7
asserts the regression-coverage invariant and would fail before fix).
Cross-cutting cover for **GROK-19103** (Step 4 derive-and-save inside
active project) and **GROK-19403** (Steps 2 + 5 right-click Share +
recipient open) is delegated here per chain
`bug_focused_candidates`. Per chain rev 3 `pyramid_layer: bug-focused`
and per `ui_coverage_plan.delegated_scenarios`, this scenario OWNS the
right-click Share dialog UI surface delegated from `share-project.md`
(steps 2 and 5 drive `pcmdShareProject`).

## Setup

1. **Fixture prerequisite:** the `multi-source-saved-projects` fixture
   produced by upstream scenarios must exist on the server: at minimum
   the `demog` project from `upload-project.md`, optionally augmented
   with `Test_Case<N>_Sync` / `_NoSync` projects from `uploading.md`.
   In `chained_tests` strategy, the Automator stage will reuse these
   from a prior chain run OR replay them via `js-api-replay` in a
   `beforeAll` block.
2. **Recipient identity for Share operations:** required for steps 2
   and 5 (right-click Share dialog). Same recipient pattern as
   `share-project.md`: the recipient is the placeholder user
   `<RECIPIENT_USERNAME_TBD>` (single test account, pending
   provisioning) AND/OR an email recipient (auto-creates a user
   account on the server). The Automator stage owns `afterAll`
   cleanup of any auto-created accounts.
3. The browser session must be authenticated as the project owner
   (otherwise the Save Copy and Save modes will be unavailable in
   the project menu).

## Scenarios

### Preview, share, open, and three-mode save sequence with GROK-19750 invariant

The following workflow runs as a chained sequence over the
`multi-source-saved-projects` fixture. Steps share state across the
chain (no `closeAll` between steps unless explicitly noted).

1. Go to **Browse** > **Dashboards**.

   > **UI-only verification moved to `projects-copy-clone-ui.md`:** the
   > preview thumbnails render-quality check (project name, picture/
   > dashboard preview, basic metadata visible across project tiles).
2. Share each previously-created project: right-click the project in
   Browse > Dashboards, select **Share**, fill recipients (registered
   user and/or email), submit. **Verify on Context Panel — Sharing
   tab:** each share completes; recipients appear on the project's
   Sharing tab.

   > Drives `pcmdShareProject` + `share-dialog-recipients` +
   > `context-panel-sharing-tab` UI flows. Per chain rev 3
   > `ui_coverage_plan.delegated_scenarios`, `share-project.md`
   > delegates the right-click Share dialog UI coverage to this
   > scenario.
3. Open each previously-created project: double-click in Browse >
   Dashboards. **Verify on open:** project loads, its tables and
   viewers render in the workspace, no errors in the browser console.
4. **Edit projects, save, and reopen** — the following four sub-flows
   run sequentially over the original project (e.g. `demog`). The
   original project state must be preserved between sub-flows for the
   GROK-19750 invariant to be testable.

   > **Spec-time restructure (2026-05-05):** the canonical scenario
   > below describes 4 sequential reopen-edit-save cycles. The spec
   > implementation collapses these into a SINGLE session
   > (1 open + cumulative addViewer + 4 sequential
   > `saveProjectWithProvenance` calls — see
   > `projects-copy-clone-spec.ts`), then verifies the
   > GROK-19750 invariant via 2 reopens (link copy + original).
   > Reason: each `closeAll` + `dapi.projects.find().open()` cycle
   > on dev under Playwright headed mode races with shell.tv
   > re-materialization; reducing reopens from 4-6 to 2 reduces
   > flake. The semantic invariant (Save Copy variants exist
   > server-side; original survives Save Copy without state leak)
   > is preserved.
   - **4a. Save the original (baseline):**
     1. Open the original project (e.g. `demog`).
     2. Add any viewer (e.g. another bar chart) to the table view.
     3. Save the project (overwrite the original — File > Save Project,
        keep the same name, OK).
     4. **Close All.**
   - **4b. Save Copy with Link:**
     1. Open the original project (e.g. `demog`).
     2. Add any viewer (e.g. another scatter plot) to the table view.
     3. Trigger **Save Copy** → choose **Link** mode → name the copy
        (e.g. `<original>-link`) → OK.
     4. **Close All.**
     5. Reopen the new copy (`<original>-link`) from Browse > Dashboards.
        **Verify on copy reopen:** the copy loads with the new viewer
        intact; tables render via the link to the original's data; no
        console errors.
     6. **Close All.**
     7. **GROK-19750 INVARIANT (mandatory regression-coverage assertion):**
        Reopen the **original** project (e.g. `demog`) from Browse >
        Dashboards. **Verify on original reopen:** the original
        project's viewers (those existing before sub-flow 4b plus the
        sub-flow 4b addition) are still present and intact —
        Save-Copy-with-Link MUST NOT have leaked state changes back
        into the source project. If the original opens without one or
        more of its viewers, this is a regression of GROK-19750.
     8. **Close All.**
   - **4c. Save Copy with Clone:**
     1. Open the original project (e.g. `demog`).
     2. Add any viewer (e.g. a line chart) to the table view.
     3. Trigger **Save Copy** → choose **Clone** mode → name the copy
        (e.g. `<original>-clone`) → OK.
     4. **Close All.**
     5. Reopen the new copy (`<original>-clone`) from Browse >
        Dashboards. **Verify on copy reopen:** the copy loads with
        the new viewer intact; tables are independent copies (not
        linked to the original); no console errors.
     6. **Close All.**
   - **4d. Save with Personal View Customizations:**
     1. Open the original project (e.g. `demog`).
     2. Add any viewer (e.g. a histogram) or apply a personal view
        customization (filter, sort, column show/hide, etc.) to the
        table view.
     3. Trigger **Save** with **personal view customizations** mode →
        name the variant (e.g. `<original>-personal-view-customizations`)
        → OK.
     4. **Close All.**
     5. Reopen the new variant (`<original>-personal-view-customizations`)
        from Browse > Dashboards. **Verify on variant reopen:** the
        variant loads (tables present, no console errors); the
        underlying data tables are shared with the original (similar
        to Link mode — assertable via JS API: post-modify source
        table, reopen variant, check propagation as expected).

        > **UI-only verification moved to `projects-copy-clone-ui.md`:**
        > confirming the view-state customization (filter / sort /
        > column show-hide / layout positioning) is preserved on the
        > variant — visual judgment.

        *(Harmonization addition: close-and-reopen verification was
        missing from this sub-step in the source `.md`; added at
        migration time for symmetry with sub-flows 4b and 4c. See
        migration report.)*
     6. **Close All.**
5. Re-share each newly-created variant project (`<original>-link`,
   `<original>-clone`, `<original>-personal-view-customizations`):
   right-click the variant in Browse > Dashboards, select **Share**,
   fill recipients, submit. **Verify on share + recipient open:**
   each variant is shareable; recipient access succeeds (the
   recipient session can open the shared variant; tables and viewers
   render; no console errors).

   > Drives `pcmdShareProject` + `share-dialog-recipients` for the
   > variant projects produced in sub-flows 4b/4c/4d. The recipient-
   > open path is also the candidate cross-cutting cover for
   > GROK-19403 (share + recipient open with potentially un-shared
   > dependencies — see chain rev 3 `bug_focused_candidates`).

## Notes

- **Original `order: 5`** — runs after `upload-project.md` /
  `uploading.md` (`order: 1`), `share-project.md` (`order: 2`),
  `opening.md` (`order: 3`), `project-url.md` (`order: 4`). Captured
  in `scenario-chains/projects.yaml` rev 3.
- **GROK-19750 invariant addition (sub-flow 4b step 7).** The original
  source `.md` does NOT contain an explicit "reopen the original after
  Save-Copy-with-Link and verify original's viewers are intact" step.
  This assertion was added at migration time per the per-cycle prompt's
  pre-known context: "MUST add an explicit assertion ... mandatory
  migration-time addition, not a candidate". Without this assertion,
  GROK-19750 regression coverage is missed by the chain. The bug's
  reproduction (per `bug-library/projects.yaml` rev 2 entry GROK-19750)
  is exactly the sub-flow 4b sequence; verifying the original's
  viewer-intactness on reopen IS the regression test.
- **Harmonization gap addition (sub-flow 4d steps 4-6).** The original
  source's 4th sub-bullet ("save personal view customizations")
  lacked the "close all, reopen it. Close all" verification that the
  preceding sub-bullets (link, clone) include. Sibling sub-flows 4b
  and 4c BOTH have the close-and-reopen verification; sub-flow 4d
  did not. This was harmonized at migration time — not a behavior
  change, just symmetry with the sibling pattern. See migration
  report for the source-text-correction history (`mig-2026-04-29-source-
  text-correction` recorded the "with layout" → "personal view
  customizations" terminology change but did NOT add the missing
  close-and-reopen verification — that gap is closed here).
- **Pyramid layer (chain rev 3): `bug-focused`.** Discriminator: the
  scenario was authored to walk GROK-19750's reproduction path
  (sub-flow 4b: open original → add viewer → Save Copy with Link →
  close all → reopen original → assert viewers intact). Sub-flow 4b
  step 7 IS the GROK-19750 regression test; without the bug fix, this
  scenario fails. Bug-focused dominates over Rule 4 multi-source per
  the heuristic priority 3 → 4.
- **UI coverage (chain rev 3): owns 6 flows.** This scenario owns the
  full Save-Copy mode dialog surface (`save-copy-with-link-dialog`,
  `save-copy-with-clone-dialog`,
  `save-personal-view-customizations-dialog`) plus the right-click
  Share dialog (`pcmdShareProject`, `share-dialog-recipients`,
  `context-panel-sharing-tab`). Per
  `ui_coverage_plan.delegated_scenarios`, `share-project.md`
  delegates UI coverage of the right-click Share dialog
  (`pcmdShareProject`) to this scenario — share-project.md uses
  JS API (`grok.dapi.permissions.grant`) for the registered-user
  share path, and the right-click Share dialog UI surface lives
  here in steps 2 and 5.
- **Cross-cutting bug awareness (chain rev 3 `bug_focused_candidates`):**
  - **GROK-19750**: cross-cutting candidate spec
    `projects-grok-19750-spec.ts` proposed (spans
    `upload-project.md:Step 1` + `projects-copy-clone.md:Step 4`).
    F-BUG-COVERAGE-01 at section-complete is the authoritative gate
    for cross-cutting bug coverage; this scenario is the primary
    bug-focused witness inside the chain.
  - **GROK-19103**: cross-cutting candidate spec
    `projects-grok-19103-spec.ts` proposed (spans include
    `projects-copy-clone.md:Step 4`). The Save Copy + viewer-add flow
    here intersects the bug's `affects: [projects.api.save,
    projects.add-relation]` surface.
  - **GROK-19403**: cross-cutting candidate spec
    `projects-grok-19403-spec.ts` proposed (spans include
    `projects-copy-clone.md:Step 5`). Step 5 right-click Share +
    recipient open is the share-with-deps reproduction path; if the
    variant has un-shared script/query deps (out of scope here) the
    recipient sees a silent null per the bug.
  - All three bugs are listed in `related_bugs:` frontmatter as
    chain-cross-cutting watch-list. Per Edit 5 Design point A, the
    citation is RECOMMENDED (early visibility) — the authoritative
    coverage gate is F-BUG-COVERAGE-01 at section-complete.
- **Invariant 3 (atlas-aware sub_features_covered for share) APPLIES.**
  Steps 2 and 5 both exercise the right-click Share dialog
  (`pcmdShareProject` per atlas rev 11 sub_feature
  `projects.shell.share-via-context-menu`). The sub_feature is
  correctly INCLUDED in `sub_features_covered`.
- **Cleanup responsibility.** Steps 4b, 4c, 4d each produce a NEW
  variant project (`<original>-link`, `<original>-clone`,
  `<original>-personal-view-customizations`). These ARE the
  consumed-by-`project-url.md` fixture (`copy-clone-customizations-
  variants` per scenario-chains rev 3). They are NOT cleaned up by
  this scenario — `deleting.md` (terminal, must_run_last) owns
  terminal cleanup. Step 5's email-invite-creates-account side
  effect IS owned by Automator's `afterAll` per per-cycle convention
  (same as `share-project.md`).
- **Sibling spec convention.** `projects-copy-clone-spec.ts` exists
  from a prior cycle; per-cycle Invariant 2 (existing -spec.ts /
  -api.ts READ-ONLY) applies — Migrator does NOT modify it. Adjacent
  specs (`opening-spec.ts`, `uploading-spec.ts`, `share-project-spec.ts`,
  `deleting-spec.ts`) all use the `loginToDatagrok`, `softStep`,
  `evalJs`, `closeAll` convention and `Date.now()` suffix naming
  pattern; consulted READ-ONLY for house style.
- **"Add any viewer" cross-feature touch.** Sub-flows 4a/4b/4c/4d each
  begin with "add any viewer". Adding a viewer to a TableView is a
  cross-feature operation (viewer system, not projects.* directly).
  The closest projects.* sub_features touched are
  `projects.shell.add-dataframe` (when the viewer addition triggers
  a project state mutation) and `projects.api.save` (when the project
  is saved with the new viewer). Listed in `sub_features_covered`
  for the saving aspect; viewer addition itself is rendered in the
  D4 viewer system, out of projects.* scope.
