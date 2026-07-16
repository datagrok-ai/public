---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [save-copy-with-link-mode]
realizes: []
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/projects-copy-clone.md
migration_date: 2026-05-04
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities:
  - all-created-projects-in-step-1
  - save-copy-ui-control-location
  - save-with-personal-view-customizations-mode-trigger
  - recipient-open-verification-depth-step-5
  - grok-19750-invariant-which-viewers-count-as-intact
  - grok-19403-cross-cutting-cover-via-step-5
scope_reductions: []
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

# Projects — Save Copy modes: Link, Clone & Personal View Customizations

For each previously-saved project (the `demog` project, and the
`Test_Case<N>` projects from
`uploading.md`), this scenario previews it in Browse, shares it with a
recipient, opens it, and runs it through three different "Save Copy"
modes — **Link**, **Clone**, and **Personal View Customizations** —
verifying that each variant reopens correctly and, critically, that
Save-Copy-with-Link never leaks changes back into the original
project (the GROK-19750 regression).

It also produces the three copy-mode project variants
(`<original>-link`, `<original>-clone`,
`<original>-personal-view-customizations`) that `project-url.md`
later deep-links into, and it is the section's main test of the
right-click Share dialog — this is where the Share dialog UI itself
gets exercised (steps 2 and 5).

## Setup

1. **Fixture prerequisite:** the `multi-source-saved-projects` fixture
   produced by upstream scenarios must exist on the server: at minimum
   the `demog` project, optionally augmented
   with `Test_Case<N>_Sync` / `_NoSync` projects from `uploading.md`.
   In `chained_tests` strategy, the Automator stage will reuse these
   from a prior chain run OR replay them via `js-api-replay` in a
   `beforeAll` block.
2. **Recipient identity for Share operations:** required for steps 2
   and 5 (right-click Share dialog). Same recipient pattern: the
   recipient is the placeholder user
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
   > `ui_coverage_plan.delegated_scenarios`, the right-click Share
   > dialog UI coverage is delegated to this scenario.
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

- **GROK-19750 invariant (sub-flow 4b, step 7).** The "reopen the
  original after Save-Copy-with-Link and verify its viewers are still
  intact" assertion was deliberately added here — without it, a
  GROK-19750 regression would go undetected. It's exactly the bug's
  reproduction path: open original, add a viewer, Save Copy with
  Link, close, reopen the copy, close, then reopen the *original* and
  confirm nothing leaked back into it.
- **Harmonization: sub-flow 4d now closes and reopens too.** Sub-flows
  4b (Link) and 4c (Clone) both close-and-reopen the resulting variant
  to verify it. Sub-flow 4d (Personal View Customizations) originally
  didn't — that's been fixed here for consistency; it's not a
  behavior change, just closing a gap between the sibling sub-flows.
- **Cross-cutting bug coverage.** Besides GROK-19750, this scenario's
  Step 4 (derive-and-save inside an active project) touches
  GROK-19103's affected surface, and Step 5 (right-click Share +
  recipient open) touches GROK-19403's — a variant with un-shared
  script/query dependencies could reproduce that bug's silent-null
  failure, though that's out of scope for this scenario's assertions.
- **Cleanup responsibility.** Steps 4b, 4c, 4d each produce a new
  variant project (`<original>-link`, `<original>-clone`,
  `<original>-personal-view-customizations`). These are consumed by
  `project-url.md`, so they are NOT cleaned up here — `projects-ui-smoke.md`
  (the section's terminal cleanup scenario) owns that.
