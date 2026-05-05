---
feature: projects
sub_features_covered:
  - projects.upload
  - projects.api.save
  - projects.api.files.sync
  - projects.api.namespaces
  - projects.api.relations.list
  - projects.add-relation
  - projects.shell.share-via-context-menu
  - projects.api.get-by-id
target_layer: playwright
coverage_type: regression
pyramid_layer: proactive-lifecycle
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/projects-lifecycle-spaces.md
migration_date: 2026-05-04
migration_report: projects-lifecycle-spaces-migration-report.md
related_bugs:
  - GROK-18345
---

# Projects — Spaces-source lifecycle

Chained lifecycle for projects sourced from a **Space** (SPGIs Space
populated via JS API prelude). Exercises proactive coverage cells
`source_class=spaces × dep_lifecycle_op=rename_external_dep` AND
`source_class=spaces × dep_lifecycle_op=share_with_recipient_open`
per chain rev 3 `proactive_lifecycle_specs[3]`. Full reproduction of
GROK-18345 (recipient cannot open shared project that uses Spaces
dataset saved with data sync — triple cross-feature interaction:
ownership transfer + datasync under different user context + binary
dataframe resolution).

UI coverage delegated to `projects-ui-smoke.md`.

## Setup

1. Authenticate as test user.
2. Project name: `lifecycle-spaces-${Date.now()}`.
3. Recipient placeholder: `<RECIPIENT_USERNAME_TBD>`.
4. Helper 3 dependency: `helpers.playwright.session.logoutAndLoginAs`
   (NOT YET REGISTERED).
5. **Environment dependencies:**
   - SPGIs Space accessible to `qa-pw` test account. Skip-with-
     logged-warning if not provisioned.
   - Recipient must have access to the Space (or the share path
     fails on a different bug — un-shared deps).
6. **Spaces prelude (JS API).** Per
   `sa-2026-05-03-spaces-inline-prelude-pattern`, create the
   Space via JS API before Step 1:
   ```js
   const space = await grok.dapi.spaces.createRoot('test-projects-demo-spaces-${Date.now()}');
   const client = await grok.dapi.spaces.spaceClient(space.id);
   const file = (await grok.dapi.files.list('System:DemoFiles', false, 'demog.csv'))[0];
   await client.addEntity(file, /*link=*/true);
   ```
7. Cleanup: delete project; revoke permissions; delete the Space
   via `grok.dapi.spaces.delete(space)` (postlude).

## Scenarios

### Main flow — spaces-source lifecycle (GROK-18345 full repro)

1. **Open table from Space.** Navigate `Browse > Spaces >
   test-projects-demo-spaces-* > demog.csv`. Double-click. Verify
   the table loads.
2. **Save project with Data Sync ON.** Save Project, name from
   Setup, Data Sync **ON**, OK. Cancel auto-share.
3. **Share project with second user (View-and-Use + Full).** Use
   `grok.dapi.permissions.grant`. Verify Sharing tab.
   - Original-user assertion: project reopens; Spaces table
     loads via Data Sync.
   - **GROK-18345 invariant assertion (recipient-side, Helper 3
     — deferred):** logout + login as
     `<RECIPIENT_USERNAME_TBD>`; open shared project; verify the
     Spaces-sourced table loads correctly under the recipient's
     identity. The bug is that this fails — the recipient sees
     a missing-data error or silent null because of the triple
     interaction (ownership transfer + datasync execution as
     recipient + binary dataframe resolution from Space).
4. **Rename external dependency — Space rename.** Rename the
   Space via JS API:
   ```js
   space.name = 'test-projects-demo-spaces-renamed-${Date.now()}';
   await grok.dapi.spaces.save(space);
   ```
   - Original-user assertion: project still opens; Space relation
     resolves (auto-resolution preferred; explicit error
     acceptable; silent null is the bug).
   - Recipient-side (Helper 3 — deferred): same.
5. **Rename project itself.** Via JS API. Same dual assertion.
6. **Cleanup.** Delete project. Revoke permissions. Delete the
   Space (postlude per Setup point 7).

### Expected results

- Save / reopen works for Spaces-sourced projects.
- **GROK-18345 invariant:** recipient can open the shared
  Spaces-sourced project after data sync transfer — full
  recipient-side validation under data sync.
- Space rename + project rename stack correctly.

## Notes

- **Origin: chain rev 3 proactive_lifecycle_specs[3]** with
  `bugs_reinforcing: [GROK-18345]` and
  `dep_lifecycle_ops_covered: [rename_external_dep,
  share_with_recipient_open]`. Authored in Phase A.
- **GROK-18345 full reproduction (Trigger 2 motivating
  example).** This bug is the primary motivating example for
  chain-analyzer-prompt.md Edit 5 Trigger 2 (single-scenario
  non-adjacent steps). In `complex.md` it spans Step 1 + Step 2
  + Step 12 + Step 13 — non-adjacent. This scenario realizes
  the equivalent in a single linear flow: open from Space + save
  with sync + share + recipient-open. The cross-cutting bug spec
  `projects-grok-18345-spec.ts` (chain rev 3
  bug_focused_candidates with spans `uploading.md:Step 7 +
  complex.md:Step 12`) targets the cross-scenario invariant; this
  scenario reinforces with full Spaces lifecycle.
- **`projects.api.namespaces` in sub_features.** Spaces are
  namespaces in atlas terms. Step 1 navigates `Browse > Spaces`
  which uses this endpoint family.
- **Spaces prelude / postlude.** Per
  `sa-2026-05-03-spaces-inline-prelude-pattern`, the Space is
  created and torn down within the scenario's lifetime. Same
  pattern as `complex.md` Setup point 4.
- **Recipient must have Space access.** Precondition for the
  share path to succeed — recipient lacking Space access is a
  different bug (GROK-19403-style un-shared-deps).
- **UI coverage delegated.** All UI surfaces are owned by
  `projects-ui-smoke.md`. JS API path used here.
- **Helper 3 deferral.** Recipient-side assertions blocked.
- **Self-cleaning.** Step 6 deletes project + Space.
- **Sequencing within Wave 2C.** Fifth lifecycle scenario per
  Plan line 220: needs SPGIs Space commit. Sequencing flexible
  with specs 3, 4, 6.
