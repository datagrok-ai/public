---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p0
realizes: [share-spaces-datasync, rename_external_dep, share_with_recipient_open, rename_project]
pyramid_layer: proactive-lifecycle
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/projects-lifecycle-spaces.md
migration_date: 2026-05-04
related_bugs:
  - GROK-18345
---

# Projects — Spaces-source lifecycle

Covers the lifecycle of a project sourced from a Space (created via
JS API for the test): save with Data Sync on, share with a second
user, rename the Space, and rename the project. This is the full
reproduction of GROK-18345 — a recipient can't open a shared project
whose table comes from a Space and was saved with Data Sync, because
resolving that data under the recipient's own identity involves a
three-way interaction of ownership transfer, data-sync execution as a
different user, and binary dataframe resolution from the Space.

UI coverage delegated to `projects-ui-smoke.md`.

## Setup

1. Authenticate as test user.
2. Project name: `lifecycle-spaces-${Date.now()}`.
3. Recipient placeholder: `<RECIPIENT_USERNAME_TBD>`.
4. Helper 3 dependency: `helpers.playwright.session.logoutAndLoginAs`
   (token-based; needs DATAGROK_AUTH_TOKEN_2).
5. **Environment dependencies:**
   - SPGIs Space accessible to the test account. Skip-with-
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

- **GROK-18345 full reproduction.** In `complex.md`, this bug's
  reproduction path is spread across four non-adjacent steps
  (Step 1 + Step 2 + Step 12 + Step 13). This scenario realizes the
  equivalent in one linear flow: open from a Space, save with sync,
  share, then verify the recipient can open it. A separate
  cross-cutting spec also targets the cross-scenario version of this
  invariant (spanning a step in `uploading.md` and a step in
  `complex.md`); this scenario reinforces it with the full Spaces
  lifecycle.
- **Spaces prelude / postlude.** The Space used by this scenario is
  created via JS API at the start and deleted at the end — same
  pattern as `complex.md`'s Setup.
- **Recipient must have Space access.** This is a precondition for
  the share path to succeed — a recipient lacking Space access would
  hit a different bug (the GROK-19403-style un-shared-dependency
  failure).
- **UI coverage delegated.** All UI surfaces are owned by
  `projects-ui-smoke.md`. This scenario uses the JS API path.
- **Deferred.** Recipient-side assertions are blocked on a
  not-yet-registered login-as-another-user test helper.
- **Self-cleaning.** Step 6 deletes the project and the Space.
