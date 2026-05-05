---
feature: projects
sub_features_covered:
  - projects.upload
  - projects.api.save
  - projects.api.files.sync
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
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/projects-lifecycle-script.md
migration_date: 2026-05-04
migration_report: projects-lifecycle-script-migration-report.md
related_bugs:
  - GROK-19403
  - GROK-19728
---

# Projects — Script-source lifecycle

Chained lifecycle for projects sourced from a **Script** output
(`Samples:Cars`). Exercises proactive coverage cells
`source_class=script × dep_lifecycle_op=rename_external_dep` AND
`source_class=script × dep_lifecycle_op=share_with_recipient_open`
per chain rev 3 `proactive_lifecycle_specs[2]`. Targets
GROK-19403 (recipient cannot open shared project when underlying
script not also shared) AND GROK-19728 (view-and-use users edit
creation script under failure state — broken-script variant).

UI coverage delegated to `projects-ui-smoke.md`.

## Setup

1. Authenticate as test user.
2. Project name: `lifecycle-script-${Date.now()}`.
3. Recipient placeholder: `<RECIPIENT_USERNAME_TBD>`.
4. Helper 3 dependency: `helpers.playwright.session.logoutAndLoginAs`
   (NOT YET REGISTERED).
5. **Environment dependencies (env-provisioning required):**
   - Script: `Samples:Cars` provisioned for `qa-pw` test account.
     Skip-with-logged-warning if not provisioned.
   - For GROK-19403 sub-flow: the script must be **un-shared
     with the recipient** initially — otherwise the bug doesn't
     reproduce.
   - For GROK-19728 sub-flow: a **broken** variant of the script
     (intentionally introduce a syntax error or missing
     dependency) is needed. Surfaced as deferred item — broken-
     script provisioning may need a separate test fixture.
6. Cleanup: delete project; revoke permissions; restore Script
   name.

## Scenarios

### Main flow — script-source lifecycle (happy path + GROK-19403)

1. **Open table from Script.** Run `Samples:Cars` script (via
   `Browse > Platform > Functions > Scripts > Samples > Cars`
   double-click OR
   `await grok.functions.eval('Samples:Cars()')`). Verify the
   resulting table is loaded.
2. **Save project with Data Sync ON.** Save Project, name from
   Setup, Data Sync **ON**, OK. Cancel auto-share.
3. **Share project with second user (View-and-Use + Full) —
   script NOT also shared (GROK-19403 setup).**
   - Verify the script `Samples:Cars` is NOT shared with the
     recipient before this step (precondition for GROK-19403).
   - Use `grok.dapi.permissions.grant(project, recipient, ...)`
     for project-level grants.
   - **GROK-19403 invariant assertion:** when recipient opens
     the shared project (Step 3.3 below), they should see an
     **explicit permission failure** (e.g. "Could not access
     script Samples:Cars — you don't have permission"), NOT a
     silent null table or success-with-empty-data.
   - Original-user assertion: project still opens for original
     user (their access is fine; only recipient lacks script
     access).
   - **Recipient-side assertion (Helper 3 — deferred):** logout +
     login as recipient; open shared project; verify the
     **explicit error** is shown (NOT a silent null per
     GROK-19403). If recipient sees an empty table without
     error, that IS the GROK-19403 failure mode.
4. **Rename external dependency — Script rename.** Via right-
   click `Rename` on `Samples:Cars` OR via JS API:
   ```js
   const s = await grok.dapi.scripts.find('Samples:Cars');
   s.name = 'Cars_renamed';
   await grok.dapi.scripts.save(s);
   ```
   Skip-with-logged-warning if rename rejects on Samples-
   namespace permission.
   - Original-user assertion: project still opens; relation to
     the renamed script either auto-resolves OR fails with
     explicit error.
5. **Rename project itself.** Via JS API.
   - Original-user assertion: opens under new name.
   - Recipient-side (Helper 3 — deferred): opens under new name
     (still with the GROK-19403 explicit-error behavior on
     recipient side).
6. **GROK-19728 sub-flow — broken creation script + view-and-
   use access (DEFERRED — broken-script fixture not yet
   provisioned).**
   - Setup: change the script content to introduce a runtime
     error (e.g. `throw new Error('intentional break')`).
     **Provisioning blocker** — needs a separate broken-script
     fixture; deferred per Setup point 5.
   - Recipient with View-and-Use access opens the project.
   - **GROK-19728 invariant assertion:** recipient (view-and-
     use level) should NOT be able to edit the creation
     script — but the bug is that they CAN edit it under the
     failure state. Test asserts the edit is rejected (script
     editing UI is read-only for view-and-use users) OR is
     gracefully blocked.
   - Realized coverage: skipped pending broken-script fixture
     provisioning. share-side assertion (project opens with
     view-and-use access) runs unconditionally.
7. **Cleanup.** Delete project. Revoke permissions. Restore
   script name (rename back). Restore script content to working
   state if Step 6 ran.

### Expected results

- Save / reopen works for Script-sourced projects.
- **GROK-19403 invariant:** un-shared script causes explicit
  permission failure on recipient side, NOT silent null.
- **GROK-19728 invariant:** view-and-use users CANNOT edit the
  creation script even under failure state.
- Script rename + project rename stack correctly.

## Notes

- **Origin: chain rev 3 proactive_lifecycle_specs[2]** with
  `bugs_reinforcing: [GROK-19403]` and
  `dep_lifecycle_ops_covered: [rename_external_dep,
  share_with_recipient_open]`. GROK-19728 added as a
  reinforcement (per chain rev 3
  `bug_focused_candidates[GROK-19728]` semantic match). Authored
  in Phase A.
- **GROK-19403 full reproduction.** Step 3 walks the exact bug
  path: un-shared script as a project dependency + recipient
  share + recipient-open expectation. The cross-cutting bug
  spec `projects-grok-19403-spec.ts` (chain rev 3
  bug_focused_candidates) targets the share+recipient-open
  invariant; this scenario reinforces with full Script-source
  lifecycle.
- **GROK-19728 sub-flow deferred.** Step 6 requires a broken-
  script fixture which is NOT YET PROVISIONED. Surfaced as
  deferred item; share-side coverage runs unconditionally.
  Recipient-side broken-state assertion blocked by Helper 3
  + broken-script fixture (compound dependency).
- **Env dependency.** Script provisioned for qa-pw +
  recipient must have read access ONLY to project, NOT
  script (precondition for GROK-19403 reproduction).
- **`projects.api.relations.list` in sub_features.** Used in
  Step 4 to assert on relation resolution after Script
  rename.
- **UI coverage delegated.** All UI surfaces (right-click
  rename, Save dialog, Sharing tab) are owned by
  `projects-ui-smoke.md`. JS API path used here.
- **Helper 3 deferral.** Recipient-side assertions in Step
  3.3 + Step 5 + Step 6 require Helper 3.
- **Self-cleaning.** Step 7 deletes project + restores script.
- **Sequencing within Wave 2C.** Fourth lifecycle scenario per
  Plan line 219: needs env Script commit + broken-script
  fixture (the latter is GROK-19728 specific).
