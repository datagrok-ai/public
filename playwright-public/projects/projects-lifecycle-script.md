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
related_bugs:
  - GROK-19403
  - GROK-19728
---

# Projects — Script-source lifecycle

Chained lifecycle for projects sourced from a **provisioned dataframe-
output Script** (created in-test via
`helpers/openers.ts:provisionDataframeScript`). Exercises proactive
coverage cells `source_class=script ×
dep_lifecycle_op=rename_external_dep` AND `source_class=script ×
dep_lifecycle_op=share_with_recipient_open` per chain rev 3
`proactive_lifecycle_specs[2]`. Targets GROK-19403 (recipient cannot
open shared project when underlying script not also shared) AND
GROK-19728 (view-and-use users edit creation script under failure
state — broken-script variant).

UI coverage delegated to `projects-ui-smoke.md`.

## Setup

1. Authenticate as test user.
2. Project name: `lifecycle-script-${Date.now()}`.
3. Recipient placeholder: `<RECIPIENT_USERNAME_TBD>`.
4. Helper 3 dependency: `helpers.playwright.session.logoutAndLoginAs`
   (NOT YET REGISTERED).
5. Cleanup: delete project; revoke permissions; delete the
   provisioned script (via `helpers/openers.ts:deleteProvisionedScript`).

## Scenarios

### Main flow — script-source lifecycle (happy path + GROK-19403)

0. **Provision JS script with `output: dataframe`.** Use
   `helpers/openers.ts:provisionDataframeScript({name:
   'lifecycleScript${stamp}', body: "df = await
   grok.data.getDemoTable('demog.csv');"})`. The helper creates
   the script via `DG.Script.create` + `grok.dapi.scripts.save`,
   waits for `DG.Func` registration, and returns `{scriptId,
   resolvedName, resolvedNqName}`. The script is namespaced
   under the test user's login — full edit/rename/delete rights.
   (This replaces the previous Samples:Cars env-dependency,
   which was scalar-output and silently skipped.)
1. **Open table from Script.** Use
   `helpers/openers.ts:openTableFromScript(page,
   provisioned.resolvedNqName, {idx: 0})`. Verify the resulting
   table is loaded and `df.tags['.script']` matches
   `<Var> = <resolvedName>(idx=0)`.
2. **Save project with Data Sync ON.** Use
   `helpers/projects.ts:saveProjectWithProvenance(page,
   projectName)`.
3. **Share project with second user (View-and-Use + Full) —
   script NOT also shared (GROK-19403 setup).**
   - The provisioned script is owned by the test user and not
     shared with anyone else by default — precondition for
     GROK-19403 reproduction holds automatically.
   - Use `grok.dapi.permissions.grant(project, recipient, false)`
     for project-level grants. Defensive skip if no second user
     exists (Helper 3 deferred).
   - **GROK-19403 invariant assertion:** when recipient opens
     the shared project, they should see an **explicit
     permission failure** ("Could not access script
     <resolvedName> — you don't have permission"), NOT a silent
     null table.
   - **Recipient-side assertion (Helper 3 — deferred):** logout +
     login as recipient; open shared project; verify the
     **explicit error** is shown.
4. **Rename external dependency — Script rename.** Via JS API:
   ```js
   const s = await grok.dapi.scripts.find(provisioned.scriptId);
   s.name = `${provisioned.resolvedName}_renamed`;
   await grok.dapi.scripts.save(s);
   ```
   The test owns the script — rename always succeeds.
   - Original-user assertion: project still opens; relation to
     the renamed script either auto-resolves OR fails with
     explicit error (mirrors github-3550 invariant for queries).
5. **Rename project itself.** Via JS API.
   - Original-user assertion: opens under new name.
   - Recipient-side (Helper 3 — deferred): opens under new name
     (still with the GROK-19403 explicit-error behavior).
6. **GROK-19728 sub-flow — broken creation script + view-and-
   use access (DEFERRED).**
   - Setup: re-provision (or update) the script content to
     introduce a runtime error (e.g. `throw new Error('intentional
     break')`). With our provisioning helper this is now
     trivially achievable — but the recipient-side assertion
     still needs Helper 3.
   - Recipient with View-and-Use access opens the project.
   - **GROK-19728 invariant assertion:** recipient (view-and-
     use level) should NOT be able to edit the creation
     script — but the bug is that they CAN edit it under the
     failure state. Test asserts the edit is rejected (script
     editing UI is read-only for view-and-use users) OR is
     gracefully blocked.
   - Realized coverage: deferred until Helper 3 lands.
7. **Cleanup.** Delete project. Revoke permissions. Delete the
   provisioned script via `deleteProvisionedScript(page,
   provisioned.scriptId)` (id is stable across rename).

### Expected results

- Save / reopen works for Script-sourced projects.
- **GROK-19403 invariant:** un-shared script causes explicit
  permission failure on recipient side, NOT silent null.
- **GROK-19728 invariant:** view-and-use users CANNOT edit the
  creation script even under failure state.
- Script rename + project rename stack correctly.

## Notes

- **Self-contained source provisioning.** This spec creates and
  deletes its own dataframe-output script — no Samples package,
  no env-provisioned fixture. The previous `Samples:Cars`
  reference was scalar-output (rejected no-input call) and the
  spec silently skipped; the in-test provisioned script wraps
  `grok.data.getDemoTable('demog.csv')` and returns a real
  dataframe.
- **Origin: chain rev 3 proactive_lifecycle_specs[2]** with
  `bugs_reinforcing: [GROK-19403]` and
  `dep_lifecycle_ops_covered: [rename_external_dep,
  share_with_recipient_open]`. GROK-19728 added as a
  reinforcement.
- **GROK-19403 full reproduction.** Step 3 walks the exact bug
  path: un-shared script as a project dependency + recipient
  share + recipient-open expectation. The provisioned script is
  owned by the test user and not auto-shared, so the GROK-19403
  precondition holds without extra setup.
- **GROK-19728 sub-flow deferred.** Step 6 needs Helper 3
  (logoutAndLoginAs) to assert recipient-side behavior. Broken-
  script provisioning is now trivial via the helper.
- **`projects.api.relations.list` in sub_features.** Used in
  Step 4 to assert on relation resolution after Script
  rename.
- **UI coverage delegated.** All UI surfaces (right-click
  rename, Save dialog, Sharing tab) are owned by
  `projects-ui-smoke.md`. JS API path used here.
- **Helper 3 deferral.** Recipient-side assertions in Step
  3.3 + Step 5 + Step 6 require Helper 3.
- **Self-cleaning.** Step 7 deletes project + provisioned
  script.
