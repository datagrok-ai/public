---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [share_with_recipient_open, rename_project]
realizes: [views.projects]
realized_as:
  - projects-lifecycle-files-spec.ts
pyramid_layer: proactive-lifecycle
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/projects-lifecycle-files.md
migration_date: 2026-05-04
related_bugs: []
---

# Projects — File-source lifecycle

Covers the lifecycle of a project sourced from a file share
(`System:DemoFiles/demog.csv`): save, share with a second user at two
access levels, and project rename — verifying the project stays
accessible to both the owner and the recipient throughout. Files have
no separate entity to rename (the file just lives at a fixed path in
the share), so this scenario's rename coverage is limited to renaming
the project itself.

UI coverage for share / delete / right-click flows lives in
`projects-ui-smoke.md`; this scenario drives sharing and permission
grants through the JS API instead, since repeatedly UI-driving those
flows here would be fragile and isn't this scenario's focus.

## Setup

1. Authenticate to Datagrok as the test user. Session is Playwright;
   share + permission grants use `grok.dapi.permissions.grant` (JS
   API substitute for right-click Share — UI is owned by
   `projects-ui-smoke.md`).
2. Project name: `lifecycle-files-${Date.now()}`.
3. Recipient placeholder: `<RECIPIENT_USERNAME_TBD>` (single user;
   resolved at Automator stage).
4. Helper dependency: this scenario's Step 3 recipient-open
   verification requires
   `helpers.playwright.session.logoutAndLoginAs` (Helper 3 from
   helpers-batch-1). NOT YET REGISTERED in `helpers-registry.yaml`.
   Step 3 recipient-side assertion is deferred pending registration;
   share-side assertion (permission grant + Sharing-tab listing) is
   NOT deferred and runs unconditionally.
5. Cleanup: delete the project at the end (Step 6); revoke the
   permission grant; log out the second-user session if Helper 3 is
   registered.

## Scenarios

### Main flow — files lifecycle

1. **Open table from File share.** Open `System:DemoFiles/demog.csv`
   via `grok.data.loadTable` (or via Browse > Files double-click —
   either path is acceptable; the scenario does not own the
   `pcmdOpen` UI surface, that's `projects-ui-smoke.md`'s job).
   Verify `grok.shell.tables.length > 0` and the `demog` table is
   loaded.
2. **Save project with Data Sync ON.** Trigger Save Project via the
   ribbon SAVE button (NOT Ctrl+S; per
   `feedback_no_ctrlS_for_layouts`). In the Save dialog: project
   name from Setup, **Data Sync** toggle **ON**, click **OK**.
   Cancel the auto-share dialog. Verify `POST /projects` succeeds
   and the project appears in `grok.dapi.projects.find(<id>)`.
3. **Share project with second user, recipient opens.**
   - Share via JS API:
     `grok.dapi.permissions.grant(project, recipient, /*edit=*/false)`
     for View-and-Use, then a second grant with `/*edit=*/true` for
     Full. Verify the Sharing tab on the Context Panel lists the
     recipient at the granted access levels (UI verification of the
     LIST is OK; UI driving of the GRANT is owned by
     `projects-ui-smoke.md`).
   - **Original-user assertion:** verify the project still opens
     for the original user (close-and-reopen via
     `grok.dapi.projects.find(<id>)` then load tables).
   - **Recipient-side assertion (deferred — Helper 3):** logout
     and login as `<RECIPIENT_USERNAME_TBD>` via
     `helpers.playwright.session.logoutAndLoginAs`; navigate to
     Browse > Dashboards > Shared with me; double-click the
     project tile; verify it opens; verify the `demog` table is
     accessible. **If Helper 3 is not registered, this sub-step
     skips with a logged warning — share-side assertion above is
     the realized coverage.**
4. **Rename external dependency — N/A for files.** File-source
   projects have no separate entity to rename (the file lives in
   the share at a fixed path). This step is a no-op for this
   source class. Cite chain rev 3
   `proactive_lifecycle_specs[0].dep_lifecycle_ops_covered:
   [share_with_recipient_open]` — `rename_external_dep` is NOT in
   this entry's coverage list.
5. **Rename project itself.** Trigger via JS API
   `project.name = '<original-name>_renamed'; await
   grok.dapi.projects.save(project)`. Verify the rename persists:
   `(await grok.dapi.projects.find(project.id)).name` returns
   the new name.
   - **Original-user assertion:** project still opens under the
     new name (`grok.dapi.projects.find(<id>)` then load tables).
   - **Recipient-side assertion (Helper 3 — deferred):** project
     still opens for the second user under the new name.
6. **Cleanup.** Delete the project via
   `grok.dapi.projects.delete(project)`. Revoke the permission
   grant via
   `grok.dapi.permissions.revoke(project, recipient)`. If Helper 3
   ran, log out the second-user session and restore the original-
   user session.

### Expected results

- Share + recipient-open works for files-source projects (the
  baseline cell).
- Project rename does NOT break the share — recipient still opens
  the project under its new name.
- No silent persistence drops on rename.

## Notes

- **No related bug.** Files-source has no associated GROK bug — this
  is pure proactive coverage.
- **UI coverage delegated.** All UI surfaces touched in this scenario
  (Save dialog, auto-share dialog dismiss, Sharing tab listing,
  Browse > Dashboards opens) are owned by `projects-ui-smoke.md`.
  The flow here uses the JS API where possible; UI is incidental at
  most.
- **Deferred.** Recipient-side assertions (logout, log in as the
  second user, verify the project opens) require a
  login-as-another-user test helper that isn't registered yet — same
  deferral as `complex-share-second-user-spec.ts` Step 13. The
  scenario degrades gracefully: the share-side assertions run
  unconditionally, and the recipient-side check only runs once the
  helper lands.
- **Self-cleaning.** Step 6 deletes the project.
