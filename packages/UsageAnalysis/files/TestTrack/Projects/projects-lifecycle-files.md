---
feature: projects
sub_features_covered:
  - projects.upload
  - projects.api.save
  - projects.api.files.sync
  - projects.shell.share-via-context-menu
  - projects.api.get-by-id
target_layer: playwright
coverage_type: regression
pyramid_layer: proactive-lifecycle
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/projects-lifecycle-files.md
migration_date: 2026-05-04
migration_report: projects-lifecycle-files-migration-report.md
related_bugs: []
---

# Projects — File-source lifecycle

Chained lifecycle for a project sourced from a **File share**
(`System:DemoFiles/demog.csv`). Exercises the proactive coverage cell
`source_class=files × dep_lifecycle_op=share_with_recipient_open` per
chain rev 3 `proactive_lifecycle_specs[0]`. Adds source-agnostic
`rename_project` op interleaved per Variant C. No external dependency
to rename (Files have no separate entity to rename — the file lives in
the share).

UI coverage for share / delete / right-click flows lives in
`projects-ui-smoke.md`. This scenario uses `grok.dapi.*` for share +
permission grants by design — UI driving for those flows is fragile
and unrelated to this scenario's regression invariant.

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

- **Origin: chain rev 3 proactive_lifecycle_specs[0].** This
  scenario realizes the `source_class=files` × `dep_lifecycle_op=
  share_with_recipient_open` proactive cell per Edit 7. Authored
  in Phase A of `projects migrate --force` cycle per Olena's
  Option γ.
- **Pyramid layer: integration (per chain rev 3 schema).** Strict
  pyramid_layer enum only allows ui-smoke / source-matrix /
  bug-focused / integration / manual. Per the Plan, this is
  conceptually a "proactive-lifecycle" scenario; the closest
  enum value is `integration`. Surfacing for retro: chain
  analyzer schema may want a dedicated `proactive-lifecycle`
  enum value.
- **bugs_reinforcing: [] (per chain rev 3
  proactive_lifecycle_specs[0]).** Files-source has no
  associated GROK bug — this is pure proactive coverage.
- **No external dependency rename.** Files have no separate
  entity to rename. `dep_lifecycle_ops_covered` is
  `[share_with_recipient_open]` only; `rename_external_dep` is
  NOT in scope for this entry.
- **UI coverage delegated.** All UI surfaces touched in this
  scenario (Save dialog, auto-share dialog dismiss, Sharing tab
  listing, Browse > Dashboards opens) are owned by
  `projects-ui-smoke.md`. The flow here is JS API where
  possible; UI is incidental at most.
- **Helper 3 deferral.** Recipient-side assertions (logout +
  login as second user + verify open) require
  `helpers.playwright.session.logoutAndLoginAs`. NOT YET
  REGISTERED. Same deferral as `complex-share-second-user-spec.ts`
  Step 13. Scenario degrades gracefully — share-side assertions
  run unconditionally; recipient-side runs only when Helper 3
  lands.
- **Self-cleaning.** Step 6 deletes the project. NOT in
  `deleting.md.depends_on`.
- **Sequencing within Wave 2C.** First lifecycle scenario
  per PROJECTS-SPLIT-COMPLETION-PLAN line 216-217: only blocker
  is Helper 3. No env provisioning required.
