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
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/projects-lifecycle-query.md
migration_date: 2026-05-04
migration_report: projects-lifecycle-query-migration-report.md
related_bugs:
  - github-3550
---

# Projects — Query-source lifecycle

Chained lifecycle for projects sourced from a **provisioned saved
query** on `System:Datagrok` (the platform's built-in metadata
Postgres connection — always present). Exercises proactive coverage
cells `source_class=query × dep_lifecycle_op=rename_external_dep`
AND `source_class=query × dep_lifecycle_op=share_with_recipient_open`
per chain rev 3 `proactive_lifecycle_specs[1]`. Full reproduction
of github-3550 (external-Query rename invalidation, sister bug to
GROK-19212 for tables).

UI coverage delegated to `projects-ui-smoke.md`.

## Setup

1. Authenticate as test user.
2. Project name: `lifecycle-query-${Date.now()}`.
3. Recipient placeholder: `<RECIPIENT_USERNAME_TBD>`.
4. Helper 3 dependency: `helpers.playwright.session.logoutAndLoginAs`
   (NOT YET REGISTERED).
5. Cleanup: delete project; revoke permissions; delete the
   provisioned saved query (via `provisioned.cleanup()`).

## Scenarios

### Main flow — query-source lifecycle

0. **Provision saved query on `System:Datagrok`.** Use
   `helpers/openers.ts:provisionSystemDatagrokQuery({nameStem:
   'lifecycle_query', sql: SYSTEM_DATAGROK_QUERIES.GROUPS_SAMPLE})`.
   The helper creates a Postgres query against the `groups` table
   and returns `{queryId, queryNqName, resolvedName, cleanup}`.
   The query is namespaced under the test user's login — so the
   user has full edit/rename/delete rights. (This replaces the
   previous Samples-package env-dependency.)
1. **Open table from Query.** Use
   `helpers/openers.ts:openTableFromDbQuery(page,
   provisioned.queryNqName)`. Verify the resulting table is loaded
   and `df.tags['.script']` matches `<Var> = <queryNqName>()`.
2. **Save project with Data Sync ON.** Use
   `helpers/projects.ts:saveProjectWithProvenance(page,
   projectName)`.
3. **Share project with second user (View-and-Use).** Use
   `grok.dapi.permissions.grant(project, recipient, false)`.
   Defensive skip if no second user exists (Helper 3 deferred).
   - Original-user assertion: project reopens; query result loads
     via Data Sync.
   - Recipient-side assertion (Helper 3 — deferred): second user
     opens; query result loads. Recipient inherits access to
     `System:Datagrok` via `allUsers` group — no per-user DB
     credential dance required.
4. **Rename external dependency — Query rename (github-3550 full
   reproduction).**
   - Rename the provisioned query via JS API:
     ```js
     const q = await grok.dapi.queries.find(provisioned.queryId);
     q.name = `${provisioned.resolvedName}_renamed`;
     await grok.dapi.queries.save(q);
     ```
   - The test owns the query (it's namespaced under the test
     user's own login), so the rename always succeeds — no
     permission fallback.
   - **github-3550 invariant assertion:** the project's
     `relations.list` should now reference the renamed query
     (or the relation should be invalidated with an explicit
     error, not a silent null). Specifically: reopen the project;
     verify the table either loads (relation auto-resolved
     post-rename — happy path) OR fails with an explicit
     "Could not resolve query <queryNqName>" error (graceful
     failure — not silent corruption).
5. **Rename project itself.** Via JS API. Verify rename persists.
6. **Cleanup.** Delete project. Invoke `provisioned.cleanup()` —
   deletes the saved query (with its renamed name; id is stable
   across rename).

### Expected results

- Save / reopen works for Query-sourced projects.
- Share + recipient-open works (recipient inherits access to
  `System:Datagrok` via the `allUsers` group).
- **github-3550 invariant:** Query rename does NOT silently
  break the project's relations — either auto-resolution
  (preferred) or explicit failure (acceptable). Silent null
  references are the bug.
- Project rename + Query rename stack: both renames persist.

## Notes

- **Self-contained source provisioning.** This spec creates and
  deletes its own saved query — no Samples package, no
  env-provisioned DB. The `System:Datagrok` connection is
  built-in (created by `ServiceConnectionsMigration` on every
  deploy) and shared with `allUsers`.
- **github-3550 full reproduction.** Step 4 walks the bug's
  exact reproduction path: external Query rename + project
  reopen + relation-resolution check. Now that the test owns
  the query, the rename always executes — no permission-based
  skip path remains.
- **`projects.api.relations.list` in sub_features.** Step 4
  asserts on `(await grok.dapi.projects.find(<id>))`'s
  relations list to verify rename impact. This sub_feature is
  a primary observation surface for github-3550.
- **UI coverage delegated.** All UI surfaces touched here
  (right-click rename, Save dialog, Sharing tab) are owned
  by `projects-ui-smoke.md`. JS API path used here.
- **Helper 3 deferral.** Recipient-side assertions blocked.
- **Self-cleaning.** Cleanup deletes project + invokes
  `provisioned.cleanup()` to delete the saved query.
