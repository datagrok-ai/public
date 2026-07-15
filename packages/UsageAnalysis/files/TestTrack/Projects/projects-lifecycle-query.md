---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [rename-dependent-entity-reopen, rename_external_dep, share_with_recipient_open, rename_project]
realizes: []
pyramid_layer: proactive-lifecycle
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/projects-lifecycle-query.md
migration_date: 2026-05-04
related_bugs:
  - github-3550
---

# Projects — Query-source lifecycle

Covers the lifecycle of a project sourced from a saved query on
`System:Datagrok` (the platform's built-in metadata Postgres
connection, always present): save, share with a second user, rename
the underlying query, and rename the project itself. This is the full
reproduction of github-3550 — renaming a query that a project depends
on must not silently break the project's reference to it (the table
should either auto-resolve or fail with an explicit error, never go
silently null). github-3550 is the query-side sibling of GROK-19212,
which covers the same failure mode for renamed tables.

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
- **github-3550 full reproduction.** Step 4 walks the bug's exact
  reproduction path: rename the external query, reopen the project,
  and check relation resolution. Because the test owns the query,
  the rename always succeeds — there's no permission-based skip.
- **UI coverage delegated.** All UI surfaces touched here
  (right-click rename, Save dialog, Sharing tab) are owned
  by `projects-ui-smoke.md`. This scenario uses the JS API path.
- **Deferred.** Recipient-side assertions are blocked on a
  not-yet-registered login-as-another-user test helper.
- **Self-cleaning.** Cleanup deletes the project and invokes
  `provisioned.cleanup()` to delete the saved query.
