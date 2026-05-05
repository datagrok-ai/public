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

Chained lifecycle for projects sourced from a **Query**
(`Samples:PostgresCustomers`). Exercises proactive coverage cells
`source_class=query × dep_lifecycle_op=rename_external_dep` AND
`source_class=query × dep_lifecycle_op=share_with_recipient_open`
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
5. **Environment dependencies (env-provisioning required):**
   - DB connection: `Postgres:NorthwindTest` accessible to test
     user. Skip-with-logged-warning if not present in environment
     (per `sa-2026-05-03-postgres-queries-public-data-substitution`
     pattern).
   - Query: `Samples:PostgresCustomers` provisioned for `qa-pw`
     test account. Skip-with-logged-warning if not provisioned.
   - Recipient must have access to the same DB connection (or
     the share path will fail with a permission error per
     GROK-19403; this scenario tests the Query rename path, NOT
     the un-shared-deps path — recipient access to the
     Connection is a precondition).
6. Cleanup: delete project; revoke permissions; restore Query
   name (rename it back to `PostgresCustomers` if Step 5
   succeeded).

## Scenarios

### Main flow — query-source lifecycle

1. **Open table from Query.** Run
   `Samples:PostgresCustomers` query (via `Browse > Platform >
   Functions > Queries` double-click OR
   `await grok.functions.eval('Samples:PostgresCustomers()')`).
   Verify the resulting table is loaded.
2. **Save project with Data Sync ON.** Save Project (ribbon
   button), name from Setup, Data Sync **ON**, OK. Cancel auto-
   share.
3. **Share project with second user (View-and-Use + Full).** Use
   `grok.dapi.permissions.grant`. Verify Sharing tab lists
   recipient.
   - Original-user assertion: project reopens; query result
     loads via Data Sync.
   - Recipient-side assertion (Helper 3 — deferred): second
     user opens; query result loads (recipient's session uses
     their own DB credentials via Spawner per platform's
     per-user query execution model).
4. **Rename external dependency — Query rename (github-3550 full
   repro).**
   - Locate the Query in `Browse > Platform > Functions >
     Queries > Samples > PostgresCustomers`.
   - Right-click the Query → **Rename** (if test account has
     edit permission on Samples namespace) OR rename via JS API:
     ```js
     const q = await grok.dapi.queries.find('Samples:PostgresCustomers');
     q.name = 'PostgresCustomers_renamed';
     await grok.dapi.queries.save(q);
     ```
   - **github-3550 invariant assertion:** the project's
     `relations.list` should now reference the renamed query
     (or the relation should be invalidated with an explicit
     error, not a silent null). Specifically: reopen the
     project; verify the table either loads (relation auto-
     resolved post-rename — happy path) OR fails with an
     explicit "Could not resolve query Samples:PostgresCustomers"
     error (graceful failure — not silent corruption).
   - Skip-with-logged-warning if the rename rejects on
     permissions (Samples namespace edit permission missing).
5. **Rename project itself.** Via JS API. Verify rename
   persists.
   - Original-user assertion: project opens under new name;
     query relation still resolves (post-rename, post-rename —
     i.e. both renames stack correctly).
   - Recipient-side assertion (Helper 3 — deferred): same.
6. **Cleanup.** Delete project. Revoke permissions. Rename
   the Query back to `PostgresCustomers` (if rename succeeded
   in Step 4) — environmental hygiene to avoid leaving
   `PostgresCustomers_renamed` for downstream tests.

### Expected results

- Save / reopen works for Query-sourced projects.
- Share + recipient-open works (recipient executes query
  under their own credentials).
- **github-3550 invariant:** Query rename does NOT silently
  break the project's relations — either auto-resolution
  (preferred) or explicit failure (acceptable). Silent null
  references are the bug.
- Project rename + Query rename stack: both renames persist.

## Notes

- **Origin: chain rev 3 proactive_lifecycle_specs[1]** with
  `bugs_reinforcing: [github-3550]` and
  `dep_lifecycle_ops_covered: [rename_external_dep,
  share_with_recipient_open]`. Authored in Phase A.
- **github-3550 full reproduction.** Step 4 walks the bug's
  exact reproduction path: external Query rename + project
  reopen + relation-resolution check. The cross-cutting bug
  spec `projects-github-3550-spec.ts` (chain rev 3
  bug_focused_candidates) targets the broader span (chain
  rev 3 spans `complex.md:Step 9` only, single-scenario);
  this scenario reinforces with full lifecycle coverage.
- **Env dependency: DB + Query + recipient DB-access.**
  This scenario CANNOT degrade gracefully without all 3 —
  if the Postgres connection is not provisioned, skip the
  entire scenario with a logged warning. Per
  `sa-2026-05-03-postgres-queries-public-data-substitution`,
  Samples:PostgresCustomers IS the canonical public-data
  substitute — but the underlying Postgres connection still
  needs provisioning.
- **`projects.api.relations.list` in sub_features.** Step 4
  asserts on `(await grok.dapi.projects.find(<id>))`'s
  relations list to verify rename impact. This sub_feature is
  a primary observation surface for github-3550.
- **Recipient must have DB connection access.** Per the
  platform's per-user query execution model, the recipient
  uses their own DB credentials via Spawner to re-execute the
  Query. If the recipient lacks access, the share path fails
  with a different bug (GROK-19403-style un-shared-deps —
  covered by `projects-lifecycle-script.md`, not this one).
- **UI coverage delegated.** All UI surfaces touched here
  (right-click rename, Save dialog, Sharing tab) are owned
  by `projects-ui-smoke.md`. JS API path used here.
- **Helper 3 deferral.** Same pattern.
- **Self-cleaning.** Step 6 deletes project + restores Query
  name. NOT in `deleting.md.depends_on`.
- **Sequencing within Wave 2C.** Third lifecycle scenario
  per Plan line 218: needs env Query commit; can be authored
  in any order with specs 4-6 once envs land.
