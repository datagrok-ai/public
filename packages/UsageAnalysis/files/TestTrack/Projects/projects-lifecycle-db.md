---
feature: projects
sub_features_covered:
  - projects.upload
  - projects.api.save
  - projects.api.files.sync
  - projects.add-relation
  - projects.shell.share-via-context-menu
  - projects.api.get-by-id
target_layer: playwright
coverage_type: regression
pyramid_layer: proactive-lifecycle
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/projects-lifecycle-db.md
migration_date: 2026-05-04
migration_report: projects-lifecycle-db-migration-report.md
related_bugs: []
---

# Projects — DB-source lifecycle (proactive)

Chained lifecycle for projects sourced from a **DB table** on
`System:Datagrok` (the platform's built-in metadata Postgres connection,
created by `ServiceConnectionsMigration` on every deploy — always
present). Exercises proactive coverage cell `source_class=db_table ×
dep_lifecycle_op=share_with_recipient_open` per chain rev 3
`proactive_lifecycle_specs[4]`. Pure proactive coverage — no GROK
ticket targets DB-source share+recipient-open today; this scenario
exists to catch regressions in DB-connection-permission gating
proactively.

UI coverage delegated to `projects-ui-smoke.md`.

Note: chain rev 3 names the spec `projects-lifecycle-db-table-spec.ts`
with `kebab(source_class.id) = "db-table"`. Per Olena's plan line 324
the canonical filename uses `db` shorthand (target file:
`projects-lifecycle-db.md` / `projects-lifecycle-db-spec.ts`) for
authoring brevity. This is a deviation from chain rev 3's strict
canonical naming convention — surfaced for B14 retro / chain rev 4
reconciliation.

## Setup

1. Authenticate as test user.
2. Project name: `lifecycle-db-${Date.now()}`.
3. Recipient placeholder: `<RECIPIENT_USERNAME_TBD>`.
4. Helper 3 dependency: `helpers.playwright.session.logoutAndLoginAs`
   (NOT YET REGISTERED).
5. Cleanup: delete project; revoke project permissions; delete the
   provisioned saved query (if Test 1 path) — all in `finally`.

## Scenarios

### Test 1 — DB query source via provisioned saved query

0. **Provision saved query on `System:Datagrok`.** Use
   `helpers/openers.ts:provisionSystemDatagrokQuery({nameStem:
   'lifecycle_db_query', sql: SYSTEM_DATAGROK_QUERIES.GROUPS_SAMPLE})`.
   The helper creates a Postgres query against the `groups` table in
   the platform's metadata DB and returns `{queryId, queryNqName,
   resolvedName, cleanup}`. The query is namespaced under the test
   user's login, so the user has full edit/rename/delete rights.
1. **Open table from saved query.** Use
   `helpers/openers.ts:openTableFromDbQuery(page,
   provisioned.queryNqName)`. Verify `df.tags['.script']` matches
   `<Var> = <queryNqName>()`.
2. **Save project with Data Sync ON** via `saveProjectWithProvenance`.
3. **Reopen** via `reopenAndAssertProvenance` — verify the query
   re-executes server-side and `.script` references the provisioned
   query's resolved name.
4. **Share project with second user (View-and-Use).** Use
   `grok.dapi.permissions.grant(project, recipient, false)`.
   Defensive skip if no second user exists (Helper 3 deferred).
5. **Cleanup.** Delete project + invoke `provisioned.cleanup()`.

### Test 2 — DB table source via ad-hoc DbQuery double-click

1. **Open `public.groups` ad-hoc.** Use
   `helpers/openers.ts:openTableFromDbTable(page, {connectionNqName:
   SYSTEM_DATAGROK_NQNAME, schemaName: 'public', tableName:
   'groups'})`. Verify `df.tags['.script']` matches
   `<Var> = DbQuery(System:Datagrok, "groups", ...)`.
2. **Save project with Data Sync ON** via `saveProjectWithProvenance`.
3. **Reopen** via `reopenAndAssertProvenance` — verify DbQuery
   re-runs server-side.
4. **Cleanup.** Delete project. (No external query to clean up —
   ad-hoc DbQuery is part of project relations only.)

### Original scope notes (deferred)

- **Rename external dependency — N/A for db_table per chain rev 3.**
  Per `proactive_lifecycle_specs[4].dep_lifecycle_ops_covered:
  [share_with_recipient_open]` (NOT `rename_external_dep`).
  Atlas rationale: rename of an external DB table renames at
  the connection layer, not the project-relation layer; the
  project's relation refers to the table by qualified name and
  is not affected by table-name churn. This step is a no-op.
- **Rename project itself.** Via JS API. Verify rename
  persists; recipient still opens under new name (Helper 3 —
  deferred).
- **Recipient-side assertion (Helper 3 — deferred):** logout +
  login as recipient; open shared project; verify the DB table
  loads correctly under recipient's session.

### Expected results

- Save / reopen works for DB-sourced projects.
- Share + recipient-open works when project is shared (recipient
  inherits access to `System:Datagrok` via `allUsers` — built-in
  permission model, no per-user DB credential dance required).
- Project rename persists; share survives rename.

## Notes

- **Self-contained source provisioning.** Both tests create their
  own DB-source prerequisites within the spec — Test 1 provisions
  a saved query against `System:Datagrok` via
  `helpers/openers.ts:provisionSystemDatagrokQuery`; Test 2 uses
  an ad-hoc DbQuery on the same connection. No external package
  (Samples) or env-provisioned DB connection is required.
- **`System:Datagrok` rationale.** Built-in read-only Postgres
  connection to the platform's metadata DB; created by
  `ServiceConnectionsMigration` on every deploy and shared with
  `allUsers` (or `admins` on public). The `groups` table is part
  of the metadata schema and always non-empty.
- **Origin: chain rev 3 proactive_lifecycle_specs[4]** with
  `bugs_reinforcing: []` and
  `dep_lifecycle_ops_covered: [share_with_recipient_open]`.
  Pure proactive — no GROK ticket targets this cell. Authored
  in Phase A.
- **Filename deviation from canonical `db-table` →
  `db`.** Per Plan line 324, authoring shorthand. Chain rev
  4 should reconcile (either rename canonical to `db` or
  rename file to `projects-lifecycle-db-table.md`).
- **No external rename for db_table.** Per chain rev 3,
  external DB table rename is a connection-layer concern,
  not project-relation; not in scope for this entry.
- **UI coverage delegated.** All UI surfaces are owned by
  `projects-ui-smoke.md`. JS API path used here.
- **Helper 3 deferral.** Recipient-side assertions blocked.
- **Self-cleaning.** Cleanup deletes project + invokes
  `provisioned.cleanup()` for the saved query.
