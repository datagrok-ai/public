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

Chained lifecycle for projects sourced from a **DB table**
(`Postgres:NorthwindTest:public:orders` via Browse > Databases double-
click). Exercises proactive coverage cell `source_class=db_table ×
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
5. **Environment dependencies:**
   - DB connection: `Postgres:NorthwindTest` provisioned for `qa-pw`.
   - Recipient must have **per-user permissions on the Connection**
     (this is the proactive coverage axis — DB connection share +
     per-user execution model).
   - Skip-with-logged-warning if any of the above is not present
     (per `sa-2026-05-03-postgres-queries-public-data-substitution`
     pattern).
6. Cleanup: delete project; revoke project permissions; revoke
   connection permissions if granted as part of the test.

## Scenarios

### Main flow — DB-source lifecycle (proactive)

1. **Open table from DB via double-click.** Navigate
   `Browse > Databases > Postgres > NorthwindTest > Schema >
   public > orders`. Double-click. Verify the table loads.
2. **Save project with Data Sync ON.** Save Project, name from
   Setup, Data Sync **ON**, OK. Cancel auto-share.
3. **Share project with second user (View-and-Use + Full).** Use
   `grok.dapi.permissions.grant`. Verify Sharing tab.
   - Verify the DB connection IS also shared with recipient (or
     grant connection-level permission as part of this step:
     `grok.dapi.permissions.grant(connection, recipient, false)`).
   - Original-user assertion: project reopens; DB table loads via
     Data Sync.
   - **Recipient-side assertion (Helper 3 — deferred):** logout +
     login as recipient; open shared project; verify the DB
     table loads correctly under recipient's session (the
     recipient's per-user DB credentials via Spawner re-execute
     the query against `public.orders`).
4. **Rename external dependency — N/A for db_table per chain rev 3.**
   Per `proactive_lifecycle_specs[4].dep_lifecycle_ops_covered:
   [share_with_recipient_open]` (NOT `rename_external_dep`).
   Atlas rationale: rename of an external DB table renames at
   the connection layer, not the project-relation layer; the
   project's relation refers to the table by qualified name and
   is not affected by table-name churn. This step is a no-op.
5. **Rename project itself.** Via JS API. Verify rename
   persists; recipient still opens under new name (Helper 3 —
   deferred).
6. **Cleanup.** Delete project. Revoke project permissions.
   Revoke connection permission grant (if granted in Step 3).

### Expected results

- Save / reopen works for DB-sourced projects.
- Share + recipient-open works when both project AND
  connection are shared.
- Project rename persists; share survives rename.
- Per-user query execution model works under recipient's
  session.

## Notes

- **Origin: chain rev 3 proactive_lifecycle_specs[4]** with
  `bugs_reinforcing: []` and
  `dep_lifecycle_ops_covered: [share_with_recipient_open]`.
  Pure proactive — no GROK ticket targets this cell. Authored
  in Phase A.
- **bugs_reinforcing: [].** Proactive coverage rationale per
  PROJECTS-SPLIT-COMPLETION-PLAN bug-traceability line 311:
  the right policy is per-source bundle covers all plausible
  interaction cells regardless of whether a GROK ticket exists.
  If this cell never has a regression, the cost is bounded
  (~7 min in CI).
- **Connection-share is the proactive surface.** This
  scenario specifically tests the coordination of project-
  share + connection-share + per-user query execution. The
  connection-share step (Step 3) is the regression candidate
  — if any future connection-share refactor breaks this path,
  THIS scenario catches it before a customer hits it.
- **Filename deviation from canonical `db-table` →
  `db`.** Per Plan line 324, authoring shorthand. Chain rev
  4 should reconcile (either rename canonical to `db` or
  rename file to `projects-lifecycle-db-table.md`).
- **No external rename for db_table.** Per chain rev 3,
  external DB table rename is a connection-layer concern,
  not project-relation; not in scope for this entry.
- **`projects.api.namespaces` NOT in sub_features.** DB
  connections live under `Postgres:` namespace but the
  project's relation is to the connection + table, not the
  namespace itself. Atlas confirms.
- **UI coverage delegated.** All UI surfaces are owned by
  `projects-ui-smoke.md`. JS API path used here.
- **Helper 3 deferral.** Recipient-side assertions blocked.
- **Self-cleaning.** Step 6 deletes project + revokes
  permissions.
- **Sequencing within Wave 2C.** Sixth lifecycle scenario per
  Plan line 221: needs DB connection per-user perms commit.
  Sequencing flexible with specs 3, 4, 5.
