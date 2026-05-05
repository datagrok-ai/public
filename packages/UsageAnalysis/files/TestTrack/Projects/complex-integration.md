---
feature: projects
sub_features_covered:
  - projects.upload
  - projects.api.save
  - projects.api.files.sync
  - projects.api.namespaces
  - projects.add-relation
  - projects.add-link
target_layer: playwright
coverage_type: regression
pyramid_layer: integration
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: decomposed
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/complex.md
migration_date: 2026-05-04
migration_report: complex-integration-migration-report.md
related_bugs: []
---

# Complex — Multi-source integration

Source-agnostic op (Wave 2B). The only intentionally-long Wave 2B
spec — value is multi-source co-existence verification. Decomposition
of `complex.md` Step 1 in full + Step 2 (Save with Sync ON + reopen)
per Plan line 170.

**8 sources** opened in one workspace + Pivot + Aggregate + Join + Clone
+ Save Sync ON + reopen. Tests that the platform correctly handles
projects with heterogeneous source coexistence.

UI coverage delegated to `projects-ui-smoke.md`.

## Setup

1. Authenticate as test user.
2. Project name: `integration-test-${Date.now()}`.
3. **Environment dependencies (multiple):**
   - File share: `System:AppData/Chem/tests/spgi-100.csv`,
     `System:DemoFiles/demog.csv`.
   - Space: inline-created `test-projects-demo-integration-*`
     populated with `demog.csv` (Spaces prelude pattern).
   - DB: `Postgres:NorthwindTest:public:orders` (env-provisioned).
   - Query: `Samples:PostgresCustomers` (env-provisioned).
   - Script: `Samples:Cars` (env-provisioned).
   - Skip-with-logged-warning if any of the above is not
     present. Per
     `sa-2026-05-03-postgres-queries-public-data-substitution`,
     `sa-2026-05-03-spgi-file-public-data-substitution`,
     `sa-2026-05-03-spaces-inline-prelude-pattern` patterns.
4. Cleanup: delete project; delete inline Space.

## Scenarios

### Main flow — multi-source integration

1. **Open 8 tables from heterogeneous sources.**
   - File share: `System:AppData/Chem/tests/spgi-100.csv`.
   - Space: `test-projects-demo-integration-* > demog.csv`
     (prelude-created).
   - DB: `Browse > Databases > Postgres > NorthwindTest > Schema >
     public > orders` (double-click).
   - Query: run `Samples:PostgresCustomers`.
   - Script: run `Samples:Cars`.
   - Pivot: configure Pivot Table on the demog table > **Add**.
   - Aggregate: configure Aggregate Rows on the cars
     (Samples:Cars) table > **Add**.
   - Join: Join the spgi-100 table with the orders table on
     a common-name key column (or any plausible join
     condition); produce as new tab.
   - Clone: right-click any open table > **Clone** to produce
     an independent copy.
   - Verify `grok.shell.tables.length >= 8` after all 8 sources
     are open.
2. **Save project with Data Sync ON (the integration save).**
   Save Project, name from Setup, Data Sync **ON** for every
   applicable table (some tables — derived Pivot/Aggregate/Join,
   Clone — do not have Data Sync configurable at that level;
   their parents do). Click **OK**. Cancel auto-share.
3. **Verify save persisted with all 8 tables.** Verify
   `(await grok.dapi.projects.find(<id>)).relations.length` is
   the expected count (depends on parent vs derived counting
   semantics; assert `>= 5` for parent tables at minimum:
   spgi, demog, orders, postgres-customers, cars).
4. **Close all views and reopen the integration project.** Close
   via `grok.shell.closeAll()`. Reopen via
   `Browse > Dashboards`. Verify all parent tables load (5+);
   verify derivations (Pivot, Aggregate, Join, Clone) restore
   correctly. **Multi-source co-existence assertion:** all 8
   tabs are present in the workspace; no missing-source errors.
5. **Cleanup.** Delete the project. Delete the inline Space.
   Optionally restore environment state (no rename ops in
   this scenario, so no restoration needed).

### Expected results

- 8 heterogeneous sources can coexist in one project.
- Save Project with Data Sync ON persists all sources +
  derivations.
- Reopen restores all 8 tables / tabs without errors.
- No source-class interaction conflicts (e.g. Spaces +
  Postgres + Files all coexist; Pivot + Aggregate + Join all
  derive correctly; Clone is independent).

## Notes

- **Origin: `complex.md` Step 1 in full + Step 2 (Wave 2B
  re-framing per Plan line 170).** Plan explicitly notes:
  "The only intentionally-long spec — value is multi-source
  co-existence." Re-framing of Migrator's Split A as a
  standalone integration test.
- **Why this is NOT bug-focused.** None of the GROK bugs
  target multi-source co-existence specifically. Each bug
  is source-specific or op-specific. This scenario covers
  the gap: regressions in multi-source coordination would
  be caught here proactively.
- **`projects.add-link` in sub_features.** Cloning
  exercises `addLink` (vs `addRelation` for fresh sources).
- **Long runtime acknowledged.** Plan explicitly accepts
  this is the longest Wave 2B spec (~10-15 min in CI).
  Value is in multi-source coverage; signal/noise improves
  vs running 8 separate single-source specs.
- **UI coverage delegated.** No UI surface beyond Save
  Project. Pivot/Aggregate/Join/Clone UI surfaces are owned
  by `uploading.md` source-matrix Cases 8/9.
- **Env-dependent (skip if missing).** All 5 env
  dependencies (File, Space, DB, Query, Script) must be
  present — defensive skip-with-logged-warning otherwise.
- **No `related_bugs`.** Pure proactive multi-source
  coverage.
- **Self-cleaning.** Step 5 deletes project + Space.
- **Sequencing within Wave 2B.** Last in Wave 2B per Plan
  execution order — heaviest, requires all 5 env
  provisioning items.
